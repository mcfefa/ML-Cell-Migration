#@ File(label='Choose a directory of videos', style='directory') import_dir
#@ String(label='File types', value='tif') file_types
#@ File (label = "Cellpose .exe", style = "file") cpExe
#@ File (label = "Cellpose model", style = "file") cpModel
#@ String (visibility=MESSAGE, value="Leave blank to analyze all files without filtering for keywords", required=false) msg
#@ String (label='Filter', value='') filters
#@ String (visibility=MESSAGE, value="Leave the following fields empty or as 0 if time slice units associated with each image stack is accurate.", required=false) timeMsg
#@ String(Label="Time units (sec, msec, usec, min, etc.)", value="") timeUnit
#@ Double(Label="Frame interval (time between two frames in above time units)", value=0.0) frameInterval
#@ Boolean(label='Process Subfolders', value=False) process_sub
#@ String (visibility=MESSAGE, value="Input comma-separated values (with no spaces) for the following tracking input variable ranges", required=false) msg1
#@ String (visibility=MESSAGE, value="Max Frame Gap values are integers", required=false) msg2
#@ String (Label="Range of Max Frame Gap values to test", value='3,4,5') maxFrameGapString
#@ String (visibility=MESSAGE, value="Alt Link Cost, Link Max, and Gap Close Max values are doubles", required=false) msg3
#@ String(Label="Range of Alternate Link Cost Factor values to test", value='1.05') altLinkCostString
#@ String(Label="Range of Linking Max Distance values to test", value='10.0,20.0,30.0') linkMaxString
#@ String(Label="Range of Gap Closing Max Distance values to test", value='10.0,20.0,30.0') gapCloseMaxString
#@ Double(Label="Weight for number of tracks kept", value=0.2) tracksKeptWeight
#@ Double(Label="Weight for average track quality", value=0.4) avgTrackQualityWeight
#@ Double(Label="Weight for average track duration", value=0.4) avgTrackDurationWeight

from fiji.plugin.trackmate.visualization.hyperstack import HyperStackDisplayer
from fiji.plugin.trackmate.io import TmXmlWriter
from fiji.plugin.trackmate import TrackMate, Logger, Dimension, Settings, SelectionModel, Spot, TrackModel, SpotCollection, Model
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettingsIO, DisplaySettings
from fiji.plugin.trackmate.features.track import TrackAnalyzer, TrackSpotQualityFeatureAnalyzer, TrackDurationAnalyzer, AbstractTrackAnalyzer
from java.io import File
#from fiji.util.gui import GenericDialogPlus
from ij import IJ, WindowManager
from ij.measure import ResultsTable
from ij.text import TextWindow
from fiji.plugin.trackmate.visualization.table import TrackTableView, AllSpotsTableView
from fiji.plugin.trackmate.tracking import SpotTrackerFactory
from fiji.plugin.trackmate.tracking.jaqaman import LAPUtils, SparseLAPTrackerFactory
from fiji.plugin.trackmate.providers import SpotAnalyzerProvider, EdgeAnalyzerProvider, TrackAnalyzerProvider
from fiji.plugin.trackmate.cellpose.CellposeSettings import PretrainedModel
from fiji.plugin.trackmate.cellpose import CellposeDetectorFactory
from fiji.plugin.trackmate.action import CaptureOverlayAction
from fiji.plugin.trackmate.util import TMUtils, LogRecorder
import fiji.plugin.trackmate.action.ExportTracksToXML as ExportTracksToXML
import math, os, sys

# We have to do the following to avoid errors with UTF8 chars generated in 
# TrackMate that will mess with our Fiji Jython.
reload(sys)
sys.setdefaultencoding('utf-8')
    
# adapted from https://imagej.net/scripting/jython/
def batch_open_images(path, file_type=None, name_filter=None, process_sub=False):
    '''Open all files in the given folder.
    :param path: The path from were to open the images. String and java.io.File are allowed.
    :param file_type: Only accept files with the given extension (default: None).
    :param name_filter: Only accept files that contain the given string (default: None).
    :param recursive: Process subfolders (default: False).
    '''
    # Converting a File object to a string.
    if isinstance(path, File):
        path = path.getAbsolutePath()
 
    def check_type(string):
        '''This function is used to check the file type.
        It is possible to use a single string or a list/tuple of strings as filter.
        This function can access the variables of the surrounding function.
        :param string: The filename to perform the check on.
        '''
        if file_type:
            # The first branch is used if file_type is a list or a tuple.
            if isinstance(file_type, (list, tuple)):
                for file_type_ in file_type:
                    if string.endswith(file_type_):
                        # Exit the function with True.
                        return True
                    else:
                        # Next iteration of the for loop.
                        continue
            # The second branch is used if file_type is a string.
            elif isinstance(file_type, string):
                if string.endswith(file_type):
                    return True
                else:
                    return False
            return False
        # Accept all files if file_type is None.
        else:
            return True
 
    def check_filter(string):
        '''This function is used to check for a given filter.
        It is possible to use a single string or a list/tuple of strings as filter.
        This function can access the variables of the surrounding function.
        :param string: The filename to perform the filtering on.
        '''
        if name_filter:
            # The first branch is used if name_filter is a list or a tuple.
            if isinstance(name_filter, (list, tuple)):
                for name_filter_ in name_filter:
                    if name_filter_ in string:
                        # Exit the function with True.
                        return True
                    else:
                        # Next iteration of the for loop.
                        continue
            # The second branch is used if name_filter is a string.
            elif isinstance(name_filter, string):
                if name_filter in string:
                    return True
                else:
                    return False
            return False
        else:
        # Accept all files if name_filter is None.
            return True
 
    # We collect all files to open in a list.
    path_to_images = []
    # Replacing some abbreviations (e.g. $HOME on Linux).
    path = os.path.expanduser(path)
    path = os.path.expandvars(path)
    # If we don't want to process subfolders, we can use os.listdir().
    if not process_sub:
        for file_name in os.listdir(path):
            full_path = os.path.join(path, file_name)
            if os.path.isfile(full_path):
                if check_type(file_name):
                    if check_filter(file_name):
                        path_to_images.append(full_path)
    # To process subfolders os.walk() is used.
    else:
        # os.walk() is iterable.
        # Each iteration of the for loop processes a different directory.
        # the first return value represents the current directory.
        # The second return value is a list of included directories.
        # The third return value is a list of included files.
        for directory, dir_names, file_names in os.walk(path):
            # We are only interested in files.
            for file_name in file_names:
                # The list contains only the file names.
                # The full path needs to be reconstructed.
                full_path = os.path.join(directory, file_name)
                # Both checks are performed to filter the files.
                if check_type(file_name):
                    if check_filter(file_name):
                        # Add the file to the list of images to open.
                        path_to_images.append(full_path)
    # Output full path of each image
    return path_to_images
 
def split_string(input_string):
    '''Split a string to a list and strip it
    :param input_string: A string that contains semicolons as separators.
    '''
    string_splitted = input_string.split(';')
    # Remove whitespace at the beginning and end of each string
    strings_striped = [string.strip() for string in string_splitted]
    return strings_striped

def displayTracks(model, imp):
    sm = SelectionModel(model)
    ds = DisplaySettingsIO.readUserDefault()
    ds.spotDisplayedAsRoi = True
    displayer =  HyperStackDisplayer( model, sm, imp, ds )
    displayer.render()
    displayer.refresh()

def findUnderSD(model, settings, featureName):
    trackIDsUnder = []
    feature = 0
    fm = model.getFeatureModel()
    numTracks = model.getTrackModel().trackIDs(True).size()
    for id in model.getTrackModel().trackIDs(True):
        featureVal = fm.getTrackFeature(id, featureName)
        feature = feature + featureVal
    meanFeature = feature/numTracks
    SS = 0
    for id in model.getTrackModel().trackIDs(True):
        featureVal = fm.getTrackFeature(id, featureName)
        SS = SS + (featureVal-meanFeature)*(featureVal-meanFeature)
    sd = math.sqrt(abs(SS)/(numTracks-1))
    for id in model.getTrackModel().trackIDs(True):
        featureVal = fm.getTrackFeature(id, featureName)
        if(featureVal < meanFeature-2*sd):  
            trackIDsUnder.append(id)
    return trackIDsUnder
    
def findTrackingQualityBatch(model, settings, maxFrameGap, altLinkCost, linkMax, gapCloseMax):
    # Configure trackerFactory, determine tracking settings
    settings.trackerSettings["MAX_FRAME_GAP"] = maxFrameGap
    settings.trackerSettings["ALTERNATIVE_LINKING_COST_FACTOR"] = altLinkCost
    settings.trackerSettings["LINKING_MAX_DISTANCE"] = linkMax
    settings.trackerSettings["GAP_CLOSING_MAX_DISTANCE"] = gapCloseMax
    
    trackmate = TrackMate(model, settings)
    
    print('Begin Tracking')
    model.getLogger().log('Begin Tracking')
    ok = trackmate.execTracking()
    if( ok ):
        print('TrackMate completed successfully.' )
        model.getLogger().log('TrackMate completed successfully.')
    else:
        print(str(trackmate.getErrorMessage()))
        model.getLogger().log(str(trackmate.getErrorMessage()))
    print('Found %d spots in %d tracks.' % ( model.getSpots().getNSpots( True ) , model.getTrackModel().nTracks( True )))
    model.getLogger().log('Found %d spots in %d tracks.' % ( model.getSpots().getNSpots( True ) , model.getTrackModel().nTracks( True )))
    print('Computing track features')
    trackmate.computeTrackFeatures(True)
    model.getLogger().log("FINISHED track features")
    
    tracksUnderDisp = findUnderSD(model, settings, "TRACK_DISPLACEMENT")
    tracksUnderQuality = findUnderSD(model, settings, "TRACK_MEAN_QUALITY")
    tracksUnderLen = findUnderSD(model, settings, "TRACK_DURATION")
    
    tracksUnder = list(set(tracksUnderDisp + tracksUnderQuality + tracksUnderLen))
    for id in tracksUnder:
        count = 0
        if id in tracksUnderDisp:
            count = count + 1
        if id in tracksUnderQuality:
            count = count + 1
        if id in tracksUnderLen:
            count = count + 1
        if count >= 2:
            model.setTrackVisibility(id, False)
    
    tracksKept = 0
    avgTrackQuality = 0
    avgTrackDuration = 0
    for id in model.getTrackModel().trackIDs(True):
        tracksKept = tracksKept + 1
        avgTrackQuality = avgTrackQuality + model.getFeatureModel().getTrackFeature(id, 'TRACK_MEAN_QUALITY')
        avgTrackDuration = avgTrackDuration + model.getFeatureModel().getTrackFeature(id, 'TRACK_DURATION')
    avgTrackQuality = avgTrackQuality/tracksKept
    avgTrackDuration = avgTrackDuration/tracksKept
    model.getLogger().log('# of tracks kept: ' + str(tracksKept))
    model.getLogger().log('Average track quality: ' + str(avgTrackQuality))
    model.getLogger().log('Average track duration: ' + str(avgTrackDuration) + ' ' + timeUnit)
    
    for id in model.getTrackModel().trackIDs(False):
        model.setTrackVisibility(id, True)
    print([maxFrameGap, altLinkCost, linkMax, gapCloseMax])
    model.getLogger().log("[" + str(maxFrameGap) + ", " + str(altLinkCost) + ", " + str(linkMax) + ", " + str(gapCloseMax) + "]")
    trackingQualityTable.addRow()
    trackingQualityTable.addValue(0, maxFrameGap)
    trackingQualityTable.addValue(1, altLinkCost)
    trackingQualityTable.addValue(2, linkMax)
    trackingQualityTable.addValue(3, gapCloseMax)
    trackingQualityTable.addValue(4, tracksKept)
    trackingQualityTable.addValue(5, avgTrackQuality)
    trackingQualityTable.addValue(6, avgTrackDuration)
    trackingQualityTable.show("Tracking_Quality")
    return trackmate

def findMiddle(model):
    bins = 25
    colCount = []
    for i in range(1, bins + 1):
        count = 0
        for spot in model.getSpots().iterable(1, True):
            X = spot.getFeature("POSITION_X")
            if (X<=i*imp.getHeight()/bins) and (X>=(i-1)*imp.getHeight()/bins):
                count = count + 1
        colCount.append(count)
    zeros = [index for index, value in enumerate(colCount) if value <= 10] 
    middle = (imp.getHeight()/bins*(zeros[0]-1)+imp.getHeight()/bins*zeros[len(zeros)-1])/2
    return(middle)

# Assigns whether the first cell in a track is to the left or right of the middle
def pathDir(model):
    middle = findMiddle(model)
    for id in model.getTrackModel().trackIDs(True):
        track = model.getTrackModel().trackSpots(id)
        spot = track.iterator().next()
        if spot.getFeature("POSITION_X") < middle:
            spot.putFeature('LR', 0)
        else:
            spot.putFeature('LR', 1)
        
# Calculates directed migration and directed migration ratio for each visible track
def calcDirMig(model):
    fm = model.getFeatureModel()
    for id in model.getTrackModel().trackIDs(True):
        # Get all spots in current track
        track = model.getTrackModel().trackSpots(id).iterator()
        # Find whole track displacement
        disp = fm.getTrackFeature(id, "TRACK_DISPLACEMENT")
        # Create vector from first spot position to last spot position
        # Find the position of the last spot, center axis on starting cell
        firstCell = track.next()
        cellPos = [0, 0]
        while(track.hasNext()):
            lastCell = track.next()
        cellPos[0] = lastCell.diffTo(firstCell, "POSITION_X")
        cellPos[1] = lastCell.diffTo(firstCell, "POSITION_Y")
        #print("Method: X: %.5f, Y: %.5f" % (lastCell.diffTo(firstCell, "POSITION_X"), lastCell.diffTo(firstCell, "POSITION_Y")))
        # Find direction vector from first spot LR
        dirVec = [0, 0]
        if firstCell.getFeature("LR") == 0:
            dirVec[0] = 1
        else:
            dirVec[0] = -1
        # Calculate projection vector and magnitude if in same direction as dirVec
        dirMigVec = [0, 0]
        dirMig = 0
        dotProd = cellPos[0]*dirVec[0] + cellPos[1]*dirVec[1]
        dirVecMag = math.sqrt(dirVec[0]**2+dirVec[1]**2)
        dirMigVec = [(dim*dotProd)/dirVecMag**2 for dim in dirVec]
        dirMig = math.sqrt(dirMigVec[0]**2+dirMigVec[1]**2)
        if ((cellPos[0]>0) and (dirVec[0]>0)):
            if dirMig < 0:
                dirMig = dirMig * -1
        elif ((cellPos[0]<0) and (dirVec[0]>0)):
            if dirMig > 0:
                dirMig = dirMig * -1
        elif ((cellPos[0]<0) and (dirVec[0]<0)):
            if dirMig < 0:
                dirMig = dirMig * -1
        elif ((cellPos[0]>0) and (dirVec[0]<0)):
            if dirMig > 0:
                dirMig = dirMig * -1
        fm.putTrackFeature(id, "DirectedMigration", dirMig)
        #print(fm.getTrackFeature(id, "DirectedMigration"))
        # Calculate directed migration ratio
        fm.putTrackFeature(id, "DirectedMigrationRatio", dirMig/disp)
        #print(fm.getTrackFeature(id, "DirectedMigrationRatio"))    
                
def outputTrackData(model, sm, ds, wkdir):
    # Save spot, edge, and track statistics
    trackTableView = TrackTableView(model, sm, ds, fullPath)
    trackTableView.getSpotTable().getPanel()
    trackTableView.getSpotTable().exportToCsv(File(os.path.join(dirName,fileName + "_spotData.csv")))
    trackTableView.getEdgeTable().exportToCsv(File(os.path.join(dirName,fileName + "_edgeData.csv")))
    trackTableView.getTrackTable().exportToCsv(File(os.path.join(dirName,fileName + "_trackData.csv")))
    IJ.log("Spot, edge, and track data saved to " + dirName + " for " + fileName)

# adapted from https://forum.image.sc/t/importing-cellpose-via-scyjava-for-trackmate-automation-in-python-3/84176/6
def cpProcess(imp, cpExe, cpModel, bestTrackerVal):
    cal = imp.getCalibration()
    
    model = Model()
    model.setPhysicalUnits(cal.getUnit(), cal.getTimeUnit())
    logger = LogRecorder(Logger.VOID_LOGGER)
    
    model.setLogger(Logger.IJ_LOGGER)
    settings = Settings(imp)
    
    # Configure Cellpose detectorFactory
    settings.detectorFactory = CellposeDetectorFactory()
    settings.detectorSettings = {
        'TARGET_CHANNEL' : 0,
        'OPTIONAL_CHANNEL_2' : 0,
        'CELLPOSE_PYTHON_FILEPATH' : cpExe,
        'CELLPOSE_MODEL' : PretrainedModel.CUSTOM ,   
        'CELLPOSE_MODEL_FILEPATH' : cpModel,
        'CELL_DIAMETER' : 0.0,
        'USE_GPU' : True,
        'SIMPLIFY_CONTOURS' : True,
    }   
    
    settings.trackerFactory = SparseLAPTrackerFactory()
    settings.trackerSettings = settings.trackerFactory.getDefaultSettings()
    settings.trackerSettings["ALLOW_GAP_CLOSING"] = True
    settings.addTrackAnalyzer(TrackSpotQualityFeatureAnalyzer())
    settings.addTrackAnalyzer(TrackDurationAnalyzer())
    
    spotAnalyzerProvider = SpotAnalyzerProvider(1)
    for key in spotAnalyzerProvider.getKeys():
        #print(key)
        settings.addSpotAnalyzerFactory(spotAnalyzerProvider.getFactory(key))
    edgeAnalyzerProvider = EdgeAnalyzerProvider()
    for  key in edgeAnalyzerProvider.getKeys():
        #print(key)
        settings.addEdgeAnalyzer(edgeAnalyzerProvider.getFactory(key))  
        
    global trackmate
    trackmate = TrackMate(model, settings)
    
    # Detection
    print('Spot detection')
    model.getLogger().log("Executing spot detection")
    trackmate.execDetection()
    model.getLogger().log("FINISHED spot detection")
    # Initial filtering, not in use, permanently deletes spots, requires re-detection to recover
    #print('Spot initial filtering')
    #ok = ok and trackmate.execInitialSpotFiltering()
    #model.getLogger().log("FINISHED spot initial filtering")
    # Compute spot features.
    print('Computing spot features')
    model.getLogger().log("Executing spot features")
    trackmate.computeSpotFeatures(True) 
    model.getLogger().log("FINISHED spot features")
    # Post filte    ring, not in use
    print('Filtering spots')
    trackmate.execSpotFiltering(True) 
    
    # Determine tracker settings and instantiate trackmate
    # Test all different combinations of input variables and measure tracking quality
    if bestTrackerVal[0] == 0:
        for maxFrameGapVal in maxFrameGap:
            for altLinkCostVal in altLinkCost:
                for linkMaxVal in linkMax:
                    for gapCloseMaxVal in gapCloseMax:
                        trackmate = findTrackingQualityBatch(model, settings, maxFrameGapVal, altLinkCostVal, linkMaxVal, gapCloseMaxVal)
                        displayTracks(model, imp)
                        
        colToCalc = trackingQualityTable.getHeadings()[4:7]
        if trackingQualityTable.size()>1:
            for col in colToCalc:
                col = str(col)
                #print("CURRENT COLUMN: " + col)
                mean = 0
                param = 0
                for i in range(trackingQualityTable.size()):
                    paramVal = trackingQualityTable.getValueAsDouble(trackingQualityTable.getColumnIndex(col), i)
                    param = param + paramVal
                mean = param/trackingQualityTable.size()
                SS = 0
                for i in range(trackingQualityTable.size()):
                    paramVal = trackingQualityTable.getValueAsDouble(trackingQualityTable.getColumnIndex(col), i)
                    SS = SS + (paramVal-mean)*(paramVal-mean)
                sd = math.sqrt(abs(SS)/(trackingQualityTable.size()))
                for row in range(trackingQualityTable.size()):
                    trackingQualityTable.setValue(trackingQualityTable.getColumnIndex(col)+3, row, (trackingQualityTable.getValueAsDouble(trackingQualityTable.getColumnIndex(col), row)-mean)/sd)
        trackingQualityTable.show("Tracking_Quality")

        tableLastCol = trackingQualityTable.getLastColumn()
        #print("PRINTING TO COL AFTER THIS: "+str(tableLastCol))
        bestTrackingQuality = float('-inf')
        for row in range(trackingQualityTable.size()):
                trackingQuality = tracksKeptWeight*trackingQualityTable.getValueAsDouble(tableLastCol-2, row) + avgTrackQualityWeight*trackingQualityTable.getValueAsDouble(tableLastCol-1, row) + avgTrackDurationWeight*trackingQualityTable.getValueAsDouble(tableLastCol, row)
                trackingQualityTable.setValue(tableLastCol+1, row, trackingQuality)
                if trackingQuality > bestTrackingQuality:
                        bestTrackerVal = [trackingQualityTable.getValueAsDouble(0, row), trackingQualityTable.getValueAsDouble(1, row), trackingQualityTable.getValueAsDouble(2, row), trackingQualityTable.getValueAsDouble(3, row)]
                        bestTrackingQuality = trackingQuality
    
    print("Begin tracking with optimized input variables on " + fileName)
    print("Max frame gap: " + str(int(bestTrackerVal[0])))
    print("Alternative linking cost factor: " + str(bestTrackerVal[1]))
    print("Linking max distance: " + str(bestTrackerVal[2]))
    print("Gap closing max distance: " + str(bestTrackerVal[3]))
    model.getLogger().log("Begin tracking with optimized input variables on " + fileName)
    model.getLogger().log("Max frame gap: " + str(int(bestTrackerVal[0])))
    model.getLogger().log("Alternative linking cost factor: " + str(bestTrackerVal[1]))
    model.getLogger().log("Linking max distance: " + str(bestTrackerVal[2]))
    model.getLogger().log("Gap closing max distance: " + str(bestTrackerVal[3]))
    
    settings.trackerSettings["MAX_FRAME_GAP"] = int(bestTrackerVal[0])
    settings.trackerSettings["ALTERNATIVE_LINKING_COST_FACTOR"] = bestTrackerVal[1]
    settings.trackerSettings["LINKING_MAX_DISTANCE"] = bestTrackerVal[2]
    settings.trackerSettings["GAP_CLOSING_MAX_DISTANCE"] = bestTrackerVal[3]
    
    ok = trackmate.execTracking()
    if( ok ):
        print('TrackMate completed successfully.')
        model.getLogger().log('TrackMate completed successfully.')
    else:
        print(str(trackmate.getErrorMessage()))
        model.getLogger().log(str(trackmate.getErrorMessage()))
    print( 'Found %d spots in %d tracks.' % ( model.getSpots().getNSpots( True ) , model.getTrackModel().nTracks( True ) ) )
    model.getLogger().log('Found %d spots in %d tracks.' % ( model.getSpots().getNSpots( True ) , model.getTrackModel().nTracks( True )))
    
    settings.addAllAnalyzers()
    print('Computing spot features')
    trackmate.computeSpotFeatures(True)
    model.getLogger().log("FINISHED spot features")
    print('Computing edge features')
    trackmate.computeEdgeFeatures(True)
    trackmate.computeTrackFeatures(True)
    
    print("Calculating additional directed migration features")
    model.getLogger().log("Calculating additional directed migration features")
    # Altered from https://github.com/TASBE/TASBEImageAnalytics/blob/5d6fc1a64b4c17e263451fa4252c94dc86193d14/scripts/cellStatsTracking.py
    spotFeat = ["LR"]
    trackFeat = ["DirectedMigration", "DirectedMigrationRatio"]                                      
    
    featureNames = {}
    featureShortNames = {}
    featureDimensions = {}
    isInt = {}
    for feature in spotFeat:
        featureNames[feature] = feature
        featureShortNames[feature] = feature
        featureDimensions[feature] = Dimension.STRING
        isInt[feature] = False
    
    model.getFeatureModel().declareSpotFeatures(spotFeat, featureNames, featureShortNames, featureDimensions, isInt)
    
    featureNames = {}
    featureShortNames = {}
    featureDimensions = {}
    isInt = {}
    for feature in trackFeat:
        featureNames[feature] = feature
        featureShortNames[feature] = feature
        featureDimensions[feature] = Dimension.STRING
        isInt[feature] = False
        
    model.getFeatureModel().declareTrackFeatures(trackFeat, featureNames, featureShortNames, featureDimensions, isInt)

    pathDir(model)
    calcDirMig(model)
    
    # Save results
    sm = SelectionModel( model )
    outputTrackData(model, SelectionModel(model), DisplaySettingsIO.readUserDefault(), dirName)
    saveFile = TMUtils.proposeTrackMateSaveFile( settings, logger )

    writer = TmXmlWriter( saveFile, logger )
    writer.appendLog( logger.toString() )
    writer.appendModel( trackmate.getModel() )
    writer.appendSettings( trackmate.getSettings() )
    writer.writeToFile();
    print( "TrackMate model saved to: " + saveFile.toString() + '\n' );
    model.getLogger().log("TrackMate model saved to: " + saveFile.toString() + '\n')
    
    saveFile = File(dirName, fileName + "_Tracks.xml")
    ExportTracksToXML.export(trackmate.getModel(), trackmate.getSettings(), saveFile)
    print( "Tracks saved to: " + saveFile.toString() + '\n' );
    model.getLogger().log("Tracks saved to: " + saveFile.toString() + '\n')
    
    # Display results
    displayTracks(model, imp)
    
    # capture overlay - RGB file
    image = trackmate.getSettings().imp
    capture = CaptureOverlayAction.capture(image, -1, imp.getNFrames(), logger)
    capture.setTitle("TracksOverlay")
    capture.show()
    
    return bestTrackerVal

if __name__ in ['__builtin__','__main__']:
    # Run the batch_open_images() function using the Scripting Parameters.
    path_to_images = batch_open_images(import_dir,
                               split_string(file_types),
                               split_string(filters),
                               process_sub
                              )  
    global bestTrackerVal
    bestTrackerVal = [0, 0.0, 0.0, 0.0]
    try:
        [int(x) for x in maxFrameGapString.split(',')]
    except ValueError, e1:
        IJ.log("One of the values in Max Frame Gap is not an integer!")
        print("Error: ", e1)
        sys.exit()
    global maxFrameGap
    global altLinkCost
    global linkMax
    global gapCloseMax
    maxFrameGap = [int(x) for x in maxFrameGapString.split(',')]
    altLinkCost = [float(x) for x in altLinkCostString.split(',')]
    linkMax = [float(x) for x in linkMaxString.split(',')]
    gapCloseMax = [float(x) for x in gapCloseMaxString.split(',')]
    IJ.log("Determining best tracking variable inputs on " + path_to_images[0])
    IJ.log("Max Frame Gap: " + maxFrameGapString.replace(",", ", "))
    IJ.log("Alt Link Cost: " + altLinkCostString.replace(",", ", "))
    IJ.log("Link Max: " + linkMaxString.replace(",", ", "))
    IJ.log("Gap Close Max: " + gapCloseMaxString.replace(",", ", "))
    # Process each image
    for imgPath in path_to_images:
        print imgPath
        global fullPath
        fullPath = imgPath
        global imp
        imp = IJ.openImage(imgPath)
        if imp:
            imp.show()
            #IJ.run("Set Scale...", "distance=1.1249 known=1 unit=um global")
            cal = imp.getCalibration()
            if timeUnit != "":
                cal.setTimeUnit(timeUnit)
            if frameInterval != 0:
                cal.frameInterval = frameInterval
            imp.setCalibration(cal)
            global fileName
            fileName = os.path.splitext(os.path.basename(imgPath))[0]
            global dirName 
            dirName = os.path.dirname(imgPath)
            if imp.getDimensions()[4] > 1:
                global trackingQualityTable
                rt_exist = WindowManager.getWindow("Tracking_Quality")
                if rt_exist == None or not isinstance(rt_exist, TextWindow):
                    trackingQualityTable = ResultsTable()
                    trackingQualityTable.show("Tracking_Quality")
                else:
                    trackingQualityTable = rt_exist.getTextPanel().getOrCreateResultsTable() 
                trackingQualityTable.setNaNEmptyCells(True)
                bestTrackerVal = cpProcess(imp, str(cpExe), str(cpModel), bestTrackerVal)
                imp.close()
                if trackingQualityTable.getHeadings()[0] == 'C1':
                    trackingQualityTable.renameColumn('C1', 'Max_Frame_Gap')
                    trackingQualityTable.renameColumn('C2', 'Alt_Link_Cost')
                    trackingQualityTable.renameColumn('C3', 'Link_Max')
                    trackingQualityTable.renameColumn('C4', 'Gap_Close_Max')
                    trackingQualityTable.renameColumn('C5', 'Num_Tracks_Kept')
                    trackingQualityTable.renameColumn('C6', 'Avg_Track_Quality')
                    trackingQualityTable.renameColumn('C7', 'Avg_Track_Length')
                    if trackingQualityTable.size()>1:
                        trackingQualityTable.renameColumn('C8', 'Norm_Num_Tracks_Kept')
                        trackingQualityTable.renameColumn('C9', 'Norm_Avg_Track_Quality')
                        trackingQualityTable.renameColumn('C10', 'Norm_Avg_Track_Length')
                        trackingQualityTable.renameColumn('C11', 'Tracking_Quality')
                    else:
                        trackingQualityTable.renameColumn('C8', 'Tracking_Quality')
                    trackingQualityTable.show("Tracking_Quality")
                if imgPath == path_to_images[0]:
                    trackingQualityTable.save(os.path.join(dirName,fileName + "_trackingQuality.csv"))
                    IJ.log("Tracking Quality calculation history saved to " + os.path.join(dirName,fileName + "_trackingQuality.csv"))
            else:
                IJ.log("Error processing " + dirName + "/" + fileName)
                IJ.log("Single time slice. Check if Z and T dimensions need to be swapped?")
    IJ.log("Batch processing for " + str(import_dir) + " complete!")
    if process_sub:
        IJ.log("All subfolders processed.")
    if filters != "":
        IJ.log("Filtered for images containing the word " + str(filters) + ".")
    if timeUnit != "":
        IJ.log("All images processed with time units: " + str(timeUnit) + " and frame interval: " + str(frameInterval))
    IJ.log("Tested the following ranges for MAX FRAME GAP: " + maxFrameGapString)
    IJ.log("Tested the following ranges for ALTERNATIVE LINK COST FACTOR: " + altLinkCostString)
    IJ.log("Tested the following ranges for LINKING MAX DISTANCE: " + linkMaxString)
    IJ.log("Tested the following ranges for GAP CLOSING MAX DISTANCE: " + gapCloseMaxString)
    IJ.log("*Optimized input variables*"
    IJ.log("Max frame gap: " + str(int(bestTrackerVal[0])))
    IJ.log("Alternative linking cost factor: " + str(bestTrackerVal[1]))
    IJ.log("Linking max distance: " + str(bestTrackerVal[2]))
    IJ.log("Gap closing max distance: " + str(bestTrackerVal[3]))
    IJ.log("Tracking quality was calculated with the following equation:")
    IJ.log(str(tracksKeptWeight) + " * # of tracks kept")
    IJ.log("+ " + str(avgTrackQualityWeight) + " * average track quality")
    IJ.log("+ " + str(avgTrackDurationWeight) + " * average track duration = tracking quality")

