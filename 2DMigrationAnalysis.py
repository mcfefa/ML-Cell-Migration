#@ File(label='Choose a directory of videos', style='directory') import_dir
#@ String(label='File types', value='tif') file_types
#@ File (label = "Cellpose .exe", style = "file") cpExe
#@ File (label = "Cellpose model", style = "file") cpModel
#@ String(label='Filter', value='') filters
#@ String (visibility=MESSAGE, value="Leave blank to analyze all files without filtering for keywords", required=false) msg
#@ Boolean(label='Process Subfolders', value=False) process_sub

from fiji.plugin.trackmate.visualization.hyperstack import HyperStackDisplayer
from fiji.plugin.trackmate.io import TmXmlReader, TmXmlWriter
from fiji.plugin.trackmate import TrackMate, Logger, Dimension, Settings, SelectionModel, Spot, TrackModel, SpotCollection, Model
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettingsIO, DisplaySettings
from fiji.plugin.trackmate.features.track import TrackAnalyzer, TrackSpotQualityFeatureAnalyzer, TrackDurationAnalyzer, AbstractTrackAnalyzer
from java.io import File
from fiji.util.gui import GenericDialogPlus
from ij import IJ, WindowManager
from ij.gui import Overlay, Line, Plot
from fiji.plugin.trackmate.visualization.table import TrackTableView, AllSpotsTableView
from jarray import array
from java.awt import Color
from fiji.plugin.trackmate.tracking import SpotTrackerFactory
from fiji.plugin.trackmate.tracking.jaqaman import LAPUtils, SparseLAPTrackerFactory
from fiji.plugin.trackmate.providers import SpotAnalyzerProvider, EdgeAnalyzerProvider, TrackAnalyzerProvider
from fiji.plugin.trackmate.cellpose.CellposeSettings import PretrainedModel
from fiji.plugin.trackmate.cellpose import CellposeDetectorFactory
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettings
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettings
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettingsIO
from fiji.plugin.trackmate.action import CaptureOverlayAction
from fiji.plugin.trackmate.cellpose.CellposeSettings import PretrainedModel
from fiji.plugin.trackmate.util import TMUtils, LogRecorder
import math, os, sys, array, time
from java.util import HashSet
from fiji.plugin.trackmate.action import TrimNotVisibleAction

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
    # Create the list that will be returned by this function.
    images = []
    for img_path in path_to_images:
        # IJ.openImage() returns an ImagePlus object or None.
        imp = IJ.openImage(img_path)
        # An object equals True and None equals False.
        if imp:
            images.append(imp)
    return images
 
def split_string(input_string):
    '''Split a string to a list and strip it
    :param input_string: A string that contains semicolons as separators.
    '''
    string_splitted = input_string.split(';')
    # Remove whitespace at the beginning and end of each string
    strings_striped = [string.strip() for string in string_splitted]
    return strings_striped
 
if __name__ in ['__builtin__','__main__']:
    # Run the batch_open_images() function using the Scripting Parameters.
    images = batch_open_images(import_dir,
                               split_string(file_types),
                               split_string(filters),
                               process_sub
                              )
    for image in images:
        # Call the toString() method of each ImagePlus object.
        print(image)

# Hard-coded, need to change later
wkdir = "E:/ML_Migration/tm_test"
file = File("/Users/ktrang/Documents/tm_test/A2_1.xml")
imp = IJ.openImage("E:/ML_Migration/tm_test/A2_1_short.tif")
imp.show()
IJ.run("Set Scale...", "distance=2016 known=1791 unit=um global")

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
		# SS = SS + (getResult("Mean", i)-meanIntensity)*(getResult("Mean", i)-meanIntensity)
	sd = math.sqrt(abs(SS)/(numTracks-1))
	for id in model.getTrackModel().trackIDs(True):
		featureVal = fm.getTrackFeature(id, featureName)
		if(featureVal < meanFeature-2*sd):	
			trackIDsUnder.append(id)
	return trackIDsUnder
	

def findTrackingQuality(model, settings, maxFrameGap, altLinkCost, linkMax, gapCloseMax):
	# Configure trackerFactory
	fm = model.getFeatureModel()
	trackingQuality = 0
	settings.trackerFactory = SparseLAPTrackerFactory()
	settings.trackerSettings = settings.trackerFactory.getDefaultSettings()
	settings.trackerSettings["MAX_FRAME_GAP"] = maxFrameGap
	settings.trackerSettings["ALTERNATIVE_LINKING_COST_FACTOR"] = altLinkCost
	settings.trackerSettings["LINKING_MAX_DISTANCE"] = linkMax
	settings.trackerSettings["GAP_CLOSING_MAX_DISTANCE"] = gapCloseMax
	settings.trackerSettings["ALLOW_GAP_CLOSING"] = True
	settings.addTrackAnalyzer(TrackSpotQualityFeatureAnalyzer())
	settings.addTrackAnalyzer(TrackDurationAnalyzer())
   	
   	trackmate = TrackMate(model, settings)
   	
   	spotAnalyzerProvider = SpotAnalyzerProvider(1)
	for key in spotAnalyzerProvider.getKeys():
		print( key )
		settings.addSpotAnalyzerFactory( spotAnalyzerProvider.getFactory( key ) )
	edgeAnalyzerProvider = EdgeAnalyzerProvider()
	for  key in edgeAnalyzerProvider.getKeys():
		print( key )
		settings.addEdgeAnalyzer( edgeAnalyzerProvider.getFactory( key ) )	
	
	trackmate.getModel().getLogger().log( settings.toStringFeatureAnalyzersInfo() )
	if model.getSpots().getNSpots( True ) == 0:
	   	# Detection
		print('Spot detection')
		trackmate.execDetection()
		model.getLogger().log("FINISHED spot detection")
		# Initial filtering, not necessary
		#print('Spot initial filtering')
		#ok = ok and trackmate.execInitialSpotFiltering()
		#model.getLogger().log("FINISHED spot initial filtering")
		# Compute spot features.
		print('Computing spot features')
		trackmate.computeSpotFeatures(True) 
		model.getLogger().log("FINISHED spot features")
		print('Filtering spots')
		trackmate.execSpotFiltering(True) 
		model.getLogger().log("FINISHED spot filtering")
	
	print('Tracking')
	ok = trackmate.execTracking()
	if( ok ):
		print('TrackMate completed successfully.' )
	else:
		print(str(trackmate.getErrorMessage()))
	print( 'Found %d spots in %d tracks.' % ( model.getSpots().getNSpots( True ) , model.getTrackModel().nTracks( True ) ) )
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
	avgTrackLength = 0
	for id in model.getTrackModel().trackIDs(True):
		tracksKept = tracksKept + 1
		avgTrackQuality= avgTrackQuality + fm.getTrackFeature(id, 'TRACK_MEAN_QUALITY')
		avgTrackLength= avgTrackLength + fm.getTrackFeature(id, 'TRACK_DISPLACEMENT')
	avgTrackQuality = avgTrackQuality/tracksKept
	avgTrackLength = avgTrackLength/tracksKept
	model.getLogger().log('# of tracks kept: ' + str(tracksKept))
	model.getLogger().log('Average track quality: ' + str(avgTrackQuality))
	model.getLogger().log('Average track length: ' + str(avgTrackLength) + ' frames')
	trackingQuality = 0.2*tracksKept + 0.4*avgTrackQuality + 0.4*avgTrackLength
	model.getLogger().log('Tracking quality: ' + str(trackingQuality))
	model.beginUpdate()
	for id in model.getTrackModel().trackIDs(False):
		track = model.getTrackModel().trackSpots(id).iterator()
		while track.hasNext():
			spot = track.next()
			for edge in model.getTrackModel().edgesOf(spot):
				model.removeEdge(edge)
	model.endUpdate()
	for id in model.getTrackModel().trackIDs(False):
		model.setTrackVisibility(id, True)
	print(trackingQuality)
	return trackmate, trackingQuality

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
	# TRY THIS: https://forum.image.sc/t/trackmate-calculate-a-new-custom-spot-feature-and-update-the-spots-table-to-export-it-as-csv/87014/3

def testNewVar(model, var):
	middle = findMiddle(model)
	for id in model.getTrackModel().trackIDs(True):
		track = model.getTrackModel().trackSpots(id)
		spot = track.iterator().next()
		print("Frame: %d, X: %.5f, %s: %d" % (spot.getFeature("FRAME"), spot.getFeature("POSITION_X"), var, spot.getFeature(var)))
		
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
		if dotProd >= 0:
			dirMigVec = [dim*dotProd for dim in dirVec]
			dirMig = math.sqrt(dirMigVec[0]**2+dirMigVec[1]**2)
			fm.putTrackFeature(id, "DirectedMigration", dirMig)
		#print(fm.getTrackFeature(id, "DirectedMigration"))
		# Calculate directed migration ratio
		fm.putTrackFeature(id, "DirectedMigrationRatio", dirMig/disp)
		#print(fm.getTrackFeature(id, "DirectedMigrationRatio"))	
				
def outputTrackData(model, sm, ds, wkdir):
	# Save spot, edge, and track statistics
	trackTableView = TrackTableView(model, sm, ds, "/Users/ktrang/Documents/tm_test/A2_1_short.tif")
	trackTableView.getSpotTable().getPanel()
	trackTableView.getSpotTable().exportToCsv(File(os.path.join(wkdir,"spotData.csv")))
	trackTableView.getEdgeTable().exportToCsv(File(os.path.join(wkdir,"edgeData.csv")))
	trackTableView.getTrackTable().exportToCsv(File(os.path.join(wkdir,"trackData.csv")))
	
def outputTrackDataNew(model, sm, ds, wkdir):
	# Save spot, edge, and track statistics
	trackTableView = TrackTableView(model, sm, ds, "/Users/ktrang/Documents/tm_test/A2_1_short.tif")
	trackTableView.getSpotTable().getPanel()
	trackTableView.getSpotTable().exportToCsv(File(os.path.join(wkdir,"spotData_NEW.csv")))
	trackTableView.getEdgeTable().exportToCsv(File(os.path.join(wkdir,"edgeData_NEW.csv")))
	trackTableView.getTrackTable().exportToCsv(File(os.path.join(wkdir,"trackData_NEW.csv")))
	
# adapted from https://forum.image.sc/t/importing-cellpose-via-scyjava-for-trackmate-automation-in-python-3/84176/6
def cpProcess(imp, cpExe, cpModel):
	model = Model()
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
   	    'SIMPLIFY_CONTOURS' : False,
   	}   	
	
	#maxFrameGap = [2, 4, 6, 8, 10]
	#altLinkCost = [0.2, 0.4, 0.6, 0.8, 1.0]
	#linkMax = [5.0, 10.0, 15.0, 20.0]
	#gapCloseMax = [20.0, 40.0, 60.0, 80.0, 100.0, 120.0]
	# Test values
	maxFrameGap = [2]
	altLinkCost = [0.2, 0.4]
	linkMax = [5.0]
	gapCloseMax = [20.0]
	
	# Determine tracker settings and instantiate trackmate
	bestTrackingQuality = 0
	bestTrackerVal = [0, 0.0, 0.0, 0.0]
	# Test all different combinations of input variables and measure tracking quality
	for maxFrameGapVal in maxFrameGap:
		for altLinkCostVal in altLinkCost:
			for linkMaxVal in linkMax:
				for gapCloseMaxVal in gapCloseMax:
					[trackmate, trackingQuality] = findTrackingQuality(model, settings, maxFrameGapVal, altLinkCostVal, linkMaxVal, gapCloseMaxVal)
					#findTrackingQuality(model, settings, maxFrameGapVal, altLinkCostVal, linkMaxVal, gapCloseMaxVal)
					#trackingQuality = 0
					if trackingQuality > bestTrackingQuality:
						bestTrackingQuality = trackingQuality
						bestTrackerVal = [maxFrameGapVal, altLinkCostVal, linkMaxVal, gapCloseMaxVal]
					print [maxFrameGapVal, altLinkCostVal, linkMaxVal, gapCloseMaxVal]
	
	settings.addAllAnalyzers()
	trackmate.computeSpotFeatures(True)
	trackmate.computeEdgeFeatures(True)
	trackmate.computeTrackFeatures(True)
	
	trim = TrimNotVisibleAction()
	trim.execute(trackmate, SelectionModel(model), DisplaySettingsIO.readUserDefault(), WindowManager.getCurrentImage())
	
	# Display results
	model = trackmate.getModel()
	sm = SelectionModel( model )
	ds = DisplaySettings()
	ds = DisplaySettingsIO.readUserDefault()
	ds.spotDisplayedAsRoi = True
	displayer =  HyperStackDisplayer( model, sm, imp, ds )
	displayer.render()
	displayer.refresh()
	
	#https://github.com/TASBE/TASBEImageAnalytics/blob/5d6fc1a64b4c17e263451fa4252c94dc86193d14/scripts/cellStatsTracking.py
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
#testNewVar(model, "LR")
# See what features we calculate for/can output
#print fm.getTrackFeatureNames()
#print fm.getSpotFeatureNames()
	calcDirMig(model)
#display_results_in_GUI(trackmate)
#outputTrackData(model, sm, ds, wkdir)
	# Save results
	outputTrackDataNew(model, sm, DisplaySettingsIO.readUserDefault(), wkdir)
	saveFile = TMUtils.proposeTrackMateSaveFile( settings, logger )

	writer = TmXmlWriter( saveFile, logger )
	writer.appendLog( logger.toString() )
	writer.appendModel( trackmate.getModel() )
	writer.appendSettings( trackmate.getSettings() )
	writer.writeToFile();
	print( "Results saved to: " + saveFile.toString() + '\n' );
	
	# Display results
	model = trackmate.getModel()
	sm = SelectionModel( model )
	#ds = DisplaySettings()
	ds = DisplaySettingsIO.readUserDefault()
	ds.spotDisplayedAsRoi = True
	displayer =  HyperStackDisplayer( model, sm, imp, ds )
	displayer.render()
	displayer.refresh()
	
	# capture overlay - RGB file
	image = trackmate.getSettings().imp
	capture = CaptureOverlayAction.capture(image, -1, imp.getNFrames(), logger)
	capture.setTitle("TracksOverlay")
	capture.show()

cpProcess(imp, str(cpExe), str(cpModel))