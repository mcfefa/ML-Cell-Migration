# ML-Cell-Migration
pipeline for an automated machine learning-driven assay for robust and accessible analysis of cell migration

## Overview
This repository contains code for an [ImageJ](https://imagej.net/software/imagej/) Jython (.py) plugin for automated analysis of 2D rectangular migration assays for robust and accessible analysis of cell migration. The code is written in [Jython](https://imagej.net/scripting/jython/) and has been tested in [FIJI (FIJI Is Just ImageJ)](https://imagej.net/software/fiji/downloads) with ImageJ version 1.54j. The plugin is dependent on an installation of [Cellpose](https://github.com/MouseLand/cellpose) for segmentation before beginning the tracking process with [TrackMate](https://imagej.net/plugins/trackmate/) (tested on v7.2.1+).

Example data used for analysis is located in the **examples** folder. Shorter image sets, denoted as _small or _tiny are also provided.

## Requirements
- FIJI distribution of ImageJ
- [TrackMate](https://imagej.net/plugins/trackmate/)

## Contents
- [2DMigrationAnalysis](2DMigrationAnalysis.py):
