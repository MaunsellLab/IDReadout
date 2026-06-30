IDReadout To Do:

Kernel average is not getting staleDirections when regressions are remain  Is a problem for plotting across kernel directions

plotSideTypeKernelAverage failed to find the summary files
Trouble finding in function probeDirs = findProbeDirsWithAverageData(sideType)
Search involves:
dataFolder = fullfile(domainFolder(mfilename('fullpath')), 'Data');
probeFolders = dir(fullfile(dataFolder, 'Probe*'));
  probeDir = str2double(token{1});
  dataFile = fullfile(dataFolder, probeName, 'AverageKernels', sideType, 'AverageKernelPlotData.mat');
  AverageKernels have not been produced, let alone sideType averageKernels
  dailyUpdate seems to be out of order

Kernel plot titles are screwed up


fitAcrossOffsetBetaMeasurements('NBoot', nBoot, 'Animal', animal); ???
saving to:
Batch summary saved: /Users/Shared/Data/IDReadout/IDR/Data/scalarNoiseRegression_batchSummary.mat
Does this need to be identified by animal?

Revise:

Getting different number of trials for beta and scale fits

  makeKernels
    strip parsing out of compileKernelSummary
    eliminate the stub: recomputeSessionKernelStruct
  
Eliminate: 
  updateAcrossOffsetSummaries -- remove selectCompStatsEntry

IDReadout Analysis Overview:
============================

IDReadout (IDR) and IDReadout2 (IDT) are the Knot plugins that created the
.dat files for this project. IDR was a single probe direction per file. IDT
allowed for randomly interleaved probe directions

Matlab Project:

All the analysis files are maintained in a Matlab Project

Folder Hierarchy:
=================

All analyses files access data files that lie within subfolder of a folder
whose path is returned by calling domainFolder(mfilename('fullpath')). This makes it 
simple to create parallel test folder hierarchies by changing the path returned by 
folderPath.m

By convention, all the Matlab files specific to IDReadout are in
$PATH/Code, and all the data are in $PATH/Data.  .dat files are in
$PATH/Data/DatFiles.  Derived .mat files are stored in assigned locations. 
Plot output from the analyses are placed in $PATH/Plots.

Matlab Project:
=================

All the analysis files are maintained in a Matlab Project

Standard Analysis:
==================

convertIDRData.m 
----------------
Converts the .dat files into a standard .mat file.  By convention, that
file contains a header structure and a trials structure.  convertIDRData
also does a bit of clean up correcting some mis-assignment of events that
occurred in early IDR files.  The converted .mat files are placed int 
$PATH/Data/FullSessions, with the same name as the base .dat file (e.g.,
baseName.mat). An additional information files is also spun off in that 
sub-folder, with the name baseName_info.mat, but it is not used in 
analysis.   

Conversion from .dat to .mat is slow, so we try to do that once for each
.dat file. 

makeKernels.m
-------------
makeKernels is a preliminary stage of analysis that breaks session files into
sub-files that are associated with individual probe directions. This parts are
stored in sub-folders named $PATH/Data/probe45, in which the final numberical 
characters specify the probe direction (±45° in this example). makeKernels 
creates three types of .mat files that are stored in different sub-folders. 
Each of these .mat files has a name that includes "_probeXX" for quick 
confirmation of their contents. 

The Kernels subfolder contains extracted kernels. There is a kernel for the
preferred direction and a kernel for the probed direction.  Ancillary data are
included.

The NoiseMatrices subfolder contains files with individual noise traces on which the
kernels were constructed. These are needed for analyses that do bootstraps across
multiple kernel to establish confidence intervals. 

In its default mode, makeKernels looks for missing derived files and creates only those.
A "replace" flag forces makeKernels to recreate each of the derived files.

When makeKernels produces output files, it also creates pdf plots showing kernels for
each probe direction within a session. These are placed in $PATH/Plots/probeXX/Kernels.







