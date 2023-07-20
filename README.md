# TreeAnalysis
Package for making and analyzing TTrees for STAR anisotropic flow analyses.

This is a **mostly** complete and standalone package to create TTrees from STAR picoDsts and also perform an anisotropic flow analysis on those trees. It will need to be run on RCF so it can access STAR libraries.

This README file is optimized for viewing in GitHub and can be found here: https://github.com/cracz/TreeAnalysis/

## Dependencies - Important!

To complete this package, use CVS to checkout the following StRoot libraries. Place them within the `StRoot/` directory of this package.

`StarClassLibrary`

`StBichsel`

`StEpdUtil`

`StEvent`

`StPicoDstMaker`

`StPicoEvent`


This package also uses a class called `ConfigReader` located in the `StRoot/` directory. This class is actually its own repository so it has been attached to this repository as a **submodule**. If this repo has been checked out from GitHub, you may need to issue the following command to also checkout the submodule:

`git submodule update --init --recursive`

The `ConfigReader` is used to read text files that contain the values of cuts that are important and/or often changed and supply them to the analysis programs. The necessary config file for reproducing results at 3.0 GeV is included in this package and called `config_3p0GeV.txt`.

## Compiling

1) Run the command `starver SL20d`.
2) In the `TreeAnalysis/` directory, run the command `cons` to perform necessary compilation of everything in the StRoot directory.
3) Go to the directory `StRoot/StEpdUtil/`, use a text editor to open the file `StEpdFastSim/StEpdFastSim.cxx`, and switch to the commented out line that includes StPicoEpdHit.h. You want to change it to look like this:

```c++
#include "../../StPicoEvent/StPicoEpdHit.h"  <---- on your laptop, need explicit path.
//#include "StRoot/StPicoEvent/StPicoEpdHit.h"
```
4) Run these commands to create the shared libraries `libStEpdUtil.so` and `libStPicoDst.so` and get back to the `TreeAnalysis/` directory:

```c++
make
cd ../StPicoEvent/
make
cd ../../
```

5) Run the following

```bash
cp StRoot/StEpdUtil/libStEpdUtil.so source/libs/
cp StRoot/StPicoEvent/libStPicoDst.so source/libs/
cp .sl73_gcc485/obj/StRoot/StBichsel/StBichsel.so source/libs/libStBichsel.so
mkdir CorrectionFiles/
```

4) Finally, go to the directory `source/` and run the command `make` to execute the makefile and compile `TreeAnalyzer.cxx` into an executable.

If nothing went wrong, then everything is all set!


## Making Trees

Before making or analyzing trees, navigate to the directory `submit/` and make sure the .xml files `treeMake.xml` and `treeAnalyze.xml` will be sending their output files into your own directories where you want them to go! In `treeMake.xml`, change every instance of `/star/data01/pwg/cracz/Data_3p0GeV_FXT/` to the location that you want the trees made from picoDsts to go, and change `/star/data01/pwg/cracz/Data_3p0GeV_FXT/out/` to the location that you want the corresponding .err and .out files to go. In `treeAnalyze.xml`, change `/star/data01/pwg/cracz/flowResults/` to the location that you want results from analyzing the trees to go. And in both xml files, change every instance of `/star/u/cracz/TreeAnalysis/` to the `TreeAnalysis/` directory that you checked out.

To start, go to the `submit/` directory and run the command 

```bash
star-submit treeMake.xml
```

The trees are made with the ROOT macro `MakeTrees.C`, which uses the `StRoot/TreeMaker/` class. If you need to modify what is saved within the produced trees, go to `StRoot/TreeMaker/`. Once the trees are created, you need to make a file list of all of those trees and supply that list to the `TreeAnalyzer` with `treeAnalyze.xml`. You can make use of the `getFileList.sh` script to quickly make a file list by going to the directory with all of the trees and running the command

```bash
getFileList.sh [jobIDofTheTrees]
```

If you need to rerun any jobs that failed to make a tree, use the script `findHolesInTrees.sh` in a similar way by supplying the job ID as the first argument and the number of total jobs submitted as the second argument. This will produce the command necessary to resubmit the failed jobs when you return to the directory where you submitted the original jobs.

## Analyzing Trees

An important part about analyzing the trees you make is that you can only analyze one file per job rather than a whole list like you could with picoDsts.

To start, go to the `submit/` directory and run the command 

```bash
star-submit treeAnalyze.xml
```

The `TreeAnalyzer` is an executable, so you only invoke it by name and supply the following arguments:

1) Input file name that has the tree.
2) Job ID given by the scheduler.
3) Name of the config file.
4) Name of the file with correction info (for re-centering and Fourier shifting of event planes) (produced by TreeAnalyzer if not present).
5) Name of the file with event plane resolutions (produced by you after 3 iterations of TreeAnalyzer are done, see below).

If you need to rerun any jobs that failed, use the script `Scripts/findHoles.sh` by supplying the job ID as the first argument and the number of total jobs submitted as the second argument. This will produce the command necessary to resubmit the failed jobs at the directory where you originally submitted them.

The various histograms and other plots will all be within files called `*.picoDst.result.root`, and the information necessary for event plane recentering and Fourier shifting will be in files called `correctionInfo_OUTPUT_*.root`. `TreeAnalyzer` takes 3 iterations to produce flat event plane distributions, and a 4th iteration to get correct flow results. After each of the first 3 iterations, go to the directory you have set in `treeAnalyze.xml` that holds the results and run

```bash
hadd correctionInfo_OUTPUT_*.root correctionInfo_INPUT.root
mv correctionInfo_INPUT.root [...]/TreeAnalysis/CorrectionFiles/3p#GeV/
rm *.root
```

where `[...]` is the path to your checked out copy of `TreeAnalysis/`, and `#` depends on what energy you're analyzing. This will set you up for the next iteration and remove the unneccessary output files from the iteration you just finished.

The program will determine what iteration you're on based on the status of `correctionInfo_INPUT.root`, so you don't have to do anything between iterations aside from combining and moving the correction info.

* Correction file doesn't exist = iteration 0 (Get re-centering info)
* Correction file exists, but TProfiles are empty = iteration 1 (Apply re-centering, get shifting info)
* Correction file exists and TProfiles are not empty = iteration 2 (Apply re-centering and shifting, get all correlations for event plane resolution calculation)
* Correction file exists, TProfiles not empty, and resolution file exists = flow can now be calculated

After iteration 3, if you want to recreate the results of the 3.0 GeV paper, use the file `resolutionInfo_INPUT_3p0GeV_averagedRes.root` for your resolutions. The details of where this is from is in the analysis note. Move or copy it to the `CorrectionFiles/` directory, change it's name to `resolutionInfo_INPUT.root`, and run the 4th iteration. Otherwise if you want to recreate a different set of resolution values, the info for event plane resolution calculations will be present in the following TProfiles of correlations between the various subevents, located in the main output ROOT files.

`p_TpcAB`

`p_TpcAEpdA`

`p_TpcAEpdB`

`p_TpcBEpdA`

`p_TpcBEpdB`

`p_EpdAEpdB`

Go to the directory you have set in `treeAnalyze.xml` that holds the results and run

```bash
hadd Results.picoDst.result.combined.root *.picoDst.result.root
mv Results.picoDst.result.combined.root [...]TreeAnalysis/Plotting/
cd [...]TreeAnalysis/Plotting/
root -l -b -q resolutions.cxx\(\"Results\",\"3\"\)
mv resolutionInfo_INPUT.root ../Correctionfiles/
```

This will combine your results from all of the jobs, move them to the `Plotting/` directory and use them to create the event plane resolutions for the inner EPD subevent, and then move those resolutions to the proper directory so that you're ready to initiate the 4th and final iteration. 

Once you run the fourth iteration, go back to the results directory and `hadd` the results together again with this command and you're ready to start plotting the results.

```bash
hadd Results.picoDst.result.combined.root *.picoDst.result.root
```

## Recreating figures from the paper

Move your main results file (`Results.picoDst.result.combined.root`) and `eventPlaneSystematics.root` to the `Plotting/` directory and run the following commands to recreate the plots shown in the 3.0 GeV paper (without systematic uncertainties).

```bash
root -l -b -q resolutionPlot.cxx
root -l -b -q prelimCentralityPlots.cxx\(\"Results\"\)
root -l -b -q prelimRapidityPlots.cxx\(\"Results\"\)
root -l -b -q prelimPtPlots.cxx\(\"Results\"\)
root -l -b -q acceptanceCombined.cxx\(\"Results\"\)
```



