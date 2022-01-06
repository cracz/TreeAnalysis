# TreeAnalysis
Package for making and analyzing TTrees for STAR anisotropic flow analyses.

This is a **mostly** complete and standalone package to create TTrees from STAR picoDsts and also perform an anisotropic flow analysis on those trees. It will need to be run on RCF so it can access STAR libraries.

## Dependencies - Important!

To complete this package, use CVS to checkout the following StRoot libraries. Place them within the `StRoot/` directory of this package.
`StarClassLibrary`
`StBichsel`
`StEpdUtil`
`StEvent`
`StPicoDstMaker`
`StPicoEvent`

This package also uses a class called `ConfigReader` located in the `StRoot/` directory. This class is actually its own repository so it has been attached to this repository as a **submodule**. Most things about the `TreeAnalysis` repo should be the same as other git repos, but a few things will be a little different when dealing with `ConfigReader`. For example, you can make changes to `ConfigReader` from within this repo and push the changes back to the main `ConfigReader` repo, but there's just some extra steps. Check out Git's info on submodules for more info.

The `ConfigReader` is used to read text files that contain the values of cuts that are important and/or often changed and supply them to the analysis programs. Before running things, you'll need to create your own config file. You can find an example of one within the `StRoot/ConfigReader/` directory, so copy that one, modify the values how you need, and place it in the main `TreeAnalysis/` directory.

## Compiling

1) Set your version of STAR libraries with the command `starver SL20d` (or whatever library version you need). This package was developed with `SL20d` so use that if there's problems with others.
2) In the `TreeAnalysis/` directory, run the command `make` to execute the makefile to compile `TreeAnalyzer.cxx` into an executable.
3) Same directory, run the command `cons` to produce the necessary \*.so files of everything in the StRoot directory. These are required by the `TreeMaker`.

If nothing went wrong, then everything is all set!

## Making Trees

The trees are made with the ROOT macro `MakeTrees.C`, which uses the `StRoot/TreeMaker/` class. If you need to modify what is saved within the produced trees, go to `StRoot/TreeMaker/`. When submitting jobs to make trees, you need to send `MakeTrees.C`, your config file, and the whole `StRoot/` directory. An example command section of the .xml file is shown here:

```xml
<command>
  starver SL20d
  cons
  root4star -q -b -l MakeTrees.C\(0,$INPUTFILECOUNT,\"$FILELIST\",\"$SCRATCH\",\"$JOBID\",\"config_3p0GeV.txt\",0\)
</command>
```

## Analyzing Trees

An important part about analyzing the trees you make is that you can only analyze one file per job rather than a whole list like you could with picoDsts.

The `TreeAnalyzer` will be an executable, so you only invoke it by name and supply the following arguments:

1) Input file name that has the tree.
2) Job ID given by the scheduler.
3) Name of the config file.
4) Name of the file with correction info (for re-centering and Fourier shifting) (produced by TreeAnalyzer if not present).
5) Name of the file with event plane resolutions (produced by you, see below).

For your .xml file, send `TreeAnalyzer`, the correction and resolution files, the config file, and the `libs/` directory. An example command section would look like this: 

```xml
<command>
  starver SL20d
  TreeAnalyzer $INPUTFILE0 $JOBID config_3p0GeV.txt correctionInfo_INPUT.root resolutionInfo_INPUT.root
</command>
```

`TreeAnalyzer` uses the event plane method with re-centering and Fourier shifting, so it takes 3 iterations to produce flow results. Currently it also applies event plane resolutions particle-by-particle, so it actually takes 4 iterations. For re-centering and Fourier shifting, the program will save the necessary information to files called `correctionInfo_OUTPUT_*.root`. You'll need to `hadd` these into one file and supply them back to the program for the next iteration. I usually just call it `correctionInfo_INPUT.root`, but the new name could be whatever you want. 

The program will determine what iteration you're on based on the status of `correctionInfo_INPUT.root`, so you don't have to do anything between iterations aside from combining and supplying the correction info. You should always supply the names of the correction and resolution files to the program since it won't be a problem if they don't exist.

* Correction file doesn't exist = iteration 0 (Get re-centering info)
* Correction file exists, but TProfiles are empty = iteration 1 (Apply re-centering, get shifting info)
* Correction file exists and TProfiles are not empty = iteration 2 (Apply re-centering and shifting, get all correlations)
* Correction file exists, TProfiles not empty, and resolution file exists = flow can now be calculated

For event plane resolutions, this info will be present in the following TProfiles of correlations between the various subevents, located in the main output ROOT files. These will be filled after the 3rd iteration and then they can be used by the analyzer to calculate the resolutions.
`p_TpcAB`
`p_TpcAEpdA`
`p_TpcAEpdB`
`p_TpcBEpdA`
`p_TpcBEpdB`
`p_EpdAEpdB`

Currently, the analyzer requires that the resolution file contains a `TH1D` of resolution vs centrality ID called `h_resolutions`. The x-axis should NOT be in percentages! It should be following the normal ID that can be seen in the section of `TreeAnalyzer.cxx` where the `centrality` variable is assigned. So essentially `h_resolutions` will have an x-axis from 0 to 16, with 16 bins, and will look like a backwards version of the resolutions vs centrality percentages.
