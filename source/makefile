SRCS = ../StRoot/ConfigReader/ConfigReader.cxx TreeAnalyzer.cxx SetupAttributes.cxx
OBJS = $(SRCS:.cxx=.o)
DEPS = ../include/FlowUtils.h ../include/SetupAttributes.h ../include/NSigmaCorrectionUtils.h
TARGET = TreeAnalyzer
TARGET2 = nSigmaInfo

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
ROOTLIBS     += -lEG

#INCFLAGS = -I$(ROOTSYS)/include -I/afs/rhic.bnl.gov/star/packages/SL19b/StRoot -I/afs/rhic.bnl.gov/star/packages/SL19b/StRoot/StPicoEvent -I/afs/rhic.bnl.gov/star/packages/SL19b/StRoot/StEpdUtil
INCFLAGS = -I$(ROOTSYS)/include -I../include/ -I./ -I../ -I../StRoot -I../StRoot/StEvent -I../StRoot/StBichsel
LIBFLAGS = -L$(ROOTSYS)/lib -L./libs -Wl,-R./libs
SOLIBS = -lStPicoDst -lStEpdUtil -lStBichsel
#for each non-standard dynamic library location -L a corresponding -Wl,-R should be specified
#the -L's to StPicoEvent, StEpdUtil, and StBichsel are for locating libStPicoEvent.so, libStEpdUtil.so, and libStBichsel.so
#the -l's are to link these .so files
# The relative path between the TreeAnalysis executable and the .so libraries needs to be maintained in order 
# for things to work using the scheduler!

CC = g++ #-m32
FLAGS = -Wall -g -fPIC -std=c++11

COMPILE = $(CC) $(FLAGS)


all: $(TARGET) $(TARGET2)

$(TARGET): $(OBJS)
	$(COMPILE) -o $@ $^ $(SOLIBS) $(ROOTCFLAGS) $(ROOTLIBS) $(INCFLAGS) $(LIBFLAGS)

$(TARGET2): nSigmaInfo.o ../StRoot/ConfigReader/ConfigReader.o
	$(COMPILE) -o $@ $^ $(SOLIBS) $(ROOTCFLAGS) $(ROOTLIBS) $(INCFLAGS) $(LIBFLAGS)

TreeAnalyzer.o: $(DEPS) ../StRoot/ConfigReader/ConfigReader.h TreeAnalyzer.cxx
	$(COMPILE) -o $@ -c TreeAnalyzer.cxx $(SOLIBS) $(ROOTCFLAGS) $(ROOTLIBS) $(INCFLAGS) $(LIBFLAGS)

../StRoot/ConfigReader/ConfigReader.o: ../StRoot/ConfigReader/ConfigReader.h ../StRoot/ConfigReader/ConfigReader.cxx
	$(COMPILE) -o $@ -c ../StRoot/ConfigReader/ConfigReader.cxx $(ROOTCFLAGS) $(ROOTLIBS)

nSigmaInfo.o: nSigmaInfo.cxx ../StRoot/ConfigReader/ConfigReader.cxx
	$(COMPILE) -o $@ -c nSigmaInfo.cxx $(SOLIBS) $(ROOTCFLAGS) $(ROOTLIBS) $(INCFLAGS) $(LIBFLAGS)

SetupAttributes.o: SetupAttributes.cxx
	$(COMPILE) -o $@ -c SetupAttributes.cxx $(SOLIBS) $(ROOTCFLAGS) $(ROOTLIBS) $(INCFLAGS) $(LIBFLAGS)


.PHONY: clean

clean:
	rm -f $(TARGET) $(OBJS) nSigmaInfo.o nSigmaInfo



# target : dependencies
#	action

# $@ = current target
# $^ = current dependencies
# $< = name of the related file that caused the action
