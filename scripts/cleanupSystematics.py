import subprocess
import numpy as np

print(
'''
"What energy is this?"
"[0] Cancel"
"[1] 3.0 GeV"
"[2] 3.2 GeV"
"[3] 3.5 GeV"
"[4] 3.9 GeV"
"[5] 4.5 GeV"
''')

choice1 = input()

resultsDirectory = "$PWG/flowResults/"
submitDirectory  = "/star/u/cracz/work/flow/TreeAnalysis/submit/systematics/"
submitSubDirectory = "fxt_3p0GeV/"

if choice1 != 0 and choice1 != 1 and choice1 != 2 and choice1 != 3 and choice1 != 4 and choice1 != 5:
    print("Please input only an integer 0 - 5.")
    exit(0)
elif choice1 == 0:
    print("Cancelling.")
    exit(0)
elif choice1 == 1:
    resultsDirectory = "flowResults/"
    submitSubDirectory = "fxt_3p0GeV/"
elif choice1 == 2:
    resultsDirectory = "flowResults_3p2GeV/"
    submitSubDirectory = "fxt_3p2GeV/"
elif choice1 == 3:
    resultsDirectory = "flowResults_3p5GeV/"
    submitSubDirectory = "fxt_3p5GeV/"
elif choice1 == 4:
    resultsDirectory = "flowResults_3p9GeV/"
    submitSubDirectory = "fxt_3p9GeV/"
elif choice1 == 5:
    resultsDirectory = "flowResults_4p5GeV/"
    submitSubDirectory = "fxt_4p5GeV/"



print(
'''
"What size of variation?"
"[0] Cancel"
"[1] 20 percent"
"[2] 30 percent"
''')

choice2 = input()

variationDirectory = "20percentVariations/"

if choice2 != 0 and choice2 != 1 and choice2 != 2:
    print("Please input only an integer 0 - 2.")
    exit(0)
elif choice2 == 0:
    print("Cancelling.")
    exit(0)
elif choice2 == 1:
    variationDirectory = "20percentVariations/"
elif choice2 == 2:
    variationDirectory = "30percentVariations/"




systemDict = np.load(submitDirectory+submitSubDirectory+"dict_systematics.npy", allow_pickle='TRUE').item()

keys = systemDict.keys()

for key in keys:
    cmd1 = "rm -r "+submitDirectory+submitSubDirectory+variationDirectory+"sched"+systemDict[key]+"*"
    cmd2 = "rm -r "+submitDirectory+submitSubDirectory+variationDirectory+"submitOut_"+key+".txt"
    cmd3 = "rm $PWG/"+resultsDirectory+systemDict[key]+"_*.root"
    cmd4 = "rm $PWG/"+resultsDirectory+"correctionInfo_OUTPUT_"+systemDict[key]+"*.root"
    cmd5 = "rm $PWG/"+resultsDirectory+"out/"+systemDict[key]+"*.err"
    cmd6 = "rm $PWG/"+resultsDirectory+"out/"+systemDict[key]+"*.out"

    p1 = subprocess.Popen(cmd1, shell=True)
    p1.communicate()
    p2 = subprocess.Popen(cmd2, shell=True)
    p2.communicate()
    p3 = subprocess.Popen(cmd3, shell=True)
    p3.communicate()
    p4 = subprocess.Popen(cmd4, shell=True)
    p4.communicate()
    p5 = subprocess.Popen(cmd5, shell=True)
    p5.communicate()
    p6 = subprocess.Popen(cmd6, shell=True)
    p6.communicate()

    print(key+" done!")

print("Cleanup is finished!")
