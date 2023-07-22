import subprocess
import numpy as np
import os.path

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

choice = input()

directory = '/star/u/cracz/work/flow/TreeAnalysis/submit/systematics/'
subDirectory = 'fxt_3p0GeV/'
fileName = 'dict_systematics.npy'

if choice != 0 and choice != 1 and choice != 2 and choice != 3 and choice != 4 and choice != 5:
    print("Please input only an integer 0 - 5.")
    exit(0)
elif choice == 0:
    print("Cancelling.")
    exit(0)
elif choice == 1:
    subDirectory="fxt_3p0GeV/"
elif choice == 2:
    subDirectory="fxt_3p2GeV/"
elif choice == 3:
    subDirectory="fxt_3p5GeV/"
elif choice == 4:
    subDirectory="fxt_3p9GeV/"
elif choice == 5:
    subDirectory="fxt_4p5GeV/"

systemDict = np.load(directory+subDirectory+fileName, allow_pickle='TRUE').item()

for key in systemDict.keys():
    if os.path.isfile("correctionInfo_OUTPUT_"+systemDict[key]+"_0.root"):
        cmd = "hadd correctionInfo_INPUT_"+key+".root correctionInfo_OUTPUT_"+systemDict[key]+"*"
        sp = subprocess.Popen(cmd, shell=True)
        sp.communicate()

print("Done!")
