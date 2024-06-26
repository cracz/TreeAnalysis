import sys
import subprocess
import numpy as np
import glob

directory = '/star/u/cracz/work/flow/TreeAnalysis/submit/systematics/'
subDirectory = 'fxt_3p0GeV/'
variationDirectory = '20percentVariations/'
fileName = 'dict_systematics.npy'
nFilesInput = 0
nFiles = '0'

if len(sys.argv) == 1:
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

    nFilesInput = input("Enter the number of files expected per job: ")

    if not isinstance(nFilesInput, int):
        print("Input was not an integer.")
        exit(0)
    else:
        nFiles = str(nFilesInput)

    print('''
    "What size of variation?"
    "[0] Cancel"
    "[1] 20 percent"
    "[2] 30 percent"
    ''')

    choice2 = input()

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

else:
    print("Incorrect number of arguments. Use zero arguments.")
    exit(0)



systemDict = np.load(directory+subDirectory+fileName, allow_pickle='TRUE').item()

for key in systemDict.keys():

    foundFiles = glob.glob(systemDict[key]+"_*.root")
    if len(foundFiles) == 0:
        print("No files found with ID \""+systemDict[key]+"\". Check if the dictionary of job IDs was updated!")
        break

    cmd1 = "~/Scripts/findHolesSyst.sh "+systemDict[key]+" "+nFiles+" "+directory+"/"+subDirectory+"/"+variationDirectory
    sp = subprocess.Popen(cmd1, shell=True)
    sp.communicate()

print("Done!")
