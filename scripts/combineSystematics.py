import subprocess
import numpy as np
import glob

systemDict = np.load('/star/u/cracz/work/flow/dict_systematics.npy', allow_pickle='TRUE').item()

#nFiles = str(raw_input("Enter the number of files expected per job: "))
if len(sys.argv) == 1:
    print("Input number of files expected per job as the first argument.")
    return

nFiles = sys.argv[1]

for key in systemDict.keys():

    foundFiles = glob.glob(systemDict[key]+"_*.root")
    if len(foundFiles) == 0:
        print("No files found with ID \""+systemDict[key]+"\". Check if the dictionary of job IDs was updated!")
        break

    #if key == "7p2":
        #nFiles = "674"

    cmd1 = "root -l -b -q ~/Scripts/combine.cxx\(\\\""+systemDict[key]+"\\\","+nFiles+"\)"
    sp = subprocess.Popen(cmd1, shell=True)
    sp.communicate()

    cmd2 = "mv "+systemDict[key]+".picoDst.result.combined.root "+key+".picoDst.result.combined.root"
    sp = subprocess.Popen(cmd2, shell=True)
    sp.communicate()

print("Done!")
