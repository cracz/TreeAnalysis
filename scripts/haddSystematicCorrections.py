import subprocess
import numpy as np
import os.path

systemDict = np.load('/star/u/cracz/work/flow/dict_systematics.npy', allow_pickle='TRUE').item()

for key in systemDict.keys():
    if os.path.isfile("correctionInfo_OUTPUT_"+systemDict[key]+"_0.root"):
        cmd = "hadd correctionInfo_INPUT_"+key+".root correctionInfo_OUTPUT_"+systemDict[key]+"*"
        sp = subprocess.Popen(cmd, shell=True)
        sp.communicate()

print("Done!")
