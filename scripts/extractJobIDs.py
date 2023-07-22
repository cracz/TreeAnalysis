# This script pieces together a proper "saveSystematics.py" script
# based on your last submission of systematics jobs by reading all 
# of the "submitOut_*.txt" files by extracting the job IDs and the
# corresponding variation code (i.e. nSigPr_high, m2Ka_low, etc.).

import subprocess
import numpy as np
import glob

foundFiles = glob.glob("submitOut_*.txt")

initialLines = ['import numpy as np', '', 'systemDict = {', ''] #line1, line2, line3
finalLines = ['}', ' ', 'np.save(\'../dict_systematics.npy\', systemDict)', ' ', 'print("Dictionary of current job IDs for systematics saved to ../dict_systematics.npy.")']

if len(foundFiles) == 0:
    print("No submitOut files found!")
else:
    lastFile = foundFiles[-1]

    newFile = open("saveSystematics.py", 'w')
    newFile.write('\n'.join(initialLines))

    for ithFileName in foundFiles:
        periodPosition = ithFileName.find(".")
        variationID = ithFileName[10: periodPosition]
        ithFile = open(ithFileName, 'r')
        jobIDFound = False;
        
        for line in ithFile:
            if 'Could not find per-packed sandbox.' in line:
                jobIDFound = True;
                underscorePosition = line.find("_")
                jobID = line[underscorePosition-32:underscorePosition]

                if ithFileName != lastFile:
                    newFile.write('\t\"'+variationID+'\": \"'+jobID+'\",\n') #comma at the end
                    #print(variationID + ": " + jobID)
                else:
                    newFile.write('\t\"'+variationID+'\": \"'+jobID+'\"\n') #no comma at the end
                break
        
        if jobIDFound == False:
            print("ERROR - No job ID found for "+variationID)

        ithFile.close()

    newFile.write('\n'.join(finalLines))
    newFile.close()
        

print("saveSystematics.py has been updated with the most recent job IDs.")
print("Done!")
