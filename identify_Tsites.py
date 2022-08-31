#Takes 3 variables as input, in the following order:
#1. name of the TMO reactivity file
#2. name of the DMS reactivity file
#3. name of the output file
#If the two profiles don't have the same number of nucleotides, the script will fail.

#Typically you don't need to edit the profile files output from shapemapper as long as you aligned to the same sequence for TMO and DMS.

# example input: python identify_Tsites.py *TMO profile name* *DMS profile name* *Name of output file*

import sys
import argparse as ap

assert len(sys.argv) >= 2, "Not enough arguments"

fileName = sys.argv[1] #"RCW_TMO_RCW_profile.txt"
fileName2 = sys.argv[2]

try:
    newFileName = sys.argv[3]
except (IndexError):
    print("No new file name detected. Outputting into file Identified_Tsites.txt")
    newFileName = "Identified_Tsites.txt"

file1 = open(fileName, "r")
file2 = open(fileName2, "r")
lineCounter = 0

newFile = open(newFileName, "w")

newFile.write("Nuc    TMO-DMS    TMO/DMS    Read Depth(TMO)    Read Depth (DMS)    Read Depth (Control)    TMO-Rate    DMS-Rate    Control-Rate")
newFile.write("\n")
numSites = 0
for line1, line2 in zip(file1, file2):
    lineCounter +=1
    if (lineCounter > 1):
        curLine1 = line1.strip().split()
        curLine2 = line2.strip().split()
        if (curLine1[1] == "A" and curLine2[1] == "A"): # check to ensure that the nucleotide being looked at is an adenine
            difference = 0
            dividend = 0
            if (int(curLine1[3]) >= 200 and int(curLine1[10]) >= 200 and int(curLine2[3]) >= 200 and int(curLine2[10]) >= 200): #check to ensure correct read depth
                tReact = curLine1[23]
                dReact = curLine2[23]
                if (tReact != 'nan' and dReact != 'nan' and float(dReact) >= -0.02): #check to ensure that the DMS reactivity isn't some wacky value far in the negative
                    if(float(dReact) <= 0):
                        dReact = 0.01 #change to 0.01 to avoid having negative dividend values or dividing by 0
                    difference  = float(tReact) - float(dReact) 
                    dividend = float(tReact)/float(dReact)
                
                if (difference >= 0.20 and dividend >= 3): #check if TMO-DMS >= 0.2 and TMO/DMS  >= 3
                    numSites +=1
                    newFile.write(curLine1[0] + "    " + str(difference) + "    " + str(dividend) + "    " + curLine1[3] + "    " + curLine2[3] + "    " + curLine1[10] + "    " + curLine1[5] + "    " + curLine2[5] + "    " + curLine1[12])
                    newFile.write("\n")
                

newFile.write("\n")    
newFile.write("Overall Summary")  
newFile.write("\n")        
newFile.write("Number of T-sites found: " + str(numSites))        
newFile.close()
            
        