#How to use this script:
#First make sure this script and ReactivityProfile.py are in the same directory as the profile you are trying to convert
#Only requires one argument: The file name of the shapemapper _profile.txt output from shapemapper
#python Convert_DMS_Single.py example_profile.txt

import sys
from ReactivityProfile import ReactivityProfile

target = sys.argv[1]

curProfile = ReactivityProfile(target, bg=0.02, ignorents="")
curProfile.normalize(DMS=True)
# write normed reactivities to file
curProfile.writeReactivity('{}.dms'.format(target[:-12])) 
