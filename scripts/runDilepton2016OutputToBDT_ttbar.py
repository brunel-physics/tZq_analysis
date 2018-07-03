#Quick script to loop over all channels to make fitting plots

#from ROOT import *

import subprocess
import sys

useSideBandRegion = ""
dirExt = "all/"
if len(sys.argv) > 3 and sys.argv[3] == "-s" :
    useSideBandRegion = "-s"
    dirExt = "sigCtrl/"

#Make the mvaInput directory
#subprocess.call("mkdir /scratch/data/TopPhysics/mvaDirs/inputs/2016/emu",shell=True)
#subprocess.call("rm /scratch/data/TopPhysics/mvaDirs/inputs/2016/emu/*",shell=True)

print "python scripts/makeMVAInputDilepton.py [\\\"emu\\\"] /scratch/data/TopPhysics/mvaDirs/skims/2016/emu/ /scratch/data/TopPhysics/mvaDirs/inputs/2016/emu_test/ --2016 "+useSideBandRegion
subprocess.call("python scripts/makeMVAInputDilepton.py [\\\"emu\\\"] /scratch/data/TopPhysics/mvaDirs/skims/2016/emu/ /scratch/data/TopPhysics/mvaDirs/inputs/2016/emu_test/ --2016 "+useSideBandRegion,shell=True)
