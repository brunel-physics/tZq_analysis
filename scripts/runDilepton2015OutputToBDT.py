#Quick script to loop over all channels to make fitting plots

#from ROOT import *

import subprocess
import sys

mzCut = sys.argv[1]
mzStr = mzCut.split(".")[0]

mwCut = sys.argv[2]
mwStr = mwCut.split(".")[0]

useSideBandRegion = ""
dirExt = "all/"
if len(sys.argv) > 3 and sys.argv[3] == "-s" :
    useSideBandRegion = "-s"
    dirExt = "sigCtrl/"

#Make the mvaInput directory
subprocess.call("mkdir /scratch/data/TopPhysics/mvaDirs/inputs/2015/"+dirExt+"mz"+mzStr+"mw"+mwStr,shell=True)

print "python scripts/makeMVAInputDilepton.py [\\\"ee\\\",\\\"mumu\\\"] /scratch/data/TopPhysics/mvaDirs/skims/2015/mz"+mzStr+"mw"+mwStr+"/ /scratch/data/TopPhysics/mvaDirs/inputs/2015/"+dirExt+"mz"+mzStr+"mw"+mwStr+"/ --2015 "+useSideBandRegion
subprocess.call("python scripts/makeMVAInputDilepton.py [\\\"ee\\\",\\\"mumu\\\"] /scratch/data/TopPhysics/mvaDirs/skims/2015/mz"+mzStr+"mw"+mwStr+"/ /scratch/data/TopPhysics/mvaDirs/inputs/2015/"+dirExt+"mz"+mzStr+"mw"+mwStr+"/ --2015 "+useSideBandRegion,shell=True)
