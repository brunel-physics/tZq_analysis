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
subprocess.call("mkdir mvaDirs/inputs/2016/"+dirExt+"mz"+mzStr+"mw"+mwStr,shell=True)

print "python scripts/makeMVAInputDilepton.py [\\\"ee\\\",\\\"mumu\\\"] mvaDirs/skims/2016/mz"+mzStr+"mw"+mwStr+"/ mvaDirs/inputs/2016/"+dirExt+"mz"+mzStr+"mw"+mwStr+"/ --2016 "+useSideBandRegion
subprocess.call("python scripts/makeMVAInputDilepton.py [\\\"ee\\\",\\\"mumu\\\"] mvaDirs/skims/2016/mz"+mzStr+"mw"+mwStr+"/ mvaDirs/inputs/2016/"+dirExt+"mz"+mzStr+"mw"+mwStr+"/ --2016 "+useSideBandRegion,shell=True)
