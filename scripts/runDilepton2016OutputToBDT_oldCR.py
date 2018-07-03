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
#subprocess.call("mkdir /scratch/data/TopPhysics/mvaDirs/inputs/2016/"+dirExt+"mz"+mzStr+"mw"+mwStr+"_zPlusJets_oldCR_test",shell=True)
#subprocess.call("rm /scratch/data/TopPhysics/mvaDirs/inputs/2016/"+dirExt+"mz"+mzStr+"mw"+mwStr+"_zPlusJets_oldCR_test/*",shell=True)

print "python scripts/makeMVAInputDilepton_oldZplusCR.py [\\\"ee\\\",\\\"mumu\\\"] /scratch/data/TopPhysics/mvaDirs/skims/2016/mz20mw20_zPlusJets_oldCR/ /scratch/data/TopPhysics/mvaDirs/inputs/2016/all/mz20mw20_zPlusJets_oldCR_test/ --2016 "+useSideBandRegion
subprocess.call("python scripts/makeMVAInputDilepton_oldZplusCR.py [\\\"ee\\\",\\\"mumu\\\"] /scratch/data/TopPhysics/mvaDirs/skims/2016/mz20mw20_zPlusJets_oldCR/ /scratch/data/TopPhysics/mvaDirs/inputs/2016/all/mz20mw20_zPlusJets_oldCR_test/ --2016 "+useSideBandRegion,shell=True)
