#Quick script to loop over all channels to make fitting plots

#from ROOT import *

import subprocess
import sys

mzCut = sys.argv[1]
mzStr = mzCut.split(".")[0]

mwCut = sys.argv[2]
mwStr = mwCut.split(".")[0]

#Make the skim directory
subprocess.call("mkdir mvaDirs/skims/mz"+mzStr+"mw"+mwStr,shell=True)

#Make the mvaInput directory
subprocess.call("mkdir mvaDirs/inputs/mz"+mzStr+"mw"+mwStr,shell=True)

print "python scripts/makeMVAInputDilepton.py [\\\"ee\\\",\\\"mumu\\\"] mvaDirs/skims/mz"+mzStr+"mw"+mwStr+"/ mvaDirs/inputs/mz"+mzStr+"mw"+mwStr+"/"
subprocess.call("python scripts/makeMVAInputDilepton.py [\\\"ee\\\",\\\"mumu\\\"] mvaDirs/skims/mz"+mzStr+"mw"+mwStr+"/ mvaDirs/inputs/mz"+mzStr+"mw"+mwStr+"/",shell=True)
