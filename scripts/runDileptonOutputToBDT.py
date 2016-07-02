#Quick script to loop over all channels to make fitting plots

#from ROOT import *

import subprocess
import sys

metCut = sys.argv[1]
metStr = metCut.split(".")[0]

mtwCut = sys.argv[2]
mtwStr = mtwCut.split(".")[0]

#Make the skim directory
subprocess.call("mkdir mvaDirs/skims/met"+metStr+"mtw"+mtwStr,shell=True)

#Make the mvaInput directory
subprocess.call("mkdir mvaDirs/inputs/met"+metStr+"mtw"+mtwStr,shell=True)

print "python scripts/makeMVAInputDilepton.py [\\\"ee\\\",\\\"mumu\\\"] mvaDirs/skims/met"+metStr+"mtw"+mtwStr+"/ mvaDirs/inputs/met"+metStr+"mtw"+mtwStr+"/"
subprocess.call("python scripts/makeMVAInputDilepton.py [\\\"ee\\\",\\\"mumu\\\"] mvaDirs/skims/met"+metStr+"mtw"+mtwStr+"/ mvaDirs/inputs/met"+metStr+"mtw"+mtwStr+"/",shell=True)

