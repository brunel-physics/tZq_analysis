#Quick script to loop over all channels to make fitting plots

#from ROOT import *

import subprocess
import sys

channels = {"1":"eee","2":"eemu","4":"emumu","8":"mumumu","16":"eeeInv","32":"eemuInv","64":"emumuInv","128":"mumumuInv"}

channelList = {"eee":"17","eemu":"34","emumu":"68","mumumu":"136"}

metCut = sys.argv[1]
metStr = metCut.split(".")[0]

mtwCut = sys.argv[2]
mtwStr = mtwCut.split(".")[0]

#Make the skim directory
subprocess.call("mkdir mvaDirs/skims/trilepton/met"+metStr+"mtw"+mtwStr,shell=True)

#Make the mvaInput directory
subprocess.call("mkdir mvaDirs/inputs/trilepton/met"+metStr+"mtw"+mtwStr,shell=True)

print "python scripts/makeMVAInput.py [\\\"eee\\\",\\\"eemu\\\",\\\"emumu\\\",\\\"mumumu\\\"] mvaDirs/skims/trilepton/met"+metStr+"mtw"+mtwStr+"/ mvaDirs/inputs/trilepton/met"+metStr+"mtw"+mtwStr+"/"
subprocess.call("python scripts/makeMVAInput.py [\\\"eee\\\",\\\"eemu\\\",\\\"emumu\\\",\\\"mumumu\\\"] mvaDirs/skims/trilepton/met"+metStr+"mtw"+mtwStr+"/ mvaDirs/inputs/trilepton/met"+metStr+"mtw"+mtwStr+"/",shell=True)

