#Quick script to loop over all channels to make fitting plots

#from ROOT import *

import subprocess
import sys

channels = {"ee":"1","mumu":"2"}

channelList = {"ee":"1","mumu":"2"}

metCut = sys.argv[1]
metStr = metCut.split(".")[0]

mtwCut = sys.argv[2]
mtwStr = mtwCut.split(".")[0]

#Make the skim directory
subprocess.call("mkdir mvaDirs/skims/met"+metStr+"mtw"+mtwStr,shell=True)

for chanName in channelList.keys():
    print "bin/analysisMain.exe -c configs/"+chanName+"MCdataConf.txt -u -t -k "+str(channelList[chanName])+" --jetRegion 2,1,3,4 -v 16383 -z --mvaDir mvaDirs/skims/met"+metStr+"mtw"+mtwStr + "/ --metCut " + str(metCut) + " --mtwCut " + str(mtwCut)
    subprocess.call("bin/analysisMain.exe -c configs/"+chanName+"MCdataConf.txt -u -t -k "+str(channelList[chanName])+" --jetRegion 2,1,3,4 -v 16383 -z --mvaDir mvaDirs/skims/met"+metStr+"mtw"+mtwStr + "/ --metCut " + str(metCut) + " --mtwCut " + str(mtwCut),shell=True)

#Make the mvaInput directory
subprocess.call("mkdir mvaDirs/inputs/met"+metStr+"mtw"+mtwStr,shell=True)

print "python scripts/makeMVAInput.py [\\\"ee\\\",\\\"mumu\\\"] mvaDirs/skims/met"+metStr+"mtw"+mtwStr+"/ mvaDirs/inputs/met"+metStr+"mtw"+mtwStr+"/"
subprocess.call("python scripts/makeMVAInput.py [\\\"ee\\\",\\\"mumu\\\"] mvaDirs/skims/met"+metStr+"mtw"+mtwStr+"/ mvaDirs/inputs/met"+metStr+"mtw"+mtwStr+"/",shell=True)

