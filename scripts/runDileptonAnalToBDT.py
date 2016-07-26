#Quick script to loop over all channels to make fitting plots

#from ROOT import *

import subprocess
import sys

channels = {"ee":"1","mumu":"2"}

channelList = {"ee":"1","mumu":"2"}

mzCut = sys.argv[1]
mzStr = mzCut.split(".")[0]

mwCut = sys.argv[2]
mwStr = mwCut.split(".")[0]

#Make the skim directory
subprocess.call("mkdir mvaDirs/skims/mz"+mzStr+"mw"+mwStr,shell=True)

for chanName in channelList.keys():
    print "bin/analysisMain.exe -c configs/"+chanName+"MCdataConf.txt -u -t -k "+str(channelList[chanName])+" -v 16383 -z --dilepton --mvaDir mvaDirs/skims/mz"+mzStr+"mw"+mwStr + "/ --mzCut " + str(mzCut) + " --mwCut " + str(mwCut)
    subprocess.call("bin/analysisMain.exe -c configs/"+chanName+"MCdataConf.txt -u -t -k "+str(channelList[chanName])+" -v 16383 -z --dilepton --mvaDir mvaDirs/skims/mz"+mzStr+"mw"+mwStr + "/ --mzCut " + str(mzCut) + " --mwCut " + str(mwCut),shell=True)

#Make the mvaInput directory
subprocess.call("mkdir mvaDirs/inputs/mz"+mzStr+"mw"+mwStr,shell=True)

print "python scripts/makeMVAInputDilepton.py [\\\"ee\\\",\\\"mumu\\\"] mvaDirs/skims/mz"+mzStr+"mw"+mwStr+"/ mvaDirs/inputs/mz"+mzStr+"mw"+mwStr+"/"
subprocess.call("python scripts/makeMVAInputDilepton.py [\\\"ee\\\",\\\"mumu\\\"] mvaDirs/skims/mz"+mzStr+"mw"+mwStr+"/ mvaDirs/inputs/mz"+mzStr+"mw"+mwStr+"/",shell=True)

