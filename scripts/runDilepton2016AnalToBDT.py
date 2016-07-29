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

useSideBandRegion = ""
if len(sys.argv) > 3 and sys.argv[3] == "-s" : 
    useSideBandRegion = "-s"

#Make the skim directory
subprocess.call("mkdir mvaDirs/skims/2016/mz"+mzStr+"mw"+mwStr,shell=True)

for chanName in channelList.keys():
    print "bin/analysisMain.exe -c configs/2016"+chanName+"MCdataConf.txt -u -t -k "+str(channelList[chanName])+" -v 16383 -z --dilepton --mvaDir mvaDirs/skims/2016/mz"+mzStr+"mw"+mwStr + "/ --mzCut " + str(mzCut) + " --mwCut " + str(mwCut) + " --2016"
    subprocess.call("bin/analysisMain.exe -c configs/2016/"+chanName+"MCdataConf.txt -u -t -k "+str(channelList[chanName])+" -v 16383 -z --dilepton --mvaDir mvaDirs/skims/2016/mz"+mzStr+"mw"+mwStr + "/ --mzCut " + str(mzCut) + " --mwCut " + str(mwCut) + " --2016",shell=True)

#Make the mvaInput directory
subprocess.call("mkdir mvaDirs/inputs/2016/mz"+mzStr+"mw"+mwStr,shell=True)

print "python scripts/make2016MVAInputDilepton.py [\\\"ee\\\",\\\"mumu\\\"] mvaDirs/skims/2016/mz"+mzStr+"mw"+mwStr+"/ mvaDirs/inputs/2016/mz"+mzStr+"mw"+mwStr+"/ --2016 "+useSideBandRegion
subprocess.call("python scripts/2016makeMVAInputDilepton.py [\\\"ee\\\",\\\"mumu\\\"] mvaDirs/skims/2016/mz"+mzStr+"mw"+mwStr+"/ mvaDirs/inputs/2016/mz"+mzStr+"mw"+mwStr+"/ --2016 "+useSideBandRegion ,shell=True)

