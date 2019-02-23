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
dirExt = "all/"
if len(sys.argv) > 3 and sys.argv[3] == "-s" :
    useSideBandRegion = "-s"
    dirExt = "sigCtrl/"

#Make the skim directory
subprocess.call("mkdir /scratch/data/TopPhysics/mvaDirs/skims/2015/mz"+mzStr+"mw"+mwStr,shell=True)

for chanName in channelList.keys():
    print "bin/analysisMain.exe -c configs/2015/"+chanName+"MCdataConf.txt -u -t -k "+str(channelList[chanName])+" -v 16383 -z --dilepton --mvaDir mvaDirs/skims/2015/mz"+mzStr+"mw"+mwStr + "/ --mzCut " + str(mzCut) + " --mwCut " + str(mwCut)
    subprocess.call("bin/analysisMain.exe -c configs/2015/"+chanName+"MCdataConf.txt -u -t -k "+str(channelList[chanName])+" -v 16383 -z --dilepton --mvaDir /scratch/data/TopPhysics/mvaDirs/skims/2015/mz"+mzStr+"mw"+mwStr + "/ --mzCut " + str(mzCut) + " --mwCut " + str(mwCut),shell=True)

#Make the mvaInput directory
subprocess.call("mkdir /scratch/data/TopPhysics/mvaDirs/inputs/2015/"+dirExt+"mz"+mzStr+"mw"+mwStr,shell=True)
subprocess.call("rm /scratch/data/TopPhysics/mvaDirs/inputs/2015/"+dirExt+"mz"+mzStr+"mw"+mwStr+"/*",shell=True)

print "python scripts/makeMVAInputDilepton.py [\\\"ee\\\",\\\"mumu\\\"] /scratch/data/TopPhysics/mvaDirs/skims/2015/mz"+mzStr+"mw"+mwStr+"/ /scratch/data/TopPhysics/mvaDirs/inputs/2015/"+dirExt+"mz"+mzStr+"mw"+mwStr+"/ --2015 "+useSideBandRegion
subprocess.call("python scripts/makeMVAInputDilepton.py [\\\"ee\\\",\\\"mumu\\\"] /scratch/data/TopPhysics/mvaDirs/skims/2015/mz"+mzStr+"mw"+mwStr+"/ /scratch/data/TopPhysics/mvaDirs/inputs/2015/"+dirExt+"mz"+mzStr+"mw"+mwStr+"/ --2015 "+useSideBandRegion ,shell=True)


