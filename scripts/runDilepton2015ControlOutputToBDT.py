#Quick script to loop over all channels to make fitting plots

#from ROOT import *

import subprocess
import sys

channels = {"emu":"4"}

channelList = {"emu":"4"}

#Make the mvaInput directory
subprocess.call("mkdir mvaDirs/inputs/2015/ttbarControl",shell=True)

print "python scripts/makeMVAInputDilepton.py [\\\"ee\\\",\\\"mumu\\\"] mvaDirs/skims/2015/ttbarControl/ mvaDirs/inputs/2015/ttbarControl --2015 --ttbar"
subprocess.call("python scripts/makeMVAInputDilepton.py [\\\"ee\\\",\\\"mumu\\\"] mvaDirs/skims/2015/ttbarControl mvaDirs/inputs/2015/ttbarControl/ --2015 --ttbar",shell=True)

