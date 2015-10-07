# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific aliases and functions


HISTSIZE=1000     # retain the last 1000 commands


source '/home/eepgadm/root/bin/thisroot.sh'
export ROOTSYS='/home/eepgadm/root/'

export SCRAM_ARCH=slc6_amd64_gcc491
export VO_CMS_SW_DIR='/cms/cmssw'

#source '/afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.sh' #Crab2
#source '/cms/crab3/3.3.16.rc2/slc6_amd64_gcc481/cms/crabclient/3.3.16.rc2/etc/profile.d/init.sh' #Crab3

alias InitCmsEnv='. /cms/cmssw/cmsset_default.sh'

alias  home='cd ~/;pwd'
alias  up='cd ..;pwd'
alias  ll='ls -l'
alias  lla='ls -la'
alias  l='ls -lh'

