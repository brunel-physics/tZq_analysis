# Point LHAPDF to pdf sets
export LHAPATH=/cvmfs/cms.cern.ch/lhapdf/pdfsets/6.1.6/

export TQZ_TOOLS_PATH='.'

# Enable required afs software
source /afs/cern.ch/sw/lcg/contrib/gcc/6.2.0abi4/x86_64-slc6-gcc62-opt/setup.sh
export PATH=/afs/cern.ch/sw/lcg/external/Python/2.7.4/x86_64-slc6-gcc48-opt/bin:$PATH
export PYTHONDIR=/afs/cern.ch/sw/lcg/external/Python/2.7.4/x86_64-slc6-gcc48-opt
export PYTHONPATH=/afs/cern.ch/sw/lcg/external/Python/2.7.4/x86_64-slc6-gcc48-opt
export PYTHONPATH=$PYTHONPATH:$ROOTSYS/lib
export LD_LIBRARY_PATH=$ROOTSYS/lib:$PYTHONDIR/lib:$LD_LIBRARY_PATH:/opt/rh/python27/root/usr/lib64
