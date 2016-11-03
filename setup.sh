# Point LHAPDF to pdf sets
export LHAPATH=/cvmfs/cms.cern.ch/lhapdf/pdfsets/6.1.6/

export TQZ_TOOLS_PATH='.'

# Enable newrt gcc and python
source /scratch/shared/sw/gcc/x86_64-slc6-gcc62-opt/setup.sh
export PATH=/scratch/shared/sw/Python/2.7.4/x86_64-slc6-gcc48-opt/bin:$PATH
export PYTHONDIR=/scratch/shared/sw/Python/2.7.4/x86_64-slc6-gcc48-opt
export PYTHONPATH=/scratch/shared/sw/Python/2.7.4/x86_64-slc6-gcc48-opt
export PYTHONPATH=$PYTHONPATH:$ROOTSYS/lib
export LD_LIBRARY_PATH=$ROOTSYS/lib:$PYTHONDIR/lib:$LD_LIBRARY_PATH
