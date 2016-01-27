#!/bin/bash          
echo Updating filelists for datasets used ...
echo First deleting old filelists ...

rm $TQZ_TOOLS_PATH/configs/datasets/fileLists/DYJets10To50Files.txt
rm $TQZ_TOOLS_PATH/configs/datasets/fileLists/DYJets50Files.txt
rm $TQZ_TOOLS_PATH/configs/datasets/fileLists/eeRun2015CFiles.txt
rm $TQZ_TOOLS_PATH/configs/datasets/fileLists/eeRun2015DFiles.txt
rm $TQZ_TOOLS_PATH/configs/datasets/fileLists/emuRun2015CFiles.txt
rm $TQZ_TOOLS_PATH/configs/datasets/fileLists/emuRun2015DFiles.txt
rm $TQZ_TOOLS_PATH/configs/datasets/fileLists/mumuRun2015CFiles.txt
rm $TQZ_TOOLS_PATH/configs/datasets/fileLists/mumuRun2015DFiles.txt
rm $TQZ_TOOLS_PATH/configs/datasets/fileLists/sChannelFiles.txt
rm $TQZ_TOOLS_PATH/configs/datasets/fileLists/synchFiles.txt
rm $TQZ_TOOLS_PATH/configs/datasets/fileLists/tbarChannelFiles.txt
rm $TQZ_TOOLS_PATH/configs/datasets/fileLists/tbarWInclusiveFiles.txt
rm $TQZ_TOOLS_PATH/configs/datasets/fileLists/tChannelFiles.txt
rm $TQZ_TOOLS_PATH/configs/datasets/fileLists/ttbarDileptonFiles.txt
rm $TQZ_TOOLS_PATH/configs/datasets/fileLists/ttWFiles.txt
rm $TQZ_TOOLS_PATH/configs/datasets/fileLists/ttZFiles.txt
rm $TQZ_TOOLS_PATH/configs/datasets/fileLists/tWInclusiveFiles.txt
rm $TQZ_TOOLS_PATH/configs/datasets/fileLists/tZqFilesRun1.txt
rm $TQZ_TOOLS_PATH/configs/datasets/fileLists/tZqFiles.txt
rm $TQZ_TOOLS_PATH/configs/datasets/fileLists/wPlusJetsFiles.txt
rm $TQZ_TOOLS_PATH/configs/datasets/fileLists/WW2l2nuFiles.txt
rm $TQZ_TOOLS_PATH/configs/datasets/fileLists/WWinclusive.txt
rm $TQZ_TOOLS_PATH/configs/datasets/fileLists/WZinclusive.txt
rm $TQZ_TOOLS_PATH/configs/datasets/fileLists/ZZ2l2nuFiles.txt
rm $TQZ_TOOLS_PATH/configs/datasets/fileLists/ZZ2l2qFiles.txt
rm $TQZ_TOOLS_PATH/configs/datasets/fileLists/ZZ4lFiles.txt

echo Done!
echo Now outputting the lists of the dataset directories into their relevant files ...

ls /scratch/data/tZqSkimsRun2/DYJetsToLL_M-10to50/* -1d >> $TQZ_TOOLS_PATH/configs/datasets/fileLists/DYJets10To50Files.txt
ls /scratch/data/tZqSkimsRun2/DYJetsToLL_M-50/* -1d >> $TQZ_TOOLS_PATH/configs/datasets/fileLists/DYJets50Files.txt
ls /scratch/data/tZqSkimsRun2/eeRun2015C/* -1d >> $TQZ_TOOLS_PATH/configs/datasets/fileLists/eeRun2015CFiles.txt
ls /scratch/data/tZqSkimsRun2/eeRun2015D/* -1d >> $TQZ_TOOLS_PATH/configs/datasets/fileLists/eeRun2015DFiles.txt
ls /scratch/data/tZqSkimsRun2/emuRun2015C/* -1d >> $TQZ_TOOLS_PATH/configs/datasets/fileLists/emuRun2015CFiles.txt
ls /scratch/data/tZqSkimsRun2/emuRun2015D/* -1d >> $TQZ_TOOLS_PATH/configs/datasets/fileLists/emuRun2015DFiles.txt
ls /scratch/data/tZqSkimsRun2/mumuRun2015C/* -1d >> $TQZ_TOOLS_PATH/configs/datasets/fileLists/mumuRun2015CFiles.txt
ls /scratch/data/tZqSkimsRun2/mumuRun2015D/* -1d >> $TQZ_TOOLS_PATH/configs/datasets/fileLists/mumuRun2015DFiles.txt
ls /scratch/data/tZqSkimsRun2/sChannel_4f/* -1d >> $TQZ_TOOLS_PATH/configs/datasets/fileLists/sChannelFiles.txt
ls /scratch/data/tZqSkimsRun2/tZq_synchTest/* -1d  >> $TQZ_TOOLS_PATH/configs/datasets/fileLists/synchFiles.txt
ls /scratch/data/tZqSkimsRun2/tbarChannel_4f/* -1d >> $TQZ_TOOLS_PATH/configs/datasets/fileLists/tbarChannelFiles.txt
ls /scratch/data/tZqSkimsRun2/tW_antitop_5f/* -1d  >> $TQZ_TOOLS_PATH/configs/datasets/fileLists/tbarWInclusiveFiles.txt
ls /scratch/data/tZqSkimsRun2/tChannel_4f/* -1d >> $TQZ_TOOLS_PATH/configs/datasets/fileLists/tChannelFiles.txt
ls /scratch/data/tZqSkimsRun2/ttbarDilepton/* -1d >> $TQZ_TOOLS_PATH/configs/datasets/fileLists/ttbarDileptonFiles.txt
ls /scratch/data/tZqSkimsRun2/ttW/* -1d >> $TQZ_TOOLS_PATH/configs/datasets/fileLists/ttWFiles.txt
ls /scratch/data/tZqSkimsRun2/ttZ/* -1d >> $TQZ_TOOLS_PATH/configs/datasets/fileLists/ttZFiles.txt
ls /scratch/data/tZqSkimsRun2/tW_top_5f/* -1d >> $TQZ_TOOLS_PATH/configs/datasets/fileLists/tWInclusiveFiles.txt
ls /nfs/data/tZqSkims/tZq4Flavour3Lepton/* -1d >> $TQZ_TOOLS_PATH/configs/datasets/fileLists/tZqFilesRun1.txt
ls /scratch/data/tZqSkimsRun2/tZq_ll_4Flavour3Lepton/* -1d >> $TQZ_TOOLS_PATH/configs/datasets/fileLists/tZqFiles.txt
ls /scratch/data/tZqSkimsRun2/wPlusJets/* -1d >> $TQZ_TOOLS_PATH/configs/datasets/fileLists/wPlusJetsFiles.txt
ls /scratch/data/tZqSkimsRun2/WW2l2nu/* -1d >> $TQZ_TOOLS_PATH/configs/datasets/fileLists/WW2l2nuFiles.txt
ls /scratch/data/tZqSkimsRun2/WWinclusive/* -1d >> $TQZ_TOOLS_PATH/configs/datasets/fileLists/WWinclusive.txt
ls /scratch/data/tZqSkimsRun2/WZJets/* -1d >> $TQZ_TOOLS_PATH/configs/datasets/fileLists/WZinclusive.txt
ls /scratch/data/tZqSkimsRun2/ZZ2l2nu/* -1d >> $TQZ_TOOLS_PATH/configs/datasets/fileLists/ZZ2l2nuFiles.txt
ls /scratch/data/tZqSkimsRun2/ZZ2l2q/* -1d >> $TQZ_TOOLS_PATH/configs/datasets/fileLists/ZZ2l2qFiles.txt
ls /scratch/data/tZqSkimsRun2/ZZ4l/* -1d >> $TQZ_TOOLS_PATH/configs/datasets/fileLists/ZZ4lFiles.txt

echo Done!
echo Filelists have been updated.
