#!/bin/bash
echo Updating filelists for datasets used ...
echo First deleting old filelists ...

# Normal datasets

rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/eeRun2016Files.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/eeRun2016BFiles.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/emuRun2016Files.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/mumuRun2016Files.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/metRun2016Files.txt

rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets10To50Files.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets50Files.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets50MadgraphFiles.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets50MadgraphUnsmearedFiles.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/sChannelFiles.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tbarChannelFiles.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tbarWInclusiveFiles.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tChannelFiles.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbarDileptonFiles.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbarInclusivePowerhegFiles.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbarInclusivePowerhegExt3Files.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttWlnuFiles.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttW2qFiles.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttZ2l2nuFiles.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttZ2qFiles.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tWInclusiveFiles.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tWZ_ll_Files.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tZqFilesRun1.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tZqFiles.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tHqFiles.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/wPlusJetsFiles.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/WW1l1nu2qFiles.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/WW2l2nuFiles.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/WZ1l1nu2q.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/WZ2l2q.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/WZjets.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ZZ2l2nuFiles.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ZZ2l2qFiles.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ZZ4lFiles.txt

# Synchronisation files

rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tZqSynchFiles.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/WZsynchFiles.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttZsynchFiles.txt
rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbarSynchFiles.txt


echo Done!
echo Now outputting the lists of the dataset directories into their relevant files ...

# Normal datasets

ls /scratch/data/tZqSkimsRun2016/eeRun2016B/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/eeRun2016BFiles.txt
ls /scratch/data/tZqSkimsRun2016/eeRun2016B/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/eeRun2016Files.txt
ls /scratch/data/tZqSkimsRun2016/eeRun2016C/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/eeRun2016Files.txt
ls /scratch/data/tZqSkimsRun2016/eeRun2016D/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/eeRun2016Files.txt
ls /scratch/data/tZqSkimsRun2016/emuRun2016B/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/emuRun2016Files.txt
ls /scratch/data/tZqSkimsRun2016/emuRun2016C/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/emuRun2016Files.txt
ls /scratch/data/tZqSkimsRun2016/emuRun2016D/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/emuRun2016Files.txt
ls /scratch/data/tZqSkimsRun2016/mumuRun2016B/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/mumuRun2016Files.txt
ls /scratch/data/tZqSkimsRun2016/mumuRun2016C/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/mumuRun2016Files.txt
ls /scratch/data/tZqSkimsRun2016/mumuRun2016D/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/mumuRun2016Files.txt
ls /scratch/data/tZqSkimsRun2016/metRun2016B/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/metRun2016Files.txt
ls /scratch/data/tZqSkimsRun2016/metRun2016C/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/metRun2016Files.txt
ls /scratch/data/tZqSkimsRun2016/metRun2016D/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/metRun2016Files.txt

ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_M-10to50/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets10To50Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_M-50_amcatnlo/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets50Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_M-50_madgraph/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets50MadgraphFiles.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_M-50_madgraph_unsmeared/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets50MadgraphUnsmearedFiles.txt
ls /scratch/data/tZqSkimsRun2016/sChannel_4f/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/sChannelFiles.txt
ls /scratch/data/tZqSkimsRun2016/tbarChannel_4f/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tbarChannelFiles.txt
ls /scratch/data/tZqSkimsRun2016/tW_antitop_5f/* -1d  >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tbarWInclusiveFiles.txt
ls /scratch/data/tZqSkimsRun2016/tWZ_5f/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tWZ_ll_Files.txt
ls /scratch/data/tZqSkimsRun2016/tChannel_4f/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tChannelFiles.txt
ls /scratch/data/tZqSkimsRun2016/ttbarDilepton/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbarDileptonFiles.txt
ls /scratch/data/tZqSkimsRun2016/ttbarInclusive_powerheg_ext3/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbarInclusivePowerhegFiles.txt
ls /scratch/data/tZqSkimsRun2016/ttbarInclusive_powerheg_ext4/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbarInclusivePowerhegFiles.txt
ls /scratch/data/tZqSkimsRun2016/ttbarInclusive_powerheg_ext3/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbarInclusivePowerhegExt3Files.txt
ls /scratch/data/tZqSkimsRun2016/ttWlnu/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttWlnuFiles.txt
ls /scratch/data/tZqSkimsRun2016/ttW2q/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttW2qFiles.txt
ls /scratch/data/tZqSkimsRun2016/ttZ2l2nu/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttZ2l2nuFiles.txt
ls /scratch/data/tZqSkimsRun2016/ttZ2q/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttZ2qFiles.txt
ls /scratch/data/tZqSkimsRun2016/tW_top_5f/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tWInclusiveFiles.txt
ls /scratch/data/tZqSkimsRun2016/tZq_ll_4Flavour3Lepton/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tZqFiles.txt
ls /scratch/data/tZqSkimsRun2016/tHq/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tHqFiles.txt
ls /scratch/data/tZqSkimsRun2016/wPlusJets/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/wPlusJetsFiles.txt
ls /scratch/data/tZqSkimsRun2016/WW1l1nu2q/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/WW1l1nu2qFiles.txt
ls /scratch/data/tZqSkimsRun2016/WW2l2nu/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/WW2l2nuFiles.txt
ls /scratch/data/tZqSkimsRun2016/WZ1l1nu2q/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/WZ1l1nu2q.txt
ls /scratch/data/tZqSkimsRun2016/WZ2l2q/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/WZ2l2q.txt
ls /scratch/data/tZqSkimsRun2016/WZJets/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/WZjets.txt
ls /scratch/data/tZqSkimsRun2016/ZZ2l2nu/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ZZ2l2nuFiles.txt
ls /scratch/data/tZqSkimsRun2016/ZZ2l2q/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ZZ2l2qFiles.txt
ls /scratch/data/tZqSkimsRun2016/ZZ4l/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ZZ4lFiles.txt

# Synchronisation files

echo Done!
echo Filelists have been updated.
