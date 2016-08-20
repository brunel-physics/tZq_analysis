#!/bin/bash          
echo Updating filelists for datasets used ...
echo First deleting old filelists ...

# Normal datasets

rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/eeRun2015CFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/eeRun2015DFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/emuRun2015CFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/emuRun2015DFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/mumuRun2015CFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/mumuRun2015DFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/metRun2015CFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/metRun2015DFiles.txt

rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/DYJets10To50Files.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/DYJets50Files.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/DYJets50v2Files.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/DYJets50MadgraphFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/sChannelFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/tbarChannelFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/tbarWInclusiveFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/tbarWInclusiveScaleUpFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/tbarWInclusiveScaleDownFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/tChannelFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/ttbarDileptonFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/ttbarInclusivePowerhegFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/ttbarInclusivePowerhegExt3Files.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/ttWlnuFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/ttW2qFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/ttZ2l2nuFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/ttZ2qFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/tWInclusiveFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/tWInclusiveScaleUpFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/tWInclusiveScaleDownFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/tZqFilesRun1.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/tZqFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/tHqFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/wPlusJetsFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/WW1l1nu2qFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/WW2l2nuFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/WZ1l1nu2q.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/WZ2l2q.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/WZjets.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/ZZ2l2nuFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/ZZ2l2qFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/ZZ4lFiles.txt

# Synchronisation files

rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/tZqSynchFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/WZsynchFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/ttZsynchFiles.txt
rm $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/ttbarSynchFiles.txt


echo Done!
echo Now outputting the lists of the dataset directories into their relevant files ...

# Normal datasets

ls /scratch/data/tZqSkimsRun2015/eeRun2015C/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/eeRun2015CFiles.txt
ls /scratch/data/tZqSkimsRun2015/eeRun2015D/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/eeRun2015DFiles.txt
ls /scratch/data/tZqSkimsRun2015/emuRun2015C/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/emuRun2015CFiles.txt
ls /scratch/data/tZqSkimsRun2015/emuRun2015D/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/emuRun2015DFiles.txt
ls /scratch/data/tZqSkimsRun2015/mumuRun2015C/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/mumuRun2015CFiles.txt
ls /scratch/data/tZqSkimsRun2015/mumuRun2015D/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/mumuRun2015DFiles.txt
ls /scratch/data/tZqSkimsRun2015/metRun2015C/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/metRun2015CFiles.txt
ls /scratch/data/tZqSkimsRun2015/metRun2015D/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/metRun2015DFiles.txt

ls /scratch/data/tZqSkimsRun2015/DYJetsToLL_M-10to50/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/DYJets10To50Files.txt
ls /scratch/data/tZqSkimsRun2015/DYJetsToLL_M-50/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/DYJets50Files.txt
ls /scratch/data/tZqSkimsRun2015/DYJetsToLL_M-50_amcatnlo_v2/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/DYJets50v2Files.txt
ls /scratch/data/tZqSkimsRun2015/DYJetsToLL_M-50_madgraph/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/DYJets50MadgraphFiles.txt
ls /scratch/data/tZqSkimsRun2015/sChannel_4f/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/sChannelFiles.txt
ls /scratch/data/tZqSkimsRun2015/tbarChannel_4f/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/tbarChannelFiles.txt
ls /scratch/data/tZqSkimsRun2015/tW_antitop_5f/* -1d  >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/tbarWInclusiveFiles.txt
ls /scratch/data/tZqSkimsRun2015/tChannel_4f/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/tChannelFiles.txt
ls /scratch/data/tZqSkimsRun2015/ttbarDilepton/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/ttbarDileptonFiles.txt
ls /scratch/data/tZqSkimsRun2015/ttbarInclusive_powerheg_ext3/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/ttbarInclusivePowerhegFiles.txt
ls /scratch/data/tZqSkimsRun2015/ttbarInclusive_powerheg_ext4/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/ttbarInclusivePowerhegFiles.txt
ls /scratch/data/tZqSkimsRun2015/ttbarInclusive_powerheg_ext3/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/ttbarInclusivePowerhegExt3Files.txt
ls /scratch/data/tZqSkimsRun2015/ttWlnu/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/ttWlnuFiles.txt
ls /scratch/data/tZqSkimsRun2015/ttW2q/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/ttW2qFiles.txt
ls /scratch/data/tZqSkimsRun2015/ttZ2l2nu/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/ttZ2l2nuFiles.txt
ls /scratch/data/tZqSkimsRun2015/ttZ2q/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/ttZ2qFiles.txt
ls /scratch/data/tZqSkimsRun2015/tW_top_5f/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/tWInclusiveFiles.txt
ls /nfs/data/tZqSkims/tZq4Flavour3Lepton/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/tZqFilesRun1.txt
ls /scratch/data/tZqSkimsRun2015/tZq_ll_4Flavour3Lepton/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/tZqFiles.txt
ls /scratch/data/tZqSkimsRun2015/tHq/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/tHqFiles.txt
ls /scratch/data/tZqSkimsRun2015/wPlusJets/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/wPlusJetsFiles.txt
ls /scratch/data/tZqSkimsRun2015/WW1l1nu2q/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/WW1l1nu2qFiles.txt
ls /scratch/data/tZqSkimsRun2015/WW2l2nu/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/WW2l2nuFiles.txt
ls /scratch/data/tZqSkimsRun2015/WZ1l1nu2q/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/WZ1l1nu2q.txt
ls /scratch/data/tZqSkimsRun2015/WZ2l2q/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/WZ2l2q.txt
ls /scratch/data/tZqSkimsRun2015/WZJets/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/WZjets.txt
ls /scratch/data/tZqSkimsRun2015/ZZ2l2nu/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/ZZ2l2nuFiles.txt
ls /scratch/data/tZqSkimsRun2015/ZZ2l2q/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/ZZ2l2qFiles.txt
ls /scratch/data/tZqSkimsRun2015/ZZ4l/* -1d >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/ZZ4lFiles.txt

# Scaling files
ls /scratch/data/tZqSkimsRun2015/tW_antitop_5f_scaleup/* >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/tbarWInclusiveScaleUpFiles.txt
ls /scratch/data/tZqSkimsRun2015/tW_antitop_5f_scaledown/* >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/tbarWInclusiveScaleDownFiles.txt
ls /scratch/data/tZqSkimsRun2015/tW_top_5f_scaleup/* >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/tWInclusiveScaleUpFiles.txt
ls /scratch/data/tZqSkimsRun2015/tW_top_5f_scaleup/* >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/tWInclusiveScaleDownFiles.txt

# Synchronisation files
ls /scratch/data/tZqSkimsRun2015/synch/tZq/* -1d  >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/tZqSynchFiles.txt
ls /scratch/data/tZqSkimsRun2015/synch/WZ/* -1d  >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/WZsynchFiles.txt
ls /scratch/data/tZqSkimsRun2015/synch/ttZ/* -1d  >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/ttZsynchFiles.txt
ls /scratch/data/tZqSkimsRun2015/synch/ttbar/* -1d  >> $TQZ_TOOLS_PATH/configs/2015/datasets/fileLists/ttbarSynchFiles.txt

echo Done!
echo Filelists have been updated.
