#!/bin/bash
echo Updating filelists for datasets used ...
echo First deleting old filelists ...

if [ -z "$TQZ_TOOLS_PATH" ]; then
    echo '$TQZ_TOOLS_PATH not set, setting to .'
    export TQZ_TOOLS_PATH='.'
fi

rm $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/*

echo Done!
echo Now outputting the lists of the dataset directories into their relevant files ...


# Normal datasets

ls /data0/data/TopPhysics/postTriggerSkims2016/eeRun2016B/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/eeRun2016Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/eeRun2016C/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/eeRun2016Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/eeRun2016D/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/eeRun2016Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/eeRun2016E/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/eeRun2016Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/eeRun2016F/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/eeRun2016Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/eeRun2016G/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/eeRun2016Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/eeRun2016H/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/eeRun2016Files.txt

ls /data0/data/TopPhysics/postTriggerSkims2016/emuRun2016B/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/emuRun2016Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/emuRun2016C/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/emuRun2016Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/emuRun2016D/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/emuRun2016Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/emuRun2016E/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/emuRun2016Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/emuRun2016F/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/emuRun2016Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/emuRun2016G/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/emuRun2016Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/emuRun2016H/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/emuRun2016Files.txt

ls /data0/data/TopPhysics/postTriggerSkims2016/mumuRun2016B/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/mumuRun2016Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/mumuRun2016C/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/mumuRun2016Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/mumuRun2016D/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/mumuRun2016Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/mumuRun2016E/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/mumuRun2016Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/mumuRun2016F/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/mumuRun2016Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/mumuRun2016G/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/mumuRun2016Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/mumuRun2016H/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/mumuRun2016Files.txt

ls /data0/data/TopPhysics/postTriggerSkims2016/mumuRun2016B/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/mumuRun2016FilesPart1.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/mumuRun2016C/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/mumuRun2016FilesPart1.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/mumuRun2016D/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/mumuRun2016FilesPart1.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/mumuRun2016E/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/mumuRun2016FilesPart1.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/mumuRun2016F/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/mumuRun2016FilesPart1.txt

ls /data0/data/TopPhysics/postTriggerSkims2016/mumuRun2016G/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/mumuRun2016FilesPart2.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/mumuRun2016H/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/mumuRun2016FilesPart2.txt

ls /data0/data/TopPhysics/postTriggerSkims2016/mumuRun2016B/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/mumuRun2016BFiles.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/mumuRun2016C/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/mumuRun2016CFiles.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/mumuRun2016D/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/mumuRun2016DFiles.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/mumuRun2016E/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/mumuRun2016EFiles.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/mumuRun2016F/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/mumuRun2016FFiles.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/mumuRun2016G/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/mumuRun2016GFiles.txt
ls /data0/data/TopPhysics/postTriggerSkims2016/mumuRun2016H/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/mumuRun2016HFiles.txt

ls /scratch/data/tZqSkimsRun2016/metRun2016B/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/metRun2016Files.txt
ls /scratch/data/tZqSkimsRun2016/metRun2016C/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/metRun2016Files.txt
ls /scratch/data/tZqSkimsRun2016/metRun2016D/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/metRun2016Files.txt
ls /scratch/data/tZqSkimsRun2016/metRun2016E/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/metRun2016Files.txt
ls /scratch/data/tZqSkimsRun2016/metRun2016F/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/metRun2016Files.txt
ls /scratch/data/tZqSkimsRun2016/metRun2016G/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/metRun2016Files.txt
ls /scratch/data/tZqSkimsRun2016/metRun2016H/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/metRun2016Files.txt

ls /scratch/data/tZqSkimsRun2016/metRun2016B/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/metRun2016Part1Files.txt
ls /scratch/data/tZqSkimsRun2016/metRun2016C/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/metRun2016Part1Files.txt
ls /scratch/data/tZqSkimsRun2016/metRun2016D/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/metRun2016Part1Files.txt
ls /scratch/data/tZqSkimsRun2016/metRun2016E/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/metRun2016Part1Files.txt
ls /scratch/data/tZqSkimsRun2016/metRun2016F/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/metRun2016Part1Files.txt
ls /scratch/data/tZqSkimsRun2016/metRun2016G/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/metRun2016Part2Files.txt
ls /scratch/data/tZqSkimsRun2016/metRun2016H/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/metRun2016Part2Files.txt

ls /scratch/data/tZqSkimsRun2016/metRun2016B/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/metRun2016BFiles.txt
ls /scratch/data/tZqSkimsRun2016/metRun2016C/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/metRun2016CFiles.txt
ls /scratch/data/tZqSkimsRun2016/metRun2016D/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/metRun2016DFiles.txt
ls /scratch/data/tZqSkimsRun2016/metRun2016E/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/metRun2016EFiles.txt
ls /scratch/data/tZqSkimsRun2016/metRun2016F/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/metRun2016FFiles.txt
ls /scratch/data/tZqSkimsRun2016/metRun2016G/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/metRun2016GFiles.txt
ls /scratch/data/tZqSkimsRun2016/metRun2016H/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/metRun2016HFiles.txt

ls /scratch/data/tZqSkimsRun2016/TGJets/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/TGJetsFiles.txt
ls /scratch/data/tZqSkimsRun2016/TTGJets/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/TTGJetsFiles.txt
ls /scratch/data/tZqSkimsRun2016/TTGJets_ext/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/TTGJetsFiles.txt

ls /scratch/data/tZqSkimsRun2016/GJet_Pt-20to40_MGG-80toInf/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/GJetPt20to40_MGG80toInfFiles.txt
ls /scratch/data/tZqSkimsRun2016/GJet_Pt-20toInf_MGG-40to80/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/GJetPt20toInf_MGG40to80Files.txt
ls /scratch/data/tZqSkimsRun2016/GJet_Pt-40toInf_MGG-80toInf/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/GJetPt40toInf_MGG80toInfFiles.txt

ls /scratch/data/tZqSkimsRun2016/ttbarttbar/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbarttbarFiles.txt
ls /scratch/data/tZqSkimsRun2016/WWG/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/WWGFiles.txt
ls /scratch/data/tZqSkimsRun2016/WZG/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/WZGFiles.txt
ls /scratch/data/tZqSkimsRun2016/WGGJets/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/WGGJetsFiles.txt
ls /scratch/data/tZqSkimsRun2016/WW4q/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/WW4qFiles.txt

ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_M-10to50_amcatnlo/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets10To50Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_M-10to50_amcatnlo_ext/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets10To50Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_M-10to50_madgraph/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets10To50MadgraphFiles.txt

ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_M-50_amcatnlo_ext2/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets50Files.txt

ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_M-50_amcatnlo_ext2/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets50Part2Files.txt

ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_M-50_madgraph/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets50MadgraphFiles.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_M-50_madgraph_ext2/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets50MadgraphFiles.txt

ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_M-50_madgraph/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets50MadgraphPart1Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_M-50_madgraph_ext2/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets50MadgraphPart2Files.txt

ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_Pt-0To50/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJetsPt0To50Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_Pt-50To100/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJetsPt50To100Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_Pt-50To100_ext3/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJetsPt50To100Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_Pt-100To250/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJetsPt100To250Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_Pt-100To250_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJetsPt100To250Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_Pt-100To250_ext2/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJetsPt100To250Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_Pt-100To250_ext5/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJetsPt100To250Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_Pt-250To400/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJetsPt250To400Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_Pt-250To400_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJetsPt250To400Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_Pt-250To400_ext2/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJetsPt250To400Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_Pt-250To400_ext5/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJetsPt250To400Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_Pt-400To650/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJetsPt400To650Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_Pt-400To650_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJetsPt400To650Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_Pt-400To650_ext2/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJetsPt400To650Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_Pt-650ToInf/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJetsPt650ToInfFiles.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_Pt-650ToInf_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJetsPt650ToInfFiles.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_Pt-650ToInf_ext2/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJetsPt650ToInfFiles.txt

ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_M-100to200_amcatnlo_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets100To200Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_M-100to200_amcatnlo_ext3/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets100To200Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_M-200to400_amcatnlo_ext2/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets200To400Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_M-200to400_amcatnlo_ext3/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets200To400Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_M-400to500_amcatnlo_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets400To500Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_M-500to700_amcatnlo_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets500To700Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_M-700to800_amcatnlo_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets700To800Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_M-800to1000_amcatnlo_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets800To1000Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_M-1000to1500_amcatnlo_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets1000To1500Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_M-1500to2000_amcatnlo_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets1500To2000Files.txt
ls /scratch/data/tZqSkimsRun2016/DYJetsToLL_M-2000to3000_amcatnlo_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DYJets2000To3000Files.txt

ls /scratch/data/tZqSkimsRun2016/DYToLL_0J_amcatnlo/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DY0JetsFiles.txt
ls /scratch/data/tZqSkimsRun2016/DYToLL_1J_amcatnlo/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DY1JetsFiles.txt
ls /scratch/data/tZqSkimsRun2016/DYToLL_2J_amcatnlo/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DY2JetsFiles.txt
ls /scratch/data/tZqSkimsRun2016/DYToLL_2J_amcatnlo_ext/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DY2JetsFiles.txt

ls /scratch/data/tZqSkimsRun2016/DY1JetsToLL_M-10to50_madgraph/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DY1Jets10to50MadgraphFiles.txt
ls /scratch/data/tZqSkimsRun2016/DY2JetsToLL_M-10to50_madgraph/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DY2Jets10to50MadgraphFiles.txt
ls /scratch/data/tZqSkimsRun2016/DY3JetsToLL_M-10to50_madgraph/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DY3Jets10to50MadgraphFiles.txt
ls /scratch/data/tZqSkimsRun2016/DY4JetsToLL_M-10to50_madgraph/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DY4Jets10to50MadgraphFiles.txt

ls /scratch/data/tZqSkimsRun2016/DY1JetsToLL_M-50_madgraph/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DY1Jets50MadgraphFiles.txt
ls /scratch/data/tZqSkimsRun2016/DY2JetsToLL_M-50_madgraph/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DY2Jets50MadgraphFiles.txt
ls /scratch/data/tZqSkimsRun2016/DY3JetsToLL_M-50_madgraph/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DY3Jets50MadgraphFiles.txt
ls /scratch/data/tZqSkimsRun2016/DY4JetsToLL_M-50_madgraph/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/DY4Jets50MadgraphFiles.txt

ls /scratch/data/tZqSkimsRun2016/sChannel_4f/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/sChannelFiles.txt
ls /scratch/data/tZqSkimsRun2016/tWZ_5f/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tWZ_ll_Files.txt

ls /scratch/data/tZqSkimsRun2016/ttbarInclusive_powerheg/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbarInclusivePowerhegFiles.txt

ls /scratch/data/tZqSkimsRun2016/ttbarInclusive_powerheg_colourFlip/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbarInclusivePowerhegColourFlipFiles.txt
ls /scratch/data/tZqSkimsRun2016/ttbarInclusive_powerheg_hdampUP/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbarInclusivePowerhegHdampUpFiles.txt
ls /scratch/data/tZqSkimsRun2016/ttbarInclusive_powerheg_hdampUP_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbarInclusivePowerhegHdampUpFiles.txt
ls /scratch/data/tZqSkimsRun2016/ttbarInclusive_powerheg_hdampDOWN/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbarInclusivePowerhegHdampDownFiles.txt
ls /scratch/data/tZqSkimsRun2016/ttbarInclusive_powerheg_hdampDOWN_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbarInclusivePowerhegHdampDownFiles.txt
ls /scratch/data/tZqSkimsRun2016/ttbarInclusive_powerheg_fsrup/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbarInclusivePowerhegFsrUpFiles.txt
ls /scratch/data/tZqSkimsRun2016/ttbarInclusive_powerheg_fsrup_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbarInclusivePowerhegFsrUpFiles.txt
ls /scratch/data/tZqSkimsRun2016/ttbarInclusive_powerheg_fsrup_ext2/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbarInclusivePowerhegFsrUpFiles.txt
ls /scratch/data/tZqSkimsRun2016/ttbarInclusive_powerheg_fsrdown/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbarInclusivePowerhegFsrDownFiles.txt
ls /scratch/data/tZqSkimsRun2016/ttbarInclusive_powerheg_fsrdown_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbarInclusivePowerhegFsrDownFiles.txt
ls /scratch/data/tZqSkimsRun2016/ttbarInclusive_powerheg_fsrdown_ext2/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbarInclusivePowerhegFsrDownFiles.txt
ls /scratch/data/tZqSkimsRun2016/ttbarInclusive_powerheg_isrup_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbarInclusivePowerhegIsrUpFiles.txt
ls /scratch/data/tZqSkimsRun2016/ttbarInclusive_powerheg_isrup_ext2/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbarInclusivePowerhegIsrUpFiles.txt
ls /scratch/data/tZqSkimsRun2016/ttbarInclusive_powerheg_isrdown/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbarInclusivePowerhegIsrDownFiles.txt
ls /scratch/data/tZqSkimsRun2016/ttbarInclusive_powerheg_isrdown_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbarInclusivePowerhegIsrDownFiles.txt
ls /scratch/data/tZqSkimsRun2016/ttbarInclusive_powerheg_isrdown_ext2/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbarInclusivePowerhegIsrDownFiles.txt

ls /scratch/data/tZqSkimsRun2016/ttHTobb/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttHTobbFiles.txt
ls /scratch/data/tZqSkimsRun2016/ttHToNonbb/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttHToNonbbFiles.txt

ls /scratch/data/tZqSkimsRun2016/ttWlnu/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttWlnuFiles.txt
ls /scratch/data/tZqSkimsRun2016/ttW2q/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttW2qFiles.txt
ls /scratch/data/tZqSkimsRun2016/ttZ2l2nu/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttZ2l2nuFiles.txt
ls /scratch/data/tZqSkimsRun2016/ttZ2q/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttZ2qFiles.txt

ls /scratch/data/tZqSkimsRun2016/tChannel_4f/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tChannelFiles.txt
ls /scratch/data/tZqSkimsRun2016/tChannel_4f_scaleup/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tChannelScaleUpFiles.txt
ls /scratch/data/tZqSkimsRun2016/tChannel_4f_scaledown/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tChannelScaleDownFiles.txt
ls /scratch/data/tZqSkimsRun2016/tChannel_4f_hdampup/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tChannelHdampUpFiles.txt
ls /scratch/data/tZqSkimsRun2016/tChannel_4f_hdampdown/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tChannelHdampDownFiles.txt
ls /scratch/data/tZqSkimsRun2016/tbarChannel_4f/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tbarChannelFiles.txt
ls /scratch/data/tZqSkimsRun2016/tbarChannel_4f_scaleup/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tbarChannelScaleUpFiles.txt
ls /scratch/data/tZqSkimsRun2016/tbarChannel_4f_scaledown/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tbarChannelScaleDownFiles.txt
ls /scratch/data/tZqSkimsRun2016/tbarChannel_4f_hdampup/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tbarChannelHdampUpFiles.txt
ls /scratch/data/tZqSkimsRun2016/tbarChannel_4f_hdampdown/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tbarChannelHdampDownFiles.txt

ls /scratch/data/tZqSkimsRun2016/tW_antitop_5f/* -1d  >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tbarWInclusiveFiles.txt
ls /scratch/data/tZqSkimsRun2016/tW_antitop_5f_scaledown/* -1d  >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tbarWInclusiveScaleDownFiles.txt
ls /scratch/data/tZqSkimsRun2016/tW_antitop_5f_scaleup/* -1d  >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tbarWInclusiveScaleUpFiles.txt
ls /scratch/data/tZqSkimsRun2016/tW_top_5f/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tWInclusiveFiles.txt
ls /scratch/data/tZqSkimsRun2016/tW_top_5f_scaledown/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tWInclusiveScaleDownFiles.txt
ls /scratch/data/tZqSkimsRun2016/tW_top_5f_scaleup/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tWInclusiveScaleUpFiles.txt

ls /scratch/data/tZqSkimsRun2016/tZq_ll_4Flavour3Lepton/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tZqFiles.txt
ls /scratch/data/tZqSkimsRun2016/tZq_ll_4Flavour3Lepton_scaleup/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tZqScaleUpFiles.txt
ls /scratch/data/tZqSkimsRun2016/tZq_ll_4Flavour3Lepton_scaledown/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tZqScaleDownFiles.txt

ls /scratch/data/tZqSkimsRun2016/tHq/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tHqFiles.txt
ls /scratch/data/tZqSkimsRun2016/wPlusJets/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/wPlusJetsFiles.txt
ls /scratch/data/tZqSkimsRun2016/wPlusJets_ext2/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/wPlusJetsFiles.txt
ls /scratch/data/tZqSkimsRun2016/WW1l1nu2q/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/WW1l1nu2qFiles.txt
ls /scratch/data/tZqSkimsRun2016/WW2l2nu/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/WW2l2nuFiles.txt
ls /scratch/data/tZqSkimsRun2016/WZ1l1nu2q/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/WZ1l1nu2q.txt
ls /scratch/data/tZqSkimsRun2016/WZ2l2q/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/WZ2l2q.txt
ls /scratch/data/tZqSkimsRun2016/WZJets/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/WZjets.txt
ls /scratch/data/tZqSkimsRun2016/ZZ2l2nu/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ZZ2l2nuFiles.txt
ls /scratch/data/tZqSkimsRun2016/ZZ2l2q/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ZZ2l2qFiles.txt
ls /scratch/data/tZqSkimsRun2016/ZZ4l/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ZZ4lFiles.txt

ls /scratch/data/tZqSkimsRun2016/WWW/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/WWWFiles.txt
ls /scratch/data/tZqSkimsRun2016/WWZ/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/WWZFiles.txt
ls /scratch/data/tZqSkimsRun2016/WZZ/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/WZZFiles.txt
ls /scratch/data/tZqSkimsRun2016/ZZZ/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ZZZFiles.txt

# FCNC Datasets

ls /scratch/data/tZqSkimsRun2016/tZq_lll_Kappa_Zct/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tZq_lll_Kappa_Zct.txt
ls /scratch/data/tZqSkimsRun2016/tZq_lll_Kappa_Zut/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tZq_lll_Kappa_Zut.txt
ls /scratch/data/tZqSkimsRun2016/tZq_lll_Zeta_Zct/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tZq_lll_Zeta_Zct.txt

ls /scratch/data/tZqSkimsRun2016/tZq_ll_Kappa_Zct/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tZq_ll_Kappa_Zct.txt
ls /scratch/data/tZqSkimsRun2016/tZq_ll_Kappa_Zut/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tZq_ll_Kappa_Zut.txt
ls /scratch/data/tZqSkimsRun2016/tZq_ll_Zeta_Zct/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tZq_ll_Zeta_Zct.txt
ls /scratch/data/tZqSkimsRun2016/tZq_ll_Zeta_Zut/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tZq_ll_Zeta_Zut.txt

ls /scratch/data/tZqSkimsRun2016/ttbar_tLeptonic_FCNC_Kappa_Zct/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbar_tLeptonic_FCNC_Kappa_Zct.txt
ls /scratch/data/tZqSkimsRun2016/ttbar_tLeptonic_FCNC_Kappa_Zut/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbar_tLeptonic_FCNC_Kappa_Zut.txt
ls /scratch/data/tZqSkimsRun2016/ttbar_tbarLeptonic_FCNC_Kappa_Zct/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbar_tbarLeptonic_FCNC_Kappa_Zct.txt
ls /scratch/data/tZqSkimsRun2016/ttbar_tbarLeptonic_FCNC_Kappa_Zut/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbar_tbarLeptonic_FCNC_Kappa_Zut.txt

ls /scratch/data/tZqSkimsRun2016/ttbar_tLeptonic_FCNC_Zeta_Zct/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbar_tLeptonic_FCNC_Zeta_Zct.txt
ls /scratch/data/tZqSkimsRun2016/ttbar_tLeptonic_FCNC_Zeta_Zut/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbar_tLeptonic_FCNC_Zeta_Zut.txt
ls /scratch/data/tZqSkimsRun2016/ttbar_tbarLeptonic_FCNC_Zeta_Zct/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbar_tbarLeptonic_FCNC_Zeta_Zct.txt
ls /scratch/data/tZqSkimsRun2016/ttbar_tbarLeptonic_FCNC_Zeta_Zut/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbar_tbarLeptonic_FCNC_Zeta_Zut.txt

ls /scratch/data/tZqSkimsRun2016/ttbar_tHadronic_FCNC_Kappa_Zct/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbar_tHadronic_FCNC_Kappa_Zct.txt
ls /scratch/data/tZqSkimsRun2016/ttbar_tHadronic_FCNC_Kappa_Zut/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbar_tHadronic_FCNC_Kappa_Zut.txt
ls /scratch/data/tZqSkimsRun2016/ttbar_tbarHadronic_FCNC_Kappa_Zct/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbar_tbarHadronic_FCNC_Kappa_Zct.txt
ls /scratch/data/tZqSkimsRun2016/ttbar_tbarHadronic_FCNC_Kappa_Zut/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbar_tbarHadronic_FCNC_Kappa_Zut.txt

ls /scratch/data/tZqSkimsRun2016/ttbar_tHadronic_FCNC_Zeta_Zct/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbar_tHadronic_FCNC_Zeta_Zct.txt
ls /scratch/data/tZqSkimsRun2016/ttbar_tHadronic_FCNC_Zeta_Zut/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbar_tHadronic_FCNC_Zeta_Zut.txt
ls /scratch/data/tZqSkimsRun2016/ttbar_tbarHadronic_FCNC_Zeta_Zct/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbar_tbarHadronic_FCNC_Zeta_Zct.txt
ls /scratch/data/tZqSkimsRun2016/ttbar_tbarHadronic_FCNC_Zeta_Zut/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbar_tbarHadronic_FCNC_Zeta_Zut.txt

# Synchronisation files
ls /scratch/data/tZqSkimsRun2016/synch/tZq/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/tZqSynchFiles.txt
ls /scratch/data/tZqSkimsRun2016/synch/data1/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/data1SynchFiles.txt
ls /scratch/data/tZqSkimsRun2016/synch/data2/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/data2SynchFiles.txt

ls /scratch/data/tZqSkimsRun2016/synch/tW/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/twSynchFiles.txt
ls /scratch/data/tZqSkimsRun2016/synch/ttbar/* -1d >> $TQZ_TOOLS_PATH/configs/2016/datasets/fileLists/ttbarSynchFiles.txt

echo Done!
echo Filelists have been updated.
