#!/bin/bash
echo Updating filelists for datasets used ...
echo First deleting old filelists ...

if [ -z "$TQZ_TOOLS_PATH" ]; then
    echo '$TQZ_TOOLS_PATH not set, setting to .'
    export TQZ_TOOLS_PATH='.'
fi

rm $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/*

echo Done!
echo Now outputting the lists of the dataset directories into their relevant files ...


# Normal datasets

ls /data0/data/TopPhysics/postTriggerSkims2017/eeRun2017B/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/eeRun2017Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2017/eeRun2017C/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/eeRun2017Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2017/eeRun2017D/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/eeRun2017Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2017/eeRun2017E/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/eeRun2017Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2017/eeRun2017F/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/eeRun2017Files.txt

ls /data0/data/TopPhysics/postTriggerSkims2017/emuRun2017B/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/emuRun2017Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2017/emuRun2017C/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/emuRun2017Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2017/emuRun2017D/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/emuRun2017Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2017/emuRun2017E/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/emuRun2017Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2017/emuRun2017F/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/emuRun2017Files.txt

ls /data0/data/TopPhysics/postTriggerSkims2017/mumuRun2017B/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/mumuRun2017Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2017/mumuRun2017C/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/mumuRun2017Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2017/mumuRun2017D/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/mumuRun2017Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2017/mumuRun2017E/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/mumuRun2017Files.txt
ls /data0/data/TopPhysics/postTriggerSkims2017/mumuRun2017F/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/mumuRun2017Files.txt

ls /data0/data/TopPhysics/postTriggerSkims2017/mumuRun2017B/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/mumuRun2017FilesPart1.txt
ls /data0/data/TopPhysics/postTriggerSkims2017/mumuRun2017C/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/mumuRun2017FilesPart1.txt
ls /data0/data/TopPhysics/postTriggerSkims2017/mumuRun2017D/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/mumuRun2017FilesPart1.txt
ls /data0/data/TopPhysics/postTriggerSkims2017/mumuRun2017E/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/mumuRun2017FilesPart1.txt
ls /data0/data/TopPhysics/postTriggerSkims2017/mumuRun2017F/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/mumuRun2017FilesPart1.txt


ls /data0/data/TopPhysics/postTriggerSkims2017/mumuRun2017B/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/mumuRun2017BFiles.txt
ls /data0/data/TopPhysics/postTriggerSkims2017/mumuRun2017C/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/mumuRun2017CFiles.txt
ls /data0/data/TopPhysics/postTriggerSkims2017/mumuRun2017D/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/mumuRun2017DFiles.txt
ls /data0/data/TopPhysics/postTriggerSkims2017/mumuRun2017E/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/mumuRun2017EFiles.txt
ls /data0/data/TopPhysics/postTriggerSkims2017/mumuRun2017F/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/mumuRun2017FFiles.txt

ls /scratch/data/tZqSkimsRun2017/metRun2017B/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/metRun2017Files.txt
ls /scratch/data/tZqSkimsRun2017/metRun2017C/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/metRun2017Files.txt
ls /scratch/data/tZqSkimsRun2017/metRun2017D/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/metRun2017Files.txt
ls /scratch/data/tZqSkimsRun2017/metRun2017E/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/metRun2017Files.txt
ls /scratch/data/tZqSkimsRun2017/metRun2017F/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/metRun2017Files.txt

ls /scratch/data/tZqSkimsRun2017/metRun2017B/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/metRun2017Part1Files.txt
ls /scratch/data/tZqSkimsRun2017/metRun2017C/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/metRun2017Part1Files.txt
ls /scratch/data/tZqSkimsRun2017/metRun2017D/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/metRun2017Part1Files.txt
ls /scratch/data/tZqSkimsRun2017/metRun2017E/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/metRun2017Part1Files.txt
ls /scratch/data/tZqSkimsRun2017/metRun2017F/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/metRun2017Part1Files.txt

ls /scratch/data/tZqSkimsRun2017/metRun2017B/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/metRun2017BFiles.txt
ls /scratch/data/tZqSkimsRun2017/metRun2017C/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/metRun2017CFiles.txt
ls /scratch/data/tZqSkimsRun2017/metRun2017D/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/metRun2017DFiles.txt
ls /scratch/data/tZqSkimsRun2017/metRun2017E/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/metRun2017EFiles.txt
ls /scratch/data/tZqSkimsRun2017/metRun2017F/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/metRun2017FFiles.txt

ls /scratch/data/tZqSkimsRun2017/TGJets/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/TGJetsFiles.txt
ls /scratch/data/tZqSkimsRun2017/TTGJets/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/TTGJetsFiles.txt
ls /scratch/data/tZqSkimsRun2017/TTGJets_ext/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/TTGJetsFiles.txt

ls /scratch/data/tZqSkimsRun2017/GJet_Pt-20to40_MGG-80toInf/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/GJetPt20to40_MGG80toInfFiles.txt
ls /scratch/data/tZqSkimsRun2017/GJet_Pt-20toInf_MGG-40to80/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/GJetPt20toInf_MGG40to80Files.txt
ls /scratch/data/tZqSkimsRun2017/GJet_Pt-40toInf_MGG-80toInf/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/GJetPt40toInf_MGG80toInfFiles.txt

ls /scratch/data/tZqSkimsRun2017/ttbarttbar/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttbarttbarFiles.txt
ls /scratch/data/tZqSkimsRun2017/WWG/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/WWGFiles.txt
ls /scratch/data/tZqSkimsRun2017/WZG/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/WZGFiles.txt
ls /scratch/data/tZqSkimsRun2017/WGGJets/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/WGGJetsFiles.txt
ls /scratch/data/tZqSkimsRun2017/WW4q/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/WW4qFiles.txt

ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_M-10to50_amcatnlo/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJets10To50Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_M-10to50_amcatnlo_ext/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJets10To50Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_M-10to50_madgraph/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJets10To50MadgraphFiles.txt

ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_M-50_amcatnlo/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJets50Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_M-50_amcatnlo_ext2/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJets50Files.txt

ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_M-50_amcatnlo/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJets50Part1Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_M-50_amcatnlo_ext2/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJets50Part2Files.txt

ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_M-50_madgraph/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJets50MadgraphFiles.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_M-50_madgraph_ext2/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJets50MadgraphFiles.txt

ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_M-50_madgraph/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJets50MadgraphPart1Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_M-50_madgraph_ext2/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJets50MadgraphPart2Files.txt

ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_Pt-0To50/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJetsPt0To50Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_Pt-50To100/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJetsPt50To100Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_Pt-50To100_ext3/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJetsPt50To100Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_Pt-100To250/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJetsPt100To250Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_Pt-100To250_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJetsPt100To250Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_Pt-100To250_ext2/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJetsPt100To250Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_Pt-100To250_ext5/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJetsPt100To250Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_Pt-250To400/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJetsPt250To400Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_Pt-250To400_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJetsPt250To400Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_Pt-250To400_ext2/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJetsPt250To400Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_Pt-250To400_ext5/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJetsPt250To400Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_Pt-400To650/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJetsPt400To650Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_Pt-400To650_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJetsPt400To650Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_Pt-400To650_ext2/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJetsPt400To650Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_Pt-650ToInf/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJetsPt650ToInfFiles.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_Pt-650ToInf_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJetsPt650ToInfFiles.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_Pt-650ToInf_ext2/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJetsPt650ToInfFiles.txt

ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_M-100to200_amcatnlo_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJets100To200Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_M-100to200_amcatnlo_ext3/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJets100To200Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_M-200to400_amcatnlo_ext2/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJets200To400Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_M-200to400_amcatnlo_ext3/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJets200To400Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_M-400to500_amcatnlo_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJets400To500Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_M-500to700_amcatnlo_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJets500To700Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_M-700to800_amcatnlo_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJets700To800Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_M-800to1000_amcatnlo_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJets800To1000Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_M-1000to1500_amcatnlo_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJets1000To1500Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_M-1500to2000_amcatnlo_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJets1500To2000Files.txt
ls /scratch/data/tZqSkimsRun2017/DYJetsToLL_M-2000to3000_amcatnlo_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DYJets2000To3000Files.txt

ls /scratch/data/tZqSkimsRun2017/DYToLL_0J_amcatnlo/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DY0JetsFiles.txt
ls /scratch/data/tZqSkimsRun2017/DYToLL_1J_amcatnlo/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DY1JetsFiles.txt
ls /scratch/data/tZqSkimsRun2017/DYToLL_2J_amcatnlo/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DY2JetsFiles.txt
ls /scratch/data/tZqSkimsRun2017/DYToLL_2J_amcatnlo_ext/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DY2JetsFiles.txt

ls /scratch/data/tZqSkimsRun2017/DY1JetsToLL_M-10to50_madgraph/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DY1Jets10to50MadgraphFiles.txt
ls /scratch/data/tZqSkimsRun2017/DY2JetsToLL_M-10to50_madgraph/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DY2Jets10to50MadgraphFiles.txt
ls /scratch/data/tZqSkimsRun2017/DY3JetsToLL_M-10to50_madgraph/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DY3Jets10to50MadgraphFiles.txt
ls /scratch/data/tZqSkimsRun2017/DY4JetsToLL_M-10to50_madgraph/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DY4Jets10to50MadgraphFiles.txt

ls /scratch/data/tZqSkimsRun2017/DY1JetsToLL_M-50_madgraph/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DY1Jets50MadgraphFiles.txt
ls /scratch/data/tZqSkimsRun2017/DY2JetsToLL_M-50_madgraph/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DY2Jets50MadgraphFiles.txt
ls /scratch/data/tZqSkimsRun2017/DY3JetsToLL_M-50_madgraph/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DY3Jets50MadgraphFiles.txt
ls /scratch/data/tZqSkimsRun2017/DY4JetsToLL_M-50_madgraph/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/DY4Jets50MadgraphFiles.txt

ls /scratch/data/tZqSkimsRun2017/sChannel_4f/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/sChannelFiles.txt
ls /scratch/data/tZqSkimsRun2017/tWZ_5f/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/tWZ_ll_Files.txt
ls /scratch/data/tZqSkimsRun2017/ttbarDilepton_madgraph/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttbarDileptonFiles.txt
ls /scratch/data/tZqSkimsRun2017/ttbarDilepton_madgraph_ext/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttbarDileptonFiles.txt
ls /scratch/data/tZqSkimsRun2017/ttbarDilepton_aMCatNLO/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttbarDileptonAmcAtNloFiles.txt
ls /scratch/data/tZqSkimsRun2017/ttbarDilepton_aMCatNLO_ext/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttbarDileptonAmcAtNloFiles.txt

ls /scratch/data/tZqSkimsRun2017/ttbarInclusive_powerheg/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttbarInclusivePowerhegFiles.txt

ls /scratch/data/tZqSkimsRun2017/ttbarInclusive_powerheg_colourFlip/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttbarInclusivePowerhegColourFlipFiles.txt
ls /scratch/data/tZqSkimsRun2017/ttbarInclusive_powerheg_hdampUP/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttbarInclusivePowerhegHdampUpFiles.txt
ls /scratch/data/tZqSkimsRun2017/ttbarInclusive_powerheg_hdampUP_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttbarInclusivePowerhegHdampUpFiles.txt
ls /scratch/data/tZqSkimsRun2017/ttbarInclusive_powerheg_hdampDOWN/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttbarInclusivePowerhegHdampDownFiles.txt
ls /scratch/data/tZqSkimsRun2017/ttbarInclusive_powerheg_hdampDOWN_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttbarInclusivePowerhegHdampDownFiles.txt
ls /scratch/data/tZqSkimsRun2017/ttbarInclusive_powerheg_fsrup/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttbarInclusivePowerhegFsrUpFiles.txt
ls /scratch/data/tZqSkimsRun2017/ttbarInclusive_powerheg_fsrup_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttbarInclusivePowerhegFsrUpFiles.txt
ls /scratch/data/tZqSkimsRun2017/ttbarInclusive_powerheg_fsrup_ext2/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttbarInclusivePowerhegFsrUpFiles.txt
ls /scratch/data/tZqSkimsRun2017/ttbarInclusive_powerheg_fsrdown/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttbarInclusivePowerhegFsrDownFiles.txt
ls /scratch/data/tZqSkimsRun2017/ttbarInclusive_powerheg_fsrdown_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttbarInclusivePowerhegFsrDownFiles.txt
ls /scratch/data/tZqSkimsRun2017/ttbarInclusive_powerheg_fsrdown_ext2/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttbarInclusivePowerhegFsrDownFiles.txt
ls /scratch/data/tZqSkimsRun2017/ttbarInclusive_powerheg_isrup_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttbarInclusivePowerhegIsrUpFiles.txt
ls /scratch/data/tZqSkimsRun2017/ttbarInclusive_powerheg_isrup_ext2/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttbarInclusivePowerhegIsrUpFiles.txt
ls /scratch/data/tZqSkimsRun2017/ttbarInclusive_powerheg_isrdown/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttbarInclusivePowerhegIsrDownFiles.txt
ls /scratch/data/tZqSkimsRun2017/ttbarInclusive_powerheg_isrdown_ext1/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttbarInclusivePowerhegIsrDownFiles.txt
ls /scratch/data/tZqSkimsRun2017/ttbarInclusive_powerheg_isrdown_ext2/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttbarInclusivePowerhegIsrDownFiles.txt

ls /scratch/data/tZqSkimsRun2017/ttHTobb/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttHTobbFiles.txt
ls /scratch/data/tZqSkimsRun2017/ttHToNonbb/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttHToNonbbFiles.txt

ls /scratch/data/tZqSkimsRun2017/ttWlnu/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttWlnuFiles.txt
ls /scratch/data/tZqSkimsRun2017/ttW2q/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttW2qFiles.txt
ls /scratch/data/tZqSkimsRun2017/ttZ2l2nu/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttZ2l2nuFiles.txt
ls /scratch/data/tZqSkimsRun2017/ttZ2q/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttZ2qFiles.txt

ls /scratch/data/tZqSkimsRun2017/tChannel_4f/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/tChannelFiles.txt
ls /scratch/data/tZqSkimsRun2017/tChannel_4f_scaleup/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/tChannelScaleUpFiles.txt
ls /scratch/data/tZqSkimsRun2017/tChannel_4f_scaledown/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/tChannelScaleDownFiles.txt
ls /scratch/data/tZqSkimsRun2017/tChannel_4f_hdampup/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/tChannelHdampUpFiles.txt
ls /scratch/data/tZqSkimsRun2017/tChannel_4f_hdampdown/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/tChannelHdampDownFiles.txt
ls /scratch/data/tZqSkimsRun2017/tbarChannel_4f/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/tbarChannelFiles.txt
ls /scratch/data/tZqSkimsRun2017/tbarChannel_4f_scaleup/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/tbarChannelScaleUpFiles.txt
ls /scratch/data/tZqSkimsRun2017/tbarChannel_4f_scaledown/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/tbarChannelScaleDownFiles.txt
ls /scratch/data/tZqSkimsRun2017/tbarChannel_4f_hdampup/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/tbarChannelHdampUpFiles.txt
ls /scratch/data/tZqSkimsRun2017/tbarChannel_4f_hdampdown/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/tbarChannelHdampDownFiles.txt

ls /scratch/data/tZqSkimsRun2017/tW_antitop_5f/* -1d  >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/tbarWInclusiveFiles.txt
ls /scratch/data/tZqSkimsRun2017/tW_antitop_5f_scaledown/* -1d  >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/tbarWInclusiveScaleDownFiles.txt
ls /scratch/data/tZqSkimsRun2017/tW_antitop_5f_scaleup/* -1d  >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/tbarWInclusiveScaleUpFiles.txt
ls /scratch/data/tZqSkimsRun2017/tW_top_5f/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/tWInclusiveFiles.txt
ls /scratch/data/tZqSkimsRun2017/tW_top_5f_scaledown/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/tWInclusiveScaleDownFiles.txt
ls /scratch/data/tZqSkimsRun2017/tW_top_5f_scaleup/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/tWInclusiveScaleUpFiles.txt

ls /scratch/data/tZqSkimsRun2017/tZq_ll_4Flavour3Lepton/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/tZqFiles.txt
ls /scratch/data/tZqSkimsRun2017/tZq_ll_4Flavour3Lepton_scaleup/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/tZqScaleUpFiles.txt
ls /scratch/data/tZqSkimsRun2017/tZq_ll_4Flavour3Lepton_scaledown/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/tZqScaleDownFiles.txt

ls /scratch/data/tZqSkimsRun2017/tHq/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/tHqFiles.txt
ls /scratch/data/tZqSkimsRun2017/wPlusJets/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/wPlusJetsFiles.txt
ls /scratch/data/tZqSkimsRun2017/WW1l1nu2q/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/WW1l1nu2qFiles.txt
ls /scratch/data/tZqSkimsRun2017/WW2l2nu/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/WW2l2nuFiles.txt
ls /scratch/data/tZqSkimsRun2017/WZ1l1nu2q/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/WZ1l1nu2q.txt
ls /scratch/data/tZqSkimsRun2017/WZ2l2q/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/WZ2l2q.txt
ls /scratch/data/tZqSkimsRun2017/WZJets/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/WZjets.txt
ls /scratch/data/tZqSkimsRun2017/ZZ2l2nu/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ZZ2l2nuFiles.txt
ls /scratch/data/tZqSkimsRun2017/ZZ2l2q/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ZZ2l2qFiles.txt
ls /scratch/data/tZqSkimsRun2017/ZZ4l/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ZZ4lFiles.txt

ls /scratch/data/tZqSkimsRun2017/WWW/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/WWWFiles.txt
ls /scratch/data/tZqSkimsRun2017/WWZ/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/WWZFiles.txt
ls /scratch/data/tZqSkimsRun2017/WZZ/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/WZZFiles.txt
ls /scratch/data/tZqSkimsRun2017/ZZZ/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ZZZFiles.txt

# QCD samples

ls /scratch/data/tZqSkimsRun2017/QCD_EMEnriched/Pt-20to30/* -1d >> configs/2017/datasets/fileLists/QCD_EMEnriched_Pt-20to30Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_EMEnriched/Pt-30to50/* -1d >> configs/2017/datasets/fileLists/QCD_EMEnriched_Pt-30to50Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_EMEnriched/Pt-30to50_ext/* -1d >> configs/2017/datasets/fileLists/QCD_EMEnriched_Pt-30to50Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_EMEnriched/Pt-50to80/* -1d >> configs/2017/datasets/fileLists/QCD_EMEnriched_Pt-50to80Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_EMEnriched/Pt-50to80_ext/* -1d >> configs/2017/datasets/fileLists/QCD_EMEnriched_Pt-50to80Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_EMEnriched/Pt-80to120/* -1d >> configs/2017/datasets/fileLists/QCD_EMEnriched_Pt-80to120Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_EMEnriched/Pt-80to120_ext/* -1d >> configs/2017/datasets/fileLists/QCD_EMEnriched_Pt-80to120Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_EMEnriched/Pt-120to170/* -1d >> configs/2017/datasets/fileLists/QCD_EMEnriched_Pt-120to170Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_EMEnriched/Pt-120to170_ext/* -1d >> configs/2017/datasets/fileLists/QCD_EMEnriched_Pt-120to170Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_EMEnriched/Pt-170to300/* -1d >> configs/2017/datasets/fileLists/QCD_EMEnriched_Pt-170to300Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_EMEnriched/Pt-300toInf/* -1d >> configs/2017/datasets/fileLists/QCD_EMEnriched_Pt-300toInfFiles.txt

ls /scratch/data/tZqSkimsRun2017/QCD_MuEnrichedPt5/Pt-15to20/* -1d >> configs/2017/datasets/fileLists/QCD_MuEnriched_Pt-15to20Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_MuEnrichedPt5/Pt-20to30/* -1d >> configs/2017/datasets/fileLists/QCD_MuEnriched_Pt-20to30Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_MuEnrichedPt5/Pt-30to50/* -1d >> configs/2017/datasets/fileLists/QCD_MuEnriched_Pt-30to50Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_MuEnrichedPt5/Pt-50to80/* -1d >> configs/2017/datasets/fileLists/QCD_MuEnriched_Pt-50to80Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_MuEnrichedPt5/Pt-80to120/* -1d >> configs/2017/datasets/fileLists/QCD_MuEnriched_Pt-80to120Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_MuEnrichedPt5/Pt-80to120_ext/* -1d >> configs/2017/datasets/fileLists/QCD_MuEnriched_Pt-80to120Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_MuEnrichedPt5/Pt-120to170/* -1d >> configs/2017/datasets/fileLists/QCD_MuEnriched_Pt-120to170Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_MuEnrichedPt5/Pt-170to300/* -1d >> configs/2017/datasets/fileLists/QCD_MuEnriched_Pt-170to300Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_MuEnrichedPt5/Pt-170to300_ext/* -1d >> configs/2017/datasets/fileLists/QCD_MuEnriched_Pt-170to300Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_MuEnrichedPt5/Pt-300to470/* -1d >> configs/2017/datasets/fileLists/QCD_MuEnriched_Pt-300to470Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_MuEnrichedPt5/Pt-300to470_ext1/* -1d >> configs/2017/datasets/fileLists/QCD_MuEnriched_Pt-300to470Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_MuEnrichedPt5/Pt-300to470_ext2/* -1d >> configs/2017/datasets/fileLists/QCD_MuEnriched_Pt-300to470Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_MuEnrichedPt5/Pt-470to600/* -1d >> configs/2017/datasets/fileLists/QCD_MuEnriched_Pt-470to600Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_MuEnrichedPt5/Pt-470to600_ext1/* -1d >> configs/2017/datasets/fileLists/QCD_MuEnriched_Pt-470to600Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_MuEnrichedPt5/Pt-470to600_ext2/* -1d >> configs/2017/datasets/fileLists/QCD_MuEnriched_Pt-470to600Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_MuEnrichedPt5/Pt-600to800/* -1d >> configs/2017/datasets/fileLists/QCD_MuEnriched_Pt-600to800Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_MuEnrichedPt5/Pt-600to800_ext/* -1d >> configs/2017/datasets/fileLists/QCD_MuEnriched_Pt-600to800Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_MuEnrichedPt5/Pt-800to1000/* -1d >> configs/2017/datasets/fileLists/QCD_MuEnriched_Pt-800to1000Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_MuEnrichedPt5/Pt-800to1000_ext1/* -1d >> configs/2017/datasets/fileLists/QCD_MuEnriched_Pt-800to1000Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_MuEnrichedPt5/Pt-800to1000_ext2/* -1d >> configs/2017/datasets/fileLists/QCD_MuEnriched_Pt-800to1000Files.txt
ls /scratch/data/tZqSkimsRun2017/QCD_MuEnrichedPt5/Pt-1000toInf/* -1d >> configs/2017/datasets/fileLists/QCD_MuEnriched_Pt-1000toInfFiles.txt
ls /scratch/data/tZqSkimsRun2017/QCD_MuEnrichedPt5/Pt-1000toInf_ext/* -1d >> configs/2017/datasets/fileLists/QCD_MuEnriched_Pt-1000toInfFiles.txt

# Synchronisation files
ls /scratch/data/tZqSkimsRun2017/synch/tZq/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/tZqSynchFiles.txt
ls /scratch/data/tZqSkimsRun2017/synch/data1/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/data1SynchFiles.txt
ls /scratch/data/tZqSkimsRun2017/synch/data2/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/data2SynchFiles.txt

ls /scratch/data/tZqSkimsRun2017/synch/tW/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/twSynchFiles.txt
ls /scratch/data/tZqSkimsRun2017/synch/ttbar/* -1d >> $TQZ_TOOLS_PATH/configs/2017/datasets/fileLists/ttbarSynchFiles.txt

echo Done!
echo Filelists have been updated.
