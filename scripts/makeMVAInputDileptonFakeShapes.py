#A tool to pull a load of information from the mva trees. Should be all easy jazz...

from ROOT import *

import math
import random
import sys
import os
from array import array
from jetCorrectionUncertainty import JetCorrectionUncertainty

def deltaR(eta1, phi1, eta2, phi2):
    ###Returns the delta R from the eta and phi of two particles.
    dEta = eta1-eta2
    dPhi = phi1-phi2
    while ( abs(dPhi) > math.pi ):
	if (dPhi > 0.0):
	    dPhi += -2*math.pi
	else:
	    dPhi += 2*math.pi
    return math.sqrt( dEta*dEta + dPhi*dPhi )

def sortOutLeptons(tree,channel):
    ###Returns two LorentzVectors containing the two z leptons. This will be VERY useful for making all of the plots.
    #Reads the position of the z leptons from variables stored at mvaTree making time, because I'm great and finally got around to doing it.
    zLep1,zLep2 = 0,0
    #Let's try commenting this out and see if everything breaks? Hopefully it won't do...
    #if tree.numElePF2PAT < 3:
    if channel == "ee":
        zLep1 = TLorentzVector(tree.elePF2PATGsfPx[tree.zLep1Index],tree.elePF2PATGsfPy[tree.zLep1Index],tree.elePF2PATGsfPz[tree.zLep1Index],tree.elePF2PATGsfE[tree.zLep1Index])
        zLep2 = TLorentzVector(tree.elePF2PATGsfPx[tree.zLep2Index],tree.elePF2PATGsfPy[tree.zLep2Index],tree.elePF2PATGsfPz[tree.zLep2Index],tree.elePF2PATGsfE[tree.zLep2Index])
    if channel == "mumu":
        zLep1 = TLorentzVector(tree.muonPF2PATPx[tree.zLep1Index],tree.muonPF2PATPy[tree.zLep1Index],tree.muonPF2PATPz[tree.zLep1Index],tree.muonPF2PATE[tree.zLep1Index])
        zLep2 = TLorentzVector(tree.muonPF2PATPx[tree.zLep2Index],tree.muonPF2PATPy[tree.zLep2Index],tree.muonPF2PATPz[tree.zLep2Index],tree.muonPF2PATE[tree.zLep2Index])
    return (zLep1,zLep2)

def sortOutHadronicW(tree,channel):

    ###Returns two LorentzVectors containing the two w quarks. This will be VERY useful for making all of the plots.
    #Reads the position of the w quarks from variables stored at mvaTree making time, because I'm great and finally got around to doing it.
    wQuark1,wQuark2 = 0,0

    wQuark1 = TLorentzVector(tree.jetPF2PATPx[tree.wQuark1Index],tree.jetPF2PATPy[tree.wQuark1Index],tree.jetPF2PATPz[tree.wQuark1Index ],tree.jetPF2PATE[tree.wQuark1Index])
    wQuark2 = TLorentzVector(tree.jetPF2PATPx[tree.wQuark2Index],tree.jetPF2PATPy[tree.wQuark2Index],tree.jetPF2PATPz[tree.wQuark2Index ],tree.jetPF2PATE[tree.wQuark2Index])
    return (wQuark1,wQuark2)

def getJets(tree,syst,jetUnc,met,is2016):
    #Makes a short list of indices of the jets in the event
    jetList = []
    jetVecList = []
    for i in range(15):
        if tree.jetInd[i] > -.5:
            jetList.append(tree.jetInd[i])
            jetVecList.append(getJetVec(tree,tree.jetInd[i],met,is2016,syst,True))
        else: continue
    return (jetList,jetVecList)

def getBjets(tree,syst,jetUnc,met,jets,is2016):
    #Return a list of the indices of the b-jets in the event
    bJetList = []
    bJetVecList = []
    for i in range(10):
        if tree.bJetInd[i] > -0.5:
            bJetList.append(tree.bJetInd[i])
            bJetVecList.append(getJetVec(tree,jets[tree.bJetInd[i]],met,is2016,syst,False))
        else:continue
#    print len(bJetList)
    return (bJetList,bJetVecList)

def getJetVec(tree, index, metVec, is2016, syst, doMetSmear):
    #Gets a vector for a jet with corrections already applied


    newSmearValue = tree.jetSmearValue[index];
#    newSmearValue = 1.0;
    returnJet = TLorentzVector();

    returnJet.SetPxPyPzE(newSmearValue*tree.jetPF2PATPx[index],newSmearValue*tree.jetPF2PATPy[index],newSmearValue*tree.jetPF2PATPz[index],newSmearValue*tree.jetPF2PATE[index]);

    if syst == 16:
        returnJet *= 1+ jetUnc.getUncertainty(returnJet.Pt(), returnJet.Eta(),1)
    elif syst == 32:
        returnJet *= 1+ jetUnc.getUncertainty(returnJet.Pt(), returnJet.Eta(),2)

    if ( doMetSmear and newSmearValue > 0.01 ) :
    #Propogate through the met. But only do it if the smear jet isn't 0.
        metVec.SetPx(metVec.Px()+tree.jetPF2PATPx[index])
        metVec.SetPy(metVec.Py()+tree.jetPF2PATPy[index])

        metVec.SetPx(metVec.Px()-returnJet.Px())
        metVec.SetPy(metVec.Py()-returnJet.Py())

    return returnJet


def doUncMet(tree,met,zLep1,zLep2,jetVecs,syst):

    uncMetX = met.Px() + zLep1.Px() + zLep2.Px()
    uncMetY = met.Py() + zLep1.Py() + zLep2.Py()

    for i in range(len(jetVecs)):
        uncMetX += jetVecs[i].Px()
        uncMetY += jetVecs[i].Py()

    if syst == 1024:
        met.SetPx(met.Px() + 0.1*uncMetX)
        met.SetPy(met.Py() + 0.1*uncMetY)
    else:
        met.SetPx(met.Px() - 0.1*uncMetX)
        met.SetPy(met.Py() - 0.1*uncMetY)

    return met

def setupInputVars():
    #Make the variables we want to save
    inputVars = {}
    inputVars["eventWeight"] = array('f',[0.])
    inputVars["eventNumber"] = array('f',[0.])
    inputVars["mTW"] = array('f',[0.])
    inputVars["wQuark1Pt"] = array('f',[0.])
    inputVars["wQuark1Eta"] = array('f',[0.])
    inputVars["wQuark1Phi"] = array('f',[0.])
    inputVars["wQuark2Pt"] = array('f',[0.])
    inputVars["wQuark2Eta"] = array('f',[0.])
    inputVars["wQuark2Phi"] = array('f',[0.])
    inputVars["wPairMass"] = array('f',[0.])
    inputVars["wPairPt"] = array('f',[0.])
    inputVars["wPairEta"] = array('f',[0.])
    inputVars["wPairPhi"] = array('f',[0.])
    inputVars["met"] = array('f',[0.])
    inputVars["nJets"] = array('f',[0.])
    inputVars["leadJetPt"] = array('f',[0.])
    inputVars["leadJetPhi"] = array('f',[0.])
    inputVars["leadJetEta"] = array('f',[0.])
    inputVars["leadJetbTag"] = array('f',[0.])
    inputVars["secJetPt"] = array('f',[0.])
    inputVars["secJetPhi"] = array('f',[0.])
    inputVars["secJetEta"] = array('f',[0.])
    inputVars["secJetbTag"] = array('f',[0.])
    inputVars["thirdJetPt"] = array('f',[0.])
    inputVars["thirdJetPhi"] = array('f',[0.])
    inputVars["thirdJetEta"] = array('f',[0.])
    inputVars["thirdJetbTag"] = array('f',[0.])
    inputVars["fourthJetPt"] = array('f',[0.])
    inputVars["fourthJetPhi"] = array('f',[0.])
    inputVars["fourthJetEta"] = array('f',[0.])
    inputVars["fourthJetbTag"] = array('f',[0.])
    inputVars["nBjets"] = array('f',[0.])
    inputVars["bTagDisc"] = array('f',[0.])
    inputVars["lep1Pt"] = array('f',[0.])
    inputVars["lep1Eta"] = array('f',[0.])
    inputVars["lep1Phi"] = array('f',[0.])
    inputVars["lep1RelIso"] = array('f',[0.])
    inputVars["lep1D0"] = array('f',[0.])
    inputVars["lep2Pt"] = array('f',[0.])
    inputVars["lep2Eta"] = array('f',[0.])
    inputVars["lep2Phi"] = array('f',[0.])
    inputVars["lep2RelIso"] = array('f',[0.])
    inputVars["lep2D0"] = array('f',[0.])
    inputVars["lepMass"] = array('f',[0.])
    inputVars["lepPt"] = array('f',[0.])
    inputVars["lepEta"] = array('f',[0.])
    inputVars["lepPhi"] = array('f',[0.])
    inputVars["zMass"] = array('f',[0.])
    inputVars["zPt"] = array('f',[0.])
    inputVars["zEta"] = array('f',[0.])
    inputVars["zPhi"] = array('f',[0.])
    inputVars["topMass"] = array('f',[0.])
    inputVars["topPt"] = array('f',[0.])
    inputVars["topEta"] = array('f',[0.])
    inputVars["topPhi"] = array('f',[0.])
    inputVars["j1j2delR"] = array('f',[0.])
    inputVars["j1j2delPhi"] = array('f',[0.])
    inputVars["w1w2delR"] = array('f',[0.])
    inputVars["w1w2delPhi"] = array('f',[0.])
    inputVars["zLepdelR"] = array('f',[0.])
    inputVars["zLepdelPhi"] = array('f',[0.])
    inputVars["zl1Quark1DelR"] = array('f',[0.])
    inputVars["zl1Quark1DelPhi"] = array('f',[0.])
    inputVars["zl1Quark2DelR"] = array('f',[0.])
    inputVars["zl1Quark2DelPhi"] = array('f',[0.])
    inputVars["zl2Quark1DelR"] = array('f',[0.])
    inputVars["zl2Quark1DelPhi"] = array('f',[0.])
    inputVars["zl2Quark2DelR"] = array('f',[0.])
    inputVars["zl2Quark2DelPhi"] = array('f',[0.])
    inputVars["zlb1DelR"] = array('f',[0.])
    inputVars["zlb1DelPhi"] = array('f',[0.])
    inputVars["zlb2DelR"] = array('f',[0.])
    inputVars["zlb2DelPhi"] = array('f',[0.])
    inputVars["lepHt"] = array('f',[0.])
    inputVars["wQuarkHt"] = array('f',[0.])
    inputVars["totPt"] = array('f',[0.])
    inputVars["totEta"] = array('f',[0.])
    inputVars["totPtVec"] = array('f',[0.])
    inputVars["totVecM"] = array('f',[0.])
    inputVars["chan"] = array('i',[0])
    inputVars["totPt2Jet"] = array('f',[0.])
    inputVars["wZdelR"] = array('f',[0.])
    inputVars["wZdelPhi"] = array('f',[0.])
    inputVars["zQuark1DelR"] = array('f',[0.])
    inputVars["zQuark1DelPhi"] = array('f',[0.])
    inputVars["zQuark2DelR"] = array('f',[0.])
    inputVars["zQuark2DelPhi"] = array('f',[0.])
    inputVars["zTopDelR"] = array('f',[0.])
    inputVars["zTopDelPhi"] = array('f',[0.])
    inputVars["zl1TopDelR"] = array('f',[0.])
    inputVars["zl1TopDelPhi"] = array('f',[0.])
    inputVars["zl2TopDelR"] = array('f',[0.])
    inputVars["zl2TopDelPhi"] = array('f',[0.])
    inputVars["wTopDelR"] = array('f',[0.])
    inputVars["wTopDelPhi"] = array('f',[0.])
    inputVars["w1TopDelR"] = array('f',[0.])
    inputVars["w1TopDelPhi"] = array('f',[0.])
    inputVars["w2TopDelR"] = array('f',[0.])
    inputVars["w2TopDelPhi"] = array('f',[0.])
    inputVars["minZJetR"] = array('f',[0.])
    inputVars["minZJetPhi"] = array('f',[0.])
    inputVars["totHt"] = array('f',[0.])
    inputVars["jetHt"] = array('f',[0.])
    inputVars["jetMass"] = array('f',[0.])
    inputVars["jetMass3"] = array('f',[0.])
    inputVars["totHtOverPt"] = array('f',[0.])
    inputVars["chi2"] = array('f',[0.])
    return inputVars

def setupBranches(tree,varMap):
    tree.Branch("EvtWeight", varMap["eventWeight"], "EvtWeight/F")
    tree.Branch("EvtNumber", varMap["eventNumber"], "EvtNumber/F")
    tree.Branch("mTW", varMap["mTW"], "mTW/F")
    tree.Branch("wQuark1Pt", varMap["wQuark1Pt"], "wQuark1Pt/F")
    tree.Branch("wQuark1Eta", varMap["wQuark1Eta"], "wQuark1Eta/F")
    tree.Branch("wQuark1Phi", varMap["wQuark1Phi"], "wQuark1Phi/F")
    tree.Branch("wQuark2Pt", varMap["wQuark2Pt"], "wQuark2Pt/F")
    tree.Branch("wQuark2Eta", varMap["wQuark2Eta"], "wQuark2Eta/F")
    tree.Branch("wQuark2Phi", varMap["wQuark2Phi"], "wQuark2Phi/F")
    tree.Branch("wPairMass", varMap["wPairMass"], "wPairMass/F")
    tree.Branch("wPairPt", varMap["wPairPt"], "wPairPt/F")
    tree.Branch("wPairEta", varMap["wPairEta"], "wPairEta/F")
    tree.Branch("wPairPhi", varMap["wPairPhi"], "wPairPhi/F")
    tree.Branch("met",varMap["met"],"met/F")
    tree.Branch("nJets",varMap["nJets"],"nJets/F")
    tree.Branch("leadJetPt",varMap["leadJetPt"],"leadJetPt/F")
    tree.Branch("leadJetEta",varMap["leadJetEta"],"leadJetEta/F")
    tree.Branch("leadJetPhi",varMap["leadJetPhi"],"leadJetPhi/F")
    tree.Branch("leadJetbTag",varMap["leadJetbTag"],"leadJetbTag/F")
    tree.Branch("secJetPt",varMap["secJetPt"],"secJetPt/F")
    tree.Branch("secJetEta",varMap["secJetEta"],"secJetEta/F")
    tree.Branch("secJetPhi",varMap["secJetPhi"],"secJetPhi/F")
    tree.Branch("secJetbTag",varMap["secJetbTag"],"secJetbTag/F")
    tree.Branch("thirdJetPt",varMap["thirdJetPt"],"thirdJetPt/F")
    tree.Branch("thirdJetEta",varMap["thirdJetEta"],"thirdJetEta/F")
    tree.Branch("thirdJetPhi",varMap["thirdJetPhi"],"thirdJetPhi/F")
    tree.Branch("thirdJetbTag",varMap["thirdJetbTag"],"thirdJetbTag/F")
    tree.Branch("fourthJetPt",varMap["fourthJetPt"],"fourthJetPt/F")
    tree.Branch("fourthJetEta",varMap["fourthJetEta"],"fourthJetEta/F")
    tree.Branch("fourthJetPhi",varMap["fourthJetPhi"],"fourthJetPhi/F")
    tree.Branch("fourthJetbTag",varMap["fourthJetbTag"],"fourthJetbTag/F")
    tree.Branch("nBjets",varMap["nBjets"],"nBjets/F")
    tree.Branch("bTagDisc",varMap["bTagDisc"],"bTagDisc/F")
    tree.Branch("lep1Pt",varMap["lep1Pt"],"lep1Pt/F")
    tree.Branch("lep1Eta",varMap["lep1Eta"],"lep1Eta/F")
    tree.Branch("lep1Phi",varMap["lep1Phi"],"lep1Phi/F")
    tree.Branch("lep1RelIso",varMap["lep1RelIso"],"lep1RelIso/F")
    tree.Branch("lep1D0",varMap["lep1D0"],"lep1D0/F")
    tree.Branch("lep2Pt",varMap["lep2Pt"],"lep2Pt/F")
    tree.Branch("lep2Eta",varMap["lep2Eta"],"lep2Eta/F")
    tree.Branch("lep2Phi",varMap["lep2Phi"],"lep2Phi/F")
    tree.Branch("lep2RelIso",varMap["lep2RelIso"],"lep2RelIso/F")
    tree.Branch("lep2D0",varMap["lep2D0"],"lep2D0/F")
    tree.Branch("lepMass",varMap["lepMass"],"lepMass/F")
    tree.Branch("lepPt",varMap["lepPt"],"lepPt/F")
    tree.Branch("lepEta",varMap["lepEta"],"lepEta/F")
    tree.Branch("lepPhi",varMap["lepPhi"],"lepPhi/F")
    tree.Branch("zMass",varMap["zMass"],"zMass/F")
    tree.Branch("zPt",varMap["zPt"],"zPt/F")
    tree.Branch("zEta",varMap["zEta"],"zEta/F")
    tree.Branch("zPhi",varMap["zPhi"],"zPhi/F")
    tree.Branch("topMass",varMap["topMass"],"topMass/F")
    tree.Branch("topPt",varMap["topPt"],"topPt/F")
    tree.Branch("topEta",varMap["topEta"],"topEta/F")
    tree.Branch("topPhi",varMap["topPhi"],"topPhi/F")
    tree.Branch("jjdelR",varMap["j1j2delR"],"jjdelR/F")
    tree.Branch("jjdelPhi",varMap["j1j2delPhi"],"jjdelPhi/F")
    tree.Branch("wwdelR",varMap["w1w2delR"],"wwdelR/F")
    tree.Branch("wwdelPhi",varMap["w1w2delPhi"],"wwdelPhi/F")
    tree.Branch("zLepdelR",varMap["zLepdelR"],"zLepdelR/F")
    tree.Branch("zLepdelPhi",varMap["zLepdelPhi"],"zLepdelPhi/F")
    tree.Branch("zl1Quark1DelR",varMap["zl1Quark1DelR"],"zl1Quark1DelR/F")
    tree.Branch("zl1Quark1DelPhi",varMap["zl1Quark1DelPhi"],"zl1Quark1DelPhi/F")
    tree.Branch("zl1Quark2DelR",varMap["zl1Quark2DelR"],"zl1Quark2DelR/F")
    tree.Branch("zl1Quark2DelPhi",varMap["zl1Quark2DelPhi"],"zl1Quark2DelPhi/F")
    tree.Branch("zl2Quark1DelR",varMap["zl2Quark1DelR"],"zl2Quark1DelR/F")
    tree.Branch("zl2Quark1DelPhi",varMap["zl2Quark1DelPhi"],"zl2Quark1DelPhi/F")
    tree.Branch("zl2Quark2DelR",varMap["zl2Quark2DelR"],"zl2Quark2DelR/F")
    tree.Branch("zl2Quark2DelPhi",varMap["zl2Quark2DelPhi"],"zl2Quark2DelPhi/F")
    tree.Branch("zlb1DelR",varMap["zlb1DelR"],"zlb1DelR/F")
    tree.Branch("zlb1DelPhi",varMap["zlb1DelPhi"],"zlb1DelPhi/F")
    tree.Branch("zlb2DelR",varMap["zlb2DelR"],"zlb2DelR/F")
    tree.Branch("zlb2DelPhi",varMap["zlb2DelPhi"],"zlb2DelPhi/F")
    tree.Branch("lepHt",varMap["lepHt"],"lepHt/F")
    tree.Branch("wQuarkHt",varMap["wQuarkHt"],"wQuarkHt/F")
    tree.Branch("totPt",varMap["totPt"],"totPt/F")
    tree.Branch("totEta",varMap["totEta"],"totEta/F")
    tree.Branch("totPtVec",varMap["totPtVec"],"totPtVec/F")
    tree.Branch("totVecM",varMap["totVecM"],"totVecM/F")
    tree.Branch("Channel",varMap["chan"],"Channel/I")
    tree.Branch("totPt2Jet",varMap["totPt2Jet"],"totPt2Jet/F")
    tree.Branch("wzdelR",varMap["wZdelR"],"wzdelR/F")
    tree.Branch("wzdelPhi",varMap["wZdelPhi"],"wzdelPhi/F")
    tree.Branch("zQuark1DelR",varMap["zQuark1DelR"],"zQuark1DelR/F")
    tree.Branch("zQuark1DelPhi",varMap["zQuark1DelPhi"],"zQuark1DelPhi/F")
    tree.Branch("zQuark2DelR",varMap["zQuark2DelR"],"zQuark2DelR/F")
    tree.Branch("zQuark2DelPhi",varMap["zQuark2DelPhi"],"zQuark2DelPhi/F")
    tree.Branch("zTopDelR",varMap["zTopDelR"],"zTopDelR/F")
    tree.Branch("zTopDelPhi",varMap["zTopDelPhi"],"zTopDelPhi/F")
    tree.Branch("zl1TopDelR",varMap["zl1TopDelR"],"zl1TopDelR/F")
    tree.Branch("zl1TopDelPhi",varMap["zl1TopDelPhi"],"zl1TopDelPhi/F")
    tree.Branch("zl2TopDelR",varMap["zl2TopDelR"],"zl2TopDelR/F")
    tree.Branch("zl2TopDelPhi",varMap["zl2TopDelPhi"],"zl2TopDelPhi/F")
    tree.Branch("wTopDelR",varMap["wTopDelR"],"wTopDelR/F")
    tree.Branch("wTopDelPhi",varMap["wTopDelPhi"],"wTopDelPhi/F")
    tree.Branch("w1TopDelR",varMap["w1TopDelR"],"w1TopDelR/F")
    tree.Branch("w1TopDelPhi",varMap["w1TopDelPhi"],"w1TopDelPhi/F")
    tree.Branch("w2TopDelR",varMap["w2TopDelR"],"w2TopDelR/F")
    tree.Branch("w2TopDelPhi",varMap["w2TopDelPhi"],"w2TopDelPhi/F")
    tree.Branch("zjminR",varMap["minZJetR"],"zjminR/F")
    tree.Branch("zjminPhi",varMap["minZJetPhi"],"zjminPhi/F")
    tree.Branch("totHt",varMap["totHt"],"totHt/F")
    tree.Branch("jetHt",varMap["jetHt"],"jetHt/F")
    tree.Branch("jetMass",varMap["jetMass"],"jetMass/F")
    tree.Branch("jetMass3",varMap["jetMass3"],"jetMass3/F")
    tree.Branch("totHtOverPt",varMap["totHtOverPt"],"totHtOverPt/F")
    tree.Branch("chi2",varMap["chi2"],"chi2/F")

def fillTree(outTreeSig, outTreeSdBnd, varMap, tree, label, jetUnc, channel, is2016, SameSignMC = False):
    #Fills the output tree. This is a new function because I want to access data and MC in different ways but do the same thing to them in the end.

    syst = 0

    if "__met__plus" in label:
        syst = 1024
    if "__met__minus" in label:
        syst = 2048

    if channel == "ee":
        varMap["chan"][0] = 1
    if channel == "mumu":
        varMap["chan"][0] = 0

    #Set up those variables as branches
    for event in range(tree.GetEntries()):
            #Fill some plots here. Let's make an example mTW plot.
            #Make a config that'll do this for me? I've done these before so should be easy. Fill expressions could be a pain?
        tree.GetEntry(event)

        ##Save event number for debugging
        varMap["eventNumber"][0] = tree.eventNum

        ##Now the real stuff!
        (zLep1,zLep2) = sortOutLeptons(tree,channel)
        metVec = TLorentzVector(tree.metPF2PATPx,tree.metPF2PATPy,0,tree.metPF2PATEt)
        (jets,jetVecs) = getJets(tree,syst,jetUnc,metVec,is2016)
        (bJets,bJetVecs) = getBjets(tree,syst,jetUnc,metVec,jets,is2016)
        (wQuark1,wQuark2) = sortOutHadronicW(tree,channel)
        #Do unclustered met stuff here now that we have all of the objects, all corrected for their various SFs etc.
        if syst == 1024 or syst == 2048:
            metVec = doUncMet(tree,metVec,zLep1,zLep2,jetVecs,syst)
        if ( SameSignMC == True and channel == "ee" ) : varMap["eventWeight"][0] = tree.eventWeight * 1.5530 # SF we weight fake ee shape by
        elif ( SameSignMC == True and channel == "mumu" ) : varMap["eventWeight"][0] = tree.eventWeight * 1.2248 # SF we weight fake ee shape by
        else : varMap["eventWeight"][0] = tree.eventWeight
        varMap["leadJetPt"][0] = jetVecs[0].Pt()
        varMap["leadJetEta"][0] = jetVecs[0].Eta()
        varMap["leadJetPhi"][0] = jetVecs[0].Phi()
        #Make all the random Pt variables I'm saving for some reason
        totPx,totPy = 0,0
        totPx += zLep1.Px() + zLep2.Px()
        totPy += zLep1.Py() + zLep2.Py()
	varMap["lep1Pt"][0] = zLep1.Pt()
	varMap["lep1Eta"][0] = zLep1.Eta()
	varMap["lep1Phi"][0] = zLep1.Phi()
        if channel == "ee":
            varMap["lep1RelIso"][0] = tree.elePF2PATComRelIsoRho[tree.zLep1Index]
	    varMap["lep1D0"][0] = tree.elePF2PATD0PV[tree.zLep1Index]
	    varMap["lep2RelIso"][0] = tree.elePF2PATComRelIsoRho[tree.zLep2Index]
	    varMap["lep2D0"][0] = tree.elePF2PATD0PV[tree.zLep2Index]
	if channel == "mumu":
            varMap["lep1RelIso"][0] = tree.muonPF2PATComRelIsodBeta[tree.zLep1Index]
	    varMap["lep1D0"][0] = tree.muonPF2PATDBPV[tree.zLep1Index]
	    varMap["lep2RelIso"][0] = tree.muonPF2PATComRelIsodBeta[tree.zLep2Index]
	    varMap["lep2D0"][0] = tree.muonPF2PATDBPV[tree.zLep2Index]

	varMap["lep2Pt"][0] = zLep2.Pt()
	varMap["lep2Eta"][0] = zLep2.Eta()
	varMap["lep2Phi"][0] = zLep2.Phi()
	varMap["lepMass"][0] = ( zLep1 + zLep2 ).M()
        varMap["lepPt"][0] = math.sqrt(totPx * totPx + totPy * totPy)
	varMap["lepEta"][0] = ( zLep1 + zLep2 ).Eta()
	varMap["lepPhi"][0] = ( zLep1 + zLep2 ).Phi()
	varMap["wQuark1Pt"][0] = wQuark1.Pt()
	varMap["wQuark1Eta"][0] = wQuark1.Eta()
	varMap["wQuark1Phi"][0] = wQuark1.Phi()
	varMap["wQuark2Pt"][0] = wQuark2.Pt()
	varMap["wQuark2Eta"][0] = wQuark2.Eta()
	varMap["wQuark2Phi"][0] = wQuark2.Phi()
	wPairMass = ( wQuark1 + wQuark2 ).M()
	varMap["wPairMass"][0] = wPairMass
	varMap["wPairPt"][0] = ( wQuark1 + wQuark2 ).Pt()
	varMap["wPairEta"][0] = ( wQuark1 + wQuark2 ).Eta()
	varMap["wPairPhi"][0] = ( wQuark1 + wQuark2 ).Phi()
        totPx += jetVecs[0].Px()
        totPy += jetVecs[0].Py()
        if len(jetVecs) > 1:
            totPx += jetVecs[1].Px()
            totPy += jetVecs[1].Py()
        varMap["totPt2Jet"][0] = math.sqrt(totPx * totPx + totPy * totPy)
        for i in range(2,len(jets)):
            totPx+=jetVecs[i].Px()
            totPy+=jetVecs[i].Py()
        varMap["totPt"][0] = math.sqrt(totPx * totPx + totPy * totPy)
        totVec = (zLep1+zLep2)
        for i in range(len(jetVecs)):
            totVec += jetVecs[i]
        varMap["totEta"][0] = totVec.Eta()
        varMap["totPtVec"][0] = totVec.Pt()
        varMap["totVecM"][0] = totVec.M()
        varMap["mTW"][0] = math.sqrt(2*tree.jetPF2PATPt[tree.wQuark1Index]*tree.jetPF2PATPt[tree.wQuark2Index] * (1-math.cos(tree.jetPF2PATPhi[tree.wQuark1Index] - tree.jetPF2PATPhi[tree.wQuark2Index])))
        varMap["nJets"][0] = float(len(jets))
        varMap["nBjets"][0] = float(len(bJets))
        varMap["met"][0] = tree.metPF2PATEt
        varMap["bTagDisc"][0] = tree.jetPF2PATBDiscriminator[jets[bJets[0]]]
        varMap["leadJetbTag"][0] = tree.jetPF2PATBDiscriminator[jets[0]]
        varMap["secJetbTag"][0] = -1.
        varMap["secJetPt"][0] = -1.
        varMap["secJetEta"][0] = -500.
        varMap["secJetPhi"][0] = -500.
        varMap["thirdJetbTag"][0] = -1.
        varMap["thirdJetPt"][0] = -1.
        varMap["thirdJetEta"][0] = -500.
        varMap["thirdJetPhi"][0] = -500.
        varMap["fourthJetbTag"][0] = -1.
        varMap["fourthJetPt"][0] = -1.
        varMap["fourthJetEta"][0] = -500.
        varMap["fourthJetPhi"][0] = -500.
        if len(jetVecs) > 1:
            varMap["secJetPt"][0] = jetVecs[1].Pt()
            varMap["secJetEta"][0] = jetVecs[1].Eta()
            varMap["secJetPhi"][0] = jetVecs[1].Phi()
            varMap["secJetbTag"][0] = tree.jetPF2PATBDiscriminator[jets[1]]
        if len(jetVecs) > 2:
            varMap["thirdJetPt"][0] = jetVecs[2].Pt()
            varMap["thirdJetEta"][0] = jetVecs[2].Eta()
            varMap["thirdJetPhi"][0] = jetVecs[2].Phi()
            varMap["thirdJetbTag"][0] = tree.jetPF2PATBDiscriminator[jets[2]]
        if len(jetVecs) > 3:
            varMap["fourthJetPt"][0] = jetVecs[3].Pt()
            varMap["fourthJetEta"][0] = jetVecs[3].Eta()
            varMap["fourthJetPhi"][0] = jetVecs[3].Phi()
            varMap["fourthJetbTag"][0] = tree.jetPF2PATBDiscriminator[jets[3]]

#        print bTagDisc[0], bJets[0], tree.jetPF2PATBDiscriminator[jets[bJets[0]]], len(bJets), nBjets[0]
	topMass = (bJetVecs[0] + wQuark1 + wQuark2).M()
        varMap["topMass"][0] = topMass
        varMap["topPt"][0] = (bJetVecs[0] + wQuark1 + wQuark2).Pt()
        varMap["topEta"][0] = (bJetVecs[0] + wQuark1 + wQuark2).Eta()
        varMap["topPhi"][0] = (bJetVecs[0] + wQuark1 + wQuark2).Phi()
        varMap["wZdelR"][0] = (zLep2 + zLep1).DeltaR(wQuark1 + wQuark2)
        varMap["wZdelPhi"][0] = (zLep2 + zLep1).DeltaPhi(wQuark1 + wQuark2)

	varMap["zQuark1DelR"][0] = (zLep2 + zLep1).DeltaR(wQuark1)
	varMap["zQuark1DelPhi"][0] = (zLep2 + zLep1).DeltaPhi(wQuark1)
	varMap["zQuark2DelR"][0] = (zLep2 + zLep1).DeltaR(wQuark2)
	varMap["zQuark2DelPhi"][0] = (zLep2 + zLep1).DeltaPhi(wQuark2)

	varMap["zTopDelR"][0] = (zLep2 + zLep1).DeltaR(bJetVecs[0] + wQuark1 + wQuark2)
	varMap["zTopDelPhi"][0] = (zLep2 + zLep1).DeltaPhi(bJetVecs[0] + wQuark1 + wQuark2)
	varMap["zl1TopDelR"][0] = (zLep1).DeltaR(bJetVecs[0] + wQuark1 + wQuark2)
	varMap["zl1TopDelPhi"][0] = (zLep1).DeltaPhi(bJetVecs[0] + wQuark1 + wQuark2)
	varMap["zl2TopDelR"][0] = (zLep2).DeltaR(bJetVecs[0] + wQuark1 + wQuark2)
	varMap["zl2TopDelPhi"][0] = (zLep2).DeltaPhi(bJetVecs[0] + wQuark1 + wQuark2)

	varMap["wTopDelR"][0] = (wQuark1 + wQuark2).DeltaR(bJetVecs[0] + wQuark1 + wQuark2)
	varMap["wTopDelPhi"][0] = (wQuark1 + wQuark2).DeltaPhi(bJetVecs[0] + wQuark1 + wQuark2)
	varMap["w1TopDelR"][0] = (wQuark1).DeltaR(bJetVecs[0] + wQuark1 + wQuark2)
	varMap["w1TopDelR"][0] = (wQuark1).DeltaR(bJetVecs[0] + wQuark1 + wQuark2)
	varMap["w1TopDelPhi"][0] = (wQuark1).DeltaPhi(bJetVecs[0] + wQuark1 + wQuark2)
	varMap["w2TopDelR"][0] = (wQuark2).DeltaR(bJetVecs[0] + wQuark1 + wQuark2)
	varMap["w2TopDelPhi"][0] = (wQuark2).DeltaPhi(bJetVecs[0] + wQuark1 + wQuark2)

        varMap["j1j2delR"][0] = -1.
        varMap["j1j2delPhi"][0] = -10.
        if len(jetVecs) > 1:
            varMap["j1j2delR"][0] = jetVecs[0].DeltaR(jetVecs[1])
            varMap["j1j2delPhi"][0] = jetVecs[0].DeltaPhi(jetVecs[1])
	varMap["w1w2delR"][0] = (wQuark1).DeltaR(wQuark2)
	varMap["w1w2delPhi"][0] = (wQuark1).DeltaPhi(wQuark2)
	varMap["zLepdelR"][0] = (zLep1).DeltaR(zLep2)
	varMap["zLepdelPhi"][0] = (zLep1).DeltaPhi(zLep2)
	varMap["zl1Quark1DelR"][0] = (zLep1).DeltaR(wQuark1)
	varMap["zl1Quark1DelPhi"][0] = (zLep1).DeltaPhi(wQuark1)
	varMap["zl1Quark2DelR"][0] = (zLep1).DeltaR(wQuark2)
	varMap["zl1Quark2DelPhi"][0] = (zLep1).DeltaPhi(wQuark2)
	varMap["zl2Quark1DelR"][0] = (zLep2).DeltaR(wQuark1)
	varMap["zl2Quark1DelPhi"][0] = (zLep2).DeltaPhi(wQuark1)
	varMap["zl2Quark2DelR"][0] = (zLep2).DeltaR(wQuark2)
	varMap["zl2Quark2DelPhi"][0] = (zLep2).DeltaPhi(wQuark2)
        varMap["minZJetR"][0] = 3.0
        jetHt = 0.
        for i in range(len(jetVecs)):
            jetHt += jetVecs[i].Pt()
            jetVector = TLorentzVector();
            if jetVecs[i].DeltaR(zLep2 + zLep1) < varMap["minZJetR"][0]:
                varMap["minZJetR"][0] = jetVecs[i].DeltaR(zLep2 + zLep1)
            if jetVecs[i].DeltaPhi(zLep2 + zLep1) < varMap["minZJetR"][0]:
                varMap["minZJetPhi"][0] = jetVecs[i].DeltaPhi(zLep2 + zLep1)
        varMap["zlb1DelR"][0] = zLep1.DeltaR(jetVecs[bJets[0]])
        varMap["zlb1DelPhi"][0] = zLep1.DeltaPhi(jetVecs[bJets[0]])
        varMap["zlb2DelR"][0] = zLep2.DeltaR(jetVecs[bJets[0]])
        varMap["zlb2DelPhi"][0] = zLep2.DeltaPhi(jetVecs[bJets[0]])
        ht = 0.
        ht += zLep1.Pt() + zLep2.Pt()
        varMap["lepHt"][0] = ht
        varMap["jetHt"][0] = jetHt
        varMap["jetMass"][0] = jetVector.M()
        varMap["jetMass3"][0] = (jetVecs[0] + jetVecs[1] + jetVecs[2]).M()
	varMap["wQuarkHt"][0] = wQuark1.Pt()+wQuark2.Pt()
        ht += jetHt
        varMap["totHt"][0] = ht
        varMap["totHtOverPt"][0] = ht / math.sqrt(totPx * totPx + totPy * totPy)
        varMap["zMass"][0] = (zLep1+zLep2).M()
        varMap["zPt"][0] = (zLep2 + zLep1).Pt()
        varMap["zEta"][0] = (zLep2 + zLep1).Eta()
        varMap["zPhi"][0] = (zLep2 + zLep1).Phi()

	wChi2Term = (wPairMass - 80.3585)/8.0
	topChi2Term = (topMass - 173.21)/30.0
	varMap["chi2"][0] = wChi2Term*wChi2Term + topChi2Term*topChi2Term

	if outTreeSdBnd :
 	    if varMap["chi2"][0] >= 40. and varMap["chi2"][0] < 150. :
                outTreeSdBnd.Fill()
            if varMap["chi2"][0] < 40. :
                 outTreeSig.Fill()
        else :
            outTreeSig.Fill()

def main():

    #Mapping of our mc names to IPHC names
#    listOfMCs = {"ttHTobb" : "ttH", "ttHToNonbb" : "ttH", "WWW" : "WWW", "WWZ" : "WWZ", "WZZ" : "WZZ", "ZZZ" : "ZZZ", "WW1l1nu2q" : "WW", "WW2l2nu":"WW","ZZ4l":"ZZ","ZZ2l2nu":"ZZ","ZZ2l2q":"ZZ","WZjets":"WZ","WZ2l2q":"WZ","WZ1l1nu2q":"WZ","sChannel":"TsChan","tChannel":"TtChan","tbarChannel":"TbartChan","tWInclusive":"TtW","tbarWInclusive":"TbartW","tZq":"tZq","tHq":"THQ","ttWlnu":"TTW","ttW2q":"TTW","ttZ2l2nu":"TTZ","ttZ2q":"TTZ","ttbarInclusivePowerheg":"TT","tWZ":"TWZ","wPlusJets":"Wjets","DYJetsToLL_M-50":"DYToLL_M50","DYJetsToLL_M-10To50":"DYToLL_M10To50"}
    listOfMCs = {}

    #jetUnc = JetCorrectionUncertainty("../scaleFactors/2015/Fall15_25nsV2_MC_Uncertainty_AK4PFchs.txt")
    #if (is2016)
    jetUnc = JetCorrectionUncertainty("scaleFactors/2016/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt")

    #mapping of channels to dataTypes
    channelToDataset = {"ee":"DataEG","mumu":"DataMu"}

    #systematics list
    systs = ["","__trig__plus","__trig__minus","__jer__plus","__jer__minus","__jes__plus","__jes__minus","__pileup__plus","__pileup__minus","__bTag__plus","__bTag__minus","__met__plus","__met__minus","__pdf__plus","__pdf__minus","__ME_PS__plus","__ME_PS__minus","__alphaS__plus","__alphaS__minus"]

    #read what channel we're using here - changing this so that multiple things can be stored in the same file. i.e. should now be a list of channels to run over
    channels = eval(sys.argv[1])

    #Might make this customisable later, but don't really feel like I need to yet
    inputDir = "mvaTest/"
    if len(sys.argv) > 2:
        inputDir = sys.argv[2]

    outputDir = "mvaInputs/"
    if len(sys.argv) > 3:
        outputDir = sys.argv[3]
    inputVars = setupInputVars()

    is2016 = False
    if len(sys.argv) > 4 and sys.argv[4] == "--2016":
        is2016 = True
    elif len(sys.argv) > 4 and sys.argv[4] == "--2015":
        is2016 = False

    useSidebandRegion = False
    if len(sys.argv) > 5 and sys.argv[5] == "-s":
        useSidebandRegion = True
    treeNamePostfixSig = ""
    treeNamePostfixSB = ""
    if useSidebandRegion:
        print "Using control region"
        treeNamePostfixSig = "sig_"
        treeNamePostfixSB = "ctrl_"

    chanMap = {}
    if is2016 :
        chanMap = {"ee":"eeRun2016","mumu":"mumuRun2016"}
    else :
        chanMap = {"ee":"eeRun2015","mumu":"mumuRun2015"}

    outFakeChannels = ["FakeEG","FakeMu"]
    outFakeChanToData = {}
    outFakeChanToData["FakeEG"] = ["ee"]
    outFakeChanToData["FakeMu"] = ["mumu"]

    listOfMCs = {"ttHTobb" : "ttH", "ttHToNonbb" : "ttH", "WWW" : "WWW", "WWZ" : "WWZ", "WZZ" : "WZZ", "ZZZ" : "ZZZ", "WW1l1nu2q" : "WW", "WW2l2nu":"WW","ZZ4l":"ZZ","ZZ2l2nu":"ZZ","ZZ2l2q":"ZZ","WZjets":"WZ","WZ2l2q":"WZ","WZ1l1nu2q":"WZ","sChannel":"TsChan","tChannel":"TtChan","tbarChannel":"TbartChan","tWInclusive":"TtW","tbarWInclusive":"TbartW","tZq":"tZq","tHq":"THQ","ttWlnu":"TTW","ttW2q":"TTW","ttZ2l2nu":"TTZ","ttZ2q":"TTZ","ttbarInclusivePowerheg":"TT","tWZ":"TWZ","wPlusJets":"Wjets","DYJetsToLL_M-50":"DYToLL_M50","DYJetsToLL_M-10To50":"DYToLL_M10To50"}

    #Loop over opposite sign samples to create fake shape
    for outChan in outFakeChannels:
        print "And finally fake (non-prompt) lepton shapes estimated from data ",outChan
        outTreeSig = TTree("Ttree_"+treeNamePostfixSig+outChan,"Ttree_"+treeNamePostfixSig+outChan)
        setupBranches(outTreeSig,inputVars)
        outTreeSdBnd = 0
        if useSidebandRegion:
            outTreeSdBnd = TTree("Ttree_"+treeNamePostfixSB+outChan,"Ttree_"+treeNamePostfixSB+outChan)
            setupBranches(outTreeSdBnd,inputVars)
        outFile = TFile(outputDir+"histofile_"+outChan+".root","RECREATE")

        # Get same sign data
        for chan in outFakeChanToData[outChan]:
            dataChain = TChain("tree")
            if is2016 :
                 dataChain.Add(inputDir+chanMap[chan]+chan+"invLepmvaOut.root")
            else :
                for run in ["C","D"]:
                    dataChain.Add(inputDir+chanMap[chan]+run+chan+"mvaOut.root")

   	    # Get expected real SS events from MC
            for sample in listOfMCs.keys():
               print "Doing SS fakes " + sample + "\n",
               sys.stdout.flush()
               dataChain.Add(inputDir+sample+chan+"invLepmvaOut.root")
            try:
               fillTree(outTreeSig, outTreeSdBnd, inputVars, dataChain, outChan, 0, chan, is2016, True)
            except AttributeError:
                print "\nAttribute Error \n"

        outFile.cd()
        outFile.Write()
        outTreeSig.Write()
        if useSidebandRegion:
            outTreeSdBnd.Write()
        outFile.Close()

if __name__ == "__main__":
    main()

