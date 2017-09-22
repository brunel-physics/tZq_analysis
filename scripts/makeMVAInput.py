#A tool to pull a load of information from the mva trees. Should be all easy jazz...

from ROOT import *

import math
import sys
import os
from array import array
from jetCorrectionUncertainty import JetCorrectionUncertainty

def sortOutLeptons(tree,channel):
    ###Returns three TLorentzVectors containing the two z leptons and the w lepton. This will be VERY useful for making all of the plots.
    #Reads the position of the w and z leptons from variables stored at mvaTree making time, because I'm great and finally got around to doing it.
    zMass = 100
    zLep1,zLep2,wLep = 0,0,0
    #Let's try commenting this out and see if everything breaks? Hopefully it won't do...
    #if tree.numElePF2PAT < 3:
    #    return (0,0,0)
    if channel == "eee":
        wLep = TLorentzVector(tree.elePF2PATGsfPx[tree.wLepIndex],tree.elePF2PATGsfPy[tree.wLepIndex],tree.elePF2PATGsfPz[tree.wLepIndex],tree.elePF2PATGsfE[tree.wLepIndex])
        zLep1 = TLorentzVector(tree.elePF2PATGsfPx[tree.zLep1Index],tree.elePF2PATGsfPy[tree.zLep1Index],tree.elePF2PATGsfPz[tree.zLep1Index],tree.elePF2PATGsfE[tree.zLep1Index])
        zLep2 = TLorentzVector(tree.elePF2PATGsfPx[tree.zLep2Index],tree.elePF2PATGsfPy[tree.zLep2Index],tree.elePF2PATGsfPz[tree.zLep2Index],tree.elePF2PATGsfE[tree.zLep2Index])
    if channel == "eemu":
        zLep1 = TLorentzVector(tree.elePF2PATGsfPx[tree.zLep1Index],tree.elePF2PATGsfPy[tree.zLep1Index],tree.elePF2PATGsfPz[tree.zLep1Index],tree.elePF2PATGsfE[tree.zLep1Index])
        zLep2 = TLorentzVector(tree.elePF2PATGsfPx[tree.zLep2Index],tree.elePF2PATGsfPy[tree.zLep2Index],tree.elePF2PATGsfPz[tree.zLep2Index],tree.elePF2PATGsfE[tree.zLep2Index])
        wLep = TLorentzVector(tree.muonPF2PATPx[tree.wLepIndex],tree.muonPF2PATPy[tree.wLepIndex],tree.muonPF2PATPz[tree.wLepIndex],tree.muonPF2PATE[tree.wLepIndex])
    if channel == "emumu":
        zLep1 = TLorentzVector(tree.muonPF2PATPx[tree.zLep1Index],tree.muonPF2PATPy[tree.zLep1Index],tree.muonPF2PATPz[tree.zLep1Index],tree.muonPF2PATE[tree.zLep1Index])
        zLep2 = TLorentzVector(tree.muonPF2PATPx[tree.zLep2Index],tree.muonPF2PATPy[tree.zLep2Index],tree.muonPF2PATPz[tree.zLep2Index],tree.muonPF2PATE[tree.zLep2Index])
        wLep = TLorentzVector(tree.elePF2PATGsfPx[tree.wLepIndex],tree.elePF2PATGsfPy[tree.wLepIndex],tree.elePF2PATGsfPz[tree.wLepIndex],tree.elePF2PATGsfE[tree.wLepIndex])
    if channel == "mumumu":
        zLep1 = TLorentzVector(tree.muonPF2PATPx[tree.zLep1Index],tree.muonPF2PATPy[tree.zLep1Index],tree.muonPF2PATPz[tree.zLep1Index],tree.muonPF2PATE[tree.zLep1Index])
        zLep2 = TLorentzVector(tree.muonPF2PATPx[tree.zLep2Index],tree.muonPF2PATPy[tree.zLep2Index],tree.muonPF2PATPz[tree.zLep2Index],tree.muonPF2PATE[tree.zLep2Index])
        wLep = TLorentzVector(tree.muonPF2PATPx[tree.wLepIndex],tree.muonPF2PATPy[tree.wLepIndex],tree.muonPF2PATPz[tree.wLepIndex],tree.muonPF2PATE[tree.wLepIndex])
    return (zLep1,zLep2,wLep)

def getJets(tree,syst,jetUnc,met):
    #Makes a short list of indices of the jets in the event
    jetList = []
    jetVecList = []
    for i in range(15):
        if tree.jetInd[i] > -.5:
            jetList.append(tree.jetInd[i])
            jetVecList.append(getJetVec(tree,tree.jetInd[i],tree.jetSmearValue[i],met,syst,True))
        else: continue
    return (jetList,jetVecList)

def getBjets(tree,syst,jetUnc,met,jets):
    #Return a list of the indices of the b-jets in the event
    bJetList = []
    bJetVecList = []
    for i in range(10):
        if tree.bJetInd[i] > -0.5:
            bJetList.append(tree.bJetInd[i])
            bJetVecList.append(getJetVec(tree,jets[tree.bJetInd[i]],tree.jetSmearValue[tree.bJetInd[i]],met,systFalse))
        else:continue
#    print len(bJetList)
    return (bJetList,bJetVecList)

def getJetVec(tree, index, smearValue, metVec, syst, doMetSmear):
    #Gets a vector for a jet with corrections already applied

    returnJet = TLorentzVector();
    returnJet.SetPxPyPzE(tree.jetPF2PATPx[index],tree.jetPF2PATPy[index],tree.jetPF2PATPz[index],tree.jetPF2PATE[index]);
    returnJet *= smearValue;

    if syst == 16:
        returnJet *= 1+ jetUnc.getUncertainty(returnJet.Pt(), returnJet.Eta(),1)
    elif syst == 32:
        returnJet *= 1+ jetUnc.getUncertainty(returnJet.Pt(), returnJet.Eta(),2)

    if ( doMetSmear and smearValue > 0.01 ) :
    #Propogate through the met. But only do it if the smear jet isn't 0.
        metVec.SetPx(metVec.Px()+tree.jetPF2PATPx[index])
        metVec.SetPy(metVec.Py()+tree.jetPF2PATPy[index])

        metVec.SetPx(metVec.Px()-returnJet.Px())
        metVec.SetPy(metVec.Py()-returnJet.Py())

    return returnJet

def doUncMet(tree,met,zLep1,zLep2,wLep,jetVecs,syst):
    #Subtracts all items from met and then varies what's left by 10% for systematic purposes.
    uncMetX = met.Px() + zLep1.Px() + zLep2.Px() + wLep.Px()
    uncMetY = met.Py() + zLep1.Py() + zLep2.Py() + wLep.Py()
    for i in range(len(jetVecs)):
        uncMetX += jetVecs[i].Px()
        uncMetY += jetVecs[i].Py()
    if syst == 1028:
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
    inputVars["mTW"] = array('f',[0.])
    inputVars["leptWPt"] = array('f',[0.])
    inputVars["leadJetPt"] = array('f',[0.])
    inputVars["totPt"] = array('f',[0.])
    inputVars["totEta"] = array('f',[0.])
    inputVars["totPhi"] = array('f',[0.])
    inputVars["totPtVec"] = array('f',[0.])
    inputVars["totVecM"] = array('f',[0.])
    inputVars["chan"] = array('i',[0])
    inputVars["nJets"] = array('f',[0.])
    inputVars["nBjets"] = array('f',[0.])
    inputVars["met"] = array('f',[0.])
    inputVars["lepPt"] = array('f',[0.])
    inputVars["lepMetPt"] = array('f',[0.])
    inputVars["totPt2Jet"] = array('f',[0.])
    inputVars["leadJetbTag"] = array('f',[0.])
    inputVars["leadJetEta"] = array('f',[0.])
    inputVars["secJetbTag"] = array('f',[0.])
    inputVars["secJetPt"] = array('f',[0.])
    inputVars["secJetEta"] = array('f',[0.])
    inputVars["bTagDisc"] = array('f',[0.])
    inputVars["topMass"] = array('f',[0.])
    inputVars["topPt"] = array('f',[0.])
    inputVars["topEta"] = array('f',[0.])
    inputVars["zPt"] = array('f',[0.])
    inputVars["zEta"] = array('f',[0.])
    inputVars["wLepEta"] = array('f',[0.])
    inputVars["wZdelR"] = array('f',[0.])
    inputVars["j1j2delR"] = array('f',[0.])
    inputVars["minZJetR"] = array('f',[0.])
    inputVars["zWLepdelR"] = array('f',[0.])
    inputVars["zmetdelPhi"] = array('f',[0.])
    inputVars["zWLepdelPhi"] = array('f',[0.])
    inputVars["lbDelR"] = array('f',[0.])
    inputVars["lbDelPhi"] = array('f',[0.])
    inputVars["zlb1DelR"] = array('f',[0.])
    inputVars["zlb1DelPhi"] = array('f',[0.])
    inputVars["zlb2DelR"] = array('f',[0.])
    inputVars["zlb2DelPhi"] = array('f',[0.])
    inputVars["totHt"] = array('f',[0.])
    inputVars["lepHt"] = array('f',[0.])
    inputVars["jetHt"] = array('f',[0.])
    inputVars["jetMass"] = array('f',[0.])
    inputVars["jetPt"] = array('f',[0.])
    inputVars["jetEta"] = array('f',[0.])
    inputVars["jetPhi"] = array('f',[0.])
    inputVars["jetMass3"] = array('f',[0.])
    inputVars["lepMetHt"] = array('f',[0.])
    inputVars["totHtOverPt"] = array('f',[0.])
    inputVars["zMass"] = array('f',[0.])
    return inputVars

def setupBranches(tree,varMap):
    tree.Branch("tree_EvtWeight", varMap["eventWeight"], "tree_EvtWeight/F")
    tree.Branch("tree_mTW", varMap["mTW"], "tree_mTW/F")
    tree.Branch("tree_leptWPt", varMap["leptWPt"], "tree_leptWPt/F")
    tree.Branch("tree_leadJetPt",varMap["leadJetPt"],"tree_leadJetPt/F")
    tree.Branch("tree_totPt",varMap["totPt"],"tree_totPt/F")
    tree.Branch("tree_totEta",varMap["totEta"],"tree_totEta/F")
    tree.Branch("tree_totPhi",varMap["totPhi"],"tree_totPhi/F")
    tree.Branch("tree_totPtVec",varMap["totPtVec"],"tree_totPtVec/F")
    tree.Branch("tree_totVecM",varMap["totVecM"],"tree_totVecM/F")
    tree.Branch("tree_Channel",varMap["chan"],"tree_Channel/I")
    tree.Branch("tree_NJets",varMap["nJets"],"tree_NJets/F")
    tree.Branch("tree_NBJets",varMap["nBjets"],"tree_NBJets/F")
    tree.Branch("tree_met",varMap["met"],"tree_met/F")
    tree.Branch("tree_lepPt",varMap["lepPt"],"tree_lepPt/F")
    tree.Branch("tree_lepMetPt",varMap["lepMetPt"],"tree_lepMetPt/F")
    tree.Branch("tree_totPt2Jet",varMap["totPt2Jet"],"tree_totPt2Jet/F")
    tree.Branch("tree_btagDiscri",varMap["bTagDisc"],"btagDiscri/F")
    tree.Branch("tree_leadJetbTag",varMap["leadJetbTag"],"leadJetbTag/F")
    tree.Branch("tree_leadJetEta",varMap["leadJetEta"],"leadJetEta/F")
    tree.Branch("tree_secJetbTag",varMap["secJetbTag"],"secJetbTag/F")
    tree.Branch("tree_secJetPt",varMap["secJetPt"],"secJetPt/F")
    tree.Branch("tree_secJetEta",varMap["secJetEta"],"secJetEta/F")
    tree.Branch("tree_topMass",varMap["topMass"],"tree_topMass/F")
    tree.Branch("tree_topPt",varMap["topPt"],"tree_topPt/F")
    tree.Branch("tree_topEta",varMap["topEta"],"tree_topEta/F")
    tree.Branch("tree_Zpt",varMap["zPt"],"tree_Zpt/F")
    tree.Branch("tree_Zeta",varMap["zEta"],"tree_Zeta/F")
    tree.Branch("tree_leptWEta",varMap["wLepEta"],"tree_leptWEta/F")
    tree.Branch("tree_wzdelR",varMap["wZdelR"],"tree_wzdelR/F")
    tree.Branch("tree_jjdelR",varMap["j1j2delR"],"tree_jjdelR/F")
    tree.Branch("tree_zjminR",varMap["minZJetR"],"tree_zjminR/F")
    tree.Branch("tree_ZlepWdelPhi",varMap["zWLepdelPhi"],"tree_ZlepWdelPhi/F")
    tree.Branch("tree_ZmetdelPhi",varMap["zmetdelPhi"],"tree_ZmetdelPhi/F")
    tree.Branch("tree_ZlepWdelR",varMap["zWLepdelR"],"tree_ZlepWdelR/F")
    tree.Branch("tree_lbDelR",varMap["lbDelR"],"tree_lbDelR/F")
    tree.Branch("tree_lbDelPhi",varMap["lbDelPhi"],"tree_lbDelPhi/F")
    tree.Branch("tree_zlb1DelR",varMap["zlb1DelR"],"tree_zlb1DelR/F")
    tree.Branch("tree_zlb1DelPhi",varMap["zlb1DelPhi"],"tree_zlb1DelPhi/F")
    tree.Branch("tree_zlb2DelR",varMap["zlb2DelR"],"tree_zlb2DelR/F")
    tree.Branch("tree_zlb2DelPhi",varMap["zlb2DelPhi"],"tree_zlb2DelPhi/F")
    tree.Branch("tree_totHt",varMap["totHt"],"tree_totHt/F")
    tree.Branch("tree_lepHt",varMap["lepHt"],"tree_lepHt/F")
    tree.Branch("tree_jetHt",varMap["jetHt"],"tree_jetHt/F")
    tree.Branch("tree_jetMass",varMap["jetMass"],"tree_jetMass/F")
    tree.Branch("tree_jetPt",varMap["jetPt"],"tree_jetPt/F")
    tree.Branch("tree_jetEta",varMap["jetEta"],"tree_jetEta/F")
    tree.Branch("tree_jetPhi",varMap["jetPhi"],"tree_jetPhi/F")
    tree.Branch("tree_jetMass3",varMap["jetMass3"],"tree_jetMass3/F")
    tree.Branch("tree_lepMetHt",varMap["lepMetHt"],"tree_lepMetHt/F")
    tree.Branch("tree_totHtOverPt",varMap["totHtOverPt"],"tree_totHtOverPt/F")
    tree.Branch("tree_zMass",varMap["zMass"],"tree_zMass/F")



def fillTree(outTree, varMap, tree, label, channel, jetUnc, overRideWeight = -1., zPtEventWeight = 0.):
    #Fills the output tree. This is a new function because I want to access data and MC in different ways but do the same thing to them in the end.

    syst = 0
    if "__jer__plus" in label:
        syst = 4
    if "__jer__minus" in label:
        syst = 8
    if "__jes__plus" in label:
        syst = 16
    if "__jes__minus" in label:
        syst = 32
    if "__met__plus" in label:
        syst = 1024
    if "__met__minus" in label:
        syst = 2048
    if channel == "eee":
        varMap["chan"][0] = 3
    if channel == "eemu":
        varMap["chan"][0] = 2
    if channel == "emumu":
        varMap["chan"][0] = 1
    if channel == "mumumu":
        varMap["chan"][0] = 0
    #        topMass
    #Set up those variables as branches
    for event in range(tree.GetEntries()):
            #Fill some plots here. Let's make an example mTW plot.
            #Make a config that'll do this for me? I've done these before so should be easy. Fill expressions could be a pain?
        tree.GetEntry(event)
        (zLep1,zLep2,wLep) = sortOutLeptons(tree,channel)
        metVec = TLorentzVector(tree.metPF2PATPx,tree.metPF2PATPy,0,tree.metPF2PATEt)
        (jets,jetVecs) = getJets(tree,syst,jetUnc,metVec)
        (bJets,bJetVecs) = getBjets(tree,syst,jetUnc,metVec,jets)
        #Do unclustered met stuff here now that we have all of the objects, all corrected for their various SFs etc.
        if syst == 1024 or syst == 2048:
            metVec = doUncMet(tree,metVec,zLep1,zLep2,wLep,jetVecs,syst)
        if wLep == 0:
            continue
        if tree.eventWeight == tree.eventWeight:
            varMap["eventWeight"][0] = tree.eventWeight
            if overRideWeight > 0.:
                varMap["eventWeight"][0] = tree.eventWeight * overRideWeight
            if zPtEventWeight > 0.1:
                varMap["eventWeight"][0] *= tree.eventWeight
            if zPtEventWeight < -0.1:
                varMap["eventWeight"][0] = 1.
        else:
            varMap["eventWeight"][0] = 0.

        if varMap["eventWeight"][0] < 0.:
            varMap["eventWeight"][0] = 0.
        varMap["leptWPt"][0] = wLep.Pt()
        varMap["wLepEta"][0] = wLep.Eta()
        varMap["leadJetPt"][0] = jetVecs[0].Pt()
        varMap["leadJetEta"][0] = jetVecs[0].Eta()
        #Make all the random Pt variables I'm saving for some reason
        totPx,totPy = 0,0
        totPx += zLep1.Px() + zLep2.Px() + wLep.Px()
        totPy += zLep1.Py() + zLep2.Py() + wLep.Py()
        varMap["lepPt"][0] = math.sqrt(totPx * totPx + totPy * totPy)
        totPx += metVec.Px()
        totPy += metVec.Py()
        varMap["lepMetPt"][0] = math.sqrt(totPx * totPx + totPy * totPy)
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
        totVec = (wLep + zLep1+zLep2)
        for i in range(len(jetVecs)):
            totVec += jetVecs[i]
        varMap["totEta"][0] = totVec.Eta()
        varMap["totPtVec"][0] = totVec.Pt()
        varMap["totVecM"][0] = totVec.M()
        varMap["mTW"][0] = math.sqrt(2*metVec.Pt()*wLep.Pt() * (1-math.cos(metVec.Phi() - wLep.Phi())))
        varMap["nJets"][0] = float(len(jets))
        varMap["nBjets"][0] = float(len(bJets))
        varMap["met"][0] = metVec.Pt()
        varMap["bTagDisc"][0] = tree.jetPF2PATBDiscriminator[jets[bJets[0]]]
        varMap["leadJetbTag"][0] = tree.jetPF2PATBDiscriminator[jets[0]]
        varMap["secJetbTag"][0] = -10.
        varMap["secJetPt"][0] = -1.
        varMap["secJetEta"][0] = -500.
        if len(jetVecs) > 1:
            varMap["secJetbTag"][0] = tree.jetPF2PATBDiscriminator[jets[1]]
            varMap["secJetPt"][0] = jetVecs[1].Pt()
            varMap["secJetEta"][0] = jetVecs[1].Eta()

#        print bTagDisc[0], bJets[0], tree.jetPF2PATBDiscriminator[jets[bJets[0]]], len(bJets), nBjets[0]
        varMap["topMass"][0] = (bJetVecs[0] + metVec + wLep).M()
        varMap["topPt"][0] = (bJetVecs[0] + metVec + wLep).Pt()
        varMap["topEta"][0] = (bJetVecs[0] + metVec + wLep).Eta()
        varMap["zPt"][0] = (zLep2 + zLep1).Pt()
        varMap["zEta"][0] = (zLep2 + zLep1).Eta()
        varMap["wZdelR"][0] = (zLep2 + zLep1).DeltaR(metVec + wLep)
        varMap["j1j2delR"][0] = -1.
        if len(jetVecs) > 1:
            varMap["j1j2delR"][0] = jetVecs[0].DeltaR(jetVecs[1])
        varMap["minZJetR"][0] = 3.0
        jetHt = 0.
        jetVector = TLorentzVector();
        for i in range(len(jetVecs)):
            jetHt += jetVecs[i].Pt()
            jetVector += jetVecs[i]
            if jetVecs[i].DeltaR(zLep2 + zLep1) < varMap["minZJetR"][0]:
                varMap["minZJetR"][0] = jetVecs[i].DeltaR(zLep2 + zLep1)
        outTree.Fill()
        varMap["zWLepdelR"][0] = (zLep2 + zLep1).DeltaR(wLep)
        varMap["zmetdelPhi"][0] = (zLep2+zLep1).DeltaPhi(metVec)
        varMap["zWLepdelPhi"][0] = (zLep2 + zLep1).DeltaPhi(wLep)
        varMap["lbDelR"][0] = wLep.DeltaR(jetVecs[bJets[0]])
        varMap["lbDelPhi"][0] = wLep.DeltaPhi(jetVecs[bJets[0]])
        varMap["zlb1DelR"][0] = zLep1.DeltaR(jetVecs[bJets[0]])
        varMap["zlb1DelPhi"][0] = zLep1.DeltaPhi(jetVecs[bJets[0]])
        varMap["zlb2DelR"][0] = zLep2.DeltaR(jetVecs[bJets[0]])
        varMap["zlb2DelPhi"][0] = zLep2.DeltaPhi(jetVecs[bJets[0]])
        ht = 0.
        ht += zLep1.Pt() + zLep2.Pt() + wLep.Pt()
        varMap["lepHt"][0] = ht
        ht += metVec.Pt()
        varMap["jetHt"][0] = jetHt
        varMap["jetMass"][0] = jetVector.M()
        varMap["jetPt"][0] = jetVector.Pt()
        varMap["jetEta"][0] = jetVector.Eta()
        varMap["jetPhi"][0] = jetVector.Phi()
        varMap["jetMass3"][0] = (jetVecs[0] + jetVecs[1] + jetVecs[2]).M()
        varMap["lepMetHt"][0] = ht
        ht += jetHt
        varMap["totHt"][0] = ht
        varMap["totHtOverPt"][0] = ht / math.sqrt(totPx * totPx + totPy * totPy)
        varMap["zMass"][0] = (zLep1+zLep2).M()



def main():

    zEnrichWeights = {"mumumu":-1.,"emumu":-1.,"eemu":-1.,"eee":-1.}
#    zEnrichWeights = {"mumumu":0.0001,"emumu":0.496,"eemu":0.0283,"eee":0.211}
    #lepton selection stage weights
#    zEnrichWeights = {"mumumu":0.057,"emumu":0.612,"eemu":0.0637,"eee":0.216}

    #Mapping of our mc names to IPHC names
    listOfMCs = {"WW2l2nu":"WW","ZZ4l":"ZZ","sChannel":"TsChan","tChannel":"TtChan","tbarChannel":"TbartChan","tWInclusive":"TtW","tbarWInclusive":"TbartW","tZq":"tZq","ttW":"TTW","ttZ":"TTZ","wPlusJets":"Wjets"}
#    listOfMCs = {"WW2l2nu":"WW","WZ3l1nu":"WZ","ZZ4l":"ZZ","sChannel":"TsChan","sbarChannel":"TbarsChan","tChannel":"TtChan","tbarChannel":"TbartChan","tWInclusive":"TtW","tbarWInclusive":"TbartW","tZq":"tZq","tZq4Flavour3Lepton":"tZq4f", "ttW":"TTW","ttZ":"TTZ","ttbarDilepton":"TT","wPlusJets":"Wjets", "zPlusJets10To50Filter":"DYToLL_M10-50","zPlusJetsTuneZ2Star":"Zjets"}
#    listOfMCs = {"WW2l2nu":"WW","WZ3l1nu":"WZ","WZ2l2nu":"WZ","ZZ2l2q":"ZZ","ZZ4l":"ZZ","sChannel":"TsChan","sbarChannel":"TbarsChan","tChannel":"TtChan","tbarChannel":"TbartChan","tWInclusive":"TtW","tbarWInclusive":"TbartW","ttW":"TTW","ttZ":"TTZ","ttbarDilepton":"TT","wPlusJets":"Wjets", "zPlusJets10To50Filter":"DYToLL_M10-50","zPlusJetsTuneZ2Star":"Zjets"}
#    listOfMCs = {"ZZ4l":"ZZ"}
#    listOfMCs = {"ttW":"TTW"}
#    listOfMCs = {"WZ3l1nu":"WZ"}
#    listOfMCs = {"WZ3l1nu":"WZ","zPlusJetsTuneZ2Star":"Zjets"}
#    listOfMCs = {"zPlusJets10To50Filter":"DYToLL_M10-50"}
#    listOfMCs = {}

    #Set-up JEC corrections
    jetUnc = JetCorrectionUncertainty("../scaleFactors/2015/Fall15_25nsV2_MC_Uncertainty_AK4PFchs.txt")

    #mapping of channels to dataTypes
    channelToDataset = {"eee":"DataEG","eemu":"DataMuEG","emumu":"DataMuEG","mumumu":"DataMu"}

    #systematics list
    systs = ["","__trig__plus","__trig__minus","__jer__plus","__jer__minus","__jes__plus","__jes__minus","__pileup__plus","__pileup__minus","__met__plus","__met__minus","__bTag__plus","__bTag__minus","__pdf__plus","__pdf__minus"]

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

    #Loop over samples
    for sample in listOfMCs.keys():
        overrideWeight = -1.
        if sample == "WZ2l2nu" or sample == "ZZ2l2p":
            continue
        print "Doing " + sample + ": ",
        sys.stdout.flush()
        if "WZ" in sample:
            overrideWeight = 1286770./2133868.

        outFile = 0
        #update the appropriate root file
        outFile = TFile(outputDir+"histofile_"+listOfMCs[sample] + ".root","UPDATE")

        for syst in systs:
            #We now define the outtree out here, coz it seems like a more sensible option.
            outTree = TTree("Ttree_"+listOfMCs[sample]+syst, "Ttree_"+listOfMCs[sample]+syst)
            setupBranches(outTree,inputVars)
            for channel in channels:
                inFile = TFile(inputDir+sample+channel+"mvaOut.root","READ")
                if "met" in syst:
                    tree = inFile.Get("tree")
                else:
                    tree = inFile.Get("tree"+syst)
                    try:
                        print syst + ": " + str(tree.GetEntriesFast())
                        sys.stdout.flush()
                        fillTree(outTree, inputVars, tree, listOfMCs[sample]+syst, channel, jetUnc, overRideWeight = overrideWeight)
                    except AttributeError:
                        print syst + ": " + "0",
                        sys.stdout.flush()
                #Various stuff needs to be saved in the same trees. Create new one if it doesn't exist, open current one if it does
                inFile.Close()
            outFile.cd()
            outFile.Write()
            outTree.Write()
        #if tree exists just update that.
        #        if outFile.GetListOfKeys().Contains("Ttree_"+listOfMCs[sample]):
        #            outTree = outFile.Get("Ttree_"+listOfMCs[sample])
        #        else:
    #next do the data files
        outFile.Write()
        outFile.Close()
        print
    chanMap = {"eee":"eeRun2015","eemu":"emuRun2015","emumu":"emuRun2015","mumumu":"mumuRun2015"}

    outChannels = ["DataEG","DataMuEG","DataMu"]
    outChanToData = {}
    outChanToData["DataEG"] = ["eee"]
    outChanToData["DataMuEG"] = ["eemu","emumu"]
    outChanToData["DataMu"] = ["mumumu"]

    for outChan in outChannels:
        print "Data ",outChan
        outTree = TTree("Ttree_"+outChan,"Ttree_"+outChan)
        setupBranches(outTree,inputVars)
        outFile = TFile(outputDir+"histofile_"+outChan+".root","UPDATE")
        for chan in outChanToData[outChan]:
            dataChain = TChain("tree")
#            for run in ["A","B","C","D"]:
            for run in ["C","D"]:
                dataChain.Add(inputDir+chanMap[chan]+run+chan+"mvaOut.root")
            fillTree(outTree, inputVars, dataChain, outChan, chan, 0)
        outFile.cd()
        outFile.Write()
        outTree.Write()
        outFile.Close()

    zEnrichSyst = ["","__zPt__plus","__zPt__minus"]
    for outChan in outChannels:
        print "And finally z-enriched data ",outChan
        outFileZ = TFile(outputDir+"histofile_"+outChan+"Zenriched.root","UPDATE")
        for systPost in zEnrichSyst:
            outTreeZ = TTree("Ttree_"+outChan+"Zenriched"+systPost,"Ttree_"+outChan+"Zenriched"+systPost)
            setupBranches(outTreeZ,inputVars)

            for chan in outChanToData[outChan]:
                dataChainZ = TChain("tree")
#                for run in ["A","B","C","D"]:
                for run in ["C","D"]:
                    dataChainZ.Add(inputDir+chanMap[chan]+run+chan+"invIsomvaOut.root")
                zPtSyst = 0.
                if "plus" in systPost:
                    zPtSyst = 1.
                if "minus" in systPost:
                    zPtSyst = -1.
                fillTree(outTreeZ, inputVars, dataChainZ, outChan, chan, 0, zEnrichWeights[chan],zPtEventWeight = zPtSyst)
            outFileZ.cd()
            outFileZ.Write()
            outTreeZ.Write()
        outFileZ.Close()

if __name__ == "__main__":
    main()
