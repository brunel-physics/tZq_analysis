//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jan 23 14:04:00 2014 by ROOT version 5.32/00
// from TChain tree/
//////////////////////////////////////////////////////////

#ifndef _AnalysisEvent_hpp_
#define _AnalysisEvent_hpp_

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <string>
#include <TLorentzVector.h>
#include <iostream>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class AnalysisEvent {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           numElePF2PAT;
   Int_t           numLooseElePF2PAT;
   Float_t         elePF2PATE[20];   //[numElePF2PAT]
   Float_t         elePF2PATET[20];   //[numElePF2PAT]
   Float_t         elePF2PATPX[20];   //[numElePF2PAT]
   Float_t         elePF2PATPY[20];   //[numElePF2PAT]
   Float_t         elePF2PATPZ[20];   //[numElePF2PAT]
   Float_t         elePF2PATPhi[20];   //[numElePF2PAT]
   Float_t         elePF2PATTheta[20];   //[numElePF2PAT]
   Float_t         elePF2PATEta[20];   //[numElePF2PAT]
   Float_t         elePF2PATPT[20];   //[numElePF2PAT]
   Int_t           elePF2PATCharge[20];   //[numElePF2PAT]
   Float_t         elePF2PATMVA[20];   //[numElePF2PAT]
   Float_t         elePF2PATImpactTransDist[20];   //[numElePF2PAT]
   Float_t         elePF2PATImpactTransError[20];   //[numElePF2PAT]
   Float_t         elePF2PATImpactTransSignificance[20];   //[numElePF2PAT]
   Float_t         elePF2PATImpact3DDist[20];   //[numElePF2PAT]
   Float_t         elePF2PATImpact3DError[20];   //[numElePF2PAT]
   Float_t         elePF2PATImpact3DSignificance[20];   //[numElePF2PAT]
   Float_t         elePF2PATChargedHadronIso[20];   //[numElePF2PAT]
   Float_t         elePF2PATNeutralHadronIso[20];   //[numElePF2PAT]
   Float_t         elePF2PATPhotonIso[20];   //[numElePF2PAT]
   Float_t         elePF2PATTrackPt[20];   //[numElePF2PAT]
   Float_t         elePF2PATTrackPhi[20];   //[numElePF2PAT]
   Float_t         elePF2PATTrackEta[20];   //[numElePF2PAT]
   Float_t         elePF2PATTrackChi2[20];   //[numElePF2PAT]
   Float_t         elePF2PATTrackNDOF[20];   //[numElePF2PAT]
   Float_t         elePF2PATTrackD0[20];   //[numElePF2PAT]
   Float_t         elePF2PATTrackDBD0[20];   //[numElePF2PAT]
   Float_t         elePF2PATD0PV[20];   //[numElePF2PAT]
   Float_t         elePF2PATDZPV[20];   //[numElePF2PAT]
   Float_t         elePF2PATBeamSpotCorrectedTrackD0[20];   //[numElePF2PAT]
   Float_t         elePF2PATTrackDz[20];   //[numElePF2PAT]
   Float_t         elePF2PATVtxZ[20];   //[numElePF2PAT]
   Int_t           elePF2PATIsGsf[20];   //[numElePF2PAT]
   Float_t         elePF2PATGsfPx[20];   //[numElePF2PAT]
   Float_t         elePF2PATGsfPy[20];   //[numElePF2PAT]
   Float_t         elePF2PATGsfPz[20];   //[numElePF2PAT]
   Float_t         elePF2PATGsfE[20];   //[numElePF2PAT]
   Float_t         elePF2PATSCEta[20];   //[numElePF2PAT]
   Float_t         elePF2PATSCE[20];   //[numElePF2PAT]
   Float_t         elePF2PATSCPhi[20];   //[numElePF2PAT]
   Float_t         elePF2PATSCSigmaEtaEta[20];   //[numElePF2PAT]
   Float_t         elePF2PATSCSigmaIEtaIEta[20];   //[numElePF2PAT]
   Float_t         elePF2PATSCE1x5[20];   //[numElePF2PAT]
   Float_t         elePF2PATSCE5x5[20];   //[numElePF2PAT]
   Float_t         elePF2PATSCE2x5max[20];   //[numElePF2PAT]
   Float_t         elePF2PATTrackIso04[20];   //[numElePF2PAT]
   Float_t         elePF2PATEcalIso04[20];   //[numElePF2PAT]
   Float_t         elePF2PATHcalIso04[20];   //[numElePF2PAT]
   Float_t         elePF2PATTrackIso03[20];   //[numElePF2PAT]
   Float_t         elePF2PATEcalIso03[20];   //[numElePF2PAT]
   Float_t         elePF2PATHcalIso03[20];   //[numElePF2PAT]
   Float_t         elePF2PATdr04EcalRecHitSumEt[20];   //[numElePF2PAT]
   Float_t         elePF2PATdr03EcalRecHitSumEt[20];   //[numElePF2PAT]
   Float_t         elePF2PATEcalIsoDeposit[20];   //[numElePF2PAT]
   Float_t         elePF2PATHcalIsoDeposit[20];   //[numElePF2PAT]
   Float_t         elePF2PATComRelIso[20];   //[numElePF2PAT]
   Float_t         elePF2PATComRelIsodBeta[20];   //[numElePF2PAT]
   Float_t         elePF2PATComRelIsoRho[20];   //[numElePF2PAT]
   Float_t         elePF2PATChHadIso[20];   //[numElePF2PAT]
   Float_t         elePF2PATNtHadIso[20];   //[numElePF2PAT]
   Float_t         elePF2PATGammaIso[20];   //[numElePF2PAT]
   Float_t         elePF2PATRhoIso[20];   //[numElePF2PAT]
   Float_t         elePF2PATAEff03[20];   //[numElePF2PAT]
   Int_t           elePF2PATMissingInnerLayers[20];   //[numElePF2PAT]
   Float_t         elePF2PATHoverE[20];   //[numElePF2PAT]
   Float_t         elePF2PATDeltaPhiSC[20];   //[numElePF2PAT]
   Float_t         elePF2PATDeltaEtaSC[20];   //[numElePF2PAT]
   Int_t           elePF2PATIsBarrel[20];   //[numElePF2PAT]
   Int_t           elePF2PATPhotonConversionTag[20];   //[numElePF2PAT]
   Float_t         elePF2PATPhotonConversionDist[20];   //[numElePF2PAT]
   Float_t         elePF2PATPhotonConversionDcot[20];   //[numElePF2PAT]
   Int_t           elePF2PATPhotonConversionVeto[20];   //[numElePF2PAT]
   Int_t           elePF2PATPhotonConversionTagCustom[20];   //[numElePF2PAT]
   Float_t         elePF2PATPhotonConversionDistCustom[20];   //[numElePF2PAT]
   Float_t         elePF2PATPhotonConversionDcotCustom[20];   //[numElePF2PAT]
   Float_t         elePF2PATTriggerMatch[20];   //[numElePF2PAT]
   Float_t         elePF2PATJetOverlap[20];   //[numElePF2PAT]
   Float_t         genElePF2PATET[20];   //[numElePF2PAT]
   Float_t         genElePF2PATPX[20];   //[numElePF2PAT]
   Float_t         genElePF2PATPY[20];   //[numElePF2PAT]
   Float_t         genElePF2PATPZ[20];   //[numElePF2PAT]
   Float_t         genElePF2PATPhi[20];   //[numElePF2PAT]
   Float_t         genElePF2PATTheta[20];   //[numElePF2PAT]
   Float_t         genElePF2PATEta[20];   //[numElePF2PAT]
   Int_t           genElePF2PATCharge[20];   //[numElePF2PAT]
   Float_t         elePF2PATlooseElectronSortedEt[20];   //[numLooseElePF2PAT]
   Float_t         elePF2PATlooseElectronSortedPt[20];   //[numLooseElePF2PAT]
   Float_t         elePF2PATlooseElectronSortedEta[20];   //[numLooseElePF2PAT]
   Float_t         elePF2PATlooseElectronSortedMVA[20];   //[numLooseElePF2PAT]
   Float_t         elePF2PATlooseElectronSortedRelIso[20];   //[numLooseElePF2PAT]
   Int_t           numMuonPF2PAT;
   Int_t           numLooseMuonPF2PAT;
   Float_t         muonPF2PATE[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATET[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATPt[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATPX[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATPY[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATPZ[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATPhi[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATTheta[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATEta[20];   //[numMuonPF2PAT]
   Int_t           muonPF2PATCharge[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATGlobalID[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATTrackID[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATChi2[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATD0[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATTrackDBD0[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATDBInnerTrackD0[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATBeamSpotCorrectedD0[20];   //[numMuonPF2PAT]
   Int_t           muonPF2PATTrackNHits[20];   //[numMuonPF2PAT]
   Int_t           muonPF2PATMuonNHits[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATNDOF[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATVertX[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATVertY[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATVertZ[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATChargedHadronIso[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATNeutralHadronIso[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATPhotonIso[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATTrackIso[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATEcalIso[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATHcalIso[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATComRelIso[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATComRelIsodBeta[20];   //[numMuonPF2PAT]
   Int_t           muonPF2PATIsPFMuon[20];   //[numMuonPF2PAT]
   Int_t           muonPF2PATNChambers[20];   //[numMuonPF2PAT]
   Int_t           muonPF2PATNMatches[20];   //[numMuonPF2PAT]
   Int_t           muonPF2PATTkLysWithMeasurements[20];   //[numMuonPF2PAT]
   Int_t           muonPF2PATVldPixHits[20];   //[numMuonPF2PAT]
   Int_t           muonPF2PATMatchedStations[20];   //[numMuonPF2PAT]
   Float_t         muonPF2PATGlbTkNormChi2[20];   //[numLooseMuonPF2PAT]
   Float_t         muonPF2PATDBPV[20];   //[numLooseMuonPF2PAT]
   Float_t         muonPF2PATDZPV[20];   //[numLooseMuonPF2PAT]
   Float_t         muonPF2PATlooseMuonSortedEt[20];   //[numLooseMuonPF2PAT]
   Float_t         muonPF2PATlooseMuonSortedPt[20];   //[numLooseMuonPF2PAT]
   Float_t         muonPF2PATlooseMuonSortedEta[20];   //[numLooseMuonPF2PAT]
   Float_t         muonPF2PATlooseMuonSortedRelIso[20];   //[numLooseMuonPF2PAT]
   Float_t         muonPF2PATlooseMuonSortedisGlb[20];   //[numLooseMuonPF2PAT]
   Float_t         muonPF2PATlooseMuonSortedisTrk[20];   //[numLooseMuonPF2PAT]
   Float_t         genMuonPF2PATET[20];   //[numMuonPF2PAT]
   Float_t         genMuonPF2PATPX[20];   //[numMuonPF2PAT]
   Float_t         genMuonPF2PATPY[20];   //[numMuonPF2PAT]
   Float_t         genMuonPF2PATPZ[20];   //[numMuonPF2PAT]
   Float_t         genMuonPF2PATPhi[20];   //[numMuonPF2PAT]
   Float_t         genMuonPF2PATTheta[20];   //[numMuonPF2PAT]
   Float_t         genMuonPF2PATEta[20];   //[numMuonPF2PAT]
   Int_t           genMuonPF2PATCharge[20];   //[numMuonPF2PAT]
   Int_t           numJetPF2PAT;
   Double_t        jetPF2PATE[40];   //[numJetPF2PAT]
   Double_t        jetPF2PATEt[40];   //[numJetPF2PAT]
   Double_t        jetPF2PATPt[40];   //[numJetPF2PAT]
   Double_t        jetPF2PATPtRaw[40];   //[numJetPF2PAT]
   Double_t        jetPF2PATUnCorEt[40];   //[numJetPF2PAT]
   Double_t        jetPF2PATUnCorPt[40];   //[numJetPF2PAT]
   Double_t        jetPF2PATEta[40];   //[numJetPF2PAT]
   Double_t        jetPF2PATTheta[40];   //[numJetPF2PAT]
   Double_t        jetPF2PATPhi[40];   //[numJetPF2PAT]
   Double_t        jetPF2PATPx[40];   //[numJetPF2PAT]
   Double_t        jetPF2PATPy[40];   //[numJetPF2PAT]
   Double_t        jetPF2PATPz[40];   //[numJetPF2PAT]
   Double_t        jetPF2PATdRClosestLepton[40];   //[numJetPF2PAT]
   Int_t           jetPF2PATNtracksInJet[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATJetCharge[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATfHPD[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATBtagSoftMuonPtRel[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATBtagSoftMuonQuality[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATCorrFactor[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATCorrResidual[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATL2L3ResErr[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATCorrErrLow[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATCorrErrHi[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATN90Hits[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATTriggered[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATSVPT[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATSVL2D[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATSVL2Dxy[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATSVL2DxyErr[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATSVL2DxySig[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATSVL3D[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATSVL3DErr[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATSVL3DSig[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATSVMass[40];   //[numJetPF2PAT]
   Int_t           jetPF2PATSVNtracks[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATSVX[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATSVY[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATSVZ[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATSVDX[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATSVDY[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATSVDZ[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATBDiscriminator[40];   //[numJetPF2PAT]
   Int_t           jetPF2PATNConstituents[40];   //[numJetPF2PAT]
   Int_t           jetPF2PATPID[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATClosestBPartonDeltaR[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATClosestCPartonDeltaR[40];   //[numJetPF2PAT]
   Float_t         genJetPF2PATET[40];   //[numJetPF2PAT]
   Float_t         genJetPF2PATPT[40];   //[numJetPF2PAT]
   Float_t         genJetPF2PATPX[40];   //[numJetPF2PAT]
   Float_t         genJetPF2PATPY[40];   //[numJetPF2PAT]
   Float_t         genJetPF2PATPZ[40];   //[numJetPF2PAT]
   Float_t         genJetPF2PATPhi[40];   //[numJetPF2PAT]
   Float_t         genJetPF2PATTheta[40];   //[numJetPF2PAT]
   Float_t         genJetPF2PATEta[40];   //[numJetPF2PAT]
   Int_t           genJetPF2PATPID[40];   //[numJetPF2PAT]
   Int_t           numLooseJetPF2PAT;
   Float_t         looseJetPF2PATEt[40];   //[numLooseBJetsPF2PAT]
   Float_t         looseJetPF2PATPt[40];   //[numLooseBJetsPF2PAT]
   Float_t         looseJetPF2PATEta[40];   //[numLooseBJetsPF2PAT]
   Float_t         looseJetPF2PATBDisc[40];   //[numLooseBJetsPF2PAT]
   Float_t         jetPF2PATBtagDisc_trackCountingHighPurBJetTags[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATBtagDisc_trackCountingHighEffBJetTags[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATBtagDisc_jetProbabilityBJetTags[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATBtagDisc_jetBProbabilityBJetTags[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATBtagDisc_softElectronBJetTags[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATBtagDisc_softMuonBJetTags[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATBtagDisc_softMuonNoIPBJetTags[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATBtagDisc_simpleSecondaryVertexHighEffBJetTags[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATBtagDisc_simpleSecondaryVertexNegativeBJetTags[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATBtagDisc_combinedSecondaryVertexMVABJetTags[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATBtagDisc_simpleSecondaryVertexHighPurBJetTags[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATMuEnergy[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATMuEnergyFraction[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATNeutralHadEnergy[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATNeutralEmEnergy[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATChargedHadronEnergyFraction[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATNeutralHadronEnergyFraction[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATChargedEmEnergyFraction[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATNeutralEmEnergyFraction[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATChargedHadronEnergyFractionCorr[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATNeutralHadronEnergyFractionCorr[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATChargedEmEnergyFractionCorr[40];   //[numJetPF2PAT]
   Float_t         jetPF2PATNeutralEmEnergyFractionCorr[40];   //[numJetPF2PAT]
   Int_t           jetPF2PATNeutralMultiplicity[40];   //[numJetPF2PAT]
   Int_t           jetPF2PATChargedMultiplicity[40];   //[numJetPF2PAT]
   Double_t        metPF2PATEt;
   Double_t        metPF2PATEtRaw;
   Double_t        metPF2PATPhi;
   Double_t        metPF2PATPt;
   Double_t        metPF2PATPx;
   Double_t        metPF2PATPy;
   Float_t         metPF2PATScalarEt;
   Float_t         metPF2PATEtUncorrected;
   Float_t         metPF2PATPhiUncorrected;
   Float_t         genMetPF2PATEt;
   Float_t         genMetPF2PATPhi;
   Float_t         genMetPF2PATPt;
   Float_t         genMetPF2PATPx;
   Float_t         genMetPF2PATPy;
   Int_t           numTauPF2PAT;
   Float_t         tauPF2PATE[1];   //[numTauPF2PAT]
   Float_t         tauPF2PATPt[1];   //[numTauPF2PAT]
   Float_t         tauPF2PATPhi[1];   //[numTauPF2PAT]
   Float_t         tauPF2PATEta[1];   //[numTauPF2PAT]
   Int_t           numGeneralTracks;
   Float_t         generalTracksPt[1000];   //[numGeneralTracks]
   Float_t         generalTracksEta[1000];   //[numGeneralTracks]
   Float_t         generalTracksTheta[1000];   //[numGeneralTracks]
   Float_t         generalTracksBeamSpotCorrectedD0[1000];   //[numGeneralTracks]
   Float_t         generalTracksPhi[1000];   //[numGeneralTracks]
   Int_t           generalTracksCharge[1000];   //[numGeneralTracks]
   Int_t           isElePlusJets;
   Float_t         genPDFScale;
   Float_t         genPDFx1;
   Float_t         genPDFx2;
   Int_t           genPDFf1;
   Int_t           genPDFf2;
   Int_t           processId;
   Float_t         processPtHat;
   Double_t        processMCWeight;
   Float_t         beamSpotX;
   Float_t         beamSpotY;
   Float_t         beamSpotZ;
   Float_t         pvX;
   Float_t         pvY;
   Float_t         pvZ;
   Float_t         pvDX;
   Float_t         pvDY;
   Float_t         pvDZ;
   Float_t         pvRho;
   Int_t           pvIsFake;
   Float_t         pvNdof;
   Float_t         pvChi2;
   Float_t         mhtPt;
   Float_t         mhtPy;
   Float_t         mhtPx;
   Float_t         mhtPhi;
   Float_t         mhtSumEt;
   Float_t         mhtSignif;
   Int_t           nTriggerBits;
   Int_t           TriggerBits[1];   //[nTriggerBits]
   Double_t        PileUpWeightRunA;
   Double_t        PileUpWeightRunB;
   Double_t        PileUpWeightRunC;
   //   Int_t           numVert;

   //2015 Data Triggers
   Int_t           HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2;
   Int_t           HLT_IsoMu20_v2;
   Int_t           HLT_IsoMu20_eta2p1_v2;
   Int_t           HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2;
   Int_t           HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2;
   Int_t           HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2;
   Int_t           HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v2;
   Int_t           HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3;
   Int_t           HLT_IsoMu18_v1;
   Int_t           HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3;
   Int_t           HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3;

   //MC Triggers
   Int_t           HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1;
   Int_t           HLT_IsoMu20_v1;
   Int_t           HLT_IsoMu20_eta2p1_v1;
   Int_t           HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1;
   Int_t           HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1;

   //Gen info
   Int_t           nGenPar;
   Float_t         genParEta[30];   //[nGenPar]
   Float_t         genParPhi[30];   //[nGenPar]
   Float_t         genParE[30];   //[nGenPar]
   Float_t         genParPt[30];   //[nGenPar]
   Int_t           genParId[30];   //[nGenPar]
   Int_t           genParCharge[30];   //[nGenPar]
   Int_t           eventRun;
   Int_t           eventNum;
   Float_t         eventLumiblock;

   // List of branches
   TBranch        *b_numElePF2PAT;   //!
   TBranch        *b_numLooseElePF2PAT;   //!
   TBranch        *b_elePF2PATE;   //!
   TBranch        *b_elePF2PATET;   //!
   TBranch        *b_elePF2PATPX;   //!
   TBranch        *b_elePF2PATPY;   //!
   TBranch        *b_elePF2PATPZ;   //!
   TBranch        *b_elePF2PATPhi;   //!
   TBranch        *b_elePF2PATTheta;   //!
   TBranch        *b_elePF2PATEta;   //!
   TBranch        *b_elePF2PATPT;   //!
   TBranch        *b_elePF2PATCharge;   //!
   TBranch        *b_elePF2PATMVA;   //!
   TBranch        *b_elePF2PATImpactTransDist;   //!
   TBranch        *b_elePF2PATImpactTransError;   //!
   TBranch        *b_elePF2PATImpactTransSignificance;   //!
   TBranch        *b_elePF2PATImpact3DDist;   //!
   TBranch        *b_elePF2PATImpact3DError;   //!
   TBranch        *b_elePF2PATImpact3DSignificance;   //!
   TBranch        *b_elePF2PATChargedHadronIso;   //!
   TBranch        *b_elePF2PATNeutralHadronIso;   //!
   TBranch        *b_elePF2PATPhotonIso;   //!
   TBranch        *b_elePF2PATTrackPt;   //!
   TBranch        *b_elePF2PATTrackPhi;   //!
   TBranch        *b_elePF2PATTrackEta;   //!
   TBranch        *b_elePF2PATTrackChi2;   //!
   TBranch        *b_elePF2PATTrackNDOF;   //!
   TBranch        *b_elePF2PATTrackD0;   //!
   TBranch        *b_elePF2PATTrackDBD0;   //!
   TBranch        *b_elePF2PATD0PV;   //!
   TBranch        *b_elePF2PATDZPV;   //!
   TBranch        *b_elePF2PATBeamSpotCorrectedTrackD0;   //!
   TBranch        *b_elePF2PATTrackDz;   //!
   TBranch        *b_elePF2PATVtxZ;   //!
   TBranch        *b_elePF2PATIsGsf;   //!
   TBranch        *b_elePF2PATGsfPx;   //!
   TBranch        *b_elePF2PATGsfPy;   //!
   TBranch        *b_elePF2PATGsfPz;   //!
   TBranch        *b_elePF2PATGsfE;   //!
   TBranch        *b_elePF2PATSCEta;   //!
   TBranch        *b_elePF2PATSCE;   //!
   TBranch        *b_elePF2PATSCPhi;   //!
   TBranch        *b_elePF2PATSCSigmaEtaEta;   //!
   TBranch        *b_elePF2PATSCSigmaIEtaIEta;   //!
   TBranch        *b_elePF2PATSCE1x5;   //!
   TBranch        *b_elePF2PATSCE5x5;   //!
   TBranch        *b_elePF2PATSCE2x5max;   //!
   TBranch        *b_elePF2PATTrackIso04;   //!
   TBranch        *b_elePF2PATEcalIso04;   //!
   TBranch        *b_elePF2PATHcalIso04;   //!
   TBranch        *b_elePF2PATTrackIso03;   //!
   TBranch        *b_elePF2PATEcalIso03;   //!
   TBranch        *b_elePF2PATHcalIso03;   //!
   TBranch        *b_elePF2PATdr04EcalRecHitSumEt;   //!
   TBranch        *b_elePF2PATdr03EcalRecHitSumEt;   //!
   TBranch        *b_elePF2PATEcalIsoDeposit;   //!
   TBranch        *b_elePF2PATHcalIsoDeposit;   //!
   TBranch        *b_elePF2PATComRelIso;   //!
   TBranch        *b_elePF2PATComRelIsodBeta;   //!
   TBranch        *b_elePF2PATComRelIsoRho;   //!
   TBranch        *b_elePF2PATChHadIso;   //!
   TBranch        *b_elePF2PATNtHadIso;   //!
   TBranch        *b_elePF2PATGammaIso;   //!
   TBranch        *b_elePF2PATRhoIso;   //!
   TBranch        *b_elePF2PATAEff03;   //!
   TBranch        *b_elePF2PATMissingInnerLayers;   //!
   TBranch        *b_elePF2PATHoverE;   //!
   TBranch        *b_elePF2PATDeltaPhiSC;   //!
   TBranch        *b_elePF2PATDeltaEtaSC;   //!
   TBranch        *b_elePF2PATIsBarrel;   //!
   TBranch        *b_elePF2PATPhotonConversionTag;   //!
   TBranch        *b_elePF2PATPhotonConversionDist;   //!
   TBranch        *b_elePF2PATPhotonConversionDcot;   //!
   TBranch        *b_elePF2PATPhotonConversionVeto;   //!
   TBranch        *b_elePF2PATPhotonConversionTagCustom;   //!
   TBranch        *b_elePF2PATPhotonConversionDistCustom;   //!
   TBranch        *b_elePF2PATPhotonConversionDcotCustom;   //!
   TBranch        *b_elePF2PATTriggerMatch;   //!
   TBranch        *b_elePF2PATJetOverlap;   //!
   TBranch        *b_genElePF2PATET;   //!
   TBranch        *b_genElePF2PATPX;   //!
   TBranch        *b_genElePF2PATPY;   //!
   TBranch        *b_genElePF2PATPZ;   //!
   TBranch        *b_genElePF2PATPhi;   //!
   TBranch        *b_genElePF2PATTheta;   //!
   TBranch        *b_genElePF2PATEta;   //!
   TBranch        *b_genElePF2PATCharge;   //!
   TBranch        *b_elePF2PATlooseElectronSortedEt;   //!
   TBranch        *b_elePF2PATlooseElectronSortedPt;   //!
   TBranch        *b_elePF2PATlooseElectronSortedEta;   //!
   TBranch        *b_elePF2PATlooseElectronSortedMVA;   //!
   TBranch        *b_elePF2PATlooseElectronSortedRelIso;   //!
   TBranch        *b_numMuonPF2PAT;   //!
   TBranch        *b_numLooseMuonPF2PAT;   //!
   TBranch        *b_muonPF2PATE;   //!
   TBranch        *b_muonPF2PATET;   //!
   TBranch        *b_muonPF2PATPt;   //!
   TBranch        *b_muonPF2PATPX;   //!
   TBranch        *b_muonPF2PATPY;   //!
   TBranch        *b_muonPF2PATPZ;   //!
   TBranch        *b_muonPF2PATPhi;   //!
   TBranch        *b_muonPF2PATTheta;   //!
   TBranch        *b_muonPF2PATEta;   //!
   TBranch        *b_muonPF2PATCharge;   //!
   TBranch        *b_muonPF2PATGlobalID;   //!
   TBranch        *b_muonPF2PATTrackID;   //!
   TBranch        *b_muonPF2PATChi2;   //!
   TBranch        *b_muonPF2PATD0;   //!
   TBranch        *b_muonPF2PATTrackDBD0;   //!
   TBranch        *b_muonPF2PATDBInnerTrackD0;   //!
   TBranch        *b_muonPF2PATBeamSpotCorrectedD0;   //!
   TBranch        *b_muonPF2PATTrackNHits;   //!
   TBranch        *b_muonPF2PATMuonNHits;   //!
   TBranch        *b_muonPF2PATNDOF;   //!
   TBranch        *b_muonPF2PATVertX;   //!
   TBranch        *b_muonPF2PATVertY;   //!
   TBranch        *b_muonPF2PATVertZ;   //!
   TBranch        *b_muonPF2PATChargedHadronIso;   //!
   TBranch        *b_muonPF2PATNeutralHadronIso;   //!
   TBranch        *b_muonPF2PATPhotonIso;   //!
   TBranch        *b_muonPF2PATTrackIso;   //!
   TBranch        *b_muonPF2PATEcalIso;   //!
   TBranch        *b_muonPF2PATHcalIso;   //!
   TBranch        *b_muonPF2PATComRelIso;   //!
   TBranch        *b_muonPF2PATComRelIsodBeta;   //!
   TBranch        *b_muonPF2PATIsPFMuon;   //!
   TBranch        *b_muonPF2PATNChambers;   //!
   TBranch        *b_muonPF2PATNMatches;   //!
   TBranch        *b_muonPF2PATTkLysWithMeasurements;   //!
   TBranch        *b_muonPF2PATVldPixHits;   //!
   TBranch        *b_muonPF2PATMatchedStations;   //!
   TBranch        *b_muonPF2PATGlbTkNormChi2;   //!
   TBranch        *b_muonPF2PATDBPV;   //!
   TBranch        *b_muonPF2PATDZPV;   //!
   TBranch        *b_muonPF2PATlooseMuonSortedEt;   //!
   TBranch        *b_muonPF2PATlooseMuonSortedPt;   //!
   TBranch        *b_muonPF2PATlooseMuonSortedEta;   //!
   TBranch        *b_muonPF2PATlooseMuonSortedRelIso;   //!
   TBranch        *b_muonPF2PATlooseMuonSortedisGlb;   //!
   TBranch        *b_muonPF2PATlooseMuonSortedisTrk;   //!
   TBranch        *b_genMuonPF2PATET;   //!
   TBranch        *b_genMuonPF2PATPX;   //!
   TBranch        *b_genMuonPF2PATPY;   //!
   TBranch        *b_genMuonPF2PATPZ;   //!
   TBranch        *b_genMuonPF2PATPhi;   //!
   TBranch        *b_genMuonPF2PATTheta;   //!
   TBranch        *b_genMuonPF2PATEta;   //!
   TBranch        *b_genMuonPF2PATCharge;   //!
   TBranch        *b_numJetPF2PAT;   //!
   TBranch        *b_jetPF2PATE;   //!
   TBranch        *b_jetPF2PATEt;   //!
   TBranch        *b_jetPF2PATPt;   //!
   TBranch        *b_jetPF2PATPtRaw;   //!
   TBranch        *b_jetPF2PATUnCorEt;   //!
   TBranch        *b_jetPF2PATUnCorPt;   //!
   TBranch        *b_jetPF2PATEta;   //!
   TBranch        *b_jetPF2PATTheta;   //!
   TBranch        *b_jetPF2PATPhi;   //!
   TBranch        *b_jetPF2PATPx;   //!
   TBranch        *b_jetPF2PATPy;   //!
   TBranch        *b_jetPF2PATPz;   //!
   TBranch        *b_jetPF2PATdRClosestLepton;   //!
   TBranch        *b_jetPF2PATNtracksInJet;   //!
   TBranch        *b_jetPF2PATJetCharge;   //!
   TBranch        *b_jetPF2PATfHPD;   //!
   TBranch        *b_jetPF2PATBtagSoftMuonPtRel;   //!
   TBranch        *b_jetPF2PATBtagSoftMuonQuality;   //!
   TBranch        *b_jetPF2PATCorrFactor;   //!
   TBranch        *b_jetPF2PATCorrResidual;   //!
   TBranch        *b_jetPF2PATL2L3ResErr;   //!
   TBranch        *b_jetPF2PATCorrErrLow;   //!
   TBranch        *b_jetPF2PATCorrErrHi;   //!
   TBranch        *b_jetPF2PATN90Hits;   //!
   TBranch        *b_jetPF2PATTriggered;   //!
   TBranch        *b_jetPF2PATSVPT;   //!
   TBranch        *b_jetPF2PATSVL2D;   //!
   TBranch        *b_jetPF2PATSVL2Dxy;   //!
   TBranch        *b_jetPF2PATSVL2DxyErr;   //!
   TBranch        *b_jetPF2PATSVL2DxySig;   //!
   TBranch        *b_jetPF2PATSVL3D;   //!
   TBranch        *b_jetPF2PATSVL3DErr;   //!
   TBranch        *b_jetPF2PATSVL3DSig;   //!
   TBranch        *b_jetPF2PATSVMass;   //!
   TBranch        *b_jetPF2PATSVNtracks;   //!
   TBranch        *b_jetPF2PATSVX;   //!
   TBranch        *b_jetPF2PATSVY;   //!
   TBranch        *b_jetPF2PATSVZ;   //!
   TBranch        *b_jetPF2PATSVDX;   //!
   TBranch        *b_jetPF2PATSVDY;   //!
   TBranch        *b_jetPF2PATSVDZ;   //!
   TBranch        *b_jetPF2PATBDiscriminator;   //!
   TBranch        *b_jetPF2PATNConstituents;   //!
   TBranch        *b_jetPF2PATPID;   //!
   TBranch        *b_jetPF2PATClosestBPartonDeltaR;   //!
   TBranch        *b_jetPF2PATClosestCPartonDeltaR;   //!
   TBranch        *b_genJetPF2PATET;   //!
   TBranch        *b_genJetPF2PATPT;   //!
   TBranch        *b_genJetPF2PATPX;   //!
   TBranch        *b_genJetPF2PATPY;   //!
   TBranch        *b_genJetPF2PATPZ;   //!
   TBranch        *b_genJetPF2PATPhi;   //!
   TBranch        *b_genJetPF2PATTheta;   //!
   TBranch        *b_genJetPF2PATEta;   //!
   TBranch        *b_genJetPF2PATPID;   //!
   TBranch        *b_numLooseBJetsPF2PAT;   //!
   TBranch        *b_looseJetPF2PATEt;   //!
   TBranch        *b_looseJetPF2PATPt;   //!
   TBranch        *b_looseJetPF2PATEta;   //!
   TBranch        *b_looseJetPF2PATBDisc;   //!
   TBranch        *b_jetPF2PATBtagDisc_trackCountingHighPurBJetTags;   //!
   TBranch        *b_jetPF2PATBtagDisc_trackCountingHighEffBJetTags;   //!
   TBranch        *b_jetPF2PATBtagDisc_jetProbabilityBJetTags;   //!
   TBranch        *b_jetPF2PATBtagDisc_jetBProbabilityBJetTags;   //!
   TBranch        *b_jetPF2PATBtagDisc_softElectronBJetTags;   //!
   TBranch        *b_jetPF2PATBtagDisc_softMuonBJetTags;   //!
   TBranch        *b_jetPF2PATBtagDisc_softMuonNoIPBJetTags;   //!
   TBranch        *b_jetPF2PATBtagDisc_simpleSecondaryVertexHighEffBJetTags;   //!
   TBranch        *b_jetPF2PATBtagDisc_simpleSecondaryVertexNegativeBJetTags;   //!
   TBranch        *b_jetPF2PATBtagDisc_combinedSecondaryVertexMVABJetTags;   //!
   TBranch        *b_jetPF2PATBtagDisc_simpleSecondaryVertexHighPurBJetTags;   //!
   TBranch        *b_jetPF2PATMuEnergy;   //!
   TBranch        *b_jetPF2PATMuEnergyFraction;   //!
   TBranch        *b_jetPF2PATNeutralHadEnergy;   //!
   TBranch        *b_jetPF2PATNeutralEmEnergy;   //!
   TBranch        *b_jetPF2PATChargedHadronEnergyFraction;   //!
   TBranch        *b_jetPF2PATNeutralHadronEnergyFraction;   //!
   TBranch        *b_jetPF2PATChargedEmEnergyFraction;   //!
   TBranch        *b_jetPF2PATNeutralEmEnergyFraction;   //!
   TBranch        *b_jetPF2PATChargedHadronEnergyFractionCorr;   //!
   TBranch        *b_jetPF2PATNeutralHadronEnergyFractionCorr;   //!
   TBranch        *b_jetPF2PATChargedEmEnergyFractionCorr;   //!
   TBranch        *b_jetPF2PATNeutralEmEnergyFractionCorr;   //!
   TBranch        *b_jetPF2PATNeutralMultiplicity;   //!
   TBranch        *b_jetPF2PATChargedMultiplicity;   //!
   TBranch        *b_metPF2PATEt;   //!
   TBranch        *b_metPF2PATEtRaw;   //!
   TBranch        *b_metPF2PATPhi;   //!
   TBranch        *b_metPF2PATPt;   //!
   TBranch        *b_metPF2PATPx;   //!
   TBranch        *b_metPF2PATPy;   //!
   TBranch        *b_metPF2PATScalarEt;   //!
   TBranch        *b_metPF2PATEtUncorrected;   //!
   TBranch        *b_metPF2PATPhiUncorrected;   //!
   TBranch        *b_genMetPF2PATEt;   //!
   TBranch        *b_genMetPF2PATPhi;   //!
   TBranch        *b_genMetPF2PATPt;   //!
   TBranch        *b_genMetPF2PATPx;   //!
   TBranch        *b_genMetPF2PATPy;   //!
   TBranch        *b_numTauPF2PAT;   //!
   TBranch        *b_tauPF2PATE;   //!
   TBranch        *b_tauPF2PATPt;   //!
   TBranch        *b_tauPF2PATPhi;   //!
   TBranch        *b_tauPF2PATEta;   //!
   TBranch        *b_numGeneralTracks;   //!
   TBranch        *b_generalTracksPt;   //!
   TBranch        *b_generalTracksEta;   //!
   TBranch        *b_generalTracksTheta;   //!
   TBranch        *b_generalTracksBeamSpotCorrectedD0;   //!
   TBranch        *b_generalTracksPhi;   //!
   TBranch        *b_generalTracksCharge;   //!
   TBranch        *b_isElePlusJets;   //!
   TBranch        *b_genPDFScale;   //!
   TBranch        *b_genPDFx1;   //!
   TBranch        *b_genPDFx2;   //!
   TBranch        *b_genPDFf1;   //!
   TBranch        *b_genPDFf2;   //!
   TBranch        *b_processId;   //!
   TBranch        *b_processPtHat;   //!
   TBranch        *b_processMCWeight;   //!
   TBranch        *b_beamSpotX;   //!
   TBranch        *b_beamSpotY;   //!
   TBranch        *b_beamSpotZ;   //!
   TBranch        *b_pvX;   //!
   TBranch        *b_pvY;   //!
   TBranch        *b_pvZ;   //!
   TBranch        *b_pvDX;   //!
   TBranch        *b_pvDY;   //!
   TBranch        *b_pvDZ;   //!
   TBranch        *b_pvRho;   //!
   TBranch        *b_pvIsFake;   //!
   TBranch        *b_pvNdof;   //!
   TBranch        *b_pvChi2;   //!
   TBranch        *b_mhtPt;   //!
   TBranch        *b_mhtPy;   //!
   TBranch        *b_mhtPx;   //!
   TBranch        *b_mhtPhi;   //!
   TBranch        *b_mhtSumEt;   //!
   TBranch        *b_mhtSignif;   //!
   TBranch        *b_nTriggerBits;   //!
   TBranch        *b_TriggerBits;   //!
   TBranch        *b_pileUpWeight;   //!
   //   TBranch        *b_numVert;    //!
   TBranch        *b_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2;
   TBranch        *b_HLT_IsoMu20_v2;
   TBranch        *b_HLT_IsoMu20_eta2p1_v2;
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2;
   TBranch        *b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2;
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2;
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v2;
   TBranch        *b_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3;
   TBranch        *b_HLT_IsoMu18_v1;
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3;
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3;
   TBranch        *b_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1;
   TBranch        *b_HLT_IsoMu20_v1;
   TBranch        *b_HLT_IsoMu20_eta2p1_v1;
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1;
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1;
   TBranch        *b_nGenPar;   //!
   TBranch        *b_genParEta;   //!
   TBranch        *b_genParPhi;   //!
   TBranch        *b_genParE;   //!
   TBranch        *b_genParPt;   //!
   TBranch        *b_genParId;   //!
   TBranch        *b_genParCharge;   //!
   TBranch        *b_eventRun;   //!
   TBranch        *b_eventNum;   //!
   TBranch        *b_eventLumiblock;   //!

   std::vector<int> electronIndexTight;
   std::vector<int> electronIndexLoose;
   std::vector<int> muonIndexTight;
   std::vector<int> muonIndexLoose;
   std::vector<int> jetIndex;
   std::vector<int> bTagIndex; 
   std::pair<TLorentzVector,TLorentzVector> zPairLeptons;
   std::pair<float,float> zPairRelIso;
   std::pair<int,int> zPairIndex; 
   TLorentzVector wLepton;
   int wLepIndex;
   float wLeptonRelIso;
   Int_t numVert;
   TBranch * b_numVert;

   AnalysisEvent(bool isMC = true, std::string triggerFlag = "", TTree *tree=0);
   virtual ~AnalysisEvent();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(bool isMC, std::string triggerFlag, TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   float getEventWeight(Long64_t entry);
};

#endif

#ifdef AnalysisEvent_cxx
AnalysisEvent::AnalysisEvent(bool isMC, std::string triggerFlag, TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_9_0_MJP_skim.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_9_0_MJP_skim.root");
      }
      f->GetObject("tree",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("tree","");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_100_0_Gu6_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_101_0_CCe_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_102_0_10I_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_103_0_OgH_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_104_0_Nzk_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_105_0_iJR_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_106_0_tN9_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_107_0_h2Y_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_108_0_hBb_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_109_0_UEL_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_10_0_Ib1_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_110_0_V1q_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_111_0_n74_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_112_0_hBx_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_113_0_WeN_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_114_0_KaM_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_115_0_aRw_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_116_0_3pE_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_117_0_gJ7_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_118_0_mdM_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_119_0_YT2_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_11_0_bm4_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_120_0_E3n_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_121_0_XdE_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_122_0_uZf_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_123_0_4GT_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_124_0_sDL_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_125_0_XCl_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_126_0_Hy1_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_127_0_b4J_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_128_0_o4T_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_129_0_N37_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_12_0_z2F_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_130_0_aU1_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_131_0_PdR_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_132_0_xFc_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_133_0_ne2_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_134_0_E9S_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_135_0_I03_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_136_0_IN2_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_137_0_Uqj_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_138_0_kas_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_139_0_Tcm_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_13_0_QjC_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_140_0_nU2_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_14_0_HCk_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_15_0_Tgc_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_16_0_0LI_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_17_1_Kxc_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_18_0_8Y8_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_19_0_EVu_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_1_1_e4c_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_20_0_pdV_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_21_1_hD4_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_22_0_KYU_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_23_0_yFZ_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_24_0_82B_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_25_4_iUG_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_26_0_pQ6_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_27_1_P51_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_28_0_K6r_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_29_0_Rkn_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_2_1_L0i_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_30_0_a9x_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_31_0_LG5_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_32_1_uRU_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_33_0_xP0_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_34_0_hCV_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_35_0_9n9_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_36_0_FAp_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_37_0_ST3_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_38_0_JeX_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_39_0_3aN_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_3_1_2ED_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_40_1_Dlo_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_41_0_kJO_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_42_0_NPy_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_43_0_pSx_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_44_0_efb_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_45_0_rj2_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_46_0_t0y_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_48_0_bFb_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_49_0_hFD_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_4_1_gY0_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_50_0_mOO_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_51_0_Tvz_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_52_0_Qyu_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_53_0_tnv_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_54_1_BWs_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_55_0_5I2_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_56_1_VTy_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_57_0_JtC_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_58_0_V1s_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_59_1_kJJ_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_5_0_ARX_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_60_0_d9h_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_61_1_4bi_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_62_0_Awp_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_63_0_NZJ_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_64_0_l4R_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_65_0_Xcp_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_66_0_msr_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_67_1_qQe_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_68_0_NNV_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_69_0_MmH_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_6_0_eIW_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_70_0_xdw_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_71_0_m7u_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_72_0_rIz_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_73_0_nnp_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_74_0_ai2_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_75_0_qBL_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_76_0_jS5_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_77_0_LlI_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_78_0_iIw_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_79_0_8BG_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_7_1_Tzw_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_80_0_sbp_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_81_0_kfW_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_82_0_OL9_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_83_0_4xd_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_84_0_92G_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_85_0_PbM_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_86_0_c6m_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_87_0_WGZ_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_88_0_8HI_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_89_0_Hjw_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_8_0_gZh_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_90_0_Aed_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_91_0_hOn_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_92_0_4zD_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_93_0_DHy_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_94_0_4l8_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_95_0_t8O_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_96_0_As6_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_97_0_zCu_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_98_0_mn7_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_99_0_Q1y_skim.root/tree");
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_9_0_MJP_skim.root/tree");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(isMC,triggerFlag,tree);
}

AnalysisEvent::~AnalysisEvent()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnalysisEvent::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AnalysisEvent::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void AnalysisEvent::Init(bool isMC, std::string triggerFlag, TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("numElePF2PAT", &numElePF2PAT, &b_numElePF2PAT);
   fChain->SetBranchAddress("numLooseElePF2PAT", &numLooseElePF2PAT, &b_numLooseElePF2PAT);
   fChain->SetBranchAddress("elePF2PATE", elePF2PATE, &b_elePF2PATE);
   fChain->SetBranchAddress("elePF2PATET", elePF2PATET, &b_elePF2PATET);
   fChain->SetBranchAddress("elePF2PATPX", elePF2PATPX, &b_elePF2PATPX);
   fChain->SetBranchAddress("elePF2PATPY", elePF2PATPY, &b_elePF2PATPY);
   fChain->SetBranchAddress("elePF2PATPZ", elePF2PATPZ, &b_elePF2PATPZ);
   fChain->SetBranchAddress("elePF2PATPhi", elePF2PATPhi, &b_elePF2PATPhi);
   fChain->SetBranchAddress("elePF2PATTheta", elePF2PATTheta, &b_elePF2PATTheta);
   fChain->SetBranchAddress("elePF2PATEta", elePF2PATEta, &b_elePF2PATEta);
   fChain->SetBranchAddress("elePF2PATPT", elePF2PATPT, &b_elePF2PATPT);
   fChain->SetBranchAddress("elePF2PATCharge", elePF2PATCharge, &b_elePF2PATCharge);
   fChain->SetBranchAddress("elePF2PATMVA", elePF2PATMVA, &b_elePF2PATMVA);
   fChain->SetBranchAddress("elePF2PATImpactTransDist", elePF2PATImpactTransDist, &b_elePF2PATImpactTransDist);
   fChain->SetBranchAddress("elePF2PATImpactTransError", elePF2PATImpactTransError, &b_elePF2PATImpactTransError);
   fChain->SetBranchAddress("elePF2PATImpactTransSignificance", elePF2PATImpactTransSignificance, &b_elePF2PATImpactTransSignificance);
   fChain->SetBranchAddress("elePF2PATImpact3DDist", elePF2PATImpact3DDist, &b_elePF2PATImpact3DDist);
   fChain->SetBranchAddress("elePF2PATImpact3DError", elePF2PATImpact3DError, &b_elePF2PATImpact3DError);
   fChain->SetBranchAddress("elePF2PATImpact3DSignificance", elePF2PATImpact3DSignificance, &b_elePF2PATImpact3DSignificance);
   fChain->SetBranchAddress("elePF2PATChargedHadronIso", elePF2PATChargedHadronIso, &b_elePF2PATChargedHadronIso);
   fChain->SetBranchAddress("elePF2PATNeutralHadronIso", elePF2PATNeutralHadronIso, &b_elePF2PATNeutralHadronIso);
   fChain->SetBranchAddress("elePF2PATPhotonIso", elePF2PATPhotonIso, &b_elePF2PATPhotonIso);
   fChain->SetBranchAddress("elePF2PATTrackPt", elePF2PATTrackPt, &b_elePF2PATTrackPt);
   fChain->SetBranchAddress("elePF2PATTrackPhi", elePF2PATTrackPhi, &b_elePF2PATTrackPhi);
   fChain->SetBranchAddress("elePF2PATTrackEta", elePF2PATTrackEta, &b_elePF2PATTrackEta);
   fChain->SetBranchAddress("elePF2PATTrackChi2", elePF2PATTrackChi2, &b_elePF2PATTrackChi2);
   fChain->SetBranchAddress("elePF2PATTrackNDOF", elePF2PATTrackNDOF, &b_elePF2PATTrackNDOF);
   fChain->SetBranchAddress("elePF2PATTrackD0", elePF2PATTrackD0, &b_elePF2PATTrackD0);
   fChain->SetBranchAddress("elePF2PATTrackDBD0", elePF2PATTrackDBD0, &b_elePF2PATTrackDBD0);
   fChain->SetBranchAddress("elePF2PATD0PV", elePF2PATD0PV, &b_elePF2PATD0PV);
   fChain->SetBranchAddress("elePF2PATDZPV", elePF2PATDZPV, &b_elePF2PATDZPV);
   fChain->SetBranchAddress("elePF2PATBeamSpotCorrectedTrackD0", elePF2PATBeamSpotCorrectedTrackD0, &b_elePF2PATBeamSpotCorrectedTrackD0);
   fChain->SetBranchAddress("elePF2PATTrackDz", elePF2PATTrackDz, &b_elePF2PATTrackDz);
   fChain->SetBranchAddress("elePF2PATVtxZ", elePF2PATVtxZ, &b_elePF2PATVtxZ);
   fChain->SetBranchAddress("elePF2PATIsGsf", elePF2PATIsGsf, &b_elePF2PATIsGsf);
   fChain->SetBranchAddress("elePF2PATGsfPx", elePF2PATGsfPx, &b_elePF2PATGsfPx);
   fChain->SetBranchAddress("elePF2PATGsfPy", elePF2PATGsfPy, &b_elePF2PATGsfPy);
   fChain->SetBranchAddress("elePF2PATGsfPz", elePF2PATGsfPz, &b_elePF2PATGsfPz);
   fChain->SetBranchAddress("elePF2PATGsfE", elePF2PATGsfE, &b_elePF2PATGsfE);
   fChain->SetBranchAddress("elePF2PATSCEta", elePF2PATSCEta, &b_elePF2PATSCEta);
   fChain->SetBranchAddress("elePF2PATSCE", elePF2PATSCE, &b_elePF2PATSCE);
   fChain->SetBranchAddress("elePF2PATSCPhi", elePF2PATSCPhi, &b_elePF2PATSCPhi);
   fChain->SetBranchAddress("elePF2PATSCSigmaEtaEta", elePF2PATSCSigmaEtaEta, &b_elePF2PATSCSigmaEtaEta);
   fChain->SetBranchAddress("elePF2PATSCSigmaIEtaIEta", elePF2PATSCSigmaIEtaIEta, &b_elePF2PATSCSigmaIEtaIEta);
   fChain->SetBranchAddress("elePF2PATSCE1x5", elePF2PATSCE1x5, &b_elePF2PATSCE1x5);
   fChain->SetBranchAddress("elePF2PATSCE5x5", elePF2PATSCE5x5, &b_elePF2PATSCE5x5);
   fChain->SetBranchAddress("elePF2PATSCE2x5max", elePF2PATSCE2x5max, &b_elePF2PATSCE2x5max);
   fChain->SetBranchAddress("elePF2PATTrackIso04", elePF2PATTrackIso04, &b_elePF2PATTrackIso04);
   fChain->SetBranchAddress("elePF2PATEcalIso04", elePF2PATEcalIso04, &b_elePF2PATEcalIso04);
   fChain->SetBranchAddress("elePF2PATHcalIso04", elePF2PATHcalIso04, &b_elePF2PATHcalIso04);
   fChain->SetBranchAddress("elePF2PATTrackIso03", elePF2PATTrackIso03, &b_elePF2PATTrackIso03);
   fChain->SetBranchAddress("elePF2PATEcalIso03", elePF2PATEcalIso03, &b_elePF2PATEcalIso03);
   fChain->SetBranchAddress("elePF2PATHcalIso03", elePF2PATHcalIso03, &b_elePF2PATHcalIso03);
   fChain->SetBranchAddress("elePF2PATdr04EcalRecHitSumEt", elePF2PATdr04EcalRecHitSumEt, &b_elePF2PATdr04EcalRecHitSumEt);
   fChain->SetBranchAddress("elePF2PATdr03EcalRecHitSumEt", elePF2PATdr03EcalRecHitSumEt, &b_elePF2PATdr03EcalRecHitSumEt);
   fChain->SetBranchAddress("elePF2PATEcalIsoDeposit", elePF2PATEcalIsoDeposit, &b_elePF2PATEcalIsoDeposit);
   fChain->SetBranchAddress("elePF2PATHcalIsoDeposit", elePF2PATHcalIsoDeposit, &b_elePF2PATHcalIsoDeposit);
   fChain->SetBranchAddress("elePF2PATComRelIso", elePF2PATComRelIso, &b_elePF2PATComRelIso);
   fChain->SetBranchAddress("elePF2PATComRelIsodBeta", elePF2PATComRelIsodBeta, &b_elePF2PATComRelIsodBeta);
   fChain->SetBranchAddress("elePF2PATComRelIsoRho", elePF2PATComRelIsoRho, &b_elePF2PATComRelIsoRho);
   fChain->SetBranchAddress("elePF2PATChHadIso", elePF2PATChHadIso, &b_elePF2PATChHadIso);
   fChain->SetBranchAddress("elePF2PATNtHadIso", elePF2PATNtHadIso, &b_elePF2PATNtHadIso);
   fChain->SetBranchAddress("elePF2PATGammaIso", elePF2PATGammaIso, &b_elePF2PATGammaIso);
   fChain->SetBranchAddress("elePF2PATRhoIso", elePF2PATRhoIso, &b_elePF2PATRhoIso);
   fChain->SetBranchAddress("elePF2PATAEff03", elePF2PATAEff03, &b_elePF2PATAEff03);
   fChain->SetBranchAddress("elePF2PATMissingInnerLayers", elePF2PATMissingInnerLayers, &b_elePF2PATMissingInnerLayers);
   fChain->SetBranchAddress("elePF2PATHoverE", elePF2PATHoverE, &b_elePF2PATHoverE);
   fChain->SetBranchAddress("elePF2PATDeltaPhiSC", elePF2PATDeltaPhiSC, &b_elePF2PATDeltaPhiSC);
   fChain->SetBranchAddress("elePF2PATDeltaEtaSC", elePF2PATDeltaEtaSC, &b_elePF2PATDeltaEtaSC);
   fChain->SetBranchAddress("elePF2PATIsBarrel", elePF2PATIsBarrel, &b_elePF2PATIsBarrel);
   fChain->SetBranchAddress("elePF2PATPhotonConversionTag", elePF2PATPhotonConversionTag, &b_elePF2PATPhotonConversionTag);
   fChain->SetBranchAddress("elePF2PATPhotonConversionDist", elePF2PATPhotonConversionDist, &b_elePF2PATPhotonConversionDist);
   fChain->SetBranchAddress("elePF2PATPhotonConversionDcot", elePF2PATPhotonConversionDcot, &b_elePF2PATPhotonConversionDcot);
   fChain->SetBranchAddress("elePF2PATPhotonConversionVeto", elePF2PATPhotonConversionVeto, &b_elePF2PATPhotonConversionVeto);
   fChain->SetBranchAddress("elePF2PATPhotonConversionTagCustom", elePF2PATPhotonConversionTagCustom, &b_elePF2PATPhotonConversionTagCustom);
   fChain->SetBranchAddress("elePF2PATPhotonConversionDistCustom", elePF2PATPhotonConversionDistCustom, &b_elePF2PATPhotonConversionDistCustom);
   fChain->SetBranchAddress("elePF2PATPhotonConversionDcotCustom", elePF2PATPhotonConversionDcotCustom, &b_elePF2PATPhotonConversionDcotCustom);
   fChain->SetBranchAddress("elePF2PATTriggerMatch", elePF2PATTriggerMatch, &b_elePF2PATTriggerMatch);
   fChain->SetBranchAddress("elePF2PATJetOverlap", elePF2PATJetOverlap, &b_elePF2PATJetOverlap);
   if (isMC) { 
     fChain->SetBranchAddress("genElePF2PATET", genElePF2PATET, &b_genElePF2PATET);
     fChain->SetBranchAddress("genElePF2PATPX", genElePF2PATPX, &b_genElePF2PATPX);
     fChain->SetBranchAddress("genElePF2PATPY", genElePF2PATPY, &b_genElePF2PATPY);
     fChain->SetBranchAddress("genElePF2PATPZ", genElePF2PATPZ, &b_genElePF2PATPZ);
     fChain->SetBranchAddress("genElePF2PATPhi", genElePF2PATPhi, &b_genElePF2PATPhi);
     fChain->SetBranchAddress("genElePF2PATTheta", genElePF2PATTheta, &b_genElePF2PATTheta);
     fChain->SetBranchAddress("genElePF2PATEta", genElePF2PATEta, &b_genElePF2PATEta);
     fChain->SetBranchAddress("genElePF2PATCharge", genElePF2PATCharge, &b_genElePF2PATCharge);
   }
   fChain->SetBranchAddress("elePF2PATlooseElectronSortedEt", elePF2PATlooseElectronSortedEt, &b_elePF2PATlooseElectronSortedEt);
   fChain->SetBranchAddress("elePF2PATlooseElectronSortedPt", elePF2PATlooseElectronSortedPt, &b_elePF2PATlooseElectronSortedPt);
   fChain->SetBranchAddress("elePF2PATlooseElectronSortedEta", elePF2PATlooseElectronSortedEta, &b_elePF2PATlooseElectronSortedEta);
   fChain->SetBranchAddress("elePF2PATlooseElectronSortedMVA", elePF2PATlooseElectronSortedMVA, &b_elePF2PATlooseElectronSortedMVA);
   fChain->SetBranchAddress("elePF2PATlooseElectronSortedRelIso", elePF2PATlooseElectronSortedRelIso, &b_elePF2PATlooseElectronSortedRelIso);
   fChain->SetBranchAddress("numMuonPF2PAT", &numMuonPF2PAT, &b_numMuonPF2PAT);
   fChain->SetBranchAddress("numLooseMuonPF2PAT", &numLooseMuonPF2PAT, &b_numLooseMuonPF2PAT);
   fChain->SetBranchAddress("muonPF2PATE", muonPF2PATE, &b_muonPF2PATE);
   fChain->SetBranchAddress("muonPF2PATET", muonPF2PATET, &b_muonPF2PATET);
   fChain->SetBranchAddress("muonPF2PATPt", muonPF2PATPt, &b_muonPF2PATPt);
   fChain->SetBranchAddress("muonPF2PATPX", muonPF2PATPX, &b_muonPF2PATPX);
   fChain->SetBranchAddress("muonPF2PATPY", muonPF2PATPY, &b_muonPF2PATPY);
   fChain->SetBranchAddress("muonPF2PATPZ", muonPF2PATPZ, &b_muonPF2PATPZ);
   fChain->SetBranchAddress("muonPF2PATPhi", muonPF2PATPhi, &b_muonPF2PATPhi);
   fChain->SetBranchAddress("muonPF2PATTheta", muonPF2PATTheta, &b_muonPF2PATTheta);
   fChain->SetBranchAddress("muonPF2PATEta", muonPF2PATEta, &b_muonPF2PATEta);
   fChain->SetBranchAddress("muonPF2PATCharge", muonPF2PATCharge, &b_muonPF2PATCharge);
   fChain->SetBranchAddress("muonPF2PATGlobalID", muonPF2PATGlobalID, &b_muonPF2PATGlobalID);
   fChain->SetBranchAddress("muonPF2PATTrackID", muonPF2PATTrackID, &b_muonPF2PATTrackID);
   fChain->SetBranchAddress("muonPF2PATChi2", muonPF2PATChi2, &b_muonPF2PATChi2);
   fChain->SetBranchAddress("muonPF2PATD0", muonPF2PATD0, &b_muonPF2PATD0);
   fChain->SetBranchAddress("muonPF2PATTrackDBD0", muonPF2PATTrackDBD0, &b_muonPF2PATTrackDBD0);
   fChain->SetBranchAddress("muonPF2PATDBInnerTrackD0", muonPF2PATDBInnerTrackD0, &b_muonPF2PATDBInnerTrackD0);
   fChain->SetBranchAddress("muonPF2PATBeamSpotCorrectedD0", muonPF2PATBeamSpotCorrectedD0, &b_muonPF2PATBeamSpotCorrectedD0);
   fChain->SetBranchAddress("muonPF2PATTrackNHits", muonPF2PATTrackNHits, &b_muonPF2PATTrackNHits);
   fChain->SetBranchAddress("muonPF2PATMuonNHits", muonPF2PATMuonNHits, &b_muonPF2PATMuonNHits);
   fChain->SetBranchAddress("muonPF2PATNDOF", muonPF2PATNDOF, &b_muonPF2PATNDOF);
   fChain->SetBranchAddress("muonPF2PATVertX", muonPF2PATVertX, &b_muonPF2PATVertX);
   fChain->SetBranchAddress("muonPF2PATVertY", muonPF2PATVertY, &b_muonPF2PATVertY);
   fChain->SetBranchAddress("muonPF2PATVertZ", muonPF2PATVertZ, &b_muonPF2PATVertZ);
   fChain->SetBranchAddress("muonPF2PATChargedHadronIso", muonPF2PATChargedHadronIso, &b_muonPF2PATChargedHadronIso);
   fChain->SetBranchAddress("muonPF2PATNeutralHadronIso", muonPF2PATNeutralHadronIso, &b_muonPF2PATNeutralHadronIso);
   fChain->SetBranchAddress("muonPF2PATPhotonIso", muonPF2PATPhotonIso, &b_muonPF2PATPhotonIso);
   fChain->SetBranchAddress("muonPF2PATTrackIso", muonPF2PATTrackIso, &b_muonPF2PATTrackIso);
   fChain->SetBranchAddress("muonPF2PATEcalIso", muonPF2PATEcalIso, &b_muonPF2PATEcalIso);
   fChain->SetBranchAddress("muonPF2PATHcalIso", muonPF2PATHcalIso, &b_muonPF2PATHcalIso);
   fChain->SetBranchAddress("muonPF2PATComRelIso", muonPF2PATComRelIso, &b_muonPF2PATComRelIso);
   fChain->SetBranchAddress("muonPF2PATComRelIsodBeta", muonPF2PATComRelIsodBeta, &b_muonPF2PATComRelIsodBeta);
   fChain->SetBranchAddress("muonPF2PATIsPFMuon", muonPF2PATIsPFMuon, &b_muonPF2PATIsPFMuon);
   fChain->SetBranchAddress("muonPF2PATNChambers", muonPF2PATNChambers, &b_muonPF2PATNChambers);
   fChain->SetBranchAddress("muonPF2PATNMatches", muonPF2PATNMatches, &b_muonPF2PATNMatches);
   fChain->SetBranchAddress("muonPF2PATTkLysWithMeasurements", muonPF2PATTkLysWithMeasurements, &b_muonPF2PATTkLysWithMeasurements);
   fChain->SetBranchAddress("muonPF2PATGlbTkNormChi2", muonPF2PATGlbTkNormChi2, &b_muonPF2PATGlbTkNormChi2);
   fChain->SetBranchAddress("muonPF2PATDBPV", muonPF2PATDBPV, &b_muonPF2PATDBPV);
   fChain->SetBranchAddress("muonPF2PATDZPV", muonPF2PATDZPV, &b_muonPF2PATDZPV);
   fChain->SetBranchAddress("muonPF2PATVldPixHits", muonPF2PATVldPixHits, &b_muonPF2PATVldPixHits);
   fChain->SetBranchAddress("muonPF2PATMatchedStations", muonPF2PATMatchedStations, &b_muonPF2PATMatchedStations);
   fChain->SetBranchAddress("muonPF2PATlooseMuonSortedEt", muonPF2PATlooseMuonSortedEt, &b_muonPF2PATlooseMuonSortedEt);
   fChain->SetBranchAddress("muonPF2PATlooseMuonSortedPt", muonPF2PATlooseMuonSortedPt, &b_muonPF2PATlooseMuonSortedPt);
   fChain->SetBranchAddress("muonPF2PATlooseMuonSortedEta", muonPF2PATlooseMuonSortedEta, &b_muonPF2PATlooseMuonSortedEta);
   fChain->SetBranchAddress("muonPF2PATlooseMuonSortedRelIso", muonPF2PATlooseMuonSortedRelIso, &b_muonPF2PATlooseMuonSortedRelIso);
   fChain->SetBranchAddress("muonPF2PATlooseMuonSortedisGlb", muonPF2PATlooseMuonSortedisGlb, &b_muonPF2PATlooseMuonSortedisGlb);
   fChain->SetBranchAddress("muonPF2PATlooseMuonSortedisTrk", muonPF2PATlooseMuonSortedisTrk, &b_muonPF2PATlooseMuonSortedisTrk);
   if (isMC) {
     fChain->SetBranchAddress("genMuonPF2PATET", genMuonPF2PATET, &b_genMuonPF2PATET);
     fChain->SetBranchAddress("genMuonPF2PATPX", genMuonPF2PATPX, &b_genMuonPF2PATPX);
     fChain->SetBranchAddress("genMuonPF2PATPY", genMuonPF2PATPY, &b_genMuonPF2PATPY);
     fChain->SetBranchAddress("genMuonPF2PATPZ", genMuonPF2PATPZ, &b_genMuonPF2PATPZ);
     fChain->SetBranchAddress("genMuonPF2PATPhi", genMuonPF2PATPhi, &b_genMuonPF2PATPhi);
     fChain->SetBranchAddress("genMuonPF2PATTheta", genMuonPF2PATTheta, &b_genMuonPF2PATTheta);
     fChain->SetBranchAddress("genMuonPF2PATEta", genMuonPF2PATEta, &b_genMuonPF2PATEta);
     fChain->SetBranchAddress("genMuonPF2PATCharge", genMuonPF2PATCharge, &b_genMuonPF2PATCharge);
   }
   fChain->SetBranchAddress("numJetPF2PAT", &numJetPF2PAT, &b_numJetPF2PAT);
   fChain->SetBranchAddress("jetPF2PATE", jetPF2PATE, &b_jetPF2PATE);
   fChain->SetBranchAddress("jetPF2PATEt", jetPF2PATEt, &b_jetPF2PATEt);
   fChain->SetBranchAddress("jetPF2PATPt", jetPF2PATPt, &b_jetPF2PATPt);
   fChain->SetBranchAddress("jetPF2PATPtRaw", jetPF2PATPtRaw, &b_jetPF2PATPtRaw);
   fChain->SetBranchAddress("jetPF2PATUnCorEt", jetPF2PATUnCorEt, &b_jetPF2PATUnCorEt);
   fChain->SetBranchAddress("jetPF2PATUnCorPt", jetPF2PATUnCorPt, &b_jetPF2PATUnCorPt);
   fChain->SetBranchAddress("jetPF2PATEta", jetPF2PATEta, &b_jetPF2PATEta);
   fChain->SetBranchAddress("jetPF2PATTheta", jetPF2PATTheta, &b_jetPF2PATTheta);
   fChain->SetBranchAddress("jetPF2PATPhi", jetPF2PATPhi, &b_jetPF2PATPhi);
   fChain->SetBranchAddress("jetPF2PATPx", jetPF2PATPx, &b_jetPF2PATPx);
   fChain->SetBranchAddress("jetPF2PATPy", jetPF2PATPy, &b_jetPF2PATPy);
   fChain->SetBranchAddress("jetPF2PATPz", jetPF2PATPz, &b_jetPF2PATPz);
   fChain->SetBranchAddress("jetPF2PATdRClosestLepton", jetPF2PATdRClosestLepton, &b_jetPF2PATdRClosestLepton);
   fChain->SetBranchAddress("jetPF2PATNtracksInJet", jetPF2PATNtracksInJet, &b_jetPF2PATNtracksInJet);
   fChain->SetBranchAddress("jetPF2PATJetCharge", jetPF2PATJetCharge, &b_jetPF2PATJetCharge);
   fChain->SetBranchAddress("jetPF2PATfHPD", jetPF2PATfHPD, &b_jetPF2PATfHPD);
   fChain->SetBranchAddress("jetPF2PATBtagSoftMuonPtRel", jetPF2PATBtagSoftMuonPtRel, &b_jetPF2PATBtagSoftMuonPtRel);
   fChain->SetBranchAddress("jetPF2PATBtagSoftMuonQuality", jetPF2PATBtagSoftMuonQuality, &b_jetPF2PATBtagSoftMuonQuality);
   fChain->SetBranchAddress("jetPF2PATCorrFactor", jetPF2PATCorrFactor, &b_jetPF2PATCorrFactor);
   fChain->SetBranchAddress("jetPF2PATCorrResidual", jetPF2PATCorrResidual, &b_jetPF2PATCorrResidual);
   fChain->SetBranchAddress("jetPF2PATL2L3ResErr", jetPF2PATL2L3ResErr, &b_jetPF2PATL2L3ResErr);
   fChain->SetBranchAddress("jetPF2PATCorrErrLow", jetPF2PATCorrErrLow, &b_jetPF2PATCorrErrLow);
   fChain->SetBranchAddress("jetPF2PATCorrErrHi", jetPF2PATCorrErrHi, &b_jetPF2PATCorrErrHi);
   fChain->SetBranchAddress("jetPF2PATN90Hits", jetPF2PATN90Hits, &b_jetPF2PATN90Hits);
   fChain->SetBranchAddress("jetPF2PATTriggered", jetPF2PATTriggered, &b_jetPF2PATTriggered);
   fChain->SetBranchAddress("jetPF2PATSVPT", jetPF2PATSVPT, &b_jetPF2PATSVPT);
   fChain->SetBranchAddress("jetPF2PATSVL2D", jetPF2PATSVL2D, &b_jetPF2PATSVL2D);
   fChain->SetBranchAddress("jetPF2PATSVL2Dxy", jetPF2PATSVL2Dxy, &b_jetPF2PATSVL2Dxy);
   fChain->SetBranchAddress("jetPF2PATSVL2DxyErr", jetPF2PATSVL2DxyErr, &b_jetPF2PATSVL2DxyErr);
   fChain->SetBranchAddress("jetPF2PATSVL2DxySig", jetPF2PATSVL2DxySig, &b_jetPF2PATSVL2DxySig);
   fChain->SetBranchAddress("jetPF2PATSVL3D", jetPF2PATSVL3D, &b_jetPF2PATSVL3D);
   fChain->SetBranchAddress("jetPF2PATSVL3DErr", jetPF2PATSVL3DErr, &b_jetPF2PATSVL3DErr);
   fChain->SetBranchAddress("jetPF2PATSVL3DSig", jetPF2PATSVL3DSig, &b_jetPF2PATSVL3DSig);
   fChain->SetBranchAddress("jetPF2PATSVMass", jetPF2PATSVMass, &b_jetPF2PATSVMass);
   fChain->SetBranchAddress("jetPF2PATSVNtracks", jetPF2PATSVNtracks, &b_jetPF2PATSVNtracks);
   fChain->SetBranchAddress("jetPF2PATSVX", jetPF2PATSVX, &b_jetPF2PATSVX);
   fChain->SetBranchAddress("jetPF2PATSVY", jetPF2PATSVY, &b_jetPF2PATSVY);
   fChain->SetBranchAddress("jetPF2PATSVZ", jetPF2PATSVZ, &b_jetPF2PATSVZ);
   fChain->SetBranchAddress("jetPF2PATSVDX", jetPF2PATSVDX, &b_jetPF2PATSVDX);
   fChain->SetBranchAddress("jetPF2PATSVDY", jetPF2PATSVDY, &b_jetPF2PATSVDY);
   fChain->SetBranchAddress("jetPF2PATSVDZ", jetPF2PATSVDZ, &b_jetPF2PATSVDZ);
   fChain->SetBranchAddress("jetPF2PATBDiscriminator", jetPF2PATBDiscriminator, &b_jetPF2PATBDiscriminator);
   fChain->SetBranchAddress("jetPF2PATNConstituents", jetPF2PATNConstituents, &b_jetPF2PATNConstituents);
   fChain->SetBranchAddress("jetPF2PATPID", jetPF2PATPID, &b_jetPF2PATPID);
   fChain->SetBranchAddress("jetPF2PATClosestBPartonDeltaR", jetPF2PATClosestBPartonDeltaR, &b_jetPF2PATClosestBPartonDeltaR);
   fChain->SetBranchAddress("jetPF2PATClosestCPartonDeltaR", jetPF2PATClosestCPartonDeltaR, &b_jetPF2PATClosestCPartonDeltaR);
   if (isMC) {
     fChain->SetBranchAddress("genJetPF2PATET", genJetPF2PATET, &b_genJetPF2PATET);
     fChain->SetBranchAddress("genJetPF2PATPT", genJetPF2PATPT, &b_genJetPF2PATPT);
     fChain->SetBranchAddress("genJetPF2PATPX", genJetPF2PATPX, &b_genJetPF2PATPX);
     fChain->SetBranchAddress("genJetPF2PATPY", genJetPF2PATPY, &b_genJetPF2PATPY);
     fChain->SetBranchAddress("genJetPF2PATPZ", genJetPF2PATPZ, &b_genJetPF2PATPZ);
     fChain->SetBranchAddress("genJetPF2PATPhi", genJetPF2PATPhi, &b_genJetPF2PATPhi);
     fChain->SetBranchAddress("genJetPF2PATTheta", genJetPF2PATTheta, &b_genJetPF2PATTheta);
     fChain->SetBranchAddress("genJetPF2PATEta", genJetPF2PATEta, &b_genJetPF2PATEta);
     fChain->SetBranchAddress("genJetPF2PATPID", genJetPF2PATPID, &b_genJetPF2PATPID);
   }
   fChain->SetBranchAddress("numLooseJetPF2PAT", &numLooseJetPF2PAT, &b_numLooseBJetsPF2PAT);
   fChain->SetBranchAddress("looseJetPF2PATEt", looseJetPF2PATEt, &b_looseJetPF2PATEt);
   fChain->SetBranchAddress("looseJetPF2PATPt", looseJetPF2PATPt, &b_looseJetPF2PATPt);
   fChain->SetBranchAddress("looseJetPF2PATEta", looseJetPF2PATEta, &b_looseJetPF2PATEta);
   fChain->SetBranchAddress("looseJetPF2PATBDisc", looseJetPF2PATBDisc, &b_looseJetPF2PATBDisc);
   fChain->SetBranchAddress("jetPF2PATBtagDisc_trackCountingHighPurBJetTags", jetPF2PATBtagDisc_trackCountingHighPurBJetTags, &b_jetPF2PATBtagDisc_trackCountingHighPurBJetTags);
   fChain->SetBranchAddress("jetPF2PATBtagDisc_trackCountingHighEffBJetTags", jetPF2PATBtagDisc_trackCountingHighEffBJetTags, &b_jetPF2PATBtagDisc_trackCountingHighEffBJetTags);
   fChain->SetBranchAddress("jetPF2PATBtagDisc_jetProbabilityBJetTags", jetPF2PATBtagDisc_jetProbabilityBJetTags, &b_jetPF2PATBtagDisc_jetProbabilityBJetTags);
   fChain->SetBranchAddress("jetPF2PATBtagDisc_jetBProbabilityBJetTags", jetPF2PATBtagDisc_jetBProbabilityBJetTags, &b_jetPF2PATBtagDisc_jetBProbabilityBJetTags);
   fChain->SetBranchAddress("jetPF2PATBtagDisc_softElectronBJetTags", jetPF2PATBtagDisc_softElectronBJetTags, &b_jetPF2PATBtagDisc_softElectronBJetTags);
   fChain->SetBranchAddress("jetPF2PATBtagDisc_softMuonBJetTags", jetPF2PATBtagDisc_softMuonBJetTags, &b_jetPF2PATBtagDisc_softMuonBJetTags);
   fChain->SetBranchAddress("jetPF2PATBtagDisc_softMuonNoIPBJetTags", jetPF2PATBtagDisc_softMuonNoIPBJetTags, &b_jetPF2PATBtagDisc_softMuonNoIPBJetTags);
   fChain->SetBranchAddress("jetPF2PATBtagDisc_simpleSecondaryVertexHighEffBJetTags", jetPF2PATBtagDisc_simpleSecondaryVertexHighEffBJetTags, &b_jetPF2PATBtagDisc_simpleSecondaryVertexHighEffBJetTags);
   fChain->SetBranchAddress("jetPF2PATBtagDisc_simpleSecondaryVertexNegativeBJetTags", jetPF2PATBtagDisc_simpleSecondaryVertexNegativeBJetTags, &b_jetPF2PATBtagDisc_simpleSecondaryVertexNegativeBJetTags);
   fChain->SetBranchAddress("jetPF2PATBtagDisc_combinedSecondaryVertexMVABJetTags", jetPF2PATBtagDisc_combinedSecondaryVertexMVABJetTags, &b_jetPF2PATBtagDisc_combinedSecondaryVertexMVABJetTags);
   fChain->SetBranchAddress("jetPF2PATBtagDisc_simpleSecondaryVertexHighPurBJetTags", jetPF2PATBtagDisc_simpleSecondaryVertexHighPurBJetTags, &b_jetPF2PATBtagDisc_simpleSecondaryVertexHighPurBJetTags);
   fChain->SetBranchAddress("jetPF2PATMuEnergy", jetPF2PATMuEnergy, &b_jetPF2PATMuEnergy);
   fChain->SetBranchAddress("jetPF2PATMuEnergyFraction", jetPF2PATMuEnergyFraction, &b_jetPF2PATMuEnergyFraction);
   fChain->SetBranchAddress("jetPF2PATNeutralHadEnergy", jetPF2PATNeutralHadEnergy, &b_jetPF2PATNeutralHadEnergy);
   fChain->SetBranchAddress("jetPF2PATNeutralEmEnergy", jetPF2PATNeutralEmEnergy, &b_jetPF2PATNeutralEmEnergy);
   fChain->SetBranchAddress("jetPF2PATChargedHadronEnergyFraction", jetPF2PATChargedHadronEnergyFraction, &b_jetPF2PATChargedHadronEnergyFraction);
   fChain->SetBranchAddress("jetPF2PATNeutralHadronEnergyFraction", jetPF2PATNeutralHadronEnergyFraction, &b_jetPF2PATNeutralHadronEnergyFraction);
   fChain->SetBranchAddress("jetPF2PATChargedEmEnergyFraction", jetPF2PATChargedEmEnergyFraction, &b_jetPF2PATChargedEmEnergyFraction);
   fChain->SetBranchAddress("jetPF2PATNeutralEmEnergyFraction", jetPF2PATNeutralEmEnergyFraction, &b_jetPF2PATNeutralEmEnergyFraction);
   fChain->SetBranchAddress("jetPF2PATChargedHadronEnergyFractionCorr", jetPF2PATChargedHadronEnergyFractionCorr, &b_jetPF2PATChargedHadronEnergyFractionCorr);
   fChain->SetBranchAddress("jetPF2PATNeutralHadronEnergyFractionCorr", jetPF2PATNeutralHadronEnergyFractionCorr, &b_jetPF2PATNeutralHadronEnergyFractionCorr);
   fChain->SetBranchAddress("jetPF2PATChargedEmEnergyFractionCorr", jetPF2PATChargedEmEnergyFractionCorr, &b_jetPF2PATChargedEmEnergyFractionCorr);
   fChain->SetBranchAddress("jetPF2PATNeutralEmEnergyFractionCorr", jetPF2PATNeutralEmEnergyFractionCorr, &b_jetPF2PATNeutralEmEnergyFractionCorr);
   fChain->SetBranchAddress("jetPF2PATNeutralMultiplicity", jetPF2PATNeutralMultiplicity, &b_jetPF2PATNeutralMultiplicity);
   fChain->SetBranchAddress("jetPF2PATChargedMultiplicity", jetPF2PATChargedMultiplicity, &b_jetPF2PATChargedMultiplicity);
   fChain->SetBranchAddress("metPF2PATEt", &metPF2PATEt, &b_metPF2PATEt);
   fChain->SetBranchAddress("metPF2PATEtRaw", &metPF2PATEtRaw, &b_metPF2PATEtRaw);
   fChain->SetBranchAddress("metPF2PATPhi", &metPF2PATPhi, &b_metPF2PATPhi);
   fChain->SetBranchAddress("metPF2PATPt", &metPF2PATPt, &b_metPF2PATPt);
   fChain->SetBranchAddress("metPF2PATPx", &metPF2PATPx, &b_metPF2PATPx);
   fChain->SetBranchAddress("metPF2PATPy", &metPF2PATPy, &b_metPF2PATPy);
   fChain->SetBranchAddress("metPF2PATScalarEt", &metPF2PATScalarEt, &b_metPF2PATScalarEt);
   fChain->SetBranchAddress("metPF2PATEtUncorrected", &metPF2PATEtUncorrected, &b_metPF2PATEtUncorrected);
   fChain->SetBranchAddress("metPF2PATPhiUncorrected", &metPF2PATPhiUncorrected, &b_metPF2PATPhiUncorrected);
   if (isMC) {
     fChain->SetBranchAddress("genMetPF2PATEt", &genMetPF2PATEt, &b_genMetPF2PATEt);
     fChain->SetBranchAddress("genMetPF2PATPhi", &genMetPF2PATPhi, &b_genMetPF2PATPhi);
     fChain->SetBranchAddress("genMetPF2PATPt", &genMetPF2PATPt, &b_genMetPF2PATPt);
     fChain->SetBranchAddress("genMetPF2PATPx", &genMetPF2PATPx, &b_genMetPF2PATPx);
     fChain->SetBranchAddress("genMetPF2PATPy", &genMetPF2PATPy, &b_genMetPF2PATPy);
   }
   fChain->SetBranchAddress("numTauPF2PAT", &numTauPF2PAT, &b_numTauPF2PAT);
   fChain->SetBranchAddress("tauPF2PATE", &tauPF2PATE, &b_tauPF2PATE);
   fChain->SetBranchAddress("tauPF2PATPt", &tauPF2PATPt, &b_tauPF2PATPt);
   fChain->SetBranchAddress("tauPF2PATPhi", &tauPF2PATPhi, &b_tauPF2PATPhi);
   fChain->SetBranchAddress("tauPF2PATEta", &tauPF2PATEta, &b_tauPF2PATEta);
   fChain->SetBranchAddress("numGeneralTracks", &numGeneralTracks, &b_numGeneralTracks);
   fChain->SetBranchAddress("generalTracksPt", generalTracksPt, &b_generalTracksPt);
   fChain->SetBranchAddress("generalTracksEta", generalTracksEta, &b_generalTracksEta);
   fChain->SetBranchAddress("generalTracksTheta", generalTracksTheta, &b_generalTracksTheta);
   fChain->SetBranchAddress("generalTracksBeamSpotCorrectedD0", generalTracksBeamSpotCorrectedD0, &b_generalTracksBeamSpotCorrectedD0);
   fChain->SetBranchAddress("generalTracksPhi", generalTracksPhi, &b_generalTracksPhi);
   fChain->SetBranchAddress("generalTracksCharge", generalTracksCharge, &b_generalTracksCharge);
   if (isMC) {
     fChain->SetBranchAddress("isElePlusJets", &isElePlusJets, &b_isElePlusJets);
     fChain->SetBranchAddress("genPDFScale", &genPDFScale, &b_genPDFScale);
     fChain->SetBranchAddress("genPDFx1", &genPDFx1, &b_genPDFx1);
     fChain->SetBranchAddress("genPDFx2", &genPDFx2, &b_genPDFx2);
     fChain->SetBranchAddress("genPDFf1", &genPDFf1, &b_genPDFf1);
     fChain->SetBranchAddress("genPDFf2", &genPDFf2, &b_genPDFf2);
   }
   fChain->SetBranchAddress("processId", &processId, &b_processId);
   fChain->SetBranchAddress("processPtHat", &processPtHat, &b_processPtHat);
   fChain->SetBranchAddress("processMCWeight", &processMCWeight, &b_processMCWeight);
   fChain->SetBranchAddress("beamSpotX", &beamSpotX, &b_beamSpotX);
   fChain->SetBranchAddress("beamSpotY", &beamSpotY, &b_beamSpotY);
   fChain->SetBranchAddress("beamSpotZ", &beamSpotZ, &b_beamSpotZ);
   fChain->SetBranchAddress("pvX", &pvX, &b_pvX);
   fChain->SetBranchAddress("pvY", &pvY, &b_pvY);
   fChain->SetBranchAddress("pvZ", &pvZ, &b_pvZ);
   fChain->SetBranchAddress("pvDX", &pvDX, &b_pvDX);
   fChain->SetBranchAddress("pvDY", &pvDY, &b_pvDY);
   fChain->SetBranchAddress("pvDZ", &pvDZ, &b_pvDZ);
   fChain->SetBranchAddress("pvRho", &pvRho, &b_pvRho);
   fChain->SetBranchAddress("pvIsFake", &pvIsFake, &b_pvIsFake);
   fChain->SetBranchAddress("pvNdof", &pvNdof, &b_pvNdof);
   fChain->SetBranchAddress("pvChi2", &pvChi2, &b_pvChi2);
   fChain->SetBranchAddress("mhtPt", &mhtPt, &b_mhtPt);
   fChain->SetBranchAddress("mhtPy", &mhtPy, &b_mhtPy);
   fChain->SetBranchAddress("mhtPx", &mhtPx, &b_mhtPx);
   fChain->SetBranchAddress("mhtPhi", &mhtPhi, &b_mhtPhi);
   fChain->SetBranchAddress("mhtSumEt", &mhtSumEt, &b_mhtSumEt);
   fChain->SetBranchAddress("mhtSignif", &mhtSignif, &b_mhtSignif);
   fChain->SetBranchAddress("nTriggerBits", &nTriggerBits, &b_nTriggerBits);
   fChain->SetBranchAddress("TriggerBits", &TriggerBits, &b_TriggerBits);
   fChain->SetBranchAddress("PileUpWeightRunA", &PileUpWeightRunA, &b_pileUpWeight);
   fChain->SetBranchAddress("PileUpWeightRunB", &PileUpWeightRunB, &b_pileUpWeight);
   fChain->SetBranchAddress("PileUpWeightRunC", &PileUpWeightRunC, &b_pileUpWeight);
   //Data trigger branches
   fChain->SetBranchAddress("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2", &HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2, &b_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2);
   fChain->SetBranchAddress("HLT_IsoMu20_v2", &HLT_IsoMu20_v2, &b_HLT_IsoMu20_v2);
   fChain->SetBranchAddress("HLT_IsoMu20_eta2p1_v2", &HLT_IsoMu20_eta2p1_v2, &b_HLT_IsoMu20_eta2p1_v2);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2", &HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2, &b_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v2", &HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v2, &b_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v2);
   fChain->SetBranchAddress("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3", &HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3, &b_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3);
   fChain->SetBranchAddress("HLT_IsoMu18_v1", &HLT_IsoMu18_v1, &b_HLT_IsoMu18_v1);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3", &HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3, &b_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3", &HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3, &b_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3);
   //MC trigger branches
   fChain->SetBranchAddress("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1", &HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1, &b_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1);
   fChain->SetBranchAddress("HLT_IsoMu20_v1", &HLT_IsoMu20_v1, &b_HLT_IsoMu20_v1);
   fChain->SetBranchAddress("HLT_IsoMu20_eta2p1_v1", &HLT_IsoMu20_eta2p1_v1, &b_HLT_IsoMu20_eta2p1_v1);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1", &HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1, &b_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1", &HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1, &b_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1);
   if (isMC) {
     fChain->SetBranchAddress("nGenPar", &nGenPar, &b_nGenPar);
     fChain->SetBranchAddress("genParEta", genParEta, &b_genParEta);
     fChain->SetBranchAddress("genParPhi", genParPhi, &b_genParPhi);
     fChain->SetBranchAddress("genParE", genParE, &b_genParE);
     fChain->SetBranchAddress("genParPt", genParPt, &b_genParPt);
     fChain->SetBranchAddress("genParId", genParId, &b_genParId);
     fChain->SetBranchAddress("genParCharge", genParCharge, &b_genParCharge);
   }
   fChain->SetBranchAddress("eventRun", &eventRun, &b_eventRun);
   fChain->SetBranchAddress("eventNum", &eventNum, &b_eventNum);
   fChain->SetBranchAddress("eventLumiblock", &eventLumiblock, &b_eventLumiblock);
   fChain->SetBranchAddress("numVert", &numVert, &b_numVert);
   Notify();
}

Bool_t AnalysisEvent::Notify()
{
  //  std::cout << "Does the notify." << std::endl;
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AnalysisEvent::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AnalysisEvent::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef AnalysisEvent_cxx

