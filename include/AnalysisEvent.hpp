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

#include <array>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class AnalysisEvent {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   static constexpr size_t NELECTRONSMAX{20};
   Int_t numElePF2PAT;
   std::array<Float_t, NELECTRONSMAX> elePF2PATE;
   std::array<Float_t, NELECTRONSMAX> elePF2PATET;
   std::array<Float_t, NELECTRONSMAX> elePF2PATPX;
   std::array<Float_t, NELECTRONSMAX> elePF2PATPY;
   std::array<Float_t, NELECTRONSMAX> elePF2PATPZ;
   std::array<Float_t, NELECTRONSMAX> elePF2PATPhi;
   std::array<Float_t, NELECTRONSMAX> elePF2PATTheta;
   std::array<Float_t, NELECTRONSMAX> elePF2PATEta;
   std::array<Float_t, NELECTRONSMAX> elePF2PATPT;
   std::array<Int_t,   NELECTRONSMAX> elePF2PATCharge;
   std::array<Float_t, NELECTRONSMAX> elePF2PATMVA;
   std::array<Int_t,   NELECTRONSMAX> elePF2PATCutIdVeto;
   std::array<Int_t,   NELECTRONSMAX> elePF2PATCutIdLoose;
   std::array<Int_t,   NELECTRONSMAX> elePF2PATCutIdMedium;
   std::array<Int_t,   NELECTRONSMAX> elePF2PATCutIdTight;
   std::array<Float_t, NELECTRONSMAX> elePF2PATImpactTransDist;
   std::array<Float_t, NELECTRONSMAX> elePF2PATImpactTransError;
   std::array<Float_t, NELECTRONSMAX> elePF2PATImpactTransSignificance;
   std::array<Float_t, NELECTRONSMAX> elePF2PATImpact3DDist;
   std::array<Float_t, NELECTRONSMAX> elePF2PATImpact3DError;
   std::array<Float_t, NELECTRONSMAX> elePF2PATImpact3DSignificance;
   std::array<Float_t, NELECTRONSMAX> elePF2PATChargedHadronIso;
   std::array<Float_t, NELECTRONSMAX> elePF2PATNeutralHadronIso;
   std::array<Float_t, NELECTRONSMAX> elePF2PATPhotonIso;
   std::array<Float_t, NELECTRONSMAX> elePF2PATTrackPt;
   std::array<Float_t, NELECTRONSMAX> elePF2PATTrackPhi;
   std::array<Float_t, NELECTRONSMAX> elePF2PATTrackEta;
   std::array<Float_t, NELECTRONSMAX> elePF2PATTrackChi2;
   std::array<Float_t, NELECTRONSMAX> elePF2PATTrackNDOF;
   std::array<Float_t, NELECTRONSMAX> elePF2PATTrackD0;
   std::array<Float_t, NELECTRONSMAX> elePF2PATTrackDBD0;
   std::array<Float_t, NELECTRONSMAX> elePF2PATD0PV;
   std::array<Float_t, NELECTRONSMAX> elePF2PATDZPV;
   std::array<Float_t, NELECTRONSMAX> elePF2PATBeamSpotCorrectedTrackD0;
   std::array<Float_t, NELECTRONSMAX> elePF2PATTrackDz;
   std::array<Float_t, NELECTRONSMAX> elePF2PATVtxZ;
   std::array<Int_t,   NELECTRONSMAX> elePF2PATIsGsf;
   std::array<Float_t, NELECTRONSMAX> elePF2PATGsfPx;
   std::array<Float_t, NELECTRONSMAX> elePF2PATGsfPy;
   std::array<Float_t, NELECTRONSMAX> elePF2PATGsfPz;
   std::array<Float_t, NELECTRONSMAX> elePF2PATGsfE;
   std::array<Float_t, NELECTRONSMAX> elePF2PATEcalEnergy;
   std::array<Float_t, NELECTRONSMAX> elePF2PATSCEta;
   std::array<Float_t, NELECTRONSMAX> elePF2PATSCE;
   std::array<Float_t, NELECTRONSMAX> elePF2PATSCPhi;
   std::array<Float_t, NELECTRONSMAX> elePF2PATSCEoverP;
   std::array<Float_t, NELECTRONSMAX> elePF2PATSCSigmaEtaEta;
   std::array<Float_t, NELECTRONSMAX> elePF2PATSCSigmaIEtaIEta;
   std::array<Float_t, NELECTRONSMAX> elePF2PATSCSigmaIEtaIEta5x5;
   std::array<Float_t, NELECTRONSMAX> elePF2PATSCE1x5;
   std::array<Float_t, NELECTRONSMAX> elePF2PATSCE5x5;
   std::array<Float_t, NELECTRONSMAX> elePF2PATSCE2x5max;
   std::array<Float_t, NELECTRONSMAX> elePF2PATTrackIso04;
   std::array<Float_t, NELECTRONSMAX> elePF2PATEcalIso04;
   std::array<Float_t, NELECTRONSMAX> elePF2PATHcalIso04;
   std::array<Float_t, NELECTRONSMAX> elePF2PATTrackIso03;
   std::array<Float_t, NELECTRONSMAX> elePF2PATEcalIso03;
   std::array<Float_t, NELECTRONSMAX> elePF2PATHcalIso03;
   std::array<Float_t, NELECTRONSMAX> elePF2PATdr04EcalRecHitSumEt;
   std::array<Float_t, NELECTRONSMAX> elePF2PATdr03EcalRecHitSumEt;
   std::array<Float_t, NELECTRONSMAX> elePF2PATEcalIsoDeposit;
   std::array<Float_t, NELECTRONSMAX> elePF2PATHcalIsoDeposit;
   std::array<Float_t, NELECTRONSMAX> elePF2PATComRelIso;
   std::array<Float_t, NELECTRONSMAX> elePF2PATComRelIsodBeta;
   std::array<Float_t, NELECTRONSMAX> elePF2PATComRelIsoRho;
   std::array<Float_t, NELECTRONSMAX> elePF2PATChHadIso;
   std::array<Float_t, NELECTRONSMAX> elePF2PATNtHadIso;
   std::array<Float_t, NELECTRONSMAX> elePF2PATGammaIso;
   std::array<Float_t, NELECTRONSMAX> elePF2PATRhoIso;
   std::array<Float_t, NELECTRONSMAX> elePF2PATAEff03;
   std::array<Int_t,   NELECTRONSMAX> elePF2PATMissingInnerLayers;
   std::array<Float_t, NELECTRONSMAX> elePF2PATHoverE;
   std::array<Float_t, NELECTRONSMAX> elePF2PATDeltaPhiSC;
   std::array<Float_t, NELECTRONSMAX> elePF2PATDeltaEtaSC;
   std::array<Float_t, NELECTRONSMAX> elePF2PATDeltaEtaSeedSC;
   std::array<Int_t,   NELECTRONSMAX> elePF2PATIsBarrel;
   std::array<Int_t,   NELECTRONSMAX> elePF2PATPhotonConversionTag;
   std::array<Float_t, NELECTRONSMAX> elePF2PATPhotonConversionDist;
   std::array<Float_t, NELECTRONSMAX> elePF2PATPhotonConversionDcot;
   std::array<Int_t,   NELECTRONSMAX> elePF2PATPhotonConversionVeto;
   std::array<Int_t,   NELECTRONSMAX> elePF2PATPhotonConversionTagCustom;
   std::array<Float_t, NELECTRONSMAX> elePF2PATPhotonConversionDistCustom;
   std::array<Float_t, NELECTRONSMAX> elePF2PATPhotonConversionDcotCustom;
   std::array<Float_t, NELECTRONSMAX> elePF2PATTriggerMatch;
   std::array<Float_t, NELECTRONSMAX> elePF2PATJetOverlap;
   std::array<Float_t, NELECTRONSMAX> genElePF2PATPT;
   std::array<Float_t, NELECTRONSMAX> genElePF2PATET;
   std::array<Float_t, NELECTRONSMAX> genElePF2PATPX;
   std::array<Float_t, NELECTRONSMAX> genElePF2PATPY;
   std::array<Float_t, NELECTRONSMAX> genElePF2PATPZ;
   std::array<Float_t, NELECTRONSMAX> genElePF2PATPhi;
   std::array<Float_t, NELECTRONSMAX> genElePF2PATTheta;
   std::array<Float_t, NELECTRONSMAX> genElePF2PATEta;
   std::array<Int_t,   NELECTRONSMAX> genElePF2PATCharge;
   std::array<Int_t,   NELECTRONSMAX> genElePF2PATPdgId;
   std::array<Int_t,   NELECTRONSMAX> genElePF2PATMotherId;
   std::array<Int_t,   NELECTRONSMAX> genElePF2PATPromptDecayed;
   std::array<Int_t,   NELECTRONSMAX> genElePF2PATPromptFinalState;
   std::array<Int_t,   NELECTRONSMAX> genElePF2PATHardProcess;
   static constexpr size_t NMUONSMAX{20};
   Int_t numMuonPF2PAT;
   std::array<Float_t, NMUONSMAX> muonPF2PATE;
   std::array<Float_t, NMUONSMAX> muonPF2PATET;
   std::array<Float_t, NMUONSMAX> muonPF2PATPt;
   std::array<Float_t, NMUONSMAX> muonPF2PATPX;
   std::array<Float_t, NMUONSMAX> muonPF2PATPY;
   std::array<Float_t, NMUONSMAX> muonPF2PATPZ;
   std::array<Float_t, NMUONSMAX> muonPF2PATPhi;
   std::array<Float_t, NMUONSMAX> muonPF2PATTheta;
   std::array<Float_t, NMUONSMAX> muonPF2PATEta;
   std::array<Int_t,   NMUONSMAX> muonPF2PATCharge;
   std::array<Float_t, NMUONSMAX> muonPF2PATGlobalID;
   std::array<Float_t, NMUONSMAX> muonPF2PATTrackID;
   std::array<Float_t, NMUONSMAX> muonPF2PATChi2;
   std::array<Float_t, NMUONSMAX> muonPF2PATD0;
   std::array<Float_t, NMUONSMAX> muonPF2PATTrackDBD0;
   std::array<Float_t, NMUONSMAX> muonPF2PATDBInnerTrackD0;
   std::array<Float_t, NMUONSMAX> muonPF2PATBeamSpotCorrectedD0;
   std::array<Int_t,   NMUONSMAX> muonPF2PATTrackNHits;
   std::array<Int_t,   NMUONSMAX> muonPF2PATMuonNHits;
   std::array<Float_t, NMUONSMAX> muonPF2PATNDOF;
   std::array<Float_t, NMUONSMAX> muonPF2PATVertX;
   std::array<Float_t, NMUONSMAX> muonPF2PATVertY;
   std::array<Float_t, NMUONSMAX> muonPF2PATVertZ;
   std::array<Float_t, NMUONSMAX> muonPF2PATChargedHadronIso;
   std::array<Float_t, NMUONSMAX> muonPF2PATNeutralHadronIso;
   std::array<Float_t, NMUONSMAX> muonPF2PATPhotonIso;
   std::array<Float_t, NMUONSMAX> muonPF2PATTrackIso;
   std::array<Float_t, NMUONSMAX> muonPF2PATEcalIso;
   std::array<Float_t, NMUONSMAX> muonPF2PATHcalIso;
   std::array<Float_t, NMUONSMAX> muonPF2PATComRelIso;
   std::array<Float_t, NMUONSMAX> muonPF2PATComRelIsodBeta;
   std::array<Int_t,   NMUONSMAX> muonPF2PATIsPFMuon;
   std::array<Int_t,   NMUONSMAX> muonPF2PATNChambers;
   std::array<Int_t,   NMUONSMAX> muonPF2PATNMatches;
   std::array<Int_t,   NMUONSMAX> muonPF2PATTkLysWithMeasurements;
   std::array<Int_t,   NMUONSMAX> muonPF2PATVldPixHits;
   std::array<Int_t,   NMUONSMAX> muonPF2PATMatchedStations;
   std::array<Float_t, NMUONSMAX> muonPF2PATGlbTkNormChi2;
   std::array<Float_t, NMUONSMAX> muonPF2PATValidFraction;
   std::array<Float_t, NMUONSMAX> muonPF2PATChi2LocalPosition;
   std::array<Float_t, NMUONSMAX> muonPF2PATTrkKick;
   std::array<Float_t, NMUONSMAX> muonPF2PATSegmentCompatibility;
   std::array<Float_t, NMUONSMAX> muonPF2PATDBPV;
   std::array<Float_t, NMUONSMAX> muonPF2PATDZPV;
   std::array<Float_t, NMUONSMAX> genMuonPF2PATPT;
   std::array<Float_t, NMUONSMAX> genMuonPF2PATET;
   std::array<Float_t, NMUONSMAX> genMuonPF2PATPX;
   std::array<Float_t, NMUONSMAX> genMuonPF2PATPY;
   std::array<Float_t, NMUONSMAX> genMuonPF2PATPZ;
   std::array<Float_t, NMUONSMAX> genMuonPF2PATPhi;
   std::array<Float_t, NMUONSMAX> genMuonPF2PATTheta;
   std::array<Float_t, NMUONSMAX> genMuonPF2PATEta;
   std::array<Int_t,   NMUONSMAX> genMuonPF2PATCharge;
   std::array<Int_t,   NMUONSMAX> genMuonPF2PATPdgId;
   std::array<Int_t,   NMUONSMAX> genMuonPF2PATMotherId;
   std::array<Int_t,   NMUONSMAX> genMuonPF2PATPromptDecayed;
   std::array<Int_t,   NMUONSMAX> genMuonPF2PATPromptFinalState;
   std::array<Int_t,   NMUONSMAX> genMuonPF2PATHardProcess;
   static constexpr size_t NJETSMAX{40};
   Int_t                  numJetPF2PAT;
   std::array<Double_t, NJETSMAX> jetPF2PATE;
   std::array<Double_t, NJETSMAX> jetPF2PATEt;
   std::array<Double_t, NJETSMAX> jetPF2PATPt;
   std::array<Double_t, NJETSMAX> jetPF2PATPtRaw;
   std::array<Double_t, NJETSMAX> jetPF2PATUnCorEt;
   std::array<Double_t, NJETSMAX> jetPF2PATUnCorPt;
   std::array<Double_t, NJETSMAX> jetPF2PATEta;
   std::array<Double_t, NJETSMAX> jetPF2PATTheta;
   std::array<Double_t, NJETSMAX> jetPF2PATPhi;
   std::array<Double_t, NJETSMAX> jetPF2PATPx;
   std::array<Double_t, NJETSMAX> jetPF2PATPy;
   std::array<Double_t, NJETSMAX> jetPF2PATPz;
   std::array<Double_t, NJETSMAX> jetPF2PATdRClosestLepton;
   std::array<Int_t,    NJETSMAX> jetPF2PATNtracksInJet;
   std::array<Float_t,  NJETSMAX> jetPF2PATJetCharge;
   std::array<Float_t,  NJETSMAX> jetPF2PATfHPD;
   std::array<Float_t,  NJETSMAX> jetPF2PATBtagSoftMuonPtRel;
   std::array<Float_t,  NJETSMAX> jetPF2PATBtagSoftMuonQuality;
   std::array<Float_t,  NJETSMAX> jetPF2PATCorrFactor;
   std::array<Float_t,  NJETSMAX> jetPF2PATCorrResidual;
   std::array<Float_t,  NJETSMAX> jetPF2PATL2L3ResErr;
   std::array<Float_t,  NJETSMAX> jetPF2PATCorrErrLow;
   std::array<Float_t,  NJETSMAX> jetPF2PATCorrErrHi;
   std::array<Float_t,  NJETSMAX> jetPF2PATN90Hits;
   std::array<Float_t,  NJETSMAX> jetPF2PATTriggered;
   std::array<Float_t,  NJETSMAX> jetPF2PATSVX;
   std::array<Float_t,  NJETSMAX> jetPF2PATSVY;
   std::array<Float_t,  NJETSMAX> jetPF2PATSVZ;
   std::array<Float_t,  NJETSMAX> jetPF2PATSVDX;
   std::array<Float_t,  NJETSMAX> jetPF2PATSVDY;
   std::array<Float_t,  NJETSMAX> jetPF2PATSVDZ;
   std::array<Float_t,  NJETSMAX> jetPF2PATBDiscriminator;
//   std::array<Float_t,  NJETSMAX> jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags;
//   std::array<Float_t,  NJETSMAX> jetPF2PATpfCombinedMVAV2BJetTags;
//   std::array<Float_t,  NJETSMAX> jetPF2PATpfCombinedCvsLJetTags;
//   std::array<Float_t,  NJETSMAX> jetPF2PATpfCombinedCvsBJetTags;
   std::array<Int_t,    NJETSMAX> jetPF2PATNConstituents;
   std::array<Int_t,    NJETSMAX> jetPF2PATPID;
   std::array<Float_t,  NJETSMAX> jetPF2PATClosestBPartonDeltaR;
   std::array<Float_t,  NJETSMAX> jetPF2PATClosestCPartonDeltaR;
   std::array<Float_t,  NJETSMAX> genJetPF2PATET;
   std::array<Float_t,  NJETSMAX> genJetPF2PATPT;
   std::array<Float_t,  NJETSMAX> genJetPF2PATPX;
   std::array<Float_t,  NJETSMAX> genJetPF2PATPY;
   std::array<Float_t,  NJETSMAX> genJetPF2PATPZ;
   std::array<Float_t,  NJETSMAX> genJetPF2PATPhi;
   std::array<Float_t,  NJETSMAX> genJetPF2PATTheta;
   std::array<Float_t,  NJETSMAX> genJetPF2PATEta;
   std::array<Int_t,    NJETSMAX> genJetPF2PATPID;
   std::array<Float_t,  NJETSMAX> jetPF2PATMuEnergy;
   std::array<Float_t,  NJETSMAX> jetPF2PATMuEnergyFraction;
   std::array<Float_t,  NJETSMAX> jetPF2PATNeutralHadEnergy;
   std::array<Float_t,  NJETSMAX> jetPF2PATNeutralEmEnergy;
   std::array<Float_t,  NJETSMAX> jetPF2PATChargedHadronEnergyFraction;
   std::array<Float_t,  NJETSMAX> jetPF2PATNeutralHadronEnergyFraction;
   std::array<Float_t,  NJETSMAX> jetPF2PATChargedEmEnergyFraction;
   std::array<Float_t,  NJETSMAX> jetPF2PATNeutralEmEnergyFraction;
   std::array<Float_t,  NJETSMAX> jetPF2PATChargedHadronEnergyFractionCorr;
   std::array<Float_t,  NJETSMAX> jetPF2PATNeutralHadronEnergyFractionCorr;
   std::array<Float_t,  NJETSMAX> jetPF2PATChargedEmEnergyFractionCorr;
   std::array<Float_t,  NJETSMAX> jetPF2PATNeutralEmEnergyFractionCorr;
   std::array<Int_t,    NJETSMAX> jetPF2PATNeutralMultiplicity;
   std::array<Int_t,    NJETSMAX> jetPF2PATChargedMultiplicity;
   Double_t metPF2PATE;
   Double_t metPF2PATEt;
   Double_t metPF2PATEtRaw;
   Double_t metPF2PATPhi;
   Double_t metPF2PATPt;
   Double_t metPF2PATPx;
   Double_t metPF2PATPy;
   Double_t metPF2PATPz;
   Float_t  metPF2PATScalarEt;
   Float_t  metPF2PATEtUncorrected;
   Float_t  metPF2PATPhiUncorrected;
   Float_t  metPF2PATUnclusteredEnUp;
   Float_t  metPF2PATUnclusteredEnDown;
   Float_t  genMetPF2PATEt;
   Float_t  genMetPF2PATPhi;
   Float_t  genMetPF2PATPt;
   Float_t  genMetPF2PATPx;
   Float_t  genMetPF2PATPy;
   Float_t  genMetPF2PATPz;
   static constexpr size_t NTAUSMAX{1};
   Int_t numTauPF2PAT;
   std::array<Float_t, NTAUSMAX> tauPF2PATE;
   std::array<Float_t, NTAUSMAX> tauPF2PATPt;
   std::array<Float_t, NTAUSMAX> tauPF2PATPhi;
   std::array<Float_t, NTAUSMAX> tauPF2PATEta;
   static constexpr size_t NTRACKSMAX{1000};
   Int_t numGeneralTracks;
   std::array<Float_t, NTRACKSMAX> generalTracksPt;
   std::array<Float_t, NTRACKSMAX> generalTracksEta;
   std::array<Float_t, NTRACKSMAX> generalTracksTheta;
   std::array<Float_t, NTRACKSMAX> generalTracksBeamSpotCorrectedD0;
   std::array<Float_t, NTRACKSMAX> generalTracksPhi;
   std::array<Int_t,   NTRACKSMAX> generalTracksCharge;
   Int_t    isElePlusJets;
   Float_t  genPDFScale;
   Float_t  genPDFx1;
   Float_t  genPDFx2;
   Int_t    genPDFf1;
   Int_t    genPDFf2;
   Double_t  topPtReweight;
   Int_t    processId;
   Float_t  processPtHat;
   Double_t processMCWeight;
   Float_t  beamSpotX;
   Float_t  beamSpotY;
   Float_t  beamSpotZ;
   Float_t  pvX;
   Float_t  pvY;
   Float_t  pvZ;
   Float_t  pvDX;
   Float_t  pvDY;
   Float_t  pvDZ;
   Float_t  pvRho;
   Int_t    pvIsFake;
   Float_t  pvNdof;
   Float_t  pvChi2;
   Float_t  mhtPt;
   Float_t  mhtPy;
   Float_t  mhtPx;
   Float_t  mhtPhi;
   Float_t  mhtSumEt;
   Float_t  mhtSignif;
   static constexpr size_t NTRIGGERBITSMAX{1};
   Int_t nTriggerBits;
   std::array<Int_t, NTRIGGERBITSMAX> TriggerBits;
   Double_t weight_muF0p5;
   Double_t weight_muF2;
   Double_t weight_muR0p5;
   Double_t weight_muR2;
   Double_t weight_muF0p5muR0p5;
   Double_t weight_muF2muR2;
   Double_t origWeightForNorm;
   Double_t weight_pdfMax;
   Double_t weight_pdfMin;
   Double_t weight_alphaMax;
   Double_t weight_alphaMin;
   // Int_t    numVert;

   //2015 Data Triggers
   Int_t           HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2;
   Int_t           HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2;
   Int_t           HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v2;
   Int_t           HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3;
   Int_t           HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3;
   Int_t           HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3;

   //2015 MC Triggers
   Int_t           HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1;
   Int_t	   HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1;
   Int_t	   HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1;
   Int_t           HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1;
   Int_t           HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1;

   //2015 MET Triggers
   Int_t           HLT_PFMET170_JetIdCleaned_v2;
   Int_t           HLT_PFHT350_PFMET100_v1;

   //2015 MET Filters
   Int_t	   Flag_CSCTightHalo2015Filter;

   //2016 Triggers
   Int_t           HLT_Ele25_eta2p1_WPTight_Gsf_v1;
   Int_t           HLT_Ele25_eta2p1_WPTight_Gsf_v2;
   Int_t           HLT_Ele25_eta2p1_WPTight_Gsf_v3;
   Int_t           HLT_Ele25_eta2p1_WPTight_Gsf_v4;
   Int_t           HLT_Ele25_eta2p1_WPTight_Gsf_v5;
   Int_t           HLT_Ele25_eta2p1_WPTight_Gsf_v6;
   Int_t           HLT_Ele25_eta2p1_WPTight_Gsf_v7;
   Int_t           HLT_Ele27_WPTight_Gsf_v1;
   Int_t           HLT_Ele27_WPTight_Gsf_v2;
   Int_t           HLT_Ele27_WPTight_Gsf_v3;
   Int_t           HLT_Ele27_WPTight_Gsf_v4;
   Int_t           HLT_Ele27_WPTight_Gsf_v5;
   Int_t           HLT_Ele27_WPTight_Gsf_v6;
   Int_t           HLT_Ele27_WPTight_Gsf_v7;
   Int_t	   HLT_Ele32_eta2p1_WPTight_Gsf_v2;
   Int_t	   HLT_Ele32_eta2p1_WPTight_Gsf_v3;
   Int_t	   HLT_Ele32_eta2p1_WPTight_Gsf_v4;
   Int_t	   HLT_Ele32_eta2p1_WPTight_Gsf_v5;
   Int_t	   HLT_Ele32_eta2p1_WPTight_Gsf_v6;
   Int_t	   HLT_Ele32_eta2p1_WPTight_Gsf_v7;
   Int_t	   HLT_Ele32_eta2p1_WPTight_Gsf_v8;
   Int_t           HLT_IsoMu24_v1;
   Int_t           HLT_IsoMu24_v2;
   Int_t           HLT_IsoMu24_v3;
   Int_t           HLT_IsoMu24_v4;
   Int_t           HLT_IsoTkMu24_v1;
   Int_t           HLT_IsoTkMu24_v2;
   Int_t           HLT_IsoTkMu24_v3;
   Int_t           HLT_IsoTkMu24_v4;
   Int_t	   HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3;
   Int_t	   HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4;
   Int_t	   HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5;
   Int_t	   HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v6;
   Int_t	   HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v7;
   Int_t	   HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v8;
   Int_t	   HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v9;


   Int_t           HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v2;
   Int_t           HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v3;
   Int_t           HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v4;
   Int_t           HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v6;
   Int_t           HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v2;
   Int_t           HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v3;
   Int_t           HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v5;

//   Int_t	   HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2; // delcared elsewhere
   Int_t	   HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v3;
   Int_t	   HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4;
   Int_t	   HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7;
//   Int_t	   HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2; // delcared elsewhere
   Int_t	   HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3;
   Int_t	   HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6;

   Int_t	   HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3;
   Int_t	   HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v4;
   Int_t	   HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v5;
   Int_t	   HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v6;
   Int_t	   HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v7;
   Int_t	   HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v8;
   Int_t	   HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v9;
   Int_t	   HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3;
   Int_t	   HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v4;
   Int_t	   HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v5;
   Int_t	   HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v6;
   Int_t	   HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v7;
   Int_t	   HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v8;
   Int_t	   HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v9;

   Int_t	   HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v1;
   Int_t	   HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3;

   Int_t	   HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1;
   Int_t	   HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2;
   Int_t	   HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3;
   Int_t	   HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4;
   Int_t	   HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v1;
   Int_t	   HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v2;
   Int_t	   HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v3;
   Int_t	   HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v4;

   //2016 MET Triggers
   Int_t	   HLT_MET200_v1;
   Int_t	   HLT_MET200_v2;
   Int_t	   HLT_MET200_v3;
   Int_t	   HLT_MET200_v4;
   Int_t	   HLT_MET200_v5;
   Int_t	   HLT_MET250_v2;
   Int_t	   HLT_MET250_v3;
   Int_t	   HLT_MET250_v4;
   Int_t	   HLT_MET250_v5;
   Int_t	   HLT_PFMET120_PFMHT120_IDTight_v3;
   Int_t	   HLT_PFMET120_PFMHT120_IDTight_v4;
   Int_t	   HLT_PFMET120_PFMHT120_IDTight_v5;
   Int_t	   HLT_PFMET120_PFMHT120_IDTight_v6;
   Int_t	   HLT_PFMET120_PFMHT120_IDTight_v7;
   Int_t	   HLT_PFMET120_PFMHT120_IDTight_v8;
   Int_t	   HLT_PFMET170_HBHECleaned_v3;
   Int_t	   HLT_PFMET170_HBHECleaned_v4;
   Int_t	   HLT_PFMET170_HBHECleaned_v5;
   Int_t	   HLT_PFMET170_HBHECleaned_v6;
   Int_t	   HLT_PFMET170_HBHECleaned_v7;
   Int_t	   HLT_PFMET170_HBHECleaned_v8;
   Int_t	   HLT_PFMET170_HBHECleaned_v9;
   Int_t	   HLT_PFHT800_v3;
   Int_t	   HLT_PFHT800_v4;
   Int_t	   HLT_PFHT800_v5;
   Int_t	   HLT_PFHT900_v4;
   Int_t	   HLT_PFHT900_v5;
   Int_t	   HLT_PFHT900_v6;
   Int_t	   HLT_PFHT750_4JetPt50_v4;
   Int_t	   HLT_PFHT750_4JetPt50_v5;
   Int_t	   HLT_PFHT750_4JetPt50_v6;
   Int_t	   HLT_PFHT750_4JetPt70_v1;
   Int_t	   HLT_PFHT750_4JetPt70_v2;
   Int_t	   HLT_PFHT750_4JetPt80_v2;
   Int_t	   HLT_PFHT300_PFMET100_v1;
   Int_t	   HLT_PFHT300_PFMET100_v2;
   Int_t	   HLT_PFHT300_PFMET100_v3;
   Int_t	   HLT_PFHT300_PFMET100_v4;
   Int_t	   HLT_PFHT300_PFMET110_v4;
   Int_t	   HLT_PFHT300_PFMET110_v5;
   Int_t	   HLT_PFHT300_PFMET110_v6;

   //2016 MET Filters
   Int_t	   Flag_globalTightHalo2016Filter;
//   Int_t           Flag_BadChargedCandidateFilter;
   Int_t           Flag_chargedHadronTrackResolutionFilter;
//   Int_t	   Flag_BadPFMuonFilter;
   Int_t           Flag_muonBadTrackFilter;
   Int_t	   Flag_ecalLaserCorrFilter;

   //Feb2017 rereco Filters
   Int_t	   Flag_badMuons;
   Int_t	   Flag_duplicateMuons;
   Int_t	   Flag_noBadMuons;

   //2015 and 2016 Triggers
   Int_t           HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2;
   Int_t           HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2;

   //2015 and 2016 MET Triggers
   Int_t           HLT_PFMET120_PFMHT120_IDTight_v2;
   Int_t	   HLT_PFMET170_HBHECleaned_v2;
   Int_t	   HLT_PFHT800_v2;
   Int_t           HLT_MET250_v1;
   Int_t           HLT_PFHT750_4JetPt50_v3;

   //2015 and 2016 MET Filters
   Int_t	   Flag_HBHENoiseFilter;
   Int_t	   Flag_HBHENoiseIsoFilter;
   Int_t	   Flag_EcalDeadCellTriggerPrimitiveFilter;
   Int_t	   Flag_goodVertices;
   Int_t	   Flag_eeBadScFilter;

   //Gen info
   static const size_t NGENPARMAX{50};
   Int_t nGenPar;
   std::array<Float_t, NGENPARMAX> genParEta;
   std::array<Float_t, NGENPARMAX> genParPhi;
   std::array<Float_t, NGENPARMAX> genParE;
   std::array<Float_t, NGENPARMAX> genParPt;
   std::array<Int_t,   NGENPARMAX> genParId;
   std::array<Int_t,   NGENPARMAX> genParMotherId;
   std::array<Int_t,   NGENPARMAX> genParCharge;
   Int_t   eventRun;
   Int_t   eventNum;
   Float_t eventLumiblock;

   // List of branches
   TBranch        *b_numElePF2PAT;   //!
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
   TBranch        *b_elePF2PATCutIdVeto;   //!
   TBranch        *b_elePF2PATCutIdLoose;   //!
   TBranch        *b_elePF2PATCutIdMedium;   //!
   TBranch        *b_elePF2PATCutIdTight;   //!
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
   TBranch        *b_elePF2PATEcalEnergy;   //!
   TBranch        *b_elePF2PATSCEta;   //!
   TBranch        *b_elePF2PATSCE;   //!
   TBranch        *b_elePF2PATSCPhi;   //!
   TBranch        *b_elePF2PATSCEoverP;   //!
   TBranch        *b_elePF2PATSCSigmaEtaEta;   //!
   TBranch        *b_elePF2PATSCSigmaIEtaIEta;   //!
   TBranch        *b_elePF2PATSCSigmaIEtaIEta5x5;   //!
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
   TBranch        *b_elePF2PATDeltaEtaSeedSC;   //!
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
   TBranch        *b_genElePF2PATPT;   //!
   TBranch        *b_genElePF2PATET;   //!
   TBranch        *b_genElePF2PATPX;   //!
   TBranch        *b_genElePF2PATPY;   //!
   TBranch        *b_genElePF2PATPZ;   //!
   TBranch        *b_genElePF2PATPhi;   //!
   TBranch        *b_genElePF2PATTheta;   //!
   TBranch        *b_genElePF2PATEta;   //!
   TBranch        *b_genElePF2PATCharge;   //!
   TBranch        *b_genElePF2PATPdgId;   //!
   TBranch        *b_genElePF2PATMotherId;   //!
   TBranch        *b_genElePF2PATPromptDecayed;   //!
   TBranch        *b_genElePF2PATPromptFinalState;   //!
   TBranch        *b_genElePF2PATHardProcess;   //!
   TBranch        *b_numMuonPF2PAT;   //!
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
   TBranch        *b_muonPF2PATValidFraction;   //!
   TBranch        *b_muonPF2PATChi2LocalPosition;   //!
   TBranch        *b_muonPF2PATTrkKick;   //!
   TBranch        *b_muonPF2PATSegmentCompatibility;   //!
   TBranch        *b_muonPF2PATDBPV;   //!
   TBranch        *b_muonPF2PATDZPV;   //!
   TBranch        *b_genMuonPF2PATPT;   //!
   TBranch        *b_genMuonPF2PATET;   //!
   TBranch        *b_genMuonPF2PATPX;   //!
   TBranch        *b_genMuonPF2PATPY;   //!
   TBranch        *b_genMuonPF2PATPZ;   //!
   TBranch        *b_genMuonPF2PATPhi;   //!
   TBranch        *b_genMuonPF2PATTheta;   //!
   TBranch        *b_genMuonPF2PATEta;   //!
   TBranch        *b_genMuonPF2PATCharge;   //!
   TBranch        *b_genMuonPF2PATPdgId;   //!
   TBranch        *b_genMuonPF2PATMotherId;   //!
   TBranch        *b_genMuonPF2PATPromptDecayed;   //!
   TBranch        *b_genMuonPF2PATPromptFinalState;   //!
   TBranch        *b_genMuonPF2PATHardProcess;   //!
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
   TBranch        *b_jetPF2PATSVX;   //!
   TBranch        *b_jetPF2PATSVY;   //!
   TBranch        *b_jetPF2PATSVZ;   //!
   TBranch        *b_jetPF2PATSVDX;   //!
   TBranch        *b_jetPF2PATSVDY;   //!
   TBranch        *b_jetPF2PATSVDZ;   //!
   TBranch        *b_jetPF2PATBDiscriminator;   //!
//   TBranch        *b_jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags;   //!
//   TBranch        *b_jetPF2PATpfCombinedMVAV2BJetTags;   //!
//   TBranch        *b_jetPF2PATpfCombinedCvsLJetTags;   //!
//   TBranch        *b_jetPF2PATpfCombinedCvsBJetTags;   //!
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
   TBranch        *b_metPF2PATE;   //!
   TBranch        *b_metPF2PATEt;   //!
   TBranch        *b_metPF2PATEtRaw;   //!
   TBranch        *b_metPF2PATPhi;   //!
   TBranch        *b_metPF2PATPt;   //!
   TBranch        *b_metPF2PATPx;   //!
   TBranch        *b_metPF2PATPy;   //!
   TBranch        *b_metPF2PATPz;   //!
   TBranch        *b_metPF2PATScalarEt;   //!
   TBranch        *b_metPF2PATEtUncorrected;   //!
   TBranch        *b_metPF2PATPhiUncorrected;   //!
   TBranch        *b_metPF2PATUnclusteredEnUp;   //!
   TBranch        *b_metPF2PATUnclusteredEnDown;   //!
   TBranch        *b_genMetPF2PATEt;   //!
   TBranch        *b_genMetPF2PATPhi;   //!
   TBranch        *b_genMetPF2PATPt;   //!
   TBranch        *b_genMetPF2PATPx;   //!
   TBranch        *b_genMetPF2PATPy;   //!
   TBranch        *b_genMetPF2PATPz;   //!
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
   TBranch        *b_topPtReweight;   //!
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
   TBranch	  *b_weight_muF0p5;   //!
   TBranch	  *b_weight_muF2;   //!
   TBranch	  *b_weight_muR0p5;   //!
   TBranch	  *b_weight_muR2;   //!
   TBranch	  *b_weight_muF0p5muR0p5;   //!
   TBranch	  *b_weight_muF2muR2;   //!
   TBranch	  *b_origWeightForNorm;   //!
   TBranch	  *b_weight_pdfMax;   //!
   TBranch	  *b_weight_pdfMin;   //!
   TBranch	  *b_weight_alphaMax;   //!
   TBranch	  *b_weight_alphaMin;   //!
   //   TBranch        *b_numVert;    //!

   //2015 Lepton and MET Triggers
   TBranch        *b_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2;
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2;
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v2;
   TBranch        *b_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3;
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3;
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3;
   TBranch        *b_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1;
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1;
   TBranch        *b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1;
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1;
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1;
   TBranch        *b_HLT_PFMET170_JetIdCleaned_v2;
   TBranch        *b_HLT_PFHT350_PFMET100_v1;

   //2015 MET Filters
   TBranch        *b_Flag_CSCTightHalo2015Filter;

   //2016 Lepton Triggers
   TBranch 	   *b_HLT_Ele25_eta2p1_WPTight_Gsf_v1;
   TBranch 	   *b_HLT_Ele25_eta2p1_WPTight_Gsf_v2;
   TBranch	   *b_HLT_Ele25_eta2p1_WPTight_Gsf_v3;
   TBranch	   *b_HLT_Ele25_eta2p1_WPTight_Gsf_v4;
   TBranch	   *b_HLT_Ele25_eta2p1_WPTight_Gsf_v5;
   TBranch	   *b_HLT_Ele25_eta2p1_WPTight_Gsf_v6;
   TBranch	   *b_HLT_Ele25_eta2p1_WPTight_Gsf_v7;
   TBranch         *b_HLT_Ele27_WPTight_Gsf_v1;
   TBranch         *b_HLT_Ele27_WPTight_Gsf_v2;
   TBranch         *b_HLT_Ele27_WPTight_Gsf_v3;
   TBranch         *b_HLT_Ele27_WPTight_Gsf_v4;
   TBranch         *b_HLT_Ele27_WPTight_Gsf_v5;
   TBranch         *b_HLT_Ele27_WPTight_Gsf_v6;
   TBranch         *b_HLT_Ele27_WPTight_Gsf_v7;
   TBranch         *b_HLT_Ele32_eta2p1_WPTight_Gsf_v2;
   TBranch         *b_HLT_Ele32_eta2p1_WPTight_Gsf_v3;
   TBranch         *b_HLT_Ele32_eta2p1_WPTight_Gsf_v4;
   TBranch         *b_HLT_Ele32_eta2p1_WPTight_Gsf_v5;
   TBranch         *b_HLT_Ele32_eta2p1_WPTight_Gsf_v6;
   TBranch         *b_HLT_Ele32_eta2p1_WPTight_Gsf_v7;
   TBranch         *b_HLT_Ele32_eta2p1_WPTight_Gsf_v8;
   TBranch         *b_HLT_IsoMu24_v1;
   TBranch         *b_HLT_IsoMu24_v2;
   TBranch         *b_HLT_IsoMu24_v3;
   TBranch         *b_HLT_IsoMu24_v4;
   TBranch         *b_HLT_IsoTkMu24_v1;
   TBranch         *b_HLT_IsoTkMu24_v2;
   TBranch         *b_HLT_IsoTkMu24_v3;
   TBranch         *b_HLT_IsoTkMu24_v4;
   TBranch	   *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3;
   TBranch	   *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4;
   TBranch	   *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5;
   TBranch	   *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v6;
   TBranch	   *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v7;
   TBranch	   *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v8;
   TBranch	   *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v9;

   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v2;
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v3;
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v4;
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v6;
   TBranch        *b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v2;
   TBranch        *b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v3;
   TBranch        *b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v5;
//   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2; // declared elsewhere
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v3;
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4;
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7;
//   TBranch        *b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2; // declared elsewhere
   TBranch        *b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3;
   TBranch        *b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6;

   TBranch	   *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3;
   TBranch	   *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v4;
   TBranch	   *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v5;
   TBranch	   *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v6;
   TBranch	   *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v7;
   TBranch	   *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v8;
   TBranch	   *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v9;
   TBranch	   *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3;
   TBranch	   *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v4;
   TBranch	   *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v5;
   TBranch	   *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v6;
   TBranch	   *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v7;
   TBranch	   *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v8;
   TBranch	   *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v9;

   TBranch	   *b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v1;
   TBranch	   *b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3;

   TBranch	   *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1;
   TBranch	   *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2;
   TBranch	   *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3;
   TBranch	   *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4;
   TBranch	   *b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v1;
   TBranch	   *b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v2;
   TBranch	   *b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v3;
   TBranch	   *b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v4;

   //2016 MET Triggers
   TBranch	   *b_HLT_MET200_v1;
   TBranch	   *b_HLT_MET200_v2;
   TBranch	   *b_HLT_MET200_v3;
   TBranch	   *b_HLT_MET200_v4;
   TBranch	   *b_HLT_MET200_v5;
   TBranch	   *b_HLT_MET250_v2;
   TBranch	   *b_HLT_MET250_v3;
   TBranch	   *b_HLT_MET250_v4;
   TBranch	   *b_HLT_MET250_v5;
   TBranch	   *b_HLT_PFMET120_PFMHT120_IDTight_v3;
   TBranch	   *b_HLT_PFMET120_PFMHT120_IDTight_v4;
   TBranch	   *b_HLT_PFMET120_PFMHT120_IDTight_v5;
   TBranch	   *b_HLT_PFMET120_PFMHT120_IDTight_v6;
   TBranch	   *b_HLT_PFMET120_PFMHT120_IDTight_v7;
   TBranch	   *b_HLT_PFMET120_PFMHT120_IDTight_v8;
   TBranch	   *b_HLT_PFMET170_HBHECleaned_v3;
   TBranch	   *b_HLT_PFMET170_HBHECleaned_v4;
   TBranch	   *b_HLT_PFMET170_HBHECleaned_v5;
   TBranch	   *b_HLT_PFMET170_HBHECleaned_v6;
   TBranch	   *b_HLT_PFMET170_HBHECleaned_v7;
   TBranch	   *b_HLT_PFMET170_HBHECleaned_v8;
   TBranch	   *b_HLT_PFMET170_HBHECleaned_v9;
   TBranch	   *b_HLT_PFHT800_v3;
   TBranch	   *b_HLT_PFHT800_v4;
   TBranch	   *b_HLT_PFHT800_v5;
   TBranch	   *b_HLT_PFHT900_v4;
   TBranch	   *b_HLT_PFHT900_v5;
   TBranch	   *b_HLT_PFHT900_v6;
   TBranch	   *b_HLT_PFHT750_4JetPt50_v4;
   TBranch	   *b_HLT_PFHT750_4JetPt50_v5;
   TBranch	   *b_HLT_PFHT750_4JetPt50_v6;
   TBranch	   *b_HLT_PFHT750_4JetPt70_v1;
   TBranch	   *b_HLT_PFHT750_4JetPt70_v2;
   TBranch	   *b_HLT_PFHT750_4JetPt80_v2;
   TBranch	   *b_HLT_PFHT300_PFMET100_v1;
   TBranch	   *b_HLT_PFHT300_PFMET100_v2;
   TBranch	   *b_HLT_PFHT300_PFMET100_v3;
   TBranch	   *b_HLT_PFHT300_PFMET100_v4;
   TBranch	   *b_HLT_PFHT300_PFMET110_v4;
   TBranch	   *b_HLT_PFHT300_PFMET110_v5;
   TBranch	   *b_HLT_PFHT300_PFMET110_v6;

   //2016 MET Filters
   TBranch        *b_Flag_globalTightHalo2016Filter;
//   TBranch	  *b_Flag_BadChargedCandidateFilter;
   TBranch        *b_Flag_chargedHadronTrackResolutionFilter;
//   TBranch	  *b_Flag_BadPFMuonFilter;
   TBranch        *b_Flag_muonBadTrackFilter;
   TBranch        *b_Flag_ecalLaserCorrFilter;

   //Feb2017 rereco Filters
   TBranch        *b_Flag_badMuons;
   TBranch        *b_Flag_duplicateMuons;
   TBranch        *b_Flag_noBadMuons;

   //2015 and 2016 MET Triggers
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2;
   TBranch        *b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2;

   //2015 and 2016 MET Filters
   TBranch        *b_HLT_PFMET120_PFMHT120_IDTight_v2;
   TBranch        *b_HLT_PFMET170_HBHECleaned_v2;
   TBranch        *b_HLT_PFHT800_v2;
   TBranch        *b_HLT_MET250_v1;
   TBranch        *b_HLT_PFHT750_4JetPt50_v3;

   //2015 and 2016 MET Filters
   TBranch        *b_Flag_HBHENoiseFilter;
   TBranch        *b_Flag_HBHENoiseIsoFilter;
   TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;
   TBranch        *b_Flag_goodVertices;
   TBranch        *b_Flag_eeBadScFilter;

   TBranch        *b_nGenPar;   //!
   TBranch        *b_genParEta;   //!
   TBranch        *b_genParPhi;   //!
   TBranch        *b_genParE;   //!
   TBranch        *b_genParPt;   //!
   TBranch        *b_genParId;   //!
   TBranch        *b_genParMotherId;   //!
   TBranch        *b_genParCharge;   //!
   TBranch        *b_eventRun;   //!
   TBranch        *b_eventNum;   //!
   TBranch        *b_eventLumiblock;   //!

   bool isMC;

   std::vector<int> electronIndexTight;
   std::vector<int> electronIndexLoose;
   std::vector<int> muonIndexTight;
   std::vector<int> muonIndexLoose;
   std::vector<int> jetIndex;
   std::vector<int> bTagIndex;
   std::vector<int> bTagLooseIndex;
   std::vector<int> cTagIndex;

   std::pair<TLorentzVector,TLorentzVector> zPairLeptons;
   std::pair<float,float> zPairRelIso;
   std::pair<int,int> zPairIndex;

   std::vector<float> jetSmearValue;
   std::vector<float> muonMomentumSF;

   std::pair<TLorentzVector,TLorentzVector> wPairQuarks;
   std::pair<int,int> wPairIndex;

   TLorentzVector wLepton;
   int wLepIndex;
   float wLeptonRelIso;

   Int_t numVert;
   TBranch * b_numVert;

   AnalysisEvent(bool isMC = true, std::string triggerFlag = "", TTree *tree=nullptr, bool is2016 = false, bool hasMetTriggers = false);
   virtual ~AnalysisEvent();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(bool isMC, std::string triggerFlag, TTree *tree, bool is2016, bool hasMetTriggers);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   float getEventWeight(Long64_t entry);
};

#endif

#ifdef AnalysisEvent_cxx
AnalysisEvent::AnalysisEvent(bool isMC, std::string triggerFlag, TTree *tree, bool is2016, bool hasMetTriggers) : fChain(nullptr)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == nullptr) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
       TFile *f{(TFile*)gROOT->GetListOfFiles()->FindObject("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_9_0_MJP_skim.root")};
      if (!f || !f->IsOpen()) {
        f = new TFile{"/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_9_0_MJP_skim.root"};
      }
      f->GetObject("tree",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain{new TChain{"tree",""}};
      chain->Add("/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_100_0_Gu6_skim.root/tree");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(isMC,triggerFlag,tree, is2016, hasMetTriggers);
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
   Long64_t centry{fChain->LoadTree(entry)};
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void AnalysisEvent::Init(bool isMC, std::string triggerFlag, TTree *tree, bool is2016, bool hasMetTriggers)
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
   fChain->SetBranchAddress("elePF2PATE", elePF2PATE.data(), &b_elePF2PATE);
   fChain->SetBranchAddress("elePF2PATET", elePF2PATET.data(), &b_elePF2PATET);
   fChain->SetBranchAddress("elePF2PATPX", elePF2PATPX.data(), &b_elePF2PATPX);
   fChain->SetBranchAddress("elePF2PATPY", elePF2PATPY.data(), &b_elePF2PATPY);
   fChain->SetBranchAddress("elePF2PATPZ", elePF2PATPZ.data(), &b_elePF2PATPZ);
   fChain->SetBranchAddress("elePF2PATPhi", elePF2PATPhi.data(), &b_elePF2PATPhi);
   fChain->SetBranchAddress("elePF2PATTheta", elePF2PATTheta.data(), &b_elePF2PATTheta);
   fChain->SetBranchAddress("elePF2PATEta", elePF2PATEta.data(), &b_elePF2PATEta);
   fChain->SetBranchAddress("elePF2PATPT", elePF2PATPT.data(), &b_elePF2PATPT);
   fChain->SetBranchAddress("elePF2PATCharge", elePF2PATCharge.data(), &b_elePF2PATCharge);
   if ( !is2016 ) fChain->SetBranchAddress("elePF2PATMVA", elePF2PATMVA.data(), &b_elePF2PATMVA);
   fChain->SetBranchAddress("elePF2PATCutIdVeto", elePF2PATCutIdVeto.data(), &b_elePF2PATCutIdVeto);
   fChain->SetBranchAddress("elePF2PATCutIdLoose", elePF2PATCutIdLoose.data(), &b_elePF2PATCutIdLoose);
   fChain->SetBranchAddress("elePF2PATCutIdMedium", elePF2PATCutIdMedium.data(), &b_elePF2PATCutIdMedium);
   fChain->SetBranchAddress("elePF2PATCutIdTight", elePF2PATCutIdTight.data(), &b_elePF2PATCutIdTight);
   fChain->SetBranchAddress("elePF2PATImpactTransDist", elePF2PATImpactTransDist.data(), &b_elePF2PATImpactTransDist);
   fChain->SetBranchAddress("elePF2PATImpactTransError", elePF2PATImpactTransError.data(), &b_elePF2PATImpactTransError);
   fChain->SetBranchAddress("elePF2PATImpactTransSignificance", elePF2PATImpactTransSignificance.data(), &b_elePF2PATImpactTransSignificance);
   fChain->SetBranchAddress("elePF2PATImpact3DDist", elePF2PATImpact3DDist.data(), &b_elePF2PATImpact3DDist);
   fChain->SetBranchAddress("elePF2PATImpact3DError", elePF2PATImpact3DError.data(), &b_elePF2PATImpact3DError);
   fChain->SetBranchAddress("elePF2PATImpact3DSignificance", elePF2PATImpact3DSignificance.data(), &b_elePF2PATImpact3DSignificance);
   fChain->SetBranchAddress("elePF2PATChargedHadronIso", elePF2PATChargedHadronIso.data(), &b_elePF2PATChargedHadronIso);
   fChain->SetBranchAddress("elePF2PATNeutralHadronIso", elePF2PATNeutralHadronIso.data(), &b_elePF2PATNeutralHadronIso);
   fChain->SetBranchAddress("elePF2PATPhotonIso", elePF2PATPhotonIso.data(), &b_elePF2PATPhotonIso);
   fChain->SetBranchAddress("elePF2PATTrackPt", elePF2PATTrackPt.data(), &b_elePF2PATTrackPt);
   fChain->SetBranchAddress("elePF2PATTrackPhi", elePF2PATTrackPhi.data(), &b_elePF2PATTrackPhi);
   fChain->SetBranchAddress("elePF2PATTrackEta", elePF2PATTrackEta.data(), &b_elePF2PATTrackEta);
   fChain->SetBranchAddress("elePF2PATTrackChi2", elePF2PATTrackChi2.data(), &b_elePF2PATTrackChi2);
   fChain->SetBranchAddress("elePF2PATTrackNDOF", elePF2PATTrackNDOF.data(), &b_elePF2PATTrackNDOF);
   fChain->SetBranchAddress("elePF2PATTrackD0", elePF2PATTrackD0.data(), &b_elePF2PATTrackD0);
   fChain->SetBranchAddress("elePF2PATTrackDBD0", elePF2PATTrackDBD0.data(), &b_elePF2PATTrackDBD0);
   fChain->SetBranchAddress("elePF2PATD0PV", elePF2PATD0PV.data(), &b_elePF2PATD0PV);
   fChain->SetBranchAddress("elePF2PATDZPV", elePF2PATDZPV.data(), &b_elePF2PATDZPV);
   fChain->SetBranchAddress("elePF2PATBeamSpotCorrectedTrackD0", elePF2PATBeamSpotCorrectedTrackD0.data(), &b_elePF2PATBeamSpotCorrectedTrackD0);
   fChain->SetBranchAddress("elePF2PATTrackDz", elePF2PATTrackDz.data(), &b_elePF2PATTrackDz);
   fChain->SetBranchAddress("elePF2PATVtxZ", elePF2PATVtxZ.data(), &b_elePF2PATVtxZ);
   fChain->SetBranchAddress("elePF2PATIsGsf", elePF2PATIsGsf.data(), &b_elePF2PATIsGsf);
   fChain->SetBranchAddress("elePF2PATGsfPx", elePF2PATGsfPx.data(), &b_elePF2PATGsfPx);
   fChain->SetBranchAddress("elePF2PATGsfPy", elePF2PATGsfPy.data(), &b_elePF2PATGsfPy);
   fChain->SetBranchAddress("elePF2PATGsfPz", elePF2PATGsfPz.data(), &b_elePF2PATGsfPz);
   fChain->SetBranchAddress("elePF2PATGsfE", elePF2PATGsfE.data(), &b_elePF2PATGsfE);
   fChain->SetBranchAddress("elePF2PATEcalEnergy", elePF2PATEcalEnergy.data(), &b_elePF2PATEcalEnergy);
   fChain->SetBranchAddress("elePF2PATSCEta", elePF2PATSCEta.data(), &b_elePF2PATSCEta);
   fChain->SetBranchAddress("elePF2PATSCE", elePF2PATSCE.data(), &b_elePF2PATSCE);
   fChain->SetBranchAddress("elePF2PATSCPhi", elePF2PATSCPhi.data(), &b_elePF2PATSCPhi);
   fChain->SetBranchAddress("elePF2PATSCEoverP", elePF2PATSCEoverP.data(), &b_elePF2PATSCEoverP);
   fChain->SetBranchAddress("elePF2PATSCSigmaEtaEta", elePF2PATSCSigmaEtaEta.data(), &b_elePF2PATSCSigmaEtaEta);
   fChain->SetBranchAddress("elePF2PATSCSigmaIEtaIEta", elePF2PATSCSigmaIEtaIEta.data(), &b_elePF2PATSCSigmaIEtaIEta);
   fChain->SetBranchAddress("elePF2PATSCSigmaIEtaIEta5x5", elePF2PATSCSigmaIEtaIEta5x5.data(), &b_elePF2PATSCSigmaIEtaIEta5x5);
   fChain->SetBranchAddress("elePF2PATSCE1x5", elePF2PATSCE1x5.data(), &b_elePF2PATSCE1x5);
   fChain->SetBranchAddress("elePF2PATSCE5x5", elePF2PATSCE5x5.data(), &b_elePF2PATSCE5x5);
   fChain->SetBranchAddress("elePF2PATSCE2x5max", elePF2PATSCE2x5max.data(), &b_elePF2PATSCE2x5max);
   fChain->SetBranchAddress("elePF2PATTrackIso04", elePF2PATTrackIso04.data(), &b_elePF2PATTrackIso04);
   fChain->SetBranchAddress("elePF2PATEcalIso04", elePF2PATEcalIso04.data(), &b_elePF2PATEcalIso04);
   fChain->SetBranchAddress("elePF2PATHcalIso04", elePF2PATHcalIso04.data(), &b_elePF2PATHcalIso04);
   fChain->SetBranchAddress("elePF2PATTrackIso03", elePF2PATTrackIso03.data(), &b_elePF2PATTrackIso03);
   fChain->SetBranchAddress("elePF2PATEcalIso03", elePF2PATEcalIso03.data(), &b_elePF2PATEcalIso03);
   fChain->SetBranchAddress("elePF2PATHcalIso03", elePF2PATHcalIso03.data(), &b_elePF2PATHcalIso03);
   fChain->SetBranchAddress("elePF2PATdr04EcalRecHitSumEt", elePF2PATdr04EcalRecHitSumEt.data(), &b_elePF2PATdr04EcalRecHitSumEt);
   fChain->SetBranchAddress("elePF2PATdr03EcalRecHitSumEt", elePF2PATdr03EcalRecHitSumEt.data(), &b_elePF2PATdr03EcalRecHitSumEt);
   fChain->SetBranchAddress("elePF2PATEcalIsoDeposit", elePF2PATEcalIsoDeposit.data(), &b_elePF2PATEcalIsoDeposit);
   fChain->SetBranchAddress("elePF2PATHcalIsoDeposit", elePF2PATHcalIsoDeposit.data(), &b_elePF2PATHcalIsoDeposit);
   fChain->SetBranchAddress("elePF2PATComRelIso", elePF2PATComRelIso.data(), &b_elePF2PATComRelIso);
   fChain->SetBranchAddress("elePF2PATComRelIsodBeta", elePF2PATComRelIsodBeta.data(), &b_elePF2PATComRelIsodBeta);
   fChain->SetBranchAddress("elePF2PATComRelIsoRho", elePF2PATComRelIsoRho.data(), &b_elePF2PATComRelIsoRho);
   fChain->SetBranchAddress("elePF2PATChHadIso", elePF2PATChHadIso.data(), &b_elePF2PATChHadIso);
   fChain->SetBranchAddress("elePF2PATNtHadIso", elePF2PATNtHadIso.data(), &b_elePF2PATNtHadIso);
   fChain->SetBranchAddress("elePF2PATGammaIso", elePF2PATGammaIso.data(), &b_elePF2PATGammaIso);
   fChain->SetBranchAddress("elePF2PATRhoIso", elePF2PATRhoIso.data(), &b_elePF2PATRhoIso);
   fChain->SetBranchAddress("elePF2PATAEff03", elePF2PATAEff03.data(), &b_elePF2PATAEff03);
   fChain->SetBranchAddress("elePF2PATMissingInnerLayers", elePF2PATMissingInnerLayers.data(), &b_elePF2PATMissingInnerLayers);
   fChain->SetBranchAddress("elePF2PATHoverE", elePF2PATHoverE.data(), &b_elePF2PATHoverE);
   fChain->SetBranchAddress("elePF2PATDeltaPhiSC", elePF2PATDeltaPhiSC.data(), &b_elePF2PATDeltaPhiSC);
   fChain->SetBranchAddress("elePF2PATDeltaEtaSC", elePF2PATDeltaEtaSC.data(), &b_elePF2PATDeltaEtaSC);
   fChain->SetBranchAddress("elePF2PATDeltaEtaSeedSC", elePF2PATDeltaEtaSeedSC.data(), &b_elePF2PATDeltaEtaSeedSC);
   fChain->SetBranchAddress("elePF2PATIsBarrel", elePF2PATIsBarrel.data(), &b_elePF2PATIsBarrel);
   fChain->SetBranchAddress("elePF2PATPhotonConversionTag", elePF2PATPhotonConversionTag.data(), &b_elePF2PATPhotonConversionTag);
   fChain->SetBranchAddress("elePF2PATPhotonConversionDist", elePF2PATPhotonConversionDist.data(), &b_elePF2PATPhotonConversionDist);
   fChain->SetBranchAddress("elePF2PATPhotonConversionDcot", elePF2PATPhotonConversionDcot.data(), &b_elePF2PATPhotonConversionDcot);
   fChain->SetBranchAddress("elePF2PATPhotonConversionVeto", elePF2PATPhotonConversionVeto.data(), &b_elePF2PATPhotonConversionVeto);
   fChain->SetBranchAddress("elePF2PATPhotonConversionTagCustom", elePF2PATPhotonConversionTagCustom.data(), &b_elePF2PATPhotonConversionTagCustom);
   fChain->SetBranchAddress("elePF2PATPhotonConversionDistCustom", elePF2PATPhotonConversionDistCustom.data(), &b_elePF2PATPhotonConversionDistCustom);
   fChain->SetBranchAddress("elePF2PATPhotonConversionDcotCustom", elePF2PATPhotonConversionDcotCustom.data(), &b_elePF2PATPhotonConversionDcotCustom);
   fChain->SetBranchAddress("elePF2PATTriggerMatch", elePF2PATTriggerMatch.data(), &b_elePF2PATTriggerMatch);
   fChain->SetBranchAddress("elePF2PATJetOverlap", elePF2PATJetOverlap.data(), &b_elePF2PATJetOverlap);
   if (isMC) {
     fChain->SetBranchAddress("genElePF2PATPT", genElePF2PATPT.data(), &b_genElePF2PATPT);
     fChain->SetBranchAddress("genElePF2PATET", genElePF2PATET.data(), &b_genElePF2PATET);
     fChain->SetBranchAddress("genElePF2PATPX", genElePF2PATPX.data(), &b_genElePF2PATPX);
     fChain->SetBranchAddress("genElePF2PATPY", genElePF2PATPY.data(), &b_genElePF2PATPY);
     fChain->SetBranchAddress("genElePF2PATPZ", genElePF2PATPZ.data(), &b_genElePF2PATPZ);
     fChain->SetBranchAddress("genElePF2PATPhi", genElePF2PATPhi.data(), &b_genElePF2PATPhi);
     fChain->SetBranchAddress("genElePF2PATTheta", genElePF2PATTheta.data(), &b_genElePF2PATTheta);
     fChain->SetBranchAddress("genElePF2PATEta", genElePF2PATEta.data(), &b_genElePF2PATEta);
     fChain->SetBranchAddress("genElePF2PATCharge", genElePF2PATCharge.data(), &b_genElePF2PATCharge);
     fChain->SetBranchAddress("genElePF2PATPdgId", genElePF2PATPdgId.data(), &b_genElePF2PATPdgId);
     fChain->SetBranchAddress("genElePF2PATMotherId", genElePF2PATMotherId.data(), &b_genElePF2PATMotherId);
     fChain->SetBranchAddress("genElePF2PATPromptDecayed", genElePF2PATPromptDecayed.data(), &b_genElePF2PATPromptDecayed);
     fChain->SetBranchAddress("genElePF2PATPromptFinalState", genElePF2PATPromptFinalState.data(), &b_genElePF2PATPromptFinalState);
     fChain->SetBranchAddress("genElePF2PATHardProcess", genElePF2PATHardProcess.data(), &b_genElePF2PATHardProcess);
   }
   fChain->SetBranchAddress("numMuonPF2PAT", &numMuonPF2PAT, &b_numMuonPF2PAT);
   fChain->SetBranchAddress("muonPF2PATE", muonPF2PATE.data(), &b_muonPF2PATE);
   fChain->SetBranchAddress("muonPF2PATET", muonPF2PATET.data(), &b_muonPF2PATET);
   fChain->SetBranchAddress("muonPF2PATPt", muonPF2PATPt.data(), &b_muonPF2PATPt);
   fChain->SetBranchAddress("muonPF2PATPX", muonPF2PATPX.data(), &b_muonPF2PATPX);
   fChain->SetBranchAddress("muonPF2PATPY", muonPF2PATPY.data(), &b_muonPF2PATPY);
   fChain->SetBranchAddress("muonPF2PATPZ", muonPF2PATPZ.data(), &b_muonPF2PATPZ);
   fChain->SetBranchAddress("muonPF2PATPhi", muonPF2PATPhi.data(), &b_muonPF2PATPhi);
   fChain->SetBranchAddress("muonPF2PATTheta", muonPF2PATTheta.data(), &b_muonPF2PATTheta);
   fChain->SetBranchAddress("muonPF2PATEta", muonPF2PATEta.data(), &b_muonPF2PATEta);
   fChain->SetBranchAddress("muonPF2PATCharge", muonPF2PATCharge.data(), &b_muonPF2PATCharge);
   fChain->SetBranchAddress("muonPF2PATGlobalID", muonPF2PATGlobalID.data(), &b_muonPF2PATGlobalID);
   fChain->SetBranchAddress("muonPF2PATTrackID", muonPF2PATTrackID.data(), &b_muonPF2PATTrackID);
   fChain->SetBranchAddress("muonPF2PATChi2", muonPF2PATChi2.data(), &b_muonPF2PATChi2);
   fChain->SetBranchAddress("muonPF2PATD0", muonPF2PATD0.data(), &b_muonPF2PATD0);
   fChain->SetBranchAddress("muonPF2PATTrackDBD0", muonPF2PATTrackDBD0.data(), &b_muonPF2PATTrackDBD0);
   fChain->SetBranchAddress("muonPF2PATDBInnerTrackD0", muonPF2PATDBInnerTrackD0.data(), &b_muonPF2PATDBInnerTrackD0);
   fChain->SetBranchAddress("muonPF2PATBeamSpotCorrectedD0", muonPF2PATBeamSpotCorrectedD0.data(), &b_muonPF2PATBeamSpotCorrectedD0);
   fChain->SetBranchAddress("muonPF2PATTrackNHits", muonPF2PATTrackNHits.data(), &b_muonPF2PATTrackNHits);
   fChain->SetBranchAddress("muonPF2PATMuonNHits", muonPF2PATMuonNHits.data(), &b_muonPF2PATMuonNHits);
   fChain->SetBranchAddress("muonPF2PATNDOF", muonPF2PATNDOF.data(), &b_muonPF2PATNDOF);
   fChain->SetBranchAddress("muonPF2PATVertX", muonPF2PATVertX.data(), &b_muonPF2PATVertX);
   fChain->SetBranchAddress("muonPF2PATVertY", muonPF2PATVertY.data(), &b_muonPF2PATVertY);
   fChain->SetBranchAddress("muonPF2PATVertZ", muonPF2PATVertZ.data(), &b_muonPF2PATVertZ);
   fChain->SetBranchAddress("muonPF2PATChargedHadronIso", muonPF2PATChargedHadronIso.data(), &b_muonPF2PATChargedHadronIso);
   fChain->SetBranchAddress("muonPF2PATNeutralHadronIso", muonPF2PATNeutralHadronIso.data(), &b_muonPF2PATNeutralHadronIso);
   fChain->SetBranchAddress("muonPF2PATPhotonIso", muonPF2PATPhotonIso.data(), &b_muonPF2PATPhotonIso);
   fChain->SetBranchAddress("muonPF2PATTrackIso", muonPF2PATTrackIso.data(), &b_muonPF2PATTrackIso);
   fChain->SetBranchAddress("muonPF2PATEcalIso", muonPF2PATEcalIso.data(), &b_muonPF2PATEcalIso);
   fChain->SetBranchAddress("muonPF2PATHcalIso", muonPF2PATHcalIso.data(), &b_muonPF2PATHcalIso);
   fChain->SetBranchAddress("muonPF2PATComRelIso", muonPF2PATComRelIso.data(), &b_muonPF2PATComRelIso);
   fChain->SetBranchAddress("muonPF2PATComRelIsodBeta", muonPF2PATComRelIsodBeta.data(), &b_muonPF2PATComRelIsodBeta);
   fChain->SetBranchAddress("muonPF2PATIsPFMuon", muonPF2PATIsPFMuon.data(), &b_muonPF2PATIsPFMuon);
   fChain->SetBranchAddress("muonPF2PATNChambers", muonPF2PATNChambers.data(), &b_muonPF2PATNChambers);
   fChain->SetBranchAddress("muonPF2PATNMatches", muonPF2PATNMatches.data(), &b_muonPF2PATNMatches);
   fChain->SetBranchAddress("muonPF2PATTkLysWithMeasurements", muonPF2PATTkLysWithMeasurements.data(), &b_muonPF2PATTkLysWithMeasurements);
   fChain->SetBranchAddress("muonPF2PATGlbTkNormChi2", muonPF2PATGlbTkNormChi2.data(), &b_muonPF2PATGlbTkNormChi2);
   fChain->SetBranchAddress("muonPF2PATValidFraction", muonPF2PATValidFraction.data(), &b_muonPF2PATValidFraction);
   fChain->SetBranchAddress("muonPF2PATChi2LocalPosition", muonPF2PATChi2LocalPosition.data(), &b_muonPF2PATChi2LocalPosition);
   fChain->SetBranchAddress("muonPF2PATTrkKick", muonPF2PATTrkKick.data(), &b_muonPF2PATTrkKick);
   fChain->SetBranchAddress("muonPF2PATSegmentCompatibility", muonPF2PATSegmentCompatibility.data(), &b_muonPF2PATSegmentCompatibility);
   fChain->SetBranchAddress("muonPF2PATDBPV", muonPF2PATDBPV.data(), &b_muonPF2PATDBPV);
   fChain->SetBranchAddress("muonPF2PATDZPV", muonPF2PATDZPV.data(), &b_muonPF2PATDZPV);
   fChain->SetBranchAddress("muonPF2PATVldPixHits", muonPF2PATVldPixHits.data(), &b_muonPF2PATVldPixHits);
   fChain->SetBranchAddress("muonPF2PATMatchedStations", muonPF2PATMatchedStations.data(), &b_muonPF2PATMatchedStations);
   if (isMC) {
     fChain->SetBranchAddress("genMuonPF2PATPT", genMuonPF2PATPT.data(), &b_genMuonPF2PATPT);
     fChain->SetBranchAddress("genMuonPF2PATET", genMuonPF2PATET.data(), &b_genMuonPF2PATET);
     fChain->SetBranchAddress("genMuonPF2PATPX", genMuonPF2PATPX.data(), &b_genMuonPF2PATPX);
     fChain->SetBranchAddress("genMuonPF2PATPY", genMuonPF2PATPY.data(), &b_genMuonPF2PATPY);
     fChain->SetBranchAddress("genMuonPF2PATPZ", genMuonPF2PATPZ.data(), &b_genMuonPF2PATPZ);
     fChain->SetBranchAddress("genMuonPF2PATPhi", genMuonPF2PATPhi.data(), &b_genMuonPF2PATPhi);
     fChain->SetBranchAddress("genMuonPF2PATTheta", genMuonPF2PATTheta.data(), &b_genMuonPF2PATTheta);
     fChain->SetBranchAddress("genMuonPF2PATEta", genMuonPF2PATEta.data(), &b_genMuonPF2PATEta);
     fChain->SetBranchAddress("genMuonPF2PATCharge", genMuonPF2PATCharge.data(), &b_genMuonPF2PATCharge);
     fChain->SetBranchAddress("genMuonPF2PATPdgId", genMuonPF2PATPdgId.data(), &b_genMuonPF2PATPdgId);
     fChain->SetBranchAddress("genMuonPF2PATMotherId", genMuonPF2PATMotherId.data(), &b_genMuonPF2PATMotherId);
     fChain->SetBranchAddress("genMuonPF2PATPromptDecayed", genMuonPF2PATPromptDecayed.data(), &b_genMuonPF2PATPromptDecayed);
     fChain->SetBranchAddress("genMuonPF2PATPromptFinalState", genMuonPF2PATPromptFinalState.data(), &b_genMuonPF2PATPromptFinalState);
     fChain->SetBranchAddress("genMuonPF2PATHardProcess", genMuonPF2PATHardProcess.data(), &b_genMuonPF2PATHardProcess);
   }
   fChain->SetBranchAddress("numJetPF2PAT", &numJetPF2PAT, &b_numJetPF2PAT);
   fChain->SetBranchAddress("jetPF2PATE", jetPF2PATE.data(), &b_jetPF2PATE);
   fChain->SetBranchAddress("jetPF2PATEt", jetPF2PATEt.data(), &b_jetPF2PATEt);
   fChain->SetBranchAddress("jetPF2PATPt", jetPF2PATPt.data(), &b_jetPF2PATPt);
   fChain->SetBranchAddress("jetPF2PATPtRaw", jetPF2PATPtRaw.data(), &b_jetPF2PATPtRaw);
   fChain->SetBranchAddress("jetPF2PATUnCorEt", jetPF2PATUnCorEt.data(), &b_jetPF2PATUnCorEt);
   fChain->SetBranchAddress("jetPF2PATUnCorPt", jetPF2PATUnCorPt.data(), &b_jetPF2PATUnCorPt);
   fChain->SetBranchAddress("jetPF2PATEta", jetPF2PATEta.data(), &b_jetPF2PATEta);
   fChain->SetBranchAddress("jetPF2PATTheta", jetPF2PATTheta.data(), &b_jetPF2PATTheta);
   fChain->SetBranchAddress("jetPF2PATPhi", jetPF2PATPhi.data(), &b_jetPF2PATPhi);
   fChain->SetBranchAddress("jetPF2PATPx", jetPF2PATPx.data(), &b_jetPF2PATPx);
   fChain->SetBranchAddress("jetPF2PATPy", jetPF2PATPy.data(), &b_jetPF2PATPy);
   fChain->SetBranchAddress("jetPF2PATPz", jetPF2PATPz.data(), &b_jetPF2PATPz);
   fChain->SetBranchAddress("jetPF2PATdRClosestLepton", jetPF2PATdRClosestLepton.data(), &b_jetPF2PATdRClosestLepton);
   fChain->SetBranchAddress("jetPF2PATNtracksInJet", jetPF2PATNtracksInJet.data(), &b_jetPF2PATNtracksInJet);
   fChain->SetBranchAddress("jetPF2PATJetCharge", jetPF2PATJetCharge.data(), &b_jetPF2PATJetCharge);
   fChain->SetBranchAddress("jetPF2PATfHPD", jetPF2PATfHPD.data(), &b_jetPF2PATfHPD);
   fChain->SetBranchAddress("jetPF2PATBtagSoftMuonPtRel", jetPF2PATBtagSoftMuonPtRel.data(), &b_jetPF2PATBtagSoftMuonPtRel);
   fChain->SetBranchAddress("jetPF2PATBtagSoftMuonQuality", jetPF2PATBtagSoftMuonQuality.data(), &b_jetPF2PATBtagSoftMuonQuality);
   fChain->SetBranchAddress("jetPF2PATCorrFactor", jetPF2PATCorrFactor.data(), &b_jetPF2PATCorrFactor);
   fChain->SetBranchAddress("jetPF2PATCorrResidual", jetPF2PATCorrResidual.data(), &b_jetPF2PATCorrResidual);
   fChain->SetBranchAddress("jetPF2PATL2L3ResErr", jetPF2PATL2L3ResErr.data(), &b_jetPF2PATL2L3ResErr);
   fChain->SetBranchAddress("jetPF2PATCorrErrLow", jetPF2PATCorrErrLow.data(), &b_jetPF2PATCorrErrLow);
   fChain->SetBranchAddress("jetPF2PATCorrErrHi", jetPF2PATCorrErrHi.data(), &b_jetPF2PATCorrErrHi);
   fChain->SetBranchAddress("jetPF2PATN90Hits", jetPF2PATN90Hits.data(), &b_jetPF2PATN90Hits);
   fChain->SetBranchAddress("jetPF2PATTriggered", jetPF2PATTriggered.data(), &b_jetPF2PATTriggered);
   fChain->SetBranchAddress("jetPF2PATSVX", jetPF2PATSVX.data(), &b_jetPF2PATSVX);
   fChain->SetBranchAddress("jetPF2PATSVY", jetPF2PATSVY.data(), &b_jetPF2PATSVY);
   fChain->SetBranchAddress("jetPF2PATSVZ", jetPF2PATSVZ.data(), &b_jetPF2PATSVZ);
   fChain->SetBranchAddress("jetPF2PATSVDX", jetPF2PATSVDX.data(), &b_jetPF2PATSVDX);
   fChain->SetBranchAddress("jetPF2PATSVDY", jetPF2PATSVDY.data(), &b_jetPF2PATSVDY);
   fChain->SetBranchAddress("jetPF2PATSVDZ", jetPF2PATSVDZ.data(), &b_jetPF2PATSVDZ);
   fChain->SetBranchAddress("jetPF2PATBDiscriminator", jetPF2PATBDiscriminator.data(), &b_jetPF2PATBDiscriminator);
//   fChain->SetBranchAddress("jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags", jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags.data(), &b_jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags);
//   fChain->SetBranchAddress("jetPF2PATpfCombinedMVAV2BJetTags", jetPF2PATpfCombinedMVAV2BJetTags.data(), &b_jetPF2PATpfCombinedMVAV2BJetTags);
//   fChain->SetBranchAddress("jetPF2PATpfCombinedCvsLJetTags", jetPF2PATpfCombinedCvsLJetTags.data(), &b_jetPF2PATpfCombinedCvsLJetTags);
//   fChain->SetBranchAddress("jetPF2PATpfCombinedCvsBJetTags", jetPF2PATpfCombinedCvsBJetTags.data(), &b_jetPF2PATpfCombinedCvsBJetTags);
   fChain->SetBranchAddress("jetPF2PATNConstituents", jetPF2PATNConstituents.data(), &b_jetPF2PATNConstituents);
   fChain->SetBranchAddress("jetPF2PATPID", jetPF2PATPID.data(), &b_jetPF2PATPID);
   fChain->SetBranchAddress("jetPF2PATClosestBPartonDeltaR", jetPF2PATClosestBPartonDeltaR.data(), &b_jetPF2PATClosestBPartonDeltaR);
   fChain->SetBranchAddress("jetPF2PATClosestCPartonDeltaR", jetPF2PATClosestCPartonDeltaR.data(), &b_jetPF2PATClosestCPartonDeltaR);
   if (isMC) {
     fChain->SetBranchAddress("genJetPF2PATET", genJetPF2PATET.data(), &b_genJetPF2PATET);
     fChain->SetBranchAddress("genJetPF2PATPT", genJetPF2PATPT.data(), &b_genJetPF2PATPT);
     fChain->SetBranchAddress("genJetPF2PATPX", genJetPF2PATPX.data(), &b_genJetPF2PATPX);
     fChain->SetBranchAddress("genJetPF2PATPY", genJetPF2PATPY.data(), &b_genJetPF2PATPY);
     fChain->SetBranchAddress("genJetPF2PATPZ", genJetPF2PATPZ.data(), &b_genJetPF2PATPZ);
     fChain->SetBranchAddress("genJetPF2PATPhi", genJetPF2PATPhi.data(), &b_genJetPF2PATPhi);
     fChain->SetBranchAddress("genJetPF2PATTheta", genJetPF2PATTheta.data(), &b_genJetPF2PATTheta);
     fChain->SetBranchAddress("genJetPF2PATEta", genJetPF2PATEta.data(), &b_genJetPF2PATEta);
     fChain->SetBranchAddress("genJetPF2PATPID", genJetPF2PATPID.data(), &b_genJetPF2PATPID);
   }
   fChain->SetBranchAddress("jetPF2PATMuEnergy", jetPF2PATMuEnergy.data(), &b_jetPF2PATMuEnergy);
   fChain->SetBranchAddress("jetPF2PATMuEnergyFraction", jetPF2PATMuEnergyFraction.data(), &b_jetPF2PATMuEnergyFraction);
   fChain->SetBranchAddress("jetPF2PATNeutralHadEnergy", jetPF2PATNeutralHadEnergy.data(), &b_jetPF2PATNeutralHadEnergy);
   fChain->SetBranchAddress("jetPF2PATNeutralEmEnergy", jetPF2PATNeutralEmEnergy.data(), &b_jetPF2PATNeutralEmEnergy);
   fChain->SetBranchAddress("jetPF2PATChargedHadronEnergyFraction", jetPF2PATChargedHadronEnergyFraction.data(), &b_jetPF2PATChargedHadronEnergyFraction);
   fChain->SetBranchAddress("jetPF2PATNeutralHadronEnergyFraction", jetPF2PATNeutralHadronEnergyFraction.data(), &b_jetPF2PATNeutralHadronEnergyFraction);
   fChain->SetBranchAddress("jetPF2PATChargedEmEnergyFraction", jetPF2PATChargedEmEnergyFraction.data(), &b_jetPF2PATChargedEmEnergyFraction);
   fChain->SetBranchAddress("jetPF2PATNeutralEmEnergyFraction", jetPF2PATNeutralEmEnergyFraction.data(), &b_jetPF2PATNeutralEmEnergyFraction);
   fChain->SetBranchAddress("jetPF2PATChargedHadronEnergyFractionCorr", jetPF2PATChargedHadronEnergyFractionCorr.data(), &b_jetPF2PATChargedHadronEnergyFractionCorr);
   fChain->SetBranchAddress("jetPF2PATNeutralHadronEnergyFractionCorr", jetPF2PATNeutralHadronEnergyFractionCorr.data(), &b_jetPF2PATNeutralHadronEnergyFractionCorr);
   fChain->SetBranchAddress("jetPF2PATChargedEmEnergyFractionCorr", jetPF2PATChargedEmEnergyFractionCorr.data(), &b_jetPF2PATChargedEmEnergyFractionCorr);
   fChain->SetBranchAddress("jetPF2PATNeutralEmEnergyFractionCorr", jetPF2PATNeutralEmEnergyFractionCorr.data(), &b_jetPF2PATNeutralEmEnergyFractionCorr);
   fChain->SetBranchAddress("jetPF2PATNeutralMultiplicity", jetPF2PATNeutralMultiplicity.data(), &b_jetPF2PATNeutralMultiplicity);
   fChain->SetBranchAddress("jetPF2PATChargedMultiplicity", jetPF2PATChargedMultiplicity.data(), &b_jetPF2PATChargedMultiplicity);
   fChain->SetBranchAddress("metPF2PATE", &metPF2PATE, &b_metPF2PATE);
   fChain->SetBranchAddress("metPF2PATEt", &metPF2PATEt, &b_metPF2PATEt);
   fChain->SetBranchAddress("metPF2PATEtRaw", &metPF2PATEtRaw, &b_metPF2PATEtRaw);
   fChain->SetBranchAddress("metPF2PATPhi", &metPF2PATPhi, &b_metPF2PATPhi);
   fChain->SetBranchAddress("metPF2PATPt", &metPF2PATPt, &b_metPF2PATPt);
   fChain->SetBranchAddress("metPF2PATPx", &metPF2PATPx, &b_metPF2PATPx);
   fChain->SetBranchAddress("metPF2PATPy", &metPF2PATPy, &b_metPF2PATPy);
   fChain->SetBranchAddress("metPF2PATPz", &metPF2PATPz, &b_metPF2PATPz);
   fChain->SetBranchAddress("metPF2PATScalarEt", &metPF2PATScalarEt, &b_metPF2PATScalarEt);
   fChain->SetBranchAddress("metPF2PATEtUncorrected", &metPF2PATEtUncorrected, &b_metPF2PATEtUncorrected);
   fChain->SetBranchAddress("metPF2PATPhiUncorrected", &metPF2PATPhiUncorrected, &b_metPF2PATPhiUncorrected);
   fChain->SetBranchAddress("metPF2PATUnclusteredEnUp", &metPF2PATUnclusteredEnUp, &b_metPF2PATUnclusteredEnUp);
   fChain->SetBranchAddress("metPF2PATUnclusteredEnDown", &metPF2PATUnclusteredEnDown, &b_metPF2PATUnclusteredEnDown);
   if (isMC) {
     fChain->SetBranchAddress("genMetPF2PATEt", &genMetPF2PATEt, &b_genMetPF2PATEt);
     fChain->SetBranchAddress("genMetPF2PATPhi", &genMetPF2PATPhi, &b_genMetPF2PATPhi);
     fChain->SetBranchAddress("genMetPF2PATPt", &genMetPF2PATPt, &b_genMetPF2PATPt);
     fChain->SetBranchAddress("genMetPF2PATPx", &genMetPF2PATPx, &b_genMetPF2PATPx);
     fChain->SetBranchAddress("genMetPF2PATPy", &genMetPF2PATPy, &b_genMetPF2PATPy);
     fChain->SetBranchAddress("genMetPF2PATPz", &genMetPF2PATPz, &b_genMetPF2PATPz);
   }
   fChain->SetBranchAddress("numTauPF2PAT", &numTauPF2PAT, &b_numTauPF2PAT);
   fChain->SetBranchAddress("tauPF2PATE", tauPF2PATE.data(), &b_tauPF2PATE);
   fChain->SetBranchAddress("tauPF2PATPt", tauPF2PATPt.data(), &b_tauPF2PATPt);
   fChain->SetBranchAddress("tauPF2PATPhi", tauPF2PATPhi.data(), &b_tauPF2PATPhi);
   fChain->SetBranchAddress("tauPF2PATEta", tauPF2PATEta.data(), &b_tauPF2PATEta);
   fChain->SetBranchAddress("numGeneralTracks", &numGeneralTracks, &b_numGeneralTracks);
   fChain->SetBranchAddress("generalTracksPt", generalTracksPt.data(), &b_generalTracksPt);
   fChain->SetBranchAddress("generalTracksEta", generalTracksEta.data(), &b_generalTracksEta);
   fChain->SetBranchAddress("generalTracksTheta", generalTracksTheta.data(), &b_generalTracksTheta);
   fChain->SetBranchAddress("generalTracksBeamSpotCorrectedD0", generalTracksBeamSpotCorrectedD0.data(), &b_generalTracksBeamSpotCorrectedD0);
   fChain->SetBranchAddress("generalTracksPhi", generalTracksPhi.data(), &b_generalTracksPhi);
   fChain->SetBranchAddress("generalTracksCharge", generalTracksCharge.data(), &b_generalTracksCharge);
   if (isMC) {
     fChain->SetBranchAddress("isElePlusJets", &isElePlusJets, &b_isElePlusJets);
     fChain->SetBranchAddress("genPDFScale", &genPDFScale, &b_genPDFScale);
     fChain->SetBranchAddress("genPDFx1", &genPDFx1, &b_genPDFx1);
     fChain->SetBranchAddress("genPDFx2", &genPDFx2, &b_genPDFx2);
     fChain->SetBranchAddress("genPDFf1", &genPDFf1, &b_genPDFf1);
     fChain->SetBranchAddress("genPDFf2", &genPDFf2, &b_genPDFf2);
     fChain->SetBranchAddress("topPtReweight", &topPtReweight, &b_topPtReweight);
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
   fChain->SetBranchAddress("TriggerBits", TriggerBits.data(), &b_TriggerBits);
   if (isMC) {
     fChain->SetBranchAddress("weight_muF0p5", &weight_muF0p5, &b_weight_muF0p5);
     fChain->SetBranchAddress("weight_muF2", &weight_muF2, &b_weight_muF2);
     fChain->SetBranchAddress("weight_muR0p5", &weight_muR0p5, &b_weight_muR0p5);
     fChain->SetBranchAddress("weight_muR2", &weight_muR2, &b_weight_muR2);
     fChain->SetBranchAddress("weight_muF0p5muR0p5", &weight_muF0p5muR0p5, &b_weight_muF0p5muR0p5);
     fChain->SetBranchAddress("weight_muF2muR2", &weight_muF2muR2, &b_weight_muF2muR2);
     fChain->SetBranchAddress("origWeightForNorm", &origWeightForNorm, &b_origWeightForNorm);
     fChain->SetBranchAddress("weight_pdfMax", &weight_pdfMax, &b_weight_pdfMax);
     fChain->SetBranchAddress("weight_pdfMin", &weight_pdfMin, &b_weight_pdfMin);
     fChain->SetBranchAddress("weight_alphaMax", &weight_alphaMax, &b_weight_alphaMax);
     fChain->SetBranchAddress("weight_alphaMin", &weight_alphaMin, &b_weight_alphaMin);
   }
   // 2015 triggers and MET filters
   if( !is2016 ) {
     //Data trigger branches
     fChain->SetBranchAddress("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2", &HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2, &b_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2);
     fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2", &HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2, &b_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2);
     fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v2", &HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v2, &b_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v2);
     fChain->SetBranchAddress("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3", &HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3, &b_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3);
     fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3", &HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3, &b_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3);
     fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3", &HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3, &b_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3);
     //MC trigger branches
     fChain->SetBranchAddress("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1", &HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1, &b_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1);
     fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1);
     fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1);
     fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1", &HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1, &b_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1);
     fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1", &HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1, &b_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1);
     //MET trigger branches
     if (hasMetTriggers) {
       fChain->SetBranchAddress("HLT_PFMET170_JetIdCleaned_v2", &HLT_PFMET170_JetIdCleaned_v2, &b_HLT_PFMET170_JetIdCleaned_v2);
       fChain->SetBranchAddress("HLT_PFHT350_PFMET100_v1", &HLT_PFHT350_PFMET100_v1, &b_HLT_PFHT350_PFMET100_v1);
     }
     //MET filter branches
     fChain->SetBranchAddress("Flag_CSCTightHalo2015Filter", &Flag_CSCTightHalo2015Filter, &b_Flag_CSCTightHalo2015Filter);
   }
   // 2016 triggers and MET filters
   else {
     //Lepton trigger branches
     fChain->SetBranchAddress("HLT_Ele25_eta2p1_WPTight_Gsf_v1", &HLT_Ele25_eta2p1_WPTight_Gsf_v1, &b_HLT_Ele25_eta2p1_WPTight_Gsf_v1);
     fChain->SetBranchAddress("HLT_Ele25_eta2p1_WPTight_Gsf_v2", &HLT_Ele25_eta2p1_WPTight_Gsf_v2, &b_HLT_Ele25_eta2p1_WPTight_Gsf_v2);
     fChain->SetBranchAddress("HLT_Ele25_eta2p1_WPTight_Gsf_v3", &HLT_Ele25_eta2p1_WPTight_Gsf_v3, &b_HLT_Ele25_eta2p1_WPTight_Gsf_v3);
     fChain->SetBranchAddress("HLT_Ele25_eta2p1_WPTight_Gsf_v4", &HLT_Ele25_eta2p1_WPTight_Gsf_v4, &b_HLT_Ele25_eta2p1_WPTight_Gsf_v4);
     fChain->SetBranchAddress("HLT_Ele25_eta2p1_WPTight_Gsf_v5", &HLT_Ele25_eta2p1_WPTight_Gsf_v5, &b_HLT_Ele25_eta2p1_WPTight_Gsf_v5);
     fChain->SetBranchAddress("HLT_Ele25_eta2p1_WPTight_Gsf_v6", &HLT_Ele25_eta2p1_WPTight_Gsf_v6, &b_HLT_Ele25_eta2p1_WPTight_Gsf_v6);
     fChain->SetBranchAddress("HLT_Ele25_eta2p1_WPTight_Gsf_v7", &HLT_Ele25_eta2p1_WPTight_Gsf_v7, &b_HLT_Ele25_eta2p1_WPTight_Gsf_v7);
     fChain->SetBranchAddress("HLT_Ele27_WPTight_Gsf_v1", &HLT_Ele27_WPTight_Gsf_v1, &b_HLT_Ele27_WPTight_Gsf_v1);
     fChain->SetBranchAddress("HLT_Ele27_WPTight_Gsf_v2", &HLT_Ele27_WPTight_Gsf_v2, &b_HLT_Ele27_WPTight_Gsf_v2);
     fChain->SetBranchAddress("HLT_Ele27_WPTight_Gsf_v3", &HLT_Ele27_WPTight_Gsf_v3, &b_HLT_Ele27_WPTight_Gsf_v3);
     fChain->SetBranchAddress("HLT_Ele27_WPTight_Gsf_v4", &HLT_Ele27_WPTight_Gsf_v4, &b_HLT_Ele27_WPTight_Gsf_v4);
     fChain->SetBranchAddress("HLT_Ele27_WPTight_Gsf_v5", &HLT_Ele27_WPTight_Gsf_v5, &b_HLT_Ele27_WPTight_Gsf_v5);
     fChain->SetBranchAddress("HLT_Ele27_WPTight_Gsf_v6", &HLT_Ele27_WPTight_Gsf_v6, &b_HLT_Ele27_WPTight_Gsf_v6);
     fChain->SetBranchAddress("HLT_Ele27_WPTight_Gsf_v7", &HLT_Ele27_WPTight_Gsf_v7, &b_HLT_Ele27_WPTight_Gsf_v7);
     fChain->SetBranchAddress("HLT_Ele32_eta2p1_WPTight_Gsf_v2", &HLT_Ele32_eta2p1_WPTight_Gsf_v2, &b_HLT_Ele32_eta2p1_WPTight_Gsf_v2);
     fChain->SetBranchAddress("HLT_Ele32_eta2p1_WPTight_Gsf_v3", &HLT_Ele32_eta2p1_WPTight_Gsf_v3, &b_HLT_Ele32_eta2p1_WPTight_Gsf_v3);
     fChain->SetBranchAddress("HLT_Ele32_eta2p1_WPTight_Gsf_v4", &HLT_Ele32_eta2p1_WPTight_Gsf_v4, &b_HLT_Ele32_eta2p1_WPTight_Gsf_v4);
     fChain->SetBranchAddress("HLT_Ele32_eta2p1_WPTight_Gsf_v5", &HLT_Ele32_eta2p1_WPTight_Gsf_v5, &b_HLT_Ele32_eta2p1_WPTight_Gsf_v5);
     fChain->SetBranchAddress("HLT_Ele32_eta2p1_WPTight_Gsf_v6", &HLT_Ele32_eta2p1_WPTight_Gsf_v6, &b_HLT_Ele32_eta2p1_WPTight_Gsf_v6);
     fChain->SetBranchAddress("HLT_Ele32_eta2p1_WPTight_Gsf_v7", &HLT_Ele32_eta2p1_WPTight_Gsf_v7, &b_HLT_Ele32_eta2p1_WPTight_Gsf_v7);
     fChain->SetBranchAddress("HLT_Ele32_eta2p1_WPTight_Gsf_v8", &HLT_Ele32_eta2p1_WPTight_Gsf_v8, &b_HLT_Ele32_eta2p1_WPTight_Gsf_v8);
     fChain->SetBranchAddress("HLT_IsoMu24_v1", &HLT_IsoMu24_v1, &b_HLT_IsoMu24_v1);
     fChain->SetBranchAddress("HLT_IsoMu24_v2", &HLT_IsoMu24_v2, &b_HLT_IsoMu24_v2);
     fChain->SetBranchAddress("HLT_IsoMu24_v3", &HLT_IsoMu24_v3, &b_HLT_IsoMu24_v3);
     fChain->SetBranchAddress("HLT_IsoMu24_v4", &HLT_IsoMu24_v4, &b_HLT_IsoMu24_v4);
     fChain->SetBranchAddress("HLT_IsoTkMu24_v1", &HLT_IsoTkMu24_v1, &b_HLT_IsoTkMu24_v1);
     fChain->SetBranchAddress("HLT_IsoTkMu24_v2", &HLT_IsoTkMu24_v2, &b_HLT_IsoTkMu24_v2);
     fChain->SetBranchAddress("HLT_IsoTkMu24_v3", &HLT_IsoTkMu24_v3, &b_HLT_IsoTkMu24_v3);
     fChain->SetBranchAddress("HLT_IsoTkMu24_v4", &HLT_IsoTkMu24_v4, &b_HLT_IsoTkMu24_v4);
     fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3);
     fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4);
     fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5);
     fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v6", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v6, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v6);
     fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v7", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v7, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v7);
     fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v8", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v8, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v8);
     fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v9", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v9, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v9);

     fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v2", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v2, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v2);
     fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v3", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v3, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v3);
     fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v4", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v4, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v4);
     fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v6", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v6, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v6);
     fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v2", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v2, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v2);
     fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v3", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v3, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v3);
     fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v5", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v5, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v5);
     fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2);
     fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v3", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v3, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v3);
     fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4);
     fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7);
     fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2);
     fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3);
     fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6);

     fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3);
     fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v4", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v4, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v4);
     fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v5", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v5, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v5);
     fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v6", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v6, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v6);
     fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v7", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v7, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v7);
     fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v8", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v8, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v8);
     fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v9", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v9, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v9);
     fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3);
     fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v4", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v4, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v4);
     fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v5", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v5, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v5);
     fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v6", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v6, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v6);
     fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v7", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v7, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v7);
     fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v8", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v8, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v8);
     fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v9", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v9, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v9);

//     fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v1", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v1, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v1);
//     fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3);

     fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1);
     fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2);
     fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3);
     fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4);
     fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v1", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v1, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v1);
     fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v2", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v2, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v2);
     fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v3", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v3, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v3);
     fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v4", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v4, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v4);
     //MET trigger branches
     if (hasMetTriggers) {
       fChain->SetBranchAddress("HLT_MET200_v1", &HLT_MET200_v1, &b_HLT_MET200_v1);
       fChain->SetBranchAddress("HLT_MET200_v2", &HLT_MET200_v2, &b_HLT_MET200_v2);
       fChain->SetBranchAddress("HLT_MET200_v3", &HLT_MET200_v3, &b_HLT_MET200_v3);
       fChain->SetBranchAddress("HLT_MET200_v4", &HLT_MET200_v4, &b_HLT_MET200_v4);
       fChain->SetBranchAddress("HLT_MET200_v5", &HLT_MET200_v5, &b_HLT_MET200_v5);
       fChain->SetBranchAddress("HLT_MET250_v2", &HLT_MET250_v2, &b_HLT_MET250_v2);
       fChain->SetBranchAddress("HLT_MET250_v3", &HLT_MET250_v3, &b_HLT_MET250_v3);
       fChain->SetBranchAddress("HLT_MET250_v4", &HLT_MET250_v4, &b_HLT_MET250_v4);
       fChain->SetBranchAddress("HLT_MET250_v5", &HLT_MET250_v5, &b_HLT_MET250_v5);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_v3", &HLT_PFMET120_PFMHT120_IDTight_v3, &b_HLT_PFMET120_PFMHT120_IDTight_v3);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_v4", &HLT_PFMET120_PFMHT120_IDTight_v4, &b_HLT_PFMET120_PFMHT120_IDTight_v4);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_v5", &HLT_PFMET120_PFMHT120_IDTight_v5, &b_HLT_PFMET120_PFMHT120_IDTight_v5);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_v6", &HLT_PFMET120_PFMHT120_IDTight_v6, &b_HLT_PFMET120_PFMHT120_IDTight_v6);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_v7", &HLT_PFMET120_PFMHT120_IDTight_v7, &b_HLT_PFMET120_PFMHT120_IDTight_v7);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_v8", &HLT_PFMET120_PFMHT120_IDTight_v8, &b_HLT_PFMET120_PFMHT120_IDTight_v8);
       fChain->SetBranchAddress("HLT_PFMET170_HBHECleaned_v3", &HLT_PFMET170_HBHECleaned_v3, &b_HLT_PFMET170_HBHECleaned_v3);
       fChain->SetBranchAddress("HLT_PFMET170_HBHECleaned_v4", &HLT_PFMET170_HBHECleaned_v4, &b_HLT_PFMET170_HBHECleaned_v4);
       fChain->SetBranchAddress("HLT_PFMET170_HBHECleaned_v5", &HLT_PFMET170_HBHECleaned_v5, &b_HLT_PFMET170_HBHECleaned_v5);
       fChain->SetBranchAddress("HLT_PFMET170_HBHECleaned_v6", &HLT_PFMET170_HBHECleaned_v6, &b_HLT_PFMET170_HBHECleaned_v6);
       fChain->SetBranchAddress("HLT_PFMET170_HBHECleaned_v7", &HLT_PFMET170_HBHECleaned_v7, &b_HLT_PFMET170_HBHECleaned_v7);
       fChain->SetBranchAddress("HLT_PFMET170_HBHECleaned_v8", &HLT_PFMET170_HBHECleaned_v8, &b_HLT_PFMET170_HBHECleaned_v8);
       fChain->SetBranchAddress("HLT_PFMET170_HBHECleaned_v9", &HLT_PFMET170_HBHECleaned_v9, &b_HLT_PFMET170_HBHECleaned_v9);
       fChain->SetBranchAddress("HLT_PFHT800_v3", &HLT_PFHT800_v3, &b_HLT_PFHT800_v3);
       fChain->SetBranchAddress("HLT_PFHT800_v4", &HLT_PFHT800_v4, &b_HLT_PFHT800_v4);
       fChain->SetBranchAddress("HLT_PFHT800_v5", &HLT_PFHT800_v5, &b_HLT_PFHT800_v5);
       fChain->SetBranchAddress("HLT_PFHT900_v4", &HLT_PFHT900_v4, &b_HLT_PFHT900_v4);
       fChain->SetBranchAddress("HLT_PFHT900_v5", &HLT_PFHT900_v5, &b_HLT_PFHT900_v5);
       fChain->SetBranchAddress("HLT_PFHT900_v6", &HLT_PFHT900_v6, &b_HLT_PFHT900_v6);
       fChain->SetBranchAddress("HLT_PFHT750_4JetPt50_v4", &HLT_PFHT750_4JetPt50_v4, &b_HLT_PFHT750_4JetPt50_v4);
       fChain->SetBranchAddress("HLT_PFHT750_4JetPt50_v5", &HLT_PFHT750_4JetPt50_v5, &b_HLT_PFHT750_4JetPt50_v5);
       fChain->SetBranchAddress("HLT_PFHT750_4JetPt50_v6", &HLT_PFHT750_4JetPt50_v6, &b_HLT_PFHT750_4JetPt50_v6);
       fChain->SetBranchAddress("HLT_PFHT750_4JetPt70_v1", &HLT_PFHT750_4JetPt70_v1, &b_HLT_PFHT750_4JetPt70_v1);
       fChain->SetBranchAddress("HLT_PFHT750_4JetPt70_v2", &HLT_PFHT750_4JetPt70_v2, &b_HLT_PFHT750_4JetPt70_v2);
       fChain->SetBranchAddress("HLT_PFHT750_4JetPt80_v2", &HLT_PFHT750_4JetPt80_v2, &b_HLT_PFHT750_4JetPt80_v2);
       fChain->SetBranchAddress("HLT_PFHT300_PFMET100_v1", &HLT_PFHT300_PFMET100_v1, &b_HLT_PFHT300_PFMET100_v1);
       fChain->SetBranchAddress("HLT_PFHT300_PFMET100_v2", &HLT_PFHT300_PFMET100_v2, &b_HLT_PFHT300_PFMET100_v2);
       fChain->SetBranchAddress("HLT_PFHT300_PFMET100_v3", &HLT_PFHT300_PFMET100_v3, &b_HLT_PFHT300_PFMET100_v3);
       fChain->SetBranchAddress("HLT_PFHT300_PFMET100_v4", &HLT_PFHT300_PFMET100_v4, &b_HLT_PFHT300_PFMET100_v4);
       fChain->SetBranchAddress("HLT_PFHT300_PFMET110_v4", &HLT_PFHT300_PFMET110_v4, &b_HLT_PFHT300_PFMET110_v4);
       fChain->SetBranchAddress("HLT_PFHT300_PFMET110_v5", &HLT_PFHT300_PFMET110_v5, &b_HLT_PFHT300_PFMET110_v5);
       fChain->SetBranchAddress("HLT_PFHT300_PFMET110_v6", &HLT_PFHT300_PFMET110_v6, &b_HLT_PFHT300_PFMET110_v6);
     }
     //MET filter branches
     fChain->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter, &b_Flag_globalTightHalo2016Filter);
//     fChain->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, &b_Flag_BadChargedCandidateFilter);
     fChain->SetBranchAddress("Flag_chargedHadronTrackResolutionFilter", &Flag_chargedHadronTrackResolutionFilter, &b_Flag_chargedHadronTrackResolutionFilter);
//     fChain->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
     fChain->SetBranchAddress("Flag_muonBadTrackFilter", &Flag_muonBadTrackFilter, &b_Flag_muonBadTrackFilter);
     fChain->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, &b_Flag_ecalLaserCorrFilter);
     //Feb2017 rereco filter branches
     fChain->SetBranchAddress("Flag_badMuons", &Flag_badMuons, &b_Flag_badMuons);
     fChain->SetBranchAddress("Flag_duplicateMuons", &Flag_duplicateMuons, &b_Flag_duplicateMuons);
     fChain->SetBranchAddress("Flag_noBadMuons", &Flag_noBadMuons, &b_Flag_noBadMuons);
   }

   // 2015 and 2016 triggers and MET filters
   //Lepton trigger branches
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2);
   // 2015 and 2016 MET triggers
   if (hasMetTriggers) {
     fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_v2", &HLT_PFMET120_PFMHT120_IDTight_v2, &b_HLT_PFMET120_PFMHT120_IDTight_v2);
     fChain->SetBranchAddress("HLT_PFMET170_HBHECleaned_v2", &HLT_PFMET170_HBHECleaned_v2, &b_HLT_PFMET170_HBHECleaned_v2);
     fChain->SetBranchAddress("HLT_PFHT800_v2", &HLT_PFHT800_v2, &b_HLT_PFHT800_v2);
     fChain->SetBranchAddress("HLT_MET250_v1", &HLT_MET250_v1, &b_HLT_MET250_v1);
     fChain->SetBranchAddress("HLT_PFHT750_4JetPt50_v3", &HLT_PFHT750_4JetPt50_v3, &b_HLT_PFHT750_4JetPt50_v3);
   }
   // 2015 and 2016 MET filter branches
   fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);

   if (isMC) {
     fChain->SetBranchAddress("nGenPar", &nGenPar, &b_nGenPar);
     fChain->SetBranchAddress("genParEta", genParEta.data(), &b_genParEta);
     fChain->SetBranchAddress("genParPhi", genParPhi.data(), &b_genParPhi);
     fChain->SetBranchAddress("genParE", genParE.data(), &b_genParE);
     fChain->SetBranchAddress("genParPt", genParPt.data(), &b_genParPt);
     fChain->SetBranchAddress("genParId", genParId.data(), &b_genParId);
     fChain->SetBranchAddress("genParMotherId", genParMotherId.data(), &b_genParMotherId);
     fChain->SetBranchAddress("genParCharge", genParCharge.data(), &b_genParCharge);
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
