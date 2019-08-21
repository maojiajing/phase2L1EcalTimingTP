// -*- C++ -*-
//
// Package:    L1Trigger/phase2L1EcalTimingAnalyzer
// Class:      phase2L1EcalTimingAnalyzer
// 
/**\class phase2L1EcalTimingAnalyzer phase2L1EcalTimingAnalyzer.cc L1Trigger/phase2L1EcalTimingTP/plugins/phase2L1EcalTimingAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jiajing Mao
//
//
#define GENPARTICLEARRAYSIZE 2500
#define EBCRYSTALARRAYSIZE 70000
#define GENJETARRAYSIZE 2000

// system include files
#include <memory>

// Math Include
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include <TLorentzVector.h>
#include <memory>
#include <math.h>
#include <vector>
#include <list>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

//#include "L1Trigger/phase2Demonstrator/interface/triggerGeometryTools.hh"
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"

//Vertex and gen particle
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"


#include "L1Trigger/phase2L1EcalTimingTP/plugins/helpers.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/L1Trigger/interface/L1PFTau.h"

#include "DataFormats/Candidate/interface/Candidate.h"

//#include "DataFormats/L1Trigger/interface/L1PFObject.h"

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/L1TrackTrigger/interface/L1TkPrimaryVertex.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

using namespace l1t;
using namespace std;

class phase2L1EcalTimingAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit phase2L1EcalTimingAnalyzer(const edm::ParameterSet&);
  ~phase2L1EcalTimingAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

// ------------ load event  ------------
  void loadEvent(const edm::Event& iEvent);

// ------------ enable branches  ------------
  virtual void setBranches();
  void enableEventInfoBranches();
  void enableEBCrystalBranches();
  void enableGenParticleBranches();
  void enableGenak4JetBranches();
  void enableGenak4JetNoNuBranches();
  void enableGenak8JetBranches();
  void enableGenak8JetNoNuBranches();

// ------------ reset branches  ------------
  virtual void resetBranches();
  void resetEventInfoBranches();
  void resetEBCrystalBranches();
  void resetGenParticleBranches();
  void resetGenak4JetBranches();
  void resetGenak4JetNoNuBranches();
  void resetGenak8JetBranches();
  void resetGenak8JetNoNuBranches();

// ------------ fill branches  ------------
  bool fillEventInfoBranches(const edm::Event& iEvent);
  bool fillEBCrystalBranches(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  bool fillGenParticleBranches();
  //std::vector<reco::GenParticle> GetGenParticles();
  //bool fillGenParticleBasicBranches(std::vector<reco::GenParticle> genParticles);
  //bool fillGenParticleMotherBranches(std::vector<reco::GenParticle> genParticles);
  //bool fillGenParticleGrandMotherBranches(std::vector<reco::GenParticle> genParticles);
  //bool fillGenParticleSiblingBranches(std::vector<reco::GenParticle> genParticles);
  //bool fillGenParticleTPBranches(std::vector<reco::GenParticle> genParticles);

  bool fillGenak4JetBranches();
  bool fillGenak4JetNoNuBranches();
  bool fillGenak8JetBranches();
  bool fillGenak8JetNoNuBranches();

// ------------corr eta phi  ------------
  vector<float> EtaPhi_Corr_EB(float X, float Y, float Z, reco::GenParticle gen);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------
  typedef std::vector<reco::GenParticle> GenParticleCollectionType;
  typedef std::vector<reco::GenJet> GenJetCollectionType;


  const CaloSubdetectorGeometry * ebGeometry;
  edm::ESHandle<CaloGeometry> caloGeometry_;

  edm::EDGetTokenT<std::vector<reco::GenJet> > genak4JetToken_;
  edm::EDGetTokenT<std::vector<reco::GenJet> > genak4JetNoNuToken_;
  edm::EDGetTokenT<std::vector<reco::GenJet> > genak8JetToken_;
  edm::EDGetTokenT<std::vector<reco::GenJet> > genak8JetNoNuToken_;
  edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken_;
  edm::EDGetTokenT<float> genTokenT_;
  edm::EDGetTokenT<EcalEBTrigPrimDigiCollection> ecalTPGBToken_;

  edm::Handle<EcalEBTrigPrimDigiCollection> ecaltpgCollection;
  edm::Handle<GenJetCollectionType> genak4JetHandle;
  edm::Handle<GenJetCollectionType> genak4JetNoNuHandle;
  edm::Handle<GenJetCollectionType> genak8JetHandle;
  edm::Handle<GenJetCollectionType> genak8JetNoNuHandle;
  edm::Handle<GenParticleCollectionType> genParticleHandle;
  edm::Handle<float> genVertexTHandle;

  edm::InputTag genSrcak4J_;
  edm::InputTag genSrcak4JN_;
  edm::InputTag genSrcak8J_;
  edm::InputTag genSrcak8JN_;
  edm::InputTag genSrc_;
  edm::InputTag genSrcT_;

  TTree* ecalTPTree;

  unsigned int runNum;
  unsigned int lumiNum;
  unsigned int eventNum;

 //ecal crystal info
 int nCrystals; 

 float eb_Et[EBCRYSTALARRAYSIZE];
 float eb_Edep[EBCRYSTALARRAYSIZE];
 float eb_time[EBCRYSTALARRAYSIZE];
 float eb_sigmat[EBCRYSTALARRAYSIZE];
 int eb_id[EBCRYSTALARRAYSIZE];
 int eb_ieta[EBCRYSTALARRAYSIZE];
 int eb_iphi[EBCRYSTALARRAYSIZE];
 int eb_ism[EBCRYSTALARRAYSIZE];
 int eb_ic[EBCRYSTALARRAYSIZE];
 float eb_cell_Eta[EBCRYSTALARRAYSIZE];
 float eb_cell_Phi[EBCRYSTALARRAYSIZE];

 //ak4 jet info
 int nGenak4Jets;

 float gak4JetMass[GENJETARRAYSIZE];
 float gak4JetE[GENJETARRAYSIZE];
 float gak4JetEt[GENJETARRAYSIZE];
 float gak4JetPt[GENJETARRAYSIZE];
 float gak4JetPx[GENJETARRAYSIZE];
 float gak4JetPy[GENJETARRAYSIZE];
 float gak4JetPz[GENJETARRAYSIZE];
 float gak4JetEta[GENJETARRAYSIZE];
 float gak4JetPhi[GENJETARRAYSIZE];

 float gak4JetArea[GENJETARRAYSIZE];

 float gak4JetPileupE[GENJETARRAYSIZE];
 int gak4JetPileupIdFlag[GENJETARRAYSIZE];

 bool gak4JetPassIdLoose[GENJETARRAYSIZE];
 bool gak4JetPassIdTight[GENJETARRAYSIZE];

 float gak4JetMuEnergy[GENJETARRAYSIZE];
 float gak4JetEmEnergy[GENJETARRAYSIZE];
 float gak4JetChargedEmEnergy[GENJETARRAYSIZE];
 float gak4JetNeutralEmEnergy[GENJETARRAYSIZE];
 float gak4JetHadronEnergy[GENJETARRAYSIZE];
 float gak4JetChargedHadronEnergy[GENJETARRAYSIZE];
 float gak4JetNeutralHadronEnergy[GENJETARRAYSIZE];

 // ak4 jet nonu info
 int nGenak4JetNoNus;

 float gak4JetNoNuMass[GENJETARRAYSIZE];
 float gak4JetNoNuE[GENJETARRAYSIZE];
 float gak4JetNoNuEt[GENJETARRAYSIZE];
 float gak4JetNoNuPt[GENJETARRAYSIZE];
 float gak4JetNoNuPx[GENJETARRAYSIZE];
 float gak4JetNoNuPy[GENJETARRAYSIZE];
 float gak4JetNoNuPz[GENJETARRAYSIZE];
 float gak4JetNoNuEta[GENJETARRAYSIZE];
 float gak4JetNoNuPhi[GENJETARRAYSIZE];

 float gak4JetNoNuArea[GENJETARRAYSIZE];

 float gak4JetNoNuPileupE[GENJETARRAYSIZE];
 int gak4JetNoNuPileupIdFlag[GENJETARRAYSIZE];

 bool gak4JetNoNuPassIdLoose[GENJETARRAYSIZE];
 bool gak4JetNoNuPassIdTight[GENJETARRAYSIZE];

 float gak4JetNoNuMuEnergy[GENJETARRAYSIZE];
 float gak4JetNoNuEmEnergy[GENJETARRAYSIZE];
 float gak4JetNoNuChargedEmEnergy[GENJETARRAYSIZE];
 float gak4JetNoNuNeutralEmEnergy[GENJETARRAYSIZE];
 float gak4JetNoNuHadronEnergy[GENJETARRAYSIZE];
 float gak4JetNoNuChargedHadronEnergy[GENJETARRAYSIZE];
 float gak4JetNoNuNeutralHadronEnergy[GENJETARRAYSIZE];

 // ak8 jet info
 int nGenak8Jets;

 float gak8JetMass[GENJETARRAYSIZE];
 float gak8JetE[GENJETARRAYSIZE];
 float gak8JetEt[GENJETARRAYSIZE];
 float gak8JetPt[GENJETARRAYSIZE];
 float gak8JetPx[GENJETARRAYSIZE];
 float gak8JetPy[GENJETARRAYSIZE];
 float gak8JetPz[GENJETARRAYSIZE];
 float gak8JetEta[GENJETARRAYSIZE];
 float gak8JetPhi[GENJETARRAYSIZE];

 float gak8JetArea[GENJETARRAYSIZE];

 float gak8JetPileupE[GENJETARRAYSIZE];
 int gak8JetPileupIdFlag[GENJETARRAYSIZE];

 bool gak8JetPassIdLoose[GENJETARRAYSIZE];
 bool gak8JetPassIdTight[GENJETARRAYSIZE];

 float gak8JetMuEnergy[GENJETARRAYSIZE];
 float gak8JetEmEnergy[GENJETARRAYSIZE];
 float gak8JetChargedEmEnergy[GENJETARRAYSIZE];
 float gak8JetNeutralEmEnergy[GENJETARRAYSIZE];
 float gak8JetHadronEnergy[GENJETARRAYSIZE];
 float gak8JetChargedHadronEnergy[GENJETARRAYSIZE];
 float gak8JetNeutralHadronEnergy[GENJETARRAYSIZE];

 // ak8 jet nonu info
 int nGenak8JetNoNus;

 float gak8JetNoNuMass[GENJETARRAYSIZE];
 float gak8JetNoNuE[GENJETARRAYSIZE];
 float gak8JetNoNuEt[GENJETARRAYSIZE];
 float gak8JetNoNuPt[GENJETARRAYSIZE];
 float gak8JetNoNuPx[GENJETARRAYSIZE];
 float gak8JetNoNuPy[GENJETARRAYSIZE];
 float gak8JetNoNuPz[GENJETARRAYSIZE];
 float gak8JetNoNuEta[GENJETARRAYSIZE];
 float gak8JetNoNuPhi[GENJETARRAYSIZE];

 float gak8JetNoNuArea[GENJETARRAYSIZE];

 float gak8JetNoNuPileupE[GENJETARRAYSIZE];
 int gak8JetNoNuPileupIdFlag[GENJETARRAYSIZE];

 bool gak8JetNoNuPassIdLoose[GENJETARRAYSIZE];
 bool gak8JetNoNuPassIdTight[GENJETARRAYSIZE];

 float gak8JetNoNuMuEnergy[GENJETARRAYSIZE];
 float gak8JetNoNuEmEnergy[GENJETARRAYSIZE];
 float gak8JetNoNuChargedEmEnergy[GENJETARRAYSIZE];
 float gak8JetNoNuNeutralEmEnergy[GENJETARRAYSIZE];
 float gak8JetNoNuHadronEnergy[GENJETARRAYSIZE];
 float gak8JetNoNuChargedHadronEnergy[GENJETARRAYSIZE];
 float gak8JetNoNuNeutralHadronEnergy[GENJETARRAYSIZE];

 //gen info
 int nGenParticles;

 float genVertexX;
 float genVertexY;
 float genVertexZ;
 float genVertexT;

 int gParticleGrandMotherId[GENPARTICLEARRAYSIZE];
 int gParticleGrandMotherIndex[GENPARTICLEARRAYSIZE];

 int gParticleMotherId[GENPARTICLEARRAYSIZE];
 int gParticleMotherIndex[GENPARTICLEARRAYSIZE];

 int gParticleSiblingId[GENPARTICLEARRAYSIZE];
 int gParticleSiblingIndex[GENPARTICLEARRAYSIZE];

 int gParticleId[GENPARTICLEARRAYSIZE];
 int gParticleStatus[GENPARTICLEARRAYSIZE];

 float gParticleE[GENPARTICLEARRAYSIZE];
 float gParticlePt[GENPARTICLEARRAYSIZE];
 float gParticlePx[GENPARTICLEARRAYSIZE];
 float gParticlePy[GENPARTICLEARRAYSIZE];
 float gParticlePz[GENPARTICLEARRAYSIZE];
 float gParticleEta[GENPARTICLEARRAYSIZE];
 float gParticlePhi[GENPARTICLEARRAYSIZE];

 float gParticle_decay_vtx_x[GENPARTICLEARRAYSIZE];
 float gParticle_decay_vtx_y[GENPARTICLEARRAYSIZE];
 float gParticle_decay_vtx_z[GENPARTICLEARRAYSIZE];

 float gParticle_prod_vtx_x[GENPARTICLEARRAYSIZE];
 float gParticle_prod_vtx_y[GENPARTICLEARRAYSIZE];
 float gParticle_prod_vtx_z[GENPARTICLEARRAYSIZE];

 float g_eb_time_Emax_02[GENPARTICLEARRAYSIZE];
 float g_eb_sigmat_Emax_02[GENPARTICLEARRAYSIZE];
 float gEmax_02[GENPARTICLEARRAYSIZE];
 int gImax_02[GENPARTICLEARRAYSIZE];
 float g_dt_sc_02[GENPARTICLEARRAYSIZE];
 float g_tmax_sc_02[GENPARTICLEARRAYSIZE];
 int g_cnt_sc_02[GENPARTICLEARRAYSIZE];
 float g_tmin_sc_02[GENPARTICLEARRAYSIZE];
 float g_dt_sm_02[GENPARTICLEARRAYSIZE];
 float g_tmax_sm_02[GENPARTICLEARRAYSIZE];
 int g_cnt_sm_02[GENPARTICLEARRAYSIZE];
 float g_tmin_sm_02[GENPARTICLEARRAYSIZE];
 float gEsc_02[GENPARTICLEARRAYSIZE];
 float genEsc_02[GENPARTICLEARRAYSIZE];
 float gE9x9_02[GENPARTICLEARRAYSIZE];
 float gE5x5_02[GENPARTICLEARRAYSIZE];
 float gE3x3_02[GENPARTICLEARRAYSIZE];

 float g_eb_time_Emax_01[GENPARTICLEARRAYSIZE];
 float g_eb_sigmat_Emax_01[GENPARTICLEARRAYSIZE];
 float gEmax_01[GENPARTICLEARRAYSIZE];
 int gImax_01[GENPARTICLEARRAYSIZE];
 float g_dt_sc_01[GENPARTICLEARRAYSIZE];
 float g_tmax_sc_01[GENPARTICLEARRAYSIZE];
 int g_cnt_sc_01[GENPARTICLEARRAYSIZE];
 float g_tmin_sc_01[GENPARTICLEARRAYSIZE];
 float g_dt_max_01[GENPARTICLEARRAYSIZE];
 float g_t_max_01[GENPARTICLEARRAYSIZE];
 float g_e_max_01[GENPARTICLEARRAYSIZE];
 float g_dt_min_01[GENPARTICLEARRAYSIZE];
 float g_t_min_01[GENPARTICLEARRAYSIZE];
 float g_e_min_01[GENPARTICLEARRAYSIZE];
 float g_dt_sm_01[GENPARTICLEARRAYSIZE];
 float g_tmax_sm_01[GENPARTICLEARRAYSIZE];
 int g_cnt_sm_01[GENPARTICLEARRAYSIZE];
 float g_tmin_sm_01[GENPARTICLEARRAYSIZE];
 float g_dt_sm2_01[GENPARTICLEARRAYSIZE];
 float g_tmax_sm2_01[GENPARTICLEARRAYSIZE];
 int g_cnt_sm2_01[GENPARTICLEARRAYSIZE];
 float g_tmin_sm2_01[GENPARTICLEARRAYSIZE];
 float g_dt_sm5_01[GENPARTICLEARRAYSIZE];
 float g_tmax_sm5_01[GENPARTICLEARRAYSIZE];
 int g_cnt_sm5_01[GENPARTICLEARRAYSIZE];
 float g_tmin_sm5_01[GENPARTICLEARRAYSIZE];
 float g_dt_sm10_01[GENPARTICLEARRAYSIZE];
 float g_tmax_sm10_01[GENPARTICLEARRAYSIZE];
 int g_cnt_sm10_01[GENPARTICLEARRAYSIZE];
 float g_tmin_sm10_01[GENPARTICLEARRAYSIZE];
 float gEsc_01[GENPARTICLEARRAYSIZE];
 float genEsc_01[GENPARTICLEARRAYSIZE];
 float gE9x9_01[GENPARTICLEARRAYSIZE];
 float gE5x5_01[GENPARTICLEARRAYSIZE];
 float gE3x3_01[GENPARTICLEARRAYSIZE];

 int gIcore[GENPARTICLEARRAYSIZE];
 float g_eb_time[GENPARTICLEARRAYSIZE];
 float g_eb_sigmat[GENPARTICLEARRAYSIZE];
 float g_tof[GENPARTICLEARRAYSIZE];
 float g_tvirtual[GENPARTICLEARRAYSIZE];
 float gen_time[GENPARTICLEARRAYSIZE];
 float gen_time_max_02[GENPARTICLEARRAYSIZE];
 float gen_time_max_01[GENPARTICLEARRAYSIZE];

 float gParticleGrandMotherE[GENPARTICLEARRAYSIZE];
 float gParticleGrandMotherPt[GENPARTICLEARRAYSIZE];
 float gParticleGrandMotherPx[GENPARTICLEARRAYSIZE];
 float gParticleGrandMotherPy[GENPARTICLEARRAYSIZE];
 float gParticleGrandMotherPz[GENPARTICLEARRAYSIZE];
 float gParticleGrandMotherEta[GENPARTICLEARRAYSIZE];
 float gParticleGrandMotherPhi[GENPARTICLEARRAYSIZE];
 float gParticleGrandMotherDR[GENPARTICLEARRAYSIZE];

 float gParticleMotherE[GENPARTICLEARRAYSIZE];
 float gParticleMotherPt[GENPARTICLEARRAYSIZE];
 float gParticleMotherPx[GENPARTICLEARRAYSIZE];
 float gParticleMotherPy[GENPARTICLEARRAYSIZE];
 float gParticleMotherPz[GENPARTICLEARRAYSIZE];
 float gParticleMotherEta[GENPARTICLEARRAYSIZE];
 float gParticleMotherPhi[GENPARTICLEARRAYSIZE];
 float gParticleMotherDR[GENPARTICLEARRAYSIZE];

 float gParticleSiblingE[GENPARTICLEARRAYSIZE];
 float gParticleSiblingPt[GENPARTICLEARRAYSIZE];
 float gParticleSiblingPx[GENPARTICLEARRAYSIZE];
 float gParticleSiblingPy[GENPARTICLEARRAYSIZE];
 float gParticleSiblingPz[GENPARTICLEARRAYSIZE];
 float gParticleSiblingEta[GENPARTICLEARRAYSIZE];
 float gParticleSiblingPhi[GENPARTICLEARRAYSIZE];
 float gParticleSiblingDR[GENPARTICLEARRAYSIZE];

};

