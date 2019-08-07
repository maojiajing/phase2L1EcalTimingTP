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
#define GENPARTICLEARRAYSIZE 2000
#define EBCRYSTALARRAYSIZE 70000

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

// ------------ reset branches  ------------
  virtual void resetBranches();
  void resetEventInfoBranches();
  void resetEBCrystalBranches();
  void resetGenParticleBranches();

// ------------corr eta phi  ------------
  vector<float> EtaPhi_Corr_EB(float X, float Y, float Z, reco::GenParticle gen);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------
  typedef std::vector<reco::GenParticle> GenParticleCollectionType;


  const CaloSubdetectorGeometry * ebGeometry;
  edm::ESHandle<CaloGeometry> caloGeometry_;

  edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken_;
  edm::EDGetTokenT<float> genTokenT_;
  edm::EDGetTokenT<EcalEBTrigPrimDigiCollection> ecalTPGBToken_;

  edm::Handle<EcalEBTrigPrimDigiCollection> ecaltpgCollection;
  edm::Handle<GenParticleCollectionType> genParticleHandle;
  edm::Handle<float> genVertexTHandle;

  edm::InputTag genSrc_;
  edm::InputTag genSrcT_;

  TTree* ecalTPTree;

  unsigned int runNum;
  unsigned int lumiNum;
  unsigned int eventNum;

 //ecal crystal info
 int nCrystals; 

 float eb_Et[EBCRYSTALARRAYSIZE];
 float eb_time[EBCRYSTALARRAYSIZE];
 int eb_id[EBCRYSTALARRAYSIZE];
 int eb_ieta[EBCRYSTALARRAYSIZE];
 int eb_iphi[EBCRYSTALARRAYSIZE];
 int eb_ism[EBCRYSTALARRAYSIZE];
 int eb_ic[EBCRYSTALARRAYSIZE];
 float eb_cell_Eta[EBCRYSTALARRAYSIZE];
 float eb_cell_Phi[EBCRYSTALARRAYSIZE];

 //gen info
 int nGenParticles;

 float genVertexX;
 float genVertexY;
 float genVertexZ;
 float genVertexT;

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

 float gE10x10[GENPARTICLEARRAYSIZE];
 float gE5x5[GENPARTICLEARRAYSIZE];
 float gE3x3[GENPARTICLEARRAYSIZE];

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

