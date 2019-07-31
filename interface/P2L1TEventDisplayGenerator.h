#ifndef P2L1TEventDisplayGenerator_H
#define P2L1TEventDisplayGenerator_H


// system include files
#include <memory>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

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
#include "FWCore/ServiceRegistry/interface/Service.h"

// GCT and RCT data formats
#include "DataFormats/L1Trigger/interface/L1PFTau.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctCollections.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"

#include "DataFormats/TauReco/interface/BaseTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauTagInfo.h"
#include "DataFormats/L1Trigger/interface/L1PFTau.h"

#include <memory>
#include <math.h>
#include <vector>
#include <list>
#include <TLorentzVector.h>

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"
#include "CondFormats/DataRecord/interface/L1CaloHcalScaleRcd.h"
#include "FWCore/Framework/interface/ESHandle.h"
//Vertex and gen particle
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"



#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
//#include "SimDataFormats/CaloTest/interface/HcalTestNumbering.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalSourcePositionData.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"

//track trigger data formats
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"

#include "L1Trigger/phase2L1EcalTimingTP/plugins/helpers.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

using std::vector;
using namespace l1t;
// class declaration
//

class P2L1TEventDisplayGenerator : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  typedef vector<reco::GenParticle> GenParticleCollectionType;

 public:
  
  // Constructor
  P2L1TEventDisplayGenerator(const edm::ParameterSet& ps);
  
  // Destructor
  virtual ~P2L1TEventDisplayGenerator();



  std::vector<double> *hcalTpgs_Pt  = new std::vector<double>; 
  std::vector<double> *hcalTpgs_Eta = new std::vector<double>; 
  std::vector<double> *hcalTpgs_Phi = new std::vector<double>; 

  std::vector<double> *ecalTpgs_Pt  = new std::vector<double>; 
  std::vector<double> *ecalTpgs_Eta = new std::vector<double>; 
  std::vector<double> *ecalTpgs_Phi = new std::vector<double>; 

  std::vector<double> *ecalCrys_Pt  = new std::vector<double>; 
  std::vector<double> *ecalCrys_Eta = new std::vector<double>; 
  std::vector<double> *ecalCrys_Phi = new std::vector<double>; 

  std::vector<double> *sumTpgs_Pt  = new std::vector<double>; 
  std::vector<double> *sumTpgs_Eta = new std::vector<double>; 
  std::vector<double> *sumTpgs_Phi = new std::vector<double>; 

  /*
  std::vector<TLorentzVector> *rlxTaus  = new std::vector<TLorentzVector>; 
  std::vector<TLorentzVector> *isoTaus  = new std::vector<TLorentzVector>; 
  std::vector<TLorentzVector> *recoTaus  = new std::vector<TLorentzVector>; 
  std::vector<TLorentzVector> *allRegions  = new std::vector<TLorentzVector>; 
  std::vector<TLorentzVector> *signalPFCands  = new std::vector<TLorentzVector>; 
  std::vector<TLorentzVector> *l1Jets  = new std::vector<TLorentzVector>; 
  std::vector<TLorentzVector> *recoJets  = new std::vector<TLorentzVector>; 
  std::vector<TLorentzVector> *caloClusters  = new std::vector<TLorentzVector>; 
  std::vector<double> *recoJetsDr  = new std::vector<double>; 
  */
  std::vector<TLorentzVector> *allEcalTPGs            = new std::vector<TLorentzVector>; 
  std::vector<TLorentzVector> *allHcalTPGs            = new std::vector<TLorentzVector>; 
  std::vector<TLorentzVector> *l1Tracks               = new std::vector<TLorentzVector>; 
  std::vector<TLorentzVector> *l1EcalClusters         = new std::vector<TLorentzVector>; 
  std::vector<TLorentzVector> *l1PFTaus               = new std::vector<TLorentzVector>; 
  std::vector<TLorentzVector> *l1EcalCrystals         = new std::vector<TLorentzVector>; 
  std::vector<TLorentzVector> *genHadronicTaus        = new std::vector<TLorentzVector>; 

  TH1F* isoTau_pt;
  TH1F* isoTau_eta;
  TH1F* isoTau_phi;

  TH1F* tau_pt;
  TH1F* tau_eta;
  TH1F* tau_phi;

  TH1F* recoTau_pt;
  TH1F* recoTau_eta;
  TH1F* recoTau_phi;
  TTree* efficiencyTree;
  TTree* efficiencyTreePiPM;
  TTree* efficiencyTreePiZero;

  int run, lumi, event;
  double isoTauPt, rlxTauPt, isoTauEta, rlxTauEta, isoTauPhi, rlxTauPhi;
  double recoPt, recoEta, recoPhi;
  double genPt, genEta, genPhi;
  double genPizeroPt, genPizeroEta, genPizeroPhi;
  int l1RlxMatched, l1IsoMatched;
  int decayMode;
  double tauEtaEcalEnt,tauPhiEcalEnt,rawEcal, rawHcal, ecal, hcal, jetEt, jetEta, jetPhi, nvtx;
  double max3ProngDeltaR, minProngPt, maxProngPt, midProngPt; int n3ProngCands;
  double pfCandsEt, signalCandsEt, isoCandsEt;
  double TPG2x2, TPGH2x2, TPGE2x2;
  double TPG5x5, TPGH5x5, TPGE5x5;
  double TPG6x6, TPGH6x6, TPGE6x6;
  double TPG7x7, TPGH7x7, TPGE7x7;



 protected:
  // Analyze
  void analyze(const edm::Event& evt, const edm::EventSetup& es);
  
  // BeginJob
  void beginJob(const edm::EventSetup &es);
  
  // EndJob
  void endJob(void);

  
 private:
  // ----------member data ---------------------------

  typedef std::vector<TTTrack<Ref_Phase2TrackerDigi_>> L1TkTrackCollectionType;

  int nev_; // Number of events processed
  bool verbose_;
  std::ofstream logFile_;

  edm::InputTag L1TrackInputTag;
  edm::InputTag L1TrackPrimaryVertexTag;

  bool debug; 
  edm::EDGetTokenT<EcalTrigPrimDigiCollection> ecalSrc_; 
  edm::EDGetTokenT<HcalTrigPrimDigiCollection> hcalSrc_;
  edm::EDGetTokenT<reco::VertexCollection> vtxLabel_;
  edm::EDGetTokenT<EcalEBTrigPrimDigiCollection> ecalTPGBToken_;
  edm::EDGetTokenT< l1t::L1PFTauCollection > L1PFTauToken_;
  edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken_;
  edm::InputTag genSrc_;
  double genMatchDeltaRcut;
  edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > ttTrackToken_;

  std::string folderName_;
  double recoPt_;

  edm::ESHandle<CaloGeometry> caloGeometry_;
  const CaloSubdetectorGeometry * ebGeometry;

  void getEcalCrystals(edm::Handle<EcalEBTrigPrimDigiCollection> ecaltpgs, 
		       vector<TLorentzVector> &ecalCrystals);

		 
 int TPGEtaRange(int ieta){
   int iEta = 0;
   // So here, -28 becomes 0.  -1 be comes 27.  +1 becomes 28. +28 becomes 55.
   // And we have mapped [-28, -1], [1, 28] onto [0, 55]   
   if(ieta < 0)
     iEta = ieta + 28;
   else if(ieta > 0)
     iEta = ieta + 27;
   if(ieta==0){
     std::cout<<"Error! ieta is 0, ieta: "<<ieta<<" iEta "<<iEta<<std::endl;
     exit(0);
   }
   return iEta;
 }

  int convertGenEta(double inputEta) {
    const double tpgEtaValues[27] = {
      0.087,      
      0.174, // HB and inner HE bins are 0.348 wide
      0.261,
      0.348,
      0.522,
      0.609,
      0.696,
      0.783,
      0.870,
      0.957,
      1.044,
      1.131,
      1.218,
      1.305,
      1.392,
      1.479,
      1.566,
      1.653,
      1.74,
      1.848,
      1.956, // Last two HE bins are 0.432 and 0.828 wide
      2.064,
      2.172,
      2.379,
      2.586,
      2.793,
      3
      //IGNORING HF
      //3.250, // HF bins are 0.5 wide
      //3.750,
      //4.250,
      //4.750
    };

    for (int n=1; n<29; n++){
      //std::cout<<"inputEta "<<inputEta<< " n "<< n <<" tpgEtaValues[n-1] "<< tpgEtaValues[n-1] << " abs(inputEta)<tpgEtaValues[n-1]"<<std::endl;
      if (std::fabs(inputEta)<tpgEtaValues[n-1]) {
	//Positive eta is >28 negative eta is 0 to 27
	if(inputEta>0){ return n + 28;}
	else{ return n;}
	break;
      }
    }
    std::cout<<"OUT OF BOUNDS!!!!  inputeta: "<<inputEta<<std::endl;
    return -9;
  }

  //-pi < phi <= +pi,
  int convertGenPhi(double inputPhi){
    double posPhi[36];
    for(int n = 0; n < 36; n++)
      posPhi[n] = (0.087) * n + 0.0435;
    double negPhi[36];
    for(int n = 0; n < 36; n++)
      negPhi[n] = -3.14159 + 0.087 * n - 0.0435;

    //1 to 36 is 0 to pi
    if( 3.1416 > inputPhi && inputPhi >= 0){

      for(int n = 1; n < 36; n++){
	//std::cout<<"inputPhi "<<inputPhi<< " posPhi[n-1] "<< posPhi[n-1] << " n "<<n<<std::endl;
	if(inputPhi <= posPhi[n-1]){
	  int tpgPhi = n;
	  return tpgPhi;
	}
      }
    }

    //37 to 72 is -pi to 0
    else if(-3.1416 < inputPhi && inputPhi < 0){
      for(int n = 1; n < 36; n++)
	if(inputPhi < negPhi[n-1]){
	  int tpgPhi = n + 36;
	  return tpgPhi;
	}
    }
    std::cout<<"OUT OF BOUNDS!!!!  inputphi: "<<inputPhi<<std::endl;
    return -9;
  }
  //map of tpg eta values by ieta
  double tpgEtaMap[57] = {0};
  //map of tpg phi values by iphi
  double tpgPhiMap[73] = {0};
  //fill gen eta map and gen phi map
  void initializeEtaPhiMaps(){
    for(int i = 0; i <57; i++){
      tpgEtaMap[i] = convertGenEta(i);
    }
    for(int i = 0; i <73; i++){
      tpgPhiMap[i] = convertGenPhi(i);
    }
  }

  void clearThePlottingVectors(){
   //Clear the vectors
    l1Tracks->clear(); 
    l1EcalClusters->clear(); 
    l1EcalCrystals->clear(); 
    l1PFTaus->clear();
    genHadronicTaus->clear(); 
    allHcalTPGs->clear();
  }


};

#endif
