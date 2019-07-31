/*
 * \file P2L1TEventDisplayGenerator.cc
 *
 * \author I. Ojalvo
 * Written for Gen
 */
#include <TLorentzVector.h>
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/Math/interface/deltaR.h"

//#include "L1Trigger/L1TCaloLayer1/src/UCTGeometry.hh"
//#include "L1Trigger/Stage3Ntuplizer/plugins/UCTRegionProcess.hh"
#include "L1Trigger/phase2L1EcalTimingTP/interface/P2L1TEventDisplayGenerator.h"
#include "L1Trigger/phase2L1EcalTimingTP/interface/triggerGeometryTools.hh"

using namespace edm;
using std::cout;
using std::endl;
using std::vector;

bool compareByPtLorentz (TLorentzVector i,TLorentzVector j) { 
  return(i.Pt() > j.Pt()); 
	 };


P2L1TEventDisplayGenerator::P2L1TEventDisplayGenerator( const ParameterSet & cfg ) :
  debug(             cfg.getUntrackedParameter<bool>("debug", false)),
  hcalSrc_(          consumes< HcalTrigPrimDigiCollection>(   cfg.getParameter<edm::InputTag>("hcalDigis"))),
  vtxLabel_(         consumes< reco::VertexCollection>(       cfg.getParameter<edm::InputTag>("vertices"))),
  ecalTPGBToken_(    consumes< EcalEBTrigPrimDigiCollection>( cfg.getParameter<edm::InputTag>("ecalTPGsBarrel"))),
  L1PFTauToken_(     consumes< L1PFTauCollection    >(        cfg.getParameter<edm::InputTag>("l1TauObjects"))),
  genSrc_ (          cfg.getParameter<edm::InputTag>( "genParticles")),
  genMatchDeltaRcut( cfg.getUntrackedParameter<double>("genMatchDeltaRcut", 0.1))
  //  tauSrc_(consumes<vector<pat::Tau> >(cfg.getParameter<edm::InputTag>("recoTau"))),
  //rctSource_(cfg.getParameter<edm::InputTag>("rctSource")),
  //packedPfCandsToken_(consumes<vector<pat::PackedCandidate> >(cfg.getParameter<edm::InputTag>("packedPfCands"))),
  //pfCandsToken_(consumes<vector<reco::PFCandidate> >(cfg.getParameter<edm::InputTag>("pfCands"))),
  //discriminatorMu_(consumes<reco::PFTauDiscriminator>(cfg.getParameter<edm::InputTag>("recoTauDiscriminatorMu"))),
  //discriminatorIso_(consumes<reco::PFTauDiscriminator>(cfg.getParameter<edm::InputTag>("recoTauDiscriminatorIso"))),
  //tauSrc_(consumes<vector<reco::PFTau> >(cfg.getParameter<edm::InputTag>("recoTaus"))),
  //slimmedTauSrc_(consumes<vector<pat::Tau> >(cfg.getParameter<edm::InputTag>("slimmedTaus"))),
  //jetSrc_(consumes<vector<pat::Jet> >(cfg.getParameter<edm::InputTag>("recoJets"))),
  //l1ExtraIsoTauSource_(consumes<vector <l1extra::L1JetParticle> >(cfg.getParameter<edm::InputTag>("l1ExtraIsoTauSource"))),
  //l1ExtraTauSource_(consumes<vector <l1extra::L1JetParticle> >(cfg.getParameter<edm::InputTag>("l1ExtraTauSource"))),
  //l1ExtraJetSource_(consumes<vector <l1extra::L1JetParticle> >(cfg.getParameter<edm::InputTag>("l1ExtraJetSource"))),
  //regionSource_(consumes<vector <L1CaloRegion> >(cfg.getParameter<edm::InputTag>("UCTRegion"))),
  //ecalCaloSrc_(consumes<vector <reco::CaloCluster> >(cfg.getParameter<edm::InputTag>("ecalCaloClusters")))
  {
    usesResource("TFileService");

    gROOT->ProcessLine(".L /cms/ojalvo/triggerPhaseII/CMSSW_9_1_0_pre2/src/L1Trigger/phase2L1EcalTimingAnalyzer/test/loader.C+");

    L1TrackInputTag = cfg.getParameter<edm::InputTag>("L1TrackInputTag");
    ttTrackToken_ = consumes< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > >(L1TrackInputTag);  
    genToken_ =     consumes<std::vector<reco::GenParticle> >(genSrc_);
    folderName_          = cfg.getUntrackedParameter<std::string>("folderName");
    recoPt_              = cfg.getParameter<double>("recoPtCut");

    edm::Service<TFileService> tfs_;
    efficiencyTree = tfs_->make<TTree>("EfficiencyTree", "Efficiency Tree");

    efficiencyTree->Branch("hcalTpgs_Pt",  &hcalTpgs_Pt); 
    efficiencyTree->Branch("hcalTpgs_Eta", &hcalTpgs_Eta); 
    efficiencyTree->Branch("hcalTpgs_Phi", &hcalTpgs_Phi); 

    efficiencyTree->Branch("ecalTpgs_Pt",  &ecalTpgs_Pt); 
    efficiencyTree->Branch("ecalTpgs_Eta", &ecalTpgs_Eta); 
    efficiencyTree->Branch("ecalTpgs_Phi", &ecalTpgs_Phi); 

    efficiencyTree->Branch("ecalCrys_Pt",  &ecalCrys_Pt); 
    efficiencyTree->Branch("ecalCrys_Eta", &ecalCrys_Eta); 
    efficiencyTree->Branch("ecalCrys_Phi", &ecalCrys_Phi); 


    efficiencyTree->Branch("sumTpgs_Pt",  &sumTpgs_Pt); 
    efficiencyTree->Branch("sumTpgs_Eta", &sumTpgs_Eta); 
    efficiencyTree->Branch("sumTpgs_Phi", &sumTpgs_Phi); 

    //putting bufsize at 32000 and changing split level to 0 so that the branch isn't split into multiple branches

    //efficiencyTree->Branch("rlxTaus", "vector<TLorentzVector>", &rlxTaus, 32000, 0); 
    //efficiencyTree->Branch("isoTaus", "vector<TLorentzVector>", &isoTaus, 32000, 0); 
    //efficiencyTree->Branch("recoTaus", "vector<TLorentzVector>", &recoTaus, 32000, 0); 
    //efficiencyTree->Branch("allRegions", "vector<TLorentzVector>", &allRegions, 32000, 0); 
    //efficiencyTree->Branch("hcalTPGs", "ROOT.std.vector(ROOT.TLorentzVector())()", &allHcalTPGs, 32000, 0); 
    efficiencyTree->Branch("hcalTPGs", "std::vector<TLorentzVector>", &allHcalTPGs, 32000, 0); 
    efficiencyTree->Branch("ecalTPGs", "std::vector<TLorentzVector>", &allEcalTPGs, 32000, 0); 
    //efficiencyTree->Branch("signalPFCands", "vector<TLorentzVector>", &signalPFCands, 32000, 0); 
    //efficiencyTree->Branch("l1Jets", "vector<TLorentzVector>", &l1Jets, 32000, 0); 
    //efficiencyTree->Branch("recoJets", "vector<TLorentzVector>", &recoJets, 32000, 0); 
    //efficiencyTree->Branch("recoJetsDR", "vector<double>", &recoJetsDr, 32000, 0); 
    //efficiencyTree->Branch("caloClusters", "vector<TLorentzVector>", &caloClusters, 32000, 0); 

    efficiencyTree->Branch("l1Tracks", "vector<TLorentzVector>", &l1Tracks, 32000, 0); 
    efficiencyTree->Branch("l1EcalClusters", "vector<TLorentzVector>", &l1EcalClusters, 32000, 0); 
    efficiencyTree->Branch("l1EcalCrystals", "vector<TLorentzVector>", &l1EcalCrystals, 32000, 0); 
    efficiencyTree->Branch("l1PFTaus",       "vector<TLorentzVector>", &l1PFTaus, 32000, 0); 
    efficiencyTree->Branch("genTaus", "vector<TLorentzVector>", &genHadronicTaus, 32000, 0); 

    efficiencyTree->Branch("run",    &run,     "run/I");
    efficiencyTree->Branch("lumi",   &lumi,    "lumi/I");
    efficiencyTree->Branch("event",  &event,   "event/I");
    efficiencyTree->Branch("nvtx",   &nvtx,         "nvtx/I");

    efficiencyTree->Branch("decayMode", &decayMode,   "decayMode/I");
    
    efficiencyTree->Branch("tauEtaEcalEnt", &tauEtaEcalEnt,"tauEtaEcalEnt/D");
    efficiencyTree->Branch("tauPhiEcalEnt", &tauPhiEcalEnt,"tauPhiEcalEnt/D");

    efficiencyTree->Branch("recoPt",        &recoPt,   "recoPt/D");
    efficiencyTree->Branch("recoEta",       &recoEta,   "recoEta/D");
    efficiencyTree->Branch("recoPhi",       &recoPhi,   "recoPhi/D");

    efficiencyTree->Branch("genPt",        &genPt,   "genPt/D");
    efficiencyTree->Branch("genEta",       &genEta,  "genEta/D");
    efficiencyTree->Branch("genPhi",       &genPhi,  "genPhi/D");

    efficiencyTree->Branch("genPizeroPt",  &genPizeroPt,   "genPizeroPt/D");
    efficiencyTree->Branch("genPizeroEta", &genPizeroEta,  "genPizeroEta/D");
    efficiencyTree->Branch("genPizeroPhi", &genPizeroPhi,  "genPizeroPhi/D");

    /////////////////////////PiZero Efficiency Tree
    efficiencyTreePiZero = tfs_->make<TTree>("EfficiencyTreePiZero", "Efficiency Tree PiZero");

    efficiencyTreePiZero->Branch("hcalTpgs_Pt",  &hcalTpgs_Pt); 
    efficiencyTreePiZero->Branch("hcalTpgs_Eta", &hcalTpgs_Eta); 
    efficiencyTreePiZero->Branch("hcalTpgs_Phi", &hcalTpgs_Phi); 

    efficiencyTreePiZero->Branch("ecalTpgs_Pt",  &ecalTpgs_Pt); 
    efficiencyTreePiZero->Branch("ecalTpgs_Eta", &ecalTpgs_Eta); 
    efficiencyTreePiZero->Branch("ecalTpgs_Phi", &ecalTpgs_Phi); 

    efficiencyTreePiZero->Branch("ecalCrys_Pt",  &ecalCrys_Pt); 
    efficiencyTreePiZero->Branch("ecalCrys_Eta", &ecalCrys_Eta); 
    efficiencyTreePiZero->Branch("ecalCrys_Phi", &ecalCrys_Phi); 


    efficiencyTreePiZero->Branch("sumTpgs_Pt",  &sumTpgs_Pt); 
    efficiencyTreePiZero->Branch("sumTpgs_Eta", &sumTpgs_Eta); 
    efficiencyTreePiZero->Branch("sumTpgs_Phi", &sumTpgs_Phi); 

    //putting bufsize at 32000 and changing split level to 0 so that the branch isn't split into multiple branches

     efficiencyTreePiZero->Branch("hcalTPGs", "std::vector<TLorentzVector>", &allHcalTPGs, 32000, 0); 
     efficiencyTreePiZero->Branch("ecalTPGs", "std::vector<TLorentzVector>", &allEcalTPGs, 32000, 0); 

    efficiencyTreePiZero->Branch("l1Tracks", "vector<TLorentzVector>", &l1Tracks, 32000, 0); 
    efficiencyTreePiZero->Branch("l1EcalClusters", "vector<TLorentzVector>", &l1EcalClusters, 32000, 0); 
    efficiencyTreePiZero->Branch("l1EcalCrystals", "vector<TLorentzVector>", &l1EcalCrystals, 32000, 0); 
    efficiencyTreePiZero->Branch("genTaus", "vector<TLorentzVector>", &genHadronicTaus, 32000, 0); 

    efficiencyTreePiZero->Branch("run",    &run,     "run/I");
    efficiencyTreePiZero->Branch("lumi",   &lumi,    "lumi/I");
    efficiencyTreePiZero->Branch("event",  &event,   "event/I");
    efficiencyTreePiZero->Branch("nvtx",   &nvtx,         "nvtx/I");

    efficiencyTreePiZero->Branch("decayMode", &decayMode,   "decayMode/I");
    
    efficiencyTreePiZero->Branch("recoPt",        &recoPt,   "recoPt/D");
    efficiencyTreePiZero->Branch("recoEta",       &recoEta,   "recoEta/D");
    efficiencyTreePiZero->Branch("recoPhi",       &recoPhi,   "recoPhi/D");

    efficiencyTreePiZero->Branch("genPt",        &genPt,   "genPt/D");
    efficiencyTreePiZero->Branch("genEta",       &genEta,  "genEta/D");
    efficiencyTreePiZero->Branch("genPhi",       &genPhi,  "genPhi/D");


    /////////////////////////PiPM Efficiency Tree
    efficiencyTreePiPM = tfs_->make<TTree>("EfficiencyTreePiPM", "Efficiency Tree PiPM");

    efficiencyTreePiPM->Branch("hcalTpgs_Pt",  &hcalTpgs_Pt); 
    efficiencyTreePiPM->Branch("hcalTpgs_Eta", &hcalTpgs_Eta); 
    efficiencyTreePiPM->Branch("hcalTpgs_Phi", &hcalTpgs_Phi); 

    efficiencyTreePiPM->Branch("ecalTpgs_Pt",  &ecalTpgs_Pt); 
    efficiencyTreePiPM->Branch("ecalTpgs_Eta", &ecalTpgs_Eta); 
    efficiencyTreePiPM->Branch("ecalTpgs_Phi", &ecalTpgs_Phi); 

    efficiencyTreePiPM->Branch("ecalCrys_Pt",  &ecalCrys_Pt); 
    efficiencyTreePiPM->Branch("ecalCrys_Eta", &ecalCrys_Eta); 
    efficiencyTreePiPM->Branch("ecalCrys_Phi", &ecalCrys_Phi); 


    efficiencyTreePiPM->Branch("sumTpgs_Pt",  &sumTpgs_Pt); 
    efficiencyTreePiPM->Branch("sumTpgs_Eta", &sumTpgs_Eta); 
    efficiencyTreePiPM->Branch("sumTpgs_Phi", &sumTpgs_Phi); 

    //putting bufsize at 32000 and changing split level to 0 so that the branch isn't split into multiple branches

     efficiencyTreePiPM->Branch("hcalTPGs", "std::vector<TLorentzVector>", &allHcalTPGs, 32000, 0); 
     efficiencyTreePiPM->Branch("ecalTPGs", "std::vector<TLorentzVector>", &allEcalTPGs, 32000, 0); 

    efficiencyTreePiPM->Branch("l1Tracks", "vector<TLorentzVector>", &l1Tracks, 32000, 0); 
    efficiencyTreePiPM->Branch("l1EcalClusters", "vector<TLorentzVector>", &l1EcalClusters, 32000, 0); 
    efficiencyTreePiPM->Branch("l1EcalCrystals", "vector<TLorentzVector>", &l1EcalCrystals, 32000, 0); 
    efficiencyTreePiPM->Branch("genTaus", "vector<TLorentzVector>", &genHadronicTaus, 32000, 0); 

    efficiencyTreePiPM->Branch("run",    &run,     "run/I");
    efficiencyTreePiPM->Branch("lumi",   &lumi,    "lumi/I");
    efficiencyTreePiPM->Branch("event",  &event,   "event/I");
    efficiencyTreePiPM->Branch("nvtx",   &nvtx,         "nvtx/I");

    efficiencyTreePiPM->Branch("decayMode", &decayMode,   "decayMode/I");
    
    efficiencyTreePiPM->Branch("recoPt",        &recoPt,   "recoPt/D");
    efficiencyTreePiPM->Branch("recoEta",       &recoEta,   "recoEta/D");
    efficiencyTreePiPM->Branch("recoPhi",       &recoPhi,   "recoPhi/D");

    efficiencyTreePiPM->Branch("genPt",        &genPt,   "genPt/D");
    efficiencyTreePiPM->Branch("genEta",       &genEta,  "genEta/D");
    efficiencyTreePiPM->Branch("genPhi",       &genPhi,  "genPhi/D");

  }

void P2L1TEventDisplayGenerator::beginJob( const EventSetup & es) {
  gROOT->ProcessLine("#include <vector>");
  gROOT->ProcessLine("#include <TLorentzVector.h>");

}

//unsigned int const P2L1TEventDisplayGenerator::N_TOWER_PHI = 72;
//unsigned int const P2L1TEventDisplayGenerator::N_TOWER_ETA = 56;

void P2L1TEventDisplayGenerator::analyze( const Event& evt, const EventSetup& es )
{
  gROOT->ProcessLine("#include <vector>");
  gROOT->ProcessLine("#include <TLorentzVector.h>");

  // Set up the ECAL Crystals
  edm::ESHandle<CaloGeometry> caloGeometryHandle;
  es.get<CaloGeometryRecord>().get(caloGeometryHandle);
  const CaloGeometry* caloGeometry_ = caloGeometryHandle.product();  
  ebGeometry = caloGeometry_->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
  
  run = evt.id().run();
  lumi = evt.id().luminosityBlock();
  event = evt.id().event();
   // L1 Tracks
  edm::Handle<L1TkTrackCollectionType> l1trackHandle;
  evt.getByLabel(L1TrackInputTag, l1trackHandle);

  edm::Handle<EcalEBTrigPrimDigiCollection> ecaltpgCollection;
  evt.getByToken( ecalTPGBToken_, ecaltpgCollection);
  vector<TLorentzVector> ecalCrystals;

  edm::Handle< std::vector<L1PFTau> > L1PFTaus;
  evt.getByToken( L1PFTauToken_, L1PFTaus);

  // Get the ECAL Crystals as ecalCrystal_t
  getEcalCrystals(ecaltpgCollection, ecalCrystals);
  
  // Get genParticles
  edm::Handle<GenParticleCollectionType> genParticleHandle;
  if(!evt.getByToken(genToken_,genParticleHandle))
    std::cout<<"No gen Particles Found "<<std::endl;
  
  vector<TTTrack<Ref_Phase2TrackerDigi_> > l1TracksRef;
  
  edm::Handle<EcalTrigPrimDigiCollection> ecalTPGs;
  edm::Handle<HcalTrigPrimDigiCollection> hcalTPGs;  
  
   //Clear the vectors
  l1Tracks->clear(); 
  l1EcalClusters->clear(); 
  l1EcalCrystals->clear(); 
  genHadronicTaus->clear(); 
  allHcalTPGs->clear();

  vector<reco::GenParticle> genTaus;
  vector<reco::GenParticle> genPiZeros;
  vector<reco::GenParticle> genPiPMs;
  vector<reco::GenParticle> genParticles;

  for(unsigned int i = 0; i< genParticleHandle->size(); i++){
    edm::Ptr<reco::GenParticle> ptr(genParticleHandle, i);
    genParticles.push_back(*ptr);
    if(abs(ptr->pdgId())==15){
      genTaus.push_back(*ptr);
    }
    if(abs(ptr->pdgId())==111 && abs(ptr->eta())<1.74){
       genPiZeros.push_back(*ptr);
       std::cout<<"Found PiZero PDGID 111 pt: "<<ptr->pt()<<" eta: "<<ptr->eta()<<" phi: "<<ptr->phi()<<std::endl;
     }

     if(abs(ptr->pdgId())==211 && abs(ptr->eta())<1.74){
       genPiPMs.push_back(*ptr);
       std::cout<<"Found PiZero PDGID 211 pt: "<<ptr->pt()<<" eta: "<<ptr->eta()<<" phi: "<<ptr->phi()<<std::endl;
     }

  }


  for(auto genPiZero: genPiZeros){
    clearThePlottingVectors();
    
    std::cout<<"got PiZero"<<std::endl;
    TLorentzVector piZero;
    //std::cout<<"Tau Decay Mode "<<GetDecayMode(&genTau)<<std::endl;
    //decayMode = GetDecayModePiZero(&genTau,piZero);
    //onlygetting the hadronic taus
    //if(decayMode<10)
    //continue;
    
    //reco::Candidate::LorentzVector visGenTau= getVisMomentum(&genPiZero, &genParticles);
    //std::cout<<"vis pi0 et: "<<visGenTau.pt()<<" eta: "<<visGenTau.eta()<<" phi: "visGenTau.phi()<<std::endl;
    TLorentzVector genPiZeroTemp;
    float et, eta, phi;
    et     = genPiZero.pt();
    genPt  = genPiZero.pt();
    
    eta    = genPiZero.eta();
    genEta = genPiZero.eta();
    
    phi    = genPiZero.phi();
    genPhi = genPiZero.phi();
    
    genPiZeroTemp.SetPtEtaPhiE(et,eta,phi,et);
    genHadronicTaus->push_back(genPiZeroTemp);
    
    
    //Find and sort the tracks
    for(size_t track_index=0; track_index<l1trackHandle->size(); ++track_index)
      {
	edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_>> ptr(l1trackHandle, track_index);
	//std::cout<<"# "<<track_index<<": "<<pt<<std::endl;
	TLorentzVector l1TrackTemp;
       float et, eta, phi;
       et  = ptr->getMomentum().perp();
       eta = ptr->getMomentum().eta();
       phi = ptr->getMomentum().phi();
       l1TrackTemp.SetPtEtaPhiE(et,eta,phi,et);
       l1Tracks->push_back(l1TrackTemp);
      }
    //std::sort(l1Tracks->begin(),l1Tracks->end(),compareByPtLorentz);   
    
    //Find and sort the crystals
    for(size_t crystal_index=0; crystal_index< ecalCrystals.size(); ++crystal_index)
      {
	//edm::Ptr<l1extra::L1EmParticle> ptr(l1EGCrystalHandle, crystal_index);
	//std::cout<<"# "<<crystal_index<<": "<<pt<<std::endl;
	TLorentzVector ecalCrystal = ecalCrystals.at(crystal_index);
	TLorentzVector l1CrystalTemp;
	float et, eta, phi;
	et  = ecalCrystal.Pt();
	eta = ecalCrystal.Eta();
	phi = ecalCrystal.Phi();//+0.00667498;
	l1CrystalTemp.SetPtEtaPhiE(et,eta,phi,et);
	l1EcalCrystals->push_back(l1CrystalTemp);
	//std::cout<<"crystal pt:"<<et<<" eta: "<<eta<<" phi: "<<phi<<std::endl;
      }
    std::sort(l1EcalCrystals->begin(),l1EcalCrystals->end(),compareByPtLorentz);   

  ESHandle<L1CaloHcalScale> hcalScale;
  es.get<L1CaloHcalScaleRcd>().get(hcalScale);

  if(!evt.getByToken(hcalSrc_, hcalTPGs))
    std::cout<<"ERROR GETTING THE HCAL TPGS"<<std::endl;
  else
    for (size_t i = 0; i < hcalTPGs->size(); ++i) {
      HcalTriggerPrimitiveDigi tpg = (*hcalTPGs)[i];
      int cal_ieta = tpg.id().ieta();
      int cal_iphi = tpg.id().iphi();
      if(cal_ieta>28)continue; 
      if(cal_ieta<-28)continue; 
      int ieta = TPGEtaRange(cal_ieta);
      short absieta = std::abs(tpg.id().ieta());
      short zside = tpg.id().zside();
      double et = hcalScale->et(tpg.SOI_compressedEt(), absieta, zside); 
      //if(et>0)
      //std::cout<<"HCAL ET "<<et<<std::endl;
      if(ieta<0){
	std::cout<<"sorry, ieta less than 1 :("<<std::endl;
	std::cout<<"cal_ieta "<<cal_ieta<<" ieta "<<ieta<<std::endl;
      }
      triggerGeometryTools trigTools;
      float eta = trigTools.getRecoEta(ieta, zside);
      float phi = trigTools.getRecoPhi(cal_iphi);    
      TLorentzVector temp ;
      temp.SetPtEtaPhiE(et,eta,phi,et);
      allHcalTPGs->push_back(temp);
    }
  
  efficiencyTreePiZero->Fill();
  } 


  ////////////////pipm
  for(auto genPiPM: genPiPMs){
    clearThePlottingVectors();
    std::cout<<"got PiPM"<<std::endl;
    TLorentzVector piZero;

    TLorentzVector genPiPMTemp;
    float et, eta, phi;
    et     = genPiPM.pt();
    genPt  = genPiPM.pt();
    
    eta    = genPiPM.eta();
    genEta = genPiPM.eta();
    
    phi    = genPiPM.phi();
    genPhi = genPiPM.phi();
    
    genPiPMTemp.SetPtEtaPhiE(et,eta,phi,et);
    genHadronicTaus->push_back(genPiPMTemp);
    
    
    //Find and sort the tracks
    for(size_t track_index=0; track_index<l1trackHandle->size(); ++track_index)
      {
	edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_>> ptr(l1trackHandle, track_index);
	//std::cout<<"# "<<track_index<<": "<<pt<<std::endl;
	TLorentzVector l1TrackTemp;
       float et, eta, phi;
       et  = ptr->getMomentum().perp();
       eta = ptr->getMomentum().eta();
       phi = ptr->getMomentum().phi();
       l1TrackTemp.SetPtEtaPhiE(et,eta,phi,et);
       l1Tracks->push_back(l1TrackTemp);
      }
    //std::sort(l1Tracks->begin(),l1Tracks->end(),compareByPtLorentz);   
    
    //Find and sort the crystals
    
    for(size_t crystal_index=0; crystal_index< ecalCrystals.size(); ++crystal_index)
      {
	//edm::Ptr<l1extra::L1EmParticle> ptr(l1EGCrystalHandle, crystal_index);
	//std::cout<<"# "<<crystal_index<<": "<<pt<<std::endl;
	TLorentzVector ecalCrystal = ecalCrystals.at(crystal_index);
	TLorentzVector l1CrystalTemp;
	float et, eta, phi;
	et  = ecalCrystal.Pt();
	eta = ecalCrystal.Eta();
	phi = ecalCrystal.Phi();
	l1CrystalTemp.SetPtEtaPhiE(et,eta,phi,et);
	l1EcalCrystals->push_back(l1CrystalTemp);
	//std::cout<<"crystal pt:"<<et<<" eta: "<<eta<<" phi: "<<phi<<std::endl;
      }
    std::sort(l1EcalCrystals->begin(),l1EcalCrystals->end(),compareByPtLorentz);   
    
    //// Put in dummy crystals
    /*
    for(float jEta = -2.5; jEta < 2.5 ; jEta += 0.0174){
      for(float jPhi = -2.5; jPhi < 2.5 ; jPhi += 0.01745){
	TLorentzVector l1CrystalTemp;
	float et, eta, phi;
	if(jEta > 0.4 && jEta < 0.52 && jPhi > -1.9 && jPhi <  -1.7 )
	  et = 3;
	else
	  et = 0;
	eta = jEta;
	phi = jPhi;
	l1CrystalTemp.SetPtEtaPhiE(et,eta,phi,et);
	l1EcalCrystals->push_back(l1CrystalTemp);
      }
    }
*/


  ESHandle<L1CaloHcalScale> hcalScale;
  es.get<L1CaloHcalScaleRcd>().get(hcalScale);

  if(!evt.getByToken(hcalSrc_, hcalTPGs))
    std::cout<<"ERROR GETTING THE HCAL TPGS"<<std::endl;
  else
    for (size_t i = 0; i < hcalTPGs->size(); ++i) {
      HcalTriggerPrimitiveDigi tpg = (*hcalTPGs)[i];
      int cal_ieta = tpg.id().ieta();
      int cal_iphi = tpg.id().iphi();
      if(cal_ieta>28)continue; 
      if(cal_ieta<-28)continue; 
      int ieta = TPGEtaRange(cal_ieta);
      short absieta = std::abs(tpg.id().ieta());
      short zside = tpg.id().zside();
      double et = hcalScale->et(tpg.SOI_compressedEt(), absieta, zside); 
      //if(et>0)
      //std::cout<<"HCAL ET "<<et<<std::endl;
      if(ieta<0){
	std::cout<<"sorry, ieta less than 1 :("<<std::endl;
	std::cout<<"cal_ieta "<<cal_ieta<<" ieta "<<ieta<<std::endl;
      }
      triggerGeometryTools trigTools;
      float eta = trigTools.getRecoEta(ieta, zside);
      float phi = trigTools.getRecoPhi(cal_iphi);    
      TLorentzVector temp ;
      temp.SetPtEtaPhiE(et,eta,phi,et);
      allHcalTPGs->push_back(temp);
    }
  
  efficiencyTreePiPM->Fill();
  } 

  ///////////
  
  for(auto genTau: genTaus){
    clearThePlottingVectors();
    std::cout<<"got tau"<<std::endl;
    TLorentzVector piZero;
    std::cout<<"Tau Decay Mode "<<GetDecayMode(&genTau)<<std::endl;
    decayMode = GetDecayModePiZero(&genTau,piZero);
    //onlygetting the hadronic taus
    if(decayMode<10)
      continue;
    
    reco::Candidate::LorentzVector visGenTau= getVisMomentum(&genTau, &genParticles);
    
    TLorentzVector genTauTemp;
    float et, eta, phi;
    et     = visGenTau.pt();
    genPt  = visGenTau.pt();
    
    eta    = visGenTau.eta();
    genEta = visGenTau.eta();
    
    phi    = visGenTau.phi();
    genPhi = visGenTau.phi();
    
    genTauTemp.SetPtEtaPhiE(et,eta,phi,et);
    genHadronicTaus->push_back(genTauTemp);
    
    
    //Find and sort the tracks
    for(size_t track_index=0; track_index<l1trackHandle->size(); ++track_index)
      {
	edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_>> ptr(l1trackHandle, track_index);
	//std::cout<<"# "<<track_index<<": "<<pt<<std::endl;
	TLorentzVector l1TrackTemp;
       float et, eta, phi;
       et  = ptr->getMomentum().perp();
       eta = ptr->getMomentum().eta();
       phi = ptr->getMomentum().phi();
       l1TrackTemp.SetPtEtaPhiE(et,eta,phi,et);
       l1Tracks->push_back(l1TrackTemp);
      }
    //std::sort(l1Tracks->begin(),l1Tracks->end(),compareByPtLorentz);   
    
    //Find and sort the tracks
    for(size_t crystal_index=0; crystal_index< ecalCrystals.size(); ++crystal_index)
      {
	//edm::Ptr<l1extra::L1EmParticle> ptr(l1EGCrystalHandle, crystal_index);
	//std::cout<<"# "<<crystal_index<<": "<<pt<<std::endl;
	TLorentzVector ecalCrystal = ecalCrystals.at(crystal_index);
	TLorentzVector l1CrystalTemp;
	float et, eta, phi;
	et  = ecalCrystal.Pt();
	eta = ecalCrystal.Eta();
	phi = ecalCrystal.Phi();
	l1CrystalTemp.SetPtEtaPhiE(et,eta,phi,et);
	l1EcalCrystals->push_back(l1CrystalTemp);
	//std::cout<<"crystal pt:"<<et<<" eta: "<<eta<<" phi: "<<phi<<std::endl;
      }
    std::sort(l1EcalCrystals->begin(),l1EcalCrystals->end(),compareByPtLorentz);   

  ESHandle<L1CaloHcalScale> hcalScale;
  es.get<L1CaloHcalScaleRcd>().get(hcalScale);

  if(!evt.getByToken(hcalSrc_, hcalTPGs))
    std::cout<<"ERROR GETTING THE HCAL TPGS"<<std::endl;
  else{
    for (size_t i = 0; i < hcalTPGs->size(); ++i) {
      HcalTriggerPrimitiveDigi tpg = (*hcalTPGs)[i];
      int cal_ieta = tpg.id().ieta();
      int cal_iphi = tpg.id().iphi();
      if(cal_ieta>28)continue; 
      if(cal_ieta<-28)continue; 
      int ieta = TPGEtaRange(cal_ieta);
      short absieta = std::abs(tpg.id().ieta());
      short zside = tpg.id().zside();
      double et = hcalScale->et(tpg.SOI_compressedEt(), absieta, zside); 
      //if(et>0)
      //std::cout<<"HCAL ET "<<et<<std::endl;
      if(ieta<0){
	std::cout<<"sorry, ieta less than 1 :("<<std::endl;
	std::cout<<"cal_ieta "<<cal_ieta<<" ieta "<<ieta<<std::endl;
      }
      triggerGeometryTools trigTools;
      float eta = trigTools.getRecoEta(ieta, zside);
      float phi = trigTools.getRecoPhi(cal_iphi);    
      TLorentzVector temp ;
      temp.SetPtEtaPhiE(et,eta,phi,et);
      allHcalTPGs->push_back(temp);
    }
  }
    for(size_t l1tau_index=0; l1tau_index< L1PFTaus->size(); ++l1tau_index)
      {
	//edm::Ptr<l1extra::L1EmParticle> ptr(l1EGCrystalHandle, crystal_index);
	//std::cout<<"# "<<crystal_index<<": "<<pt<<std::endl;
	L1PFTau l1Tau = L1PFTaus->at(l1tau_index);
	TLorentzVector l1TauTemp;
	float et, eta, phi;
	et  = l1Tau.p4().Pt();
	eta = l1Tau.p4().Eta();
	phi = l1Tau.p4().Phi();
	l1TauTemp.SetPtEtaPhiE(et,eta,phi,et);
	l1PFTaus->push_back(l1TauTemp);
	//std::cout<<"crystal pt:"<<et<<" eta: "<<eta<<" phi: "<<phi<<std::endl;
      }
    std::sort(l1PFTaus->begin(),l1PFTaus->end(),compareByPtLorentz);   

  efficiencyTree->Fill();
  }
}


void P2L1TEventDisplayGenerator::getEcalCrystals(edm::Handle<EcalEBTrigPrimDigiCollection> ecalTPGs, vector<TLorentzVector> &ecalCrystals)
{

  if(debug)
    std::cout<<"-------------- ECAL Crystals --------------"<<std::endl;
    
  for(auto& tpg : *ecalTPGs.product())
    {

      if(tpg.encodedEt() > 0) 
	{

	  GlobalVector position;
	  auto cell = ebGeometry->getGeometry(tpg.id());

	  float et = tpg.encodedEt()/8.;

	  if(et<0.1) continue;// LSB is 0.1
	  //float energy = et / sin(position.theta());
	  float eta = cell->getPosition().eta();
	  float phi = cell->getPosition().phi();
	  EBDetId detID = tpg.id();
	  float iEta = detID.ieta();
	  float iPhi = detID.iphi();

	  //DEBUG STATMENT
	  if(debug)
	    if(et>1){
	      std::cout<<"ET "<<et<<std::endl;
	      std::cout<<" eta  "<< eta<< " phi  "<< phi<<std::endl;
	      std::cout<<" iEta"<<iEta<<" iphi "<<iPhi<<std::endl;
	    }
	  TLorentzVector tempCrystal;
	  tempCrystal.SetPtEtaPhiE(et, eta, phi,et);
	  //ecalCrystal_t tempCrystal;
	  //tempCrystal.p4.SetPtEtaPhiE(et, eta, phi,et);
	  //tempCrystal.iEta = iEta;
	  //tempCrystal.iPhi = iPhi;
	  //tempCrystal.id = detID;

	  ecalCrystals.push_back(tempCrystal);
	}
    }
}


void P2L1TEventDisplayGenerator::endJob() {
}

P2L1TEventDisplayGenerator::~P2L1TEventDisplayGenerator(){
}

DEFINE_FWK_MODULE(P2L1TEventDisplayGenerator);
