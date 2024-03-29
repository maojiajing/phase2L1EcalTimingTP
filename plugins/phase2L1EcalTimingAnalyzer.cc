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
#include "phase2L1EcalTimingAnalyzer.h"
#include <random>
#define debug 0
//#include "L1Trigger/phase2L1EcalTimingTP/plugins/phase2L1EcalTimingAnalyzer.h"
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
phase2L1EcalTimingAnalyzer::phase2L1EcalTimingAnalyzer(const edm::ParameterSet& cfg):
  ecalTPGBToken_(   consumes<EcalEBTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("ecalTPGsBarrel"))),
  caloJetTag_ ((        cfg.getParameter<edm::InputTag>( "caloJets"))),
  genSrcak4J_ ((        cfg.getParameter<edm::InputTag>( "ak4GenJets"))),
  genSrcak4JN_ ((        cfg.getParameter<edm::InputTag>( "ak4GenJetsNoNu"))),
  genSrcak8J_ ((        cfg.getParameter<edm::InputTag>( "ak8GenJets"))),
  genSrcak8JN_ ((        cfg.getParameter<edm::InputTag>( "ak8GenJetsNoNu"))),
  genSrc_ ((        cfg.getParameter<edm::InputTag>( "genParticles"))),
  genSrcT_ ((        cfg.getParameter<edm::InputTag>( "genParticles_t0")))
{
  //now do what ever initialization is needed
  usesResource("TFileService");
  genak4JetToken_ =     consumes<std::vector<reco::GenJet> >(genSrcak4J_);
  genak4JetNoNuToken_ =     consumes<std::vector<reco::GenJet> >(genSrcak4JN_);
  genak8JetToken_ =     consumes<std::vector<reco::GenJet> >(genSrcak8J_);
  genak8JetNoNuToken_ =     consumes<std::vector<reco::GenJet> >(genSrcak8JN_);
  genToken_ =     consumes<std::vector<reco::GenParticle> >(genSrc_);
  genTokenT_ =     consumes<float>(genSrcT_);
  caloJetToken_ =     consumes<l1t::JetBxCollection>(caloJetTag_);

  edm::Service<TFileService> fs;
  ecalTPTree = fs->make<TTree>("ecalTPTree", "Crystal cluster individual crystal pt values");
}


phase2L1EcalTimingAnalyzer::~phase2L1EcalTimingAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)


}

// ------------ load event  ------------
void phase2L1EcalTimingAnalyzer::loadEvent(const edm::Event& iEvent){

  //control plot for ecal crystals
  iEvent.getByToken( ecalTPGBToken_, ecaltpgCollection);

  if(!iEvent.getByToken(caloJetToken_,caloJetHandle))
    std::cout<<"calo jets Found "<<std::endl;
  else{
    //std::cout<<"calo Jets size "<<caloJetHandle->size()<<std::endl;
  }

  if(!iEvent.getByToken(genak4JetToken_,genak4JetHandle))
    std::cout<<"No gen ak4 jets Found "<<std::endl;
  else{
    //std::cout<<"Gen i ak4 Jets size "<<genak4JetHandle->size()<<std::endl;
  }
	
  if(!iEvent.getByToken(genak4JetNoNuToken_,genak4JetNoNuHandle))
    std::cout<<"No gen ak4 nonu jets Found "<<std::endl;
  else{
    //std::cout<<"Gen ak4 nonu Jets size "<<genak4JetNoNuHandle->size()<<std::endl;
  }
	
  if(!iEvent.getByToken(genak8JetToken_,genak8JetHandle))
    std::cout<<"No gen ak8 jets Found "<<std::endl;
  else{
    //std::cout<<"Gen ak8 Jets size "<<genak8JetHandle->size()<<std::endl;
  }
	
  if(!iEvent.getByToken(genak8JetNoNuToken_,genak8JetNoNuHandle))
    std::cout<<"No gen ak8 nonu jets Found "<<std::endl;
  else{
    //std::cout<<"Gen ak8 nonu Jets size "<<genak8JetNoNuHandle->size()<<std::endl;
  }
	
  if(!iEvent.getByToken(genToken_,genParticleHandle))
    std::cout<<"No gen Particles Found "<<std::endl;
  else{
    //std::cout<<"Gen Particles size "<<genParticleHandle->size()<<std::endl;
  }
	
  if(!iEvent.getByToken(genTokenT_,genVertexTHandle))
    std::cout<<"No gen Particles T Found "<<std::endl;
  else{
    //std::cout<<"Gen Particles T value "<<*genVertexTHandle<<std::endl;
    //genVertexT = *genVertexTHandle;
  }

}

// ------------ enable branches  ------------
void phase2L1EcalTimingAnalyzer::setBranches(){
  enableEventInfoBranches();
  //enableEBCrystalBranches();
  enableGenParticleBranches();
  //enableGenak4JetBranches();
  //enableGenak4JetNoNuBranches();
  //enableGenak8JetBranches();
  //enableGenak8JetNoNuBranches();
  enableCaloJetBranches();

};

void phase2L1EcalTimingAnalyzer::enableEventInfoBranches(){
  ecalTPTree->Branch("runNum",    &runNum,     "runNum/I");
  ecalTPTree->Branch("lumiNum",   &lumiNum,    "lumiNum/I");
  ecalTPTree->Branch("eventNum",  &eventNum,   "eventNum/I");
};

void phase2L1EcalTimingAnalyzer::enableEBCrystalBranches(){
 //ecal crystal info
 ecalTPTree->Branch("nCrystals", &nCrystals, "nCrystals/I"); 

 ecalTPTree->Branch("eb_Et", &eb_Et, "eb_Et[nCrystals]/F");
 ecalTPTree->Branch("eb_Edep", &eb_Edep, "eb_Edep[nCrystals]/F");
 ecalTPTree->Branch("eb_time", &eb_time, "eb_time[nCrystals]/F");
 ecalTPTree->Branch("eb_time_60ps", &eb_time_60ps, "eb_time_60ps[nCrystals]/F");
 ecalTPTree->Branch("eb_time_100ps", &eb_time_100ps, "eb_time_100ps[nCrystals]/F");
 ecalTPTree->Branch("eb_time_60ps_func", &eb_time_60ps_func, "eb_time_60ps_func[nCrystals]/F");
 ecalTPTree->Branch("eb_time_100ps_func", &eb_time_100ps_func, "eb_time_100ps_func[nCrystals]/F");
 ecalTPTree->Branch("eb_sigmat", &eb_sigmat, "eb_sigmat[nCrystals]/F");
 ecalTPTree->Branch("eb_sigmat_60ps_func", &eb_sigmat_60ps_func, "eb_sigmat_60ps_func[nCrystals]/F");
 ecalTPTree->Branch("eb_sigmat_100ps_func", &eb_sigmat_100ps_func, "eb_sigmat_100ps_func[nCrystals]/F");
 ecalTPTree->Branch("eb_id", &eb_id, "eb_id[nCrystals]/I");
 ecalTPTree->Branch("eb_ieta", &eb_ieta, "eb_ieta[nCrystals]/I");
 ecalTPTree->Branch("eb_iphi", &eb_iphi, "eb_iphi[nCrystals]/I");
 ecalTPTree->Branch("eb_ism", &eb_ism, "eb_ism[nCrystals]/I");
 ecalTPTree->Branch("eb_ic", &eb_ic, "eb_ic[nCrystals]/I");
 ecalTPTree->Branch("eb_cell_Eta", &eb_cell_Eta, "eb_cell_Eta[nCrystals]/F");
 ecalTPTree->Branch("eb_cell_Phi", &eb_cell_Phi, "eb_cell_Phi[nCrystals]/F");

 ecalTPTree->Branch("gen_time_index", &gen_time_index, "gen_time_index[nCrystals]/I");
 ecalTPTree->Branch("gen_time_e", &gen_time_e, "gen_time_e[nCrystals]/F");
 ecalTPTree->Branch("gen_time_dr", &gen_time_dr, "gen_time_dr[nCrystals]/F");
 ecalTPTree->Branch("gen_time_tp", &gen_time_tp, "gen_time_tp[nCrystals]/F");
 ecalTPTree->Branch("gen_time_tp_60ps", &gen_time_tp_60ps, "gen_time_tp_60ps[nCrystals]/F");
 ecalTPTree->Branch("gen_time_tp_100ps", &gen_time_tp_100ps, "gen_time_tp_100ps[nCrystals]/F");
 ecalTPTree->Branch("gen_time_tp_60ps_func", &gen_time_tp_60ps_func, "gen_time_tp_60ps_func[nCrystals]/F");
 ecalTPTree->Branch("gen_time_tp_100ps_func", &gen_time_tp_100ps_func, "gen_time_tp_100ps_func[nCrystals]/F");
};

void phase2L1EcalTimingAnalyzer::enableGenParticleBranches(){
 //gen info
 ecalTPTree->Branch("nGenParticles", &nGenParticles, "nGenParticles/I"); 

 ecalTPTree->Branch("genVertexX", &genVertexX, "genVertexX/F"); 
 ecalTPTree->Branch("genVertexY", &genVertexY, "genVertexY/F"); 
 ecalTPTree->Branch("genVertexZ", &genVertexZ, "genVertexZ/F"); 
 ecalTPTree->Branch("genVertexT", &genVertexT, "genVertexT/F"); 

 ecalTPTree->Branch("gParticleGrandMotherId", &gParticleGrandMotherId, "gParticleGrandMotherId[nGenParticles]/I");
 ecalTPTree->Branch("gParticleGrandMotherIndex", &gParticleGrandMotherIndex, "gParticleGrandMotherIndex[nGenParticles]/I");

 ecalTPTree->Branch("gParticleMotherId", &gParticleMotherId, "gParticleMotherId[nGenParticles]/I");
 ecalTPTree->Branch("gParticleMotherIndex", &gParticleMotherIndex, "gParticleMotherIndex[nGenParticles]/I");

 ecalTPTree->Branch("gParticleSiblingId", &gParticleSiblingId, "gParticleSiblingId[nGenParticles]/I");
 ecalTPTree->Branch("gParticleSiblingIndex", &gParticleSiblingIndex, "gParticleSiblingIndex[nGenParticles]/I");

 ecalTPTree->Branch("gParticleId", &gParticleId, "gParticleId[nGenParticles]/I");
 ecalTPTree->Branch("gParticleStatus", &gParticleStatus, "gParticleStatus[nGenParticles]/I");

 ecalTPTree->Branch("gParticleE", &gParticleE, "gParticleE[nGenParticles]/F");
 ecalTPTree->Branch("gParticlePt", &gParticlePt, "gParticlePt[nGenParticles]/F");
 ecalTPTree->Branch("gParticlePx", &gParticlePx, "gParticlePx[nGenParticles]/F");
 ecalTPTree->Branch("gParticlePy", &gParticlePy, "gParticlePy[nGenParticles]/F");
 ecalTPTree->Branch("gParticlePz", &gParticlePz, "gParticlePz[nGenParticles]/F");
 ecalTPTree->Branch("gParticleEta", &gParticleEta, "gParticleEta[nGenParticles]/F");
 ecalTPTree->Branch("gParticlePhi", &gParticlePhi, "gParticlePhi[nGenParticles]/F");

 ecalTPTree->Branch("gParticle_decay_vtx_x", &gParticle_decay_vtx_x, "gParticle_decay_vtx_x[nGenParticles]/F");
 ecalTPTree->Branch("gParticle_decay_vtx_y", &gParticle_decay_vtx_y, "gParticle_decay_vtx_y[nGenParticles]/F");
 ecalTPTree->Branch("gParticle_decay_vtx_z", &gParticle_decay_vtx_z, "gParticle_decay_vtx_z[nGenParticles]/F");

 ecalTPTree->Branch("gParticle_prod_vtx_x", &gParticle_prod_vtx_x, "gParticle_prod_vtx_x[nGenParticles]/F");
 ecalTPTree->Branch("gParticle_prod_vtx_y", &gParticle_prod_vtx_y, "gParticle_prod_vtx_y[nGenParticles]/F");
 ecalTPTree->Branch("gParticle_prod_vtx_z", &gParticle_prod_vtx_z, "gParticle_prod_vtx_z[nGenParticles]/F");

 ecalTPTree->Branch("g_eb_time_Emax_02", &g_eb_time_Emax_02, "g_eb_time_Emax_02[nGenParticles]/F");
 ecalTPTree->Branch("g_eb_sigmat_Emax_02", &g_eb_sigmat_Emax_02, "g_eb_sigmat_Emax_02[nGenParticles]/F");
 ecalTPTree->Branch("gEmax_02", &gEmax_02, "gEmax_02[nGenParticles]/F");
 ecalTPTree->Branch("gImax_02", &gImax_02, "gImax_02[nGenParticles]/I");
 ecalTPTree->Branch("gEsc_02", &gEsc_02, "gEsc_02[nGenParticles]/F");
 ecalTPTree->Branch("g_dt_sc_02", &g_dt_sc_02, "g_dt_sc_02[nGenParticles]/F");
 ecalTPTree->Branch("g_cnt_sc_02", &g_cnt_sc_02, "g_cnt_sc_02[nGenParticles]/I");
 ecalTPTree->Branch("g_tmax_sc_02", &g_tmax_sc_02, "g_tmax_sc_02[nGenParticles]/F");
 ecalTPTree->Branch("g_tmin_sc_02", &g_tmin_sc_02, "g_tmin_sc_02[nGenParticles]/F");
 ecalTPTree->Branch("g_dt_sm_02", &g_dt_sm_02, "g_dt_sm_02[nGenParticles]/F");
 ecalTPTree->Branch("g_cnt_sm_02", &g_cnt_sm_02, "g_cnt_sm_02[nGenParticles]/I");
 ecalTPTree->Branch("g_tmax_sm_02", &g_tmax_sm_02, "g_tmax_sm_02[nGenParticles]/F");
 ecalTPTree->Branch("g_tmin_sm_02", &g_tmin_sm_02, "g_tmin_sm_02[nGenParticles]/F");
 ecalTPTree->Branch("genEsc_02", &genEsc_02, "genEsc_02[nGenParticles]/F");
 ecalTPTree->Branch("gE9x9_02", &gE9x9_02, "gE9x9_02[nGenParticles]/F");
 ecalTPTree->Branch("gE5x5_02", &gE5x5_02, "gE5x5_02[nGenParticles]/F");
 ecalTPTree->Branch("gE3x3_02", &gE3x3_02, "gE3x3_02[nGenParticles]/F");

 ecalTPTree->Branch("g_eb_time_Emax_01", &g_eb_time_Emax_01, "g_eb_time_Emax_01[nGenParticles]/F");
 ecalTPTree->Branch("g_eb_sigmat_Emax_01", &g_eb_sigmat_Emax_01, "g_eb_sigmat_Emax_01[nGenParticles]/F");
 ecalTPTree->Branch("gEmax_01", &gEmax_01, "gEmax_01[nGenParticles]/F");
 ecalTPTree->Branch("gImax_01", &gImax_01, "gImax_01[nGenParticles]/I");
 ecalTPTree->Branch("gEsc_01", &gEsc_01, "gEsc_01[nGenParticles]/F");
 ecalTPTree->Branch("g_dt_sc_01", &g_dt_sc_01, "g_dt_sc_01[nGenParticles]/F");
 ecalTPTree->Branch("g_tmax_sc_01", &g_tmax_sc_01, "g_tmax_sc_01[nGenParticles]/F");
 ecalTPTree->Branch("g_cnt_sc_01", &g_cnt_sc_01, "g_cnt_sc_01[nGenParticles]/I");
 ecalTPTree->Branch("g_tmin_sc_01", &g_tmin_sc_01, "g_tmin_sc_01[nGenParticles]/F");
 ecalTPTree->Branch("g_dt_sm_01", &g_dt_sm_01, "g_dt_sm_01[nGenParticles]/F");
 ecalTPTree->Branch("g_tmax_sm_01", &g_tmax_sm_01, "g_tmax_sm_01[nGenParticles]/F");
 ecalTPTree->Branch("g_cnt_sm_01", &g_cnt_sm_01, "g_cnt_sm_01[nGenParticles]/I");
 ecalTPTree->Branch("g_tmin_sm_01", &g_tmin_sm_01, "g_tmin_sm_01[nGenParticles]/F");
 ecalTPTree->Branch("g_dt_sm2_01", &g_dt_sm2_01, "g_dt_sm2_01[nGenParticles]/F");
 ecalTPTree->Branch("g_tmax_sm2_01", &g_tmax_sm2_01, "g_tmax_sm2_01[nGenParticles]/F");
 ecalTPTree->Branch("g_cnt_sm2_01", &g_cnt_sm2_01, "g_cnt_sm2_01[nGenParticles]/I");
 ecalTPTree->Branch("g_tmin_sm2_01", &g_tmin_sm2_01, "g_tmin_sm2_01[nGenParticles]/F");
 ecalTPTree->Branch("g_dt_sm5_01", &g_dt_sm5_01, "g_dt_sm5_01[nGenParticles]/F");
 ecalTPTree->Branch("g_tmax_sm5_01", &g_tmax_sm5_01, "g_tmax_sm5_01[nGenParticles]/F");
 ecalTPTree->Branch("g_cnt_sm5_01", &g_cnt_sm5_01, "g_cnt_sm5_01[nGenParticles]/I");
 ecalTPTree->Branch("g_tmin_sm5_01", &g_tmin_sm5_01, "g_tmin_sm5_01[nGenParticles]/F");
 ecalTPTree->Branch("g_dt_sm10_01", &g_dt_sm10_01, "g_dt_sm10_01[nGenParticles]/F");
 ecalTPTree->Branch("g_tmax_sm10_01", &g_tmax_sm10_01, "g_tmax_sm10_01[nGenParticles]/F");
 ecalTPTree->Branch("g_cnt_sm10_01", &g_cnt_sm10_01, "g_cnt_sm10_01[nGenParticles]/I");
 ecalTPTree->Branch("g_tmin_sm10_01", &g_tmin_sm10_01, "g_tmin_sm10_01[nGenParticles]/F");
 ecalTPTree->Branch("g_dt_max_01", &g_dt_max_01, "g_dt_max_01[nGenParticles]/F");
 ecalTPTree->Branch("g_t_max_01", &g_t_max_01, "g_t_max_01[nGenParticles]/F");
 ecalTPTree->Branch("g_e_max_01", &g_e_max_01, "g_e_max_01[nGenParticles]/F");
 ecalTPTree->Branch("g_dt_min_01", &g_dt_min_01, "g_dt_min_01[nGenParticles]/F");
 ecalTPTree->Branch("g_t_min_01", &g_t_min_01, "g_t_min_01[nGenParticles]/F");
 ecalTPTree->Branch("g_e_min_01", &g_e_min_01, "g_e_min_01[nGenParticles]/F");
 ecalTPTree->Branch("genEsc_01", &genEsc_01, "genEsc_01[nGenParticles]/F");
 ecalTPTree->Branch("gE9x9_01", &gE9x9_01, "gE9x9_01[nGenParticles]/F");
 ecalTPTree->Branch("gE5x5_01", &gE5x5_01, "gE5x5_01[nGenParticles]/F");
 ecalTPTree->Branch("gE3x3_01", &gE3x3_01, "gE3x3_01[nGenParticles]/F");

 ecalTPTree->Branch("g_mindr", &g_mindr, "g_mindr[nGenParticles]/F");
 ecalTPTree->Branch("g_e_mindr", &g_e_mindr, "g_e_mindr[nGenParticles]/F");
 ecalTPTree->Branch("g_t_mindr", &g_t_mindr, "g_t_mindr[nGenParticles]/F");

 ecalTPTree->Branch("gIcore", &gIcore, "gIcore[nGenParticles]/I");
 ecalTPTree->Branch("g_eb_time", &g_eb_time, "g_eb_time[nGenParticles]/F");
 ecalTPTree->Branch("g_eb_sigmat", &g_eb_sigmat, "g_eb_sigmat[nGenParticles]/F");
 ecalTPTree->Branch("g_tof", &g_tof, "g_tof[nGenParticles]/F");
 ecalTPTree->Branch("g_tvirtual", &g_tvirtual, "g_tvirtual[nGenParticles]/F");
 ecalTPTree->Branch("gen_time", &gen_time, "gen_time[nGenParticles]/F");
 ecalTPTree->Branch("gen_time_max_02", &gen_time_max_02, "gen_time_max_02[nGenParticles]/F");
 ecalTPTree->Branch("gen_time_max_01", &gen_time_max_01, "gen_time_max_01[nGenParticles]/F");

 ecalTPTree->Branch("gParticleGrandMotherE", &gParticleGrandMotherE, "gParticleGrandMother[nGenParticles]/F");
 ecalTPTree->Branch("gParticleGrandMotherPt", &gParticleGrandMotherPt, "gParticleGrandMotherPt[nGenParticles]/F");
 ecalTPTree->Branch("gParticleGrandMotherPx", &gParticleGrandMotherPx, "gParticleGrandMotherPx[nGenParticles]/F");
 ecalTPTree->Branch("gParticleGrandMotherPy", &gParticleGrandMotherPy, "gParticleGrandMotherPy[nGenParticles]/F");
 ecalTPTree->Branch("gParticleGrandMotherPz", &gParticleGrandMotherPz, "gParticleGrandMotherPz[nGenParticles]/F");
 ecalTPTree->Branch("gParticleGrandMotherEta", &gParticleGrandMotherEta, "gParticleGrandMotherEta[nGenParticles]/F");
 ecalTPTree->Branch("gParticleGrandMotherPhi", &gParticleGrandMotherPhi, "gParticleGrandMotherPhi[nGenParticles]/F");
 ecalTPTree->Branch("gParticleGrandMotherDR", &gParticleGrandMotherDR, "gParticleGrandMotherDR[nGenParticles]/F");

 ecalTPTree->Branch("gParticleMotherE", &gParticleMotherE, "gParticleMother[nGenParticles]/F");
 ecalTPTree->Branch("gParticleMotherPt", &gParticleMotherPt, "gParticleMotherPt[nGenParticles]/F");
 ecalTPTree->Branch("gParticleMotherPx", &gParticleMotherPx, "gParticleMotherPx[nGenParticles]/F");
 ecalTPTree->Branch("gParticleMotherPy", &gParticleMotherPy, "gParticleMotherPy[nGenParticles]/F");
 ecalTPTree->Branch("gParticleMotherPz", &gParticleMotherPz, "gParticleMotherPz[nGenParticles]/F");
 ecalTPTree->Branch("gParticleMotherEta", &gParticleMotherEta, "gParticleMotherEta[nGenParticles]/F");
 ecalTPTree->Branch("gParticleMotherPhi", &gParticleMotherPhi, "gParticleMotherPhi[nGenParticles]/F");
 ecalTPTree->Branch("gParticleMotherDR", &gParticleMotherDR, "gParticleMotherDR[nGenParticles]/F");

 ecalTPTree->Branch("gParticleSiblingE", &gParticleSiblingE, "gParticleSiblingE[nGenParticles]/F");
 ecalTPTree->Branch("gParticleSiblingPt", &gParticleSiblingPt, "gParticleSiblingPt[nGenParticles]/F");
 ecalTPTree->Branch("gParticleSiblingPx", &gParticleSiblingPx, "gParticleSiblingPx[nGenParticles]/F");
 ecalTPTree->Branch("gParticleSiblingPy", &gParticleSiblingPy, "gParticleSiblingPy[nGenParticles]/F");
 ecalTPTree->Branch("gParticleSiblingPz", &gParticleSiblingPz, "gParticleSiblingPz[nGenParticles]/F");
 ecalTPTree->Branch("gParticleSiblingEta", &gParticleSiblingEta, "gParticleSiblingEta[nGenParticles]/F");
 ecalTPTree->Branch("gParticleSiblingPhi", &gParticleSiblingPhi, "gParticleSiblingPhi[nGenParticles]/F");
 ecalTPTree->Branch("gParticleSiblingDR", &gParticleSiblingDR, "gParticleSiblingDR[nGenParticles]/F");

};

void phase2L1EcalTimingAnalyzer::enableCaloJetBranches(){
 //jet info

 ecalTPTree->Branch("nCaloJets", &nCaloJets, "nCaloJets/I");
 ecalTPTree->Branch("caloJetE", &caloJetE, "caloJetE[nCaloJets]/F");
 ecalTPTree->Branch("caloJetPt", &caloJetPt, "caloJetPt[nCaloJets]/F");
 ecalTPTree->Branch("caloJetEta", &caloJetEta, "caloJetEta[nCaloJets]/F");
 ecalTPTree->Branch("caloJetPhi", &caloJetPhi, "caloJetPhi[nCaloJets]/F");
 ecalTPTree->Branch("CaloJetTime", &caloJetTime, "CaloJetTime[nCaloJets]/F");
 ecalTPTree->Branch("CaloJetEsum", &caloJetEsum, "CaloJetEsum[nCaloJets]/F");
 ecalTPTree->Branch("CaloJetTime_t", &caloJetTime_t, "CaloJetTime_t[nCaloJets]/F");
 ecalTPTree->Branch("CaloJetTime_t_60ps", &caloJetTime_t_60ps, "CaloJetTime_t_60ps[nCaloJets]/F");
 ecalTPTree->Branch("CaloJetTime_t_100ps", &caloJetTime_t_100ps, "CaloJetTime_t_100ps[nCaloJets]/F");
 ecalTPTree->Branch("CaloJetTime_t_60ps_func", &caloJetTime_t_60ps_func, "CaloJetTime_t_60ps_func[nCaloJets]/F");
 ecalTPTree->Branch("CaloJetTime_t_100ps_func", &caloJetTime_t_100ps_func, "CaloJetTime_t_100ps_func[nCaloJets]/F");
 ecalTPTree->Branch("CaloJetEsum_t", &caloJetEsum_t, "CaloJetEsum_t[nCaloJets]/F");

};


void phase2L1EcalTimingAnalyzer::enableGenak4JetNoNuBranches(){
 //jet info
 ecalTPTree->Branch("nGenak4JetNoNus", &nGenak4JetNoNus, "nGenak4JetNoNus/I");

 ecalTPTree->Branch("gak4JetNoNuTime", &gak4JetNoNuTime, "gak4JetNoNuTime[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuEsum", &gak4JetNoNuEsum, "gak4JetNoNuEsum[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuTime_t", &gak4JetNoNuTime_t, "gak4JetNoNuTime_t[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuTime_t_60ps", &gak4JetNoNuTime_t_60ps, "gak4JetNoNuTime_t_60ps[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuTime_t_100ps", &gak4JetNoNuTime_t_100ps, "gak4JetNoNuTime_t_100ps[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuTime_t_60ps_func", &gak4JetNoNuTime_t_60ps_func, "gak4JetNoNuTime_t_60ps_func[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuTime_t_100ps_func", &gak4JetNoNuTime_t_100ps_func, "gak4JetNoNuTime_t_100ps_func[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuEsum_t", &gak4JetNoNuEsum_t, "gak4JetNoNuEsum_t[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuTime_t0p5", &gak4JetNoNuTime_t0p5, "gak4JetNoNuTime_t0p5[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuEsum_t0p5", &gak4JetNoNuEsum_t0p5, "gak4JetNoNuEsum_t0p5[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuTime_t1", &gak4JetNoNuTime_t1, "gak4JetNoNuTime_t1[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuEsum_t1", &gak4JetNoNuEsum_t1, "gak4JetNoNuEsum_t1[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuTime_t2", &gak4JetNoNuTime_t2, "gak4JetNoNuTime_t2[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuEsum_t2", &gak4JetNoNuEsum_t2, "gak4JetNoNuEsum_t2[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuTime_t3", &gak4JetNoNuTime_t3, "gak4JetNoNuTime_t3[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuEsum_t3", &gak4JetNoNuEsum_t3, "gak4JetNoNuEsum_t3[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuTime_t5", &gak4JetNoNuTime_t5, "gak4JetNoNuTime_t5[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuEsum_t5", &gak4JetNoNuEsum_t5, "gak4JetNoNuEsum_t5[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuTime_t10", &gak4JetNoNuTime_t10, "gak4JetNoNuTime_t10[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuEsum_t10", &gak4JetNoNuEsum_t10, "gak4JetNoNuEsum_t10[nGenak4JetNoNus]/F");

 ecalTPTree->Branch("gak4JetNoNuMass", &gak4JetNoNuMass, "gak4JetNoNuMass[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuE", &gak4JetNoNuE, "gak4JetNoNuE[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuEt", &gak4JetNoNuEt, "gak4JetNoNuEt[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuPt", &gak4JetNoNuPt, "gak4JetNoNuPt[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuPx", &gak4JetNoNuPx, "gak4JetNoNuPx[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuPy", &gak4JetNoNuPy, "gak4JetNoNuPy[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuPz", &gak4JetNoNuPz, "gak4JetNoNuPz[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuEta", &gak4JetNoNuEta, "gak4JetNoNuEta[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuPhi", &gak4JetNoNuPhi, "gak4JetNoNuPhi[nGenak4JetNoNus]/F");

 ecalTPTree->Branch("gak4JetNoNuArea", &gak4JetNoNuArea, "gak4JetNoNuArea[nGenak4JetNoNus]/F");

 ecalTPTree->Branch("gak4JetNoNuPileupE", &gak4JetNoNuPileupE, "gak4JetNoNuPileupE[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuPileupIdFlag", &gak4JetNoNuPileupIdFlag, "gak4JetNoNuPileupIdFlag[nGenak4JetNoNus]/I");

 ecalTPTree->Branch("gak4JetNoNuPassIdLoose", &gak4JetNoNuPassIdLoose, "gak4JetNoNuPassIdLoose[nGenak4JetNoNus]/O");
 ecalTPTree->Branch("gak4JetNoNuPassIdTight", &gak4JetNoNuPassIdTight, "gak4JetNoNuPassIdTight[nGenak4JetNoNus]/O");


 ecalTPTree->Branch("gak4JetNoNuMuEnergy", &gak4JetNoNuMuEnergy, "gak4JetNoNuMuEnergy[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuEmEnergy", &gak4JetNoNuEmEnergy, "gak4JetNoNuEmEnergy[nGenak4JetNoNus]/F");
 //ecalTPTree->Branch("gak4JetNoNuChargedEmEnergy", &gak4JetNoNuChargedEmEnergy, "gak4JetNoNuChargedEmEnergy[nGenak4JetNoNus]/F");
 //ecalTPTree->Branch("gak4JetNoNuNeutralEmEnergy", &gak4JetNoNuNeutralEmEnergy, "gak4JetNoNuNeutralEmEnergy[nGenak4JetNoNus]/F");
 ecalTPTree->Branch("gak4JetNoNuHadronEnergy", &gak4JetNoNuHadronEnergy, "gak4JetNoNuHadronEnergy[nGenak4JetNoNus]/F");
 //ecalTPTree->Branch("gak4JetNoNuChargedHadronEnergy", &gak4JetNoNuChargedHadronEnergy, "gak4JetNoNuChargedHadronEnergy[nGenak4JetNoNus]/F");
 //ecalTPTree->Branch("gak4JetNoNuNeutralHadronEnergy", &gak4JetNoNuNeutralHadronEnergy, "gak4JetNoNuNeutralHadronEnergy[nGenak4JetNoNus]/F");

};




void phase2L1EcalTimingAnalyzer::enableGenak8JetBranches(){
 //jet info
 ecalTPTree->Branch("nGenak8Jets", &nGenak8Jets, "nGenak8Jets/I");

 ecalTPTree->Branch("gak8JetMass", &gak8JetMass, "gak8JetMass[nGenak8Jets]/F");
 ecalTPTree->Branch("gak8JetE", &gak8JetE, "gak8JetE[nGenak8Jets]/F");
 ecalTPTree->Branch("gak8JetEt", &gak8JetEt, "gak8JetEt[nGenak8Jets]/F");
 ecalTPTree->Branch("gak8JetPt", &gak8JetPt, "gak8JetPt[nGenak8Jets]/F");
 ecalTPTree->Branch("gak8JetPx", &gak8JetPx, "gak8JetPx[nGenak8Jets]/F");
 ecalTPTree->Branch("gak8JetPy", &gak8JetPy, "gak8JetPy[nGenak8Jets]/F");
 ecalTPTree->Branch("gak8JetPz", &gak8JetPz, "gak8JetPz[nGenak8Jets]/F");
 ecalTPTree->Branch("gak8JetEta", &gak8JetEta, "gak8JetEta[nGenak8Jets]/F");
 ecalTPTree->Branch("gak8JetPhi", &gak8JetPhi, "gak8JetPhi[nGenak8Jets]/F");

 ecalTPTree->Branch("gak8JetArea", &gak8JetArea, "gak8JetArea[nGenak8Jets]/F");

 ecalTPTree->Branch("gak8JetPileupE", &gak8JetPileupE, "gak8JetPileupE[nGenak8Jets]/F");
 ecalTPTree->Branch("gak8JetPileupIdFlag", &gak8JetPileupIdFlag, "gak8JetPileupIdFlag[nGenak8Jets]/I");

 ecalTPTree->Branch("gak8JetPassIdLoose", &gak8JetPassIdLoose, "gak8JetPassIdLoose[nGenak8Jets]/O");
 ecalTPTree->Branch("gak8JetPassIdTight", &gak8JetPassIdTight, "gak8JetPassIdTight[nGenak8Jets]/O");


 ecalTPTree->Branch("gak8JetMuEnergy", &gak8JetMuEnergy, "gak8JetMuEnergy[nGenak8Jets]/F");
 ecalTPTree->Branch("gak8JetEmEnergy", &gak8JetEmEnergy, "gak8JetEmEnergy[nGenak8Jets]/F");
 //ecalTPTree->Branch("gak8JetChargedEmEnergy", &gak8JetChargedEmEnergy, "gak8JetChargedEmEnergy[nGenak8Jets]/F");
 //ecalTPTree->Branch("gak8JetNeutralEmEnergy", &gak8JetNeutralEmEnergy, "gak8JetNeutralEmEnergy[nGenak8Jets]/F");
 ecalTPTree->Branch("gak8JetHadronEnergy", &gak8JetHadronEnergy, "gak8JetHadronEnergy[nGenak8Jets]/F");
 //ecalTPTree->Branch("gak8JetChargedHadronEnergy", &gak8JetChargedHadronEnergy, "gak8JetChargedHadronEnergy[nGenak8Jets]/F");
 //ecalTPTree->Branch("gak8JetNeutralHadronEnergy", &gak8JetNeutralHadronEnergy, "gak8JetNeutralHadronEnergy[nGenak8Jets]/F");

};

void phase2L1EcalTimingAnalyzer::enableGenak8JetNoNuBranches(){
 //jet info
 ecalTPTree->Branch("nGenak8JetNoNus", &nGenak8JetNoNus, "nGenak8JetNoNus/I");

 ecalTPTree->Branch("gak8JetNoNuMass", &gak8JetNoNuMass, "gak8JetNoNuMass[nGenak8JetNoNus]/F");
 ecalTPTree->Branch("gak8JetNoNuE", &gak8JetNoNuE, "gak8JetNoNuE[nGenak8JetNoNus]/F");
 ecalTPTree->Branch("gak8JetNoNuEt", &gak8JetNoNuEt, "gak8JetNoNuEt[nGenak8JetNoNus]/F");
 ecalTPTree->Branch("gak8JetNoNuPt", &gak8JetNoNuPt, "gak8JetNoNuPt[nGenak8JetNoNus]/F");
 ecalTPTree->Branch("gak8JetNoNuPx", &gak8JetNoNuPx, "gak8JetNoNuPx[nGenak8JetNoNus]/F");
 ecalTPTree->Branch("gak8JetNoNuPy", &gak8JetNoNuPy, "gak8JetNoNuPy[nGenak8JetNoNus]/F");
 ecalTPTree->Branch("gak8JetNoNuPz", &gak8JetNoNuPz, "gak8JetNoNuPz[nGenak8JetNoNus]/F");
 ecalTPTree->Branch("gak8JetNoNuEta", &gak8JetNoNuEta, "gak8JetNoNuEta[nGenak8JetNoNus]/F");
 ecalTPTree->Branch("gak8JetNoNuPhi", &gak8JetNoNuPhi, "gak8JetNoNuPhi[nGenak8JetNoNus]/F");

 ecalTPTree->Branch("gak8JetNoNuArea", &gak8JetNoNuArea, "gak8JetNoNuArea[nGenak8JetNoNus]/F");

 ecalTPTree->Branch("gak8JetNoNuPileupE", &gak8JetNoNuPileupE, "gak8JetNoNuPileupE[nGenak8JetNoNus]/F");
 ecalTPTree->Branch("gak8JetNoNuPileupIdFlag", &gak8JetNoNuPileupIdFlag, "gak8JetNoNuPileupIdFlag[nGenak8JetNoNus]/I");

 ecalTPTree->Branch("gak8JetNoNuPassIdLoose", &gak8JetNoNuPassIdLoose, "gak8JetNoNuPassIdLoose[nGenak8JetNoNus]/O");
 ecalTPTree->Branch("gak8JetNoNuPassIdTight", &gak8JetNoNuPassIdTight, "gak8JetNoNuPassIdTight[nGenak8JetNoNus]/O");


 ecalTPTree->Branch("gak8JetNoNuMuEnergy", &gak8JetNoNuMuEnergy, "gak8JetNoNuMuEnergy[nGenak8JetNoNus]/F");
 ecalTPTree->Branch("gak8JetNoNuEmEnergy", &gak8JetNoNuEmEnergy, "gak8JetNoNuEmEnergy[nGenak8JetNoNus]/F");
 //ecalTPTree->Branch("gak8JetNoNuChargedEmEnergy", &gak8JetNoNuChargedEmEnergy, "gak8JetNoNuChargedEmEnergy[nGenak8JetNoNus]/F");
 //ecalTPTree->Branch("gak8JetNoNuNeutralEmEnergy", &gak8JetNoNuNeutralEmEnergy, "gak8JetNoNuNeutralEmEnergy[nGenak8JetNoNus]/F");
 ecalTPTree->Branch("gak8JetNoNuHadronEnergy", &gak8JetNoNuHadronEnergy, "gak8JetNoNuHadronEnergy[nGenak8JetNoNus]/F");
 //ecalTPTree->Branch("gak8JetNoNuChargedHadronEnergy", &gak8JetNoNuChargedHadronEnergy, "gak8JetNoNuChargedHadronEnergy[nGenak8JetNoNus]/F");
 //ecalTPTree->Branch("gak8JetNoNuNeutralHadronEnergy", &gak8JetNoNuNeutralHadronEnergy, "gak8JetNoNuNeutralHadronEnergy[nGenak8JetNoNus]/F");

};


// ------------ reset branches  ------------
void phase2L1EcalTimingAnalyzer::resetBranches(){
  resetEventInfoBranches();
  resetEBCrystalBranches();
  resetGenParticleBranches();
  resetGenak4JetBranches();
  resetGenak4JetNoNuBranches();
  resetGenak8JetBranches();
  resetGenak8JetNoNuBranches();

  nCaloJets = 0;
  for(int i=0; i<GENJETARRAYSIZE; i++){
    caloJetTime[i] = -666.;
    caloJetEsum[i] = -666.;
    caloJetTime_t[i] = -666.;
    caloJetTime_t_60ps[i] = -666.;
    caloJetTime_t_100ps[i] = -666.;
    caloJetTime_t_60ps_func[i] = -666.;
    caloJetTime_t_100ps_func[i] = -666.;
    caloJetEsum_t[i] = -666.;
    caloJetE[i]    = -666.;
    caloJetPt[i]   = -666.;
    caloJetEta[i]  = -666.;
    caloJetPhi[i]  = -666.;
  }

};

void phase2L1EcalTimingAnalyzer::resetEventInfoBranches(){
 runNum=0;
 lumiNum=0;
 eventNum=0;
};

void phase2L1EcalTimingAnalyzer::resetEBCrystalBranches(){
 //ecal crystal info
 nCrystals=0; 

 for(int i=0; i<EBCRYSTALARRAYSIZE; i++){
 eb_Et[i]       = -666.;
 eb_Edep[i]     = -666.;
 eb_time[i]     = -666.;
 eb_time_60ps[i]     = -666.;
 eb_time_100ps[i]     = -666.;
 eb_time_60ps_func[i]     = -666.;
 eb_time_100ps_func[i]     = -666.;
 eb_sigmat[i]     = -666.;
 eb_sigmat_60ps_func[i]     = -666.;
 eb_sigmat_100ps_func[i]     = -666.;
 eb_id[i]       = -666;
 eb_ieta[i]     = -666;
 eb_iphi[i]     = -666;
 eb_ism[i]      = -666;
 eb_ic[i]       = -666;
 eb_cell_Eta[i] = -666.;
 eb_cell_Phi[i] = -666.;

 gen_time_index[i]     = -666;
 gen_time_e[i]     = -666.;
 gen_time_dr[i]     = -666.;
 gen_time_tp[i]     = -666.;
 gen_time_tp_60ps[i]     = -666.;
 gen_time_tp_100ps[i]     = -666.;
 gen_time_tp_60ps_func[i]     = -666.;
 gen_time_tp_100ps_func[i]     = -666.;
 }

};

void phase2L1EcalTimingAnalyzer::resetGenParticleBranches(){
 //gen info
 nGenParticles=0;

 genVertexX=0.;
 genVertexY=0.;
 genVertexZ=0.;
 genVertexT=0.;

 for(int i=0; i<GENPARTICLEARRAYSIZE; i++){

 gParticleGrandMotherId[i]     = -666;
 gParticleGrandMotherIndex[i]  = -666;

 gParticleMotherId[i]     = -666;
 gParticleMotherIndex[i]  = -666;

 gParticleSiblingId[i]     = -666;
 gParticleSiblingIndex[i]  = -666;

 gParticleId[i]           = -666;
 gParticleStatus[i]       = -666;

 gParticleE[i]            = -666.;
 gParticlePt[i]           = -666.;
 gParticlePx[i]           = -666.;
 gParticlePy[i]           = -666.;
 gParticlePz[i]           = -666.;
 gParticleEta[i]          = -666.;
 gParticlePhi[i]          = -666.;

 gParticle_decay_vtx_x[i] = -666.;
 gParticle_decay_vtx_y[i] = -666.;
 gParticle_decay_vtx_z[i] = -666.;

 gParticle_prod_vtx_x[i] = -666.;
 gParticle_prod_vtx_y[i] = -666.;
 gParticle_prod_vtx_z[i] = -666.;

 g_eb_time_Emax_02[i] = -666.;
 g_eb_sigmat_Emax_02[i] = -666.;
 gEmax_02[i] = -666.;
 gImax_02[i] = -666;
 g_dt_sc_02[i] = -666.;
 g_tmax_sc_02[i] = -666.;
 g_cnt_sc_02[i] = -666;
 g_tmin_sc_02[i] = -666.;
 g_dt_sm_02[i] = -666.;
 g_tmax_sm_02[i] = -666.;
 g_cnt_sm_02[i] = -666;
 g_tmin_sm_02[i] = -666.;
 gEsc_02[i] = -666.;
 genEsc_02[i] = -666.;
 gE9x9_02[i] = -666.;
 gE5x5_02[i] = -666.;
 gE3x3_02[i] = -666.;

 g_eb_time_Emax_01[i] = -666.;
 g_eb_sigmat_Emax_01[i] = -666.;
 gEmax_01[i] = -666.;
 gImax_01[i] = -666;
 g_dt_sc_01[i] = -666.;
 g_tmax_sc_01[i] = -666.;
 g_cnt_sc_01[i] = -666;
 g_tmin_sc_01[i] = -666.;
 g_dt_sm_01[i] = -666.;
 g_tmax_sm_01[i] = -666.;
 g_cnt_sm_01[i] = -666;
 g_tmin_sm_01[i] = -666.;
 g_dt_sm2_01[i] = -666.;
 g_tmax_sm2_01[i] = -666.;
 g_cnt_sm2_01[i] = -666;
 g_tmin_sm2_01[i] = -666.;
 g_dt_sm5_01[i] = -666.;
 g_tmax_sm5_01[i] = -666.;
 g_cnt_sm5_01[i] = -666;
 g_tmin_sm5_01[i] = -666.;
 g_dt_sm10_01[i] = -666.;
 g_tmax_sm10_01[i] = -666.;
 g_cnt_sm10_01[i] = -666;
 g_tmin_sm10_01[i] = -666.;
 g_dt_max_01[i] = -666.;
 g_t_max_01[i] = -666.;
 g_e_max_01[i] = -666.;
 g_dt_min_01[i] = -666.;
 g_t_min_01[i] = -666.;
 g_e_min_01[i] = -666.;
 gEsc_01[i] = -666.;
 genEsc_01[i] = -666.;
 gE9x9_01[i] = -666.;
 gE5x5_01[i] = -666.;
 gE3x3_01[i] = -666.;

 g_mindr[i] = -666.;
 g_e_mindr[i] = -666.;
 g_t_mindr[i] = -666.;

 gIcore[i] = -666;
 g_eb_time[i] = -666.;
 g_eb_sigmat[i] = -666.;
 g_tof[i] = -666.;
 g_tvirtual[i] = -666.;
 gen_time[i] = -666.;
 gen_time_max_02[i] = -666.;
 gen_time_max_01[i] = -666.;

 gParticleGrandMotherE[i]            = -666.;
 gParticleGrandMotherPt[i]           = -666.;
 gParticleGrandMotherPx[i]           = -666.;
 gParticleGrandMotherPy[i]           = -666.;
 gParticleGrandMotherPz[i]           = -666.;
 gParticleGrandMotherEta[i]          = -666.;
 gParticleGrandMotherPhi[i]          = -666.;
 gParticleGrandMotherDR[i]          = -666.;

 gParticleMotherE[i]            = -666.;
 gParticleMotherPt[i]           = -666.;
 gParticleMotherPx[i]           = -666.;
 gParticleMotherPy[i]           = -666.;
 gParticleMotherPz[i]           = -666.;
 gParticleMotherEta[i]          = -666.;
 gParticleMotherPhi[i]          = -666.;
 gParticleMotherDR[i]          = -666.;

 gParticleSiblingE[i]            = -666.;
 gParticleSiblingPt[i]           = -666.;
 gParticleSiblingPx[i]           = -666.;
 gParticleSiblingPy[i]           = -666.;
 gParticleSiblingPz[i]           = -666.;
 gParticleSiblingEta[i]          = -666.;
 gParticleSiblingPhi[i]          = -666.;
 gParticleSiblingDR[i]          = -666.;


 }

};


void phase2L1EcalTimingAnalyzer::resetGenak4JetBranches(){
 // ak4 jet info
 nGenak4Jets = 0;

 for(int i=0; i<GENJETARRAYSIZE; i++){
 gak4JetMass[i] = -666.;
 gak4JetE[i]    = -666.;
 gak4JetEt[i]   = -666.;
 gak4JetPt[i]   = -666.;
 gak4JetPx[i]   = -666.;
 gak4JetPy[i]   = -666.;
 gak4JetPz[i]   = -666.;
 gak4JetEta[i]  = -666.;
 gak4JetPhi[i]  = -666.;

 gak4JetArea[i] = -666.;

 gak4JetPileupE[i] = -666. ;
 gak4JetPileupIdFlag[i] = -666;

 gak4JetPassIdLoose[i] = false;
 gak4JetPassIdTight[i] = false;

 gak4JetMuEnergy[i] = -666.;
 gak4JetEmEnergy[i] =-666.;
 gak4JetChargedEmEnergy[i] =-666.;
 gak4JetNeutralEmEnergy[i] = -666.;
 gak4JetHadronEnergy[i] = -666.;
 gak4JetChargedHadronEnergy[i] = -666.;
 gak4JetNeutralHadronEnergy[i] = -666.;
 }

};

void phase2L1EcalTimingAnalyzer::resetGenak4JetNoNuBranches(){
 // ak4 jet nonu info
 nGenak4JetNoNus = 0;

 for(int i=0; i<GENJETARRAYSIZE; i++){
 gak4JetNoNuTime[i] = -666.;
 gak4JetNoNuEsum[i] = -666.;
 gak4JetNoNuTime_t[i] = -666.;
 gak4JetNoNuTime_t_60ps[i] = -666.;
 gak4JetNoNuTime_t_100ps[i] = -666.;
 gak4JetNoNuTime_t_60ps_func[i] = -666.;
 gak4JetNoNuTime_t_100ps_func[i] = -666.;
 gak4JetNoNuEsum_t[i] = -666.;
 gak4JetNoNuTime_t0p5[i] = -666.;
 gak4JetNoNuEsum_t0p5[i] = -666.;
 gak4JetNoNuTime_t1[i] = -666.;
 gak4JetNoNuEsum_t1[i] = -666.;
 gak4JetNoNuTime_t2[i] = -666.;
 gak4JetNoNuEsum_t2[i] = -666.;
 gak4JetNoNuTime_t3[i] = -666.;
 gak4JetNoNuEsum_t3[i] = -666.;
 gak4JetNoNuTime_t5[i] = -666.;
 gak4JetNoNuEsum_t5[i] = -666.;
 gak4JetNoNuTime_t10[i] = -666.;
 gak4JetNoNuEsum_t10[i] = -666.;

 gak4JetNoNuMass[i] = -666.;
 gak4JetNoNuE[i]    = -666.;
 gak4JetNoNuEt[i]   = -666.;
 gak4JetNoNuPt[i]   = -666.;
 gak4JetNoNuPx[i]   = -666.;
 gak4JetNoNuPy[i]   = -666.;
 gak4JetNoNuPz[i]   = -666.;
 gak4JetNoNuEta[i]  = -666.;
 gak4JetNoNuPhi[i]  = -666.;

 gak4JetNoNuArea[i] = -666.;

 gak4JetNoNuPileupE[i] = -666. ;
 gak4JetNoNuPileupIdFlag[i] = -666;

 gak4JetNoNuPassIdLoose[i] = false;
 gak4JetNoNuPassIdTight[i] = false;

 gak4JetNoNuMuEnergy[i] = -666.;
 gak4JetNoNuEmEnergy[i] =-666.;
 gak4JetNoNuChargedEmEnergy[i] =-666.;
 gak4JetNoNuNeutralEmEnergy[i] = -666.;
 gak4JetNoNuHadronEnergy[i] = -666.;
 gak4JetNoNuChargedHadronEnergy[i] = -666.;
 gak4JetNoNuNeutralHadronEnergy[i] = -666.;
 }

};

void phase2L1EcalTimingAnalyzer::resetGenak8JetBranches(){
 //ak8 jet info
 nGenak8Jets = 0;

 for(int i=0; i<GENJETARRAYSIZE; i++){
 gak8JetMass[i] = -666.;
 gak8JetE[i]    = -666.;
 gak8JetEt[i]   = -666.;
 gak8JetPt[i]   = -666.;
 gak8JetPx[i]   = -666.;
 gak8JetPy[i]   = -666.;
 gak8JetPz[i]   = -666.;
 gak8JetEta[i]  = -666.;
 gak8JetPhi[i]  = -666.;

 gak8JetArea[i] = -666.;

 gak8JetPileupE[i] = -666. ;
 gak8JetPileupIdFlag[i] = -666;

 gak8JetPassIdLoose[i] = false;
 gak8JetPassIdTight[i] = false;

 gak8JetMuEnergy[i] = -666.;
 gak8JetEmEnergy[i] =-666.;
 gak8JetChargedEmEnergy[i] =-666.;
 gak8JetNeutralEmEnergy[i] = -666.;
 gak8JetHadronEnergy[i] = -666.;
 gak8JetChargedHadronEnergy[i] = -666.;
 gak8JetNeutralHadronEnergy[i] = -666.;
 }

};

void phase2L1EcalTimingAnalyzer::resetGenak8JetNoNuBranches(){
 //ak8 jet nonu info
 nGenak8JetNoNus = 0;

 for(int i=0; i<GENJETARRAYSIZE; i++){
 gak8JetNoNuMass[i] = -666.;
 gak8JetNoNuE[i]    = -666.;
 gak8JetNoNuEt[i]   = -666.;
 gak8JetNoNuPt[i]   = -666.;
 gak8JetNoNuPx[i]   = -666.;
 gak8JetNoNuPy[i]   = -666.;
 gak8JetNoNuPz[i]   = -666.;
 gak8JetNoNuEta[i]  = -666.;
 gak8JetNoNuPhi[i]  = -666.;

 gak8JetNoNuArea[i] = -666.;

 gak8JetNoNuPileupE[i] = -666. ;
 gak8JetNoNuPileupIdFlag[i] = -666;

 gak8JetNoNuPassIdLoose[i] = false;
 gak8JetNoNuPassIdTight[i] = false;

 gak8JetNoNuMuEnergy[i] = -666.;
 gak8JetNoNuEmEnergy[i] =-666.;
 gak8JetNoNuChargedEmEnergy[i] =-666.;
 gak8JetNoNuNeutralEmEnergy[i] = -666.;
 gak8JetNoNuHadronEnergy[i] = -666.;
 gak8JetNoNuChargedHadronEnergy[i] = -666.;
 gak8JetNoNuNeutralHadronEnergy[i] = -666.;
 }

};


// ------------ fill branches  ------------
bool phase2L1EcalTimingAnalyzer::fillEventInfoBranches(const edm::Event& iEvent){
  //event info
  runNum   = iEvent.id().run();
  lumiNum  = iEvent.id().luminosityBlock();
  eventNum = iEvent.id().event();

  return true;
};

bool phase2L1EcalTimingAnalyzer::fillEBCrystalBranches(const edm::Event& iEvent, const edm::EventSetup& iSetup){
//bool phase2L1EcalTimingAnalyzer::fillEBCrystalBranches(const edm::Event& iEvent, const edm::EventSetup& iSetup, std::vector<reco::GenParticle> genParticles){

  std::vector<reco::GenParticle>  genParticles = phase2L1EcalTimingAnalyzer::GetGenParticles();

  edm::ESHandle<CaloGeometry> caloGeometryHandle;
  iSetup.get<CaloGeometryRecord>().get(caloGeometryHandle);
  const CaloGeometry* caloGeometry_ = caloGeometryHandle.product();  
  ebGeometry = caloGeometry_->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);

  //ecal barrel crystals
  for(auto& tpg : *ecaltpgCollection.product()) {

    nCrystals++; 

    //std::cout<<"Et " << tpg.encodedEt()<<std::endl;
    //std::cout<<"l1aSpike " << tpg.l1aSpike()<<std::endl;
    if(tpg.time()!=0 && debug) {
	std::cout<<"time " << tpg.time()<<std::endl;
	std::cout<<"id " << tpg.id()<<std::endl;
	//std::cout<<"id iEta " << tpg.id().ieta() << " iPhi " << tpg.id().iphi()<<std::endl;
	//std::cout<<"tower iEta " << tpg.id().tower_ieta() << " iPhi " << tpg.id().tower_iphi()<<std::endl;
	std::cout<<"approx Eta " << tpg.id().approxEta() <<std::endl;
	std::cout<<"time above is not zero " <<std::endl;
    }

    eb_id[nCrystals-1] = tpg.id();
    eb_ieta[nCrystals-1] = tpg.id().ieta();
    eb_iphi[nCrystals-1] = tpg.id().iphi();
    eb_ism[nCrystals-1] = tpg.id().ism();
    eb_ic[nCrystals-1] = tpg.id().ic();
    
    GlobalVector position;
    auto cell = ebGeometry->getGeometry(tpg.id());
    float eta = cell->getPosition().eta();
    float phi = cell->getPosition().phi();
    float cellX = cell->getPosition().x();
    float cellY = cell->getPosition().y();
    float cellZ = cell->getPosition().z();
    eb_cell_Eta[nCrystals-1] = eta;
    eb_cell_Phi[nCrystals-1] = phi;
    
    //std::cout<<"id " << tpg.id().ieta() << tpg.id().iphi() << tpg.id().ism() << tpg.id().ic() <<std::endl;
    //std::cout<<"approx Eta " << tpg.id().approxEta() <<std::endl;
    //std::cout<<"cell Eta " << eta << " Phi " << phi <<std::endl;
    //std::cout<<"Et above is not zero " <<std::endl;
    
    if(tpg.encodedEt() > 0) 
      {
	float et = tpg.encodedEt()/8.; // convert ADC to GeV
	
	//if(et<0.001) continue;
	//float energy = et / sin(position.theta());
	eb_Et[nCrystals-1] = et;
	//eb_Edep[nCrystals-1] = Et_to_E(et, tpg.id().ieta());
	float e = Et_to_E(et, tpg.id().ieta());
	eb_Edep[nCrystals-1] = e;
	//eb_time[nCrystals-1] = tpg.time();
	float sigma_t = pow(e, -1)*0.540; //similar to eta=0 lumi 300/fb, e = 20 GeV, sigma_t = 43 ps
	float sigma_t_60ps_func = pow(e, -1)*1.2; // e = 20 GeV, sigma_t = 60 ps
	float sigma_t_100ps_func = pow(e, -1)*2.0; // e = 20 GeV, sigma_t = 100 ps
	std::random_device rd;
	//std::default_random_engine generator;
	std::mt19937 generator(rd());
	std::normal_distribution<float> distribution(0,sigma_t);
	std::normal_distribution<float> distribution_60ps(0, 0.06);
	std::normal_distribution<float> distribution_100ps(0, 0.1);
	std::normal_distribution<float> distribution_60ps_func(0,sigma_t_60ps_func);
	std::normal_distribution<float> distribution_100ps_func(0,sigma_t_100ps_func);
	eb_time[nCrystals-1] = distribution(generator);  // generates number
	eb_time_60ps[nCrystals-1] = distribution_60ps(generator);  // generates number
	eb_time_100ps[nCrystals-1] = distribution_100ps(generator);  // generates number
	eb_time_60ps_func[nCrystals-1] = distribution_60ps_func(generator);  // generates number
	eb_time_100ps_func[nCrystals-1] = distribution_100ps_func(generator);  // generates number
	eb_sigmat[nCrystals-1] = sigma_t;  // generates number
	eb_sigmat_60ps_func[nCrystals-1] = sigma_t_60ps_func;  // generates number
	eb_sigmat_100ps_func[nCrystals-1] = sigma_t_100ps_func;  // generates number

	//std::cout<<eb_time[nCrystals-1] <<std::endl;
	//float sigma_t = TMath::pow(et,-0.34242268082)*0.121
	//std::cout<<"time " << tpg.time()<<std::endl;

	//std::cout << "check: " << e << " : " << sigma_t << " | " << eb_time[nCrystals-1] << "\n";

      }
  	
    float unit = 2*TMath::Pi()/360;
    float time_tp = -666.;
    float time_tp_60ps = -666.;
    float time_tp_100ps = -666.;
    float time_tp_60ps_func = -666.;
    float time_tp_100ps_func = -666.;
    float deltar = 666.;
    float distance = 666.;
    int index = -666;
    float e = -666.;
    float t_crystal = 0.;
    float t_crystal_60ps = 0.;
    float t_crystal_100ps = 0.;
    float t_crystal_60ps_func = 0.;
    float t_crystal_100ps_func = 0.;
    for(unsigned int i = 0; i< genParticles.size(); i++){

      if(gParticleId[i]==111 && abs(gParticleEta[i]<1.5) )
	{
	  reco::GenParticle gen = genParticles[i];

	  vector<float> etaphi = EtaPhi_Corr_EB(gParticle_prod_vtx_x[i], gParticle_prod_vtx_y[i], gParticle_prod_vtx_z[i], gen);
	  float gen_eta = etaphi[0];
	  float gen_phi = etaphi[1];
	  float gen_tof = etaphi[2];
	  float gen_tvirtual = etaphi[3];

	  //std::cout<<gen_eta <<" , "<< eta <<std::endl;
	  //dr between pi0 and crystal
	  distance = deltaR(gen_eta, gen_phi, eta, phi);

	  if(distance<2*unit)
	    {
	      //std::cout<<distance <<" , "<< dr <<std::endl;
	      if(distance <= deltar)
		{
		  index = i;
		  e = gen.energy();
		  deltar = deltaR(gen_eta, gen_phi, eta, phi);
		  if(eb_time[nCrystals-1]==-666) t_crystal = 0.;
		  else t_crystal = eb_time[nCrystals-1];
		  time_tp = gen_tof + genVertexT + t_crystal - gen_tvirtual;

		  if(eb_time_60ps[nCrystals-1]==-666) t_crystal_60ps = 0.;
		  else t_crystal_60ps = eb_time_60ps[nCrystals-1];
		  time_tp_60ps = gen_tof + genVertexT + t_crystal_60ps - gen_tvirtual;

		  if(eb_time_100ps[nCrystals-1]==-666) t_crystal_100ps = 0.;
		  else t_crystal_100ps = eb_time_100ps[nCrystals-1];
		  time_tp_100ps = gen_tof + genVertexT + t_crystal_100ps - gen_tvirtual;

		  if(eb_time_60ps_func[nCrystals-1]==-666) t_crystal_60ps_func = 0.;
		  else t_crystal_60ps_func = eb_time_60ps_func[nCrystals-1];
		  time_tp_60ps_func = gen_tof + genVertexT + t_crystal_60ps_func - gen_tvirtual;

		  if(eb_time_100ps_func[nCrystals-1]==-666) t_crystal_100ps_func = 0.;
		  else t_crystal_100ps_func = eb_time_100ps_func[nCrystals-1];
		  time_tp_100ps_func = gen_tof + genVertexT + t_crystal_100ps_func - gen_tvirtual;
		  //std::cout<<distance <<" , "<< dr << " , time_tp "<< time_tp << " , gen_tof " << gen_tof << " , genVertexT " << genVertexT << ", eb_time[nCrystal-1] "<< eb_time[nCrystals-1] << ", gen_tvirtual " << gen_tvirtual<<std::endl;
					
		}	

	    }
			
	}
    }//match TP to close/fast pi0
    //std::cout<<distance <<" , "<< dr <<std::endl;

    //if(deltar = 666) time_tp = -666;
	
    //if we cannot find a pi0, then use time of an imaginary photon
    if (time_tp == -666 && tpg.encodedEt() > 0) {
      vector<float> timeVector = CalculateTime_EB(genVertexX, genVertexY, genVertexZ, cellX, cellY, cellZ);
      float gen_tof = timeVector[0];
      float gen_tvirtual = timeVector[1];
      t_crystal = eb_time[nCrystals-1];
      time_tp = gen_tof + genVertexT + t_crystal - gen_tvirtual;		        
    }


    gen_time_index[nCrystals-1] = index;
    gen_time_e[nCrystals-1] = e;
    gen_time_dr[nCrystals-1] = deltar;
    gen_time_tp[nCrystals-1] = time_tp;
    gen_time_tp_60ps[nCrystals-1] = time_tp_60ps;
    gen_time_tp_100ps[nCrystals-1] = time_tp_100ps;
    gen_time_tp_60ps_func[nCrystals-1] = time_tp_60ps_func;
    gen_time_tp_100ps_func[nCrystals-1] = time_tp_100ps_func;
    }
  return true;
};

bool phase2L1EcalTimingAnalyzer::fillGenParticleBranches(){ 
  std::vector<reco::GenParticle>  genParticles = phase2L1EcalTimingAnalyzer::GetGenParticles();
  phase2L1EcalTimingAnalyzer::fillGenParticleBasicBranches( genParticles);
  //  phase2L1EcalTimingAnalyzer::fillGenParticleMotherBranches( genParticles);
  //  phase2L1EcalTimingAnalyzer::fillGenParticleGrandMotherBranches( genParticles);
  //  //phase2L1EcalTimingAnalyzer::fillGenParticleSiblingBranches( genParticles);
  //  phase2L1EcalTimingAnalyzer::fillGenParticleTPBranches( genParticles);
  return true;
};
//  //fill gen info
//  
std::vector<reco::GenParticle>  phase2L1EcalTimingAnalyzer::GetGenParticles(){
  // Get genParticles
  std::vector<reco::GenParticle> genParticles;
  

  bool foundGenVertex = false;
  int pi = 0;
  for(unsigned int i = 0; i< genParticleHandle->size(); i++){

    edm::Ptr<reco::GenParticle> ptr(genParticleHandle, i);
    if(
       (abs((ptr)->pdgId()) == 2212 ) //protons
       || (abs((ptr)->pdgId()) >= 1 && abs((ptr)->pdgId()) <= 6 )//&& ( (ptr)->status() < 30 )) //quarks
       || (abs((ptr)->pdgId()) >= 11 && abs((ptr)->pdgId()) <= 16) //leptons
       || (abs((ptr)->pdgId()) == 21) // && (ptr)->status() < 30) //gluons
       || (abs((ptr)->pdgId()) >= 22 && abs((ptr)->pdgId()) <= 25) //&& ( (ptr)->status() < 30)) //gammas, z0, w, higgs
       || (abs((ptr)->pdgId()) >= 32 && abs((ptr)->pdgId()) <= 42) // other gauge and higgs bosons
       || (abs((ptr)->pdgId()) == 111 || abs((ptr)->pdgId()) == 211) // pi0 or pi+ pi-
       || (abs((ptr)->pdgId()) == 311 || abs((ptr)->pdgId()) == 321) // k0 or k+ k-
       || (abs((ptr)->pdgId()) == 130 || abs((ptr)->pdgId()) == 310) // k0l or k0s
       || (abs((ptr)->pdgId()) >= 1000001 && abs((ptr)->pdgId()) <= 1000039) //susy particles
       || (abs((ptr)->pdgId()) == 9000006 || abs((ptr)->pdgId()) == 9000007) //llp
     ){
    genParticles.push_back(*ptr);
/*
    reco::GenParticle gen = *ptr;
    gParticleId[pi] = gen.pdgId();
    pi++;
*/
    }
    //genParticles.push_back(*ptr);
  //std::cout<<"Finished push back gen particles"<<std::endl;

    int num = ptr->numberOfDaughters();
    if (!foundGenVertex)
    {
      for (int j=0; j< num; ++j)
      {
        const reco::Candidate *dau = ptr->daughter(j);
        if (dau)
        {
          genVertexX = dau->vx();
	  //std::cout<<"vx " << dau->vx() << " genVertexX " << genVertexX <<std::endl;
          genVertexY = dau->vy();
          genVertexZ = dau->vz();
          foundGenVertex = true;
          break;
        }
      }
    }//find genvertex
 

  }

  genVertexT = *genVertexTHandle;

  return genParticles;
};
//
//
bool phase2L1EcalTimingAnalyzer::fillGenParticleBasicBranches(std::vector<reco::GenParticle> genParticles){
//  //std::cout<<"Finished  fill gen vtx"<<std::endl;

  nGenParticles = genParticles.size();
  for(unsigned int i = 0; i< genParticles.size(); i++){
    reco::GenParticle gen = genParticles[i];

    gParticleId[i] = gen.pdgId();
    gParticleStatus[i] = gen.status();

    gParticleE[i] = gen.energy();
    gParticlePt[i] = gen.pt();
    gParticlePx[i] = gen.px();
    gParticlePy[i] = gen.py();
    gParticlePz[i] = gen.pz();
    gParticleEta[i] = gen.eta();
    gParticlePhi[i] = gen.phi();

    gParticle_decay_vtx_x[i] = gen.vx();
    gParticle_decay_vtx_y[i] = gen.vy();
    gParticle_decay_vtx_z[i] = gen.vz();

    gParticle_prod_vtx_x[i] = gen.vx();
    gParticle_prod_vtx_y[i] = gen.vy();
    gParticle_prod_vtx_z[i] = gen.vz();
//   }
//
//   return true;
//};
//  //std::cout<<"Finished  fill  basic gen particles"<<std::endl;
//
//bool phase2L1EcalTimingAnalyzer::fillGenParticleMotherBranches(std::vector<reco::GenParticle> genParticles){
//  for(unsigned int i = 0; i< genParticles.size(); i++){
//    reco::GenParticle gen = genParticles[i];
    // gen mother
    if(gen.numberOfMothers() > 0){

//	const reco::Candidate* mother = gen.mother(0);
//	//const reco::Candidate mother = *mo;
//
//        			gParticleMotherId[i] = mother->pdgId();
//	for(unsigned int j = 0; j < genParticles.size(); j++)
//	{
//		if(genParticles[j].pdgId() == mother->pdgId() && genParticles[j].energy() == mother->energy() && genParticles[j].eta() ==mother->eta() && genParticles[j].phi() == mother->phi() && genParticles[j].pt() == mother->pt() && genParticles[j].px() == mother->px() )
//			gParticleMotherIndex[i] = j;
//        }
//        			gParticleMotherE[i] = mother->energy();
//        			gParticleMotherPt[i] = mother->pt();
//        			gParticleMotherPx[i] = mother->px();
//        			gParticleMotherPy[i] = mother->py();
//        			gParticleMotherPz[i] = mother->pz();
//        			gParticleMotherEta[i] = mother->eta();
//        			gParticleMotherPhi[i] = mother->phi();
//        			gParticleMotherDR[i] = deltaR(mother->eta(), mother->phi(), genParticles[i].eta(), genParticles[i].phi());

//careful version
	const reco::Candidate *firstMotherWithDifferentID = findFirstMotherWithDifferentID(&gen);
        if (firstMotherWithDifferentID)
        {
        	gParticleMotherId[i] = firstMotherWithDifferentID->pdgId();
        }

	//find the mother and keep going up the mother chain if the ID's are the same
	const reco::Candidate *originalMotherWithSameID = findOriginalMotherWithSameID(&gen);
	for(unsigned int j = 0; j < genParticles.size(); j++)
	{
		//std::cout<<" Particle GEN Id " << gen.pdgId() << " MotherId " << gParticleMotherId[i] << " origin mather id  " << originalMotherWithSameID->pdgId() << " test id " << genParticles[j].pdgId() 
		//<< " pointer &i "  << &gen << " &o " << &originalMotherWithSameID << " &j " << &(genParticles[j]) 
		//<< " o  " << originalMotherWithSameID << " j " << genParticles[j] 
		//<< " *j " << *genParticles[j] 
		//<< " *o " << *originalMotherWithSameID 
		//<<std::endl;
		if(genParticles[j].pdgId() == originalMotherWithSameID->pdgId() && genParticles[j].status() == originalMotherWithSameID->status() && genParticles[j].energy() == originalMotherWithSameID->energy() && genParticles[j].eta() ==originalMotherWithSameID->eta() && genParticles[j].phi() == originalMotherWithSameID->phi() && genParticles[j].pt() == originalMotherWithSameID->pt() && genParticles[j].px() == originalMotherWithSameID->px() )
		{
			//std::cout<<" find "<<std::endl;

			gParticleMotherIndex[i] = j;

			//mother info
			int mo = gParticleMotherIndex[i];
			if(mo!=-666){
    				reco::GenParticle mother = genParticles[mo];

        			gParticleMotherE[i] = mother.energy();
        			gParticleMotherPt[i] = mother.pt();
        			gParticleMotherPx[i] = mother.px();
        			gParticleMotherPy[i] = mother.py();
        			gParticleMotherPz[i] = mother.pz();
        			gParticleMotherEta[i] = mother.eta();
        			gParticleMotherPhi[i] = mother.phi();
        			gParticleMotherDR[i] = deltaR(mother.eta(), mother.phi(), genParticles[i].eta(), genParticles[i].phi());

			}
	
			break;
		}
	}
     }// finish gen mother
//   }
//   return true;
//};
//
//  //std::cout<<"Finished  fill gen particles mother"<<std::endl;
//
//bool phase2L1EcalTimingAnalyzer::fillGenParticleGrandMotherBranches(std::vector<reco::GenParticle> genParticles){
//  for(unsigned int i = 0; i< genParticles.size(); i++){
//    reco::GenParticle gen = genParticles[i];
	  int k1 = gParticleMotherIndex[i]; //mother index

	  int k11 = gParticleMotherIndex[k1]; // grand mother index

	 //grand mother info
	 if(k1!=-666 && k11 != -666){
	 	reco::GenParticle grandmom = genParticles[k11];

		gParticleGrandMotherId[i] = genParticles[k11].pdgId();
		gParticleGrandMotherIndex[i] = k11;

        	gParticleGrandMotherE[i] = grandmom.energy();
        	gParticleGrandMotherPt[i] = grandmom.pt();
        	gParticleGrandMotherPx[i] = grandmom.px();
        	gParticleGrandMotherPy[i] = grandmom.py();
        	gParticleGrandMotherPz[i] = grandmom.pz();
        	gParticleGrandMotherEta[i] = grandmom.eta();
        	gParticleGrandMotherPhi[i] = grandmom.phi();
        	gParticleGrandMotherDR[i] = deltaR(grandmom.eta(), grandmom.phi(), genParticles[i].eta(), genParticles[i].phi());
		
	 }//finish grand mother
//   }
//   return true;
//};
//
//  //std::cout<<"Finished  fill gen particles grand mother"<<std::endl;
//
//bool phase2L1EcalTimingAnalyzer::fillGenParticleSiblingBranches(std::vector<reco::GenParticle> genParticles){
//  for(unsigned int i = 0; i< genParticles.size(); i++){
//    reco::GenParticle gen = genParticles[i];
/*
	  float mindr = 1000.;
	  //int k1 = gParticleMotherIndex[i]; //mother index
	  for(unsigned int p = 0; p < genParticles.size(); p++)
	  {	
	  	int k2 = gParticleMotherIndex[p];

		float dr = deltaR(gParticleEta[i], gParticlePhi[i], gParticleEta[p], gParticlePhi[p]);
  		//std::cout<<"Finished  dr"<<std::endl;
		if(p!=i && k1==k2){
  		//std::cout<<"Finished  if id"<<std::endl;
			if(dr <= mindr){
  		//std::cout<<"Finished  if dr"<<std::endl;
			
			mindr = dr;

			gParticleSiblingId[i] = genParticles[p].pdgId();
			gParticleSiblingIndex[i] = p;

		        gParticleSiblingE[i] = genParticles[p].energy();
		        gParticleSiblingPt[i] = genParticles[p].pt();
		        gParticleSiblingPx[i] = genParticles[p].px();
		        gParticleSiblingPy[i] = genParticles[p].py();
		        gParticleSiblingPz[i] = genParticles[p].pz();
		        gParticleSiblingEta[i] = genParticles[p].eta();
		        gParticleSiblingPhi[i] = genParticles[p].phi();
        		gParticleSiblingDR[i] = deltaR(genParticles[p].eta(), genParticles[p].phi(), genParticles[i].eta(), genParticles[i].phi());

			gParticleSiblingId[p] = genParticles[i].pdgId();
			gParticleSiblingIndex[p] = i;
		
	        	gParticleSiblingE[p] = genParticles[i].energy();
	        	gParticleSiblingPt[p] = genParticles[i].pt();
	        	gParticleSiblingPx[p] = genParticles[i].px();
	        	gParticleSiblingPy[p] = genParticles[i].py();
	        	gParticleSiblingPz[p] = genParticles[i].pz();
	        	gParticleSiblingEta[p] = genParticles[i].eta();
	        	gParticleSiblingPhi[p] = genParticles[i].phi();
        		gParticleSiblingDR[p] = deltaR(genParticles[i].eta(), genParticles[i].phi(), genParticles[p].eta(), genParticles[p].phi());

			}//mindr
			//if(gParticleId[i]==22) std::cout<<" Particle GEN Id " << gen.pdgId() << " MotherId " << gParticleMotherId[i]  << " test id " << genParticles[p].pdgId() <<std::endl;
			//if(gParticleId[i]==22) std::cout<<" mother index " << gParticleMotherIndex[i] << " test mother index " << gParticleMotherIndex[p] << " test num dau " << genParticles[p].numberOfDaughters() <<std::endl;
		}// p!=i k1=k2

	  }//finish sibling
*/
//  }
//  //std::cout<<"Finished  fill gen particles sibling"<<std::endl;
//    
//   return true;
//};
//
//  //std::cout<<"Finished  fill gen particles grand mother"<<std::endl;
//
//bool phase2L1EcalTimingAnalyzer::fillGenParticleTPBranches(std::vector<reco::GenParticle> genParticles){
//  for(unsigned int i = 0; i< genParticles.size(); i++){
//    reco::GenParticle gen = genParticles[i];
//
	//corr eta phi to origin
	//vector<float> etaphi = EtaPhi_Corr_EB(genVertexX, genVertexY, genVertexZ, gen);
	vector<float> etaphi = EtaPhi_Corr_EB(gParticle_prod_vtx_x[i], gParticle_prod_vtx_y[i], gParticle_prod_vtx_z[i], gen);
	float eta = etaphi[0];
	float phi = etaphi[1];
	float tof = etaphi[2];
	float tvirtual = etaphi[3];

	g_tof[i] = tof;
	g_tvirtual[i] = tvirtual;

  }//loop of gen

  for(unsigned int i = 0; i< genParticles.size(); i++){
    reco::GenParticle gen = genParticles[i];

    //gammas and pi0s in ecal barrel
    if(( (abs(gen.pdgId())==22 && abs(gParticleMotherId[i])==111) || abs(gen.pdgId())==111)  && abs(gen.eta()) < 1.5 ){
	//conversion between ieta, iphi and detID
	int iEta = 1;
	int iPhi = 1;
	float unit = 2*TMath::Pi()/360;

	//corr eta phi to origin
	//vector<float> etaphi = EtaPhi_Corr_EB(genVertexX, genVertexY, genVertexZ, gen);
	vector<float> etaphi = EtaPhi_Corr_EB(gParticle_prod_vtx_x[i], gParticle_prod_vtx_y[i], gParticle_prod_vtx_z[i], gen);
	float eta = etaphi[0];
	float phi = etaphi[1];
	float tof = etaphi[2];
	float tvirtual = etaphi[3];

	g_tof[i] = tof;
	g_tvirtual[i] = tvirtual;

	iEta = eta_to_iEta(eta);
	iPhi = eta_to_iEta(phi);

	int id = detID_from_iEtaiPhi(iEta, iPhi, true, false);

	//std::cout<<" Particle GEN Eta " << gen.eta() << " Phi " << gen.phi() << " ieta " << gen.eta()/unit << " iphi " << gen.phi()/unit << " iEta " << iEta << " iPhi " << iPhi << " id " << id <<std::endl;
	//std::cout<<" Particle CORR Eta " << eta << " Phi " << phi << " ieta " << eta/unit << " iphi " << phi/unit << " iEta " << iEta << " iPhi " << iPhi << " id " << id <<std::endl;
	
	//super cluster gen energy
	int cnt_sc_02 = 0;
	int cnt_sc_01 = 0;
	float max_sc_02 = 0.;
	float min_sc_02 = 666.;
	float max_sc_01 = 0.;
	float min_sc_01 = 666.;
	int cnt_sm_02 = 0;
	int cnt_sm_01 = 0;
	float max_sm_02 = 0.;
	float min_sm_02 = 666.;
	float max_sm_01 = 0.;
	float min_sm_01 = 666.;
	int cnt_sm2_01 = 0;
	float max_sm2_01 = 0.;
	float min_sm2_01 = 666.;
	int cnt_sm5_01 = 0;
	float max_sm5_01 = 0.;
	float min_sm5_01 = 666.;
	int cnt_sm10_01 = 0;
	float max_sm10_01 = 0.;
	float min_sm10_01 = 666.;
	float t_max_01 = 0.;
	float t_min_01 = 0.;
	float e_max_01 = 0.;
	float e_min_01 = 6000.;
	float t1 = 0.;
	float t2 = 0.;
	float Egen_sc_02 = gen.energy();
	float Egen_sc_01 = gen.energy();
	float dr = 666.;
	float e_mindr = 0.;
	float t_mindr = 0.;

	for(unsigned int q = 0; q < genParticles.size(); q++)
	{
    		reco::GenParticle neighbor = genParticles[q];
		if( abs(neighbor.eta()) >= 1.5) continue;
        	float distance = deltaR(gen.eta(), gen.phi(), genParticles[q].eta(), genParticles[q].phi());

		if(gen.pdgId()==neighbor.pdgId() && distance<0.2 )
		{
			if(q!=i) Egen_sc_02 += neighbor.energy();

			if(g_tof[q] == -666) t2 = 0.;	
			else t2 = g_tof[q];	

			if(t2 >= max_sc_02)
				max_sc_02 = t2;
			if(t2 <= min_sc_02 && t2>0)
				min_sc_02 = t2;

			cnt_sc_02++;

			//if(gParticleMotherIndex[q] == gParticleMotherIndex[i])
			if(neighbor.energy()>1)
			{
	
				if(t2 >= max_sm_02)
					max_sm_02 = t2;
				if(t2 <= min_sm_02 && t2>0)
					min_sm_02 = t2;
	
				cnt_sm_02++;

			}


			if(distance<0.1)
			{
				if(q!=i) Egen_sc_01 += neighbor.energy();

				if(g_tof[q] == -666) t1 = 0.;	
				else t1 = g_tof[q];	

				if(distance<unit*2)
				{
				if(q!=i && gParticle_prod_vtx_x[q]==gParticle_prod_vtx_x[i] && gParticle_prod_vtx_y[q]==gParticle_prod_vtx_y[i] && gParticle_prod_vtx_z[q]==gParticle_prod_vtx_z[i] && distance <= dr)
				{
					dr = distance;
					e_mindr = neighbor.energy();
					t_mindr = t1;
				}
				}

				if(t1 >= max_sc_01)
					max_sc_01 = t1;
				if(t1 <= min_sc_01 && t1>0)
					min_sc_01 = t1;


				if(neighbor.energy() >= e_max_01)
				{
					e_max_01 = neighbor.energy();
					t_max_01 = t1;
				}
				if(neighbor.energy() <= e_min_01 && t1>0)
				{
					e_min_01 = neighbor.energy();
					t_min_01 = t1;
				}

	//if(t1>0 && g_tof[q]>0 && tof>0 && max_sc_01>0 && min_sc_01>0) std::cout<<"index " <<i << " dr 0.1 " << t1 << " " << g_tof[q] <<" "  << g_tof[i] <<" " << " max " << max_sc_01 <<  " min " << min_sc_01 <<std::endl;
	//std::cout<<"index " <<i << " dr 0.1, index q: "<< q << " t1:  " << t1 << " t1 origin (include -666) :  " << g_tof[q] <<" compare to core: "  << g_tof[i] <<" result " << " max: " << max_sc_01 <<  " min: " << min_sc_01 <<std::endl;

				cnt_sc_01++;

				//if(gParticleMotherIndex[q] == gParticleMotherIndex[i])
				if(neighbor.energy()>1)
				{
	
					if(t1 >= max_sm_01)
						max_sm_01 = t1;
					if(t1 <= min_sm_01 && t1>0)
						min_sm_01 = t1;
	
					cnt_sm_01++;

				}

				if(neighbor.energy()>2)
				{
	
					if(t1 >= max_sm2_01)
						max_sm2_01 = t1;
					if(t1 <= min_sm2_01 && t1>0)
						min_sm2_01 = t1;
	
					cnt_sm2_01++;

				}

				if(neighbor.energy()>5)
				{
	
					if(t1 >= max_sm5_01)
						max_sm5_01 = t1;
					if(t1 <= min_sm5_01 && t1>0)
						min_sm5_01 = t1;
	
					cnt_sm5_01++;

				}

				if(neighbor.energy()>10)
				{
	
					if(t1 >= max_sm10_01)
						max_sm10_01 = t1;
					if(t1 <= min_sm10_01 && t1>0)
						min_sm10_01 = t1;
	
					cnt_sm10_01++;

				}

			}//dr 0.1

		}//dr 0.2

	}

	if(dr==666) dr = 0.;
	g_mindr[i] = dr;
	g_e_mindr[i] = e_mindr;
	g_t_mindr[i] = t_mindr;

	if(min_sc_02==666) min_sc_02 =0.;
	if(min_sc_01==666) min_sc_01 =0.;

	g_cnt_sc_02[i] = cnt_sc_02;
	g_cnt_sc_01[i] = cnt_sc_01;

	g_tmax_sc_02[i] = max_sc_02;
	g_tmax_sc_01[i] = max_sc_01;
	g_tmin_sc_02[i] = min_sc_02;
	g_tmin_sc_01[i] = min_sc_01;

	g_dt_sc_02[i] = max_sc_02 - min_sc_02;
	g_dt_sc_01[i] = max_sc_01 - min_sc_01;

	if(min_sm_02==666) min_sm_02 =0.;
	if(min_sm_01==666) min_sm_01 =0.;

	g_cnt_sm_02[i] = cnt_sm_02;
	g_cnt_sm_01[i] = cnt_sm_01;

	g_tmax_sm_02[i] = max_sm_02;
	g_tmax_sm_01[i] = max_sm_01;
	g_tmin_sm_02[i] = min_sm_02;
	g_tmin_sm_01[i] = min_sm_01;

	g_dt_sm_02[i] = max_sm_02 - min_sm_02;
	g_dt_sm_01[i] = max_sm_01 - min_sm_01;

	if(min_sm2_01==666) min_sm2_01 =0.;
	g_cnt_sm2_01[i] = cnt_sm2_01;
	g_tmax_sm2_01[i] = max_sm2_01;
	g_tmin_sm2_01[i] = min_sm2_01;
	g_dt_sm2_01[i] = max_sm2_01 - min_sm2_01;

	if(min_sm5_01==666) min_sm5_01 =0.;
	g_cnt_sm5_01[i] = cnt_sm5_01;
	g_tmax_sm5_01[i] = max_sm5_01;
	g_tmin_sm5_01[i] = min_sm5_01;
	g_dt_sm5_01[i] = max_sm5_01 - min_sm5_01;

	if(min_sm10_01==666) min_sm10_01 =0.;
	g_cnt_sm10_01[i] = cnt_sm10_01;
	g_tmax_sm10_01[i] = max_sm10_01;
	g_tmin_sm10_01[i] = min_sm10_01;
	g_dt_sm10_01[i] = max_sm10_01 - min_sm10_01;

	g_e_max_01[i] = e_max_01;
	g_t_max_01[i] = t_max_01;
	g_dt_max_01[i] = g_tof[i] - t_max_01;

	if(e_min_01==6000) e_min_01 =0.;
	g_e_min_01[i] = e_min_01;
	g_t_min_01[i] = t_min_01;
	g_dt_min_01[i] = g_tof[i] - t_min_01;

	genEsc_02[i] = Egen_sc_02;
	genEsc_01[i] = Egen_sc_01;

	//std::cout<<"Finished  fill gen particles Egen SC"<<std::endl;

	//in cone 0.1/ 0.2 find crystal with Emax
	float Edep = 0.;
	float eb_time_Emax_02 = 0.;
	float eb_sigmat_Emax_02 = 0.;
	float eb_time_Emax_01 = 0.;
	float eb_sigmat_Emax_01 = 0.;
	float Emax_02 = 0.;
	float Emax_01 = 0.;
	int Imax_02 = -666; //crystal index
	int Imax_01 = -666; //crystal index
	int Icore = -666; //crystal index
	float deltaR_eb = 1000.;

	for(int k=0; k<nCrystals; k++){
		//float eb_eta = iEta_to_eta(eb_ieta[k]);
		//float eb_phi = iEta_to_eta(eb_iphi[k]);
		float eb_eta = eb_cell_Eta[k];
		float eb_phi = eb_cell_Phi[k];
		float deltaR_eb_gen = deltaR(eb_eta, eb_phi, eta, phi); // eta, phi are corrected eta, phi wrt origin
		//if(eb_Et[k]!=-666 ) Edep = Et_to_E(eb_Et[k], eb_ieta[k]);
		if(eb_Et[k]!=-666 ) Edep = eb_Edep[k];
		else Edep = 0.;

		//find highest energy deposit within cone 0.2 or 0.1
		if(deltaR_eb_gen < 0.2 ){
			if(Edep >= Emax_02){
				//Edep_sc_02 += Edep;
				Emax_02 = Edep;
				Imax_02 = k;
				eb_time_Emax_02 = eb_time[k];
				eb_sigmat_Emax_02 = eb_sigmat[k];
			}
			if(deltaR_eb_gen < 0.1 ){
				if(Edep >= Emax_01){
					//Edep_sc_01 += Edep;
					Emax_01 = Edep;
					Imax_01 = k;
					eb_time_Emax_01 = eb_time[k];
					eb_sigmat_Emax_01 = eb_sigmat[k];
				}
				if(deltaR_eb_gen <= deltaR_eb){
					deltaR_eb = deltaR_eb_gen;
					Icore = k;
				}
			}
		}
	}//loop of crystals

	g_eb_time_Emax_02[i] = eb_time_Emax_02;
	g_eb_sigmat_Emax_02[i] = eb_sigmat_Emax_02;
	gEmax_02[i] = Emax_02;
	gImax_02[i] = Imax_02;

	g_eb_time_Emax_01[i] = eb_time_Emax_01;
	g_eb_sigmat_Emax_01[i] = eb_sigmat_Emax_01;
	gEmax_01[i] = Emax_01;
	gImax_01[i] = Imax_01;

	gIcore[i] = Icore;

	if(eb_time[Icore]>0){
	g_eb_time[i] = eb_time[Icore];
	}else{
	g_eb_time[i] = 0.;
	}

	if(eb_time[Imax_01]>0){
	g_eb_time_Emax_01[i] = eb_time[Imax_01];
	}else{
	g_eb_time_Emax_01[i] = 0.;
	}

	if(eb_time[Imax_02]>0){
	g_eb_time_Emax_02[i] = eb_time[Imax_02];
	}else{
	g_eb_time_Emax_02[i] = 0.;
	}

	//t = tof + eb_t + t_bs - tvirtual
	float mimic_gen_time = tof + g_eb_time[i] + genVertexT -tvirtual;
	float mimic_gen_time_max_02 = tof + g_eb_time_Emax_02[i] + genVertexT -tvirtual;
	float mimic_gen_time_max_01 = tof + g_eb_time_Emax_01[i] + genVertexT -tvirtual;
        
	gen_time[i] = mimic_gen_time;
	gen_time_max_02[i] = mimic_gen_time_max_02;
	gen_time_max_01[i] = mimic_gen_time_max_01;

  //std::cout<<"Finished  fill gen Emax and Imax"<<std::endl;

	int seed_ieta_02 = eb_ieta[Imax_02];
	int seed_iphi_02 = eb_iphi[Imax_02];

	int seed_ieta_01 = eb_ieta[Imax_01];
	int seed_iphi_01 = eb_iphi[Imax_01];

	float seed_eta_02 = eb_cell_Eta[Imax_02];
	float seed_phi_02 = eb_cell_Phi[Imax_02];

	float seed_eta_01 = eb_cell_Eta[Imax_01];
	float seed_phi_01 = eb_cell_Phi[Imax_01];

	float Edep9x9_02 = 0.;
	float Edep5x5_02 = 0.;
	float Edep3x3_02 = 0.;

	float Edep9x9_01 = 0.;
	float Edep5x5_01 = 0.;
	float Edep3x3_01 = 0.;

	//super cluster Edep
	float Edep_sc_02 = 0.;
	float Edep_sc_01 = 0.;
	
	for(int k=0; k<nCrystals; k++){
		int ieta = eb_ieta[k];
		int iphi = eb_iphi[k];
		
		if( eb_Edep[k]!=-666 && ieta >= seed_ieta_02-4 && ieta <= seed_ieta_02+4 && iphi >= seed_iphi_02-4 && iphi <= seed_iphi_02+4){
			Edep9x9_02 += eb_Edep[k];
			if( eb_Edep[k]!=-666 && ieta >= seed_ieta_02-2 && ieta <= seed_ieta_02+2 && iphi >= seed_iphi_02-2 && iphi <= seed_iphi_02+2){
				Edep5x5_02 += eb_Edep[k];
				if( eb_Edep[k]!=-666 && ieta >= seed_ieta_02-1 && ieta <= seed_ieta_02+1 && iphi >= seed_iphi_02-1 && iphi <= seed_iphi_02+1){
					Edep3x3_02 += eb_Edep[k];
				}
			}
		}//cone 0.2

		if( eb_Edep[k]!=-666 && ieta >= seed_ieta_01-4 && ieta <= seed_ieta_01+4 && iphi >= seed_iphi_01-4 && iphi <= seed_iphi_01+4){
			Edep9x9_01 += eb_Edep[k];
			if( eb_Edep[k]!=-666 && ieta >= seed_ieta_01-2 && ieta <= seed_ieta_01+2 && iphi >= seed_iphi_01-2 && iphi <= seed_iphi_01+2){
				Edep5x5_01 += eb_Edep[k];
				if( eb_Edep[k]!=-666 && ieta >= seed_ieta_01-1 && ieta <= seed_ieta_01+1 && iphi >= seed_iphi_01-1 && iphi <= seed_iphi_01+1){
					Edep3x3_01 += eb_Edep[k];
				}
			}
		}//cone 0.1

		float eb_eta = eb_cell_Eta[k];
		float eb_phi = eb_cell_Phi[k];

		float deltaR_eb_emax_02 = deltaR(eb_eta, eb_phi, seed_eta_02, seed_phi_02); 
		if( eb_Edep[k]!=-666 && deltaR_eb_emax_02<0.2){
			Edep_sc_02 += eb_Edep[k];
		}//cone 0.2

		float deltaR_eb_emax_01 = deltaR(eb_eta, eb_phi, seed_eta_01, seed_phi_01); 
		if( eb_Edep[k]!=-666 && deltaR_eb_emax_01<0.1){
			Edep_sc_01 += eb_Edep[k];
		}//cone 0.2

	}//loop of crystals

	gE9x9_02[i] = Edep9x9_02;
	gE5x5_02[i] = Edep5x5_02;
	gE3x3_02[i] = Edep3x3_02;

	gE9x9_01[i] = Edep9x9_02;
	gE5x5_01[i] = Edep5x5_02;
	gE3x3_01[i] = Edep3x3_02;

	gEsc_02[i] = Edep_sc_02;
	gEsc_01[i] = Edep_sc_01;

  //std::cout<<"Finished  fill gen particles clusters"<<std::endl;
	//std::cout<<"Eta " << gen.eta() << " Phi " << gen.phi() <<  " iEta " << iEta << " iPhi " << iPhi << " E " << gen.energy()  << " Pt " << gen.pt() << " Edep 9x9 " << Edep9x9_01 << " Edep 5x5 " << Edep5x5_01 << " Edep3x3 " << Edep3x3_01  <<std::endl;
	//std::cout<<" Edep SC 0.1 " << Edep_sc_01  <<std::endl;
	//if(Edep_sc_01-Edep3x3_01>0)  std::cout<<" true ?" << Edep_sc_01-Edep3x3_01  <<std::endl;
	//else  std::cout<<" false ?" << Edep_sc_01-Edep3x3_01  <<std::endl;
	

    }//photons and pi0s 

  }//loop of gen


  return true;
};

bool phase2L1EcalTimingAnalyzer::fillGenak4JetBranches(){

  //ak4 gen jet info
  //nGenak4Jets = genJetHandle->size();
  for(const reco::GenJet &jet : *genak4JetHandle){
	nGenak4Jets ++;

	gak4JetMass[nGenak4Jets-1] = jet.mass();
	gak4JetE[nGenak4Jets-1] = jet.energy();
	gak4JetEt[nGenak4Jets-1] = jet.et();
	gak4JetPt[nGenak4Jets-1] = jet.pt();
	gak4JetPx[nGenak4Jets-1] = jet.px();
	gak4JetPy[nGenak4Jets-1] = jet.py();
	gak4JetPz[nGenak4Jets-1] = jet.pz();
	gak4JetEta[nGenak4Jets-1] = jet.eta();
	gak4JetPhi[nGenak4Jets-1] = jet.phi();

	gak4JetArea[nGenak4Jets-1] = jet.jetArea();

	gak4JetPileupE[nGenak4Jets-1] = jet.pileup();
	gak4JetPileupIdFlag[nGenak4Jets-1] = 0;

	//gak4JetPassIdLoose[nGenak4Jets-1] = passJetID(&jet, 0);
	//gak4JetPassIdTight[nGenak4Jets-1] = passJetID(&jet, 1);

	gak4JetMuEnergy[nGenak4Jets-1] = jet.muonEnergy();
	gak4JetEmEnergy[nGenak4Jets-1] = jet.emEnergy();
	//gak4JetChargedEmEnergy[nGenak4Jets-1] = jet.chargedEmEnergy();
	//gak4JetNeutralEmEnergy[nGenak4Jets-1] = jet.neutralEmEnergy();
	gak4JetHadronEnergy[nGenak4Jets-1] = jet.hadEnergy();
	//gak4JetChargedHadronEnergy[nGenak4Jets-1] = jet.chargedHadronEnergy();
	//gak4JetNeutralHadronEnergy[nGenak4Jets-1] = jet.neutralHadronEnergy();

	//std::cout<<" Gen ak4Jet " << nGenak4Jets-1 << " Em energy " << jet.emEnergy()  << " charged " << jet.chargedEmEnergy() << " neutral " << jet.neutralEmEnergy() <<std::endl;
	//std::cout<<" Gen ak4Jet " << nGenak4Jets-1 << " Em energy " << jet.emEnergy()  << " charged " << jet.chargedEmMultiplicity() << " neutral " << jet.neutralEmMultiplicity() <<std::endl;
  }

  return true;
};

bool phase2L1EcalTimingAnalyzer::fillGenak4JetNoNuBranches(){

  std::vector<reco::GenParticle>  genParticles = phase2L1EcalTimingAnalyzer::GetGenParticles();

  //ak4 nonu gen jet info
  //nGenak4JetNoNus = genJetHandle->size();
  for(const reco::GenJet &jet : *genak4JetNoNuHandle){
	nGenak4JetNoNus ++;

	gak4JetNoNuMass[nGenak4JetNoNus-1] = jet.mass();
	gak4JetNoNuE[nGenak4JetNoNus-1] = jet.energy();
	gak4JetNoNuEt[nGenak4JetNoNus-1] = jet.et();
	gak4JetNoNuPt[nGenak4JetNoNus-1] = jet.pt();
	gak4JetNoNuPx[nGenak4JetNoNus-1] = jet.px();
	gak4JetNoNuPy[nGenak4JetNoNus-1] = jet.py();
	gak4JetNoNuPz[nGenak4JetNoNus-1] = jet.pz();
	gak4JetNoNuEta[nGenak4JetNoNus-1] = jet.eta();
	gak4JetNoNuPhi[nGenak4JetNoNus-1] = jet.phi();

	gak4JetNoNuArea[nGenak4JetNoNus-1] = jet.jetArea();

	gak4JetNoNuPileupE[nGenak4JetNoNus-1] = jet.pileup();
	gak4JetNoNuPileupIdFlag[nGenak4JetNoNus-1] = 0;

	//gak4JetNoNuPassIdLoose[nGenak4JetNoNus-1] = passJetID(&jet, 0);
	//gak4JetNoNuPassIdTight[nGenak4JetNoNus-1] = passJetID(&jet, 1);

	gak4JetNoNuMuEnergy[nGenak4JetNoNus-1] = jet.muonEnergy();
	gak4JetNoNuEmEnergy[nGenak4JetNoNus-1] = jet.emEnergy();
	//gak4JetNoNuChargedEmEnergy[nGenak4JetNoNus-1] = jet.chargedEmEnergy();
	//gak4JetNoNuNeutralEmEnergy[nGenak4JetNoNus-1] = jet.neutralEmEnergy();
	gak4JetNoNuHadronEnergy[nGenak4JetNoNus-1] = jet.hadEnergy();
	//gak4JetNoNuChargedHadronEnergy[nGenak4JetNoNus-1] = jet.chargedHadronEnergy();
	//gak4JetNoNuNeutralHadronEnergy[nGenak4JetNoNus-1] = jet.neutralHadronEnergy();

	float Edep = 0.;
	float Edep_sum = 0.;
	float Edep_sum_temp = 0.;
	float jet_time_temp = 0.;
	float jet_time_temp_60ps = 0.;
	float jet_time_temp_100ps = 0.;
	float jet_time_temp_60ps_func = 0.;
	float jet_time_temp_100ps_func = 0.;
	float Edep_sum_temp0p5 = 0.;
	float jet_time_temp0p5 = 0.;
	float Edep_sum_temp1 = 0.;
	float jet_time_temp1 = 0.;
	float Edep_sum_temp2 = 0.;
	float jet_time_temp2 = 0.;
	float Edep_sum_temp3 = 0.;
	float jet_time_temp3 = 0.;
	float Edep_sum_temp5 = 0.;
	float jet_time_temp5 = 0.;
	float Edep_sum_temp10 = 0.;
	float jet_time_temp10 = 0.;

	//cout << "GenJEt: " << jet.pt() << " " << jet.eta() << " " << jet.phi() << " " << jet.energy() << "\n";

	for(unsigned int k=0; k<61200; k++){
		float eb_eta = eb_cell_Eta[k];
		float eb_phi = eb_cell_Phi[k];
		float deltaR_eb_jet = deltaR(eb_eta, eb_phi, jet.eta(), jet.phi()); 

		if(eb_Et[k]!=-666 ) Edep = eb_Edep[k];
		else Edep = 0.;

		//find cells within cone 0.4
		if(deltaR_eb_jet < 0.4 ){
			Edep_sum += Edep;
			if(gen_time_tp[k]!=-666)
			{

			  //cout << "TP " << Edep << " : " << gen_time_tp[k] << "\n";
			  
				Edep_sum_temp += Edep;
				jet_time_temp += gen_time_tp[k]*Edep;
				jet_time_temp_60ps += gen_time_tp_60ps[k]*Edep;
				jet_time_temp_100ps += gen_time_tp_100ps[k]*Edep;
				jet_time_temp_60ps_func += gen_time_tp_60ps_func[k]*Edep;
				jet_time_temp_100ps_func += gen_time_tp_100ps_func[k]*Edep;
				//if(Edep!=0) std::cout<<"in loop TP sum E" << Edep_sum<<"TP time temp"<<jet_time_temp <<" Edep "<< Edep<< std::endl;
				if(Edep>0.5)
				{
					Edep_sum_temp0p5 += Edep;
					jet_time_temp0p5 += gen_time_tp[k]*Edep;
				if(Edep>1)
				{
					Edep_sum_temp1 += Edep;
					jet_time_temp1 += gen_time_tp[k]*Edep;
				if(Edep>2)
				{
					Edep_sum_temp2 += Edep;
					jet_time_temp2 += gen_time_tp[k]*Edep;
				if(Edep>3)
				{
					Edep_sum_temp3 += Edep;
					jet_time_temp3 += gen_time_tp[k]*Edep;
				if(Edep>5)
				{
					Edep_sum_temp5 += Edep;
					jet_time_temp5 += gen_time_tp[k]*Edep;
				if(Edep>10)
				{
					Edep_sum_temp10 += Edep;
					jet_time_temp10 += gen_time_tp[k]*Edep;
				}//10
				}//5
				}//3
				}//2
				}//1
				}//0.5
			}
		}//0.4
		
	}//loop of TP

	//cout << "jet time: " << jet_time_temp/Edep_sum << " | Esum = " << Edep_sum << "\n";

	float jet_time = jet_time_temp/Edep_sum;
	gak4JetNoNuTime[nGenak4JetNoNus-1] = jet_time;	
	gak4JetNoNuEsum[nGenak4JetNoNus-1] = Edep_sum;	

	float jet_time_t = jet_time_temp/Edep_sum_temp;
	gak4JetNoNuTime_t[nGenak4JetNoNus-1] = jet_time_t;	
	gak4JetNoNuEsum_t[nGenak4JetNoNus-1] = Edep_sum_temp;	

	float jet_time_t_60ps = jet_time_temp_60ps/Edep_sum_temp;
	float jet_time_t_100ps = jet_time_temp_100ps/Edep_sum_temp;
	float jet_time_t_60ps_func = jet_time_temp_60ps_func/Edep_sum_temp;
	float jet_time_t_100ps_func = jet_time_temp_100ps_func/Edep_sum_temp;
	gak4JetNoNuTime_t_60ps[nGenak4JetNoNus-1] = jet_time_t_60ps;	
	gak4JetNoNuTime_t_100ps[nGenak4JetNoNus-1] = jet_time_t_100ps;	
	gak4JetNoNuTime_t_60ps_func[nGenak4JetNoNus-1] = jet_time_t_60ps_func;	
	gak4JetNoNuTime_t_100ps_func[nGenak4JetNoNus-1] = jet_time_t_100ps_func;	

	float jet_time_t0p5 = jet_time_temp0p5/Edep_sum_temp0p5;
	gak4JetNoNuTime_t0p5[nGenak4JetNoNus-1] = jet_time_t0p5;	
	gak4JetNoNuEsum_t0p5[nGenak4JetNoNus-1] = Edep_sum_temp0p5;	

	float jet_time_t1 = jet_time_temp1/Edep_sum_temp1;
	gak4JetNoNuTime_t1[nGenak4JetNoNus-1] = jet_time_t1;	
	gak4JetNoNuEsum_t1[nGenak4JetNoNus-1] = Edep_sum_temp1;	

	float jet_time_t2 = jet_time_temp2/Edep_sum_temp2;
	gak4JetNoNuTime_t2[nGenak4JetNoNus-1] = jet_time_t2;	
	gak4JetNoNuEsum_t2[nGenak4JetNoNus-1] = Edep_sum_temp2;	

	float jet_time_t3 = jet_time_temp3/Edep_sum_temp3;
	gak4JetNoNuTime_t3[nGenak4JetNoNus-1] = jet_time_t3;	
	gak4JetNoNuEsum_t3[nGenak4JetNoNus-1] = Edep_sum_temp3;	

	float jet_time_t5 = jet_time_temp5/Edep_sum_temp5;
	gak4JetNoNuTime_t5[nGenak4JetNoNus-1] = jet_time_t5;	
	gak4JetNoNuEsum_t5[nGenak4JetNoNus-1] = Edep_sum_temp5;	

	float jet_time_t10 = jet_time_temp10/Edep_sum_temp10;
	gak4JetNoNuTime_t10[nGenak4JetNoNus-1] = jet_time_t10;	
	gak4JetNoNuEsum_t10[nGenak4JetNoNus-1] = Edep_sum_temp10;	

	//std::cout<<" TP sum E" << Edep_sum<<"TP time temp"<<jet_time_temp <<"time "<< jet_time<< std::endl;
	//std::cout<<" TP" << eb_time[0]<<eb_ieta[0] << eb_Et[0]<<eb_Edep[0]<<eb_id[0]<< std::endl;
	//std::cout<<" Gen ak4JetNoNu " << nGenak4JetNoNus-1 << " Em energy " << jet.emEnergy()  << " charged " << jet.chargedEmEnergy() << " neutral " << jet.neutralEmEnergy() <<std::endl;
	//std::cout<<" Gen ak4JetNoNu " << nGenak4JetNoNus-1 << " Em energy " << jet.emEnergy()  << " charged " << jet.chargedEmMultiplicity() << " neutral " << jet.neutralEmMultiplicity() <<std::endl;


	//*************************************************************
	// Debugging Code Below
	//*************************************************************
	// if (Edep_sum_temp <= 0 && abs(jet.eta()) < 1.5) {
	//   cout << "\n\n";
	//   cout << "GenJet has no time stamp : " << jet.energy() << " " << jet.pt() << " " << jet.eta() << " " << jet.phi() << "\n";
	//   cout << "ECAL TPs\n";
	  
	//   for(unsigned int k=0; k<61200; k++){
	//     if (eb_Edep[k] == -666) continue;

	//     float eb_eta = eb_cell_Eta[k];
	//     float eb_phi = eb_cell_Phi[k];
	//     float deltaR_eb_jet = deltaR(eb_eta, eb_phi, jet.eta(), jet.phi()); 

	//     //find cells within cone 0.4
	//     if(deltaR_eb_jet < 0.4 ){
	//       cout << "TP " << k << " : " << eb_eta << " " <<  eb_phi << " | " << deltaR_eb_jet << " | " << eb_Edep[k] << "\n";	   		

	//       for(unsigned int i = 0; i< genParticles.size(); i++){
	// 	reco::GenParticle gen = genParticles[i];
		
	// 	//gammas and pi0s in ecal barrel
	// 	if( deltaR(genParticles[i].eta(), genParticles[i].phi(), eb_eta , eb_phi ) < 0.05
	// 	    && ((abs(gen.pdgId())==22 || abs(gen.pdgId())==111) || gen.status() == 1)
	// 	    ) {
	// 	  cout << "nearby particle: " << gen.pdgId() << " | " << deltaR(genParticles[i].eta(), genParticles[i].phi(), eb_eta , eb_phi ) << " | "
	// 	       << gen.pt() << " " << gen.eta() << " " << gen.phi()
	// 	       << "\n";
	// 	}
		
	//       }

	//       cout << "\n";
	//     }
	//   }//loop of TP
	 
	//   cout << "\nPi0's nearby\n";	  
	//   for(unsigned int i = 0; i< genParticles.size(); i++){
	//     reco::GenParticle gen = genParticles[i];

	//     //gammas and pi0s in ecal barrel
	//     if( deltaR(genParticles[i].eta(), genParticles[i].phi(), jet.eta() , jet.phi() ) < 0.6
	// 	&& (abs(gen.pdgId())==22 || abs(gen.pdgId())==111)
	// 	) {
	//       cout << "particle: " << gen.pdgId() << " | " << deltaR(genParticles[i].eta(), genParticles[i].phi(), jet.eta() , jet.phi() ) << " | "
	// 	   << gen.pt() << " " << gen.eta() << " " << gen.phi()
	// 	   << "\n";
	//     }
	    
	//   } //loop over gen particles
	// } //if it's a jet with no timestamp




  } //loop over gen jets

  return true;
};



bool phase2L1EcalTimingAnalyzer::fillCaloJetBranches(){

  std::vector<reco::GenParticle>  genParticles = phase2L1EcalTimingAnalyzer::GetGenParticles();

  //ak4 nonu gen jet info
  //nGenak4JetNoNus = genJetHandle->size();
  for(const l1t::Jet &jet : *caloJetHandle){
	nCaloJets++;

	caloJetE[nCaloJets-1] = jet.energy();
	caloJetPt[nCaloJets-1] = jet.pt();
	caloJetEta[nCaloJets-1] = jet.eta();
	caloJetPhi[nCaloJets-1] = jet.phi();

	float Edep = 0.;
	float Edep_sum = 0.;
	float Edep_sum_temp = 0.;
	float jet_time_temp = 0.;
	float jet_time_temp_60ps = 0.;
	float jet_time_temp_100ps = 0.;
	float jet_time_temp_60ps_func = 0.;
	float jet_time_temp_100ps_func = 0.;
	float Edep_sum_temp0p5 = 0.;
	float jet_time_temp0p5 = 0.;
	float Edep_sum_temp1 = 0.;
	float jet_time_temp1 = 0.;
	float Edep_sum_temp2 = 0.;
	float jet_time_temp2 = 0.;
	float Edep_sum_temp3 = 0.;
	float jet_time_temp3 = 0.;
	float Edep_sum_temp5 = 0.;
	float jet_time_temp5 = 0.;
	float Edep_sum_temp10 = 0.;
	float jet_time_temp10 = 0.;

	for(unsigned int k=0; k<61200; k++){
		float eb_eta = eb_cell_Eta[k];
		float eb_phi = eb_cell_Phi[k];
		float deltaR_eb_jet = deltaR(eb_eta, eb_phi, jet.eta(), jet.phi()); 

		if(eb_Et[k]!=-666 ) Edep = eb_Edep[k];
		else Edep = 0.;

		//find cells within cone 0.4
		if(deltaR_eb_jet < 0.4 ){
			Edep_sum += Edep;
			if(gen_time_tp[k]!=-666)
			{

			  //cout << "TP " << Edep << " : " << gen_time_tp[k] << "\n";
			  
				Edep_sum_temp += Edep;
				jet_time_temp += gen_time_tp[k]*Edep;
				jet_time_temp_60ps += gen_time_tp_60ps[k]*Edep;
				jet_time_temp_100ps += gen_time_tp_100ps[k]*Edep;
				jet_time_temp_60ps_func += gen_time_tp_60ps_func[k]*Edep;
				jet_time_temp_100ps_func += gen_time_tp_100ps_func[k]*Edep;
				//if(Edep!=0) std::cout<<"in loop TP sum E" << Edep_sum<<"TP time temp"<<jet_time_temp <<" Edep "<< Edep<< std::endl;
				if(Edep>0.5)
				{
					Edep_sum_temp0p5 += Edep;
					jet_time_temp0p5 += gen_time_tp[k]*Edep;
				if(Edep>1)
				{
					Edep_sum_temp1 += Edep;
					jet_time_temp1 += gen_time_tp[k]*Edep;
				if(Edep>2)
				{
					Edep_sum_temp2 += Edep;
					jet_time_temp2 += gen_time_tp[k]*Edep;
				if(Edep>3)
				{
					Edep_sum_temp3 += Edep;
					jet_time_temp3 += gen_time_tp[k]*Edep;
				if(Edep>5)
				{
					Edep_sum_temp5 += Edep;
					jet_time_temp5 += gen_time_tp[k]*Edep;
				if(Edep>10)
				{
					Edep_sum_temp10 += Edep;
					jet_time_temp10 += gen_time_tp[k]*Edep;
				}//10
				}//5
				}//3
				}//2
				}//1
				}//0.5
			}
		}//0.4
		
	}//loop of TP

	//cout << "jet time: " << jet_time_temp/Edep_sum << " | Esum = " << Edep_sum << "\n";

	float jet_time = jet_time_temp/Edep_sum;
	caloJetTime[nCaloJets-1] = jet_time;	
	caloJetEsum[nCaloJets-1] = Edep_sum;	

	float jet_time_t = jet_time_temp/Edep_sum_temp;
	caloJetTime_t[nCaloJets-1] = jet_time_t;	
	caloJetEsum_t[nCaloJets-1] = Edep_sum_temp;	

	float jet_time_t_60ps = jet_time_temp_60ps/Edep_sum_temp;
	float jet_time_t_100ps = jet_time_temp_100ps/Edep_sum_temp;
	float jet_time_t_60ps_func = jet_time_temp_60ps_func/Edep_sum_temp;
	float jet_time_t_100ps_func = jet_time_temp_100ps_func/Edep_sum_temp;
	caloJetTime_t_60ps[nCaloJets-1] = jet_time_t_60ps;	
	caloJetTime_t_100ps[nCaloJets-1] = jet_time_t_100ps;	
	caloJetTime_t_60ps_func[nCaloJets-1] = jet_time_t_60ps_func;	
	caloJetTime_t_100ps_func[nCaloJets-1] = jet_time_t_100ps_func;	


	//*************************************************************
	// Debugging Code Below
	//*************************************************************
	if (Edep_sum_temp <= 0 && abs(jet.eta()) < 1.5) {
	  cout << "\n\n";
	  cout << "Jet has no time stamp : " << jet.energy() << " " << jet.pt() << " " << jet.eta() << " " << jet.phi() << "\n";
	  cout << "ECAL TPs\n";
	  
	  for(unsigned int k=0; k<61200; k++){
	    if (eb_Edep[k] == -666) continue;

	    float eb_eta = eb_cell_Eta[k];
	    float eb_phi = eb_cell_Phi[k];
	    float deltaR_eb_jet = deltaR(eb_eta, eb_phi, jet.eta(), jet.phi()); 

	    //find cells within cone 0.4
	    if(deltaR_eb_jet < 0.4 ){
	      cout << "TP " << k << " : " << eb_eta << " " <<  eb_phi << " | " << deltaR_eb_jet << " | " << eb_Edep[k] << "\n";	   		

	      for(unsigned int i = 0; i< genParticles.size(); i++){
		reco::GenParticle gen = genParticles[i];
		
		//gammas and pi0s in ecal barrel
		if( deltaR(genParticles[i].eta(), genParticles[i].phi(), eb_eta , eb_phi ) < 0.05
		    && ((abs(gen.pdgId())==22 || abs(gen.pdgId())==111) || gen.status() == 1)
		    ) {
		  cout << "nearby particle: " << gen.pdgId() << " | " << deltaR(genParticles[i].eta(), genParticles[i].phi(), eb_eta , eb_phi ) << " | "
		       << gen.pt() << " " << gen.eta() << " " << gen.phi()
		       << "\n";
		}
		
	      }

	      cout << "\n";
	    }
	  }//loop of TP
	 
	  cout << "\nPi0's nearby\n";	  
	  for(unsigned int i = 0; i< genParticles.size(); i++){
	    reco::GenParticle gen = genParticles[i];

	    //gammas and pi0s in ecal barrel
	    if( deltaR(genParticles[i].eta(), genParticles[i].phi(), jet.eta() , jet.phi() ) < 0.6
		&& (abs(gen.pdgId())==22 || abs(gen.pdgId())==111)
		) {
	      cout << "particle: " << gen.pdgId() << " | " << deltaR(genParticles[i].eta(), genParticles[i].phi(), jet.eta() , jet.phi() ) << " | "
		   << gen.pt() << " " << gen.eta() << " " << gen.phi()
		   << "\n";
	    }
	    
	  } //loop over gen particles
	} //if it's a jet with no timestamp




  } //loop over calojets

  return true;
};




bool phase2L1EcalTimingAnalyzer::fillGenak8JetBranches(){

  //ak8 gen jet info
  //nGenak8Jets = genJetHandle->size();
  for(const reco::GenJet &jet : *genak8JetHandle){
	nGenak8Jets ++;

	gak8JetMass[nGenak8Jets-1] = jet.mass();
	gak8JetE[nGenak8Jets-1] = jet.energy();
	gak8JetEt[nGenak8Jets-1] = jet.et();
	gak8JetPt[nGenak8Jets-1] = jet.pt();
	gak8JetPx[nGenak8Jets-1] = jet.px();
	gak8JetPy[nGenak8Jets-1] = jet.py();
	gak8JetPz[nGenak8Jets-1] = jet.pz();
	gak8JetEta[nGenak8Jets-1] = jet.eta();
	gak8JetPhi[nGenak8Jets-1] = jet.phi();

	gak8JetArea[nGenak8Jets-1] = jet.jetArea();

	gak8JetPileupE[nGenak8Jets-1] = jet.pileup();
	gak8JetPileupIdFlag[nGenak8Jets-1] = 0;

	//gak8JetPassIdLoose[nGenak8Jets-1] = passJetID(&jet, 0);
	//gak8JetPassIdTight[nGenak8Jets-1] = passJetID(&jet, 1);

	gak8JetMuEnergy[nGenak8Jets-1] = jet.muonEnergy();
	gak8JetEmEnergy[nGenak8Jets-1] = jet.emEnergy();
	//gak8JetChargedEmEnergy[nGenak8Jets-1] = jet.chargedEmEnergy();
	//gak8JetNeutralEmEnergy[nGenak8Jets-1] = jet.neutralEmEnergy();
	gak8JetHadronEnergy[nGenak8Jets-1] = jet.hadEnergy();
	//gak8JetChargedHadronEnergy[nGenak8Jets-1] = jet.chargedHadronEnergy();
	//gak8JetNeutralHadronEnergy[nGenak8Jets-1] = jet.neutralHadronEnergy();

	//std::cout<<" Gen ak8Jet " << nGenak8Jets-1 << " Em energy " << jet.emEnergy()  << " charged " << jet.chargedEmEnergy() << " neutral " << jet.neutralEmEnergy() <<std::endl;
	//std::cout<<" Gen ak8Jet " << nGenak8Jets-1 << " Em energy " << jet.emEnergy()  << " charged " << jet.chargedEmMultiplicity() << " neutral " << jet.neutralEmMultiplicity() <<std::endl;
  }

  return true;
};

bool phase2L1EcalTimingAnalyzer::fillGenak8JetNoNuBranches(){

  //ak8 nonu gen jet info
  //nGenak8JetNoNus = genJetHandle->size();
  for(const reco::GenJet &jet : *genak8JetNoNuHandle){
	nGenak8JetNoNus ++;

	gak8JetNoNuMass[nGenak8JetNoNus-1] = jet.mass();
	gak8JetNoNuE[nGenak8JetNoNus-1] = jet.energy();
	gak8JetNoNuEt[nGenak8JetNoNus-1] = jet.et();
	gak8JetNoNuPt[nGenak8JetNoNus-1] = jet.pt();
	gak8JetNoNuPx[nGenak8JetNoNus-1] = jet.px();
	gak8JetNoNuPy[nGenak8JetNoNus-1] = jet.py();
	gak8JetNoNuPz[nGenak8JetNoNus-1] = jet.pz();
	gak8JetNoNuEta[nGenak8JetNoNus-1] = jet.eta();
	gak8JetNoNuPhi[nGenak8JetNoNus-1] = jet.phi();

	gak8JetNoNuArea[nGenak8JetNoNus-1] = jet.jetArea();

	gak8JetNoNuPileupE[nGenak8JetNoNus-1] = jet.pileup();
	gak8JetNoNuPileupIdFlag[nGenak8JetNoNus-1] = 0;

	//gak8JetNoNuPassIdLoose[nGenak8JetNoNus-1] = passJetID(&jet, 0);
	//gak8JetNoNuPassIdTight[nGenak8JetNoNus-1] = passJetID(&jet, 1);

	gak8JetNoNuMuEnergy[nGenak8JetNoNus-1] = jet.muonEnergy();
	gak8JetNoNuEmEnergy[nGenak8JetNoNus-1] = jet.emEnergy();
	//gak8JetNoNuChargedEmEnergy[nGenak8JetNoNus-1] = jet.chargedEmEnergy();
	//gak8JetNoNuNeutralEmEnergy[nGenak8JetNoNus-1] = jet.neutralEmEnergy();
	gak8JetNoNuHadronEnergy[nGenak8JetNoNus-1] = jet.hadEnergy();
	//gak8JetNoNuChargedHadronEnergy[nGenak8JetNoNus-1] = jet.chargedHadronEnergy();
	//gak8JetNoNuNeutralHadronEnergy[nGenak8JetNoNus-1] = jet.neutralHadronEnergy();

	//std::cout<<" Gen ak8JetNoNu " << nGenak8JetNoNus-1 << " Em energy " << jet.emEnergy()  << " charged " << jet.chargedEmEnergy() << " neutral " << jet.neutralEmEnergy() <<std::endl;
	//std::cout<<" Gen ak8JetNoNu " << nGenak8JetNoNus-1 << " Em energy " << jet.emEnergy()  << " charged " << jet.chargedEmMultiplicity() << " neutral " << jet.neutralEmMultiplicity() <<std::endl;
  }

  return true;
};

// ------------corr eta phi  ------------
vector<float> phase2L1EcalTimingAnalyzer::EtaPhi_Corr_EB(float X, float Y, float Z, reco::GenParticle gen){
	//X, Y, Z are for PV

	//gen
	float E = gen.energy();
	float Pt = gen.pt();
	float Px = gen.px();
	float Py = gen.py();
	float Pz = gen.pz();
	float Phi = gen.phi();
	float Eta = gen.eta();

	float x_ecal = 0.;
	float y_ecal = 0.;
	float z_ecal = 0.;

	float radius = sqrt( pow(X,2) + pow(Y,2) ); 
	float radius_ecal = 129.0; //cm

	float t_ecal = (1/30.)*(radius_ecal-radius)/(Pt/E);
	
	x_ecal = X + 30.*(Px/E)*t_ecal;
	y_ecal = Y + 30.*(Py/E)*t_ecal;
	z_ecal = Z + 30.*(Pz/E)*t_ecal;

	float t_virtual = (1/30.)*sqrt(pow(radius_ecal,2)+pow(z_ecal,2));

	//corrections of phi and eta wrt origin
	float phi = atan((y_ecal-0.)/(x_ecal-X-0.));
	if(x_ecal < 0.0) phi = TMath::Pi() + phi;

	phi = dPhi(phi, 0.);

	float theta = atan(sqrt(pow(x_ecal-0.,2)+pow(y_ecal-0.,2))/abs(z_ecal-0.));	
	float eta = -1.0*TMath::Sign(1.0,z_ecal-0.)*log(tan(theta/2));

	vector<float> etaphi;
	etaphi.push_back(eta);
	etaphi.push_back(phi);
	etaphi.push_back(t_ecal);
	etaphi.push_back(t_virtual);
	return etaphi; 
};


// ------------corr eta phi  ------------
vector<float> phase2L1EcalTimingAnalyzer::CalculateTime_EB(float X, float Y, float Z, float hitX, float hitY, float hitZ) {
	//X, Y, Z are for PV

	//gen
	float dist = sqrt( pow(hitX-X,2) + pow(hitY-Y,2) + pow(hitZ-Z,2) ); 
	float t_ecal = (1/30.)*(dist)/1.0;
	
	float radius_ecal = sqrt( pow(hitX,2) + pow(hitY,2) + pow(hitZ,2));
	float t_virtual = (1/30.)*radius_ecal / 1.0;

	vector<float> timeVector;
	timeVector.push_back(t_ecal);
	timeVector.push_back(t_virtual);
	return timeVector; 
};

//
// member functions
//

// ------------ method called for each event  ------------
void
phase2L1EcalTimingAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
  //load event info
  loadEvent(iEvent);

  //reset branches
  resetBranches();

  //fill branches
  // fill event info
  fillEventInfoBranches(iEvent);
  //std::cout<<"Finished  fill event info"<<std::endl;
  //fill gen info
  fillGenParticleBranches();
  //std::cout<<"Finished  fill gen particles"<<std::endl;
  //ecal barrel crystals
  fillEBCrystalBranches(iEvent, iSetup);
  //fillEBCrystalBranches(iEvent, iSetup, genParticles);
  //std::cout<<"Finished  fill ecal barrel crystals"<<std::endl;
  //fill jet info
  //fillGenak4JetBranches();
  //fillGenak4JetNoNuBranches();
  //fillGenak8JetBranches();
  //fillGenak8JetNoNuBranches();
  fillCaloJetBranches();

  ecalTPTree->Fill();
  std::cout<<"Finished Analyzing"<<std::endl;
}


// ------------ method called once each job just before starting event loop  ------------
void 
phase2L1EcalTimingAnalyzer::beginJob()
{
  setBranches();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
phase2L1EcalTimingAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
phase2L1EcalTimingAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(phase2L1EcalTimingAnalyzer);
