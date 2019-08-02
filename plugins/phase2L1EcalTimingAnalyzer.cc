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
  genSrc_ ((        cfg.getParameter<edm::InputTag>( "genParticles"))),
  genSrcT_ ((        cfg.getParameter<edm::InputTag>( "genParticles_t0")))
{
  //now do what ever initialization is needed
  usesResource("TFileService");
  genToken_ =     consumes<std::vector<reco::GenParticle> >(genSrc_);
  genTokenT_ =     consumes<float>(genSrcT_);
  
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

  if(!iEvent.getByToken(genToken_,genParticleHandle))
    std::cout<<"No gen Particles Found "<<std::endl;
  else
    std::cout<<"Gen Particles size "<<genParticleHandle->size()<<std::endl;

  if(!iEvent.getByToken(genTokenT_,genVertexTHandle))
    std::cout<<"No gen Particles T Found "<<std::endl;
  else
    std::cout<<"Gen Particles T value "<<*genVertexTHandle<<std::endl;

}

// ------------ enable branches  ------------
void phase2L1EcalTimingAnalyzer::setBranches(){
  enableEventInfoBranches();
  enableEBCrystalBranches();
  enableGenParticleBranches();
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
 ecalTPTree->Branch("eb_time", &eb_time, "eb_time[nCrystals]/F");
 ecalTPTree->Branch("eb_ieta", &eb_ieta, "eb_ieta[nCrystals]/I");
 ecalTPTree->Branch("eb_iphi", &eb_iphi, "eb_iphi[nCrystals]/I");
 ecalTPTree->Branch("eb_ism", &eb_ism, "eb_ism[nCrystals]/I");
 ecalTPTree->Branch("eb_ic", &eb_ic, "eb_ic[nCrystals]/I");
 ecalTPTree->Branch("eb_cell_Eta", &eb_cell_Eta, "eb_cell_Eta[nCrystals]/F");
 ecalTPTree->Branch("eb_cell_Phi", &eb_cell_Phi, "eb_cell_Phi[nCrystals]/F");

};

void phase2L1EcalTimingAnalyzer::enableGenParticleBranches(){
 //gen info
 ecalTPTree->Branch("nGenParticles", &nGenParticles, "nGenParticles/I"); 

 ecalTPTree->Branch("genVertexX", &genVertexX, "genVertexX/F"); 
 ecalTPTree->Branch("genVertexY", &genVertexY, "genVertexY/F"); 
 ecalTPTree->Branch("genVertexZ", &genVertexZ, "genVertexZ/F"); 
 ecalTPTree->Branch("genVertexT", &genVertexT, "genVertexT/F"); 

 ecalTPTree->Branch("gParticleMotherId", &gParticleMotherId, "gParticleMotherId[nGenParticles]/I");
 ecalTPTree->Branch("gParticleMotherIndex", &gParticleMotherIndex, "gParticleMotherIndex[nGenParticles]/I");

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

};

// ------------ reset branches  ------------
void phase2L1EcalTimingAnalyzer::resetBranches(){
  resetEventInfoBranches();
  resetEBCrystalBranches();
  resetGenParticleBranches();
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
 eb_time[i]     = -666.;
 eb_ieta[i]     = -666;
 eb_iphi[i]     = -666;
 eb_ism[i]      = -666;
 eb_ic[i]       = -666;
 eb_cell_Eta[i] = -666.;
 eb_cell_Phi[i] = -666.;
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
 gParticleMotherId[i]     = -666;
 gParticleMotherIndex[i]  = -666;

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
 }

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

  
 
  runNum   = iEvent.id().run();
  lumiNum  = iEvent.id().luminosityBlock();
  eventNum = iEvent.id().event();
  
/*
  //edm::Handle< std::vector<L1CaloCluster> > l1Clusters;
  //iEvent.getByToken( L1ClustersToken_, l1Clusters);
  
  //  edm::Handle< std::vector<L1PFObject> > l1PFChargedCandidates;
  //  iEvent.getByToken( L1PFToken_, l1PFChargedCandidates);
  

  edm::ESHandle<CaloGeometry> caloGeometryHandle;
  iSetup.get<CaloGeometryRecord>().get(caloGeometryHandle);
  const CaloGeometry* caloGeometry_ = caloGeometryHandle.product();  
  ebGeometry = caloGeometry_->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);

  for(auto& tpg : *ecaltpgCollection.product())
    {
 	nCrystals++; 

    //std::cout<<"Et " << tpg.encodedEt()<<std::endl;
    //std::cout<<"l1aSpike " << tpg.l1aSpike()<<std::endl;
    //
    if(tpg.time()!=0) {
	std::cout<<"time " << tpg.time()<<std::endl;
	std::cout<<"id " << tpg.id()<<std::endl;
	//std::cout<<"id iEta " << tpg.id().ieta() << " iPhi " << tpg.id().iphi()<<std::endl;
	//std::cout<<"tower iEta " << tpg.id().tower_ieta() << " iPhi " << tpg.id().tower_iphi()<<std::endl;
	std::cout<<"approx Eta " << tpg.id().approxEta() <<std::endl;
	std::cout<<"time above is not zero " <<std::endl;
     }

      if(tpg.encodedEt() > 0) 
      {

	      GlobalVector position;
	      auto cell = ebGeometry->getGeometry(tpg.id());

	      float et = tpg.encodedEt()/8.;

	      if(et<0.001) continue;//
	        //float energy = et / sin(position.theta());
	      float eta = cell->getPosition().eta();
	      float phi = cell->getPosition().phi();
	      eb_Et[nCrystals-1] = et;
	      eb_time[nCrystals-1] = tpg.time();
	      eb_ieta[nCrystals-1] = tpg.id().ieta();
	      eb_iphi[nCrystals-1] = tpg.id().iphi();
	      eb_ism[nCrystals-1] = tpg.id().ism();
	      eb_ic[nCrystals-1] = tpg.id().ic();
	      eb_cell_Eta[nCrystals-1] = eta;
	      eb_cell_Phi[nCrystals-1] = phi;
		      //l1EcalCrystals->Fill(eta,phi,et);
	//std::cout<<"time " << tpg.time()<<std::endl;
	//std::cout<<"id " << tpg.id()<<std::endl;
	//std::cout<<"approx Eta " << tpg.id().approxEta() <<std::endl;
	//std::cout<<"cell Eta " << eta << " Phi " << phi <<std::endl;
	//std::cout<<"Et above is not zero " <<std::endl;
		      }
    }


  // Get genParticles
  std::vector<reco::GenParticle> genPiZeros;
  std::vector<reco::GenParticle> genParticles;

  nGenParticles = genParticleHandle->size();
  for(unsigned int i = 0; i< genParticleHandle->size(); i++){
    edm::Ptr<reco::GenParticle> ptr(genParticleHandle, i);
    genParticles.push_back(*ptr);
 
    //gParticleMotherId[i] = ptr->mother(0)->pdgId();

    gParticleId[i] = ptr->pdgId();
    if(i==0)std::cout<<"Gen Particles Id "<<gParticleId[i]<<std::endl;
    gParticleStatus[i] = ptr->status();

    gParticleE[i] = ptr->energy();
    gParticlePt[i] = ptr->pt();
    gParticlePx[i] = ptr->px();
    gParticlePy[i] = ptr->py();
    gParticlePz[i] = ptr->pz();
    gParticleEta[i] = ptr->eta();
    gParticlePhi[i] = ptr->phi();

    gParticle_decay_vtx_x[i] = ptr->vx();
    gParticle_decay_vtx_y[i] = ptr->vy();
    gParticle_decay_vtx_z[i] = ptr->vz();

    if(abs(ptr->pdgId())==111 && abs(ptr->eta()<1.5)){
      genPiZeros.push_back(*ptr);
      //std::cout<<"Found PiZero PDGID 111 pt: "<<ptr->pt()<<" eta: "<<ptr->eta()<<" phi: "<<ptr->phi()<<std::endl;
      //std::cout<<"Found PiZero PDGID 111 status: "<<ptr->status()<<" vx: "<<ptr->vx()<<" vy: "<<ptr->vy()<<" vz: "<<ptr->vz()<<std::endl;
      //std::cout<<"Found PiZero PDGID 111 daughter num: "<<ptr->numberOfDaughters()<<" mother num: "<<ptr->numberOfMothers()<<std::endl;
      //const reco::Candidate *dau0 = ptr->mother(0);
      int num = ptr->numberOfMothers();

      for(int ntr = 0; ntr < num; ntr++){
	const reco::Candidate *mo = ptr->mother(ntr);
      int num1 = mo->numberOfMothers();
	if(abs(mo->pdgId())==2212)
      		std::cout<<"Found PiZero from proton mo"<<std::endl;
      	for(int ntr1 = 0; ntr1 < num1; ntr1++){
		
	const reco::Candidate *mo1 = mo->mother(ntr1);
      int num2 = mo1->numberOfMothers();
	if(abs(mo1->pdgId())==2212)
      		std::cout<<"Found PiZero from proton mo1"<<std::endl;
      	for(int ntr2 = 0; ntr2 < num2; ntr2++){
		
	const reco::Candidate *mo2 = mo1->mother(ntr2);
      int num3 = mo2->numberOfMothers();
	if(abs(mo2->pdgId())==2212)
      		std::cout<<"Found PiZero from proton mo2"<<std::endl;
      	for(int ntr3 = 0; ntr3 < num3; ntr3++){
		
	const reco::Candidate *mo3 = mo2->mother(ntr3);
      int num4 = mo3->numberOfMothers();
	if(abs(mo3->pdgId())==2212)
      		std::cout<<"Found PiZero from proton mo3"<<std::endl;
      	for(int ntr4 = 0; ntr4 < num4; ntr4++){
		
	const reco::Candidate *mo4 = mo3->mother(ntr4);
      int num5 = mo4->numberOfMothers();
	if(abs(mo4->pdgId())==2212)
      		std::cout<<"Found PiZero from proton mo4"<<std::endl;
      	for(int ntr5 = 0; ntr5 < num5; ntr5++){
		
	const reco::Candidate *mo5 = mo4->mother(ntr4);
      int num6 = mo5->numberOfMothers();
	if(abs(mo5->pdgId())==2212)
      		std::cout<<"Found PiZero from proton mo5"<<std::endl;
	}
	}
	}
	}
	}
		
      }

      if(ptr->mother(0) )std::cout<<"Found PiZero PDGID 111 mother id: "<<ptr->mother(0)->pdgId()<<std::endl;
      if(ptr->mother(0)->mother(0) )std::cout<<"Found PiZero PDGID 111 grand mother id: "<<ptr->mother(0)->mother(0)->pdgId()<<std::endl;
      if(ptr->mother(0)->mother(0)->mother(0) )std::cout<<"Found PiZero PDGID 111 great grand mother id: "<<ptr->mother(0)->mother(0)->mother(0)->pdgId()<<std::endl;

      if(ptr->mother(0) )std::cout<<"Found PiZero PDGID 111 mother num: "<<ptr->numberOfMothers()<<std::endl;
      if(ptr->mother(0)->mother(0) )std::cout<<"Found PiZero PDGID 111 grand mother num: "<<ptr->mother(0)->numberOfMothers()<<std::endl;
      if(ptr->mother(0)->mother(0)->mother(0) )std::cout<<"Found PiZero PDGID 111 great grand mother num: "<<ptr->mother(0)->mother(0)->numberOfMothers()<<std::endl;
      //std::cout<<"Found PiZero PDGID 111 flag Photon: "<<ptr->isPhoton()<<std::endl;
      
      if(ptr->numberOfDaughters()==2){
      std::cout<<"Found PiZero PDGID 111 flag daughter 0 id: "<<ptr->daughter(0)->pdgId()<<std::endl;
      std::cout<<"Found PiZero PDGID 111 flag daughter 1 id: "<<ptr->daughter(1)->pdgId()<<std::endl;
      const reco::Candidate *dau0 = ptr->daughter(0);
      const reco::Candidate *dau1 = ptr->daughter(1);
      std::cout<<"Found PiZero dau0: "<<" vx: "<<dau0->vx()<<" vy: "<<dau0->vy()<<" vz: "<<dau0->vz()<<std::endl;
      std::cout<<"Found PiZero dau0: "<<" pt: "<<dau0->pt()<<" eta: "<<dau0->eta()<<" phi: "<<dau0->phi()<<std::endl;
      std::cout<<"Found PiZero dau1: "<<" vx: "<<dau1->vx()<<" vy: "<<dau1->vy()<<" vz: "<<dau1->vz()<<std::endl;
      std::cout<<"Found PiZero dau1: "<<" pt: "<<dau1->pt()<<" eta: "<<dau1->eta()<<" phi: "<<dau1->phi()<<std::endl;
      }

    }//genPiZero

    if(abs(ptr->pdgId())==2212 && abs(ptr->eta()<1.5)){
      //genPiPluss.push_back(*ptr);
      std::cout<<"Found proton PDGID 2212 pt: "<<ptr->pt()<<" eta: "<<ptr->eta()<<" phi: "<<ptr->phi()<<std::endl;
      std::cout<<"Found proton PDGID 2212 status: "<<ptr->status()<<" vx: "<<ptr->vx()<<" vy: "<<ptr->vy()<<" vz: "<<ptr->vz()<<std::endl;
    }

  }
  ecalTPTree->Fill();

*/   
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
