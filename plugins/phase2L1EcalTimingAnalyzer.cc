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
  else{
    std::cout<<"Gen Particles T value "<<*genVertexTHandle<<std::endl;
    //genVertexT = *genVertexTHandle;
  }

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
 ecalTPTree->Branch("eb_id", &eb_id, "eb_id[nCrystals]/I");
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

 ecalTPTree->Branch("gE10x10", &gE10x10, "gE10x10[nGenParticles]/F");
 ecalTPTree->Branch("gE5x5", &gE5x5, "gE5x5[nGenParticles]/F");
 ecalTPTree->Branch("gE3x3", &gE3x3, "gE3x3[nGenParticles]/F");
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
 eb_id[i]       = -666;
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

 gE10x10[i] = -666.;
 gE5x5[i] = -666.;
 gE3x3[i] = -666.;
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

  
  //event info
  runNum   = iEvent.id().run();
  lumiNum  = iEvent.id().luminosityBlock();
  eventNum = iEvent.id().event();
  
  edm::ESHandle<CaloGeometry> caloGeometryHandle;
  iSetup.get<CaloGeometryRecord>().get(caloGeometryHandle);
  const CaloGeometry* caloGeometry_ = caloGeometryHandle.product();  
  ebGeometry = caloGeometry_->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);

  //ecal barrel crystals
  for(auto& tpg : *ecaltpgCollection.product())
    {
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
      if(tpg.encodedEt() > 0) 
      {

	      GlobalVector position;
	      auto cell = ebGeometry->getGeometry(tpg.id());

	      int et = tpg.encodedEt();
	      //float et = tpg.encodedEt()/8.;

	      if(et<0.001) continue;//
	        //float energy = et / sin(position.theta());
	      float eta = cell->getPosition().eta();
	      float phi = cell->getPosition().phi();
	      eb_Et[nCrystals-1] = et;
	      //eb_time[nCrystals-1] = tpg.time();
	      float sigma_t = pow(et,-0.34242268082)*0.121; //similar to eta=0 lumi 300/fb
	      std::random_device rd;
	      std::default_random_engine generator;
              std::normal_distribution<float> distribution(0,sigma_t);
              eb_time[nCrystals-1] = distribution(generator);  // generates number
	      //std::cout<<"dice " << eb_time[nCrystals-1] <<std::endl;
	      //float sigma_t = TMath::pow(et,-0.34242268082)*0.121;
	      eb_id[nCrystals-1] = tpg.id();
	      eb_ieta[nCrystals-1] = tpg.id().ieta();
	      eb_iphi[nCrystals-1] = tpg.id().iphi();
	      eb_ism[nCrystals-1] = tpg.id().ism();
	      eb_ic[nCrystals-1] = tpg.id().ic();
	      eb_cell_Eta[nCrystals-1] = eta;
	      eb_cell_Phi[nCrystals-1] = phi;
		      //l1EcalCrystals->Fill(eta,phi,et);
	//std::cout<<"time " << tpg.time()<<std::endl;
	//std::cout<<"id " << tpg.id()<<std::endl;
	//std::cout<<"id " << tpg.id().ieta() << tpg.id().iphi() << tpg.id().ism() << tpg.id().ic() <<std::endl;
	//std::cout<<"approx Eta " << tpg.id().approxEta() <<std::endl;
	//std::cout<<"cell Eta " << eta << " Phi " << phi <<std::endl;
	//std::cout<<"Et above is not zero " <<std::endl;
	}
    }
  // Get genParticles
  std::vector<reco::GenParticle> genPiZeros;
  std::vector<reco::GenParticle> genParticles;
  
  bool foundGenVertex = false;
  for(unsigned int i = 0; i< genParticleHandle->size(); i++){

    edm::Ptr<reco::GenParticle> ptr(genParticleHandle, i);
    if(
       (abs((ptr)->pdgId()) == 2212 ) //protons
       || (abs((ptr)->pdgId()) >= 1 && abs((ptr)->pdgId()) <= 6 && ( (ptr)->status() < 30 )) //quarks
       || (abs((ptr)->pdgId()) >= 11 && abs((ptr)->pdgId()) <= 16) //leptons
       || (abs((ptr)->pdgId()) == 21 && (ptr)->status() < 30) //gluons
       || (abs((ptr)->pdgId()) >= 22 && abs((ptr)->pdgId()) <= 25 && ( (ptr)->status() < 30)) //gammas, z0, w, higgs
       || (abs((ptr)->pdgId()) >= 32 && abs((ptr)->pdgId()) <= 42) // other gauge and higgs bosons
       || (abs((ptr)->pdgId()) == 111 && abs((ptr)->pdgId()) == 211) // pi0 or pi+ pi-
       || (abs((ptr)->pdgId()) == 311 && abs((ptr)->pdgId()) == 321) // k0 or k+ k-
       || (abs((ptr)->pdgId()) == 130 && abs((ptr)->pdgId()) == 310) // k0l or k0s
       || (abs((ptr)->pdgId()) >= 1000001 && abs((ptr)->pdgId()) <= 1000039) //susy particles
       || (abs((ptr)->pdgId()) == 9000006 || abs((ptr)->pdgId()) == 9000007) //llp
     ){
    genParticles.push_back(*ptr);
    }

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

  //fill gen info
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

    // gen mother
    if(gen.numberOfMothers() > 0){
	const reco::Candidate* firstMotherWithDifferentID = findFirstMotherWithDifferentID(&gen);
        if (firstMotherWithDifferentID)
        {
        	gParticleMotherId[i] = firstMotherWithDifferentID->pdgId();
        }

	//find the mother and keep going up the mother chain if the ID's are the same
	const reco::Candidate* originalMotherWithSameID = findOriginalMotherWithSameID(&gen);
	for(unsigned int j = 0; j < genParticles.size(); j++)
	{
		if(gen.pdgId() == originalMotherWithSameID->pdgId() )
		{
			gParticleMotherIndex[i] = j;
			break;
		}
	}
     }// finish gen mother

    //pi0 and Egammas in ecal barrel
    //if( (abs(gen.pdgId())==111 || abs(gen.pdgId())==22 || abs(gen.pdgId())==11) && abs(gen.eta()) < 1.5 ){
    if(abs(gen.pdgId())==22 && abs(gParticleMotherId[i])==111  && gen.pt()>50 && abs(gen.eta()) < 1.5 ){
	//conversion between ieta, iphi and detID
	int iEta = 1;
	int iPhi = 1;
	float unit = 2*TMath::Pi()/360;

	iEta = eta_to_iEta(gen.eta());
	iPhi = eta_to_iEta(gen.phi());

	int id = detID_from_iEtaiPhi(iEta, iPhi, true, false);
	std::cout<<" Particle Eta " << gen.eta() << " Phi " << gen.phi() << " ieta " << gen.eta()/unit << " iphi " << gen.phi()/unit << " iEta " << iEta << " iPhi " << iPhi << " id " << id <<std::endl;

	float Edep10x10 = 0.;
	float Edep5x5 = 0.;
	float Edep3x3 = 0.;
	float Edep = 0.;
	for(int k=0; k<nCrystals; k++){
		for(int k1=iEta-5;k1<iEta+5;k1++){
			for(int k2=iPhi-5;k2<iPhi+5;k2++){
				int id_neighbbor = detID_from_iEtaiPhi(k1, k2, true, false);
				if(eb_id[k]==id_neighbbor){
					Edep = Et_to_E(eb_Et[k], eb_ieta[k]);
  					std::cout<<"Found EB crystal for particle " << gen.pdgId() << " with id " << id_neighbbor << " iEta " << k1  << " iPhi "<< k2  << " Et " << eb_Et[k] <<" E " << Edep <<std::endl;
					Edep10x10 += Edep;
					if(k1>=iEta-2 && k1<=iEta+2 && k2>=iPhi-2 && k2<=iPhi+2) Edep5x5 += Edep;
					if(k1>=iEta-1 && k1<=iEta+1 && k2>=iPhi-1 && k2<=iPhi+1) Edep3x3 += Edep;
				}
			}
		}
	}
	gE10x10[i] = Edep10x10;
	gE5x5[i] = Edep5x5;
	gE3x3[i] = Edep3x3;

	std::cout<<"Eta " << gen.eta() << " Phi " << gen.phi() <<  " iEta " << iEta << " iPhi " << iPhi << " E " << gen.energy()  << " Pt " << gen.pt() << " Edep 10x10 " << Edep10x10 << " Edep 5x5 " << Edep5x5 << " Edep3x3 " << Edep3x3  <<std::endl;
    } 


  }//loop of gen

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
