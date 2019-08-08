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
 ecalTPTree->Branch("eb_Edep", &eb_Edep, "eb_Edep[nCrystals]/F");
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

 ecalTPTree->Branch("gEmax_02", &gEmax_02, "gEmax_02[nGenParticles]/F");
 ecalTPTree->Branch("gImax_02", &gImax_02, "gImax_02[nGenParticles]/I");
 ecalTPTree->Branch("gE9x9_02", &gE9x9_02, "gE9x9_02[nGenParticles]/F");
 ecalTPTree->Branch("gE5x5_02", &gE5x5_02, "gE5x5_02[nGenParticles]/F");
 ecalTPTree->Branch("gE3x3_02", &gE3x3_02, "gE3x3_02[nGenParticles]/F");

 ecalTPTree->Branch("gEmax_01", &gEmax_01, "gEmax_01[nGenParticles]/F");
 ecalTPTree->Branch("gImax_01", &gImax_01, "gImax_01[nGenParticles]/I");
 ecalTPTree->Branch("gE9x9_01", &gE9x9_01, "gE9x9_01[nGenParticles]/F");
 ecalTPTree->Branch("gE5x5_01", &gE5x5_01, "gE5x5_01[nGenParticles]/F");
 ecalTPTree->Branch("gE3x3_01", &gE3x3_01, "gE3x3_01[nGenParticles]/F");

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
 eb_Edep[i]     = -666.;
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

 gEmax_02[i] = -666.;
 gImax_02[i] = -666;
 gE9x9_02[i] = -666.;
 gE5x5_02[i] = -666.;
 gE3x3_02[i] = -666.;

 gEmax_01[i] = -666.;
 gImax_01[i] = -666;
 gE9x9_01[i] = -666.;
 gE5x5_01[i] = -666.;
 gE3x3_01[i] = -666.;

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

	//corrections of phi and eta wrt origin
	float phi = atan((y_ecal-0.)/(x_ecal-X-0.));
	if(x_ecal < 0.0) phi = TMath::Pi() + phi;

	phi = dPhi(phi, 0.);

	float theta = atan(sqrt(pow(x_ecal-0.,2)+pow(y_ecal-0.,2))/abs(z_ecal-0.));	
	float eta = -1.0*TMath::Sign(1.0,z_ecal-0.)*log(tan(theta/2));

	vector<float> etaphi;
	etaphi.push_back(eta);
	etaphi.push_back(phi);
	return etaphi; 
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

	eb_id[nCrystals-1] = tpg.id();
	eb_ieta[nCrystals-1] = tpg.id().ieta();
	eb_iphi[nCrystals-1] = tpg.id().iphi();
	eb_ism[nCrystals-1] = tpg.id().ism();
	eb_ic[nCrystals-1] = tpg.id().ic();

	GlobalVector position;
	auto cell = ebGeometry->getGeometry(tpg.id());
	float eta = cell->getPosition().eta();
	float phi = cell->getPosition().phi();
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
	      eb_Edep[nCrystals-1] = Et_to_E(et, tpg.id().ieta());
	      //eb_time[nCrystals-1] = tpg.time();
	      float sigma_t = pow(et,-0.34242268082)*0.121; //similar to eta=0 lumi 300/fb
	      std::random_device rd;
	      std::default_random_engine generator;
              std::normal_distribution<float> distribution(0,sigma_t);
              eb_time[nCrystals-1] = distribution(generator);  // generates number
	      //std::cout<<"dice " << eb_time[nCrystals-1] <<std::endl;
	      //float sigma_t = TMath::pow(et,-0.34242268082)*0.121
	//std::cout<<"time " << tpg.time()<<std::endl;
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
       || (abs((ptr)->pdgId()) >= 1 && abs((ptr)->pdgId()) <= 6 )//&& ( (ptr)->status() < 30 )) //quarks
       || (abs((ptr)->pdgId()) >= 11 && abs((ptr)->pdgId()) <= 16) //leptons
       || (abs((ptr)->pdgId()) == 21 && (ptr)->status() < 30) //gluons
       || (abs((ptr)->pdgId()) >= 22 && abs((ptr)->pdgId()) <= 25) //&& ( (ptr)->status() < 30)) //gammas, z0, w, higgs
       || (abs((ptr)->pdgId()) >= 32 && abs((ptr)->pdgId()) <= 42) // other gauge and higgs bosons
       || (abs((ptr)->pdgId()) == 111 || abs((ptr)->pdgId()) == 211) // pi0 or pi+ pi-
       || (abs((ptr)->pdgId()) == 311 || abs((ptr)->pdgId()) == 321) // k0 or k+ k-
       || (abs((ptr)->pdgId()) == 130 || abs((ptr)->pdgId()) == 310) // k0l or k0s
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
		if(genParticles[j].pdgId() == originalMotherWithSameID->pdgId() && genParticles[j].energy() == originalMotherWithSameID->energy() && genParticles[j].eta() ==originalMotherWithSameID->eta() && genParticles[j].phi() == originalMotherWithSameID->phi() && genParticles[j].pt() == originalMotherWithSameID->pt() && genParticles[j].px() == originalMotherWithSameID->px() )
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

    //pi0 and Egammas in ecal barrel
    //if( (abs(gen.pdgId())==111 || abs(gen.pdgId())==22 || abs(gen.pdgId())==11) && abs(gen.eta()) < 1.5 ){
    //if(abs(gen.pdgId())==22 && abs(gParticleMotherId[i])==111  && gen.pt()>50 && abs(gen.eta()) < 1.5 ){
    //
    //gammas in ecal barrel
    if(abs(gen.pdgId())==22 && abs(gParticleMotherId[i])==111  && abs(gen.eta()) < 1.5 ){
	//conversion between ieta, iphi and detID
	int iEta = 1;
	int iPhi = 1;
	float unit = 2*TMath::Pi()/360;

	//corr eta phi to origin
	vector<float> etaphi = EtaPhi_Corr_EB(genVertexX, genVertexY, genVertexZ, gen);
	float eta = etaphi[0];
	float phi = etaphi[1];

	iEta = eta_to_iEta(eta);
	iPhi = eta_to_iEta(phi);

	int id = detID_from_iEtaiPhi(iEta, iPhi, true, false);

	//std::cout<<" Particle GEN Eta " << gen.eta() << " Phi " << gen.phi() << " ieta " << gen.eta()/unit << " iphi " << gen.phi()/unit << " iEta " << iEta << " iPhi " << iPhi << " id " << id <<std::endl;
	//std::cout<<" Particle CORR Eta " << eta << " Phi " << phi << " ieta " << eta/unit << " iphi " << phi/unit << " iEta " << iEta << " iPhi " << iPhi << " id " << id <<std::endl;

	float Edep = 0.;
	float Emax_02 = 0.;
	float Emax_01 = 0.;
	int Imax_02 = -666; //crystal index
	int Imax_01 = -666; //crystal index
	for(int k=0; k<nCrystals; k++){
		//float eb_eta = iEta_to_eta(eb_ieta[k]);
		//float eb_phi = iEta_to_eta(eb_iphi[k]);
		float eb_eta = eb_cell_Eta[k];
		float eb_phi = eb_cell_Phi[k];
		float deltaR_eb_gen = deltaR(eb_eta, eb_phi, eta, phi);
		//if(eb_Et[k]!=-666 ) Edep = Et_to_E(eb_Et[k], eb_ieta[k]);
		if(eb_Et[k]!=-666 ) Edep = eb_Edep[k];
		else Edep = 0.;

		//find highest energy deposit within cone 0.2 or 0.1
		if(deltaR_eb_gen < 0.2 ){
			if(Edep >= Emax_02){
				Emax_02 = Edep;
				Imax_02 = k;
			}
			if(deltaR_eb_gen < 0.1 ){
				if(Edep >= Emax_01){
					Emax_01 = Edep;
					Imax_01 = k;
				}
			}
		}
	}//loop of crystals

	gEmax_02[i] = Emax_02;
	gImax_02[i] = Imax_02;

	gEmax_01[i] = Emax_01;
	gImax_01[i] = Imax_01;

	int seed_ieta_02 = eb_ieta[Imax_02];
	int seed_iphi_02 = eb_iphi[Imax_02];

	int seed_ieta_01 = eb_ieta[Imax_01];
	int seed_iphi_01 = eb_iphi[Imax_01];

	float Edep9x9_02 = 0.;
	float Edep5x5_02 = 0.;
	float Edep3x3_02 = 0.;

	float Edep9x9_01 = 0.;
	float Edep5x5_01 = 0.;
	float Edep3x3_01 = 0.;

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
	}//loop of crystals

	gE9x9_02[i] = Edep9x9_02;
	gE5x5_02[i] = Edep5x5_02;
	gE3x3_02[i] = Edep3x3_02;

	gE9x9_01[i] = Edep9x9_02;
	gE5x5_01[i] = Edep5x5_02;
	gE3x3_01[i] = Edep3x3_02;

/*
	for(int k=0; k<nCrystals; k++){
		for(int k1=-85;k1<85;k1++){
			for(int k2=0;k2<360;k2++){
				int id_neighbbor = detID_from_iEtaiPhi(k1, k2, true, false);
				float eb_eta = iEta_to_eta(k1);
				float eb_phi = iEta_to_eta(k2);
				float deltaR_eb_gen = deltaR(eb_eta, eb_phi, eta, phi);
				Edep = Et_to_E(eb_Et[k], eb_ieta[k]);
				if(deltaR_eb_gen < 0.3 && eb_id[k]==id_neighbbor ){
  					//std::cout<<"Found EB crystal for particle " << gen.pdgId() << " with id " << id_neighbbor << " iEta " << k1  << " iPhi "<< k2  << " Et " << eb_Et[k] <<" E " << Edep <<std::endl;
					Edep9x9 += Edep;
					if(k1>=iEta-2 && k1<=iEta+2 && k2>=iPhi-2 && k2<=iPhi+2) Edep5x5 += Edep;
					if(k1>=iEta-1 && k1<=iEta+1 && k2>=iPhi-1 && k2<=iPhi+1) Edep3x3 += Edep;
				}
			}
		}
	}
*/
	//std::cout<<"Eta " << gen.eta() << " Phi " << gen.phi() <<  " iEta " << iEta << " iPhi " << iPhi << " E " << gen.energy()  << " Pt " << gen.pt() << " Edep 9x9 " << Edep9x9 << " Edep 5x5 " << Edep5x5 << " Edep3x3 " << Edep3x3  <<std::endl;
    }//photons 

  }//loop of gen

  for(unsigned int i = 0; i< genParticles.size(); i++){
	  reco::GenParticle gen = genParticles[i];

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
		
	 }

	  float mindr = 0.;
	  for(unsigned int p = 0; p < genParticles.size(); p++)
	  {	
	  	int k2 = gParticleMotherIndex[p];

		//float dr = deltaR(gParticleEta[i], gParticlePhi[i], gParticleEta[p], gParticlePhi[p]);
		if(p!=i && k1==k2){
			//if(dr > mindr){}
			gParticleSiblingId[i] = genParticles[p].pdgId();
			gParticleSiblingIndex[i] = p;

			//sibling info
			int si = gParticleSiblingIndex[i];
			if(si!=-666){
		    		reco::GenParticle sibling = genParticles[si];
		
		        	gParticleSiblingE[i] = sibling.energy();
		        	gParticleSiblingPt[i] = sibling.pt();
		        	gParticleSiblingPx[i] = sibling.px();
		        	gParticleSiblingPy[i] = sibling.py();
		        	gParticleSiblingPz[i] = sibling.pz();
		        	gParticleSiblingEta[i] = sibling.eta();
		        	gParticleSiblingPhi[i] = sibling.phi();
        			gParticleSiblingDR[i] = deltaR(sibling.eta(), sibling.phi(), genParticles[i].eta(), genParticles[i].phi());
			}
			//if(gParticleId[i]==22) std::cout<<" Particle GEN Id " << gen.pdgId() << " MotherId " << gParticleMotherId[i]  << " test id " << genParticles[p].pdgId() <<std::endl;
			//if(gParticleId[i]==22) std::cout<<" mother index " << gParticleMotherIndex[i] << " test mother index " << gParticleMotherIndex[p] << " test num dau " << genParticles[p].numberOfDaughters() <<std::endl;
		}
	  }
  }

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
