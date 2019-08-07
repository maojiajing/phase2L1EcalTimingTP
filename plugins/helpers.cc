// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
//#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "helpers.h"
#include <iostream>
#include "TMath.h"
#include <vector>
#include <tgmath.h>

//Gets visible 4-momentum of a particle from list of daughters
reco::Candidate::LorentzVector GetVisibleP4(std::vector<const reco::GenParticle*>& daughters)
{
  reco::Candidate::LorentzVector p4_vis(0,0,0,0);
  for(size_t i = 0; i < daughters.size(); ++i){
    if (!isNeutrino(daughters[i]) && daughters[i]->status() == 1){  //list should include intermediate daughters, so check for status = 1
      p4_vis += daughters[i]->p4();
    }
  }
  return p4_vis;
}

bool isNeutrino(const reco::Candidate* daughter)
{
  return (TMath::Abs(daughter->pdgId()) == 12 || TMath::Abs(daughter->pdgId()) == 14 || TMath::Abs(daughter->pdgId()) == 16 || TMath::Abs(daughter->pdgId()) == 18);
}


int GetDecayModePiZero(const reco::GenParticle* tau, TLorentzVector &pizero){
  int decayMode = -1;
  pizero.SetPtEtaPhiE(0,0,0,0);
  std::vector<int> counts(4,0);
  //std::cout << " ->";
  for(size_t j = 0; j < tau->numberOfDaughters(); ++j){  //looping through first level of decay products 
    int pdg = tau->daughter(j)->pdgId();
    if (TMath::Abs(pdg)==11) ++counts[0];
    if (TMath::Abs(pdg)==13) ++counts[1];
    if (TMath::Abs(pdg)==111) ++counts[2];
    if (TMath::Abs(pdg)==211 || TMath::Abs(pdg)==321) ++counts[3];
    if (TMath::Abs(pdg)==111){
      pizero.SetPtEtaPhiE(tau->daughter(j)->pt(),tau->daughter(j)->eta(),tau->daughter(j)->phi(),tau->daughter(j)->pt());
      std::cout<<"pizero pt: "<<tau->daughter(j)->pt()<<std::endl;
    }
    // if (tau->daughter(j)->status() == 1){  //status=1 means no further decay for this daughter
    // std::cout << " " << pdg;
    //} 
    if (tau->daughter(j)->status() == 2){  //status=2 means this daughter decays 
      //std::cout << " (" << pdg;
      GetDaughterDecayMode(tau->daughter(j), counts); //we check the daughter for decay products
      //std::cout << " )";
    }
  }
  if (counts[0] > 0) decayMode = 3;
  if (counts[1] > 0) decayMode = 4;
  if (counts[3] > 0) decayMode = 10*counts[3]+counts[2];
  return decayMode;
}


int GetDecayMode(const reco::GenParticle* tau){
  int decayMode = -1;
  std::vector<int> counts(4,0);
  //std::cout << " ->";
  for(size_t j = 0; j < tau->numberOfDaughters(); ++j){  //looping through first level of decay products 
    int pdg = tau->daughter(j)->pdgId();
    if (TMath::Abs(pdg)==11) ++counts[0];
    if (TMath::Abs(pdg)==13) ++counts[1];
    if (TMath::Abs(pdg)==111) ++counts[2];
    if (TMath::Abs(pdg)==211 || TMath::Abs(pdg)==321) ++counts[3];
    // if (tau->daughter(j)->status() == 1){  //status=1 means no further decay for this daughter
    // std::cout << " " << pdg;
    //} 
    if (tau->daughter(j)->status() == 2){  //status=2 means this daughter decays 
      //std::cout << " (" << pdg;
      GetDaughterDecayMode(tau->daughter(j), counts); //we check the daughter for decay products
      //std::cout << " )";
    }
  }
  if (counts[0] > 0) decayMode = 3;
  if (counts[1] > 0) decayMode = 4;
  if (counts[3] > 0) decayMode = 10*counts[3]+counts[2];
  return decayMode;
}

/*
const reco::GenParticle* findBestGenMatch(const reco::PFTau& tauObj,
					  std::vector<const reco::GenParticle*>& GenPart, double maxDR) {
  const reco::GenParticle* output = NULL;
  double bestDeltaR = maxDR;
  for (size_t i = 0; i < GenPart.size(); ++i) {
    double deltaR = reco::deltaR(tauObj, *GenPart[i]);
    if (deltaR < maxDR) {
      if (deltaR < bestDeltaR) {
	output = GenPart[i];
	bestDeltaR = deltaR;
      }
    }
  }
  return output;
}
*/
void GetDaughterDecayMode(const reco::Candidate* particle, std::vector<int> &counts){
  //std::cout << " ->";
  for(size_t j = 0; j < particle->numberOfDaughters(); ++j){  //looping through first level of decay products 
    int pdg = particle->daughter(j)->pdgId();
    if (TMath::Abs(pdg)==11) ++counts[0];
    if (TMath::Abs(pdg)==13) ++counts[1];
    if (TMath::Abs(pdg)==111) ++counts[2];
    if (TMath::Abs(pdg)==211 || TMath::Abs(pdg)==321) ++counts[3];

    // if (particle->daughter(j)->status() == 1){  //status=1 means no further decay for this daughter
    // std::cout << " " << pdg;
    //} 
    if (particle->daughter(j)->status() == 2){  //status=2 means this daughter decays 
      //std::cout << " (" << pdg;
      GetDaughterDecayMode(particle->daughter(j), counts); //we check the daughter for decay products
      //std::cout << " )";
    }
  }
}

//Creates a vector of all (including intermediate) daughters for a given mother particle
/*
void findDaughters(const reco::GenParticle* mother, std::vector<const reco::GenParticle*>& daughters)
{
  unsigned numDaughters = mother->numberOfDaughters();
  if (numDaughters == 0) std::cout << " none ";
  for (unsigned iDaughter = 0; iDaughter < numDaughters; ++iDaughter ) {
    const reco::GenParticle* daughter = mother->daughterRef(iDaughter).get();
    if (daughter.status() == 1){  //status = 1 is a final state daughter
      daughters.push_back(daughter); 
    }
    if (daughter->status() == 2){  //status = 2 is an intermediate daughter; will decay further
      daughters.push_back(daughter); 
      findDaughters(daughter, daughters);
    }
  }
  }*/


reco::Candidate::LorentzVector getVisMomentum(const std::vector<const reco::GenParticle*>& daughters, int status)
{
  reco::Candidate::LorentzVector p4Vis(0,0,0,0);

  for ( std::vector<const reco::GenParticle*>::const_iterator daughter = daughters.begin();
	daughter != daughters.end(); ++daughter ) {
    //std::cout << "adding daughter: status "<< (*daughter)->status() <<" pdgId = " << (*daughter)->pdgId() << ", Pt = " << (*daughter)->pt() << ","
    //<< " eta = " << (*daughter)->eta() << ", phi = " << (*daughter)->phi()*180./TMath::Pi() << std::endl;
    if ( (status == -1 || (*daughter)->status() == status) && !isNeutrino(*daughter) ) {
      p4Vis += (*daughter)->p4();
    }
  }

  //std::cout << "--> vis. Momentum: Pt = " << p4Vis.pt() << ", eta = " << p4Vis.eta() << ", phi = " << p4Vis.phi() << std::endl;

  return p4Vis;
}

reco::Candidate::LorentzVector getVisMomentum(const reco::GenParticle* genLeg, const reco::GenParticleCollection* genParticles)
{
  std::vector<const reco::GenParticle*> stableDaughters;
  findDaughters(genLeg, stableDaughters, -1);
  std::cout<<"stableDaughters Size: "<<stableDaughters.size();
  reco::Candidate::LorentzVector p4Vis = getVisMomentum(stableDaughters);

  return p4Vis;
}

//adding option for daughter is prompt
void findDaughters(const reco::GenParticle* mother, std::vector<const reco::GenParticle*>& daughters, int status)
{
  unsigned numDaughters = mother->numberOfDaughters();
  for ( unsigned iDaughter = 0; iDaughter < numDaughters; ++iDaughter ) {
    const reco::GenParticle* daughter = mother->daughterRef(iDaughter).get();

    if ( status == -1 || daughter->status() == status ) daughters.push_back(daughter);

    findDaughters(daughter, daughters, status);
  }
}

const reco::Candidate *findFirstMotherWithDifferentID(const reco::Candidate *particle){

  if( particle == 0 ){
    printf("ERROR! null candidate pointer, this should never happen\n");
    return 0;
  }

  // Is this the first parent with a different ID? If yes, return, otherwise
  // go deeper into recursion
  if (particle->numberOfMothers() > 0 && particle->pdgId() != 0) {
    if (particle->pdgId() == particle->mother(0)->pdgId()
	&& particle->mother(0)->status() != 11  // prevent infinite loop for sherpa documentation gluons
	) {
      return findFirstMotherWithDifferentID(particle->mother(0));
    } else {
      return particle->mother(0);
    }
  }

  return 0;
};

const reco::Candidate *findFirstMotherWithDifferentID(const reco::GenParticle *particle){

  if( particle == 0 ){
    printf("ERROR! null candidate pointer, this should never happen\n");
    return 0;
  }

  // Is this the first parent with a different ID? If yes, return, otherwise
  // go deeper into recursion
  if (particle->numberOfMothers() > 0 && particle->pdgId() != 0) {
    if (particle->pdgId() == particle->mother(0)->pdgId()
	&& particle->mother(0)->status() != 11  // prevent infinite loop for sherpa documentation gluons
	) {
      return findFirstMotherWithDifferentID(particle->mother(0));
    } else {
      return particle->mother(0);
    }
  }

  return 0;
};

const reco::Candidate *findOriginalMotherWithSameID(const reco::Candidate *particle){

  if( particle == 0 ){
    printf("ERROR! null candidate pointer, this should never happen\n");
    return 0;
  }

  // Is there another parent with the same ID? If yes, go deeper into recursion
    if (particle->numberOfMothers() > 0 && particle->pdgId() != 0) {
    if (particle->mother(0)->numberOfMothers() == 0 ||
	particle->mother(0)->status() == 11 ||  // prevent infinite loop for sherpa documentation gluons
	(particle->mother(0)->numberOfMothers() > 0 && particle->mother(0)->mother(0)->pdgId() != particle->mother(0)->pdgId())
	) {
      return particle->mother(0);
    } else {
      return findOriginalMotherWithSameID(particle->mother(0));
    }
  }

  return 0;
}

const reco::Candidate *findOriginalMotherWithSameID(const reco::GenParticle *particle){

  if( particle == 0 ){
    printf("ERROR! null candidate pointer, this should never happen\n");
    return 0;
  }

  // Is there another parent with the same ID? If yes, go deeper into recursion
    if (particle->numberOfMothers() > 0 && particle->pdgId() != 0) {
    if (particle->mother(0)->numberOfMothers() == 0 ||
	particle->mother(0)->status() == 11 ||  // prevent infinite loop for sherpa documentation gluons
	(particle->mother(0)->numberOfMothers() > 0 && particle->mother(0)->mother(0)->pdgId() != particle->mother(0)->pdgId())
	) {
      return particle->mother(0);
    } else {
      return findOriginalMotherWithSameID(particle->mother(0));
    }
  }

  return 0;
}

////conversion between DetId <-> ieta/ix/iphi/iy

int detID_from_iEtaiPhi(int iEta_or_iX=1, int iPhi_or_iY=1, bool isEB = true, bool isEEMinus = false)
{
	uint32_t detID = 0;
	int Ecal = 3;
	int EcalBarrel=1;
	int EcalEndcap=2;
	int iz = isEEMinus?-1:1;

	if(isEB) 
	{
		detID = ((Ecal&0xF)<<28)|((EcalBarrel&0x7)<<25);
		detID |= ((iEta_or_iX>0)?(0x10000|(iEta_or_iX<<9)):((-iEta_or_iX)<<9))|(iPhi_or_iY&0x1FF);
	}
	else
	{
		detID = ((Ecal&0xF)<<28)|((EcalEndcap&0x7)<<25);
		
		detID |=(iPhi_or_iY&0x7f)|((iEta_or_iX&0x7f)<<7)|((iz>0)?(0x4000):(0));
	}

	return int(detID);	
};


int iEta_or_iX_from_detID(int detID=1, bool isEB = true)
{
	int iEta_or_iX = 0;
	uint32_t id_ = uint32_t(detID);
	if(isEB)
	{
		int zside = (id_&0x10000)?(1):(-1);
		int ietaAbs = (id_>>9)&0x7F;
		iEta_or_iX =  zside*ietaAbs;
	}
	else
	{
		iEta_or_iX = (id_>>7)&0x7F;
	}
	return iEta_or_iX;
	
};

int iPhi_or_iY_from_detID(int detID=1, bool isEB = true)
{
	int iPhi_or_iY = 0;
	uint32_t id_ = uint32_t(detID);
	if(isEB)
	{
		iPhi_or_iY =  id_&0x1FF;
	}
	else
	{
		iPhi_or_iY = id_&0x7F;
	}
	return iPhi_or_iY;
};

int eta_to_iEta(float eta){
	float unit = 2*TMath::Pi()/360;
	float ieta = eta/unit;
	int iEta = 1;
	if(ieta >= 0){
		if(abs(ieta-int(ieta)) <= abs(ieta-int(ieta)-1)){
			iEta = int(ieta);
		}
		else{
			iEta = int(ieta)+1;
		} 
	}
	else{
		if(abs(ieta-int(ieta)) <= abs(ieta-int(ieta)+1)){
			iEta = int(ieta);
		}
		else{
			iEta = int(ieta)-1;
		} 
	}
	return iEta;
};

float iEta_to_eta(int iEta){
	float unit = 2*TMath::Pi()/360;
	float eta = iEta*unit;
	return eta;
};

float Et_to_E(float Et, int ieta){
	float E = 0.;
	float theta = 0.;
	float eta = iEta_to_eta(ieta);
	theta = 2*TMath::ATan(exp(-eta));
	E = Et/sin(theta);
  	//std::cout<< " iEta " << ieta  << " Et " << Et <<" E " << E << " theta " << theta <<std::endl;
	return E;
};

float dPhi(float phi1, float phi2){
	float dphi = phi1-phi2;
	while (dphi > TMath::Pi())
  	{	
    		dphi -= TMath::TwoPi();
  	}
  	while (dphi <= -TMath::Pi())
  	{
    	dphi += TMath::TwoPi();
  	}
  	return dphi;
};

float deltaR(float eta1, float phi1, float eta2, float phi2){
	float dphi12 = dPhi(phi1,phi2);
	float deta12 = eta1 - eta2;
	return sqrt( dphi12*dphi12 + deta12*deta12);
};
