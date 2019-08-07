/*
 * =====================================================================================
 *
 *       Filename:  Helpers.h
 *
 *    Description:  Common Gen Functions
 *
 *         Author:  M. Cepeda, S. Dasu, E. Friis
 *        Company:  UW Madison
 *
 * =====================================================================================
 */

#ifndef HELPERS_W9QK6HND
#define HELPERS_W9QK6HND
#include "DataFormats/PatCandidates/interface/Jet.h"
#include <TLorentzVector.h>
//#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
//MINIAOD

reco::Candidate::LorentzVector GetVisibleP4(std::vector<const reco::GenParticle*>& daughters);
bool isNeutrino(const reco::Candidate* daughter);
//int GetDecayMode(const reco::GenParticle* tau);
int GetDecayMode(const reco::GenParticle* tau);
int GetDecayModePiZero(const reco::GenParticle* tau, TLorentzVector &piZero);

//const reco::GenParticle* findBestGenMatch(const reco::PFTau& TagTauObj,std::vector<const reco::GenParticle*>& GenPart, double maxDR);
void GetDaughterDecayMode(const reco::Candidate* particle, std::vector<int> &counts);
//void findDaughters(const reco::GenParticle* mother, std::vector<const reco::GenParticle*>& daughters)
reco::Candidate::LorentzVector getVisMomentum(const std::vector<const reco::GenParticle*>&, int = 1);
reco::Candidate::LorentzVector getVisMomentum(const reco::GenParticle* genLeg, const reco::GenParticleCollection* genParticles);
void findDaughters(const reco::GenParticle* mother, std::vector<const reco::GenParticle*>& daughters, int status);

// ------------ help functions  ------------
const reco::Candidate* findFirstMotherWithDifferentID(const reco::Candidate *particle);
const reco::Candidate* findFirstMotherWithDifferentID(const reco::GenParticle *particle);
const reco::Candidate* findOriginalMotherWithSameID(const reco::Candidate *particle);
const reco::Candidate* findOriginalMotherWithSameID(const reco::GenParticle *particle);
//conversion between DetId <-> ieta/ix/iphi/iy
int detID_from_iEtaiPhi(int iEta_or_iX, int iPhi_or_iY, bool isEB, bool isEEMinus);
int iEta_or_iX_from_detID(int detID, bool isEB);
int iPhi_or_iY_from_detID(int detID, bool isEB);
//conversion
int eta_to_iEta(float eta);
float iEta_to_eta(int iEta);
float Et_to_E(float Et, int ieta);
//dR
float dPhi(float phi1, float phi2);
float deltaR(float eta1, float phi1, float eta2, float phi2);
#endif
