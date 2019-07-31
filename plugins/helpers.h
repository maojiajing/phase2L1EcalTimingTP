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

#endif
