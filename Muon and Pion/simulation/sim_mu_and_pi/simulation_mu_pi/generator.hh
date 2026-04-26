#ifndef GENERATOR_HH
#define GENERATOR_HH

#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"

class MyPrimaryGenerator : public G4VUserPrimaryGeneratorAction

{
public: 
	MyPrimaryGenerator();
	~MyPrimaryGenerator();
	
	virtual void GeneratePrimaries(G4Event*);
	
private: 
	G4ParticleGun *fParticleGun;
	G4ParticleDefinition *muonPlus;
   	G4ParticleDefinition *muonMinus;
   	G4ParticleDefinition *pionPlus;
   	G4ParticleDefinition *pionMinus;
	//G4ParticleDefinition *muon;
	//G4ParticleDefinition *pion;
};
#endif
