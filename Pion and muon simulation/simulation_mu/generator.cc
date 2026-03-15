#include "generator.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"

MyPrimaryGenerator::MyPrimaryGenerator()
{
	fParticleGun = new G4ParticleGun(1);
	
	G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName="mu-";
	G4ParticleDefinition *particle = particleTable->FindParticle("mu-");
	
	G4ThreeVector pos(0., 0., 0.);
	G4ThreeVector mom(0., 0., 1.);
	
	fParticleGun->SetParticlePosition(pos);
	fParticleGun->SetParticleMomentumDirection(mom);
	fParticleGun->SetParticleMomentum(100.*GeV);
	fParticleGun->SetParticleDefinition(particle);
}

MyPrimaryGenerator::~MyPrimaryGenerator()
{
	delete fParticleGun;
}

void MyPrimaryGenerator::GeneratePrimaries(G4Event *anEvent)
{	
	fParticleGun->GeneratePrimaryVertex(anEvent);
	
	G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
    
    	G4double pMin = 50. * MeV;
    	G4double pMax = 1000. * MeV;
    	
    	G4double pRandomMinus = pMin + (pMax - pMin) * G4UniformRand();
    	
    	fParticleGun->SetParticleDefinition(particleTable->FindParticle("mu-"));
    	fParticleGun->SetParticleMomentum(pRandomMinus);
    	fParticleGun->GeneratePrimaryVertex(anEvent);
    
    	G4double pRandomPlus = pMin + (pMax - pMin) * G4UniformRand();
    	
    	fParticleGun->SetParticleDefinition(particleTable->FindParticle("mu+"));
    	fParticleGun->SetParticleMomentum(pRandomPlus);
    	fParticleGun->GeneratePrimaryVertex(anEvent);
}
