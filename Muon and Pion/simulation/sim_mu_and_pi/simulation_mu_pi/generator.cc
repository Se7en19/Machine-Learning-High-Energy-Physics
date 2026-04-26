#include "generator.hh"
#include "G4SystemOfUnits.hh"

MyPrimaryGenerator::MyPrimaryGenerator()
{
	fParticleGun = new G4ParticleGun(1);

	G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
	//G4ParticleDefinition *particle = particleTable->FindParticle("mu+");
	//muon = particleTable->FindParticle("mu+");
   	//pion = particleTable->FindParticle("pi+");
   	
   	muonPlus  = particleTable->FindParticle("mu+");
   	muonMinus = particleTable->FindParticle("mu-");
   	pionPlus  = particleTable->FindParticle("pi+");
   	pionMinus = particleTable->FindParticle("pi-");
   	
	//G4ThreeVector pos(0., 0., -2.*m);
	//G4ThreeVector mom(0., 0., 1.);
	//fParticleGun->SetParticleDefinition(particle);
   	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
	//fParticleGun->SetParticlePosition(pos);
	//fParticleGun->SetParticleMomentumDirection(mom);
	//fParticleGun->SetParticleEnergy(100.*MeV);  
	// En lugar de 100*MeV fijo, hacemos un barrido aleatorio
	G4double energy = (G4UniformRand() * 490. + 10.) * MeV; 
	fParticleGun->SetParticleEnergy(energy);
	//fParticleGun->SetParticleDefinition(particle);
}

MyPrimaryGenerator::~MyPrimaryGenerator()
{
	delete fParticleGun;
}

void MyPrimaryGenerator::GeneratePrimaries(G4Event *anEvent)
{
	G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();

    // ... (Keep your position randomization code here) ...

   	G4double rand = G4UniformRand();

   	if (rand < 0.25) {
   	     fParticleGun->SetParticleDefinition(muonPlus);
   	 } 
   	else if (rand < 0.50) {
       		fParticleGun->SetParticleDefinition(muonMinus);
    	} 
   	else if (rand < 0.75) {
       		fParticleGun->SetParticleDefinition(pionPlus);
    	} 
   	else {
       		fParticleGun->SetParticleDefinition(pionMinus);
    	}

    	fParticleGun->GeneratePrimaryVertex(anEvent);

}
