#include "generator.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

MyPrimaryGenerator::MyPrimaryGenerator()
{
	fParticleGun = new G4ParticleGun(1);

	G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
	fMuon = particleTable->FindParticle("mu+");
	fPion = particleTable->FindParticle("pi+");

	fParticleGun->SetParticleEnergy(100.*GeV);   // default; sobreescrito por barrido_continuo.mac
	fParticleGun->SetParticleDefinition(fMuon);
}

MyPrimaryGenerator::~MyPrimaryGenerator()
{
	delete fParticleGun;
}

void MyPrimaryGenerator::GeneratePrimaries(G4Event *anEvent)
{
	// Cono: fuente puntual en (0,0,−2m) apuntando a punto uniforme en cara del Fe (±35 cm)
	G4double tx = (G4UniformRand() - 0.5) * 70.*cm;
	G4double ty = (G4UniformRand() - 0.5) * 70.*cm;
	G4ThreeVector source(0., 0., -2.*m);
	G4ThreeVector target(tx, ty, 0.);

	// Selección aleatoria 50/50: mu+ o pi+
	G4ParticleDefinition *particle = (G4UniformRand() < 0.5) ? fMuon : fPion;
	fParticleGun->SetParticleDefinition(particle);
	fParticleGun->SetParticlePosition(source);
	fParticleGun->SetParticleMomentumDirection((target - source).unit());
	fParticleGun->GeneratePrimaryVertex(anEvent);
}
