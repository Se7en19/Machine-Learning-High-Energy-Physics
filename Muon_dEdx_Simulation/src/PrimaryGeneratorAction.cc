//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any  representation or  warranty, express or implied,*
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and available *
// * URL above for the full disclaimer and the limitation of liability.*
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your  *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
#include "PrimaryGeneratorAction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <cmath>


namespace B2
{

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  fParticleGun=new G4ParticleGun(1);
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  fMuPlus = particleTable->FindParticle("mu+");
  fMuMinus = particleTable->FindParticle("mu-");

  fParticleGun->SetParticleDefinition(fMuMinus);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
  fParticleGun->SetParticleEnergy(1.0 * GeV);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  // mundo 
  G4double worldZHalfLength = 0;
  G4LogicalVolume* worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
  G4Box* worldBox = nullptr;
  if (worldLV) worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
  if (worldBox)
    worldZHalfLength = worldBox->GetZHalfLength();
  else {
    G4cerr << "World volume of box not found." << G4endl;
    G4cerr << "Perhaps you have changed geometry." << G4endl;
    G4cerr << "The gun will be place in the center." << G4endl;
  }
  G4bool isMuPlus = (G4UniformRand()<0.5);
  G4ParticleDefinition* particle = isMuPlus ? fMuPlus : fMuMinus;
  fParticleGun->SetParticleDefinition(particle);

  G4double logPmin = std::log(fPmin);
  G4double logPmax = std::log(fPmax);
  G4double pMag = std::exp(logPmin + G4UniformRand() * (logPmax - logPmin));


  G4double totalEnergy = std::sqrt(pMag * pMag + fMuonMass * fMuonMass);
  G4double eKin = totalEnergy - fMuonMass;  // MeV

  fParticleGun->SetParticleEnergy(eKin * MeV);

    // ── Dirección: siempre +Z ──────────────────────────────────────────
    // La rigidez p/q = signo(carga) * |p| / |q|
    // μ⁺ (q=+1): p/q > 0  →  banda derecha
    // μ⁻ (q=-1): p/q < 0  →  banda izquierda
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
  fParticleGun->SetParticlePosition(
        G4ThreeVector(0., 0., -worldZHalfLength + 1.*um));

  fParticleGun->GeneratePrimaryVertex(event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B2
