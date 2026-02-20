#include "physics.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
MyPhysicsList::MyPhysicsList()

{
	RegisterPhysics (new G4EmStandardPhysics());
	RegisterPhysics (new G4OpticalPhysics());
  RegisterPhysics (new G4HadronPhysicsFTFP_BERT());
}

MyPhysicsList::~MyPhysicsList()
{}

