#ifndef GENERATOR_HH
#define GENERATOR_HH

#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"

/**
 * @class MyPrimaryGenerator
 * @brief Configuración de la fuente primaria de partículas.
 * * Utiliza un G4ParticleGun para disparar piones o muones cargados (@c pi+ / @c mu+) con una 
 * energía cinética de 100 GeV (esto puede ser modificado) y una dirección de momento fija hacia el detector.
 */

class MyPrimaryGenerator : public G4VUserPrimaryGeneratorAction

{
public: 
	MyPrimaryGenerator();
	~MyPrimaryGenerator();
	
	virtual void GeneratePrimaries(G4Event*);
	
private: 
	G4ParticleGun *fParticleGun;
};
#endif
