#ifndef PHYSICS_HH
#define PHYSICS_HH

#include "G4VModularPhysicsList.hh"
#include "G4EmStandardPhysics.hh"
#include "G4OpticalPhysics.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"

/**
 * @class MyActionInitialization
 * @brief Clase puente para registrar las acciones del usuario en el kernel de Geant4.
 * * Instancia y registra el generador de primarias y las acciones de ejecución 
 * (RunAction) para que el RunManager las ejecute durante la simulación.
 */

class MyPhysicsList : public G4VModularPhysicsList
{
public: 
	MyPhysicsList();
	~MyPhysicsList();


};

#endif
