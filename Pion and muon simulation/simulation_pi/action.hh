#ifndef ACTION_HH
#define ACTION_HH

#include "G4VUserActionInitialization.hh"

#include "generator.hh"
#include "run.hh"
/**
 * @class MyActionInitialization
 * @brief Clase puente para registrar las acciones del usuario en el kernel de Geant4.
 * * Instancia y registra el generador de primarias y las acciones de ejecución 
 * (RunAction) para que el RunManager las ejecute durante la simulación.
 */
class MyActionInitialization : public G4VUserActionInitialization

{
public: 
	MyActionInitialization();
	~MyActionInitialization();
	
	virtual void Build() const;
};
#endif
