#ifndef RUN_HH
#define RUN_HH

#include "G4UserRunAction.hh"

#include "G4AnalysisManager.hh"
//#include "g4root.hh"

/**
 * @class MyRunAction
 * @brief Gestor de las acciones realizadas durante una corrida de simulación.
 * * Se encarga de la inicialización de G4AnalysisManager, definiendo la estructura 
 * de los archivos de salida ROOT y asegurando que los datos se guarden 
 * correctamente al finalizar el Run.
 */

class MyRunAction : public G4UserRunAction 
{
public: 
	MyRunAction();
	~MyRunAction();
	
	virtual void BeginOfRunAction(const G4Run*);
	virtual void EndOfRunAction(const G4Run*);


};

#endif
