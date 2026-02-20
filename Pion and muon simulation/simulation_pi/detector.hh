#ifndef DETECTOR_HH
#define DETECTOR_HH

#include "G4VSensitiveDetector.hh"

#include "G4AnalysisManager.hh"

/**
 * @class MySensitiveDetector
 * @brief Extensión de G4VSensitiveDetector para el procesamiento de impactos (hits).
 * * Esta clase intercepta cada paso de la partícula dentro de los volúmenes lógicos 
 * designados para calcular la energía depositada (@f$ dE @f$) y la pérdida de 
 * energía lineal (@f$ dE/dx @f$).
 */

class MySensitiveDetector : public G4VSensitiveDetector
{
public: 
	MySensitiveDetector(G4String);
	~MySensitiveDetector();
	
private: 
	virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *);
};

#endif

