#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"

#include "detector.hh"

/**   
 * @class MyDetectorConstruction
 * @brief Clase fundamental para el desarollo del mundo (entorno experimental)
 * * Gestiona la seleccion de materiales usando la base de datos NIST, la creacion
 * de volumenes logicos y fisicos, y la segmentacion del detector en un rejilla.
 * * @note Implementa la configuracion de Detecteres Sensibles (SD) en volumenes
 *  especificos.
 */

class MyDetectorConstruction : public G4VUserDetectorConstruction
{
public:
	MyDetectorConstruction();
	~MyDetectorConstruction();

	virtual G4VPhysicalVolume *Construct();
	
private:
	G4LogicalVolume *logicDetector;
	virtual void ConstructSDandField();
};

#endif
