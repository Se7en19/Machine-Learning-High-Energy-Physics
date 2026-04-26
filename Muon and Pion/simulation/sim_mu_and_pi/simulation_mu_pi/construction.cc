#include "construction.hh"
#include "detector.hh"

MyDetectorConstruction::MyDetectorConstruction()
	: logicDetector(nullptr)
{}

MyDetectorConstruction::~MyDetectorConstruction()
{}

G4VPhysicalVolume *MyDetectorConstruction::Construct()
{
	G4NistManager *nist = G4NistManager::Instance();

	G4Material *worldMat = nist->FindOrBuildMaterial("G4_Galactic");
	G4Material *hierro   = nist->FindOrBuildMaterial("G4_Fe");
	G4Material *plastico = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

	// Mundo: semi-lados 6m x 6m x 5m
	G4Box *solidWorld = new G4Box("solidWorld", 1.0*m, 1.0*m, 2.5*m);
	G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
	G4VPhysicalVolume *physWorld = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicWorld, "physWorld", 0, false, 0, true);

	// Absorbedor: cubo de hierro 5x5x5 cm (semi-lados = 0.025 m)
	// Cara delantera en z=0, cara trasera en z=0.05 m
	//G4Box *solidHierro = new G4Box("solidHierro", 0.025*m, 0.025*m, 0.025*m);
	//G4LogicalVolume *logicHierro = new G4LogicalVolume(solidHierro, hierro, "logicHierro");
	//new G4PVPlacement(0, G4ThreeVector(0., 0., 0.025*m), logicHierro, "physHierro", logicWorld, false, 0, true);

	// Detector: plástico centellador BC404, 10x10x0.1 m (semi-lados = 5.0, 5.0, 0.05 m)
	// Separado 1 m del hierro: cara delantera en z=1.05 m, cara trasera en z=1.15 m
	G4Box *solidPlastico = new G4Box("solidPlastico", 5.0*m, 5.0*m, 0.05*m);
	logicDetector = new G4LogicalVolume(solidPlastico, plastico, "logicDetector");
	new G4PVPlacement(0, G4ThreeVector(0., 0., 1.10*m), logicDetector, "physPlastico", logicWorld, false, 0, true);
	
	G4Box *solidBarraX = new G4Box("solidBarraX", 0.5*m, 2.5*cm, 0.5*m);
	G4Box *solidBarraY = new G4Box("solidBarraY", 2.5*cm, 0.5*m, 0.5*cm);
	G4LogicalVolume *logicBarraX = new G4LogicalVolume(solidBarraX, plastico, "logicBarraX");
	G4LogicalVolume *logicBarraY  = new G4LogicalVolume(solidBarraY, plastico, "logicBarraY");
	
	G4double zCapa1 = 1.0*m;
	G4double zCapa2 = zCapa1 + 1.0*cm;
	
	for(int i = 0; i < 20; i++){
	
		G4double posicion = -47.5*cm + (i * 5.0*cm);
	
		new G4PVPlacement(0, G4ThreeVector(0, posicion, zCapa1), logicBarraX, "physBarraX", logicWorld, false, i, true);
	
		new G4PVPlacement(0, G4ThreeVector(posicion, 0, zCapa2), logicBarraY, "physBarraY", logicWorld, false, i + 20, true); 
	
	}
	
	return physWorld;
}

void MyDetectorConstruction::ConstructSDandField()
{
	MySensitiveDetector *sensDet = new MySensitiveDetector("SensitiveDetector");
	SetSensitiveDetector("logicBarraX", sensDet);
	SetSensitiveDetector("logicBarraY", sensDet);

	//if(logicDetector != NULL)
		//logicDetector->SetSensitiveDetector(sensDet);
}
