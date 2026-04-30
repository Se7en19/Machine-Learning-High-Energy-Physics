#include "construction.hh"
#include "detector.hh"

MyDetectorConstruction::MyDetectorConstruction()
	: logicBarraX(nullptr), logicBarraY(nullptr)
{}

MyDetectorConstruction::~MyDetectorConstruction()
{}

G4VPhysicalVolume *MyDetectorConstruction::Construct()
{
	G4NistManager *nist = G4NistManager::Instance();

	G4Material *worldMat = nist->FindOrBuildMaterial("G4_Galactic");
	G4Material *hierro   = nist->FindOrBuildMaterial("G4_Fe");
	G4Material *plastico = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

	// ── Mundo: semi-lados 2×2×3 m ──────────────────────────────────────────
	G4Box *solidWorld = new G4Box("solidWorld", 2.0*m, 2.0*m, 3.0*m);
	G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
	G4VPhysicalVolume *physWorld = new G4PVPlacement(
		0, G4ThreeVector(0., 0., 0.), logicWorld, "physWorld", 0, false, 0, true);

	// ── Absorbedor: cubo de hierro 70×70×70 cm ─────────────────────────────
	// Cara delantera en z=0, cara trasera en z=70 cm, centro en z=35 cm
	// 70 cm de Fe ≈ 4.2 longitudes de interacción nuclear (λ_I ≈ 16.77 cm)
	G4double ironHalf = 35.0*cm;
	G4Box *solidHierro = new G4Box("solidHierro", ironHalf, ironHalf, ironHalf);
	G4LogicalVolume *logicHierro = new G4LogicalVolume(solidHierro, hierro, "logicHierro");
	new G4PVPlacement(0, G4ThreeVector(0., 0., ironHalf),
	                  logicHierro, "physHierro", logicWorld, false, 0, true);

	// ── Barras centelladora: 1 m largo × 5 cm ancho × 1 cm grosor ──────────
	// Half-lengths: (50 cm, 2.5 cm, 0.5 cm) para barra a lo largo de X
	//               (2.5 cm, 50 cm, 0.5 cm) para barra a lo largo de Y
	G4double barLong   = 50.0*cm;
	G4double barAncho  = 2.5*cm;
	G4double barGrosor = 0.5*cm;

	G4Box *solidBarraX = new G4Box("solidBarraX", barLong,  barAncho, barGrosor);
	G4Box *solidBarraY = new G4Box("solidBarraY", barAncho, barLong,  barGrosor);

	logicBarraX = new G4LogicalVolume(solidBarraX, plastico, "logicBarraX");
	logicBarraY = new G4LogicalVolume(solidBarraY, plastico, "logicBarraY");

	// ── Posicionamiento de capas ────────────────────────────────────────────
	// Cara trasera del hierro: z = 70 cm
	// Gap post-Fe: 30 cm
	// Capa 1 centro z: 70 + 30 + 0.5 = 100.5 cm
	// Capa 2 centro z: 100.5 + 3.0   = 103.5 cm (separación entre centros = 3 cm)
	G4double ironBackFace = 2.0 * ironHalf;          // 70 cm
	G4double gapPostFe    = 30.0*cm;
	G4double zCapa1 = ironBackFace + gapPostFe + barGrosor;  // 100.5 cm
	G4double zCapa2 = zCapa1 + 3.0*cm;                       // 103.5 cm

	// 20 barras por capa, juntas: paso = 5 cm = 2×barAncho
	// Posiciones del centro: -47.5, -42.5, ..., +42.5, +47.5 cm
	for(G4int i = 0; i < 20; i++)
	{
		G4double desp = -47.5*cm + i * 5.0*cm;

		// Capa 1: barras a lo largo de X, desplazadas en Y
		// Copy number 0-19
		new G4PVPlacement(0, G4ThreeVector(0., desp, zCapa1),
		                  logicBarraX, "physBarraX", logicWorld, false, i, true);

		// Capa 2: barras a lo largo de Y, desplazadas en X (perpendiculares)
		// Copy number 20-39
		new G4PVPlacement(0, G4ThreeVector(desp, 0., zCapa2),
		                  logicBarraY, "physBarraY", logicWorld, false, i + 20, true);
	}

	return physWorld;
}

void MyDetectorConstruction::ConstructSDandField()
{
	MySensitiveDetector *sensDet = new MySensitiveDetector("SensitiveDetector");

	if(logicBarraX != nullptr)
		logicBarraX->SetSensitiveDetector(sensDet);
	if(logicBarraY != nullptr)
		logicBarraY->SetSensitiveDetector(sensDet);
}
