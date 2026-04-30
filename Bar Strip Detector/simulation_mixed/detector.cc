#include "detector.hh"
#include "G4RunManager.hh"

MySensitiveDetector::MySensitiveDetector(G4String name)
	: G4VSensitiveDetector(name)
{}

MySensitiveDetector::~MySensitiveDetector()
{}

G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist)
{
	// Solo registrar la traza primaria (TrackID=1)
	G4int trackID = aStep->GetTrack()->GetTrackID();
	if(trackID != 1) return false;

	G4double edep = aStep->GetTotalEnergyDeposit();
	G4double stepLength = aStep->GetStepLength();

	G4double dEdx = 0.;
	if(stepLength > 0.)
		dEdx = edep / stepLength;

	G4double ekin        = aStep->GetPreStepPoint()->GetKineticEnergy();
	G4double tof         = aStep->GetPreStepPoint()->GetGlobalTime();
	G4double trackLength = aStep->GetTrack()->GetTrackLength();

	G4ThreeVector dirPre  = aStep->GetPreStepPoint()->GetMomentumDirection();
	G4ThreeVector dirPost = aStep->GetPostStepPoint()->GetMomentumDirection();
	G4double scatteringAngle = dirPre.angle(dirPost);

	// Posición y número de copia del volumen sensible
	const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable();
	G4VPhysicalVolume  *physVol   = touchable->GetVolume();
	G4ThreeVector posDetector     = physVol->GetTranslation();

	// Identificación de capa y barra:
	// copy number 0-19  → Capa 1 (barras a lo largo de X)
	// copy number 20-39 → Capa 2 (barras a lo largo de Y)
	G4int copyNo  = physVol->GetCopyNo();
	G4int layerID = (copyNo < 20) ? 0 : 1;
	G4int barID   = copyNo % 20;

	const G4Event *evt    = G4RunManager::GetRunManager()->GetCurrentEvent();
	G4int          eventID = evt->GetEventID();

	// particleID: 0 = mu+, 1 = pi+
	G4String pname = aStep->GetTrack()->GetDefinition()->GetParticleName();
	G4int particleID = (pname == "mu+") ? 0 : 1;

	// ConeAngle: ángulo entre la dirección inicial de la traza y el eje z del haz (rad)
	// Usa GetVertexMomentumDirection() para obtener la dirección original en el vértice,
	// antes de cualquier dispersión en el absorbedor.
	G4ThreeVector vdir      = aStep->GetTrack()->GetVertexMomentumDirection();
	G4double      coneAngle = vdir.angle(G4ThreeVector(0., 0., 1.));

	G4AnalysisManager *man = G4AnalysisManager::Instance();

	man->FillNtupleDColumn(0, posDetector[0]);
	man->FillNtupleDColumn(1, posDetector[1]);
	man->FillNtupleDColumn(2, posDetector[2]);
	man->FillNtupleDColumn(3, edep);
	man->FillNtupleDColumn(4, dEdx);
	man->FillNtupleDColumn(5, ekin);
	man->FillNtupleDColumn(6, tof);
	man->FillNtupleDColumn(7, trackLength);
	man->FillNtupleDColumn(8, scatteringAngle);
	man->FillNtupleDColumn(9, aStep->GetPreStepPoint()->GetMomentum().mag());
	man->FillNtupleIColumn(10, eventID);
	man->FillNtupleIColumn(11, layerID);
	man->FillNtupleIColumn(12, barID);
	man->FillNtupleIColumn(13, particleID);
	man->FillNtupleDColumn(14, coneAngle);
	man->AddNtupleRow(0);

	return true;
}
