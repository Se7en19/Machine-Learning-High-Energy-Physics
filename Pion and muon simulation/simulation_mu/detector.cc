#include "detector.hh"

#include "G4RunManager.hh"

MySensitiveDetector::MySensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{}

MySensitiveDetector::~MySensitiveDetector()
{}

G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist)
{
	// Solo registrar la traza primaria (TrackID=1)
	// Esto define el hit: la partícula primaria entra al centellador y deposita energía
	G4int trackID = aStep->GetTrack()->GetTrackID();
	if (trackID != 1) return false;

	G4double ekin = aStep->GetPreStepPoint()->GetKineticEnergy();

	G4double edep = aStep->GetTotalEnergyDeposit();

	G4double stepLength = aStep->GetStepLength();

	G4double dEdx = 0.;
	if (stepLength > 0.) {
		dEdx = edep / stepLength;
	}

	// Tiempo global en nanosegundos
	G4double tof = aStep->GetPreStepPoint()->GetGlobalTime();

	// Longitud total recorrida por la partícula
	G4double trackLength = aStep->GetTrack()->GetTrackLength();

	// Dirección antes y después del step (scattering)
	G4ThreeVector dirPre  = aStep->GetPreStepPoint()->GetMomentumDirection();
	G4ThreeVector dirPost = aStep->GetPostStepPoint()->GetMomentumDirection();
	G4double scatteringAngle = dirPre.angle(dirPost);

	// Posición del volumen detector
	const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable();
	G4VPhysicalVolume  *physVol   = touchable->GetVolume();
	G4ThreeVector posDetector     = physVol->GetTranslation();

	G4AnalysisManager *man = G4AnalysisManager::Instance();

	const G4Event *evt = G4RunManager::GetRunManager()->GetCurrentEvent();
	G4int eventID = evt->GetEventID();

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
	man->AddNtupleRow(0);

	return true;
}
