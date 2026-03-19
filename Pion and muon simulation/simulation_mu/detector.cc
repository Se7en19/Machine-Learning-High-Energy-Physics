#include "detector.hh"

#include "G4RunManager.hh"

MySensitiveDetector::MySensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{}

MySensitiveDetector::~MySensitiveDetector()
{}

G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist)
{
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

	//G4cout << "Detector position: " << posDetector << G4endl;

	G4AnalysisManager *man = G4AnalysisManager::Instance();

	const G4Event *evt = G4RunManager::GetRunManager()->GetCurrentEvent();
	G4int eventID = evt->GetEventID();

	man->FillNtupleIColumn(0, eventID);
	man->FillNtupleDColumn(1, posDetector[0]);
	man->FillNtupleDColumn(2, posDetector[1]);
	man->FillNtupleDColumn(3, posDetector[2]);
	man->FillNtupleDColumn(4, edep);
	man->FillNtupleDColumn(5, dEdx);
	man->FillNtupleDColumn(6, ekin);
	man->FillNtupleDColumn(7, tof);
	man->FillNtupleDColumn(8, trackLength);
	man->FillNtupleDColumn(9, scatteringAngle);
	man->FillNtupleDColumn(10, aStep->GetPreStepPoint()->GetMomentum().mag());
	man->AddNtupleRow(0);

	return true;
}
