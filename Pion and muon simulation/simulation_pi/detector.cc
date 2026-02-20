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
		
	//G4Track *track = aStep->GetTrack();
	
	//track->SetTrackStatus(fStopAndKill);
	
	// Extrae el tiempo global en nanosegundos
	G4double tof = aStep->GetPreStepPoint()->GetGlobalTime();
	
	// Extrae la longitud total recorrida por la partícula
	G4double trackLength = aStep->GetTrack()->GetTrackLength();
	
	// Dirección que traía la partícula ANTES de entrar
	G4ThreeVector dirPre = aStep->GetPreStepPoint()->GetMomentumDirection();

// Dirección que tomó DESPUÉS de chocar
	G4ThreeVector dirPost = aStep->GetPostStepPoint()->GetMomentumDirection();

// Calculamos el ángulo entre esos dos vectores (el "scattering")
	G4double scatteringAngle = dirPre.angle(dirPost);
	
	G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
	G4StepPoint *posStepPoint = aStep->GetPostStepPoint();
	
	G4ThreeVector posPhoton = preStepPoint->GetPosition();
	
	//G4cout << "Mu  position: " << posPhoton << G4endl; 
	
	const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable();
	
	G4int copyNo = touchable->GetCopyNumber();
	
	//G4cout << "Copy number: " <<copyNo << G4endl;
	
	G4VPhysicalVolume *physVol = touchable->GetVolume();
	G4ThreeVector posDetector = physVol->GetTranslation();
	
	G4cout << "Detector position: " <<posDetector << G4endl;
	
	edep = aStep->GetTotalEnergyDeposit();
	stepLength = aStep->GetStepLength();

	dEdx = 0.;
   	if (stepLength > 0.) {
        dEdx = edep / stepLength;
   	 }
   
	ekin = aStep->GetPreStepPoint()->GetKineticEnergy();
	
	G4AnalysisManager *man = G4AnalysisManager::Instance();
	
	const G4Event* evt = G4RunManager::GetRunManager()->GetCurrentEvent();
	G4int eventID = evt->GetEventID();
	//const G4Event* currentEvent = G4RunManager::GetRunManager()->GetCurrentEvent();
	//G4int evt = currentEvent->GetEventID();

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
	man->AddNtupleRow(0);
	
	return true;
}
