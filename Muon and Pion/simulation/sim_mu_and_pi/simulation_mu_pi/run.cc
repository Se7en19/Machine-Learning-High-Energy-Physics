#include "run.hh"
#include "G4Run.hh"
#include <sstream>

MyRunAction::MyRunAction()
{}

MyRunAction::~MyRunAction()
{}

void MyRunAction::BeginOfRunAction(const G4Run* aRun)
{
	G4AnalysisManager *man = G4AnalysisManager::Instance();

	// Un archivo por run, árbol siempre llamado "Hits"
	std::ostringstream oss;
	oss << "output_run" << aRun->GetRunID() << ".root";
	//oss << "../../../Classifier/data/muon/output_run" << aRun->GetRunID() << ".root";
	G4String fname = oss.str();
	G4cout << "!!! CREANDO ARCHIVO EN: " << fname << G4endl;
	man->OpenFile(fname);

	man->CreateNtuple("Hits", "Hits");
	man->CreateNtupleDColumn("fX");
	man->CreateNtupleDColumn("fY");
	man->CreateNtupleDColumn("fZ");
	man->CreateNtupleDColumn("fEdep");
	man->CreateNtupleDColumn("fdEdx");
	man->CreateNtupleDColumn("Ekin");
	man->CreateNtupleDColumn("TOF");
	man->CreateNtupleDColumn("TrackLength");
	man->CreateNtupleDColumn("ScatteringAng");
	man->CreateNtupleDColumn("Momentum");
	man->CreateNtupleIColumn("fEvent");
	man->CreateNtupleIColumn("fPDG");
	man->FinishNtuple(0);
}

void MyRunAction::EndOfRunAction(const G4Run* aRun)
{

	G4AnalysisManager *man = G4AnalysisManager::Instance();
	
	G4int nofEvents = aRun->GetNumberOfEvent();
	G4cout << "Total de eventos generados: " << nofEvents << G4endl;
	
	man->Write();
	man->CloseFile();
}
