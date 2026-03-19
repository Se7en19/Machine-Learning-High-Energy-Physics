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

	std::ostringstream oss;
	oss << "output_run" << aRun->GetRunID() << ".root";
	G4String fname = oss.str();
	man->OpenFile(fname);

	man->CreateNtuple("Hits", "Hits");
	man->CreateNtupleIColumn("fEvent");
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
	man->FinishNtuple(0);
}

void MyRunAction::EndOfRunAction(const G4Run*)
{
	G4AnalysisManager *man = G4AnalysisManager::Instance();
	man->Write();
	man->CloseFile();
}
