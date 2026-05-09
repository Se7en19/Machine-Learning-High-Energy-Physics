#ifndef DETECTOR_HH
#define DETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include "G4AnalysisManager.hh"

constexpr G4double kDEDXThreshold = 0.05;  // MeV/mm — umbral mínimo de dE/dx

class MySensitiveDetector : public G4VSensitiveDetector
{
public:
	MySensitiveDetector(G4String);
	~MySensitiveDetector();

private:
	virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *);
};

#endif
