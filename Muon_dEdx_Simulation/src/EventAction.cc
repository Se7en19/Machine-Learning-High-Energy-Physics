//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any  representation or  warranty, express or implied,*
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and available *
// * URL above for the full disclaimer and the limitation of liability.*
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your  *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#include "EventAction.hh"
#include "TrackerHit.hh"
#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4TrajectoryContainer.hh"
#include "G4ios.hh"


namespace B2
{

void EventAction::BeginOfEventAction(const G4Event*) {}

void EventAction::EndOfEventAction(const G4Event* event)
{
  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID("TrackerHitsCollection");
  auto hce=event->GetHCofThisEvent();
  if(!hce) return;

  auto hitsCollection=
    static_cast<B2::TrackerHitsCollection*>(hce->GetHC(hcID));
  if (!hitsCollection) return;

  auto analysisManager = G4AnalysisManager::Instance();
  G4int eventID = event->GetEventID();
  std::size_t nHits = hitsCollection->entries();

  const G4int nChambers = 5;
    G4double totalEdep[nChambers]   = {0.};
    G4double totalLength[nChambers] = {0.};
    G4double rigidity = 0.;
    G4String particleName = "";
    G4bool   primaryFound = false;


  for (std::size_t i=0; i<nHits; i++){
        TrackerHit* hit = (*hitsCollection)[i];

        if (hit->GetTrackID() != 1) continue;

        G4int chamb = hit->GetChamberNb();
        if (chamb < 0 || chamb >= nChambers) continue;

        totalEdep  [chamb] += hit->GetEdep();
        totalLength[chamb] += hit->GetStepLength();

        if (!primaryFound) {
            rigidity     = hit->GetRigidity();
            particleName = hit->GetParticleName();
            primaryFound = true;
        }
  }
  if (!primaryFound) return;

    // ── Llenar ntuple: un entry por cámara ────────────────────────────
    for (G4int chamb = 0; chamb < nChambers; chamb++) {
        if (totalLength[chamb] <= 0.) continue;

        G4double dedx = totalEdep[chamb] / totalLength[chamb];
        if (dedx <= 0.) continue;

        analysisManager->FillNtupleIColumn(0, 0, eventID);
        analysisManager->FillNtupleIColumn(0, 1, 1);          // trackID primaria
        analysisManager->FillNtupleIColumn(0, 2, chamb);
        analysisManager->FillNtupleDColumn(0, 3, dedx);
        analysisManager->FillNtupleDColumn(0, 4, rigidity);
        analysisManager->FillNtupleDColumn(0, 5, rigidity);   // |p| ≈ rigidity para muones
        analysisManager->FillNtupleSColumn(0, 6, particleName);
        analysisManager->AddNtupleRow(0);
    }

  if (eventID < 100 || eventID % 100 == 0) {
        G4cout << ">>> Event: " << eventID
               << "  hits: " << nHits << G4endl;
  }
}



}  // namespace B2
