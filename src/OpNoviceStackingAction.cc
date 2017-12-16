
#include "OpNoviceStackingAction.hh"

#include "G4VProcess.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceStackingAction::OpNoviceStackingAction(B4cEventAction* eventAction, G4String CutOption)
  : G4UserStackingAction(),
    fCutOption(CutOption),
    fEventAction(eventAction)   
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceStackingAction::~OpNoviceStackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
OpNoviceStackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  if(aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) // particle is optical photon
  {
    if(aTrack->GetParentID()>0)//particle is secondary
    { 
      if(fCutOption == "LightColl" && aTrack->GetTrackID() % 100 != 0)
         return fKill;
      else
         if(fCutOption == "Timing" && aTrack->GetGlobalTime() > 34.*ns)
            return fKill;
         else
            fEventAction->IncreaseTotPhot();
      /*if(aTrack->GetCreatorProcess()->GetProcessName() == "Scintillation")
      {
        fScintillationCounter++;
                //return fKill;
      }
      
      
      if(aTrack->GetCreatorProcess()->GetProcessName() == "Cerenkov")
      {
         fCerenkovCounter++;
         if(fCerenkovCounter%100!=0)//The number of cerenkov photons is reduced by a factor 100, as much as the factor for scintillation photons
         {
            OpNoviceCreateTree::Instance()->NPhot_tot--;
            OpNoviceCreateTree::Instance()->NPhot_ch--;
            return fKill;
         }
      }
      */
         
      /*
      if(aTrack->GetTotalEnergy()<375.*nm)
      {
        OpNoviceCreateTree::Instance()->NPhot_abs += 1;
        return fKill;
      }
      */
    }
  }
  
  return fUrgent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
void OpNoviceStackingAction::NewStage()
{
  G4cout << "Number of Scintillation photons produced in this event : "
         << fScintillationCounter << G4endl;
  G4cout << "Number of Cerenkov photons produced in this event : "
         << fCerenkovCounter << G4endl;
  
  //G4cout << "Number of optical photon hitting the SiPM : "<<fScoring<<G4endl;
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
void OpNoviceStackingAction::PrepareNewEvent()
{
  fScintillationCounter = 0;
  fCerenkovCounter = 0;
  fScoring=0;
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
