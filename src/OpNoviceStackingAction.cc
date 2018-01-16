
#include "OpNoviceStackingAction.hh"

#include "G4VProcess.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "ConfigFile.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceStackingAction::OpNoviceStackingAction(B4cEventAction* eventAction, std::string configFileName)
  : G4UserStackingAction(),
    fEventAction(eventAction)   
{
   ConfigFile config (configFileName) ;
   if (config.keyExists("CutOption"))
      fCutOption = config.read<std::string> ("CutOption");  
   else
      fCutOption = "";

   if(fCutOption == "Timing")
   {
      if (config.keyExists("TimeCut"))
         fTimeCut = config.read<float> ("TimeCut")*ns;  
      else
         fTimeCut = 34.*ns;
   } 
}

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
         if(fCutOption == "Timing" && aTrack->GetGlobalTime() > fTimeCut)
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
