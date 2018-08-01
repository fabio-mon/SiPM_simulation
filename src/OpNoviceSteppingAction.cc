
#include "OpNoviceSteppingAction.hh"
#include "B4cEventAction.hh"
#include "G4SystemOfUnits.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4Electron.hh"


#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "SiPM_SD.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


OpNoviceSteppingAction::OpNoviceSteppingAction(B4cEventAction* eventAction,G4String CutOption)
: G4UserSteppingAction(),
fEventAction(eventAction),
fCutOption(CutOption)
{ 
  fEventNumber = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceSteppingAction::~OpNoviceSteppingAction()
{ ; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceSteppingAction::UserSteppingAction(const G4Step* step)
{
  G4int eventNumber = G4RunManager::GetRunManager()-> GetCurrentEvent()->GetEventID();
  
  if (eventNumber != fEventNumber) 
     fEventNumber = eventNumber;
   
  G4StepPoint* point1 = step->GetPreStepPoint();
  G4StepPoint* point2 = step->GetPostStepPoint();
	
    	//get the name of the volume in which the particles is at the begin of the step
  G4TouchableHandle touch1 = point1->GetTouchableHandle();
  G4VPhysicalVolume* pre_volume = touch1->GetVolume();
  G4String pre_volName = "" ; 
  if ( pre_volume ) 
     pre_volName = pre_volume -> GetName () ;

    	//get the name of the volume in which the particles is at the end of the step
   G4TouchableHandle touch2 = point2->GetTouchableHandle();
   G4VPhysicalVolume* post_volume = touch2->GetVolume();
   G4String post_volName = "" ; 
   if ( post_volume ) 
      post_volName = post_volume -> GetName () ;

  G4Track* track = step->GetTrack();
  G4String ParticleName = track->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
  G4TrackStatus track_status = track->GetTrackStatus(); 
  
  //std::cout<<ParticleName<<"\t"<<pre_volName<<" ---> "<< post_volName<<std::endl;
  //if(track -> GetTrackStatus () == fStopAndKill || track -> GetTrackStatus () == fKillTrackAndSecondaries)
  //   std::cout<<"\tKILL"<<std::endl;

  if (ParticleName == "opticalphoton")
  {
    //if(post_volName == "Optic_window" && pre_volName == "World")
    //  G4cout<<"reflection"<<G4endl;
    /*if((pre_volName == "Optic_window" && post_volName != "Optic_Glue") || (post_volName == "Optic_window" && pre_volName != "Optic_Glue"))
     { 
        G4cout<<pre_volName<< " -> "<<post_volName<<G4endl;
        if(track -> GetTrackStatus () == fStopAndKill || track -> GetTrackStatus () == fKillTrackAndSecondaries)
           G4cout<<"\tKILL\n"<<G4endl;
     }
    */
    /*if(pre_volName == "SiPM"  && post_volName == "SiPM")
     {
       std::cout<<pre_volName <<" -> "<<post_volName<<std::endl;
        if(track -> GetTrackStatus () == fStopAndKill || track -> GetTrackStatus () == fKillTrackAndSecondaries)
           std::cout<<"\tKILL"<<std::endl;
	   }*/
     /*
     if(pre_volName == "Optic_Glue" && post_volName!="Crystal")
     {
        std::cout<<"Glue -> "<<post_volName<<std::endl;
        if(track -> GetTrackStatus () == fStopAndKill || track -> GetTrackStatus () == fKillTrackAndSecondaries)
           std::cout<<"\tKILL"<<std::endl;
     }     
     if(pre_volName == "SiPM" && post_volName=="Optic_Glue")
        std::cout<<"SiPM reflection\n"<<std::endl;
     //if(pre_volName == "SiPM")
     //   std::cout<<"SiPM -> "<<post_volName<<std::endl;*/
 /*    if(pre_volName == "SiPM_border" && post_volName!="Optic_Glue")
        assert(track -> GetTrackStatus () == fStopAndKill || track -> GetTrackStatus () == fKillTrackAndSecondaries);
 */
     if ((track_status == fStopAndKill || track -> GetTrackStatus () == fKillTrackAndSecondaries)&& pre_volName == "Optic_window" && post_volName=="SiPM")
     {
        fEventAction->IncreaseDetPhot();
        G4SDManager* SDman = G4SDManager::GetSDMpointer();
        G4String sdName="SiPM_detector";
        SiPM_SD* SD = (SiPM_SD*)SDman->FindSensitiveDetector(sdName);
        if(SD)SD->ProcessHits_manually(step,NULL);
        //track->SetTrackStatus(fStopAndKill);
/*
        if(track -> GetCreatorProcess() -> GetProcessName() == "Scintillation")
        {   
           OpNoviceCreateTree::Instance()->NPhot_SiPM += 100; 
           OpNoviceCreateTree::Instance()->NPhot_SiPM_scint += 100; 
        }
        else
        {
           OpNoviceCreateTree::Instance()->NPhot_SiPM += 1; 
           OpNoviceCreateTree::Instance()->NPhot_SiPM_ch += 1;
        } */
     }
     else 
        if ((pre_volName == "World" && post_volName=="" && track_status == fStopAndKill) || (pre_volName == "World" && post_volName=="SiPM"))
        {
           fEventAction->IncreaseLostPhot();
           //track->SetTrackStatus(fStopAndKill);
           //std::cout<<pre_volName<<" -> "<<post_volName<<"\tKILL"<<std::endl;
          /* if(track -> GetCreatorProcess() -> GetProcessName() == "Scintillation")
              OpNoviceCreateTree::Instance()->NPhot_lost += 100; 
           else
              OpNoviceCreateTree::Instance()->NPhot_lost += 1; */
        }
        else
           if (pre_volName == "Crystal" && post_volName!="SiPM" &&  post_volName!="Optic_Glue" && post_volName!="Optic_window"  && track_status == fStopAndKill)
           {
              fEventAction->IncreaseAbsPhot();
            /*  if(track -> GetCreatorProcess() -> GetProcessName() == "Scintillation")
              {
                 OpNoviceCreateTree::Instance()->NPhot_abs += 100; 
                 OpNoviceCreateTree::Instance()->NRefl_abs -> push_back(particle_refl[track->GetTrackID()]);
                 OpNoviceCreateTree::Instance()->phot_length -> push_back(track->GetTrackLength());
              }
              else
                 OpNoviceCreateTree::Instance()->NPhot_abs += 1; */
           }
     return;
  }
  
  else
    if ((ParticleName == "mu-" || ParticleName == "mu+" || ParticleName == "pi-" || ParticleName == "pi+") && pre_volName == "Crystal")
    {
       
       fEventAction->AddEdep(step->GetTotalEnergyDeposit());
       fEventAction->AddLength(step->GetStepLength());
       /*if(pre_volName == "World" && post_volName=="Crystal")
       {
         G4ThreeVector step_pos (point2 -> GetPosition());
         OpNoviceCreateTree::Instance()->mu_x_hit = step_pos.getX();
         OpNoviceCreateTree::Instance()->mu_y_hit = step_pos.getY();
       }*/
    }
    /*else 
       std::cout<<ParticleName<<"\t"<<(track->GetKineticEnergy())/CLHEP::MeV<<" MeV"<<"\t"<<pre_volName<<"->"<<post_volName<<std::endl;*/

  const std::vector<const G4Track*>* secondaries = step->GetSecondaryInCurrentStep();
/*
  if (secondaries->size()>0) 
  {
     for(unsigned int i=0; i<secondaries->size(); ++i)
     {
        if (secondaries->at(i)->GetParentID()>0) 
        {
           if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
           {
              fEventAction->IncreaseTotPhot();
           }
          /* else
              if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition() == G4Electron::ElectronDefinition() && pre_volName == "Crystal")
                 OpNoviceCreateTree::Instance() -> N_electr += 1;
        }
     }
  }*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
