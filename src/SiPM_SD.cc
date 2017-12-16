
#include "SiPM_SD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include <cmath>
#include "Randomize.hh"

double PDEKETEK3x3_15um(const double &energy)
{
   //cout<<"energy = "<<energy<<endl;
   //cout<<"PDE = "<<0.39784 * exp(-(energy - 2.9159)*(energy - 2.9159) / (2 * 0.5512 * 0.5512))<<endl;
   //cout<<"Reflec = "<<SiliconReflectivity(energy)<<endl;
   //cout<<"TOT PDE = "<<0.39784/0.48/SiliconReflectivity(energy) * exp(-(energy - 2.9159)*(energy - 2.9159) / (2 * 0.5512 * 0.5512))<<endl<<endl;
   return 0.39784/*/0.6*0.04/SiliconReflectivity(energy) */ * exp(-(energy - 2.9159)*(energy - 2.9159) / (2 * 0.5512 * 0.5512));
}

double PDEKETEKDIFFUSE_15um(const double &energy)
{
   //cout<<"energy = "<<energy<<endl;
   //cout<<"PDE = "<<0.39784 * exp(-(energy - 2.9159)*(energy - 2.9159) / (2 * 0.5512 * 0.5512))<<endl;
   //cout<<"Reflec = "<<SiliconReflectivity(energy)<<endl;
   //cout<<"TOT PDE = "<<0.39784/0.48/SiliconReflectivity(energy) * exp(-(energy - 2.9159)*(energy - 2.9159) / (2 * 0.5512 * 0.5512))<<endl<<endl;
   return 0.39784/0.7*0.16/*/SiliconReflectivity(energy) */ * exp(-(energy - 2.9159)*(energy - 2.9159) / (2 * 0.5512 * 0.5512));
}

double PDEFBKDIFFUSE_20um(const double &energy)
{
   //cout<<"energy = "<<energy<<endl;
   //cout<<"PDE = "<<0.39784 * exp(-(energy - 2.9159)*(energy - 2.9159) / (2 * 0.5512 * 0.5512))<<endl;
   //cout<<"Reflec = "<<SiliconReflectivity(energy)<<endl;
   //cout<<"TOT PDE = "<<0.39784/0.48/SiliconReflectivity(energy) * exp(-(energy - 2.9159)*(energy - 2.9159) / (2 * 0.5512 * 0.5512))<<endl<<endl;
   return 5.*5./(10.*10)*(0.985* exp(-(energy -3.633)*(energy -3.633) / (2 * 1.0795*1.0795))-0.225 +energy*(0.2389) -0.1075*energy*energy);
}

double PDEKETEK3x3_25um(const double &energy)
{
   //cout<<"energy = "<<energy<<endl;
   //cout<<"PDE = "<<0.39784 * exp(-(energy - 2.9159)*(energy - 2.9159) / (2 * 0.5512 * 0.5512))<<endl;
   //cout<<"Reflec = "<<SiliconReflectivity(energy)<<endl;
   //cout<<"TOT PDE = "<<0.39784/0.48/SiliconReflectivity(energy) * exp(-(energy - 2.9159)*(energy - 2.9159) / (2 * 0.5512 * 0.5512))<<endl<<endl;
   return 0.442608/*/0.6*0.04/SiliconReflectivity(energy) */ * exp(-(energy - 2.82216)*(energy - 2.82216) / (2 * 0.51774 * 0.51774)) + 0.0707552 + energy*(-0.0808317) + energy*energy*(0.0274066);
}

double PDEFBK5x5_20um(const double &energy)
{
   //cout<<"energy = "<<energy<<endl;
   //cout<<"PDE = "<<0.39784 * exp(-(energy - 2.9159)*(energy - 2.9159) / (2 * 0.5512 * 0.5512))<<endl;
   //cout<<"Reflec = "<<SiliconReflectivity(energy)<<endl;
   //cout<<"TOT PDE = "<<0.39784/0.48/SiliconReflectivity(energy) * exp(-(energy - 2.9159)*(energy - 2.9159) / (2 * 0.5512 * 0.5512))<<endl<<endl;
   return (0.985* exp(-(energy -3.633)*(energy -3.633) / (2 * 1.0795*1.0795))-0.225 +energy*(0.2389) -0.1075*energy*energy);
}

double PDEHPK_MPPC(const double &energy)
{
  //std::cout<<"ook";
  return 0.37043 * exp(-(energy - 2.69767)*(energy - 2.69767) / (2 * 0.692431 * 0.692431));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SiPM_SD::SiPM_SD(const G4String& name, const G4String& hitsCollectionName, G4int nofCells, G4String PDEoption)
 : G4VSensitiveDetector(name),
   fHitsCollection(0),
   fPDEoption(PDEoption),
   fNofCells(nofCells)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SiPM_SD::~SiPM_SD() 
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SiPM_SD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection = new B4cCalorHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce
  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 

  // Create hits
  // fNofCells for cells + one more for total sums 
  for (G4int i=0; i<fNofCells+1; i++ ) 
  {
    fHitsCollection->insert(new B4cCalorHit());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool SiPM_SD::ProcessHits(G4Step* step, G4TouchableHistory*)
{  
  //std::cout<<"Running default ProcessHits"<<std::endl;
  return false;
}

G4bool SiPM_SD::ProcessHits_manually(const G4Step* step, G4TouchableHistory*)
{
  //std::cout<<"Running default ProcessHits_manually"<<std::endl;
  G4Track* track = step->GetTrack();
  G4StepPoint* point1 = step->GetPreStepPoint();
  G4StepPoint* point2 = step->GetPostStepPoint();
  G4TouchableHistory* touchable = (G4TouchableHistory*)(step->GetPreStepPoint()->GetTouchable());

  G4TouchableHandle touch1 = point1->GetTouchableHandle();
  G4VPhysicalVolume* pre_volume = touch1->GetVolume();
  G4String pre_volName = "" ; 
  if ( pre_volume ) 
     pre_volName = pre_volume -> GetName () ;

  G4TouchableHandle touch2 = point2->GetTouchableHandle();
  G4VPhysicalVolume* post_volume = touch2->GetVolume();
  G4String post_volName = "" ; 
  if ( post_volume ) 
     post_volName = post_volume -> GetName () ;

  G4String ParticleName = track->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
  //std::cout<<"\nProc hit"<<ParticleName<<"\t"<<pre_volName<<std::endl;
  //if(track -> GetTrackStatus () == fStopAndKill || track -> GetTrackStatus () == fKillTrackAndSecondaries)
  //   std::cout<<"KILL\n"<<std::endl;
  //if (ParticleName != "opticalphoton" || pre_volName != "Optic_Glue" || (track->GetTrackStatus()!=fStopAndKill && track -> GetTrackStatus () != fKillTrackAndSecondaries)) return false; 
  // step length
  G4double stepLength = 0.;
  stepLength = step->GetStepLength();

  //photon energy
  G4double PhotEnergy = point2 -> GetTotalEnergy ();

  //photon position
  G4double Ph_x_hit = point2 ->GetPosition().x();
  G4double Ph_y_hit = point2 ->GetPosition().y();

  //photon creation process
  G4int CreatorProcess=-1;
  if(track -> GetCreatorProcess() -> GetProcessName() == "Scintillation")
    CreatorProcess=0;
  else
    CreatorProcess=1;

  //Arrival time
  G4double time = track->GetGlobalTime(); 

  //SiPM number
  //G4cout<<pre_volName<<"->"<<post_volName<<"\tSiPM Copy Nb = "<<post_volume->GetCopyNo()<<G4endl;
  G4int SiPMNumber = post_volume->GetCopyNo();

  //if ( stepLength == 0. ) return false;      

  // Get hit accounting data for this cell
  B4cCalorHit* hit = (*fHitsCollection)[SiPMNumber];
  if ( ! hit ) 
  {
    G4ExceptionDescription msg;
    msg << "Cannot access hit " << SiPMNumber; 
    G4Exception("SiPM_SD::ProcessHits()","MyCode0004", FatalException, msg);
  }         

  // Get hit for total accounting
  B4cCalorHit* hitTotal = (*fHitsCollection)[fHitsCollection->entries()-1];

  //Apply PDE
  G4double rndm;
  if (fPDEoption == "KETEK3x3_25um")
  {
     rndm = G4UniformRand();
     //std::cout<<"rndm = "<<rndm<<"\tPDE ( "<<PhotEnergy/CLHEP::eV <<") = "<<PDEKETEK3x3_25um(PhotEnergy/CLHEP::eV)<<std::endl;
     if(rndm<PDEKETEK3x3_25um(PhotEnergy/CLHEP::eV))
     {
        // Add values
        hit->AddPhotonInfo(time,PhotEnergy,CreatorProcess,SiPMNumber,Ph_x_hit,Ph_y_hit);
        hitTotal->AddPhotonInfo(time,PhotEnergy,CreatorProcess,SiPMNumber,Ph_x_hit,Ph_y_hit);
        //std::cout<<time<<"\t"<<PhotEnergy<<"\t"<<CreatorProcess<<"\t"<<SiPMNumber<<"\t"<<Ph_x_hit<<"\t"<<Ph_y_hit<<std::endl;    
        return true;
     }
     else return false;
  }
  else
     if(fPDEoption == "KETEK3x3_15um")
     {
        rndm = G4UniformRand();
        if(rndm<PDEKETEK3x3_15um(PhotEnergy/CLHEP::eV))
        {
           // Add values
           hit->AddPhotonInfo(time,PhotEnergy,CreatorProcess,SiPMNumber,Ph_x_hit,Ph_y_hit);
           hitTotal->AddPhotonInfo(time,PhotEnergy,CreatorProcess,SiPMNumber,Ph_x_hit,Ph_y_hit);
           //std::cout<<time<<"\t"<<PhotEnergy<<"\t"<<CreatorProcess<<"\t"<<SiPMNumber<<"\t"<<Ph_x_hit<<"\t"<<Ph_y_hit<<std::endl;    
           return true;
        }
        else return false;
     }
     else
        if(fPDEoption == "FBK5x5_20um")
        {
           rndm = G4UniformRand();
           if(rndm<PDEFBK5x5_20um(PhotEnergy/CLHEP::eV))
           {
              // Add values
              hit->AddPhotonInfo(time,PhotEnergy,CreatorProcess,SiPMNumber,Ph_x_hit,Ph_y_hit);
              hitTotal->AddPhotonInfo(time,PhotEnergy,CreatorProcess,SiPMNumber,Ph_x_hit,Ph_y_hit);
              //std::cout<<time<<"\t"<<PhotEnergy<<"\t"<<CreatorProcess<<"\t"<<SiPMNumber<<"\t"<<Ph_x_hit<<"\t"<<Ph_y_hit<<std::endl;    
              return true;
           }
           else return false;
        }
        else
           if(fPDEoption == "KETEKDIFFUSE_15um")
           {
              rndm = G4UniformRand();
              if(rndm<PDEKETEKDIFFUSE_15um(PhotEnergy/CLHEP::eV))
              {
                 // Add values
                 hit->AddPhotonInfo(time,PhotEnergy,CreatorProcess,SiPMNumber,Ph_x_hit,Ph_y_hit);
                 hitTotal->AddPhotonInfo(time,PhotEnergy,CreatorProcess,SiPMNumber,Ph_x_hit,Ph_y_hit);
                 //std::cout<<time<<"\t"<<PhotEnergy<<"\t"<<CreatorProcess<<"\t"<<SiPMNumber<<"\t"<<Ph_x_hit<<"\t"<<Ph_y_hit<<std::endl;    
                 return true;
              }
              else return false;
           }
           else
              if(fPDEoption == "HPK_MPPC")
              {
                 rndm = G4UniformRand();
                 if(rndm<PDEHPK_MPPC(PhotEnergy/CLHEP::eV))
                 {
                    // Add values
                    hit->AddPhotonInfo(time,PhotEnergy,CreatorProcess,SiPMNumber,Ph_x_hit,Ph_y_hit);
                    hitTotal->AddPhotonInfo(time,PhotEnergy,CreatorProcess,SiPMNumber,Ph_x_hit,Ph_y_hit);
                    //std::cout<<time<<"\t"<<PhotEnergy<<"\t"<<CreatorProcess<</*"\t"<<SiPMNumber<<"\t"<<Ph_x_hit<<"\t"<<Ph_y_hit<<*/std::endl;    
                    return true;
                 }
                 else return false;
              }
              else
                 if(fPDEoption == "PDEFBKDIFFUSE_20um")
                 {
                    rndm = G4UniformRand();
                    if(rndm<PDEFBKDIFFUSE_20um(PhotEnergy/CLHEP::eV))
                    {
                       // Add values
                       hit->AddPhotonInfo(time,PhotEnergy,CreatorProcess,SiPMNumber,Ph_x_hit,Ph_y_hit);
                       hitTotal->AddPhotonInfo(time,PhotEnergy,CreatorProcess,SiPMNumber,Ph_x_hit,Ph_y_hit);
                       //std::cout<<time<<"\t"<<PhotEnergy<<"\t"<<CreatorProcess<</*"\t"<<SiPMNumber<<"\t"<<Ph_x_hit<<"\t"<<Ph_y_hit<<*/std::endl;    
                       return true;
                    }
                    else return false;
                 }      
                 else
                 {
                    hit->AddPhotonInfo(time,PhotEnergy,CreatorProcess,SiPMNumber,Ph_x_hit,Ph_y_hit);
                    hitTotal->AddPhotonInfo(time,PhotEnergy,CreatorProcess,SiPMNumber,Ph_x_hit,Ph_y_hit);
                    return true;
                 }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SiPM_SD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) 
  { 
     G4int nofHits = fHitsCollection->entries();
     G4cout
       << G4endl 
       << "-------->Hits Collection: in this event they are " << nofHits 
       << " hits in the tracker chambers: " << G4endl;
     for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
