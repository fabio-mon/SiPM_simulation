
#include "B4cEventAction.hh"
#include "SiPM_SD.hh"
#include "B4cCalorHit.hh"
#include "B4Analysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
B4cEventAction::B4cEventAction()
 : G4UserEventAction(),
   fAbsHCID(-1),
   fGapHCID(-1),
   fPh_time(1, 0.),
   //fPh_energy(5000, 0.), 
   fPh_creator_proc(1,-1),
   fPh_SiPM_number(1,-1),
   //fPh_x_hit(5000,-20),
   //fPh_y_hit(5000,-20),
   fPh_detected_Nb(0),
   fPh_tot_Nb(0),
   fPh_lost_Nb(0),
   fPh_abs_Nb(0),
   fmu_Edep(0.),
   fmu_TrackLength(0.)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cEventAction::~B4cEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cEventAction::IncreaseDetPhot()
{
   fPh_detected_Nb++;
}
  
void B4cEventAction::IncreaseTotPhot()
{
   fPh_tot_Nb++;
}

void B4cEventAction::IncreaseLostPhot()
{
   fPh_lost_Nb++;
}

void B4cEventAction::IncreaseAbsPhot()
{
   fPh_abs_Nb++;
}

void B4cEventAction::AddEdep(const G4double &Edep)
{
  fmu_Edep+=Edep;
}

void B4cEventAction::AddLength(const G4double &StepLength)
{
  fmu_TrackLength+=StepLength;
}


B4cCalorHitsCollection* B4cEventAction::GetHitsCollection(G4int hcID,const G4Event* event) const
{
  B4cCalorHitsCollection* hitsCollection = static_cast<B4cCalorHitsCollection*>(event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) 
  {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID; 
    G4Exception("B4cEventAction::GetHitsCollection()","MyCode0003", FatalException, msg);
  }         

  return hitsCollection;
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cEventAction::PrintEventStatistics(const G4double &mu_x,const G4double &mu_y) const
{
  // print event statistics
  G4cout
     << "   Muon vertex position ( x , y ) = ( " 
     << mu_x/CLHEP::mm <<" , "<< mu_y/CLHEP::mm<<" ) mm"
     <<G4endl
     << "       Total track length = " 
     << fmu_TrackLength/CLHEP::mm<<" mm"
     << G4endl
     << "        Energy loss =  " 
     << fmu_Edep/CLHEP::MeV<<" MeV"
     << G4endl
     << "#Photons total = " 
     << fPh_tot_Nb
     << G4endl
     << "#Photons detected = " 
     << fPh_detected_Nb
     << G4endl
     << "Arrival time of last detected photon ("<<fPh_time.size()<<"th phot) = " 
     << fPh_time.back()/CLHEP::ns<<" ns"
     << G4endl
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cEventAction::BeginOfEventAction(const G4Event* event)
{
   //G4int eventID = event->GetEventID();
   //G4cout<<"Event "<<eventID<<"starts"<<G4endl;
   std::cout<<std::endl;
   fPh_detected_Nb=0;
   fPh_tot_Nb=0;
   fPh_lost_Nb=0;
   fPh_abs_Nb=0;
   fmu_Edep=0.;
   fmu_TrackLength=0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cEventAction::Freespace()
{
   //for(std::vector<G4float>::iterator it_time = fPh_time.begin();it_time<fPh_time.end();++it_time)
      //std::cout<<*it_time<<std::endl;
   std::cout << "initial size: " << fPh_time.size() << "\n";
   std::vector<G4float>::iterator it_time = fPh_time.end()-1;
   std::vector<G4float>::iterator it_energy = fPh_energy.end()-1;
   std::vector<G4float>::iterator it_x_hit = fPh_x_hit.end()-1;
   std::vector<G4float>::iterator it_y_hit = fPh_y_hit.end()-1;
   std::vector<G4int>::iterator it_creator_proc = fPh_creator_proc.end()-1;
   std::vector<G4int>::iterator it_SiPM_number = fPh_SiPM_number.end()-1;
   while(*it_time==0.)
   {
      fPh_time.pop_back();
      fPh_energy.pop_back();
      fPh_x_hit.pop_back();
      fPh_y_hit.pop_back();
      fPh_creator_proc.pop_back();
      fPh_SiPM_number.pop_back();
      it_time = fPh_time.end()-1;
      it_energy = fPh_energy.end()-1;
      it_x_hit = fPh_x_hit.end()-1;
      it_y_hit = fPh_y_hit.end()-1;
      it_creator_proc = fPh_creator_proc.end()-1;
      it_SiPM_number = fPh_SiPM_number.end()-1;
   }
   std::cout << "final size: " << fPh_time.size() << "\n";
   //for(std::vector<G4float>::iterator it_time = fPh_time.begin();it_time<fPh_time.end();++it_time)
      //std::cout<<*it_time<<std::endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cEventAction::EndOfEventAction(const G4Event* event)
{  
  // Get hits collections IDs (only once)
  if ( fAbsHCID == -1 ) 
  {
    fAbsHCID = G4SDManager::GetSDMpointer()->GetCollectionID("SiPMHitsCollection");
  }

  // Get hits collections
  B4cCalorHitsCollection* absoHC = GetHitsCollection(fAbsHCID, event);

  // Get hit with total values
  B4cCalorHit* absoHit = (*absoHC)[absoHC->entries()-1];
 
  G4double mu_x = event->GetPrimaryVertex()->GetPosition().x();
  G4double mu_y = event->GetPrimaryVertex()->GetPosition().y();

  // Print per event (modulo n)
  //
  G4int eventID = event->GetEventID();
  
  absoHit->GetTime(fPh_time);
  //fPh_energy = absoHit->GetEnergy();
  absoHit->GetCreatorProc(fPh_creator_proc);
  absoHit->GetSiPM_number(fPh_SiPM_number);
  //fPh_SiPM_number =  absoHit->GetSiPMNumber();
  //fPh_x_hit = absoHit->GetPh_x_hit();
  //fPh_y_hit = absoHit->GetPh_y_hit();
  //Freespace();

  //G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  //if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) 
  //{
    G4cout << "---> End of event: " << eventID << G4endl;     
    PrintEventStatistics(mu_x,mu_y);
  //}
  
  // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 
  // fill ntuple
  analysisManager->FillNtupleIColumn(0, eventID);
  analysisManager->FillNtupleFColumn(1, mu_x);
  analysisManager->FillNtupleFColumn(2, mu_y);
  analysisManager->FillNtupleFColumn(3, fmu_Edep);
  analysisManager->FillNtupleFColumn(4, fmu_TrackLength);
  analysisManager->FillNtupleIColumn(5, fPh_tot_Nb);  
  analysisManager->FillNtupleIColumn(6, fPh_detected_Nb);
  analysisManager->FillNtupleIColumn(7, fPh_lost_Nb);
  analysisManager->FillNtupleIColumn(8, fPh_abs_Nb);
//vector coloumns should be automatically filled
  analysisManager->AddNtupleRow();  
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
