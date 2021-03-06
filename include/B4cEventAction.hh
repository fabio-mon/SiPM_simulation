
#ifndef B4cEventAction_h
#define B4cEventAction_h 1

#include "G4UserEventAction.hh"

#include "B4cCalorHit.hh"

#include "globals.hh"

/// Event action class
///
/// In EndOfEventAction(), it prints the accumulated quantities of the energy 
/// deposit and track lengths of charged particles in Absober and Gap layers 
/// stored in the hits collections.

class B4cEventAction : public G4UserEventAction
{
public:
  B4cEventAction();
  virtual ~B4cEventAction();

  virtual void  BeginOfEventAction(const G4Event* event);
  virtual void    EndOfEventAction(const G4Event* event);
  void IncreaseDetPhot();
  void IncreaseTotPhot();
  void IncreaseLostPhot();
  void IncreaseAbsPhot();
  void AddEdep(const G4double &Edep);
  void AddLength(const G4double &StepLength);

  G4int GetDetPhot() { return fPh_detected_Nb; }
  G4int GetAbsPhot() { return fPh_abs_Nb; }
  G4int GetTotPhot() { return fPh_tot_Nb; }
  G4int GetLosttPhot() { return fPh_lost_Nb; }
  G4float GetEdep() { return fmu_Edep; }
  G4float GetLength() { return fmu_TrackLength; }
  std::vector<G4float>& GetTimeVec() { return fPh_time; }
  std::vector<G4float>& GetEnergyVec() { return fPh_energy; }
  std::vector<G4int>& GetCreatorProcVec() { return fPh_creator_proc; }
  std::vector<G4int>& GetSiPM_numberVec() { return fPh_SiPM_number; }
  std::vector<G4float>& GetPhXVec() { return fPh_x_hit; }
  std::vector<G4float>& GetPhYVec() { return fPh_y_hit; }
    
private:
  // methods
  void Freespace();
  B4cCalorHitsCollection* GetHitsCollection(G4int hcID, const G4Event* event) const;
  void PrintEventStatistics(const G4double &mu_x,const G4double &mu_y) const;
  
  // data members                   
  G4int  fAbsHCID;
  G4int  fGapHCID;
  G4int  fPh_detected_Nb;
  G4int  fPh_tot_Nb;
  G4int  fPh_lost_Nb;
  G4int  fPh_abs_Nb;
  G4float  fmu_Edep;
  G4float  fmu_TrackLength;
  std::vector<G4float> fPh_time;
  std::vector<G4float> fPh_energy; 
  std::vector<G4int> fPh_creator_proc;
  std::vector<G4int> fPh_SiPM_number;
  std::vector<G4float> fPh_x_hit; 
  std::vector<G4float> fPh_y_hit; 
};
                     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
