
#ifndef OpNoviceDetectorConstruction_5_h
#define OpNoviceDetectorConstruction_5_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class OpNoviceDetectorConstruction_5 : public G4VUserDetectorConstruction
{
  public:
    OpNoviceDetectorConstruction_5(std::string configFileName);
    virtual ~OpNoviceDetectorConstruction_5();

  public:
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

  private:
    G4String fPDEoption;  
    G4int fsurface_type;
    G4double fwrapping_refl;
    G4double fSigmaAlpha;

    G4double fExpHall_x;
    G4double fExpHall_y;
    G4double fExpHall_z;

    G4double fCryst_x;
    G4double fCryst_y;
    G4double fCryst_z;

    G4double fSiPM_x;
    G4double fSiPM_y;
    G4double fSiPM_z;
    
    G4double fGlue_x;
    G4double fGlue_y;
    G4double fGlue_z;

    G4double fSiPM_window_x;
    G4double fSiPM_window_y;
    G4double fSiPM_window_z;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*OpNoviceDetectorConstruction_5_h*/
