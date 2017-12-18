
#ifndef OpNoviceDetectorConstruction_1_h
#define OpNoviceDetectorConstruction_1_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "MyMaterials.hh"

class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;

class OpNoviceDetectorConstruction_1 : public G4VUserDetectorConstruction
{
  public:
    OpNoviceDetectorConstruction_1(std::string configFileName);
    virtual ~OpNoviceDetectorConstruction_1();

  public:
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();
     
  private:
    //G4VPhysicalVolume* DefineVolumes();
    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger; 
                                      // magnetic field messenger
    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps
    G4String fPDEoption;  
    G4int fsurface_type;
    G4double fwrapping_refl;
    G4double fSigmaAlpha,fSS,fSL,fBS;

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

#endif
