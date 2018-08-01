
#ifndef DetectorConstruction_planar_h
#define DetectorConstruction_planar_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "MyMaterials.hh"
#include "G4OpticalSurface.hh"

class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;

class DetectorConstruction_planar : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction_planar(std::string configFileName);
    virtual ~DetectorConstruction_planar();

  public:
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();
     
  private:
    void SetSurfaceProperties(G4OpticalSurface* opSurface, int surface_type, double SigmaAlpha, double refl, double SL, double SS, double BS);

    //G4VPhysicalVolume* DefineVolumes();
    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger; 
                                      // magnetic field messenger
    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps
    G4String fPDEoption;  
    G4double fweight;
    G4double ftilt_angle;
  //bool fGlueWindowEverywhere;
  
    G4int fsurface_type;
    G4double fwrapping_refl;
    G4double fSigmaAlpha,fSS,fSL,fBS;
    G4int ffrontsurface_type;
    G4double ffrontwrapping_refl;
    G4double ffrontSigmaAlpha,ffrontSS,ffrontSL,ffrontBS;
    G4int flateralsurface_type;
    G4double flateralwrapping_refl;
    G4double flateralSigmaAlpha,flateralSS,flateralSL,flateralBS;
    G4int fbacksurface_type;
    G4double fbackwrapping_refl;
    G4double fbackSigmaAlpha,fbackSS,fbackSL,fbackBS;
    G4int fWindowAirsurface_type;
    G4double fWindowAirwrapping_refl;
    G4double fWindowAirSigmaAlpha,fWindowAirSS,fWindowAirSL,fWindowAirBS;


  
  
    G4double fExpHall_x;
    G4double fExpHall_y;
    G4double fExpHall_z;

    G4double fCryst_x;
    G4double fCryst_y;
    G4double fCryst_z;

    G4double fSiPM_x;
    G4double fSiPM_y;
    G4double fSiPM_z;
    
    std::vector<double> fSiPM_x_pos;
    std::vector<double> fSiPM_y_pos;
    std::vector<double> fGlue_x_pos;
    std::vector<double> fGlue_y_pos;
    std::vector<double> fSiPM_window_x_pos;
    std::vector<double> fSiPM_window_y_pos;

    G4double fGlue_x;
    G4double fGlue_y;
    G4double fGlue_z;

    G4double fSiPM_window_x;
    G4double fSiPM_window_y;
    G4double fSiPM_window_z;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
