
#ifndef SiPM_SD_h
#define SiPM_SD_h 1

#include "G4VSensitiveDetector.hh"

#include "B4cCalorHit.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent;

/// Calorimeter sensitive detector class
///
/// In Initialize(), it creates one hit for each calorimeter layer and one more
/// hit for accounting the total quantities in all layers.
///
/// The values are accounted in hits in ProcessHits() function which is called
/// by Geant4 kernel at each step.

class SiPM_SD : public G4VSensitiveDetector
{
  public:
    SiPM_SD(const G4String& name, const G4String& hitsCollectionName, G4int nofCells, G4String PDEoption, G4double weight);
    virtual ~SiPM_SD();
  
    // methods from base class
    virtual void   Initialize(G4HCofThisEvent* hitCollection);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual G4bool ProcessHits_manually(const G4Step* step, G4TouchableHistory* history);
    virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

  private:
    B4cCalorHitsCollection* fHitsCollection;
    G4int     fNofCells;
    G4String fPDEoption;
    G4double fweight;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

