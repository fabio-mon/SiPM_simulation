//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef OpNoviceDetectorConstruction_4_h
#define OpNoviceDetectorConstruction_4_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class OpNoviceDetectorConstruction_4 : public G4VUserDetectorConstruction
{
  public:
    OpNoviceDetectorConstruction_4(std::string configFileName);
    virtual ~OpNoviceDetectorConstruction_4();

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

#endif /*OpNoviceDetectorConstruction_4_h*/
