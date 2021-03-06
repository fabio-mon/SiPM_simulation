
#ifndef B4PrimaryGeneratorAction_h
#define B4PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include <vector>

class G4ParticleGun;
class G4Event;

/// The primary generator action class with particle gum.
///
/// It defines a single particle which hits the calorimeter 
/// perpendicular to the input face. The type of the particle
/// can be changed via the G4 build-in commands of G4ParticleGun class 
/// (see the macros provided with this example).

class B4PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  B4PrimaryGeneratorAction(std::string configFileName);    
  virtual ~B4PrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event* event);
  
  // set methods
  void SetRandomFlag(G4bool value);
  std::vector<double>& get_fX_SourcePos()       { return fX_SourcePos; }
  std::vector<double>& get_fY_SourcePos()       { return fY_SourcePos; }

private:
  G4ParticleGun*  fParticleGun; // G4 particle gun
  std::string fconfigFileName;
  std::string fsource_distr;
  std::vector<double> fX_SourcePos;
  std::vector<double> fY_SourcePos;
  G4int fNevents;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


