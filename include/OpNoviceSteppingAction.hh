
#ifndef OpNoviceSteppingAction_h
#define OpNoviceSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include "B4cEventAction.hh"

/// Stepping action class
/// 

class OpNoviceSteppingAction : public G4UserSteppingAction
{
  public:
    OpNoviceSteppingAction(B4cEventAction* eventAction,G4String CutOption);
    virtual ~OpNoviceSteppingAction();

    // method from the base class
    virtual void UserSteppingAction(const G4Step*);

  private:
    B4cEventAction* fEventAction;
    G4int fEventNumber;
    G4String fCutOption;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
