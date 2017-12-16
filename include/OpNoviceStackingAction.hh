
#ifndef OpNoviceStackingAction_H
#define OpNoviceStackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"
#include "B4cEventAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class OpNoviceStackingAction : public G4UserStackingAction
{
  public:
    OpNoviceStackingAction(B4cEventAction* eventAction,G4String CutOption);
    virtual ~OpNoviceStackingAction();

  public:
    virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
    //virtual void NewStage();
    //virtual void PrepareNewEvent();

  private:
    G4String fCutOption;
    B4cEventAction *fEventAction;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif