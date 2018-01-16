
#ifndef OpNoviceStackingAction_H
#define OpNoviceStackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"
#include "B4cEventAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class OpNoviceStackingAction : public G4UserStackingAction
{
  public:
    OpNoviceStackingAction(B4cEventAction* eventAction,std::string configFileName);
    virtual ~OpNoviceStackingAction();

  public:
    virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
    //virtual void NewStage();
    //virtual void PrepareNewEvent();

  private:
    G4String fCutOption;
    G4double fTimeCut;
    B4cEventAction *fEventAction;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
