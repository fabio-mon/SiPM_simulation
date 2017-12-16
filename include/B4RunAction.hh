
#ifndef B4RunAction_h
#define B4RunAction_h 1

#include "G4UserRunAction.hh"
#include "B4cEventAction.hh"
#include "globals.hh"

class G4Run;
class B4RunAction : public G4UserRunAction
{
  public:
    B4RunAction(B4cEventAction* eventAction, const std::string& configFileName);
    virtual ~B4RunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);
 
  private:
    B4cEventAction* fEventAction;
    std::string fconfigFileName;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

