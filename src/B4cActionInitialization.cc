
#include "B4cActionInitialization.hh"
#include "B4PrimaryGeneratorAction.hh"
#include "B4RunAction.hh"
#include "B4cEventAction.hh"
#include "OpNoviceSteppingAction.hh"
#include "OpNoviceStackingAction.hh"
#include "ConfigFile.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cActionInitialization::B4cActionInitialization(std::string configFileName)
 : G4VUserActionInitialization(),
   fconfigFileName(configFileName)
{
   ConfigFile config (fconfigFileName) ;
  if (config.keyExists("CutOption"))
     fCutOption = config.read<std::string> ("CutOption");  
  else
     fCutOption = "";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cActionInitialization::~B4cActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cActionInitialization::BuildForMaster() const
{
   B4cEventAction* eventAction = new B4cEventAction;
   SetUserAction(new B4RunAction(eventAction,fconfigFileName));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cActionInitialization::Build() const
{
   SetUserAction(new B4PrimaryGeneratorAction(fconfigFileName));
   B4cEventAction *eventAction = new B4cEventAction;
   SetUserAction(eventAction);
   SetUserAction(new B4RunAction(eventAction,fconfigFileName));
   SetUserAction(new OpNoviceSteppingAction(eventAction,fCutOption));
   SetUserAction(new OpNoviceStackingAction(eventAction,fCutOption));
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
