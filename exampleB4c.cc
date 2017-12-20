
#include "OpNovicePhysicsList.hh"
#include "ConfigFile.hh"
#include "DetectorConstruction_planar.hh"
#include "OpNoviceDetectorConstruction_longtile.hh"

#include "B4cActionInitialization.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "FTFP_BERT.hh"

#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " ./exampleB4c [filename.cfg]" << G4endl;
    G4cerr << "   note: cfg is mandatory (default template.cfg)"
           << G4endl;
    
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Evaluate arguments
  //
  if ( argc == 1 || argc >= 3) 
  {
    PrintUsage();
    return 1;
  }

  //G4String macro;
  G4String session;

#ifdef G4MULTITHREADED
  G4int nThreads = 4;
#endif

  std::string configFileName = argv[1];
  ConfigFile config(configFileName);

  // Detect interactive mode (if no macro provided) and define UI session
  //
  G4int Nevents=0;
  if (config.keyExists("Nevents"))
     Nevents = config.read<int> ("Nevents");  

  G4long myseed = 0;//time(NULL);
  if (config.keyExists("Seed"))
     myseed= config.read<int> ("Seed");  

  G4UIExecutive* ui = 0;
  if ( Nevents<=0 /*! macro.size()*/ ) 
  {
    ui = new G4UIExecutive(argc, argv);
  }

  // Choose the Random engine
  //
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  G4Random::setTheSeed(myseed);
  
  // Construct the default run manager
  //
#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
  if ( nThreads > 0 )
  { 
    runManager->SetNumberOfThreads(nThreads);
  }  
#else
  G4RunManager * runManager = new G4RunManager;
#endif

  // Set the geometry
  //

  G4String geometry;
  if (config.keyExists("geometry"))
    geometry = config.read<std::string> ("geometry");
  else
  {
    G4cerr<<"WARNING: geometry not set, default geometry: planar"<<G4endl;
    geometry = "planar";
  }

  if(geometry == "planar" || geometry == "Planar" || geometry == "PLANAR")
    runManager-> SetUserInitialization(new DetectorConstruction_planar(configFileName));
  else
    if(geometry == "longtile" || geometry == "Longtile" || geometry == "LONGTILE")
      runManager-> SetUserInitialization(new OpNoviceDetectorConstruction_longtile(configFileName));
    else
    {
      G4cerr<<"ERROR: geometry "<<geometry<<" is NOT VALID!"<<G4endl;
      exit(EXIT_FAILURE);  
    }

//Initialize physics list
  //G4VModularPhysicsList* physicsList = new FTFP_BERT;
  runManager->SetUserInitialization(new OpNovicePhysicsList());

//Initialize ActionInizializiation which initializes the other required classes
  B4cActionInitialization* actionInitialization = new B4cActionInitialization(configFileName);
  runManager->SetUserInitialization(actionInitialization);
  


  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if ( Nevents > 0 /*macro.size()*/ ) 
  {
    // batch mode
    //G4String command = "/control/execute ";
    //UImanager->ApplyCommand(command+macro);
    runManager->Initialize();
    runManager->BeamOn(Nevents);
  }
  else  
  {  
    // Initialize visualization
    G4VisManager* visManager = new G4VisExecutive;
    // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
    // G4VisManager* visManager = new G4VisExecutive("Quiet");
    visManager->Initialize();
    // interactive mode : define UI session
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    if (ui->IsGUI()) 
    {
      UImanager->ApplyCommand("/control/execute gui.mac");
    }
    ui->SessionStart();
    delete ui;
    delete visManager;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !

  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
