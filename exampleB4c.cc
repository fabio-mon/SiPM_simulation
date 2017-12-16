
#include "OpNovicePhysicsList.hh"
#include "ConfigFile.hh"
#include "OpNoviceDetectorConstruction_1.hh"
#include "OpNoviceDetectorConstruction_longtile.hh"
#include "OpNoviceDetectorConstruction_1_tilted.hh"
#include "OpNoviceDetectorConstruction_4.hh"
#include "OpNoviceDetectorConstruction_5.hh"

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
/*
  for ( G4int i=1; i<argc; i=i+2 ) 
  {
     if      ( G4String(argv[i]) == "-m" ) macro   = argv[i+1];
     //else if ( G4String(argv[i]) == "-u" ) session = argv[i+1];
     else 
        if ( G4String(argv[i]) == "-r" )
           myseed  = atoi(argv[i+1]);
        else 
           if (G4String(argv[i]) == "-tilt")
              tilt_angle =  G4UIcommand::ConvertToDouble(argv[i+1]);
           else
              if ( G4String(argv[i]) == "-geom" )
              {
                 geom = G4UIcommand::ConvertToInt(argv[i+1]);
                 if(argc==i+4)
                 {
                    Crystal_length = G4UIcommand::ConvertToDouble(argv[i+2]);
                    SiPM_length = G4UIcommand::ConvertToDouble(argv[i+3]);
                    break;
                 }
                 else 
                 {
                   G4cout<<argc<<"\t"<<i+4<<G4endl;
                   PrintUsage();
                   return 1;
                 }
              }
           
     
#ifdef G4MULTITHREADED
     else if ( G4String(argv[i]) == "-t" ) {
                    nThreads = G4UIcommand::ConvertToInt(argv[i+1]);
    }
#endif
            else
            {
               PrintUsage();
               return 1;
            }
  }
*/
  std::string configFileName = argv[1];
  ConfigFile config(configFileName);
/*  if(argc==4)
  {
    if(G4String(argv[2])=="-m") 
      macro = argv[3];
    else
    {
      PrintUsage();
      return 1;
    }
  }
*/
/*
  G4String CutOption="";
  if (config.keyExists("CutOption"))
     CutOption = config.read<std::string> ("CutOption");  

  G4String PDEoption="";
  if (config.keyExists("PDEoption"))
     PDEoption = config.read<std::string> ("PDEoption");
*/

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

  G4int NSiPM=1;
  if (config.keyExists("NSiPM"))
     NSiPM= config.read<int> ("NSiPM");

  G4double tilt_angle=0.*CLHEP::deg;
  if (config.keyExists("tilt_angle"))
     tilt_angle= config.read<double> ("tilt_angle") * CLHEP::deg;

/*

  G4double Crystal_x=12.*mm;
  if (config.keyExists("Crystal_x"))
     Crystal_x= config.read<double> ("Crystal_x") * mm;

  G4double Crystal_y=12.*mm;
  if (config.keyExists("Crystal_y"))
     Crystal_y= config.read<double> ("Crystal_y") * mm;

  G4double Crystal_z=3.*mm;
  if (config.keyExists("Crystal_z"))
     Crystal_z= config.read<double> ("Crystal_z") * mm;

  G4double SiPM_x=5.*mm;
  if (config.keyExists("SiPM_x"))
     SiPM_x= config.read<double> ("SiPM_x") * mm;

  G4double SiPM_y=5.*mm;
  if (config.keyExists("SiPM_y"))
     SiPM_y= config.read<double> ("SiPM_y") * mm;

  G4double SiPM_z=0.8*mm;
  if (config.keyExists("SiPM_z"))
     SiPM_z= config.read<double> ("SiPM_z") * mm;

  G4int surface_type=0;
  if (config.keyExists("surface_type"))
     surface_type= config.read<int> ("surface_type");

  G4double wrapping_refl=97.;
  if (config.keyExists("wrapping_refl"))
     wrapping_refl= config.read<int> ("wrapping_refl");

  G4double tilt_angle=0.*deg;
  if (config.keyExists("tilt_angle"))
     tilt_angle= config.read<double> ("tilt_angle") * deg;


  G4double Crystal_length = std::max( Crystal_x , std::max(Crystal_y,Crystal_z) );
  G4double SiPM_length = std::max( SiPM_x , std::max(SiPM_y,SiPM_z) );
*/
  switch (NSiPM)
  {
     case 1:
/*
        if(Crystal_length>50.*CLHEP::mm || SiPM_length>50.*CLHEP::mm)
        {
           G4cout<<"ERROR: geometry not available"<<G4endl;
           exit (EXIT_FAILURE);
        }
*/
        if(tilt_angle == 0.)
           runManager-> SetUserInitialization(new OpNoviceDetectorConstruction_1(configFileName));
        else
           runManager-> SetUserInitialization(new OpNoviceDetectorConstruction_1_tilted(configFileName));
        break;

  case 2:
/*
     if(Crystal_length>50.*CLHEP::mm || SiPM_length>50.*CLHEP::mm)
     {
        G4cout<<"ERROR: geometry not available"<<G4endl;
        exit (EXIT_FAILURE);
     }
*/
     runManager-> SetUserInitialization(new OpNoviceDetectorConstruction_longtile(configFileName));
     break;

  case 4:
/*
     if(Crystal_length>50.*CLHEP::mm || 2*SiPM_length>Crystal_length)
     {
        G4cout<<"ERROR: geometry not available"<<G4endl;
        exit (EXIT_FAILURE);
     }
*/
     runManager-> SetUserInitialization(new OpNoviceDetectorConstruction_4(configFileName));
     break;

  case 5:
/*
     if(Crystal_length>50.*CLHEP::mm || 3*SiPM_length>Crystal_length)
     {
        G4cout<<"ERROR: geometry not available"<<G4endl;
        exit (EXIT_FAILURE);
     }
*/
     runManager-> SetUserInitialization(new OpNoviceDetectorConstruction_5(configFileName));
     break;

  default:
     if(tilt_angle == 0.)
        runManager-> SetUserInitialization(new OpNoviceDetectorConstruction_1(configFileName));
     else
        runManager-> SetUserInitialization(new OpNoviceDetectorConstruction_1_tilted(configFileName));
  }
 /* 
  G4cout<<"\nCrystal: "<< Crystal_x <<"x"<<Crystal_y<<"x"<<Crystal_z<<"mm^3\n"
        <<"SiPM: " << SiPM_x <<"x"<<SiPM_y<<"x"<<SiPM_z<<"mm^3"<<G4endl;
*/
  //G4VModularPhysicsList* physicsList = new FTFP_BERT;
  runManager->SetUserInitialization(new OpNovicePhysicsList());
/*
  G4String filename = "diffuse98";
  filename += "_Tile" 
           + G4UIcommand::ConvertToString(Crystal_length) + "x"
           + G4UIcommand::ConvertToString(Crystal_length) + "x4_"
           + G4UIcommand::ConvertToString(geom) + "SiPM" 
           + G4UIcommand::ConvertToString(SiPM_length) + "x"
           + G4UIcommand::ConvertToString(SiPM_length) + "_tilt"
           + G4UIcommand::ConvertToString(tilt_angle) + "_PDE"
           + PDEoption;
*/
  B4cActionInitialization* actionInitialization = new B4cActionInitialization(configFileName);
  runManager->SetUserInitialization(actionInitialization);
  
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

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
    // interactive mode : define UI session
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    if (ui->IsGUI()) {
      UImanager->ApplyCommand("/control/execute gui.mac");
    }
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !

  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
