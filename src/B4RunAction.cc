
#include "B4RunAction.hh"
#include "B4Analysis.hh"
#include "ConfigFile.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4RunAction::B4RunAction(B4cEventAction *eventAction, const std::string& configFileName)
 : G4UserRunAction(),
  fEventAction(eventAction),
  fconfigFileName(configFileName)
{ 
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in B4Analysis.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories 
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFirstHistoId(1);

  // Creating ntuple
  //
  analysisManager->CreateNtuple("ntu", "Photon properties of each event");
  analysisManager->CreateNtupleIColumn("eventNb");
  analysisManager->CreateNtupleFColumn("mu_x_hit");
  analysisManager->CreateNtupleFColumn("mu_y_hit");
  analysisManager->CreateNtupleFColumn("mu_Edep");
  analysisManager->CreateNtupleFColumn("mu_TrackLength");
  analysisManager->CreateNtupleIColumn("Ph_tot_Nb");  
  analysisManager->CreateNtupleIColumn("Ph_detected_Nb");
  analysisManager->CreateNtupleIColumn("Ph_lost_Nb");
  analysisManager->CreateNtupleIColumn("Ph_labs_Nb");
  analysisManager->CreateNtupleFColumn("Ph_time", fEventAction->GetTimeVec());
  //analysisManager->CreateNtupleFColumn("Ph_energy", fEventAction->GetEnergyVec());
  analysisManager->CreateNtupleIColumn("Ph_creator_proc", fEventAction->GetCreatorProcVec());
  analysisManager->CreateNtupleIColumn("Ph_SiPM_number", fEventAction->GetSiPM_numberVec());
  //analysisManager->CreateNtupleFColumn("Ph_x_hit", fEventAction->GetPhXVec());
  //analysisManager->CreateNtupleFColumn("Ph_y_hit", fEventAction->GetPhYVec());
  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4RunAction::~B4RunAction()
{
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4RunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  //Create the output filename
  ConfigFile config (fconfigFileName) ;
  G4String geometry="planar";
  if (config.keyExists("geometry"))
     geometry = config.read<std::string> ("geometry");

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

  G4double wrapping_refl=0.97;
  if (config.keyExists("wrapping_refl"))
     wrapping_refl= config.read<double> ("wrapping_refl");

  G4double tilt_angle=0.*CLHEP::deg;
  if (config.keyExists("tilt_angle"))
     tilt_angle= config.read<double> ("tilt_angle") * CLHEP::deg; 

  G4String PDEoption="";
  if (config.keyExists("PDEoption"))
     PDEoption = config.read<std::string> ("PDEoption");

  G4String source_distr="";
  if (config.keyExists("SourceDistribution"))
     source_distr = config.read<std::string> ("SourceDistribution");

  G4String path="";
  if (config.keyExists("output_path"))
     path = config.read<std::string> ("output_path");

  G4int Seed=0;
  if (config.keyExists("Seed"))
     Seed = config.read<int> ("Seed");

  G4String filename =/* path +*/ "surface";
  filename += G4UIcommand::ConvertToString(surface_type); 
  if(surface_type!=0)
    filename += "_refl"+ G4UIcommand::ConvertToString(wrapping_refl*100);
  filename += "_Tile" 
           + G4UIcommand::ConvertToString(Crystal_x) + "x"
           + G4UIcommand::ConvertToString(Crystal_y) + "x"
           + G4UIcommand::ConvertToString(Crystal_z) + "x_"
           + geometry + "_geometry_" 
           + G4UIcommand::ConvertToString(SiPM_x) + "x"
           + G4UIcommand::ConvertToString(SiPM_y) + "x"
           + G4UIcommand::ConvertToString(SiPM_z) + "_tilt"
           + G4UIcommand::ConvertToString(tilt_angle) + "_PDE"
           + PDEoption + "_source"
           + source_distr + "_seed"
           + G4UIcommand::ConvertToString(Seed) + ".root";

  // Open an output file
  //
  analysisManager->OpenFile(filename);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // print histogram statistics
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // save histograms & ntuple
  //
  analysisManager->Write();
  analysisManager->CloseFile();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
