
#include "DetectorConstruction_planar.hh"
#include "SiPM_SD.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "ConfigFile.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* DetectorConstruction_planar::fMagFieldMessenger = 0; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction_planar::DetectorConstruction_planar(std::string configFileName)
 : G4VUserDetectorConstruction(),
   fCheckOverlaps(true)
{
  fExpHall_x = fExpHall_y = fExpHall_z = 0.10*m;

  //get geom parameters from config file
  ConfigFile config (configFileName) ;

  //Crystal
  if (config.keyExists("Crystal_x"))
     fCryst_x= config.read<double> ("Crystal_x") * mm/2;
  else
     fCryst_x=12.*mm/2;

  if (config.keyExists("Crystal_y"))
     fCryst_y= config.read<double> ("Crystal_y") * mm/2;
  else
     fCryst_y=12.*mm/2;

  if (config.keyExists("Crystal_z"))
     fCryst_z= config.read<double> ("Crystal_z") * mm/2;
  else
     fCryst_z=12.*mm/2;

  //SiPM
  if (config.keyExists("SiPM_x"))
     fSiPM_x= config.read<double> ("SiPM_x") * mm/2;
  else
     fSiPM_x=5.*mm/2;

  if (config.keyExists("SiPM_y"))
     fSiPM_y= config.read<double> ("SiPM_y") * mm/2;
  else
     fSiPM_y=5.*mm/2;

  if (config.keyExists("SiPM_z"))
     fSiPM_z= config.read<double> ("SiPM_z") * mm/2;
  else
     fSiPM_z=0.8*mm/2;
  
  //Glue
  if (config.keyExists("Glue_x"))
     fGlue_x= config.read<double> ("Glue_x") * mm/2;
  else
     fGlue_x=fSiPM_x;

  if (config.keyExists("Glue_y"))
     fGlue_y= config.read<double> ("Glue_y") * mm/2;
  else
     fGlue_y=fSiPM_y;

  if (config.keyExists("Glue_z"))
     fGlue_z= config.read<double> ("Glue_z") * mm/2;
  else
     fGlue_z=0.05*mm;
  
  //SiPM optic window 
  if (config.keyExists("SiPM_window_x"))
     fSiPM_window_x= config.read<double> ("SiPM_window_x") * mm/2;
  else
     fSiPM_window_x=fSiPM_x;

  if (config.keyExists("SiPM_window_y"))
     fSiPM_window_y= config.read<double> ("SiPM_window_y") * mm/2;
  else
     fSiPM_window_y=fSiPM_y;

  if (config.keyExists("SiPM_window_z"))
     fSiPM_window_z= config.read<double> ("SiPM_window_z") * mm/2;
  else
     fSiPM_window_z= 0.15*mm;
  
  

  /*
  if (config.keyExists("GlueWindowEverywhere"))
     fGlueWindowEverywhere = config.read<bool> ("GlueWindowEverywhere");
  else
     fGlueWindowEverywhere = false;
  
  if(fGlueWindowEverywhere == false)
  {
    fGlue_x = fSiPM_x;
    fGlue_y = fSiPM_y;
    fGlue_z = 0.05*mm;

    fSiPM_window_x = fSiPM_x;
    fSiPM_window_y = fSiPM_y;
    fSiPM_window_z = 0.15*mm;
  }

  else
  {
    fGlue_x = fCryst_x;
    fGlue_y = fCryst_y;
    fGlue_z = 0.05*mm;

    fSiPM_window_x = fCryst_x;
    fSiPM_window_y = fCryst_y;
    fSiPM_window_z = 0.15*mm;
  }
  */
    
  //SiPM positions
  if (config.keyExists("SiPM_x_pos"))
     config.readIntoVect(fSiPM_x_pos,"SiPM_x_pos");
  else
     fSiPM_x_pos = std::vector<double> (1,0);

  if (config.keyExists("SiPM_y_pos"))
     config.readIntoVect(fSiPM_y_pos,"SiPM_y_pos");
  else
     fSiPM_y_pos = std::vector<double> (1,0);

  assert(fSiPM_x_pos.size() == fSiPM_y_pos.size() );

  //Glue positions
  if (config.keyExists("Glue_x_pos"))
     config.readIntoVect(fGlue_x_pos,"Glue_x_pos");
  else
    fGlue_x_pos = fSiPM_x_pos;

  if (config.keyExists("Glue_y_pos"))
     config.readIntoVect(fGlue_y_pos,"Glue_y_pos");
  else
    fGlue_y_pos = fSiPM_y_pos;

  assert(fGlue_x_pos.size() == fGlue_y_pos.size() );

  //optic window positions
  if (config.keyExists("SiPM_window_x_pos"))
     config.readIntoVect(fSiPM_window_x_pos,"SiPM_window_x_pos");
  else
    fSiPM_window_x_pos = fSiPM_x_pos;

  if (config.keyExists("SiPM_window_y_pos"))
     config.readIntoVect(fSiPM_window_y_pos,"SiPM_window_y_pos");
  else
    fSiPM_window_y_pos = fSiPM_y_pos;

  assert(fSiPM_window_x_pos.size() == fSiPM_window_y_pos.size() );

  //Tilt SiPM wrt to beam
  if (config.keyExists("tilt_angle"))
     ftilt_angle = config.read<double> ("tilt_angle")*deg;
  else
     ftilt_angle = 0.*deg;

  //PDE settings
  if (config.keyExists("PDEoption"))
     fPDEoption = config.read<std::string> ("PDEoption");
  else
     fPDEoption="";

  if (config.keyExists("PDEweight"))
     fweight = config.read<double> ("PDEweight");
  else
     fweight = 1.;

  //SURFACE in general
  if (config.keyExists("surface_type"))
     fsurface_type= config.read<int> ("surface_type");
  else
     fsurface_type=0;

  if (config.keyExists("wrapping_refl"))
     fwrapping_refl= config.read<double> ("wrapping_refl");
  else
     fwrapping_refl=0.97;

  if (config.keyExists("SigmaAlpha"))
     fSigmaAlpha = config.read<double> ("SigmaAlpha");
  else
     fSigmaAlpha=3.4*deg;

  if (config.keyExists("SpecularSpike") && config.keyExists("SpecularLobe") && config.keyExists("BackScatter"))
  {
     fSS = config.read<double> ("SpecularSpike");     
     fSL = config.read<double> ("SpecularLobe");
     fBS = config.read<double> ("BackScatter");
  }
  else
  {
     fSS=0.;
     fSL=1.;
     fBS=0.;
  }

  //FRONT SURFACE (opposite side wrt SiPM)
  if (config.keyExists("frontsurface_type"))
     ffrontsurface_type= config.read<int> ("frontsurface_type");
  else
     ffrontsurface_type=fsurface_type;

  if (config.keyExists("frontwrapping_refl"))
     ffrontwrapping_refl= config.read<double> ("frontwrapping_refl");
  else
     ffrontwrapping_refl=fwrapping_refl;

  if (config.keyExists("frontSigmaAlpha"))
     ffrontSigmaAlpha = config.read<double> ("frontSigmaAlpha");
  else
     ffrontSigmaAlpha=fSigmaAlpha;

  if (config.keyExists("frontSpecularSpike") && config.keyExists("frontSpecularLobe") && config.keyExists("frontBackScatter"))
  {
     ffrontSS = config.read<double> ("frontSpecularSpike");     
     ffrontSL = config.read<double> ("frontSpecularLobe");
     ffrontBS = config.read<double> ("frontBackScatter");
  }
  else
  {
     ffrontSS=fSS;
     ffrontSL=fSL;
     ffrontBS=fBS;
  }


  //LATERAL SURFACES
  if (config.keyExists("lateralsurface_type"))
     flateralsurface_type= config.read<int> ("lateralsurface_type");
  else
     flateralsurface_type=fsurface_type;

  if (config.keyExists("lateralwrapping_refl"))
     flateralwrapping_refl= config.read<double> ("lateralwrapping_refl");
  else
     flateralwrapping_refl=fwrapping_refl;

  if (config.keyExists("lateralSigmaAlpha"))
     flateralSigmaAlpha = config.read<double> ("lateralSigmaAlpha");
  else
     flateralSigmaAlpha=fSigmaAlpha;

  if (config.keyExists("lateralSpecularSpike") && config.keyExists("lateralSpecularLobe") && config.keyExists("lateralBackScatter"))
  {
     flateralSS = config.read<double> ("lateralSpecularSpike");     
     flateralSL = config.read<double> ("lateralSpecularLobe");
     flateralBS = config.read<double> ("lateralBackScatter");
  }
  else
  {
     flateralSS=fSS;
     flateralSL=fSL;
     flateralBS=fBS;
  }

  //BACK SURFACE (same side of SiPM)
  if (config.keyExists("backsurface_type"))
     fbacksurface_type= config.read<int> ("backsurface_type");
  else
     fbacksurface_type=fsurface_type;

  if (config.keyExists("backwrapping_refl"))
     fbackwrapping_refl= config.read<double> ("backwrapping_refl");
  else
     fbackwrapping_refl=fwrapping_refl;

  if (config.keyExists("backSigmaAlpha"))
     fbackSigmaAlpha = config.read<double> ("backSigmaAlpha");
  else
     fbackSigmaAlpha=fSigmaAlpha;

  if (config.keyExists("backSpecularSpike") && config.keyExists("backSpecularLobe") && config.keyExists("backBackScatter"))
  {
     fbackSS = config.read<double> ("backSpecularSpike");     
     fbackSL = config.read<double> ("backSpecularLobe");
     fbackBS = config.read<double> ("backBackScatter");
  }
  else
  {
     fbackSS=fSS;
     fbackSL=fSL;
     fbackBS=fBS;
  }

  //SURFACE opticwindow-air
  if (config.keyExists("WindowAirsurface_type"))
     fWindowAirsurface_type= config.read<int> ("WindowAirsurface_type");
  else
     fWindowAirsurface_type=0;

  if (config.keyExists("WindowAirwrapping_refl"))
     fWindowAirwrapping_refl= config.read<double> ("WindowAirwrapping_refl");
  else
     fWindowAirwrapping_refl=0.97;

  if (config.keyExists("WindowAirSigmaAlpha"))
     fWindowAirSigmaAlpha = config.read<double> ("WindowAirSigmaAlpha");
  else
     fWindowAirSigmaAlpha=0.;

  if (config.keyExists("WindowAirSpecularSpike") && config.keyExists("WindowAirSpecularLobe") && config.keyExists("WindowAirbackScatter"))
  {
     fWindowAirSS = config.read<double> ("WindowAirSpecularSpike");     
     fWindowAirSL = config.read<double> ("WindowAirSpecularLobe");
     fWindowAirBS = config.read<double> ("WindowAirbackScatter");
  }
  else
  {
     fWindowAirSS=fSS;
     fWindowAirSL=fSL;
     fWindowAirBS=fBS;
  }
  
  std::cout<<"initialization ok"<<std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction_planar::~DetectorConstruction_planar()
{ 
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction_planar::Construct()
{

  // -----------------------> Get materials

  G4Material* Air = MyMaterials::Air () ;
  G4Material* Vacuum = MyMaterials::Vacuum () ;
  G4Material *LYSO = MyMaterials::LYSO () ;
  G4Material *SiPM_mat = MyMaterials::Silicon () ;
  G4Material *Glue_mat = MyMaterials::Glue () ;
  G4Material *Window_mat = MyMaterials::OpticWindow () ;
  G4Material *PlasticScint_mat = MyMaterials::PlasticScintillator () ;
  
  if (! Vacuum || !  LYSO || !  Glue_mat || !  SiPM_mat || ! PlasticScint_mat || ! Window_mat) 
  {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("B4DetectorConstruction::DefineVolumes()", "MyCode0001", FatalException, msg);
  }  

  //------------------------> declare geom + logic + phys Volumes
  G4double posx, posy, posz;
  G4RotationMatrix *Rm = new G4RotationMatrix();
  Rm->rotateX(ftilt_angle);

// The experimental Hall

  G4Box* expHall_box = new G4Box("World",fExpHall_x,fExpHall_y,fExpHall_z);
  G4LogicalVolume* expHall_log = new G4LogicalVolume(expHall_box,Vacuum,"World",0,0,0);
  G4VPhysicalVolume* expHall_phys = new G4PVPlacement(0,G4ThreeVector(),expHall_log,"World",0,false,0);

// Small spaces around the crystal useful to define different surfaces around the crystal
  //Lateral surfaces
  G4Box* world_lateral_x_crystal_box = new G4Box("Lateral_world_x",0.2*mm,fCryst_y,fCryst_z);
  G4Box* world_lateral_y_crystal_box = new G4Box("Lateral_world_y",fCryst_x,0.2*mm,fCryst_z);
  G4Box* world_lateral_z_crystal_box = new G4Box("Lateral_world_z",fCryst_x,fCryst_y,0.2*mm);

  G4LogicalVolume* world_lateral_x_crystal_log = new G4LogicalVolume(world_lateral_x_crystal_box,Vacuum,"Lateral_world_x",0,0,0);
  G4LogicalVolume* world_lateral_y_crystal_log = new G4LogicalVolume(world_lateral_y_crystal_box,Vacuum,"Lateral_world_y",0,0,0);
  G4LogicalVolume* world_lateral_z_crystal_log = new G4LogicalVolume(world_lateral_z_crystal_box,Vacuum,"Lateral_world_z",0,0,0);

  G4ThreeVector world_lateral_x1_crystal_vector(-fCryst_x-0.2*mm,0.,0.);
  world_lateral_x1_crystal_vector.rotateX(-ftilt_angle);
  G4ThreeVector world_lateral_x2_crystal_vector(+fCryst_x+0.2*mm,0.,0.);
  world_lateral_x2_crystal_vector.rotateX(-ftilt_angle);
  G4ThreeVector world_lateral_y1_crystal_vector(0.,-fCryst_x-0.2,0.);
  world_lateral_y1_crystal_vector.rotateX(-ftilt_angle);   
  G4ThreeVector world_lateral_y2_crystal_vector(0.,+fCryst_x+0.2,0.);
  world_lateral_y2_crystal_vector.rotateX(-ftilt_angle);   
  G4ThreeVector world_lateral_z_crystal_vector(0.,0.,-fCryst_z-0.2); 
  world_lateral_z_crystal_vector.rotateX(-ftilt_angle);

  G4VPhysicalVolume* world_lateral_x1_crystal_phys = 
new G4PVPlacement(Rm , world_lateral_x1_crystal_vector , world_lateral_x_crystal_log,"Lateral_world",expHall_log,false,0);
  G4VPhysicalVolume* world_lateral_x2_crystal_phys = 
new G4PVPlacement(Rm , world_lateral_x2_crystal_vector , world_lateral_x_crystal_log,"Lateral_world",expHall_log,false,0);
  G4VPhysicalVolume* world_lateral_y1_crystal_phys = 
new G4PVPlacement(Rm , world_lateral_y1_crystal_vector , world_lateral_y_crystal_log,"Lateral_world",expHall_log,false,0);
  G4VPhysicalVolume* world_lateral_y2_crystal_phys = 
new G4PVPlacement(Rm , world_lateral_y2_crystal_vector , world_lateral_y_crystal_log,"Lateral_world",expHall_log,false,0);
  G4VPhysicalVolume* world_lateral_z_crystal_phys = 
new G4PVPlacement(Rm , world_lateral_z_crystal_vector , world_lateral_z_crystal_log,"Lateral_world",expHall_log,false,0);
   
// The LYSO crystal

  G4ThreeVector crystal_vector;
  crystal_vector.rotateX(-ftilt_angle);
  G4Box* crystal_box = new G4Box("Crystal",fCryst_x,fCryst_y,fCryst_z);
  G4LogicalVolume* crystal_log = new G4LogicalVolume(crystal_box,LYSO,"Crystal",0,0,0);
  G4VPhysicalVolume* crystal_phys = new G4PVPlacement(Rm, crystal_vector,crystal_log,"Crystal", expHall_log,false,0);
                       
//loop to place the SiPM+glue+window units

  // Optical glue             
  G4ThreeVector Glue_vector;
  G4Box* Glue_box = new G4Box("Optic_Glue",fGlue_x,fGlue_y,fGlue_z);
  G4LogicalVolume* Glue_log = new G4LogicalVolume(Glue_box,Glue_mat,"Optic_Glue",0,0,0);
  std::vector<G4VPhysicalVolume*> Glue_phys;
  for (unsigned i=0; i < fGlue_x_pos.size(); i++) 
  {    
    Glue_vector = G4ThreeVector(fGlue_x_pos.at(i)*mm,fGlue_y_pos.at(i)*mm, posz = fCryst_z + fGlue_z );
    Glue_vector.rotateX(-ftilt_angle);              
    Glue_phys.push_back( new G4PVPlacement(Rm,Glue_vector,Glue_log, "Optic_Glue",expHall_log,false,0) );
  }
  std::cout<<"glue placed"<<std::endl;
  /*
  if(fGlueWindowEverywhere==true)
  {
    Glue_vector = G4ThreeVector( 0. , 0. , posz = fCryst_z + fGlue_z );
    Glue_vector.rotateX(-ftilt_angle);              
    Glue_phys.push_back( new G4PVPlacement(Rm,Glue_vector,Glue_log, "Optic_Glue",expHall_log,false,0) );
  }
  */

  //SiPM optical grease
  G4ThreeVector Window_vector;
  G4Box* Window_box = new G4Box("Optic_window",fSiPM_window_x,fSiPM_window_y,fSiPM_window_z);
  G4LogicalVolume* Window_log = new G4LogicalVolume(Window_box,Window_mat,"Optic_window",0,0,0);
  std::vector<G4VPhysicalVolume*> Window_phys;
  for (unsigned i=0; i < fSiPM_window_x_pos.size(); i++) 
  {    
    //SiPM optical grease
    Window_vector = G4ThreeVector(fSiPM_window_x_pos.at(i)*mm,fSiPM_window_y_pos.at(i)*mm, posz = fCryst_z +2* fGlue_z +fSiPM_window_z);
    Window_vector.rotateX(-ftilt_angle);
    Window_phys.push_back( new G4PVPlacement(Rm,Window_vector,Window_log,"Optic_window", expHall_log,false,0) );
  }
  std::cout<<"window placed"<<std::endl;

  /*
  if(fGlueWindowEverywhere==true)
  {
    Window_vector = G4ThreeVector( 0. , 0. , posz = fCryst_z +2* fGlue_z +fSiPM_window_z);
    Window_vector.rotateX(-ftilt_angle);              
    Window_phys.push_back( new G4PVPlacement(Rm,Window_vector,Window_log, "Optic_window",expHall_log,false,0) );
  }
  */

  //SiPM
  G4ThreeVector SiPM_vector;
  G4Box* SiPM_box = new G4Box("SiPM",fSiPM_x,fSiPM_y,fSiPM_z);
  G4LogicalVolume* SiPM_log = new G4LogicalVolume(SiPM_box,SiPM_mat,"SiPM",0,0,0);
  std::vector<G4VPhysicalVolume*> SiPM_phys;
  for (unsigned i=0; i < fSiPM_x_pos.size(); i++) 
  {    
    //SiPM
    SiPM_vector = G4ThreeVector(fSiPM_x_pos.at(i)*mm,fSiPM_y_pos.at(i)*mm, posz = fCryst_z + 2*fGlue_z + 2*fSiPM_window_z + fSiPM_z );
    SiPM_vector.rotateX(-ftilt_angle);
    SiPM_phys.push_back( new G4PVPlacement(Rm,SiPM_vector,SiPM_log,"SiPM", expHall_log,false,0) );  
    /*
    if(fGlueWindowEverywhere==true) //in this case window and glue have been already placed
      continue;

    // Optical glue
    Glue_vector = G4ThreeVector(fSiPM_x_pos.at(i)*mm,fSiPM_y_pos.at(i)*mm, posz = fCryst_z + fGlue_z );
    Glue_vector.rotateX(-ftilt_angle);              
    Glue_phys.push_back( new G4PVPlacement(Rm,Glue_vector,Glue_log, "Optic_Glue",expHall_log,false,0) );

    //SiPM optical grease
    Window_vector = G4ThreeVector(fSiPM_x_pos.at(i)*mm,fSiPM_y_pos.at(i)*mm, posz = fCryst_z +2* fGlue_z +fSiPM_window_z);
    Window_vector.rotateX(-ftilt_angle);
    Window_phys.push_back( new G4PVPlacement(Rm,Window_vector,Window_log,"Optic_window", expHall_log,false,0) );
    */
  }


//------------------> Define surfaces

// LYSO-air front surface
  G4OpticalSurface* opLYSOFrontSurface = new G4OpticalSurface("LYSOSurface");
  new G4LogicalBorderSurface("LYSOSurface", crystal_phys,world_lateral_z_crystal_phys,opLYSOFrontSurface);
// LYSO-air lateral surface
  G4OpticalSurface* opLYSOLateralx1Surface = new G4OpticalSurface("LYSOSurface");
  new G4LogicalBorderSurface("LYSOSurface", crystal_phys,world_lateral_x1_crystal_phys,opLYSOLateralx1Surface);
  G4OpticalSurface* opLYSOLateralx2Surface = new G4OpticalSurface("LYSOSurface");
  new G4LogicalBorderSurface("LYSOSurface", crystal_phys,world_lateral_x2_crystal_phys,opLYSOLateralx2Surface);
  G4OpticalSurface* opLYSOLateraly1Surface = new G4OpticalSurface("LYSOSurface");
  new G4LogicalBorderSurface("LYSOSurface", crystal_phys,world_lateral_y1_crystal_phys,opLYSOLateraly1Surface);
  G4OpticalSurface* opLYSOLateraly2Surface = new G4OpticalSurface("LYSOSurface");
  new G4LogicalBorderSurface("LYSOSurface", crystal_phys,world_lateral_y2_crystal_phys,opLYSOLateraly2Surface);
  // LYSO-air back surface
  G4OpticalSurface* opLYSOBackSurface = new G4OpticalSurface("LYSOSurface");
  new G4LogicalBorderSurface("LYSOSurface", crystal_phys,expHall_phys,opLYSOBackSurface);

  SetSurfaceProperties(opLYSOFrontSurface,ffrontsurface_type,ffrontSigmaAlpha, ffrontwrapping_refl, ffrontSL, ffrontSS, ffrontBS);
  SetSurfaceProperties(opLYSOLateralx1Surface,flateralsurface_type,flateralSigmaAlpha, flateralwrapping_refl, flateralSL, flateralSS, flateralBS);
  SetSurfaceProperties(opLYSOLateralx2Surface,flateralsurface_type,flateralSigmaAlpha, flateralwrapping_refl, flateralSL, flateralSS, flateralBS);
  SetSurfaceProperties(opLYSOLateraly1Surface,flateralsurface_type,flateralSigmaAlpha, flateralwrapping_refl, flateralSL, flateralSS, flateralBS);
  SetSurfaceProperties(opLYSOLateraly2Surface,flateralsurface_type,flateralSigmaAlpha, flateralwrapping_refl, flateralSL, flateralSS, flateralBS);
  SetSurfaceProperties(opLYSOBackSurface,fbacksurface_type,fbackSigmaAlpha, fbackwrapping_refl, fbackSL, fbackSS, fbackBS);

  // Optic grease - air surface (trick to simulate SiPM reflective coating)
  vector<G4OpticalSurface*> opWindowAir;
  for(unsigned i=0;i<Window_phys.size();++i)
  {
    opWindowAir.push_back( new G4OpticalSurface("WindowAirSurface") );
    new G4LogicalBorderSurface("WindowAirSurface", Window_phys.at(i),expHall_phys,opWindowAir.at(i));
    SetSurfaceProperties(opWindowAir.at(i),fWindowAirsurface_type,fWindowAirSigmaAlpha, fWindowAirwrapping_refl, fWindowAirSL, fWindowAirSS, fWindowAirBS);
  }

//SiPM optical surface
//

G4MaterialPropertiesTable* myMPT5 = new G4MaterialPropertiesTable();
  G4double energySilicon[] =
            {4.96, 4.7692307692, 4.5925925926, 4.4285714286, 4.275862069, 4.1333333333, 4, 3.875, 3.7575757576, 3.6470588235, 3.5428571429, 3.4444444444, 3.3513513514, 3.2631578947, 3.1794871795, 3.1, 3.0243902439, 2.9523809524, 2.8837209302, 2.8181818182, 2.7555555556, 2.6956521739, 2.6382978723, 2.5833333333, 2.5306122449, 2.48, 2.431372549, 2.3846153846, 2.3396226415, 2.2962962963, 2.2545454545, 2.2142857143, 2.1754385965, 2.1379310345, 2.1016949153, 2.0666666667, 2.0327868852, 2, 1.9682539683, 1.9375, 1.9076923077, 1.8787878788, 1.8507462687, 1.8235294118, 1.7971014493, 1.7714285714, 1.7464788732, 1.7222222222, 1.698630137, 1.6756756757, 1.6533333333, 1.6315789474, 1.6103896104, 1.5897435897, 1.5696202532, 1.55, 1.5308641975, 1.512195122, 1.4939759036, 1.4761904762, 1.4588235294, 1.4418604651, 1.4252873563, 1.4090909091, 1.393258427, 1.3777777778, 1.3626373626, 1.347826087, 1.3333333333, 1.3191489362, 1.3052631579, 1.2916666667, 1.2783505155, 1.2653061224, 1.2525252525, 1.24};

  G4double Re_n[] = 
            {1.694, 1.8, 2.129, 3.052, 4.426, 5.055, 5.074, 5.102, 5.179, 5.293, 5.483, 6.014, 6.863, 6.548, 5.976, 5.587, 5.305, 5.091, 4.925, 4.793, 4.676, 4.577, 4.491, 4.416, 4.348, 4.293, 4.239, 4.192, 4.15, 4.11, 4.077, 4.044, 4.015, 3.986, 3.962, 3.939, 3.916, 3.895, 3.879, 3.861, 3.844, 3.83, 3.815, 3.8, 3.787, 3.774, 3.762, 3.751, 3.741, 3.732, 3.723, 3.714, 3.705, 3.696, 3.688, 3.681, 3.674, 3.668, 3.662, 3.656, 3.65, 3.644, 3.638, 3.632, 3.626, 3.62, 3.614, 3.608, 3.602, 3.597, 3.592, 3.587, 3.582, 3.578, 3.574, 3.57};

  G4double Im_n[] = 
            {3.666, 4.072, 4.69, 5.258, 5.16, 4.128, 3.559, 3.269, 3.085, 2.951, 2.904, 2.912, 2.051, 0.885, 0.465, 0.303, 0.22, 0.167, 0.134, 0.109, 0.091, 0.077, 0.064, 0.057, 0.05, 0.045, 0.039, 0.036, 0.033, 0.03, 0.028, 0.026, 0.024, 0.023, 0.021, 0.02, 0.018, 0.017, 0.016, 0.015, 0.015, 0.014, 0.013, 0.012, 0.011, 0.011, 0.011, 0.01, 0.009, 0.008, 0.008, 0.007, 0.007, 0.006, 0.006, 0.005, 0.005, 0.005, 0.004, 0.004, 0.004, 0.003, 0.003, 0.003, 0.002, 0.002, 0.002, 0.002, 0.002, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001};
  G4double Eff[] = 
            {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};

  G4double Refl[] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};/*{0.672612594, 0.7051739998, 0.7320895527, 0.7229564109, 0.6842353612, 0.623487589, 0.5904758352, 0.5741303379, 0.5656774122, 0.5617493182, 0.5653802759, 0.5829110024, 0.5842708013, 0.5465022924, 0.5109736438, 0.4860210277, 0.4668532597, 0.4515215805, 0.4391232489, 0.4289072653, 0.4195856995, 0.4114859504, 0.4042813942, 0.3978791808, 0.3919647328, 0.3871055313, 0.3822645215, 0.3779990919, 0.3741420137, 0.3704285222, 0.3673359361, 0.3642162212, 0.3614517486, 0.3586671346, 0.3563449752, 0.3541066718, 0.3518536143, 0.3497852413, 0.3482013386, 0.3464114799, 0.3447139284, 0.3433093164, 0.3417986254, 0.340281901, 0.3389624224, 0.3376390016, 0.3364132847, 0.3352856254, 0.33425759, 0.3333299988, 0.3324006727, 0.3314686517, 0.3305348304, 0.3295983594, 0.3287643932, 0.3280328678, 0.3273003059, 0.3266713015, 0.3260410062, 0.3254099758, 0.3247779306, 0.3241446496, 0.3235105692, 0.3228754694, 0.3222391901, 0.3216020461, 0.3209638781, 0.3203246844, 0.3196844635, 0.3191500638, 0.3186150452, 0.3180793106, 0.317542859, 0.3171131809, 0.316683043, 0.3162524448};*/
              
  assert(sizeof(energySilicon) == sizeof(Re_n) && sizeof(energySilicon) == sizeof(Im_n) && sizeof(energySilicon) == sizeof(Eff) );
  const G4int nEnSi = sizeof(energySilicon)/sizeof(G4double);
  //myMPT5->AddProperty("EFFICIENCY",energySilicon, Eff, nEnSi);
  myMPT5->AddProperty("REALRINDEX", energySilicon, Re_n, nEnSi);
  myMPT5->AddProperty("IMAGINARYRINDEX", energySilicon, Im_n, nEnSi);
  //myMPT5->AddProperty("REFLECTIVITY",energySilicon, Refl, nEnSi);
  G4cout << "Silicon G4MaterialPropertiesTable" << G4endl;

  G4OpticalSurface* SiPM_opsurf = new G4OpticalSurface("SiPM_opsurf",glisur,polished, dielectric_metal);
  SiPM_opsurf->SetMaterialPropertiesTable(myMPT5); 
  new G4LogicalSkinSurface("SiPM_surf",SiPM_log,SiPM_opsurf);
/*
  std::vector<G4OpticalSurface*> SiPM_opsurf;
  for(unsigned i=0;i<SiPM_log.size();i++)
  {
     SiPM_opsurf.push_back( new G4OpticalSurface("SiPM_opsurf",glisur,polished, dielectric_metal) );
     //new G4LogicalBorderSurface("SiPM_opsurf", Glue_phys,SiPM_phys,SiPM_opsurf);
     SiPM_opsurf.at(i)->SetMaterialPropertiesTable(myMPT5); 
     new G4LogicalSkinSurface("SiPM_surf",SiPM_log.at(i),SiPM_opsurf.at(i));
  }
*/
//Print parameters
  G4cout
    << G4endl 
    << "------------------------------------------------------------" << G4endl
    << "---> The Tile is " << fCryst_x*2 << " x " << fCryst_y*2 <<" x "<<fCryst_z*2<<" mm^3 "<< G4endl
    << "------------------------------------------------------------" << G4endl;
  
  //                                        
  // Visualization attributes
  //
  G4VisAttributes* simpleBoxVisAtt0= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt0->SetVisibility(false);
  expHall_log->SetVisAttributes(simpleBoxVisAtt0); //(G4VisAttributes::Invisible);

  G4VisAttributes* simpleBoxVisAtt1= new G4VisAttributes(G4Colour(1.0,0.0,0.0));//red
  simpleBoxVisAtt1->SetVisibility(true);
  crystal_log->SetVisAttributes(simpleBoxVisAtt1);

  G4VisAttributes* simpleBoxVisAtt2= new G4VisAttributes(G4Colour(0.0,1.0,0.0));//green
  simpleBoxVisAtt2->SetVisibility(true);
  SiPM_log->SetVisAttributes(simpleBoxVisAtt2);

  G4VisAttributes* simpleBoxVisAtt3= new G4VisAttributes(G4Colour(0.0,0.0,1.0));//blue
  simpleBoxVisAtt3->SetVisibility(true);
  Glue_log->SetVisAttributes(simpleBoxVisAtt3);

  G4VisAttributes* simpleBoxVisAtt4= new G4VisAttributes(G4Colour(0.8,0.4,0.0));//orange
  simpleBoxVisAtt4->SetVisibility(true);
  Window_log->SetVisAttributes(simpleBoxVisAtt4);

  //
  // Always return the physical World
  //
  return expHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......---------------------------------------
void DetectorConstruction_planar::SetSurfaceProperties
(G4OpticalSurface* opSurface, int surface_type, double SigmaAlpha, double refl, double SL, double SS, double BS)
{
  opSurface->SetType(dielectric_dielectric);
  opSurface->SetFinish(G4OpticalSurfaceFinish(surface_type));

  if(surface_type == 0 || surface_type == 1)
    opSurface->SetModel(glisur);
  else
    opSurface->SetModel(unified);

  const G4int NUM = 2;
  G4double XX[NUM] = {1.8*eV,3.4*eV} ; 
  G4double refractiveIndex[NUM] = {1.002, 1.002};
  G4double specularLobe[NUM]    = {SL, SL};
  G4double specularSpike[NUM]   = {SS, SS};
  G4double backScatter[NUM]     = {BS, BS};
  G4double REFLECT[NUM] = {refl,refl};
  G4MaterialPropertiesTable* myST1 = new G4MaterialPropertiesTable();

  if(surface_type >=2)
  {
     opSurface->SetSigmaAlpha(SigmaAlpha);
     myST1->AddProperty("SPECULARLOBECONSTANT",  XX, specularLobe,   NUM);
     myST1->AddProperty("SPECULARSPIKECONSTANT", XX, specularSpike,   NUM);
     myST1->AddProperty("BACKSCATTERCONSTANT",   XX, backScatter,     NUM);
  }

  if(surface_type == 1 || surface_type == 2 || surface_type == 4 || surface_type == 5)
     myST1->AddProperty("REFLECTIVITY", XX, REFLECT,NUM);

  if(surface_type == 2 || surface_type == 5)
     myST1->AddProperty("RINDEX", XX, refractiveIndex, NUM);

  if (surface_type >=1 && surface_type<=5)
    opSurface->SetMaterialPropertiesTable(myST1);
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......---------------------------------------

void DetectorConstruction_planar::ConstructSDandField()
{
  // G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // 
  // Sensitive detectors
  //
  SiPM_SD* SiPM = new SiPM_SD("SiPM_detector", "SiPMHitsCollection", 1,fPDEoption,fweight);
  SetSensitiveDetector("SiPM",SiPM);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
