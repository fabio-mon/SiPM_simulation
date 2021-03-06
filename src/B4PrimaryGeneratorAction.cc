
#include "B4PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "ConfigFile.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4PrimaryGeneratorAction::B4PrimaryGeneratorAction( std::string configFileName )
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0),
   fconfigFileName(configFileName)
{
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

  //get source parameters from config file
  ConfigFile config (configFileName) ;
    
  std::string Particle;
  if (config.keyExists("Particle"))
     Particle = config.read<std::string> ("Particle");
  else
     Particle = "pi+";
  G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle(Particle.c_str());
  fParticleGun->SetParticleDefinition(particleDefinition);

  G4double Energy;
  if (config.keyExists("Energy"))
     Energy = config.read<double> ("Energy") * GeV;
  else
     Energy = 100.*GeV;
  fParticleGun->SetParticleEnergy(Energy);

  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));

  if (config.keyExists("SourceDistribution"))
     fsource_distr = config.read<std::string> ("SourceDistribution");
  else
     fsource_distr = "Uniform";

  if(fsource_distr == "Point" || fsource_distr == "POINT" || fsource_distr == "point")
  {
     if (config.keyExists("X_SourcePosition"))
        config.readIntoVect(fX_SourcePos,"X_SourcePosition");
     else
        fX_SourcePos = std::vector<double> (1,0);

     if (config.keyExists("Y_SourcePosition"))
        config.readIntoVect(fY_SourcePos,"Y_SourcePosition");
     else
        fY_SourcePos = std::vector<double> (1,0);

     if(fX_SourcePos.size() != fY_SourcePos.size())
     {
        G4cerr << "ERROR: different size of X_SourcePosition and Y_SourcePosition" << G4endl;
        exit(EXIT_FAILURE);
     }

  }

  if (config.keyExists("Nevents"))
     fNevents = config.read<int> ("Nevents");
  else
     fNevents = 0;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4PrimaryGeneratorAction::~B4PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume
  // from G4LogicalVolumeStore
  //
  G4int eventNb = anEvent->GetEventID();
  G4double CrystalHalfXLength = 0;
  G4double CrystalHalfYLength = 0;
  G4LogicalVolume* crystalLV = G4LogicalVolumeStore::GetInstance()->GetVolume("Crystal");
  G4Box* crystalBox = 0;
  G4double Xmin;
  G4double Xmax;
  G4double Ymin;
  G4double Ymax;
  G4double X,Y,Z;
    
  if(fsource_distr == "Uniform" || fsource_distr == "uniform" || fsource_distr == "UNIFORM")
  {
     if ( crystalLV) 
        crystalBox = dynamic_cast< G4Box*>(crystalLV->GetSolid()); 
     if ( crystalBox ) 
     {
        CrystalHalfXLength = crystalBox->GetXHalfLength();  
        CrystalHalfYLength = crystalBox->GetYHalfLength();  
     }
     else  
     {
        G4ExceptionDescription msg;
        msg << "Crystal volume of box not found." << G4endl;
        msg << "Perhaps you have changed geometry." << G4endl;
        //msg << "The gun will be place in the center.";
        //G4Exception("B4PrimaryGeneratorAction::GeneratePrimaries()","MyCode0002", JustWarning, msg);
     } 
     // Set gun position
     Xmin = -CrystalHalfXLength;
     Xmax = +CrystalHalfXLength;
     Ymin = -CrystalHalfYLength;
     Ymax = +CrystalHalfYLength;
     X = Xmin + G4UniformRand()*(Xmax-Xmin);
     Y = Ymin + G4UniformRand()*(Ymax-Ymin);
     Z = -0.05*m;
  }
  else
  {
     if(fsource_distr == "Point" || fsource_distr == "POINT" || fsource_distr == "point")
     {
        G4double event_percentage = 1. * eventNb / fNevents;
        int vec_position = (int) (event_percentage * fX_SourcePos.size());
        X = fX_SourcePos.at(vec_position);
        Y = fY_SourcePos.at(vec_position);
        Z = -0.05*m;
     }
     else
     {
        G4cerr<<"ERROR: source distribution of type "<<fsource_distr.c_str()<<" does not exist"<<G4endl;
        exit(EXIT_FAILURE);
     }
  }
  fParticleGun ->SetParticlePosition(G4ThreeVector(X, Y, Z));
  fParticleGun->GeneratePrimaryVertex(anEvent);

//  fParticleGun ->SetParticlePosition(G4ThreeVector(0, Ymin + ( (Ymax-Ymin)/(9+1) + (Ymax-Ymin)/(9+1) * ((int)(G4UniformRand()*9.)) )*mm ,-0.05*m));
// G4cout<<"-----------------------       "<< Ymin <<" + "<< (Ymax-Ymin)/(3+1)<<" + "<< (Ymax-Ymin)/(3+1) <<"*"<< ((int)(G4UniformRand()*3.)) <<" = " << Ymin + (Ymax-Ymin)/(3+1) + (Ymax-Ymin)/(3+1) * ((int)(G4UniformRand()*3.)) <<G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

