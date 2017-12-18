#ifndef MyMaterials_hh
#define MyMaterials_hh

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "ConfigFile.hh"



class MyMaterials
{
private:
  
public:
  MyMaterials();
  ~MyMaterials();
  
  static G4Material* Vacuum();
  static G4Material* Air();
  static G4Material* LYSO();
  static G4Material* Silicon();
  static G4Material* Glue();
  static G4Material* OpticWindow();
  static G4Material* PlasticScintillator();

  static G4double fromNmToEv(G4double wavelength);
  static G4double fromEvToNm(G4double energy);

};

#endif
