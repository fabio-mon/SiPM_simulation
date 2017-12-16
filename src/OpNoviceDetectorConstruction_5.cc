/*

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
the SiPM optical window has not been implemented yet

*/

#include "OpNoviceDetectorConstruction_5.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "SiPM_SD.hh"
#include "ConfigFile.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceDetectorConstruction_5::OpNoviceDetectorConstruction_5(std::string configFileName)
 : G4VUserDetectorConstruction()
{
   
  fExpHall_x = fExpHall_y = fExpHall_z = 0.10*m;

//get geom parameters from config file
  ConfigFile config (configFileName) ;
    
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

  fGlue_x = fSiPM_x;
  fGlue_y = fSiPM_y;
  fGlue_z = 0.05*mm;

  fSiPM_window_x = fSiPM_x;
  fSiPM_window_y = fSiPM_y;
  fSiPM_window_z = 0.15*mm;

  if (config.keyExists("surface_type"))
     fsurface_type= config.read<int> ("surface_type");
  else
     fsurface_type=0;

  if (config.keyExists("wrapping_refl"))
     fwrapping_refl= config.read<double> ("wrapping_refl");
  else
     fwrapping_refl=0.97;

  if (config.keyExists("PDEoption"))
     fPDEoption = config.read<std::string> ("PDEoption");
  else
     fPDEoption="";

  if (config.keyExists("SigmaAlpha"))
     fSigmaAlpha = config.read<double> ("SigmaAlpha");
  else
     fSigmaAlpha=3.4*deg;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceDetectorConstruction_5::~OpNoviceDetectorConstruction_5(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* OpNoviceDetectorConstruction_5::Construct()
{

// ------------- Materials -------------

  G4int nel;
  G4double density; 
  G4double temperature= 300*kelvin;
  G4double pressure= 2.0e-7*bar;
  G4double massfraction,z,a;
  G4double posx,posy,posz;

// Air
//
  G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);

  G4Material* Air = new G4Material("Air", density=1.29*mg/cm3, nel=2);
  Air->AddElement(N, 70.*perCent);
  Air->AddElement(O, 30.*perCent);
/*
// Water
//
  G4Element* H = new G4Element("Hydrogen", "H", z=1 , a=1.01*g/mole);

  G4Material* water = new G4Material("Water", density= 1.0*g/cm3, nelements=2);
  water->AddElement(H, 2);
  water->AddElement(O, 1);
*/
// Vacuum
  G4Material* Vacuum = new G4Material("Vacuum", density=2.376e-15*g/cm3, nel=1,
                                        kStateGas,temperature,pressure);
  Vacuum-> AddMaterial(Air, massfraction= 1.);
  
// LYSO
  G4Material *Scint_mat = new G4Material("Scint", density=7.4*g/cm3, nel = 4);
  G4Element* Lu = new G4Element("Lutethium", "Lu", z=71 , a=174.967*g/mole);
  G4Element* Si = new G4Element("Silicium", "Si", z=14 , a=28.086*g/mole);
  G4Element* Y = new G4Element("YYYY", "Y", z=39 , a=88.906*g/mole);
  G4Element* Ce = new G4Element("Cerium"  , "Ce", z=58 , a=140.116*g/mole);
  Scint_mat->AddElement(Lu, 71*perCent);
  Scint_mat->AddElement(Si, 7*perCent);
  Scint_mat->AddElement(O, 18*perCent);   
  Scint_mat->AddElement(Y, 4*perCent);
  G4Material *LYSO = new G4Material("LYSO", density=7.4*g/cm3, nel = 2);
  LYSO->AddMaterial(Scint_mat, 99.81*perCent);
  LYSO->AddElement(Ce, 0.19*perCent); 

// SiPM
  G4Material *SiPM_mat = new G4Material("SiPM", density=2.3290*g/cm3, nel = 1);
  SiPM_mat->AddElement(Si, 100*perCent);
  
// Optical Glue see https://en.wikipedia.org/wiki/Bisphenol_A
  G4int n_atoms;
  G4Material *Glue_mat = new G4Material("Optic_Glue", density=1.2*g/cm3, nel = 3);
  G4Element* C = new G4Element("Carbon", "C", z=6 , a=12.0107*g/mole);
  G4Element* H = new G4Element("Hydrogen", "H", z=1 , a=1.008*g/mole);
  Glue_mat ->AddElement(C, n_atoms = 15);
  Glue_mat ->AddElement(H, n_atoms = 16);
  Glue_mat ->AddElement(O, n_atoms = 2);   






//
// ------------ Generate & Add Material Properties Table ------------
//
//
// LYSO
//

  G4double photonEnergy[] =
            { 3.4238821301*eV,
3.3489075539*eV,
3.3007165749*eV,
3.2789833036*eV,
3.2639051572*eV,
3.2489735602*eV,
3.2341695226*eV,
3.2122188661*eV,
3.1905641638*eV,
3.1691994704*eV,
3.1204388746*eV,
3.001694505*eV,
2.9034780858*eV,
2.8396014482*eV,
2.7784748862*eV,
2.7095415194*eV,
2.6292577012*eV,
2.5263329462*eV,
2.3782261885*eV,
2.2465192539*eV,
2.1350376903*eV,
2.0757306907*eV};
            
  const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);
G4double photonEnergy_scint[] =
{ 3.3513513514*eV, 3.34682861*eV, 3.3423180593*eV, 3.3378196501*eV, 3.3333333333*eV, 3.3288590604*eV, 3.3243967828*eV, 3.3199464525*eV, 3.3155080214*eV, 3.3110814419*eV, 3.3066666667*eV, 3.3022636485*eV, 3.2978723404*eV, 3.2934926959*eV, 3.2891246684*eV, 3.2847682119*eV, 3.2804232804*eV, 3.2760898283*eV, 3.27176781*eV, 3.2674571805*eV, 3.2631578947*eV, 3.258869908*eV, 3.2545931759*eV, 3.250327654*eV, 3.2460732984*eV, 3.2418300654*eV, 3.2375979112*eV, 3.2333767927*eV, 3.2291666667*eV, 3.2249674902*eV, 3.2207792208*eV, 3.2166018158*eV, 3.2124352332*eV, 3.2082794308*eV, 3.2041343669*eV, 3.2*eV, 3.1958762887*eV, 3.1917631918*eV, 3.1876606684*eV, 3.1835686778*eV, 3.1794871795*eV, 3.1754161332*eV, 3.1713554987*eV, 3.1673052363*eV, 3.1632653061*eV, 3.1592356688*eV, 3.155216285*eV, 3.1512071156*eV, 3.1472081218*eV, 3.1432192649*eV, 3.1392405063*eV, 3.1352718078*eV, 3.1313131313*eV, 3.1273644388*eV, 3.1234256927*eV, 3.1194968553*eV, 3.1155778894*eV, 3.1116687578*eV, 3.1077694236*eV, 3.1038798498*eV, 3.1*eV, 3.0961298377*eV, 3.0922693267*eV, 3.0884184309*eV, 3.0845771144*eV, 3.0807453416*eV, 3.0769230769*eV, 3.073110285*eV, 3.0693069307*eV, 3.065512979*eV, 3.0617283951*eV, 3.0579531443*eV, 3.0541871921*eV, 3.0504305043*eV, 3.0466830467*eV, 3.0429447853*eV, 3.0392156863*eV, 3.035495716*eV, 3.0317848411*eV, 3.0280830281*eV, 3.0243902439*eV, 3.0207064555*eV, 3.0170316302*eV, 3.0133657351*eV, 3.0097087379*eV, 3.0060606061*eV, 3.0024213075*eV, 2.9987908102*eV, 2.9951690821*eV, 2.9915560917*eV, 2.9879518072*eV, 2.9843561974*eV, 2.9807692308*eV, 2.9771908764*eV, 2.9736211031*eV, 2.9700598802*eV, 2.966507177*eV, 2.962962963*eV, 2.9594272076*eV, 2.9558998808*eV, 2.9523809524*eV, 2.9488703924*eV, 2.945368171*eV, 2.9418742586*eV, 2.9383886256*eV, 2.9349112426*eV, 2.9314420804*eV, 2.9279811098*eV, 2.9245283019*eV, 2.9210836278*eV, 2.9176470588*eV, 2.9142185664*eV, 2.9107981221*eV, 2.9073856975*eV, 2.9039812646*eV, 2.9005847953*eV, 2.8971962617*eV, 2.8938156359*eV, 2.8904428904*eV, 2.8870779977*eV, 2.8837209302*eV, 2.8803716609*eV, 2.8770301624*eV, 2.8736964079*eV, 2.8703703704*eV, 2.8670520231*eV, 2.8637413395*eV, 2.860438293*eV, 2.8571428571*eV, 2.8538550058*eV, 2.8505747126*eV, 2.8473019518*eV, 2.8440366972*eV, 2.8407789233*eV, 2.8375286041*eV, 2.8342857143*eV, 2.8310502283*eV, 2.8278221209*eV, 2.8246013667*eV, 2.8213879408*eV, 2.8181818182*eV, 2.8149829739*eV, 2.8117913832*eV, 2.8086070215*eV, 2.8054298643*eV, 2.802259887*eV, 2.7990970655*eV, 2.7959413754*eV, 2.7927927928*eV, 2.7896512936*eV, 2.7865168539*eV, 2.7833894501*eV, 2.7802690583*eV, 2.7771556551*eV, 2.774049217*eV, 2.7709497207*eV, 2.7678571429*eV, 2.7647714604*eV, 2.7616926503*eV, 2.7586206897*eV, 2.7555555556*eV, 2.7524972253*eV, 2.7494456763*eV, 2.7464008859*eV, 2.7433628319*eV, 2.7403314917*eV, 2.7373068433*eV, 2.7342888644*eV, 2.731277533*eV, 2.7282728273*eV, 2.7252747253*eV, 2.7222832053*eV, 2.7192982456*eV, 2.7163198248*eV, 2.7133479212*eV, 2.7103825137*eV, 2.7074235808*eV, 2.7044711014*eV, 2.7015250545*eV, 2.6985854189*eV, 2.6956521739*eV, 2.6927252986*eV, 2.6898047722*eV, 2.6868905742*eV, 2.683982684*eV, 2.6810810811*eV, 2.6781857451*eV, 2.6752966559*eV, 2.6724137931*eV, 2.6695371367*eV, 2.6666666667*eV, 2.6638023631*eV, 2.660944206*eV, 2.6580921758*eV, 2.6552462527*eV, 2.6524064171*eV, 2.6495726496*eV, 2.6467449306*eV, 2.6439232409*eV, 2.6411075612*eV, 2.6382978723*eV, 2.6354941552*eV, 2.6326963907*eV, 2.6299045599*eV, 2.6271186441*eV, 2.6243386243*eV, 2.621564482*eV, 2.6187961985*eV, 2.6160337553*eV, 2.6132771338*eV, 2.6105263158*eV, 2.6077812829*eV, 2.6050420168*eV, 2.6023084995*eV, 2.5995807128*eV, 2.5968586387*eV, 2.5941422594*eV, 2.5914315569*eV, 2.5887265136*eV, 2.5860271116*eV, 2.5833333333*eV, 2.5806451613*eV, 2.577962578*eV, 2.5752855659*eV, 2.5726141079*eV, 2.5699481865*eV, 2.5672877847*eV, 2.5646328852*eV, 2.5619834711*eV, 2.5593395253*eV, 2.5567010309*eV, 2.5540679712*eV, 2.5514403292*eV, 2.5488180884*eV, 2.546201232*eV, 2.5435897436*eV, 2.5409836066*eV, 2.5383828045*eV, 2.5357873211*eV, 2.5331971399*eV, 2.5306122449*eV, 2.5280326198*eV, 2.5254582485*eV, 2.522889115*eV, 2.5203252033*eV, 2.5177664975*eV, 2.5152129817*eV, 2.5126646403*eV, 2.5101214575*eV, 2.5075834176*eV, 2.5050505051*eV, 2.5025227043*eV, 2.5*eV, 2.4974823766*eV, 2.4949698189*eV, 2.4924623116*eV, 2.4899598394*eV, 2.4874623872*eV, 2.4849699399*eV, 2.4824824825*eV, 2.48*eV, 2.4775224775*eV, 2.4750499002*eV, 2.4725822532*eV, 2.4701195219*eV, 2.4676616915*eV, 2.4652087475*eV, 2.4627606753*eV, 2.4603174603*eV, 2.4578790882*eV, 2.4554455446*eV, 2.453016815*eV, 2.4505928854*eV, 2.4481737414*eV, 2.4457593688*eV, 2.4433497537*eV, 2.4409448819*eV, 2.4385447394*eV, 2.4361493124*eV, 2.4337585868*eV, 2.431372549*eV, 2.4289911851*eV, 2.4266144814*eV, 2.4242424242*eV, 2.421875*eV, 2.4195121951*eV, 2.4171539961*eV, 2.4148003895*eV, 2.4124513619*eV, 2.4101068999*eV, 2.4077669903*eV, 2.4054316198*eV, 2.4031007752*eV, 2.4007744434*eV, 2.3984526112*eV, 2.3961352657*eV, 2.3938223938*eV, 2.3915139826*eV, 2.3892100193*eV, 2.3869104909*eV, 2.3846153846*eV, 2.3823246878*eV, 2.3800383877*eV, 2.3777564717*eV, 2.3754789272*eV, 2.3732057416*eV, 2.3709369025*eV, 2.3686723973*eV, 2.3664122137*eV, 2.3641563394*eV, 2.3619047619*eV, 2.3596574691*eV, 2.3574144487*eV, 2.3551756885*eV, 2.3529411765*eV, 2.3507109005*eV, 2.3484848485*eV, 2.3462630085*eV, 2.3440453686*eV, 2.3418319169*eV, 2.3396226415*eV, 2.3374175306*eV, 2.3352165725*eV, 2.3330197554*eV, 2.3308270677*eV, 2.3286384977*eV, 2.3264540338*eV, 2.3242736645*eV, 2.3220973783*eV, 2.3199251637*eV, 2.3177570093*eV, 2.3155929038*eV, 2.3134328358*eV, 2.311276794*eV, 2.3091247672*eV, 2.3069767442*eV, 2.3048327138*eV, 2.3026926648*eV, 2.3005565863*eV, 2.2984244671*eV, 2.2962962963*eV, 2.2941720629*eV, 2.292051756*eV, 2.2899353647*eV, 2.2878228782*eV, 2.2857142857*eV, 2.2836095764*eV, 2.2815087397*eV, 2.2794117647*eV, 2.277318641*eV, 2.2752293578*eV, 2.2731439047*eV, 2.2710622711*eV, 2.2689844465*eV, 2.2669104205*eV, 2.2648401826*eV, 2.2627737226*eV, 2.2607110301*eV, 2.2586520947*eV, 2.2565969063*eV, 2.2545454545*eV, 2.2524977293*eV, 2.2504537205*eV, 2.248413418*eV, 2.2463768116*eV, 2.2443438914*eV, 2.2423146474*eV, 2.2402890696*eV, 2.238267148*eV, 2.2362488729*eV, 2.2342342342*eV, 2.2322232223*eV, 2.2302158273*eV, 2.2282120395*eV, 2.2262118492*eV, 2.2242152466*eV, 2.2222222222*eV, 2.2202327663*eV, 2.2182468694*eV, 2.2162645219*eV, 2.2142857143*eV, 2.2123104371*eV, 2.2103386809*eV, 2.2083704363*eV, 2.206405694*eV, 2.2044444444*eV, 2.2024866785*eV, 2.2005323869*eV, 2.1985815603*eV, 2.1966341895*eV, 2.1946902655*eV, 2.192749779*eV, 2.1908127208*eV, 2.1888790821*eV, 2.1869488536*eV, 2.1850220264*eV, 2.1830985915*eV, 2.18117854*eV, 2.1792618629*eV, 2.1773485514*eV, 2.1754385965*eV, 2.1735319895*eV, 2.1716287215*eV, 2.1697287839*eV, 2.1678321678*eV, 2.1659388646*eV, 2.1640488656*eV, 2.1621621622*eV, 2.1602787456*eV, 2.1583986075*eV, 2.1565217391*eV, 2.1546481321*eV, 2.1527777778*eV, 2.1509106678*eV, 2.1490467938*eV, 2.1471861472*eV, 2.1453287197*eV, 2.143474503*eV, 2.1416234888*eV, 2.1397756687*eV, 2.1379310345*eV, 2.136089578*eV, 2.1342512909*eV, 2.1324161651*eV, 2.1305841924*eV, 2.1287553648*eV, 2.1269296741*eV, 2.1251071123*eV, 2.1232876712*eV, 2.121471343*eV, 2.1196581197*eV, 2.1178479932*eV, 2.1160409556*eV, 2.1142369991*eV, 2.1124361158*eV, 2.1106382979*eV, 2.1088435374*eV, 2.1070518267*eV, 2.1052631579*eV, 2.1034775233*eV, 2.1016949153*eV, 2.099915326*eV, 2.0981387479*eV, 2.0963651733*eV, 2.0945945946*eV, 2.0928270042*eV, 2.0910623946*eV, 2.0893007582*eV, 2.0875420875*eV, 2.0857863751*eV, 2.0840336134*eV, 2.0822837951*eV, 2.0805369128*eV, 2.0787929589*eV, 2.0770519263*eV, 2.0753138075*eV, 2.0735785953*eV, 2.0718462824*eV, 2.0701168614*eV, 2.0683903253*eV, 2.0666666667*eV, 2.0649458784*eV, 2.0632279534*eV, 2.0615128845*eV, 2.0598006645*eV, 2.0580912863*eV, 2.056384743*eV, 2.0546810273*eV, 2.0529801325*eV, 2.0512820513*eV, 2.0495867769*eV, 2.0478943022*eV, 2.0462046205*eV, 2.0445177246*eV, 2.0428336079*eV, 2.0411522634*eV, 2.0394736842*eV, 2.0377978636*eV, 2.0361247947*eV, 2.0344544709*eV, 2.0327868852*eV, 2.0311220311*eV, 2.0294599018*eV, 2.0278004906*eV, 2.0261437908*eV, 2.0244897959*eV, 2.0228384992*eV, 2.0211898941*eV, 2.0195439739*eV, 2.0179007323*eV, 2.0162601626*eV, 2.0146222583*eV, 2.012987013*eV, 2.0113544201*eV, 2.0097244733*eV, 2.008097166*eV, 2.0064724919*eV, 2.0048504446*eV, 2.0032310178*eV, 2.001614205*eV, 2*eV, 1.9983883965*eV, 1.9967793881*eV, 1.9951729686*eV, 1.9935691318*eV, 1.9919678715*eV, 1.9903691814*eV, 1.9887730553*eV, 1.9871794872*eV, 1.9855884708*eV, 1.984*eV, 1.9824140687*eV, 1.9808306709*eV, 1.9792498005*eV, 1.9776714514*eV, 1.9760956175*eV, 1.974522293*eV, 1.9729514718*eV, 1.9713831479*eV, 1.9698173153*eV, 1.9682539683*eV, 1.9666931007*eV, 1.9651347068*eV, 1.9635787807*eV, 1.9620253165*eV, 1.9604743083*eV, 1.9589257504*eV, 1.9573796369*eV, 1.9558359621*eV, 1.9542947203*eV, 1.9527559055*eV, 1.9512195122*eV, 1.9496855346*eV, 1.948153967*eV, 1.9466248038*eV, 1.9450980392*eV, 1.9435736677*eV, 1.9420516836*eV, 1.9405320814*eV, 1.9390148554*eV, 1.9375*eV, 1.9359875098*eV, 1.9344773791*eV, 1.9329696025*eV, 1.9314641745*eV, 1.9299610895*eV, 1.9284603421*eV, 1.926961927*eV, 1.9254658385*eV, 1.9239720714*eV, 1.9224806202*eV, 1.9209914795*eV, 1.919504644*eV, 1.9180201083*eV, 1.9165378671*eV, 1.9150579151*eV, 1.9135802469*eV, 1.9121048574*eV, 1.9106317411*eV, 1.909160893*eV, 1.9076923077*eV};

  const G4int nEntries_scint = sizeof(photonEnergy_scint)/sizeof(G4double);

  G4double scintilFast[] =
{ 5.653, 5.4636, 5.63734, 6.16206, 6.92776, 7.66037, 8.31898, 8.96946, 9.87718, 10.9865, 11.9857, 12.8428, 13.7793, 15.0148, 16.4099, 17.8478, 19.3287, 21.0805, 22.9021, 24.5062, 25.8105, 27.2142, 28.9335, 30.8429, 32.8791, 35.0609, 37.3185, 39.3839, 41.4014, 43.5675, 45.7439, 47.9825, 50.1839, 52.7438, 55.3837, 58.2847, 61.1579, 64.0054, 67.0658, 70.4437, 74.0542, 77.5629, 81.1023, 84.6858, 88.6336, 92.7958, 97.1218, 101.169, 104.992, 109.123, 114.307, 119.915, 125.469, 130.29, 135.238, 140.709, 146.97, 153.162, 159.004, 164.085, 169.451, 175.754, 182.885, 189.379, 195.153, 201.089, 208.267, 215.836, 223.271, 230.181, 237.281, 244.55, 251.886, 258.49, 264.626, 270.962, 278.568, 286.92, 294.748, 301.371, 307.717, 314.617, 322.14, 329.431, 336.257, 343.2, 350.314, 357.603, 364.459, 370.829, 376.748, 382.991, 389.198, 395.71, 402.268, 408.612, 413.849, 417.947, 421.539, 425.68, 430.191, 434.598, 438.754, 442.234, 445.78, 448.622, 451.009, 452.786, 454.942, 457.56, 460.04, 462.057, 463.242, 463.821, 463.377, 463.276, 463.127, 462.911, 461.394, 459.735, 459.022, 458.867, 457.579, 455.133, 452.718, 451.539, 450.123, 447.287, 443.166, 439.28, 435.461, 431.953, 427.739, 423.637, 419.6, 415.43, 410.295, 404.695, 400.104, 396.864, 393.958, 389.895, 385.277, 380.083, 374.983, 369.633, 365.02, 360.168, 355.582, 349.983, 344.681, 338.988, 334.079, 328.699, 323.472, 318.217, 313.087, 307.641, 301.744, 296.691, 292.221, 287.83, 282.369, 277.093, 272.321, 268.204, 263.798, 258.906, 254.012, 249.595, 245.158, 240.661, 235.879, 231.441, 227.343, 223.734, 219.742, 215.219, 210.283, 206.047, 202.765, 199.791, 196.716, 193.255, 189.514, 185.732, 182.466, 179.522, 176.568, 173.577, 170.509, 167.446, 164.057, 161.322, 158.893, 156.702, 153.922, 151.167, 148.247, 145.554, 142.892, 140.429, 137.981, 135.234, 132.773, 130.756, 129.134, 127.345, 125.319, 123.422, 121.626, 119.627, 117.285, 115.335, 114.014, 112.921, 111.578, 110.095, 108.645, 107.122, 105.699, 104.449, 103.153, 101.588, 99.9306, 98.3639, 97.0263, 95.9876, 95.2024, 94.3203, 93.1508, 91.948, 91.0619, 90.1571, 89.157, 88.0977, 86.8273, 84.9982, 83.0915, 81.7858, 81.1847, 80.5788, 79.6912, 78.7278, 77.9198, 77.3506, 77.0474, 76.4193, 75.4407, 74.0747, 73.1508, 72.6319, 72.5521, 72.2614, 71.7871, 70.9097, 69.8691, 68.7421, 67.7007, 67.0353, 66.7051, 66.4623, 65.9513, 65.3533, 64.9251, 64.4503, 63.6395, 62.5925, 61.7203, 61.0173, 60.3999, 60.0213, 59.7588, 59.301, 58.5086, 57.8103, 57.2736, 56.8014, 56.3713, 55.7925, 55.1853, 54.5648, 54.134, 53.5281, 52.9118, 52.2439, 51.7601, 51.3232, 50.7976, 50.2158, 49.5814, 49.1203, 48.7459, 48.3676, 47.8083, 47.3026, 46.7184, 46.0388, 45.3425, 44.717, 44.1642, 43.6988, 43.2002, 42.6547, 42.0509, 41.5613, 41.1375, 40.4938, 39.8158, 39.4343, 39.3372, 39.0962, 38.5009, 37.8593, 37.197, 36.7779, 36.4319, 36.3835, 36.0534, 35.7086, 35.1736, 34.7548, 34.1971, 33.7196, 33.2726, 32.8913, 32.3675, 31.8364, 31.4165, 31.0645, 30.6849, 30.2035, 29.8484, 29.6951, 29.5738, 29.2072, 28.7467, 28.4189, 28.2671, 28.0665, 27.5595, 26.9723, 26.517, 26.2508, 25.9854, 25.6909, 25.3453, 25.0105, 24.6232, 24.2962, 24.1093, 23.8839, 23.5489, 23.1205, 22.8357, 22.57, 22.3399, 21.9981, 21.8733, 21.7177, 21.4793, 21.1315, 20.761, 20.3485, 19.8555, 19.545, 19.3563, 19.3117, 19.1074, 18.9681, 18.9401, 18.9404, 18.7391, 18.1658, 17.5783, 17.2826, 17.2077, 17.0026, 16.735, 16.6009, 16.6244, 16.5614, 16.1924, 15.7444, 15.3775, 15.2063, 15.1864, 15.1662, 15.1082, 14.9224, 14.5938, 14.2524, 14.023, 13.933, 13.8722, 13.7555, 13.5214, 13.2853, 13.0202, 12.7442, 12.4969, 12.3436, 12.2443, 12.1348, 11.966, 11.7971, 11.5861, 11.3991, 11.2389, 11.1522, 11.0773, 10.9327, 10.7311, 10.4013, 10.2201, 10.2196, 10.3596, 10.2405, 9.91198, 9.54354, 9.43766, 9.46134, 9.54266, 9.46208, 9.3142, 8.98796, 8.77055, 8.60793, 8.6928, 8.74686, 8.70336, 8.40274, 8.15511, 8.06225, 8.07562, 7.98042, 7.78679, 7.66921, 7.56145, 7.43326, 7.19748, 7.0976, 7.16639, 7.26398, 7.12839, 6.86168, 6.7198, 6.64605, 6.51279, 6.27161, 6.17325, 6.14705, 6.16148, 6.19039, 6.29996, 6.33137, 6.2492, 6.12285, 5.98873, 5.79497, 5.58503, 5.4682, 5.50877, 5.59986, 5.68007, 5.68804, 5.5964, 5.44079, 5.27271, 5.11378, 5.03778, 5.07701, 5.12873, 4.99995, 4.79876, 4.62148, 4.58934, 4.573, 4.64084, 4.70288, 4.7405, 4.64024, 4.43192, 4.21036, 4.04867, 4.00833, 4.01152, 4.01571, 3.96784, 3.96336, 3.97369, 4.01847, 3.97173, 3.86762, 3.73496, 3.66361, 3.61999, 3.57782, 3.5273, 3.50585, 3.50042, 3.44388, 3.32714, 3.28898, 3.33563, 3.35992, 3.22947, 3.05819, 2.94839, 2.92942, 2.94832, 2.92595, 2.87991, 2.7922, 2.79748, 2.82908, 2.85229, 2.75509, 2.63947, 2.57423, 2.57102, 2.58364, 2.59854, 2.61219, 2.56296, 2.47742, 2.44607, 2.50589, 2.55083, 2.50669, 2.44329, 2.41881, 2.40857, 2.36062, 2.24879, 2.19017, 2.18642, 2.20356, 2.15279, 2.11331, 2.101, 2.14267, 2.13498, 2.05203, 1.93688, 1.83428, 1.86925, 1.99454, 2.1211, 2.09036, 1.98252, 1.84871, 1.79907, 1.752, 1.70636, 1.63422, 1.63009, 1.63046, 1.65704, 1.71014, 1.828 };//La normalizzazione non ha alcun effetto??????
            
  assert(sizeof(scintilFast) == sizeof(photonEnergy_scint));


  G4double photonEnergy_rindx[] =
            //{405*nm, 420*nm, 436*nm, 461*nm, 486*nm, 516*nm, 546*nm};
              {3.4444444444*eV, 3.397260274*eV, 3.3513513514*eV, 3.3066666667*eV, 3.2631578947*eV, 3.2207792208*eV, 3.1794871795*eV, 3.1392405063*eV, 3.1*eV, 3.0617283951*eV, 3.0243902439*eV, 2.9879518072*eV, 2.9523809524*eV, 2.9176470588*eV, 2.8837209302*eV, 2.8505747126*eV, 2.8181818182*eV, 2.7865168539*eV, 2.7555555556*eV, 2.7252747253*eV, 2.6956521739*eV, 2.6666666667*eV, 2.6382978723*eV, 2.6105263158*eV, 2.5833333333*eV, 2.5567010309*eV, 2.5306122449*eV, 2.5050505051*eV, 2.48*eV, 2.4554455446*eV, 2.431372549*eV, 2.4077669903*eV, 2.3846153846*eV, 2.3619047619*eV, 2.3396226415*eV, 2.3177570093*eV, 2.2962962963*eV, 2.2752293578*eV, 2.2545454545*eV, 2.2342342342*eV, 2.2142857143*eV, 2.1946902655*eV, 2.1754385965*eV, 2.1565217391*eV, 2.1379310345*eV, 2.1196581197*eV, 2.1016949153*eV, 2.0840336134*eV, 2.0666666667*eV, 2.0495867769*eV, 2.0327868852*eV, 2.0162601626*eV, 2*eV, 1.984*eV, 1.9682539683*eV, 1.9527559055*eV, 1.9375*eV, 1.9224806202*eV, 1.9076923077*eV};
  G4double refractiveIndex1[] =
            //{1.833, 1.827, 1.822, 1.818, 1.813, 1.810, 1.806};
            {1.8476156814, 1.8454919931, 1.8434638008, 1.8415253313, 1.8396712502, 1.8378966214, 1.8361968715, 1.8345677577, 1.8330053392, 1.8315059513, 1.8300661818, 1.8286828507, 1.8273529907, 1.8260738304, 1.8248427783, 1.8236574095, 1.8225154522, 1.8214147762, 1.8203533825, 1.8193293935, 1.818341044, 1.8173866734, 1.8164647178, 1.8155737038, 1.8147122417, 1.8138790202, 1.8130728007, 1.8122924128, 1.8115367496, 1.8108047636, 1.8100954628, 1.8094079073, 1.8087412061, 1.8080945137, 1.8074670277, 1.8068579857, 1.8062666634, 1.8056923722, 1.8051344567, 1.8045922932, 1.8040652879, 1.8035528746, 1.8030545139, 1.8025696913, 1.8020979157, 1.8016387184, 1.8011916518, 1.8007562881, 1.8003322185, 1.7999190519, 1.7995164144, 1.7991239477, 1.7987413092, 1.7983681705, 1.7980042169, 1.7976491468, 1.797302671, 1.796964512, 1.7966344035};
  const G4int nEntries_rindx = sizeof(photonEnergy_rindx)/sizeof(G4double);

  assert(sizeof(refractiveIndex1) == sizeof(photonEnergy_rindx));


G4double photonEnergy_absorb[] =
{3.4444444444*eV,
3.4159779614*eV,
3.3879781421*eV,
3.3604336043*eV,
3.3513513514*eV,
3.3333333333*eV,
3.3007165749*eV,
3.2804232804*eV,
3.2654781029*eV,
3.2539014015*eV,
3.2309840952*eV,
3.2083956469*eV,
3.1972112067*eV,
3.1750827698*eV,
3.1424625249*eV,
3.1*eV,
3.0283862365*eV,
2.7722072806*eV,
2.2656802431*eV,
2.0902058684*eV,
2.0666666667*eV,
1.9076923077*eV};

  G4double absorption[] =
{42.*cm,
42.*cm,
42.*cm,
42.*cm,
42.*cm,
42.*cm,
42.*cm,
42.*cm,
42.*cm,
42.*cm,
42.*cm,
42.*cm,
42.*cm,
42.*cm,
42.*cm,
42.*cm,
42.*cm,
42.*cm,
42.*cm,
42.*cm,
42.*cm,
42.*cm};

  const G4int nEntries_absorb = sizeof(photonEnergy_absorb)/sizeof(G4double);
  assert(sizeof(absorption) == sizeof(photonEnergy_absorb));


  G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();

  myMPT1->AddProperty("RINDEX",photonEnergy_rindx, refractiveIndex1,nEntries_rindx) 
  -> SetSpline(true);//LYSO corrected: see http://www.hep.caltech.edu/~zhu/papers/N49-1.pdf
  
  myMPT1->AddProperty("FASTCOMPONENT",photonEnergy_scint, scintilFast,     nEntries_scint)
        ->SetSpline(true);//LYSO corrected  see: http://www.hep.caltech.edu/~zhu/papers/N49-1.pdf

  myMPT1->AddProperty("ABSLENGTH",photonEnergy_absorb, absorption, nEntries_absorb)->SetSpline(true);
  //myMPT1->AddConstProperty("ABSLENGTH",42.*cm); //LYSO corrected: see https:www.google.it/url?sa=t&rct=j&q=&esrc=s&source=web&cd=2&ved=0ahUKEwio_sGQhLXSAhWDORQKHSrzD68QFggnMAE&url=http%3A%2F%2Fssd-rd.web.cern.ch%2Fssd-rd%2Fpad_hpd%2Fliterature%2Fcrystal_nim_deleo.doc&usg=AFQjCNFuyfsP-3YJHQ8wyVagA68w9W-vDA&sig2=VLZV_60i3JDYVoodeLmmZA&bvm=bv.148441817,d.d24
  
  myMPT1->AddConstProperty("SCINTILLATIONYIELD",30./keV);//LYSO corrected: see wikipedia
  
  myMPT1->AddConstProperty("RESOLUTIONSCALE",1.0);//1.0 -> exactly poissonian fluctuation of the number of photon 
  
  //myMPT1->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  
  myMPT1->AddConstProperty("FASTTIMECONSTANT",40.*ns);//LYSO corrected: see http://www.hep.caltech.edu/~zhu/papers/N49-1.pdf
  
  myMPT1->AddConstProperty("YIELDRATIO",1.);// = fasttimecomponent_yield/total_yield

  G4cout << "LYSO G4MaterialPropertiesTable" << G4endl;
  myMPT1->DumpTable();

  LYSO->SetMaterialPropertiesTable(myMPT1);

  // Set the Birks Constant for the LYSO scintillator
  LYSO->GetIonisation()->SetBirksConstant(1.027*mm/MeV);//LYSO corrected: see http://www.aesj.or.jp/publication/pnst001/data/218.pdf


//
// Air
//
  G4double refractiveIndex2[] =
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00};
  assert(sizeof(refractiveIndex2) == sizeof(photonEnergy));
  
  G4MaterialPropertiesTable* myMPT2 = new G4MaterialPropertiesTable();
  myMPT2->AddProperty("RINDEX", photonEnergy, refractiveIndex2, nEntries);

  G4cout << "Air G4MaterialPropertiesTable" << G4endl;
  myMPT2->DumpTable();
  
  
  //air->SetMaterialPropertiesTable(myMPT2);
  Vacuum->SetMaterialPropertiesTable(myMPT2);

//
// Optical Glue
//

  G4MaterialPropertiesTable* myMPT3 = new G4MaterialPropertiesTable();
  G4double refractiveIndex3[] =
            { 1.57, 1.57, 1.57, 1.57, 1.57, 1.57, 1.57,
              1.57, 1.57, 1.57, 1.57, 1.57, 1.57, 1.57,
              1.57, 1.57, 1.57, 1.57, 1.57, 1.57, 1.57,
              1.57};
              
   assert(sizeof(refractiveIndex3) == sizeof(photonEnergy));
  
  myMPT3->AddProperty("RINDEX", photonEnergy, refractiveIndex3, nEntries);


  G4cout << "Optical Glue G4MaterialPropertiesTable" << G4endl;
  myMPT3->DumpTable();
  Glue_mat->SetMaterialPropertiesTable(myMPT3);

//
// ------------- Volumes --------------

// The experimental Hall
//
  G4Box* expHall_box = new G4Box("World",fExpHall_x,fExpHall_y,fExpHall_z);

  G4LogicalVolume* expHall_log
    = new G4LogicalVolume(expHall_box,Vacuum,"World",0,0,0);

  G4VPhysicalVolume* expHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),expHall_log,"World",0,false,0);

// Small spaces around the crystal useful to define different surfaces around the crystal
  //Lateral surfaces
  G4Box* world_lateral_x_crystal_box = new G4Box("Lateral_world_x",0.2*mm,fCryst_y,fCryst_z);
  G4Box* world_lateral_y_crystal_box = new G4Box("Lateral_world_y",fCryst_x,0.2*mm,fCryst_z);
  G4Box* world_lateral_z_crystal_box = new G4Box("Lateral_world_z",fCryst_x,fCryst_y,0.2*mm);

  G4LogicalVolume* world_lateral_x_crystal_log = new G4LogicalVolume(world_lateral_x_crystal_box,Vacuum,"Lateral_world_x",0,0,0);
  G4LogicalVolume* world_lateral_y_crystal_log = new G4LogicalVolume(world_lateral_y_crystal_box,Vacuum,"Lateral_world_y",0,0,0);
  G4LogicalVolume* world_lateral_z_crystal_log = new G4LogicalVolume(world_lateral_z_crystal_box,Vacuum,"Lateral_world_z",0,0,0);

  G4VPhysicalVolume* world_lateral_x1_crystal_phys = 
new G4PVPlacement(0,G4ThreeVector(-fCryst_x-0.2*mm,0.,0.),world_lateral_x_crystal_log,"Lateral_world",expHall_log,false,0);
  G4VPhysicalVolume* world_lateral_x2_crystal_phys = 
new G4PVPlacement(0,G4ThreeVector(+fCryst_x+0.2*mm,0.,0.),world_lateral_x_crystal_log,"Lateral_world",expHall_log,false,0);
  G4VPhysicalVolume* world_lateral_y1_crystal_phys = 
new G4PVPlacement(0,G4ThreeVector(0.,-fCryst_x-0.2,0.),world_lateral_y_crystal_log,"Lateral_world",expHall_log,false,0);
  G4VPhysicalVolume* world_lateral_y2_crystal_phys = 
new G4PVPlacement(0,G4ThreeVector(0.,+fCryst_x+0.2,0.),world_lateral_y_crystal_log,"Lateral_world",expHall_log,false,0);
  G4VPhysicalVolume* world_lateral_z_crystal_phys = 
new G4PVPlacement(0,G4ThreeVector(0.,0.,-fCryst_z-0.2),world_lateral_z_crystal_log,"Lateral_world",expHall_log,false,0);
                     
// The LYSO crystal
  G4Box* crystal_box = new G4Box("Crystal",fCryst_x,fCryst_y,fCryst_z);

  G4LogicalVolume* crystal_log
    = new G4LogicalVolume(crystal_box,LYSO,"Crystal",0,0,0);

  G4VPhysicalVolume* crystal_phys
    = new G4PVPlacement(0,G4ThreeVector(),crystal_log,"Crystal",
                        expHall_log,false,0);
                       
                        
// Optical glue
    G4Box* Glue_box = new G4Box("Optic_Glue",fGlue_x,fGlue_y,fGlue_z);

  G4LogicalVolume* Glue_log
    = new G4LogicalVolume(Glue_box,Glue_mat,"Optic_Glue",0,0,0);

  G4VPhysicalVolume* Glue_phys0
    = new G4PVPlacement(0,G4ThreeVector(posx=-fCryst_x*0.666667 ,posy=-fCryst_x*0.666667, posz = fCryst_z + fGlue_z ),Glue_log,"Optic_Glue",
                        expHall_log,false,0);   
  G4VPhysicalVolume* Glue_phys1
    = new G4PVPlacement(0,G4ThreeVector(posx=-fCryst_x*0.666667 ,posy=+fCryst_x*0.666667, posz = fCryst_z + fGlue_z ),Glue_log,"Optic_Glue",
                        expHall_log,false,1);  
  G4VPhysicalVolume* Glue_phys2
    = new G4PVPlacement(0,G4ThreeVector(posx=+fCryst_x*0.666667 ,posy=-fCryst_x*0.666667, posz = fCryst_z + fGlue_z ),Glue_log,"Optic_Glue",
                        expHall_log,false,2);  
  G4VPhysicalVolume* Glue_phys3
    = new G4PVPlacement(0,G4ThreeVector(posx=+fCryst_x*0.666667 ,posy=+fCryst_x*0.666667, posz = fCryst_z + fGlue_z ),Glue_log,"Optic_Glue",
                        expHall_log,false,3);    
  G4VPhysicalVolume* Glue_phys4
    = new G4PVPlacement(0,G4ThreeVector(posx=0 ,posy=0, posz = fCryst_z + fGlue_z ),Glue_log,"Optic_Glue",
                        expHall_log,false,4);   

                        
// SiPM 
  G4Box* SiPM_box = new G4Box("SiPM",fSiPM_x,fSiPM_y,fSiPM_z);

  G4LogicalVolume* SiPM_log
    = new G4LogicalVolume(SiPM_box,SiPM_mat,"SiPM",0,0,0);

  G4VPhysicalVolume* SiPM_phys0
    = new G4PVPlacement(0,G4ThreeVector(posx=-fCryst_x*0.666667 ,posy=-fCryst_x*0.666667, posz = fCryst_z + 2*fGlue_z + fSiPM_z ),SiPM_log,"SiPM",
                        expHall_log,false,0);  
  G4VPhysicalVolume* SiPM_phys1
    = new G4PVPlacement(0,G4ThreeVector(posx=-fCryst_x*0.666667 ,posy=+fCryst_x*0.666667, posz = fCryst_z + 2*fGlue_z + fSiPM_z ),SiPM_log,"SiPM",
                        expHall_log,false,1);  
  G4VPhysicalVolume* SiPM_phys2
    = new G4PVPlacement(0,G4ThreeVector(posx=+fCryst_x*0.666667 ,posy=-fCryst_x*0.666667, posz = fCryst_z + 2*fGlue_z + fSiPM_z ),SiPM_log,"SiPM",
                        expHall_log,false,2);  
  G4VPhysicalVolume* SiPM_phys3
    = new G4PVPlacement(0,G4ThreeVector(posx=+fCryst_x*0.666667 ,posy=+fCryst_x*0.666667, posz = fCryst_z + 2*fGlue_z + fSiPM_z ),SiPM_log,"SiPM",
                        expHall_log,false,3);  
  G4VPhysicalVolume* SiPM_phys4
    = new G4PVPlacement(0,G4ThreeVector(posx=0 ,posy=0, posz = fCryst_z + 2*fGlue_z + fSiPM_z ),SiPM_log,"SiPM",
                        expHall_log,false,4);  
       

// ------------- Surfaces --------------
//

// LYSO-air surface
//
// LYSO-air back surface
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

  opLYSOFrontSurface->SetType(dielectric_dielectric);
  opLYSOLateralx1Surface->SetType(dielectric_dielectric);
  opLYSOLateralx2Surface->SetType(dielectric_dielectric);
  opLYSOLateraly1Surface->SetType(dielectric_dielectric);
  opLYSOLateraly2Surface->SetType(dielectric_dielectric);
  opLYSOFrontSurface->SetFinish(G4OpticalSurfaceFinish(fsurface_type));
  opLYSOLateralx1Surface->SetFinish(G4OpticalSurfaceFinish(fsurface_type));
  opLYSOLateralx2Surface->SetFinish(G4OpticalSurfaceFinish(fsurface_type));
  opLYSOLateraly1Surface->SetFinish(G4OpticalSurfaceFinish(fsurface_type));
  opLYSOLateraly2Surface->SetFinish(G4OpticalSurfaceFinish(fsurface_type));

  if(fsurface_type == 0 || fsurface_type == 1)
  {
    opLYSOFrontSurface->SetModel(glisur);
    opLYSOLateralx1Surface->SetModel(glisur);
    opLYSOLateralx2Surface->SetModel(glisur);
    opLYSOLateraly1Surface->SetModel(glisur);
    opLYSOLateraly2Surface->SetModel(glisur);
  }
  else
  { 
    opLYSOFrontSurface->SetModel(unified);
    opLYSOLateralx1Surface->SetModel(unified);
    opLYSOLateralx2Surface->SetModel(unified);
    opLYSOLateraly1Surface->SetModel(unified);
    opLYSOLateraly2Surface->SetModel(unified);
  }

  const G4int NUM = 2;
  G4double XX[NUM] = {1.8*eV,3.4*eV} ; 
  G4double refractiveIndex[NUM] = {1.002, 1.002};
  G4double specularLobe[NUM]    = {0.5, 0.5};
  G4double specularSpike[NUM]   = {0.5, 0.5};
  G4double backScatter[NUM]     = {0., 0.};
  G4double REFLECT[NUM] = {fwrapping_refl,fwrapping_refl};
  G4MaterialPropertiesTable* myST1 = new G4MaterialPropertiesTable();

  if(fsurface_type >=2)
  {
     opLYSOFrontSurface->SetSigmaAlpha(fSigmaAlpha);
     opLYSOLateralx1Surface->SetSigmaAlpha(fSigmaAlpha);
     opLYSOLateralx2Surface->SetSigmaAlpha(fSigmaAlpha);
     opLYSOLateraly1Surface->SetSigmaAlpha(fSigmaAlpha);
     opLYSOLateraly2Surface->SetSigmaAlpha(fSigmaAlpha);

     myST1->AddProperty("SPECULARLOBECONSTANT",  XX, specularLobe,   NUM);
     myST1->AddProperty("SPECULARSPIKECONSTANT", XX, specularSpike,   NUM);
     myST1->AddProperty("BACKSCATTERCONSTANT",   XX, backScatter,     NUM);
  }

  if(fsurface_type == 1 || fsurface_type == 2 || fsurface_type == 4 || fsurface_type == 5)
  {
     myST1->AddProperty("REFLECTIVITY", XX, REFLECT,NUM);
  }

  if(fsurface_type == 2 || fsurface_type == 5)
     myST1->AddProperty("RINDEX", XX, refractiveIndex, NUM);

  if (fsurface_type >=1 && fsurface_type<=5)
  {
    opLYSOFrontSurface->SetMaterialPropertiesTable(myST1);
    opLYSOLateralx1Surface->SetMaterialPropertiesTable(myST1);
    opLYSOLateralx2Surface->SetMaterialPropertiesTable(myST1);
    opLYSOLateraly1Surface->SetMaterialPropertiesTable(myST1);
    opLYSOLateraly2Surface->SetMaterialPropertiesTable(myST1);
  }


//
//SiPM optical surface
//

G4OpticalSurface* SiPM_opsurf= new G4OpticalSurface("SiPM_opsurf",glisur,polished, dielectric_metal);
  //new G4LogicalBorderSurface("SiPM_opsurf", Glue_phys,SiPM_phys,SiPM_opsurf);

G4MaterialPropertiesTable* myMPT5 = new G4MaterialPropertiesTable();
  G4double energySilicon[] =
            {4.96, 4.7692307692, 4.5925925926, 4.4285714286, 4.275862069, 4.1333333333, 4, 3.875, 3.7575757576, 3.6470588235, 3.5428571429, 3.4444444444, 3.3513513514, 3.2631578947, 3.1794871795, 3.1, 3.0243902439, 2.9523809524, 2.8837209302, 2.8181818182, 2.7555555556, 2.6956521739, 2.6382978723, 2.5833333333, 2.5306122449, 2.48, 2.431372549, 2.3846153846, 2.3396226415, 2.2962962963, 2.2545454545, 2.2142857143, 2.1754385965, 2.1379310345, 2.1016949153, 2.0666666667, 2.0327868852, 2, 1.9682539683, 1.9375, 1.9076923077, 1.8787878788, 1.8507462687, 1.8235294118, 1.7971014493, 1.7714285714, 1.7464788732, 1.7222222222, 1.698630137, 1.6756756757, 1.6533333333, 1.6315789474, 1.6103896104, 1.5897435897, 1.5696202532, 1.55, 1.5308641975, 1.512195122, 1.4939759036, 1.4761904762, 1.4588235294, 1.4418604651, 1.4252873563, 1.4090909091, 1.393258427, 1.3777777778, 1.3626373626, 1.347826087, 1.3333333333, 1.3191489362, 1.3052631579, 1.2916666667, 1.2783505155, 1.2653061224, 1.2525252525, 1.24};

  G4double Re_n[] = 
            {1.694, 1.8, 2.129, 3.052, 4.426, 5.055, 5.074, 5.102, 5.179, 5.293, 5.483, 6.014, 6.863, 6.548, 5.976, 5.587, 5.305, 5.091, 4.925, 4.793, 4.676, 4.577, 4.491, 4.416, 4.348, 4.293, 4.239, 4.192, 4.15, 4.11, 4.077, 4.044, 4.015, 3.986, 3.962, 3.939, 3.916, 3.895, 3.879, 3.861, 3.844, 3.83, 3.815, 3.8, 3.787, 3.774, 3.762, 3.751, 3.741, 3.732, 3.723, 3.714, 3.705, 3.696, 3.688, 3.681, 3.674, 3.668, 3.662, 3.656, 3.65, 3.644, 3.638, 3.632, 3.626, 3.62, 3.614, 3.608, 3.602, 3.597, 3.592, 3.587, 3.582, 3.578, 3.574, 3.57};

  G4double Im_n[] = 
            {3.666, 4.072, 4.69, 5.258, 5.16, 4.128, 3.559, 3.269, 3.085, 2.951, 2.904, 2.912, 2.051, 0.885, 0.465, 0.303, 0.22, 0.167, 0.134, 0.109, 0.091, 0.077, 0.064, 0.057, 0.05, 0.045, 0.039, 0.036, 0.033, 0.03, 0.028, 0.026, 0.024, 0.023, 0.021, 0.02, 0.018, 0.017, 0.016, 0.015, 0.015, 0.014, 0.013, 0.012, 0.011, 0.011, 0.011, 0.01, 0.009, 0.008, 0.008, 0.007, 0.007, 0.006, 0.006, 0.005, 0.005, 0.005, 0.004, 0.004, 0.004, 0.003, 0.003, 0.003, 0.002, 0.002, 0.002, 0.002, 0.002, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001};
  G4double Eff[] = 
            {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};

  G4double Refl[] = 
            {0.672612594, 0.7051739998, 0.7320895527, 0.7229564109, 0.6842353612, 0.623487589, 0.5904758352, 0.5741303379, 0.5656774122, 0.5617493182, 0.5653802759, 0.5829110024, 0.5842708013, 0.5465022924, 0.5109736438, 0.4860210277, 0.4668532597, 0.4515215805, 0.4391232489, 0.4289072653, 0.4195856995, 0.4114859504, 0.4042813942, 0.3978791808, 0.3919647328, 0.3871055313, 0.3822645215, 0.3779990919, 0.3741420137, 0.3704285222, 0.3673359361, 0.3642162212, 0.3614517486, 0.3586671346, 0.3563449752, 0.3541066718, 0.3518536143, 0.3497852413, 0.3482013386, 0.3464114799, 0.3447139284, 0.3433093164, 0.3417986254, 0.340281901, 0.3389624224, 0.3376390016, 0.3364132847, 0.3352856254, 0.33425759, 0.3333299988, 0.3324006727, 0.3314686517, 0.3305348304, 0.3295983594, 0.3287643932, 0.3280328678, 0.3273003059, 0.3266713015, 0.3260410062, 0.3254099758, 0.3247779306, 0.3241446496, 0.3235105692, 0.3228754694, 0.3222391901, 0.3216020461, 0.3209638781, 0.3203246844, 0.3196844635, 0.3191500638, 0.3186150452, 0.3180793106, 0.317542859, 0.3171131809, 0.316683043, 0.3162524448};
              
  assert(sizeof(energySilicon) == sizeof(Re_n) && sizeof(energySilicon) == sizeof(Im_n) && sizeof(energySilicon) == sizeof(Eff) );
  const G4int nEnSi = sizeof(energySilicon)/sizeof(G4double);
  myMPT5->AddProperty("EFFICIENCY",energySilicon, Eff, nEnSi);
  myMPT5->AddProperty("REALRINDEX", energySilicon, Re_n, nEnSi);
  myMPT5->AddProperty("IMAGINARYRINDEX", energySilicon, Im_n, nEnSi);
  //myMPT5->AddProperty("REFLECTIVITY",energySilicon, Refl, nEnSi);
  G4cout << "Silicon G4MaterialPropertiesTable" << G4endl;
  SiPM_opsurf->SetMaterialPropertiesTable(myMPT5); 
  new G4LogicalSkinSurface("SiPM_surf",SiPM_log,SiPM_opsurf);

/*
  G4OpticalSurface* SiPM_opsurf=
     new G4OpticalSurface("SiPM_opsurf",glisur,polished,
                           dielectric_metal);
  G4double SiPM_EFF[]=
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00 };
              
  G4double SiPM_REFL[]=
            { 0.,0.,0.,0.,0.,0.,0.,
              0.,0.,0.,0.,0.,0.,0.,
              0.,0.,0.,0.,0.,0.,0.,
              0.};
  const G4int numentries_SiPM = sizeof(SiPM_EFF)/sizeof(G4double);
  G4MaterialPropertiesTable* SiPM_mt = new G4MaterialPropertiesTable();
  SiPM_mt->AddProperty("EFFICIENCY",photonEnergy,SiPM_EFF,numentries_SiPM);
  SiPM_mt->AddProperty("REFLECTIVITY",photonEnergy,SiPM_REFL,numentries_SiPM);
  G4cout << "SiPM Surface G4MaterialPropertiesTable" << G4endl;
  SiPM_mt->DumpTable();
  SiPM_opsurf->SetMaterialPropertiesTable(SiPM_mt);
  new G4LogicalSkinSurface("SiPM_surf",SiPM_log,SiPM_opsurf);
*/



//always return the physical World
  return expHall_phys;
}


void OpNoviceDetectorConstruction_5::ConstructSDandField()
{
  // G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // 
  // Sensitive detectors
  //
  SiPM_SD* SiPM = new SiPM_SD("SiPM_detector", "SiPMHitsCollection", 5,fPDEoption);
  SetSensitiveDetector("SiPM",SiPM);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
