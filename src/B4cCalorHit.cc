
#include "B4cCalorHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

G4ThreadLocal G4Allocator<B4cCalorHit>* B4cCalorHitAllocator = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cCalorHit::B4cCalorHit()
 : G4VHit(),
   fEdep(0.),
   fTrackLength(0.),
   fPh_time(1, 0.),
   //fPh_energy(5000, 0.), 
   fPh_SiPM_number(1,-1),
   //fPh_x_hit(5000,-20.),
   //fPh_y_hit(5000,-20.),
   fPh_creator_proc(1,-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cCalorHit::~B4cCalorHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cCalorHit::B4cCalorHit(const B4cCalorHit& right)
  : G4VHit()
{
  fEdep        = right.fEdep;
  fTrackLength = right.fTrackLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const B4cCalorHit& B4cCalorHit::operator=(const B4cCalorHit& right)
{
  fEdep        = right.fEdep;
  fTrackLength = right.fTrackLength;

  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int B4cCalorHit::operator==(const B4cCalorHit& right) const
{
  return ( this == &right ) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cCalorHit::Print()
{
  G4cout
     << "Edep: " 
     << std::setw(7) << G4BestUnit(fEdep,"Energy")
     << " track length: " 
     << std::setw(7) << G4BestUnit( fTrackLength,"Length")
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cCalorHit::GetTime(std::vector<G4float> &event_Ph_time) const
{
   event_Ph_time = fPh_time;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<G4float> B4cCalorHit::GetEnergy() const
{
   return fPh_energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cCalorHit::GetCreatorProc(std::vector<G4int> &event_Ph_creator_proc) const
{
   event_Ph_creator_proc = fPh_creator_proc;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cCalorHit::GetSiPM_number(std::vector<G4int> &event_Ph_SiPM_number) const
{
   event_Ph_SiPM_number = fPh_SiPM_number;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<G4int> B4cCalorHit::GetSiPMNumber() const
{
   return fPh_SiPM_number;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<G4float> B4cCalorHit::GetPh_x_hit() const
{
   return fPh_x_hit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<G4float> B4cCalorHit::GetPh_y_hit() const
{
   return fPh_y_hit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void B4cCalorHit::AddPhotonInfo(const G4double &time,const G4double &PhotEnergy,const G4int &CreatorProcess,const G4int &SiPMNumber,const G4double &Ph_x_hit, const G4double &Ph_y_hit)
{
   //for (unsigned i=0;i<1000;i++)
   //   std::cout<<i<<"\t"<<fPh_time.at(i)<<std::endl;

   if (fPh_time.at(0)==0.)//it's the first hit
   {
      fPh_time.at(0) = time;
      //fPh_energy.at(0) = PhotEnergy;
      fPh_creator_proc.at(0) = CreatorProcess;
      fPh_SiPM_number.at(0) = SiPMNumber;
      //fPh_x_hit.at(0) = Ph_x_hit;
      //fPh_y_hit.at(0) = Ph_y_hit;
   }
   else
   {
      std::vector<G4float>::iterator it_time = fPh_time.begin();
      //std::vector<G4float>::iterator it_energy = fPh_energy.begin();
      //std::vector<G4float>::iterator it_x_hit = fPh_x_hit.begin();
      //std::vector<G4float>::iterator it_y_hit = fPh_y_hit.begin();
      std::vector<G4int>::iterator it_creator_proc = fPh_creator_proc.begin();
      std::vector<G4int>::iterator it_SiPM_number = fPh_SiPM_number.begin();
      bool flag=0;
      while(it_time<fPh_time.end() && flag==0)
      {
         if(time<*it_time || *it_time==0)
            flag=1;
         else
         {
            it_time++;
            //it_energy++;
            //it_x_hit++;
            //it_y_hit++;
            it_creator_proc++;
            it_SiPM_number++;
         }
      }
      if(flag==1)
      {
         fPh_time.insert (it_time,time); 
         //fPh_time.pop_back();
         //fPh_energy.insert (it_energy,PhotEnergy); 
         //fPh_energy.pop_back();
         //fPh_x_hit.insert (it_x_hit,Ph_x_hit); 
         //fPh_x_hit.pop_back();
         //fPh_y_hit.insert (it_y_hit,Ph_y_hit); 
         //fPh_y_hit.pop_back();
         fPh_creator_proc.insert (it_creator_proc,CreatorProcess); 
         //fPh_creator_proc.pop_back();
         fPh_SiPM_number.insert (it_SiPM_number,SiPMNumber); 
         //fPh_SiPM_number.pop_back();
      }
      else
         if(it_time==fPh_time.end())
         {
            fPh_time.push_back(time); 
            fPh_creator_proc.push_back(CreatorProcess); 
            fPh_SiPM_number.push_back(SiPMNumber); 
         }
   }

/*   if (fPh_time.at(0)==0.)//it's the first hit
   {
      fPh_time.at(0) = time;
      fPh_energy.at(0) = PhotEnergy;
      fPh_creator_proc.at(0) = CreatorProcess;
      fPh_SiPM_number.at(0) = SiPMNumber;
      fPh_x_hit.at(0) = Ph_x_hit;
      fPh_y_hit.at(0) = Ph_y_hit;
   }
   else
   {
      unsigned it=0;
      while(time > fPh_time.at(it) && fPh_time.at(it)>0. && (it<1000))
      {
         it++;
         if(it==1000)
            break;
      }
      if(it < 999)
      {
         unsigned it2 = 998;
         for( ; it2>it ; it2--)
         {
            fPh_time.at(it2+1) = fPh_time.at(it2);
            fPh_energy.at(it2+1) = fPh_energy.at(it2);
            fPh_creator_proc.at(it2+1) = fPh_creator_proc.at(it2);
            fPh_SiPM_number.at(it2+1) = fPh_SiPM_number.at(it2);
            fPh_x_hit.at(it2+1) = fPh_x_hit.at(it2);
            fPh_y_hit.at(it2+1) = fPh_y_hit.at(it2);
         }
         fPh_time.at(it) = time;
         fPh_energy.at(it) = PhotEnergy;
         fPh_creator_proc.at(it) = CreatorProcess;
         fPh_SiPM_number.at(it) = SiPMNumber;
         fPh_x_hit.at(it) = Ph_x_hit;
         fPh_y_hit.at(it) = Ph_y_hit;
      }
      else
         if(it==999)
         {
            fPh_time.at(it) = time;
            fPh_energy.at(it) = PhotEnergy;
            fPh_creator_proc.at(it) = CreatorProcess;
            fPh_SiPM_number.at(it) = SiPMNumber;
            fPh_x_hit.at(it) = Ph_x_hit;
            fPh_y_hit.at(it) = Ph_y_hit;
         }
   }
*/
}

