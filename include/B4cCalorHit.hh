
#ifndef B4cCalorHit_h
#define B4cCalorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"
#include <vector>

/// Calorimeter hit class
///
/// It defines data members to store the the energy deposit and track lengths
/// of charged particles in a selected volume:
/// - fEdep, fTrackLength

class B4cCalorHit : public G4VHit
{
  public:
    B4cCalorHit();
    B4cCalorHit(const B4cCalorHit&);
    virtual ~B4cCalorHit();

    // operators
    const B4cCalorHit& operator=(const B4cCalorHit&);
    G4int operator==(const B4cCalorHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw() {}
    virtual void Print();

    // methods to handle data
    void Add(G4double de, G4double dl);
    void AddPhotonInfo(const G4double &time,const G4double &PhotEnergy,const G4int &CreatorProcess,const G4int &SiPMNumber,const G4double &Ph_x_hit, const G4double &Ph_y_hit);

    // get methods
    G4double GetEdep() const;
    G4double GetTrackLength() const;
    //std::vector<G4float> GetTime() const;
    void GetTime(std::vector<G4float> &event_Ph_time) const;
    std::vector<G4float> GetEnergy() const;
    //std::vector<G4int> GetCreatorProc() const; 
    void GetCreatorProc(std::vector<G4int> &event_Ph_creator_proc) const;
    void GetSiPM_number(std::vector<G4int> &event_Ph_SiPM_number) const;
    std::vector<G4int> GetSiPMNumber() const;
    std::vector<G4float> GetPh_x_hit() const;
    std::vector<G4float> GetPh_y_hit() const;

  private:
    G4double fEdep;        ///< Energy deposit in the sensitive volume
    G4double fTrackLength; ///< Track length in the  sensitive volume
    std::vector<G4float> fPh_time;
    std::vector<G4float> fPh_energy; 
    std::vector<G4int> fPh_creator_proc;
    std::vector<G4int> fPh_SiPM_number;
    std::vector<G4float> fPh_x_hit; 
    std::vector<G4float> fPh_y_hit; 

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<B4cCalorHit> B4cCalorHitsCollection;

extern G4ThreadLocal G4Allocator<B4cCalorHit>* B4cCalorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* B4cCalorHit::operator new(size_t)
{
  if(!B4cCalorHitAllocator)
    B4cCalorHitAllocator = new G4Allocator<B4cCalorHit>;
  void *hit;
  hit = (void *) B4cCalorHitAllocator->MallocSingle();
  return hit;
}

inline void B4cCalorHit::operator delete(void *hit)
{
  if(!B4cCalorHitAllocator)
    B4cCalorHitAllocator = new G4Allocator<B4cCalorHit>;
  B4cCalorHitAllocator->FreeSingle((B4cCalorHit*) hit);
}

inline void B4cCalorHit::Add(G4double de, G4double dl) {
  fEdep += de; 
  fTrackLength += dl;
}

inline G4double B4cCalorHit::GetEdep() const { 
  return fEdep; 
}

inline G4double B4cCalorHit::GetTrackLength() const { 
  return fTrackLength; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
