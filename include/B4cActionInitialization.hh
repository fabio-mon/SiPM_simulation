
#ifndef B4cActionInitialization_h
#define B4cActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "G4String.hh"
#include <string>

/// Action initialization class.
///

class B4cActionInitialization : public G4VUserActionInitialization
{
  public:
    B4cActionInitialization(std::string configFileName);
    virtual ~B4cActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;
  private:
    G4String fCutOption;
    std::string fconfigFileName;
};

#endif

    
