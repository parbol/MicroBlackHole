#ifndef MicroBlackHolePhysics_h
#define MicroBlackHolePhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class MicroBlackHolePhysicsMessenger;
class MicroBlackHole;


class MicroBlackHolePhysics : public G4VPhysicsConstructor
{
public:

  MicroBlackHolePhysics(const G4String& nam = "Monopole Physics");

  ~MicroBlackHolePhysics();

  // This method is dummy for physics
  virtual void ConstructParticle();

  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type
  virtual void ConstructProcess();

  void SetElectricCharge(G4double);
  void SetMass(G4double);

private:

  // hide assignment operator
  MicroBlackHolePhysics & operator=(const MicroBlackHolePhysics &right);
  MicroBlackHolePhysics(const MicroBlackHolePhysics&);

  G4double    fElCharge;
  G4double    fMass;

  MicroBlackHole* fBH;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

