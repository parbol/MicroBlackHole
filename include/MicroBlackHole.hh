#ifndef MicroBlackHole_h
#define MicroBlackHole_h 1

#include "G4ParticleDefinition.hh"
#include "globals.hh"

#include "CLHEP/Units/SystemOfUnits.h"


class MicroBlackHole : public G4ParticleDefinition
{
private:

  static MicroBlackHole*  theMicroBlackHole;

  MicroBlackHole(
       const G4String&     aName,        G4double            mass,
       G4double            width,        G4double            charge,   
       G4int               iSpin,        G4int               iParity,    
       G4int               iConjugation, G4int               iIsospin,   
       G4int               iIsospin3,    G4int               gParity,
       const G4String&     pType,        G4int               lepton,      
       G4int               baryon,       G4int               encoding,
       G4bool              stable,       G4double            lifetime,
       G4DecayTable        *decaytable );

  virtual ~MicroBlackHole();

public: 
  
  static MicroBlackHole* MicroBlackHoleDefinition(G4double mass = 100.*CLHEP::GeV, 
                                        G4double elCharge  = 0.0);

  static MicroBlackHole* MicroBH();

private:

};

#endif
