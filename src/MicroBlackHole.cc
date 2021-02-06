#include "MicroBlackHole.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"


MicroBlackHole* MicroBlackHole::theMicroBlackHole = 0;


MicroBlackHole::MicroBlackHole(
       const G4String&     aName,        G4double            mass,
       G4double            width,        G4double            charge,   
       G4int               iSpin,        G4int               iParity,    
       G4int               iConjugation, G4int               iIsospin,   
       G4int               iIsospin3,    G4int               gParity,
       const G4String&     pType,        G4int               lepton,      
       G4int               baryon,       G4int               encoding,
       G4bool              stable,       G4double            lifetime,
       G4DecayTable        *decaytable)
 : G4ParticleDefinition( aName, mass, width, charge, iSpin, iParity,
           iConjugation, iIsospin, iIsospin3, gParity, pType,
           lepton, baryon, encoding, stable, lifetime, decaytable )
{}


MicroBlackHole::~MicroBlackHole()
{}

//     
//    Arguments for constructor are as follows
//               name             mass          width         charge
//             2*spin           parity  C-conjugation
//          2*Isospin       2*Isospin3       G-parity
//               type    lepton number  baryon number   PDG encoding
//             stable         lifetime    decay table 
//

MicroBlackHole* MicroBlackHole::MicroBlackHoleDefinition(G4double mass, G4double eCharge)
{    
  if(!theMicroBlackHole) {
    theMicroBlackHole = new MicroBlackHole(
       "microblackhole",         mass,       0.0*MeV,       eplus*eCharge, 
                0,               0,             0,          
                0,               0,             0,             
          "boson",               0,             0,           0,
             true,            -1.0,             0);
    
    
    G4cout << "Microblackhole is created: m(GeV)= " << theMicroBlackHole->GetPDGMass()/GeV 
           << " Qel= " << theMicroBlackHole->GetPDGCharge()/eplus
           << G4endl;
  }
  return theMicroBlackHole;
}


MicroBlackHole* MicroBlackHole::MicroBH()
{    
  if(!theMicroBlackHole) { theMicroBlackHole = MicroBlackHoleDefinition(); }
  return theMicroBlackHole;
} 


