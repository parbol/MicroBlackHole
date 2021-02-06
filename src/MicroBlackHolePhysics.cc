#include "MicroBlackHolePhysics.hh"

#include "MicroBlackHole.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"

#include "G4StepLimiter.hh"
#include "G4Transportation.hh"
#include "G4hMultipleScattering.hh"
#include "G4mplIonisation.hh"
#include "G4mplIonisationWithDeltaModel.hh"
#include "G4hhIonisation.hh"
#include "G4hIonisation.hh"

#include "G4PhysicsListHelper.hh"

#include "G4BuilderType.hh"
#include "G4SystemOfUnits.hh"


MicroBlackHolePhysics::MicroBlackHolePhysics(const G4String& nam)
  : G4VPhysicsConstructor(nam),
    fBH(0)
{
  fElCharge  = 1.0;
  fMass = 100.*GeV;
  SetPhysicsType(bUnknown);
}

MicroBlackHolePhysics::~MicroBlackHolePhysics()
{
}


void MicroBlackHolePhysics::ConstructParticle()
{
  fBH = MicroBlackHole::MicroBlackHoleDefinition(fMass, fElCharge);
}


void MicroBlackHolePhysics::ConstructProcess()
{
  if(verboseLevel > 0) {
    G4cout << "MicroBlackHolePhysics::ConstructProcess" << G4endl;
  }
  
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  G4ProcessManager* pmanager = fBH->GetProcessManager();
  
  G4double emin = fMass/20000.;
  if(emin < keV) { emin = keV; }
  G4double emax = std::max(10.*TeV, fMass*100);
  G4int nbin = G4lrint(10*std::log10(emax/emin));

  if(fBH->GetPDGCharge() != 0.0) {
    G4hIonisation* hhioni = new G4hIonisation();
    hhioni->SetDEDXBinning(nbin);
    hhioni->SetMinKinEnergy(emin);
    hhioni->SetMaxKinEnergy(emax);
    ph->RegisterProcess(hhioni, fBH);
  }
  
  ph->RegisterProcess(new G4StepLimiter(), fBH);
}

void MicroBlackHolePhysics::SetElectricCharge(G4double val)
{
  fElCharge = val;
}

void MicroBlackHolePhysics::SetMass(G4double mass)
{
  fMass = mass;
}


