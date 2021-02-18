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

#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4StoppingPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4ionIonisation.hh"
#include "G4hIonisation.hh"
#include "G4HadronicAbsorptionBertini.hh"
#include "G4HadronicAbsorptionFritiof.hh"
#include "G4MuonMinusCapture.hh"


MicroBlackHolePhysics::MicroBlackHolePhysics(const G4String& nam)
  : G4VPhysicsConstructor(nam),
    fBH(0)
{
  fElCharge  = 1.0;
  fMass = 1.220910e19*GeV;
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
 

  ph->RegisterProcess(new G4ionIonisation(), fBH);
  ph->RegisterProcess(new G4hIonisation(), fBH);
  ph->RegisterProcess(new G4HadronicAbsorptionBertini(), fBH);
  ph->RegisterProcess(new G4HadronicAbsorptionFritiof(), fBH);

  //ph->RegisterProcess(new G4StepLimiter(), fBH);

}

void MicroBlackHolePhysics::SetElectricCharge(G4double val)
{
  fElCharge = val;
}

void MicroBlackHolePhysics::SetMass(G4double mass)
{
  fMass = mass;
}


