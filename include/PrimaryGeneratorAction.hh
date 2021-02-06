//----------------------------------------------------------------------//
// ___  ___                    _____           _                        //
// |  \/  |                   /  ___|         | |                       //
// | .  . |_   _  ___  _ __   \ `--. _   _ ___| |_ ___ _ __ ___  ___    //
// | |\/| | | | |/ _ \| '_ \   `--. \ | | / __| __/ _ \ '_ ` _ \/ __|   //
// | |  | | |_| | (_) | | | | /\__/ / |_| \__ \ ||  __/ | | | | \__ \   //
// \_|  |_/\__,_|\___/|_| |_| \____/ \__, |___/\__\___|_| |_| |_|___/   //
//                                    __/ |                             //
//----------------------------------------------------------------------//
// A project by: C. Diez, P. Gomez and P. Martinez                      //
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//
// PrimaryGeneratorAction.hh                                            //
//----------------------------------------------------------------------//
// This class handles the generation of events by calling CRY.          //
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//
#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "G4DataVector.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleGun.hh"
#include "Randomize.hh"
#include "globals.hh"
#include "ConfigurationGeometry.hh"
#include "CRYSetup.h"
#include "CRYGenerator.h"
#include "CRYParticle.h"
#include "CRYUtils.h"
//include "CRYGenerator.h"
//include "CRYParticle.h"
//include "CRYUtils.h"
#include "vector"
#include "RNGWrapper.hh"
#include "PrimaryGeneratorMessenger.hh"

class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction(ConfigurationGeometry *, G4String, G4long);
    ~PrimaryGeneratorAction();

public:
    void GeneratePrimaries(G4Event* anEvent);
    void InputCRY();
    void UpdateCRY(std::string* MessInput);
    void CRYFromFile(G4String newValue);
    G4double timeSimulated;

private:
    std::vector<CRYParticle*> *vect; // vector of generated particles
    G4ParticleTable* particleTable;
    G4ParticleGun* particleGun;
    CRYGenerator* gen;
    G4int InputState;
    PrimaryGeneratorMessenger* gunMessenger;
    ConfigurationGeometry *myGeom;
    G4long randomSeed;

};

#endif
