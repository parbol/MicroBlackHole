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
// DirftChamberLayer.hh                                                 //
//----------------------------------------------------------------------//
// This class implements a drift chamber layer detector.                //
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


#ifndef DriftChamberLayer_h
#define DriftChamberLayer_h 1

#include "G4VSensitiveDetector.hh"
#include "DriftChamberLayerHit.hh"
#include "Detector.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class DriftChamberLayer : public G4VSensitiveDetector
{

public:
    DriftChamberLayer(G4String name, G4String collection);
    virtual ~DriftChamberLayer();

    virtual void Initialize(G4HCofThisEvent*HCE);
    virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
    G4bool ProcessHits_old(G4Step*aStep,G4TouchableHistory*ROhist);
    virtual void SetStructure(Detector *);
    virtual void EndOfEvent(G4HCofThisEvent*HCE);
    std::vector<G4double> getWirePos(G4double, G4int);
    std::vector<G4double> getProbabilisticWirePos(G4double, G4ThreeVector, G4int);

private:
    DriftChamberLayerHitsCollection * hitsCollection;
    G4ThreeVector theLocalPosError;
    G4int HCID;
    Detector *det;
};




#endif

