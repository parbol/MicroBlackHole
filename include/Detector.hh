#ifndef Detector_h
#define Detector_h 1

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ErrorMatrix.hh"


class Detector {

public:
    Detector(G4double, G4double, G4double, G4double, G4double, G4double, G4double, G4double, G4double);
    G4int AddLayer(G4double, G4double, G4double, G4double, G4double, G4double, G4double, G4double, G4double, G4double, G4double, std::vector<G4double>);
    G4RotationMatrix *getRot();
    G4ThreeVector     getPos();
    G4ThreeVector     getSizes();
    G4int             getNLayers();
    G4ThreeVector     toGlobal(G4ThreeVector);
    G4ThreeVector     toLocal(G4ThreeVector);
    G4ThreeVector     toGlobal(G4ThreeVector, G4int);
    G4ErrorMatrix     toGlobalStateVector(G4ErrorMatrix);
    G4ErrorMatrix     toGlobalStateCov(G4ErrorMatrix);
    G4ThreeVector     toLocal(G4ThreeVector, G4int);
    G4ThreeVector     getPosLayer(G4int);
    G4ThreeVector     getSizeLayer(G4int);
    G4ThreeVector     getRotsLayer(G4int);
    G4RotationMatrix *getRotLayer(G4int);
    G4double          getEffLayer(G4int);
    G4double          getUncertaintyLayer(G4int);
    std::vector< std::vector<G4double> > getWires();
    G4VPhysicalVolume * getVolume(G4int);
    void              setVolume(G4VPhysicalVolume *);
    void              Print();

private:
    G4RotationMatrix rot, invrot;
    G4ThreeVector pos;
    G4ThreeVector sizes;
    G4ThreeVector rots;
    G4ErrorMatrix *superMatrix;
    G4ErrorMatrix *superVector;
    std::vector<G4ThreeVector> PosLayer;
    std::vector<G4ThreeVector> SizeLayer;
    std::vector<G4ThreeVector> RotsLayer;
    std::vector<G4RotationMatrix> RotLayer;
    std::vector<G4RotationMatrix> InvRotLayer;
    std::vector< std::vector<G4double> > wires;
    std::vector<G4VPhysicalVolume *> volumes;
    std::vector<G4double> efficiency;
    std::vector<G4double> uncertainty;
};



#endif

