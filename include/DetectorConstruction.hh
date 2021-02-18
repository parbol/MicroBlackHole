#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"
#include "G4FieldManager.hh"
#include "ConfigurationGeometry.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Box.hh"
#include "G4SDManager.hh"


class G4VPhysicalVolume;
class G4Material;
class G4VSensitiveDetector;
class G4VisAttributes;
class DetectorConstMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction(ConfigurationGeometry *);
    virtual ~DetectorConstruction();

public:
    virtual G4VPhysicalVolume* Construct();

private:
    void ConstructMaterials();
    void DestroyMaterials();
    void DumpGeometricalTree(G4VPhysicalVolume*, G4int);
    G4LogicalVolume* getPhysicalChamber(G4String , Detector * , G4int, G4SDManager* );
    G4ThreeVector getMaterialColor(G4String);
    
private:
    DetectorConstMessenger* messenger;

    ConfigurationGeometry *myConf;

    //In the future we will move this to a map
    std::map<G4String, G4Material*> materials;



};

#endif

