#ifndef ConfigurationGeometry_h
#define ConfigurationGeometry_h 1


#include <vector>
#include <iostream>
#include <fstream>
#include "globals.hh"

#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include "assert.h"

#include "Detector.hh"


class ConfigurationGeometry {

public:
    ConfigurationGeometry(G4String);
    bool isGood();
    G4double getSizeX();
    G4double getSizeY();
    G4double getSizeZ();
    G4double getZOffsetCRY();
    G4double getSizeBoxCRY();

    G4double getSizeXLayer();
    G4double getSizeYLayer();
    G4double getSizeZLayer();
    G4double getXOffsetLayer();
    G4double getYOffsetLayer();
    G4double getZOffsetLayer();
    G4String getOuterMaterialLayer();


    Detector *getDetector1();
    Detector *getDetector2();

    void Print();

private:

    G4double uniSizeX, uniSizeY, uniSizeZ;
    G4double zOffsetCRY, sizeBoxCRY;
    G4double xSizeLayer, ySizeLayer, zSizeLayer, xOffsetLayer, yOffsetLayer, zOffsetLayer;
    G4String outerMaterialLayer;

    Detector *detector1, *detector2;
    bool goodGeometry;

};



#endif

