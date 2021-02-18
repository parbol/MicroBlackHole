#include "ConfigurationGeometry.hh"
#include <json/json.h>
#include <json/value.h>
//#include "jsoncpp.cpp"



//----------------------------------------------------------------------//
// Constructor                                                          //
//----------------------------------------------------------------------//
ConfigurationGeometry::ConfigurationGeometry(G4String file) {

    xSizeLayer = 0;
    ySizeLayer = 0;
    zSizeLayer = 0;
    xOffsetLayer = 0;
    yOffsetLayer = 0;
    zOffsetLayer = 0;

    Json::Value root;
    Json::Reader reader;

    //We open the JSON file, read it and put it on a string
    std::stringstream  filecontent;
    std::ifstream infile(file.c_str());

    if(!(infile.good())) {
        G4cerr << "\033[1;31m" << "Error opening geometry file: " + file << "\033[0m" << G4endl;
        goodGeometry = false;
        return;
    }

    std::string currline;
    while(getline(infile,currline)) {
        filecontent<<currline;
    }

    infile.close();

    std::string FileContent=filecontent.str();

    //Parsing of the JSON file
    bool parsingSuccessful = reader.parse( FileContent, root );
    if ( !parsingSuccessful ) {
        G4cerr << "\033[1;31m" << "Error parsing file: " + file << "\033[0m" << G4endl;
        goodGeometry = false;
        return;
    }

    if( root.size() > 0 ) {

        G4double xSize = atof(root["TheWorld"]["xSizeWorld"].asString().c_str());
        G4double ySize = atof(root["TheWorld"]["ySizeWorld"].asString().c_str());
        G4double zSize = atof(root["TheWorld"]["zSizeWorld"].asString().c_str());
        G4double zOffsetCRY_ = atof(root["TheWorld"]["zOffsetCRY"].asString().c_str());
        G4double sizeBoxCRY_ = atof(root["TheWorld"]["sizeBoxCRY"].asString().c_str());

        if(xSize <= 0 || ySize <= 0|| zSize <= 0) {
            G4cerr << "\033[1;31m" << "The size of the Universe has been greater than 0" << "\033[0m" << G4endl;
            goodGeometry = false;
            return;
        }

        if(zOffsetCRY_ < 0) {
            G4cerr << "\033[1;31m" << "Cry should be producing muons above the surface" << "\033[0m" << G4endl;
            goodGeometry = false;
            return;
        }

        if(sizeBoxCRY_ <= 0) {
            G4cerr << "\033[1;31m" << "Cry should have a positive size of production" << "\033[0m" << G4endl;
            goodGeometry = false;
            return;
        }

        uniSizeX = xSize * CLHEP::cm;
        uniSizeY = ySize * CLHEP::cm;
        uniSizeZ = zSize * CLHEP::cm;
        zOffsetCRY = zOffsetCRY_ * CLHEP::cm;
        sizeBoxCRY = sizeBoxCRY_ * CLHEP::cm;

        if(!root["LeadLayer"].isNull()) {
            G4double xSizeLayer_ = atof(root["LeadLayer"]["xSize"].asString().c_str());
            G4double ySizeLayer_ = atof(root["LeadLayer"]["ySize"].asString().c_str());
            G4double zSizeLayer_ = atof(root["LeadLayer"]["zSize"].asString().c_str());
            G4String outerMaterial_ = root["LeadLayer"]["outerMaterial"].asString().c_str();
            if(xSize - xSizeLayer_ <= 0 || ySize - ySizeLayer_ <= 0|| zSize - zSizeLayer_ <= 0) {
                G4cerr << "\033[1;31m" << "The size of the layer has to fit the universe " << "\033[0m" << G4endl;
                goodGeometry = false;
                return;
            }
	    G4double x_offset_ = atof(root["LeadLayer"]["xOffset"].asString().c_str());
            G4double y_offset_ = atof(root["LeadLayer"]["yOffset"].asString().c_str());
            G4double z_offset_ = atof(root["LeadLayer"]["zOffset"].asString().c_str());
            if(xSize/2 - (abs(x_offset_) + xSizeLayer_/2) <= 0 ||ySize/2 - (abs(y_offset_) + ySizeLayer_/2) <= 0  ||zSize/2 - (abs(z_offset_) + zSizeLayer_/2) <= 0 ) {
                G4cerr << "\033[1;31m" << "The size and the offset of the layer has to fit the universe " << "\033[0m" << G4endl;
                goodGeometry = false;
                return;
            }
	    xSizeLayer = xSizeLayer_ * CLHEP::cm;
            ySizeLayer = ySizeLayer_ * CLHEP::cm;
            zSizeLayer = zSizeLayer_ * CLHEP::cm;
	    xOffsetLayer = x_offset_ * CLHEP::cm;
            yOffsetLayer = y_offset_ * CLHEP::cm;
            zOffsetLayer = z_offset_ * CLHEP::cm;
            outerMaterialLayer = outerMaterial_ ;
            
        }
        

        G4double xPos1 = atof(root["Detector1"]["xPosDetector"].asString().c_str()) * CLHEP::cm;
        G4double yPos1 = atof(root["Detector1"]["yPosDetector"].asString().c_str()) * CLHEP::cm;
        G4double zPos1 = atof(root["Detector1"]["zPosDetector"].asString().c_str()) * CLHEP::cm;
        G4double xDir1 = atof(root["Detector1"]["xDirDetector"].asString().c_str()) * CLHEP::degree;
        G4double yDir1 = atof(root["Detector1"]["yDirDetector"].asString().c_str()) * CLHEP::degree;
        G4double zDir1 = atof(root["Detector1"]["zDirDetector"].asString().c_str()) * CLHEP::degree;
        G4double xSize1 = atof(root["Detector1"]["xSizeDetector"].asString().c_str()) * CLHEP::cm;
        G4double ySize1 = atof(root["Detector1"]["ySizeDetector"].asString().c_str()) * CLHEP::cm;
        G4double zSize1 = atof(root["Detector1"]["zSizeDetector"].asString().c_str()) * CLHEP::cm;
        
        detector1 = new Detector(xPos1, yPos1, zPos1, xDir1, yDir1, zDir1, xSize1, ySize1, zSize1);

        const Json::Value Layer1 = root["Detector1"]["layers"];
        for(G4int icoll = 0; icoll < Layer1.size(); ++icoll) {
            G4double xPosLayer = atof(root["Detector1"]["layers"][icoll]["xPosLayer"].asString().c_str()) * CLHEP::cm;
            G4double yPosLayer = atof(root["Detector1"]["layers"][icoll]["yPosLayer"].asString().c_str()) * CLHEP::cm;
            G4double zPosLayer = atof(root["Detector1"]["layers"][icoll]["zPosLayer"].asString().c_str()) * CLHEP::cm;
            G4double xDirLayer = atof(root["Detector1"]["layers"][icoll]["xDirLayer"].asString().c_str()) * CLHEP::degree;
            G4double yDirLayer = atof(root["Detector1"]["layers"][icoll]["yDirLayer"].asString().c_str()) * CLHEP::degree;
            G4double zDirLayer = atof(root["Detector1"]["layers"][icoll]["zDirLayer"].asString().c_str()) * CLHEP::degree;
            G4double xSizeLayer_ = atof(root["Detector1"]["layers"][icoll]["xSizeLayer"].asString().c_str()) * CLHEP::cm;
            G4double ySizeLayer_ = atof(root["Detector1"]["layers"][icoll]["ySizeLayer"].asString().c_str()) * CLHEP::cm;
            G4double zSizeLayer_ = atof(root["Detector1"]["layers"][icoll]["zSizeLayer"].asString().c_str()) * CLHEP::cm;
            G4double efficiency = atof(root["Detector1"]["layers"][icoll]["efficiency"].asString().c_str());
            G4double uncertainty = atof(root["Detector1"]["layers"][icoll]["uncertainty"].asString().c_str()) * CLHEP::cm;
            if(efficiency < 0 || efficiency > 1) {
                G4cerr << "\033[1;31m" << "Chamber 1 efficiency is not in the range [0, 1]" << "\033[0m" << G4endl;
                goodGeometry = false;
                return;
            }
            std::vector<G4double> wires_layer;
            const Json::Value Wires = root["Detector1"]["layers"][icoll]["wires"];
            for(G4int jcoll = 0; jcoll < Wires.size(); ++jcoll) {
                wires_layer.push_back(atof(Wires[jcoll].asString().c_str()) * CLHEP::cm);
            }
            detector1->AddLayer(xPosLayer, yPosLayer, zPosLayer, xDirLayer, yDirLayer, zDirLayer, xSizeLayer_, ySizeLayer_, zSizeLayer_, efficiency, uncertainty, wires_layer);
        }

        
        G4double xPos2 = atof(root["Detector2"]["xPosDetector"].asString().c_str()) * CLHEP::cm;
        G4double yPos2 = atof(root["Detector2"]["yPosDetector"].asString().c_str()) * CLHEP::cm;
        G4double zPos2 = atof(root["Detector2"]["zPosDetector"].asString().c_str()) * CLHEP::cm;
        G4double xDir2 = atof(root["Detector2"]["xDirDetector"].asString().c_str()) * CLHEP::degree;
        G4double yDir2 = atof(root["Detector2"]["yDirDetector"].asString().c_str()) * CLHEP::degree;
        G4double zDir2 = atof(root["Detector2"]["zDirDetector"].asString().c_str()) * CLHEP::degree;
        G4double xSize2 = atof(root["Detector2"]["xSizeDetector"].asString().c_str()) * CLHEP::cm;
        G4double ySize2 = atof(root["Detector2"]["ySizeDetector"].asString().c_str()) * CLHEP::cm;
        G4double zSize2 = atof(root["Detector2"]["zSizeDetector"].asString().c_str()) * CLHEP::cm;
        
        detector2 = new Detector(xPos2, yPos2, zPos2, xDir2, yDir2, zDir2, xSize2, ySize2, zSize2);
        const Json::Value Layer2 = root["Detector2"]["layers"];
        for(G4int icoll = 0; icoll < Layer2.size(); ++icoll) {
            G4double xPosLayer = atof(root["Detector2"]["layers"][icoll]["xPosLayer"].asString().c_str()) * CLHEP::cm;
            G4double yPosLayer = atof(root["Detector2"]["layers"][icoll]["yPosLayer"].asString().c_str()) * CLHEP::cm;
            G4double zPosLayer = atof(root["Detector2"]["layers"][icoll]["zPosLayer"].asString().c_str()) * CLHEP::cm;
            G4double xDirLayer = atof(root["Detector2"]["layers"][icoll]["xDirLayer"].asString().c_str()) * CLHEP::degree;
            G4double yDirLayer = atof(root["Detector2"]["layers"][icoll]["yDirLayer"].asString().c_str()) * CLHEP::degree;
            G4double zDirLayer = atof(root["Detector2"]["layers"][icoll]["zDirLayer"].asString().c_str()) * CLHEP::degree;
            G4double xSizeLayer_ = atof(root["Detector2"]["layers"][icoll]["xSizeLayer"].asString().c_str()) * CLHEP::cm;
            G4double ySizeLayer_ = atof(root["Detector2"]["layers"][icoll]["ySizeLayer"].asString().c_str()) * CLHEP::cm;
            G4double zSizeLayer_ = atof(root["Detector2"]["layers"][icoll]["zSizeLayer"].asString().c_str()) * CLHEP::cm;
            G4double efficiency = atof(root["Detector2"]["layers"][icoll]["efficiency"].asString().c_str());
            G4double uncertainty = atof(root["Detector2"]["layers"][icoll]["uncertainty"].asString().c_str()) * CLHEP::cm;
            if(efficiency < 0 || efficiency > 1) {
                G4cerr << "\033[1;31m" << "Chamber 2 efficiency is not in the range [0, 1]" << "\033[0m" << G4endl;
                goodGeometry = false;
                return;
            }
            std::vector<G4double> wires_layer;
            const Json::Value Wires = root["Detector2"]["layers"][icoll]["wires"];
            for(G4int jcoll = 0; jcoll < Wires.size(); ++jcoll) {
                wires_layer.push_back(atof(Wires[jcoll].asString().c_str()) * CLHEP::cm);
            }
            detector2->AddLayer(xPosLayer, yPosLayer, zPosLayer, xDirLayer, yDirLayer, zDirLayer, xSizeLayer_, ySizeLayer_, zSizeLayer_, efficiency, uncertainty, wires_layer);
        }

    }

    goodGeometry = true;
    return;

}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Accesor to class information                                         //
//----------------------------------------------------------------------//
bool ConfigurationGeometry::isGood() {
    return goodGeometry;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Accesor to class information                                         //
//----------------------------------------------------------------------//
G4double ConfigurationGeometry::getSizeX() {
    return uniSizeX;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Accesor to class information                                         //
//----------------------------------------------------------------------//
G4double ConfigurationGeometry::getSizeY() {
    return uniSizeY;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Accesor to class information                                         //
//----------------------------------------------------------------------//
G4double ConfigurationGeometry::getSizeZ() {
    return uniSizeZ;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Accesor to class information                                         //
//----------------------------------------------------------------------//
G4double ConfigurationGeometry::getZOffsetCRY() {
    return zOffsetCRY;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Accesor to class information                                         //
//----------------------------------------------------------------------//
G4double ConfigurationGeometry::getSizeBoxCRY() {
    return sizeBoxCRY;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Accesor to class information                                         //
//----------------------------------------------------------------------//
G4double ConfigurationGeometry::getSizeXLayer() {
    return xSizeLayer;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Accesor to class information                                         //
//----------------------------------------------------------------------//
G4double ConfigurationGeometry::getSizeYLayer() {
    return ySizeLayer;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Accesor to class information                                         //
//----------------------------------------------------------------------//
G4double ConfigurationGeometry::getSizeZLayer() {
    return zSizeLayer;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Accesor to class information                                         //
//----------------------------------------------------------------------//
G4double ConfigurationGeometry::getXOffsetLayer() {
    return xOffsetLayer;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Accesor to class information                                         //
//----------------------------------------------------------------------//
G4double ConfigurationGeometry::getYOffsetLayer() {
    return yOffsetLayer;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Accesor to class information                                         //
//----------------------------------------------------------------------//
G4double ConfigurationGeometry::getZOffsetLayer() {
    return zOffsetLayer;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Accesor to class information                                         //
//----------------------------------------------------------------------//
G4String ConfigurationGeometry::getOuterMaterialLayer() {
    return outerMaterialLayer;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//



//----------------------------------------------------------------------//
// Accesor to class information                                         //
//----------------------------------------------------------------------//
Detector * ConfigurationGeometry::getDetector1() {
    return detector1;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Accesor to class information                                         //
//----------------------------------------------------------------------//
Detector * ConfigurationGeometry::getDetector2() {
    return detector2;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Print the class information                                          //
//----------------------------------------------------------------------//
void ConfigurationGeometry::Print() {

    G4cout << "\033[1;34m" << "---------------------------------------Geometry information-----------------------------------" << "\033[0m" << G4endl;
    G4cout << "\033[1;34m" << "The loaded geometry is as follows: " << G4endl;
    G4cout << "\033[1;34m" << "The world is contained in [" << -uniSizeX/2.0 << ", " << uniSizeX/2.0 << "]x[" << -uniSizeY/2.0 << ", " << uniSizeY/2.0 << "]x[" << -uniSizeZ/2.0 << ", " << uniSizeZ/2.0 << "]" << G4endl;
    G4cout << "\033[1;34m" << "The size of the layer is: " << xSizeLayer/CLHEP::cm << " " << ySizeLayer/CLHEP::cm << " " << zSizeLayer/CLHEP::cm << G4endl;
    G4cout << "\033[1;34m" << "The material of the layer is: " << outerMaterialLayer << G4endl;
    G4cout << "\033[1;34m" << "Detector 1 located" << G4endl;
    detector1->Print();
    G4cout << "\033[1;34m" << "Detector 2 located" << G4endl;
    detector2->Print();
    G4cout << "\033[0m" << G4endl;

}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//




