#include "DetectorConstruction.hh"
#include "DriftChamberLayer.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"

#include "G4Para.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"

#include "G4ios.hh"


#include "G4PVReplica.hh"

#include "G4SubtractionSolid.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include <algorithm>
#include <math.h>
#include <sys/time.h>


//Check volumes
G4bool checkSurface = false;
//----------------------------------------------------------------------//
// Constructor                                                          //
//----------------------------------------------------------------------//
DetectorConstruction::DetectorConstruction(ConfigurationGeometry *w) {
    myConf = w;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Destructor                                                           //
//----------------------------------------------------------------------//
DetectorConstruction::~DetectorConstruction() {

    DestroyMaterials();

}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//

//----------------------------------------------------------------------//
// Creates all the geometrical structures                               //
//----------------------------------------------------------------------//
G4VPhysicalVolume* DetectorConstruction::Construct() {

    if(!checkSurface) std::cout << " **** WARNING : Deactivate checking overlaps for volumes. Activate at DetectorConstruction.cc setting checkSurface to true, for a better solid location control." << std::endl;
    //Building the materials
    ConstructMaterials();

    //Printing the geometry
    myConf->Print();

    //Manager of objects in memory
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    G4String SDname;

    //Creating the world
    G4VSolid* worldSolid = new G4Box("worldBox", myConf->getSizeX()/2.0, myConf->getSizeY()/2.0, myConf->getSizeZ()/2.0 );
    G4LogicalVolume* worldLogical = new G4LogicalVolume(worldSolid, materials["air"], "worldLogical",0,0,0);
    G4VPhysicalVolume* worldPhysical = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), worldLogical, "worldPhysical", 0, false, 0, checkSurface);
    G4VisAttributes *worldVisAtt = new G4VisAttributes(G4Colour(255.0,6.0,0.0));
    worldVisAtt->SetVisibility(false);
    worldLogical->SetVisAttributes(worldVisAtt);


        
    G4ThreeVector color_vect = getMaterialColor(myConf->getOuterMaterialLayer());
    G4double red = color_vect.x();
    G4double green = color_vect.y();
    G4double blue = color_vect.z();
    G4double alpha = 1.;
        
    G4Colour* color1 = new G4Colour(red, green, blue, alpha);
    G4VisAttributes *attr_leadLayer = new G4VisAttributes(G4Colour::Magenta());
    attr_leadLayer->SetVisibility(true);        
    attr_leadLayer->SetColor(*color1);
        
        
    G4Box * leadLayer = new G4Box("leadLayer", myConf->getSizeXLayer()/2.0, myConf->getSizeYLayer()/2.0, myConf->getSizeZLayer()/2.0);
    G4LogicalVolume *leadLayerLogical = new G4LogicalVolume(leadLayer, materials[myConf->getOuterMaterialLayer()], "leadLayerLogical", 0, 0, 0);
    new G4PVPlacement(0, G4ThreeVector(myConf->getXOffsetLayer(), myConf->getYOffsetLayer(), myConf->getZOffsetLayer()), leadLayerLogical, "leadLayerPhysical", worldLogical, false, 0, checkSurface);
    leadLayerLogical->SetVisAttributes(attr_leadLayer);
        


    //Chamber attr
    G4VisAttributes *attr_chamber = new G4VisAttributes(G4Colour::Magenta());
    G4Colour* chamber_color = new G4Colour(0, 0, 1, 0.5);
    attr_chamber->SetVisibility(false);
    attr_chamber->SetColor(*chamber_color);
    
    //Chambers and layers definition    
    G4LogicalVolume* layers_det1 [myConf->getDetector1()->getNLayers()];
    for(G4int c_layer = 0; c_layer < myConf->getDetector1()->getNLayers(); c_layer++) {
        layers_det1[c_layer] = getPhysicalChamber("Detector1", myConf->getDetector1(), c_layer, SDman);
    }
	
    G4LogicalVolume* layers_det2 [myConf->getDetector2()->getNLayers()];
    for(G4int c_layer = 0; c_layer < myConf->getDetector2()->getNLayers(); c_layer++) {    
            layers_det2[c_layer] = getPhysicalChamber("Detector2", myConf->getDetector2(), c_layer, SDman);
    }
    
    G4VSolid* chamber1 = new G4Box("Chamber1", myConf->getDetector1()->getSizes().x()/2.0, myConf->getDetector1()->getSizes().y()/2.0, myConf->getDetector1()->getSizes().z()/2.0);
    G4LogicalVolume* chamber1Logical = new G4LogicalVolume(chamber1, materials["air"], "Chamber1Log", 0, 0, 0);
    chamber1Logical->SetVisAttributes(attr_chamber);
    
    G4VSolid* chamber2 = new G4Box("Chamber2", myConf->getDetector2()->getSizes().x()/2.0, myConf->getDetector2()->getSizes().y()/2.0, myConf->getDetector2()->getSizes().z()/2.0);
    G4LogicalVolume* chamber2Logical = new G4LogicalVolume(chamber2, materials["air"], "Chamber2Log", 0, 0, 0);
    chamber2Logical->SetVisAttributes(attr_chamber);
    
    for(G4int c_layer = 0; c_layer < myConf->getDetector1()->getNLayers(); c_layer++) {

        //G4VPhysicalVolume *vol = new G4PVPlacement(myConf->getDetector1()->getRotLayer(c_layer), myConf->getDetector1()->getPosLayer(c_layer), layerLogical1, "p_layer1", chamber1Logical, false, c_layer, true);
        G4VPhysicalVolume *vol = new G4PVPlacement(myConf->getDetector1()->getRotLayer(c_layer), myConf->getDetector1()->getPosLayer(c_layer), layers_det1[c_layer], "p_layer1"+ std::to_string(c_layer), chamber1Logical, false, c_layer, checkSurface);
        myConf->getDetector1()->setVolume(vol);

    }

    for(G4int c_layer = 0; c_layer < myConf->getDetector2()->getNLayers(); c_layer++) {

        //G4VPhysicalVolume *vol = new G4PVPlacement(myConf->getDetector2()->getRotLayer(c_layer), myConf->getDetector2()->getPosLayer(c_layer), layerLogical2, "p_layer2", chamber2Logical, false, c_layer, true);
        G4VPhysicalVolume *vol = new G4PVPlacement(myConf->getDetector2()->getRotLayer(c_layer), myConf->getDetector2()->getPosLayer(c_layer), layers_det2[c_layer], "p_layer2", chamber2Logical, false, c_layer, checkSurface);
        myConf->getDetector2()->setVolume(vol);

    }


    //Careful to the units
    new G4PVPlacement(myConf->getDetector1()->getRot(), myConf->getDetector1()->getPos(), chamber1Logical, "chamber1_phys", worldLogical, false, 0, checkSurface);
    new G4PVPlacement(myConf->getDetector2()->getRot(), myConf->getDetector2()->getPos(), chamber2Logical, "chamber2_phys", worldLogical, false, 0, checkSurface);

    DumpGeometricalTree(worldPhysical, 3);


    return worldPhysical;

}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Creates a chamber with a sensitive detection volumen
//----------------------------------------------------------------------//
G4LogicalVolume*  DetectorConstruction::getPhysicalChamber(G4String name, Detector *myDetector, G4int c_layer, G4SDManager* SDman) {

    
    //Variables de caracterización de los objetos que forman el layer    
    G4ThreeVector layer_pos = myDetector->getPosLayer(c_layer);
    G4double glassFiber_thickness = 1.6 * CLHEP::mm;
    G4double copperLayer_thickness = 0.035 * CLHEP::mm;
    G4double wire_radius = 0.025 / 2. * CLHEP::mm;

    G4double profile_thickness = 1.3 * CLHEP::mm;
    G4double profile_height = 1 * CLHEP::cm;
    G4double profile_width = 2 * CLHEP::cm;
    G4double profile_length_x = myDetector->getSizeLayer(c_layer).x() - 2 * profile_width;
    G4double profile_length_y = myDetector->getSizeLayer(c_layer).y() - 2 * profile_width;
    
    G4double cathod_gap = 10 * CLHEP::mm;
    
    G4double total_height = cathod_gap + 2*glassFiber_thickness + 2*copperLayer_thickness + profile_height;
        
    G4double delta_z_copper =  cathod_gap/2. + copperLayer_thickness/2.0;
    G4double delta_z_glassfiber = cathod_gap/2. + copperLayer_thickness + glassFiber_thickness/2.0;
    G4double delta_z_profile =  cathod_gap/2. + copperLayer_thickness + glassFiber_thickness + profile_height/2.0;
    // altura del volumen formado por la fibra de PC, el cobre y la distancia entre placas
    G4double sandwich_height = cathod_gap + 2*copperLayer_thickness + 2*glassFiber_thickness;
    
    
    
    // DEFINICION DE ATRIBUTOS
    // Atributos contenedor layer
    G4VisAttributes *attr_layer = new G4VisAttributes(G4Colour::Magenta());
    G4Colour* layer_color = new G4Colour(53./255., 81./255., 92./255., 0.25);
    attr_layer->SetVisibility(false);
    attr_layer->SetColor(*layer_color);
            
    // Atributos fibra de vidrio
    G4VisAttributes *attr_glassfiber = new G4VisAttributes(G4Colour::Magenta());
    G4Colour* glassfiber_color = new G4Colour(0./255., 187./255., 45./255., 1);
    attr_glassfiber->SetVisibility(true);
    attr_glassfiber->SetColor(*glassfiber_color);
    
    // Atributos cobre
    G4VisAttributes *attr_copper = new G4VisAttributes(G4Colour::Magenta());
    G4Colour* copper_color = new G4Colour(203./255., 109./255., 81./255., 1);
    attr_copper->SetVisibility(true);
    attr_copper->SetColor(*copper_color);
    
    // Atributos argon
    G4VisAttributes *attr_argon = new G4VisAttributes(G4Colour::Magenta());
    G4Colour* argon_color = new G4Colour(255./255., 34./255., 0./255., .5);
    attr_argon->SetVisibility(true);
    attr_argon->SetColor(*argon_color);
    
    // Atributos wolframio
    G4VisAttributes *attr_wolframio = new G4VisAttributes(G4Colour::Magenta());
    attr_wolframio->SetVisibility(true);
    
    
    // Atributos perfil aluminio
    G4VisAttributes *attr_profile = new G4VisAttributes(G4Colour::Magenta());
    G4Colour* profile_color = new G4Colour(143./255., 143./255., 143./255., 1);
    attr_profile->SetVisibility(true);
    attr_profile->SetColor(*profile_color);
    
    //____________________________________________________________
    
    
    
    // DEFINICION DE LOS ELEMENTOS FISICOS Y LOGICOS DE LA CAMARA
    // Contenedor de la cámara
    G4String myName = name + "_layer_" + std::to_string(c_layer);
    //Defino el volumen "layer" que englobará todos los elementos.
    G4VSolid* layer = new G4Box(myName, myDetector->getSizeLayer(c_layer).x()/2.0, myDetector->getSizeLayer(c_layer).y()/2.0, total_height/2.);
    G4LogicalVolume* layerLogical = new G4LogicalVolume(layer, materials["air"], myName+"_Log", 0, 0, 0);
    layerLogical->SetVisAttributes(attr_layer);
    
    // Perfiles
    myName = name + "_profile_" + std::to_string(c_layer);
    G4VSolid* layer_profile_outer = new G4Box("outer", profile_length_x/2.0, profile_width/2.0, profile_height/2.0);
    G4VSolid* layer_profile_inner = new G4Box("inner", profile_length_x/2.0, (profile_width-profile_thickness)/2.0, (profile_height-profile_thickness)/2.0);
    G4SubtractionSolid* layer_profile = new G4SubtractionSolid(myName,layer_profile_outer,layer_profile_inner);
    G4LogicalVolume* layerLogical_profile = new G4LogicalVolume(layer_profile, materials["glassfiber"], myName+"_Log", 0, 0, 0);
    layerLogical_profile->SetVisAttributes(attr_profile);
    
    // Perfiles un poco más largos que se colocan en el borde de la cámara
    myName = name + "_profileLarge_top_" + std::to_string(c_layer);
    G4VSolid* layer_profileLarge_outer = new G4Box("outer", profile_width/2.0, (profile_length_y + 2 * profile_width) / 2.0, profile_height/2.0);
    G4VSolid* layer_profileLarge_inner = new G4Box("inner", (profile_width-profile_thickness)/2.0, (profile_length_y + 2 * profile_width)/2.0, (profile_height-profile_thickness)/2.0);
    G4SubtractionSolid* layer_profileLarge = new G4SubtractionSolid(myName,layer_profileLarge_outer,layer_profileLarge_inner);
    G4LogicalVolume* layerLogical_profileLarge = new G4LogicalVolume(layer_profileLarge, materials["glassfiber"], myName+"_Log", 0, 0, 0);
    layerLogical_profileLarge->SetVisAttributes(attr_profile);
    
    
    // Capa de Fibra de vidrio
    myName = name + "_glassFiber_top_" + std::to_string(c_layer);
    G4VSolid* layer_glassFiber_top = new G4Box(myName, myDetector->getSizeLayer(c_layer).x()/2.0, myDetector->getSizeLayer(c_layer).y()/2.0, glassFiber_thickness/2.0);
    G4LogicalVolume* layerLogical_glassFiber_top = new G4LogicalVolume(layer_glassFiber_top, materials["glassfiber"], myName+"_Log", 0, 0, 0);
    layerLogical_glassFiber_top->SetVisAttributes(attr_glassfiber);
    
    // Capa de Cobre
    myName = name + "_copper_top_" + std::to_string(c_layer);
    G4VSolid* layer_copper_top = new G4Box(myName, myDetector->getSizeLayer(c_layer).x()/2.0, myDetector->getSizeLayer(c_layer).y()/2.0, copperLayer_thickness/2.0);
    G4LogicalVolume* layerLogical_copper_top = new G4LogicalVolume(layer_copper_top, materials["copper"], myName+"_Log", 0, 0, 0);
    layerLogical_copper_top->SetVisAttributes(attr_copper);
    
    
    // Gas argon
    myName = name + "_argon_" + std::to_string(c_layer);
    //Defino el volumen "layer" que englobará todos los elementos.
    G4VSolid* layer_argon = new G4Box(myName, myDetector->getSizeLayer(c_layer).x()/2.0, myDetector->getSizeLayer(c_layer).y()/2.0, cathod_gap/2.);
    G4LogicalVolume* layerLogical_argon = new G4LogicalVolume(layer_argon, materials["argon"], myName+"_Log", 0, 0, 0);
    layerLogical_argon->SetVisAttributes(attr_argon);    
    
    // Capa hilos
    std::vector< std::vector<G4double> > myWires = myDetector->getWires();
    G4RotationMatrix * rot_wires = new G4RotationMatrix;
    rot_wires->rotateX(M_PI/2*rad);
    
    myName = name + "_wire_" + std::to_string(c_layer) + "_";
    G4Tubs* wire = new G4Tubs(myName, 0.0, wire_radius, myDetector->getSizeLayer(c_layer).y()/2.0, 0 * CLHEP::deg, 360 * CLHEP::deg);
    G4LogicalVolume* layerLogical_wire = new G4LogicalVolume(wire, materials["wolframio"], myName+"_Log", 0, 0, 0);
    layerLogical_wire->SetVisAttributes(attr_wolframio);
    
    
     
            
    // POSICIONAMIENTO ESPACIAL
    
    // posciono los volumenes en las cámaras superiores del conjunto XY
    // las cámaras pares tendrán perfiles en su lado superior
    if(c_layer%2 == 0){
        // Poscionamiento de los perfiles
        
        //Perfil borde derecho
        myName = name + "_profile_top_p_rightEdge_" + std::to_string(c_layer);
        new G4PVPlacement(0, 
            G4ThreeVector(0,myDetector->getSizeLayer(c_layer).y()*(1./2.)-profile_width/2.0, sandwich_height/2.),
            layerLogical_profile, 
            myName, layerLogical, false, c_layer, checkSurface);
        //Perfil borde izquierdo
        myName = name + "_profile_top_p_leftEdge_" + std::to_string(c_layer);
        new G4PVPlacement(0, 
            G4ThreeVector(0,myDetector->getSizeLayer(c_layer).y()*(-1./2.)+profile_width/2.0, sandwich_height/2.),
            layerLogical_profile, 
            myName, layerLogical, false, c_layer, checkSurface);
        //Perfil borde superior
        myName = name + "_profile_top_p_topEdge_" + std::to_string(c_layer);
        new G4PVPlacement(0, 
            G4ThreeVector(myDetector->getSizeLayer(c_layer).x()*(1./2.)-profile_width/2.0, 0, sandwich_height/2.),
            layerLogical_profileLarge, 
            myName, layerLogical, false, c_layer, checkSurface);
        //Perfil borde inferior
        myName = name + "_profile_top_p_bottomEdge_" + std::to_string(c_layer);
        new G4PVPlacement(0, 
            G4ThreeVector(myDetector->getSizeLayer(c_layer).x()*(-1./2.)+profile_width/2.0, 0, sandwich_height/2.),
            layerLogical_profileLarge, 
            myName, layerLogical, false, c_layer, checkSurface);
        //Perfiles en el interior del marco exterior. 
        for(int prof = 1; prof<4;prof++){
            myName = name + "_profile_top_p_" + std::to_string(c_layer) + "_" + std::to_string(prof);
            new G4PVPlacement(0, 
                G4ThreeVector(0, myDetector->getSizeLayer(c_layer).y()*(-1./2. + prof * 1./4.), sandwich_height/2.),
                layerLogical_profile, 
                myName, layerLogical, false, c_layer, checkSurface);
            
        }
        
        
        // posciono la fibra de vidrio de la capa superior
        myName = name + "_glassFiber_top_p_" + std::to_string(c_layer);
        new G4PVPlacement(0, 
            G4ThreeVector(0,0,-profile_height/2.+delta_z_glassfiber),
            layerLogical_glassFiber_top, 
            myName, layerLogical, false, c_layer, checkSurface);

        // posiciono el cobre de la capa superior
        myName = name + "_copper_top_p_" + std::to_string(c_layer);
        new G4PVPlacement(0, 
            G4ThreeVector(0,0, -profile_height/2. + delta_z_copper),
            layerLogical_copper_top, 
            myName, layerLogical, false, c_layer, checkSurface);
        
              
        //Posiciono los hilos dentro del argon
        for (int w = 0;w < myWires[c_layer].size();w++){    
//        for (int w = 0;w < 5;w++){    
            myName = name + "_wire_p" + std::to_string(c_layer) + "_"+  std::to_string(w);
            new G4PVPlacement(rot_wires, 
            G4ThreeVector(myWires[c_layer][w] * CLHEP::mm, 0, 0.),
                layerLogical_wire, 
                myName, layerLogical_argon, false, c_layer, checkSurface);
        }
        
        // posiciono el argon de la cámara
        myName = name + "_argon_p_" + std::to_string(c_layer);
        new G4PVPlacement(0, 
            G4ThreeVector(0,0, -profile_height/2.),
            layerLogical_argon, 
            myName, layerLogical, false, c_layer, checkSurface);
        

        // posiciono la fibra de vidrio de la capa inferior
        myName = name + "_glassFiber_bottom_p_" + std::to_string(c_layer);
        new G4PVPlacement(0, 
            G4ThreeVector(0,0,  -profile_height/2. -  delta_z_glassfiber),
            layerLogical_glassFiber_top, 
            myName, layerLogical, false, c_layer, checkSurface);

        // posiciono el cobre de la capa inferior
        myName = name + "_copper_bottom_p_" + std::to_string(c_layer);
        new G4PVPlacement(0, 
            G4ThreeVector(0, 0, -profile_height/2. -  delta_z_copper),
            layerLogical_copper_top, 
            myName, layerLogical, false, c_layer, checkSurface);

        
    
    
    
    // posciono los volumenes en las cámaras inferiores del conjunto XY
    // las cámaras impares tendrán perfiles en su lado inferior
    }else{        
        
        // posciono la fibra de vidrio de la capa superior
        myName = name + "_glassFiber_top_p_" + std::to_string(c_layer);
        new G4PVPlacement(0, 
                G4ThreeVector(0,0,profile_height/2.+delta_z_glassfiber),
                layerLogical_glassFiber_top, 
                 myName, layerLogical, false, c_layer, checkSurface);

        // posiciono el cobre de la capa superior
        myName = name + "_copper_top_p_" + std::to_string(c_layer);
        new G4PVPlacement(0, 
            G4ThreeVector(0,0, profile_height/2. + delta_z_copper),
            layerLogical_copper_top, 
            myName, layerLogical, false, c_layer, checkSurface);
               
        //Posiciono los hilos dentro del argón       
        for (int w = 0;w < myWires[c_layer].size();w++){    
//        for (int w = 0;w < 5;w++){    
            myName = name + "_wire_p" + std::to_string(c_layer) + "_"+  std::to_string(w);
            new G4PVPlacement(rot_wires, 
                G4ThreeVector(myWires[c_layer][w] * CLHEP::mm, 0, 0.),
                layerLogical_wire, 
                myName, layerLogical_argon, false, c_layer, checkSurface);
        }
        
        // posiciono el argón de la cámara
        myName = name + "_argon_p_" + std::to_string(c_layer);
        new G4PVPlacement(0, 
            G4ThreeVector(0,0, profile_height/2.),
            layerLogical_argon, 
            myName, layerLogical, false, c_layer, checkSurface);
        

        // posiciono la fibra de vidrio de la capa inferior
        myName = name + "_glassFiber_bottom_p_" + std::to_string(c_layer);
        new G4PVPlacement(0, 
            G4ThreeVector(0,0,  profile_height/2. -  delta_z_glassfiber),
            layerLogical_glassFiber_top, 
            myName, layerLogical, false, c_layer, checkSurface);

        // posiciono el cobre de la capa inferior
        myName = name + "_copper_bottom_p_" + std::to_string(c_layer);
        new G4PVPlacement(0, 
            G4ThreeVector(0, 0, profile_height/2. -  delta_z_copper),
            layerLogical_copper_top, 
            myName, layerLogical, false, c_layer, checkSurface);

        // Poscionamiento de los perfiles
        
        //Perfil borde derecho
        myName = name + "_profile_bottom_p_rightEdge_" + std::to_string(c_layer);
        new G4PVPlacement(0, 
            G4ThreeVector(0,myDetector->getSizeLayer(c_layer).y()*(1./2.)-profile_width/2.0,-sandwich_height/2.),
            layerLogical_profile, 
            myName, layerLogical, false, c_layer, checkSurface);
        //Perfil borde izquierdo
        myName = name + "_profile_bottom_p_leftEdge_" + std::to_string(c_layer);
        new G4PVPlacement(0, 
            G4ThreeVector(0,myDetector->getSizeLayer(c_layer).y()*(-1./2.)+profile_width/2.0,-sandwich_height/2.),
            layerLogical_profile, 
            myName, layerLogical, false, c_layer, checkSurface);
        //Perfil borde superior
        myName = name + "_profile_bottom_p_topEdge_" + std::to_string(c_layer);
        new G4PVPlacement(0, 
            G4ThreeVector(myDetector->getSizeLayer(c_layer).x()*(1./2.)-profile_width/2.0, 0, -sandwich_height/2.),
            layerLogical_profileLarge, 
            myName, layerLogical, false, c_layer, checkSurface);
        //Perfil borde superior
        myName = name + "_profile_bottom_p_bottomEdge_" + std::to_string(c_layer);
        new G4PVPlacement(0, 
            G4ThreeVector(myDetector->getSizeLayer(c_layer).x()*(-1./2.)+profile_width/2.0, 0, -sandwich_height/2.),
            layerLogical_profileLarge, 
            myName, layerLogical, false, c_layer, checkSurface);
        //Perfiles en el interior del marco exterior. 
        for(int prof = 1; prof<4;prof++){
            myName = name + "_profile_bottom_p_" + std::to_string(c_layer) + "_" + std::to_string(prof);
            new G4PVPlacement(0, 
                G4ThreeVector(0,myDetector->getSizeLayer(c_layer).y()*(-1./2. + prof * 1./4.),-sandwich_height/2.),
                layerLogical_profile, 
                myName, layerLogical, false, c_layer, checkSurface);
            
        }
        
    }
    
    
    // Asociamos al argon la parte sensible del detector
    G4String SDname;
    DriftChamberLayer *theChamber = new DriftChamberLayer(SDname = name + "_Chamber_" + std::to_string(c_layer), "HitsCollection_"+name+"_"+std::to_string(c_layer));
    theChamber->SetStructure(myDetector);
    SDman->AddNewDetector(theChamber);
    layerLogical_argon->SetSensitiveDetector(theChamber);
    
    
    
    return layerLogical;
    
//    layers_det1[c_layer] = layerLogical1;
    

}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Creates a chamber with a sensitive detection volumen
//----------------------------------------------------------------------//
G4ThreeVector  DetectorConstruction::getMaterialColor(G4String material_name) {
    
    G4double red, green, blue;

    if(material_name == "steel") {
        red = 202./255.;
        green = 111./255.;
        blue = 30./255.;
    }else if(material_name == "iron") {
        red = 255./255.;
        green = 175./255.;
        blue = 51./255.;    
    }else if(material_name == "silicon") {
        red = 242./255.;
        green = 255./255.;
        blue = 51./255.;        
    }else if(material_name == "scrap") {
        red = 72./255.;
        green = 25./255.;
        blue = 5./255.;
    }else if(material_name == "air") {
        red = 52./255.;
        green = 12./255.;
        blue = 219./255.;
    }else if(material_name == "cascarilla") {
        red = 72./255.;
        green = 25./255.;
        blue = 5./255.;
    }else if(material_name == "water") {
        red = 72./255.;
        green = 25./255.;
        blue = 5./255.;
    }else if(material_name == "muscle") {
        red = 161./255.;
        green = 44./255.;
        blue = 44./255.;
    }else if(material_name == "boneCompact") {
        red = 227./255.;
        green = 218./255.;
        blue = 201./255.;
    }else if(material_name == "slag") {
        red = 168./255.;
        green = 167./255.;
        blue = 164./255.;
    }else if(material_name == "tailing_2.1") {
        red = 133./255.;
        green = 87./255.;
        blue = 35./255.;
    }else if(material_name == "tailing_1.9") {
        red = 146./255.;
        green = 104./255.;
        blue = 41./255.;
    }else if(material_name == "tailing_1.7") {
        red = 185./255.;
        green = 156./255.;
        blue = 107./255.;
    }else if(material_name == "tailing_1.5") {
        red = 234./255.;
        green = 182./255.;
        blue = 79./255.;
    }else if(material_name == "uranium") {
        red = 0./255.;
        green = 255./255.;
        blue = 0./255.;
    }else if(material_name == "aluminium") {
        red = 143./255.;
        green = 143./255.;
        blue = 143./255.;
    }else if(material_name == "argon") {
        red = 52./255.;
        green = 12./255.;
        blue = 219./255.;
    }else if(material_name == "lead") {
        red = 121./255.;
        green = 128./255.;
        blue = 129./255.;
    }else if(material_name == "copper") {
        red = 184./255.;
        green = 115./255.;
        blue = 51./255.;
    }else if(material_name == "copper") {
        red = 203./255.;
        green = 109./255.;
        blue = 81./255.;
    }else if(material_name == "concrete") {
        red = 210./255.;
        green = 209./255.;
        blue = 205./255.;
    }else if(material_name == "concrete_2.4") {
        red = 168./255.;
        green = 167./255.;
        blue = 164./255.;
    }else if(material_name == "carbonRefractory") {
        red = 25./255.;
        green = 25./255.;
        blue = 25./255.;
    }else if(material_name == "liquidSteel") {
        red = 211./255.;
        green = 34./255.;
        blue = 19./255.;
    }else if(material_name == "rockwool") {
        red = 243./255.;
        green = 255./255.;
        blue = 105./255.;
    }else if(material_name == "glassfiber") {
        red = 0./255.;
        green = 187./255.;
        blue = 45./255.;
    }else if(material_name == "carbonRefractory_electricFurnace") {
        red = 98./255.;
        green = 93./255.;
        blue = 93./255.;
    }else if(material_name == "liquidSilicon") {
        red = 211./255.;
        green = 89./255.;
        blue = 19./255.;
    }else if(material_name == "carbon_electricFurnace") {
        red = 92./255.;
        green = 97./255.;
        blue = 97./255.;
    }else if(material_name == "silicon_electricFurnace") {
        red = 172./255.;
        green = 191./255.;
        blue = 96./255.;
    }else if(material_name == "electrode") {
        red = 28./255.;
        green = 28./255.;
        blue = 28./255.;
    }else if(material_name == "quartzCarbon") {
        red = 211./255.;
        green = 34./255.;
        blue = 19./255.;
    }else {
        red = 255./255.;
        green = 0./255.;
        blue = 255./255.;
    }
    G4String electrode_material = "electrode";

    G4ThreeVector color_vect(red,green,blue);
    return color_vect;


}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//



//----------------------------------------------------------------------//
// Construct all the materials                                          //
//----------------------------------------------------------------------//
void DetectorConstruction::ConstructMaterials() {

    G4NistManager* man = G4NistManager::Instance();

    //G4double density, fractionmass;
    //G4int nElem;

    materials.insert(std::pair<G4String, G4Material *>("air", man->FindOrBuildMaterial("G4_AIR")));
    materials.insert(std::pair<G4String, G4Material *>("iron", man->FindOrBuildMaterial("G4_Fe")));
    materials.insert(std::pair<G4String, G4Material *>("uranium", man->FindOrBuildMaterial("G4_U")));
    materials.insert(std::pair<G4String, G4Material *>("aluminium", man->FindOrBuildMaterial("G4_Al")));
    materials.insert(std::pair<G4String, G4Material *>("argon", man->FindOrBuildMaterial("G4_Ar")));
    materials.insert(std::pair<G4String, G4Material *>("lead", man->FindOrBuildMaterial("G4_Pb")));
    materials.insert(std::pair<G4String, G4Material *>("silicon", man->FindOrBuildMaterial("G4_Si")));
    materials.insert(std::pair<G4String, G4Material *>("copper", man->FindOrBuildMaterial("G4_Cu")));
    materials.insert(std::pair<G4String, G4Material *>("wolframio", man->FindOrBuildMaterial("G4_W")));
    
    


    //G4Material* liquidSteel = new G4Material("liquidSteel", density = 6.00*g/cm3, nElem = 6);
    //liquidSteel->AddElement(man->FindOrBuildElement("G4_C"), fractionmass=0.1*perCent);
    //materials.insert(std::pair<G4String, G4Material *>("liquidSteel", liquidSteel));
    //


    G4double z, a, fractionmass, density;
    G4String name, symbol;
    G4int ncomponents;
    a = 14.01*g/mole;
    G4Element* elN  = new G4Element(name="Nitrogen",symbol="N", z= 7., a);

    a = 16.00*g/mole;
    G4Element* elO  = new G4Element(name="oxygen",symbol="O", z= 8., a);

    a = 1.01*g/mole;
    G4Element* elH = new G4Element(name="Hydrogen", symbol="H", z=1., a);
    a = 12.01*g/mole;
    G4Element* elC = new G4Element(name="Carbon", symbol="C", z=6., a);
    a = 22.99*g/mole;
    G4Element* elNa= new G4Element(name="Sodio", symbol="Na", z=11., a);
    a = 24.32*g/mole;
    G4Element* elMg= new G4Element(name="Magnesium", symbol="Mg", z=12., a);
    a = 26.98*g/mole;
    G4Element* elAl= new G4Element(name="aluminium", symbol="Al", z=13., a);
    a = 28*g/mole;
    G4Element* elSi= new G4Element(name="silicon", symbol="Si", z=14., a);
    a = 30.97*g/mole;
    G4Element* elP= new G4Element(name="Phosphorous", symbol="P", z=15., a);
    a = 32.06*g/mole;
    G4Element* elS= new G4Element(name="Sulphur", symbol="S", z=16., a);
    a = 35.45*g/mole;
    G4Element* elCl= new G4Element(name="Chlorine", symbol="Cl", z=17., a);
    a = 39.09*g/mole;
    G4Element* elK= new G4Element(name="Potassium", symbol="K", z=19., a);
    a = 40.08*g/mole;
    G4Element* elCa= new G4Element(name="Calcium", symbol="Ca", z=20., a);
    a = 55.85*g/mole;
    G4Element* elFe= new G4Element(name="Iron", symbol="Fe", z=26., a);
    a = 65.38*g/mole;
    G4Element* elZn= new G4Element(name="Zinc", symbol="Zn", z=26., a);
    a = 126.9*g/mole;
    G4Element* elI = new G4Element(name="Iodine", symbol="I", z=53., a);
    a = 132.9*g/mole;
    G4Element* elCs= new G4Element(name="Cesium", symbol="Cs", z=55., a);
    
    G4Element* elCr = man->FindOrBuildElement("Cr");
    G4Element* elMn = man->FindOrBuildElement("Mn");
    G4Element* elNi = man->FindOrBuildElement("Ni");
    G4Element* elTi = man->FindOrBuildElement("Ti");
    G4Element* elCu = man->FindOrBuildElement("Cu");
    G4Element* elSn = man->FindOrBuildElement("Sn");
    G4Element* elMo = man->FindOrBuildElement("Mo");
    G4Element* elV = man->FindOrBuildElement("V");
    
    //Concrete
    density = 2.3*g/cm3;
    G4Material* concrete = new G4Material(name="Concrete", density, ncomponents=10);
    concrete->AddElement( elH, 0.01    );
    concrete->AddElement( elC, 0.001   );
    concrete->AddElement( elO, 0.529107);
    concrete->AddElement( elNa, 0.016   );
    concrete->AddElement(elMg, 0.002   );
    concrete->AddElement(elAl, 0.033872);
    concrete->AddElement(elSi, 0.337021);
    concrete->AddElement(elK, 0.013   );
    concrete->AddElement(elCa, 0.044   );
    concrete->AddElement(elFe, 0.014   );
    materials.insert(std::pair<G4String, G4Material *>("concrete", concrete));
    
    //Concrete 2.4
    density = 2.4*g/cm3;
    G4Material* concrete_2_4 = new G4Material(name="Concrete_2.4", density, ncomponents=10);
    concrete_2_4->AddElement( elH, 0.01    );
    concrete_2_4->AddElement( elC, 0.001   );
    concrete_2_4->AddElement( elO, 0.529107);
    concrete_2_4->AddElement( elNa, 0.016   );
    concrete_2_4->AddElement(elMg, 0.002   );
    concrete_2_4->AddElement(elAl, 0.033872);
    concrete_2_4->AddElement(elSi, 0.337021);
    concrete_2_4->AddElement(elK, 0.013   );
    concrete_2_4->AddElement(elCa, 0.044   );
    concrete_2_4->AddElement(elFe, 0.014   );
    materials.insert(std::pair<G4String, G4Material *>("concrete_2.4", concrete_2_4));
    

    //Hematita Fe2 O3
    density = 5.27*g/cm3;
    G4Material* hematita = new G4Material(name="Hematita", density, ncomponents=2);
    hematita->AddElement(elFe, 2);
    hematita->AddElement(elO, 3);
    materials.insert(std::pair<G4String, G4Material *>("hematita", hematita));
    
    
    //Magnetita Fe3 O4
    density = 5.175*g/cm3;
    G4Material* magnetita = new G4Material(name="Magnetita", density, ncomponents=2);
    magnetita->AddElement(elFe, 3);
    magnetita->AddElement(elO, 4);
    materials.insert(std::pair<G4String, G4Material *>("magnetita", magnetita));
    
    //Wustita Fe O
    density = 5.75*g/cm3;
    G4Material* wustita = new G4Material(name="Wustita", density, ncomponents=2);
    wustita->AddElement(elFe, 1);
    wustita->AddElement(elO, 1);
    materials.insert(std::pair<G4String, G4Material *>("wustita", wustita));
    
    //Cascarilla 2% Hematita 4% Magnetita 94% Wustita
    density = 4.0*g/cm3;
    G4Material* cascarilla = new G4Material(name="Cascarilla", density, ncomponents=3);
    cascarilla->AddMaterial(hematita, fractionmass=0.02);
    cascarilla->AddMaterial(magnetita, fractionmass=0.04);
    cascarilla->AddMaterial(wustita, fractionmass=0.94);
    materials.insert(std::pair<G4String, G4Material *>("cascarilla", cascarilla));
    
    density = 3.50*g/cm3;
    G4Material* carbonRefractory = new G4Material(name="carbonRefractory",density,ncomponents = 1);
    carbonRefractory->AddElement(elC, fractionmass=100*perCent);
    materials.insert(std::pair<G4String, G4Material *>("carbonRefractory", carbonRefractory));

    density = 6.00*g/cm3;
    G4Material* liquidSteel = new G4Material("liquidSteel", density, ncomponents = 6);
    liquidSteel->AddElement(elC, fractionmass=0.001);
    liquidSteel->AddElement(elSi, fractionmass=0.007);
    liquidSteel->AddElement(elCr, fractionmass=0.18);
    liquidSteel->AddElement(elMn, fractionmass=0.01);
    liquidSteel->AddElement(elFe, fractionmass=0.712);
    liquidSteel->AddElement(elNi, fractionmass=0.09);
    materials.insert(std::pair<G4String, G4Material *>("liquidSteel", liquidSteel));
    
    density = 7.85*g/cm3;
    G4Material* steel = new G4Material("Steel", density, ncomponents = 14);
    steel->AddElement(elFe, fractionmass=0.98);
    steel->AddElement(elC,  fractionmass=0.001);
    steel->AddElement(elMn, fractionmass=0.007);
    steel->AddElement(elSi, fractionmass=0.002);
    
    steel->AddElement(elS,  fractionmass=0.001);
    steel->AddElement(elP,  fractionmass=0.001);
    steel->AddElement(elCu, fractionmass=0.001);
    steel->AddElement(elSn, fractionmass=0.001);
    steel->AddElement(elNi, fractionmass=0.001);
    steel->AddElement(elMo, fractionmass=0.001);
    steel->AddElement(elCr, fractionmass=0.001);
    steel->AddElement(elV,  fractionmass=0.001);
    steel->AddElement(elAl, fractionmass=0.001);
    steel->AddElement(elZn, fractionmass=0.001);
    materials.insert(std::pair<G4String, G4Material *>("steel", steel));
        
    density = 2.00*g/cm3;
    G4Material* slag = new G4Material("slag", density, ncomponents = 6);
    slag->AddElement(elC, fractionmass=0.001);
    slag->AddElement(elSi, fractionmass=0.007);
    slag->AddElement(elCr, fractionmass=0.18);
    slag->AddElement(elMn, fractionmass=0.01);
    slag->AddElement(elFe, fractionmass=0.712);
    slag->AddElement(elNi, fractionmass=0.09);
    materials.insert(std::pair<G4String, G4Material *>("slag", slag));

    density = 3.00*g/cm3;
    G4Material* scrap = new G4Material("scrap", density, ncomponents = 6);
    scrap->AddElement(elC, fractionmass=0.001);
    scrap->AddElement(elSi, fractionmass=0.007);
    scrap->AddElement(elCr, fractionmass=0.18);
    scrap->AddElement(elMn, fractionmass=0.01);
    scrap->AddElement(elFe, fractionmass=0.712);
    scrap->AddElement(elNi, fractionmass=0.09);
    materials.insert(std::pair<G4String, G4Material *>("scrap", scrap));

    density=70*kg/m3;
    G4Material* rockwool = new G4Material("rockwool", density, ncomponents = 8 );
    rockwool->AddElement(elO, fractionmass=0.481);
    rockwool->AddElement(elSi, fractionmass=0.277);
    rockwool->AddElement(elAl, fractionmass=0.081);
    rockwool->AddElement(elFe, fractionmass=0.05);
    rockwool->AddElement(elCa, fractionmass=0.036);
    rockwool->AddElement(elNa, fractionmass=0.028);
    rockwool->AddElement(elK, fractionmass=0.026);
    rockwool->AddElement(elMg, fractionmass=0.021);
    materials.insert(std::pair<G4String, G4Material *>("rockwool", rockwool));
    
    
    density = 1.000*g/cm3;
    G4int ncomp = 2;
    G4Material* water = new G4Material("Water",density,ncomp=2);
    G4int natoms;
    water->AddElement(elH, natoms=2);
    water->AddElement(elO, natoms=1);
    materials.insert(std::pair<G4String, G4Material *>("water", water));



    //Cuarzo
    density = 2.65*g/cm3;
    ncomp = 2;
    G4Material* SiO2 = new G4Material("SiO2",density,ncomp=2);
    SiO2->AddElement(elSi, natoms=1);
    SiO2->AddElement(elO, natoms=2);



    //Alumina
    density = 2.65*g/cm3;
    ncomp = 2;
    G4Material* Al2O3 = new G4Material("Al2O3",density,ncomp=2);
    Al2O3->AddElement(elAl, natoms=2);
    Al2O3->AddElement(elO, natoms=3);


    //Oxido de hierro III
    density = 5.242*g/cm3;
    ncomp = 2;
    G4Material* Fe2O3 = new G4Material("Fe2O3",density,ncomp=2);
    Fe2O3->AddElement(elFe, natoms=2);
    Fe2O3->AddElement(elO, natoms=3);


    //Oxido de manganeso
    density = 5.43*g/cm3;
    ncomp = 2;
    G4Material* MnO = new G4Material("MnO",density,ncomp=2);
    MnO->AddElement(elMn, natoms=1);
    MnO->AddElement(elO, natoms=1);


    //Oxido de magnesio
    density =3.6*g/cm3;
    ncomp = 2;
    G4Material* MgO = new G4Material("MgO",density,ncomp=2);
    MgO->AddElement(elMg, natoms=1);
    MgO->AddElement(elO, natoms=1);


    //Oxido de Calcioden
    density =3.3*g/cm3;
    ncomp = 2;
    G4Material* CaO = new G4Material("CaO",density,ncomp=2);
    CaO->AddElement(elCa, natoms=1);
    CaO->AddElement(elO, natoms=1);


    //Oxido de Sodio
    density =2.27*g/cm3;
    ncomp = 2;
    G4Material* Na2O = new G4Material("Na2O",density,ncomp=2);
    Na2O->AddElement(elNa, natoms=2);
    Na2O->AddElement(elO, natoms=1);


    //Oxido de Potasio
    density =2.13*g/cm3;
    ncomp = 2;
    G4Material* K2O = new G4Material("K2O",density,ncomp=2);
    K2O->AddElement(elK, natoms=2);
    K2O->AddElement(elO, natoms=1);


    //Oxido de Titanio
    density =4.2*g/cm3;
    ncomp = 2;
    G4Material* TiO2 = new G4Material("TiO2",density,ncomp=2);
    TiO2->AddElement(elTi, natoms=1);
    TiO2->AddElement(elO, natoms=2);


    //Oxido de Fosforo
    density =2.39*g/cm3;
    ncomp = 2;
    G4Material* P2O5 = new G4Material("P2O5",density,ncomp=2);
    P2O5->AddElement(elP, natoms=2);
    P2O5->AddElement(elO, natoms=5);

    //L.O.I. aproximado con agua
    density = 1.000*g/cm3;
    ncomp = 2;
    G4Material* LOI = new G4Material("LOI",density,ncomp=2);
    LOI->AddElement(elH, natoms=2);
    LOI->AddElement(elO, natoms=1);

    //Tailing 2.1
    density = 2.1*g/cm3;
    G4Material* tailing_21 = new G4Material("tailing_2.1", density, ncomponents = 11);
    tailing_21->AddMaterial(SiO2, fractionmass=0.712);
    tailing_21->AddMaterial(Al2O3,  fractionmass=0.0395);
    tailing_21->AddMaterial(Fe2O3, fractionmass=0.2032);
    tailing_21->AddMaterial(MnO, fractionmass=0.0008);
    tailing_21->AddMaterial(MgO, fractionmass=0.0072);
    tailing_21->AddMaterial(CaO, fractionmass=0.0048);
    tailing_21->AddMaterial(Na2O, fractionmass=0.0021);
    tailing_21->AddMaterial(K2O, fractionmass=0.0075);
    tailing_21->AddMaterial(TiO2, fractionmass=0.0021);
    tailing_21->AddMaterial(P2O5, fractionmass=0.0015);
    tailing_21->AddMaterial(LOI, fractionmass=0.0193);
    
    materials.insert(std::pair<G4String, G4Material *>("tailing_2.1", tailing_21));


    //Tailing 2.0
    density = 2.0*g/cm3;
    G4Material* tailing_20 = new G4Material("tailing_2.0", density, ncomponents = 11);
    tailing_20->AddMaterial(SiO2, fractionmass=0.712);
    tailing_20->AddMaterial(Al2O3,  fractionmass=0.0395);
    tailing_20->AddMaterial(Fe2O3, fractionmass=0.2032);
    tailing_20->AddMaterial(MnO, fractionmass=0.0008);
    tailing_20->AddMaterial(MgO, fractionmass=0.0072);
    tailing_20->AddMaterial(CaO, fractionmass=0.0048);
    tailing_20->AddMaterial(Na2O, fractionmass=0.0021);
    tailing_20->AddMaterial(K2O, fractionmass=0.0075);
    tailing_20->AddMaterial(TiO2, fractionmass=0.0021);
    tailing_20->AddMaterial(P2O5, fractionmass=0.0015);
    tailing_20->AddMaterial(LOI, fractionmass=0.0193);
    
    materials.insert(std::pair<G4String, G4Material *>("tailing_2.0", tailing_20));


    //Tailing 1.9
    density = 1.9*g/cm3;
    G4Material* tailing_19 = new G4Material("tailing_1.9", density, ncomponents = 11);
    tailing_19->AddMaterial(SiO2, fractionmass=0.712);
    tailing_19->AddMaterial(Al2O3,  fractionmass=0.0395);
    tailing_19->AddMaterial(Fe2O3, fractionmass=0.2032);
    tailing_19->AddMaterial(MnO, fractionmass=0.0008);
    tailing_19->AddMaterial(MgO, fractionmass=0.0072);
    tailing_19->AddMaterial(CaO, fractionmass=0.0048);
    tailing_19->AddMaterial(Na2O, fractionmass=0.0021);
    tailing_19->AddMaterial(K2O, fractionmass=0.0075);
    tailing_19->AddMaterial(TiO2, fractionmass=0.0021);
    tailing_19->AddMaterial(P2O5, fractionmass=0.0015);
    tailing_19->AddMaterial(LOI, fractionmass=0.0193);
    
    materials.insert(std::pair<G4String, G4Material *>("tailing_1.9", tailing_19));


    //Tailing 1.8
    density = 1.8*g/cm3;
    G4Material* tailing_18 = new G4Material("tailing_1.8", density, ncomponents = 11);
    tailing_18->AddMaterial(SiO2, fractionmass=0.712);
    tailing_18->AddMaterial(Al2O3,  fractionmass=0.0395);
    tailing_18->AddMaterial(Fe2O3, fractionmass=0.2032);
    tailing_18->AddMaterial(MnO, fractionmass=0.0008);
    tailing_18->AddMaterial(MgO, fractionmass=0.0072);
    tailing_18->AddMaterial(CaO, fractionmass=0.0048);
    tailing_18->AddMaterial(Na2O, fractionmass=0.0021);
    tailing_18->AddMaterial(K2O, fractionmass=0.0075);
    tailing_18->AddMaterial(TiO2, fractionmass=0.0021);
    tailing_18->AddMaterial(P2O5, fractionmass=0.0015);
    tailing_18->AddMaterial(LOI, fractionmass=0.0193);
    
    materials.insert(std::pair<G4String, G4Material *>("tailing_1.8", tailing_18));


    //Tailing 1.7
    density = 1.7*g/cm3;
    G4Material* tailing_17 = new G4Material("tailing_1.7", density, ncomponents = 11);
    tailing_17->AddMaterial(SiO2, fractionmass=0.712);
    tailing_17->AddMaterial(Al2O3,  fractionmass=0.0395);
    tailing_17->AddMaterial(Fe2O3, fractionmass=0.2032);
    tailing_17->AddMaterial(MnO, fractionmass=0.0008);
    tailing_17->AddMaterial(MgO, fractionmass=0.0072);
    tailing_17->AddMaterial(CaO, fractionmass=0.0048);
    tailing_17->AddMaterial(Na2O, fractionmass=0.0021);
    tailing_17->AddMaterial(K2O, fractionmass=0.0075);
    tailing_17->AddMaterial(TiO2, fractionmass=0.0021);
    tailing_17->AddMaterial(P2O5, fractionmass=0.0015);
    tailing_17->AddMaterial(LOI, fractionmass=0.0193);
    
    materials.insert(std::pair<G4String, G4Material *>("tailing_1.7", tailing_17));


    //Tailing 1.6
    density = 1.6*g/cm3;
    G4Material* tailing_16 = new G4Material("tailing_1.6", density, ncomponents = 11);
    tailing_16->AddMaterial(SiO2, fractionmass=0.712);
    tailing_16->AddMaterial(Al2O3,  fractionmass=0.0395);
    tailing_16->AddMaterial(Fe2O3, fractionmass=0.2032);
    tailing_16->AddMaterial(MnO, fractionmass=0.0008);
    tailing_16->AddMaterial(MgO, fractionmass=0.0072);
    tailing_16->AddMaterial(CaO, fractionmass=0.0048);
    tailing_16->AddMaterial(Na2O, fractionmass=0.0021);
    tailing_16->AddMaterial(K2O, fractionmass=0.0075);
    tailing_16->AddMaterial(TiO2, fractionmass=0.0021);
    tailing_16->AddMaterial(P2O5, fractionmass=0.0015);
    tailing_16->AddMaterial(LOI, fractionmass=0.0193);
    
    materials.insert(std::pair<G4String, G4Material *>("tailing_1.6", tailing_16));


    //Tailing 1.5
    density = 1.5*g/cm3;
    G4Material* tailing_15 = new G4Material("tailing_1.5", density, ncomponents = 11);
    tailing_15->AddMaterial(SiO2, fractionmass=0.712);
    tailing_15->AddMaterial(Al2O3,  fractionmass=0.0395);
    tailing_15->AddMaterial(Fe2O3, fractionmass=0.2032);
    tailing_15->AddMaterial(MnO, fractionmass=0.0008);
    tailing_15->AddMaterial(MgO, fractionmass=0.0072);
    tailing_15->AddMaterial(CaO, fractionmass=0.0048);
    tailing_15->AddMaterial(Na2O, fractionmass=0.0021);
    tailing_15->AddMaterial(K2O, fractionmass=0.0075);
    tailing_15->AddMaterial(TiO2, fractionmass=0.0021);
    tailing_15->AddMaterial(P2O5, fractionmass=0.0015);
    tailing_15->AddMaterial(LOI, fractionmass=0.0193);
    
    materials.insert(std::pair<G4String, G4Material *>("tailing_1.5", tailing_15));

    
    //Fiber Glass
    density = 2.5*g/cm3;
    G4Material* glassFiber = new G4Material("glassFiber", density, ncomponents = 2);
    glassFiber->AddElement(elSi, natoms=1);
    glassFiber->AddElement(elO, natoms=4);
    
    
    materials.insert(std::pair<G4String, G4Material *>("glassfiber", glassFiber));


    
    //Electrode
    density = 1.65*g/cm3;
    ncomp = 1;
    G4Material* electrode = new G4Material("electrode",density,ncomp=1);
    electrode->AddElement(elC, natoms=1);    
    materials.insert(std::pair<G4String, G4Material *>("electrode", electrode));
    
    density = 2.60*g/cm3;
    G4Material* carbonRefractory_electricFurnace = new G4Material(name="carbonRefractory_electricFurnace",density,ncomponents = 1);
    carbonRefractory_electricFurnace->AddElement(elC, fractionmass=100*perCent);
    materials.insert(std::pair<G4String, G4Material *>("carbonRefractory_electricFurnace", carbonRefractory_electricFurnace));
    
    density = 1.60*g/cm3;
    G4Material* carbon_electricFurnace = new G4Material(name="carbon_electricFurnace",density,ncomponents = 1);
    carbon_electricFurnace->AddElement(elC, fractionmass=100*perCent);
    materials.insert(std::pair<G4String, G4Material *>("carbon_electricFurnace", carbon_electricFurnace));
    
    density = 2.00*g/cm3;
    G4Material* silicon_electricFurnace = new G4Material(name="silicon_electricFurnace",density,ncomponents = 1);
    silicon_electricFurnace->AddElement(elSi, fractionmass=100*perCent);
    materials.insert(std::pair<G4String, G4Material *>("silicon_electricFurnace", silicon_electricFurnace));
    
    density = 2.50*g/cm3;
    G4Material* liquidSilicon = new G4Material(name="liquidSilicon",density,ncomponents = 1);
    liquidSilicon->AddElement(elSi, fractionmass=100*perCent);
    materials.insert(std::pair<G4String, G4Material *>("liquidSilicon", liquidSilicon));
    
    
//    density =2.39*g/cm3;
//    ncomp = 2;
//    G4Material* carbon = new G4Material("P2O5",density,ncomp=2);
//    P2O5->AddElement(elP, natoms=2);
//    P2O5->AddElement(elO, natoms=5);
    
    density = 2.00*g/cm3;
    G4Material* quartzCarbon = new G4Material(name="quartzCarbon",density,ncomponents = 2);
    quartzCarbon->AddMaterial(SiO2, fractionmass=0.75);
    quartzCarbon->AddElement(elC,  fractionmass=0.25);
    materials.insert(std::pair<G4String, G4Material *>("quartzCarbon", quartzCarbon));
    
    density = 0.4*g/cm3;
    G4Material* snow_4 = new G4Material(name="snow_400", density, ncomponents = 2);
    snow_4->AddElement(elH, 2);
    snow_4->AddElement(elO,1);
    materials.insert(std::pair<G4String, G4Material *>("snow_400", snow_4));

    density = 0.38*g/cm3;
    G4Material* snow_38 = new G4Material(name="snow_380", density, ncomponents = 2);
    snow_38->AddElement(elH, 2);
    snow_38->AddElement(elO,1);
    materials.insert(std::pair<G4String, G4Material *>("snow_380", snow_38));
    
    density = 0.36*g/cm3;
    G4Material* snow_36 = new G4Material(name="snow_360", density, ncomponents = 2);
    snow_36->AddElement(elH, 2);
    snow_36->AddElement(elO,1);
    materials.insert(std::pair<G4String, G4Material *>("snow_360", snow_36));

    density = 0.34*g/cm3;
    G4Material* snow_34 = new G4Material(name="snow_340", density, ncomponents = 2);
    snow_34->AddElement(elH, 2);
    snow_34->AddElement(elO,1);
    materials.insert(std::pair<G4String, G4Material *>("snow_340", snow_34));

    density = 0.32*g/cm3;
    G4Material* snow_32 = new G4Material(name="snow_320", density, ncomponents = 2);
    snow_32->AddElement(elH, 2);
    snow_32->AddElement(elO,1);
    materials.insert(std::pair<G4String, G4Material *>("snow_320", snow_32));

    density = 0.30*g/cm3;
    G4Material* snow_30 = new G4Material(name="snow_300", density, ncomponents = 2);
    snow_30->AddElement(elH, 2);
    snow_30->AddElement(elO,1);
    materials.insert(std::pair<G4String, G4Material *>("snow_300", snow_30));

//    // Acero de diferentes densidades
//    density = 7.85*g/cm3;
//    G4Material* steel = new G4Material("Steel", density, ncomponents = 14);
//    steel->AddElement(elFe, fractionmass=0.98);
//    steel->AddElement(elC,  fractionmass=0.001);
//    steel->AddElement(elMn, fractionmass=0.007);
//    steel->AddElement(elSi, fractionmass=0.002);
//    
//    steel->AddElement(elS,  fractionmass=0.001);
//    steel->AddElement(elP,  fractionmass=0.001);
//    steel->AddElement(elCu, fractionmass=0.001);
//    steel->AddElement(elSn, fractionmass=0.001);
//    steel->AddElement(elNi, fractionmass=0.001);
//    steel->AddElement(elMo, fractionmass=0.001);
//    steel->AddElement(elCr, fractionmass=0.001);
//    steel->AddElement(elV,  fractionmass=0.001);
//    steel->AddElement(elAl, fractionmass=0.001);
//    steel->AddElement(elZn, fractionmass=0.001);
//  
//    materials.insert(std::pair<G4String, G4Material *>("steel", steel));
//    
    
    
    
    
    
    
    
    

  ///////////////////////////////////////////
  ///          Biological materials       ///
  ///////////////////////////////////////////
  //Bone Compact
  density = 1.85* g/cm3;
  G4int nElem = 0;
  G4Material* boneCompact = new G4Material( "boneCompact", density, nElem = 8 );
  boneCompact -> AddElement ( elH, 0.063984 ); //Z=1
  boneCompact -> AddElement ( elC, 0.278 ); //Z=6
  boneCompact -> AddElement ( elN, 0.027 ); //Z=7
  boneCompact -> AddElement ( elO, 0.410016 ); //Z=8
  boneCompact -> AddElement ( elMg, 0.002 ); //Z=12
  boneCompact -> AddElement ( elP, 0.07 ); //Z=15
  boneCompact -> AddElement ( elS, 0.002 ); //Z=16
  boneCompact -> AddElement ( elCa, 0.147 ); //Z=20
  
  materials.insert(std::pair<G4String, G4Material *>("boneCompact", boneCompact));

  //Bone Cortical
  density = 1.85* g/cm3;
  G4Material* boneCortical = new G4Material( "boneCortical", density, nElem = 9 );
  boneCortical -> AddElement ( elH, 0.047234 ); //Z=1
  boneCortical -> AddElement ( elC, 0.144330 ); //Z=6
  boneCortical -> AddElement ( elN, 0.041990 ); //Z=7
  boneCortical -> AddElement ( elO, 0.446096 ); //Z=8
  boneCortical -> AddElement ( elMg, 0.002200 ); //Z=12
  boneCortical -> AddElement ( elP, 0.104970 ); //Z=15
  boneCortical -> AddElement ( elS, 0.003150 ); //Z=16
  boneCortical -> AddElement ( elCa, 0.209930 ); //Z=20
  boneCortical -> AddElement ( elZn, 0.000100 ); //Z=30

  materials.insert(std::pair<G4String, G4Material *>("boneCortical", boneCortical));

  //Muscle
  density = 1.04* g/cm3;
  G4Material* muscle = new G4Material( "muscleSkeletal", density, nElem = 13 );
  muscle -> AddElement ( elH, 0.100637 ); //Z=1
  muscle -> AddElement ( elC, 0.107830 ); //Z=6
  muscle -> AddElement ( elN, 0.027680 ); //Z=7
  muscle -> AddElement ( elO, 0.754773 ); //Z=8
  muscle -> AddElement ( elNa, 0.000750 ); //Z=11
  muscle -> AddElement ( elMg, 0.000190 ); //Z=12
  muscle -> AddElement ( elP, 0.001800 ); //Z=15
  muscle -> AddElement ( elS, 0.002410 ); //Z=16
  muscle -> AddElement ( elCl, 0.000790 ); //Z=17
  muscle -> AddElement ( elK, 0.003020 ); //Z=19
  muscle -> AddElement ( elCa, 0.000030 ); //Z=20
  muscle -> AddElement ( elFe, 0.000040 ); //Z=26
  muscle -> AddElement ( elZn, 0.000050 ); //Z=30

  materials.insert(std::pair<G4String, G4Material *>("muscle", muscle));
  
  //Blood
  density = 1.06* g/cm3;
  G4Material* blood = new G4Material( "blood", density, nElem = 14 );
  blood -> AddElement (elH, 0.101866); //Z=1
  blood -> AddElement (elC, 0.100020); //Z=6
  blood -> AddElement (elN, 0.029640); //Z=7
  blood -> AddElement (elO, 0.759414); //Z=8
  blood -> AddElement (elNa, 0.001850); //Z=11
  blood -> AddElement (elMg, 0.000040); //Z=12
  blood -> AddElement (elSi, 0.000030); //Z=14
  blood -> AddElement (elP, 0.000350); //Z=15
  blood -> AddElement (elS, 0.001850); //Z=16
  blood -> AddElement (elCl, 0.002780); //Z=17
  blood -> AddElement (elK, 0.001630); //Z=19
  blood -> AddElement (elCa, 0.000060); //Z=20
  blood -> AddElement (elFe, 0.000460); //Z=26
  blood -> AddElement (elZn, 0.000010); //Z=30
  
  materials.insert(std::pair<G4String, G4Material *>("muscle", muscle));

  //Skin
  density = 1.1* g/cm3;
  G4Material* skin = new G4Material( "skin", density, nElem = 13 );
  skin -> AddElement (elH, 0.100588); //Z=1
  skin -> AddElement (elC, 0.228250); //Z=6
  skin -> AddElement (elN, 0.046420); //Z=7
  skin -> AddElement (elO, 0.619002); //Z=8
  skin -> AddElement (elNa, 0.000070); //Z=11
  skin -> AddElement (elMg, 0.000060); //Z=12
  skin -> AddElement (elP, 0.000330); //Z=15
  skin -> AddElement (elS, 0.001590); //Z=16
  skin -> AddElement (elCl, 0.002670); //Z=17
  skin -> AddElement (elK, 0.000850); //Z=19
  skin -> AddElement (elCa, 0.000150); //Z=20
  skin -> AddElement (elFe, 0.000010); //Z=26
  skin -> AddElement (elZn, 0.000010); //Z=30

  materials.insert(std::pair<G4String, G4Material *>("skin", skin));
  
  //AdiposeTissue
  density = 0.92* g/cm3;
  G4Material* adiposeTissue = new G4Material( "adiposeTissue", density, nElem = 13 );
  adiposeTissue -> AddElement (elH, 0.119477); //Z=1
  adiposeTissue -> AddElement (elC, 0.637240); //Z=6
  adiposeTissue -> AddElement (elN, 0.007970); //Z=7
  adiposeTissue -> AddElement (elO, 0.232333); //Z=8
  adiposeTissue -> AddElement (elNa, 0.000500); //Z=11
  adiposeTissue -> AddElement (elMg, 0.000020); //Z=12
  adiposeTissue -> AddElement (elP, 0.000160); //Z=15
  adiposeTissue -> AddElement (elS, 0.000730); //Z=16
  adiposeTissue -> AddElement (elCl, 0.001190); //Z=17
  adiposeTissue -> AddElement (elK, 0.000320); //Z=19
  adiposeTissue -> AddElement (elCa, 0.000020); //Z=20
  adiposeTissue -> AddElement (elFe, 0.000020); //Z=26
  adiposeTissue -> AddElement (elZn, 0.000020); //Z=30
  
  materials.insert(std::pair<G4String, G4Material *>("adiposeTissue", adiposeTissue));






}
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Destroy all the materials                                            //
//----------------------------------------------------------------------//
void DetectorConstruction::DestroyMaterials() {
    // Destroy all allocated elements and materials
    size_t i;
    G4MaterialTable* matTable = (G4MaterialTable*)G4Material::GetMaterialTable();
    for(i=0; i<matTable->size(); i++) delete (*(matTable))[i];
    matTable->clear();
    G4ElementTable* elemTable = (G4ElementTable*)G4Element::GetElementTable();
    for(i=0; i<elemTable->size(); i++) delete (*(elemTable))[i];
    elemTable->clear();

}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//



void DetectorConstruction::DumpGeometricalTree(G4VPhysicalVolume* aVolume,G4int depth)
{

    for(int isp=0; isp<depth; isp++)
    {
        G4cout << "  ";
    }
    G4cout << aVolume->GetName() << "[" << aVolume->GetCopyNo() << "] "
           << aVolume->GetLogicalVolume()->GetName() << " "
           << aVolume->GetLogicalVolume()->GetNoDaughters() << " "
           << aVolume->GetLogicalVolume()->GetMaterial()->GetName();
    if(aVolume->GetLogicalVolume()->GetSensitiveDetector())
    {
        G4cout << " " << aVolume->GetLogicalVolume()->GetSensitiveDetector()->GetFullPathName();
    }
    G4cout << G4endl;
    for(int i=0; i<aVolume->GetLogicalVolume()->GetNoDaughters(); i++)
    {
        DumpGeometricalTree(aVolume->GetLogicalVolume()->GetDaughter(i),depth+1);
    }

}
