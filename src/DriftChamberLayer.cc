#include "DriftChamberLayer.hh"
#include "DriftChamberLayerHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4Navigator.hh"
#include "G4ios.hh"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandFlat.h"
#include "math.h"

//----------------------------------------------------------------------//
// Constructor                                                          //
//----------------------------------------------------------------------//
DriftChamberLayer::DriftChamberLayer(G4String name, G4String collection)
    :G4VSensitiveDetector(name)
{
    G4String HCname;
    collectionName.insert(collection);
    HCID = -1;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Constructor                                                          //
//----------------------------------------------------------------------//
DriftChamberLayer::~DriftChamberLayer() {
    ;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Constructor                                                          //
//----------------------------------------------------------------------//
void DriftChamberLayer::Initialize(G4HCofThisEvent*HCE)
{
    hitsCollection = new DriftChamberLayerHitsCollection(SensitiveDetectorName, collectionName[0]);
    if(HCID<0) HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection);
    HCE->AddHitsCollection(HCID,hitsCollection);
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Process the hits and fills the relevant information                  //
//----------------------------------------------------------------------//

G4bool DriftChamberLayer::ProcessHits(G4Step*aStep,G4TouchableHistory*  /*ROhist*/) {
    

    
    G4double charge = aStep->GetTrack()->GetDefinition()->GetPDGCharge();
    if(charge==0.) return true;
    G4int type = abs(aStep->GetTrack()->GetDefinition()->GetPDGEncoding());
    if(abs(type) != 13) return true;

    G4String preMaterial = aStep->GetPreStepPoint()->GetMaterial()->GetName();
    G4String postMaterial = aStep->GetPostStepPoint()->GetMaterial()->GetName();
    if(preMaterial != "G4_Ar" || postMaterial != "G4_Cu") return true;
  

    G4StepPoint* preStepPoint = aStep->GetPreStepPoint();    
    G4TouchableHistory* theTouchable = (G4TouchableHistory*)(preStepPoint->GetTouchable());
    G4VPhysicalVolume* theMotherPhysical = theTouchable->GetVolume(1); // mother
    G4VPhysicalVolume* thePhysical = theTouchable->GetVolume(0); // child
    // localizacion de preStep en el sistema de ref. del mundo
    G4ThreeVector prePos_world = preStepPoint->GetPosition();
    // localizacion de preStep en el sistema de ref. del detector
    G4ThreeVector prePos_detector = det->toLocal(prePos_world);
    
    // localizacion de postStep en el sistema de ref. del mundo
    G4ThreeVector postPos_world = aStep->GetPostStepPoint()->GetPosition();
    // localizacion de postStep en el sistema de ref. del detector
    G4ThreeVector postPos_detector = det->toLocal(postPos_world);
    
    
    G4int copyNo = -1;
    
    

    for(G4int a = 0; a < det->getNLayers(); a++) {
        G4ThreeVector poslayer = det->getPosLayer(a);
        G4ThreeVector sizelayer = det->getSizeLayer(a);
        if(prePos_detector[2] >= poslayer[2] - 1.*sizelayer[2]/2.0 && prePos_detector[2] <= poslayer[2] + 1.*sizelayer[2]/2.0) {
            copyNo = a;
            break;
        }
    }
    
    if(copyNo == -1) return true;//TODO: esto no debería ser false?????

    //Simulating efficiency
    if(CLHEP::RandFlat::shoot(0.0, 1.0) > det->getEffLayer(copyNo)) return false;
    
    
    // localizacion de preStep en el sistema de ref. del layer
    G4ThreeVector prePos_layer = det->toLocal(det->toLocal(prePos_world), copyNo); 
    G4double energy = preStepPoint->GetTotalEnergy();
    
    // localizacion de postStep en el sistema de ref. del layer
    G4ThreeVector postPos_layer = det->toLocal(det->toLocal(postPos_world), copyNo);

    // Obtengo la dirección de la particula en el sistema de referencia del detector
    G4ThreeVector particleDir_detector = postPos_detector - prePos_detector;
    // Obtengo la dirección de la particula en el sistema de referencia del layer
    G4ThreeVector particleDir_layer = postPos_layer - prePos_layer;
    
    
    // Defino la Z donde están los hilos en el sistema de referencia del layer
    //TODO: el valor esta hardcoded pero habría que leerlo del detector de alguna manera.
    G4double wires_z_layer = 1./2. * CLHEP::cm;
    if(copyNo%2 == 0) wires_z_layer = -1./2. * CLHEP::cm;
    
    // Genero un punto centrado en el origen X=0 Y=0 con la z de los hilos en el layer
     G4ThreeVector wirePos_layer(0,0,wires_z_layer);
    
    // Defino la Z de los hilos en coordenadas del detector
    G4ThreeVector wirePos_detector = det->toGlobal(wirePos_layer, copyNo);
//    std::cout << "wirePos Global: ( "<<wirePos_detector<<" )"<<std::endl;    
    G4double wires_z_detector = wirePos_detector.z();

   
    
    //Teniendo en cuenta la dirección calculo el corte con la posicion de los hilos
    G4double w_x_detector = prePos_detector.x() + (wires_z_detector - prePos_detector.z()) *  particleDir_detector.x()/particleDir_detector.z();
    G4double w_y_detector = prePos_detector.y() + (wires_z_detector - prePos_detector.z()) *  particleDir_detector.y()/particleDir_detector.z();

    // Genero el punto de corte con el eje de los hilos en el sis. ref. detector
     G4ThreeVector w_detector(w_x_detector, w_y_detector, wires_z_detector);
    G4double theta_x_detector = M_PI/2. - atan(fabs(particleDir_detector.z()/particleDir_detector.x()));
    G4double theta_y_detector = M_PI/2. - atan(fabs(particleDir_detector.z()/particleDir_detector.y()));
    
      
    
    // localizacion de preStep en el sistema de ref. del layer
    G4ThreeVector w_layer = det->toLocal(w_detector, copyNo);
            
    // Si el corte con los hilos queda fuera del área sensor cancelamos el hit. 
    // La traza debe pasar muy cerca de un hilo para que haya señal
    if(fabs(w_layer.x()) > (det->getSizeLayer(copyNo)/2.).x()) return true;
    
//    std::vector<G4double> pair = getWirePos(prePos_layer.x(), copyNo);
    std::vector<G4double> pair = getProbabilisticWirePos(w_layer.x(), particleDir_layer, copyNo);
       
    
    if(pair[1] < 0) return false;
//    G4ThreeVector LayerMeas(pair[0], 0, 0);
    G4ThreeVector LayerMeas(pair[0], 0, w_layer.z());
    G4ThreeVector ChamberMeas = det->toGlobal(LayerMeas, copyNo);
       
    //Simulating resolution
    DriftChamberLayerHit* aHit = new DriftChamberLayerHit(copyNo);
    aHit->SetLocalPos(prePos_detector);
    aHit->SetLocalMeas(ChamberMeas);
    aHit->SetTime(preStepPoint->GetGlobalTime());
    aHit->SetEnergy(energy);
    aHit->SetLayerID(copyNo);
    

    hitsCollection->insert(aHit);
    return true;
}



G4bool DriftChamberLayer::ProcessHits_old(G4Step*aStep,G4TouchableHistory*  /*ROhist*/) {


    
    G4double charge = aStep->GetTrack()->GetDefinition()->GetPDGCharge();
    if(charge==0.) return true;
    G4int type = abs(aStep->GetTrack()->GetDefinition()->GetPDGEncoding());
    if(abs(type) != 13) return true;


    G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
    G4TouchableHistory* theTouchable = (G4TouchableHistory*)(preStepPoint->GetTouchable());
    G4VPhysicalVolume* theMotherPhysical = theTouchable->GetVolume(1); // mother
    G4VPhysicalVolume* thePhysical = theTouchable->GetVolume(0); // child
    G4ThreeVector worldPos = preStepPoint->GetPosition();
    G4ThreeVector ChamberPos = det->toLocal(worldPos);
    G4int copyNo = -1;
    
    

    for(G4int a = 0; a < det->getNLayers(); a++) {
        G4ThreeVector poslayer = det->getPosLayer(a);
        G4ThreeVector sizelayer = det->getSizeLayer(a);
        if(ChamberPos[2] >= poslayer[2] - 1.*sizelayer[2]/2.0 && ChamberPos[2] <= poslayer[2] + 1.*sizelayer[2]/2.0) {
            copyNo = a;
            break;
        }
    }
    
    if(copyNo == -1) return true;//TODO: esto no debería ser false?????

    //Simulating efficiency
    if(CLHEP::RandFlat::shoot(0.0, 1.0) > det->getEffLayer(copyNo)) return false;
    
//    std::cout << "Layer Hit: " << copyNo << std::endl;
//    std::cout << "Prestep material: " << aStep->GetPreStepPoint()->GetMaterial()->GetName() << std::endl;
//    std::cout << "Poststep material: " << aStep->GetPostStepPoint()->GetMaterial()->GetName() << std::endl;
//    std::cout << "Delta step: " << aStep ->GetDeltaPosition() << std::endl;
//    std::cout << "Posicion Pre: (" << aStep->GetPreStepPoint()->GetPosition().x()<<" , "<< aStep->GetPreStepPoint()->GetPosition().y() << " , " << aStep->GetPreStepPoint()->GetPosition().z()<<" )"<<std::endl;
//    std::cout << "Posicion Post: (" << aStep->GetPostStepPoint()->GetPosition().x()<<" , "<< aStep->GetPostStepPoint()->GetPosition().y() << " , " << aStep->GetPostStepPoint()->GetPosition().z()<<" )"<<std::endl;
//    if(det->getPos().z() == -500 && copyNo==3) std::cout << ""<<std::endl;
    
    
    G4ThreeVector LayerPos = det->toLocal(ChamberPos, copyNo);
    G4double energy = preStepPoint->GetTotalEnergy();

    std::vector<G4double> pair = getWirePos(LayerPos.x(), copyNo);
    if(pair[1] < 0) return false;
    G4ThreeVector LayerMeas(pair[0], 0, 0);
    G4ThreeVector ChamberMeas = det->toGlobal(LayerMeas, copyNo);

    //Simulating resolution
    DriftChamberLayerHit* aHit = new DriftChamberLayerHit(copyNo);
    aHit->SetLocalPos(ChamberPos);
    aHit->SetLocalMeas(ChamberMeas);
    aHit->SetTime(preStepPoint->GetGlobalTime());
    aHit->SetEnergy(energy);
    aHit->SetLayerID(copyNo);
    

    hitsCollection->insert(aHit);

    return true;
}

//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Get the position of the closest wire                                 //
//----------------------------------------------------------------------//

std::vector<G4double> DriftChamberLayer::getProbabilisticWirePos(G4double position, G4ThreeVector direction, G4int a) {
    G4double w = position;
    
    //Obtengo el hilo 
    std::vector<G4double> val;
    G4int closestWire = -1;
    G4int closestWire_2 = -1;
    G4double minimumSoFar = 20*CLHEP::m;
    G4double minimumSoFar_2 = 20*CLHEP::m;
    std::vector<G4double> wires = (det->getWires())[a];
    if(wires.size() < 2) return val;
    
    G4double distance = wires[1]-wires[0];
    G4int approxi = (int)((w - wires[0])/distance);
    G4int imin = (int)(fmax(0, approxi-2));
    G4int imax = (int)(fmin(wires.size()-1, approxi+2));
    // Obtengo el hilo más cercano al punto de corte
    for(G4int i = imin; i <= imax; i++) {
        G4double dist = fabs(wires[i] - w);
        if(dist < minimumSoFar) {
            minimumSoFar = dist;
            closestWire = i;
        }
    }
        
    //Calculo el angulo de la traza en la cámara
    G4double theta_x_layer = M_PI/2. - atan(fabs(direction.z()/direction.x()));
    
    //Simulo el encendido de hilos
    //Obtengo un valor aleatorio
    G4double myRand = CLHEP::RandFlat::shoot(0.0, 1.0);
    G4int wiresOn = 1;
    G4double P1 = 0, P2 = 0, P3 = 0, P4 = 0, P5 = 0;
    
    //TODO: Añadir los porcentajes de cada ángulo a traves de un fichero.
    // introducimos los porcentajes asociados a cada cámara a través de datos reales.
    // los datos reales
    if(theta_x_layer < 5. * M_PI / 180.){
        
        if (a==0){
            P5 = 0.0215;
            P4 = 0.0216;
            P3 = 0.0816;
            P2 = 0.1543;        
        }else if (a==1){
            P5 = 0.0087;
            P4 = 0.0122;
            P3 = 0.0328;
            P2 = 0.1246;  
        }else if (a==2){
            P5 = 0.0031;
            P4 = 0.0039;
            P3 = 0.0098;
            P2 = 0.0439;
        }else if (a==3){
            P5 = 0.0174;
            P4 = 0.0254;
            P3 = 0.0614;
            P2 = 0.1455;
        }
        
    }else if(theta_x_layer < 10. * M_PI / 180.){
        if (a==0){
            P5 = 0.0206;
            P4 = 0.0283;
            P3 = 0.0882;
            P2 = 0.2473;        
        }else if (a==1){
            P5 = 0.0083;
            P4 = 0.0119;
            P3 = 0.0337;
            P2 = 0.1824;
        }else if (a==2){
            P5 = 0.0029;
            P4 = 0.0046;
            P3 = 0.0106;
            P2 = 0.0867;
        }else if (a==3){
            P5 = 0.0203;
            P4 = 0.0253;
            P3 = 0.069;
            P2 = 0.1874;
        }   
    }else if(theta_x_layer < 15. * M_PI / 180.){
        if (a==0){
            P5 = 0.0214;
            P4 = 0.0353;
            P3 = 0.0977;
            P2 = 0.3341;        
        }else if (a==1){
            P5 = 0.0099;
            P4 = 0.0144;
            P3 = 0.0371;
            P2 = 0.2694;
        }else if (a==2){
            P5 = 0.0029;
            P4 = 0.0048;
            P3 = 0.011;
            P2 = 0.1311;
        }else if (a==3){
            P5 = 0.0211;
            P4 = 0.0305;
            P3 = 0.0799;
            P2 = 0.2632;
        }       
    }else if(theta_x_layer < 20. * M_PI / 180.){
        if (a==0){
            P5 = 0.235;
            P4 = 0.416;
            P3 = 0.1124;
            P2 = 0.4184;        
        }else if (a==1){
            P5 = 0.0101;
            P4 = 0.0149;
            P3 = 0.0428;
            P2 = 0.3684;
        }else if (a==2){
            P5 = 0.0029;
            P4 = 0.0053;
            P3 = 0.0114;
            P2 = 0.1907;
        }else if (a==3){
            P5 = 0.0232;
            P4 = 0.0371;
            P3 = 0.0884;
            P2 = 0.3478;
        }   
    }else if(theta_x_layer < 25. * M_PI / 180.){
        if (a==0){
            P5 = 0.0234;
            P4 = 0.0522;
            P3 = 0.1259;
            P2 = 0.513;        
        }else if (a==1){
            P5 = 0.0101;
            P4 = 0.0176;
            P3 = 0.0494;
            P2 = 0.4712; 
        }else if (a==2){
            P5 = 0.0029;
            P4 = 0.0056;
            P3 = 0.0117;
            P2 = 0.2513;
        }else if (a==3){
            P5 = 0.0251;
            P4 = 0.0441;
            P3 = 0.1013;
            P2 = 0.4305;
        }
    }else if(theta_x_layer < 30. * M_PI / 180.){
        if (a==0){
            P5 = 0.0271;
            P4 = 0.0612;
            P3 = 0.1702;
            P2 = 0.5567;        
        }else if (a==1){
            P5 = 0.0111;
            P4 = 0.0206;
            P3 = 0.0769;
            P2 = 0.5514;
        }else if (a==2){
            P5 = 0.0029;
            P4 = 0.0058;
            P3 = 0.0168;
            P2 = 0.3199;
        }else if (a==3){
            P5 = 0.0278;
            P4 = 0.0538;
            P3 = 0.1285;
            P2 = 0.4843;
        }      
    }else if(theta_x_layer < 35. * M_PI / 180.){
        if (a==0){
            P5 = 0.0259;
            P4 = 0.0694;
            P3 = 0.257;
            P2 = 0.52;        
        }else if (a==1){
            P5 = 0.0141;
            P4 = 0.0297;
            P3 = 0.1454;
            P2 = 0.5669;
        }else if (a==2){
            P5 = 0.0027;
            P4 = 0.0059;
            P3 = 0.0315;
            P2 = 0.3779;
        }else if (a==3){
            P5 = 0.0254;
            P4 = 0.0608;
            P3 = 0.1811;
            P2 = 0.4799;
        }
    }else if(theta_x_layer < 40. * M_PI / 180.){
        if (a==0){
            P5 = 0.0295;
            P4 = 0.0943;
            P3 = 0.3709;
            P2 = 0.3965;        
        }else if (a==1){
            P5 = 0.0157;
            P4 = 0.0362;
            P3 = 0.2568;
            P2 = 0.4817;
        }else if (a==2){
            P5 = 0.0034;
            P4 = 0.0084;
            P3 = 0.0578;
            P2 = 0.3993;
        }else if (a==3){
            P5 = 0.0273;
            P4 = 0.0657;
            P3 = 0.2706;
            P2 = 0.4113;
        }     
    }else if(theta_x_layer < 45 * M_PI / 180.){
        if (a==0){
            P5 = 0.0449;
            P4 = 0.132;
            P3 = 0.4268;
            P2 = 0.2645;        
        }else if (a==1){
            P5 = 0.0216;
            P4 = 0.0528;
            P3 = 0.3638;
            P2 = 0.3542;
        }else if (a==2){
            P5 = 0.0042;
            P4 = 0.0107;
            P3 = 0.0871;
            P2 = 0.3532;
        }else if (a==3){
            P5 = 0.0377;
            P4 = 0.0893;
            P3 = 0.2998;
            P2 = 0.2852;
        }
    }
    if(myRand < P5) wiresOn = 5;
    else if(myRand < P4 + P5) wiresOn = 4;
    else if(myRand < P3 + P4 + P5) wiresOn = 3;
    else if(myRand < P2 + P3 + P4 + P5) wiresOn = 2;
    else wiresOn = 1;
    
    
    //Si se enciende un número par de hilos tomamos el valor medio
    // Obtengo el segundo hilo más cercano al punto de corte
    G4double interWire = 0;
    if(w <= wires[closestWire])  interWire = -distance / 2.;
    else interWire = distance / 2.;
    
    if(wiresOn == 2 || wiresOn == 4){
        val.push_back(wires[closestWire] + interWire);
    }else{
        val.push_back(wires[closestWire]);
    }
    
//    std::cout << "Hilos activados: "<< wiresOn << " Valor retornado: " << val[0] << std::endl; 
    
    if(closestWire == wires.size()-1) {
        val.push_back(fabs(wires[wires.size()-1]-wires[wires.size()-2]));
    } else {
        val.push_back(fabs(wires[closestWire+1]-wires[closestWire]));
    }
//    if(imin<1){        
//        std::cout << "index "<< closestWire;        
//        std::cout <<" w " << w;
//        std::cout <<" approxi " << approxi;
//        std::cout <<" imin " << imin;
//        std::cout <<" imax " << imax;
//        std::cout <<" value " << wires[closestWire];
//        std::cout << std::endl;
//    }
    return val;

}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Get the position of the closest wire                                 //
//----------------------------------------------------------------------//

std::vector<G4double> DriftChamberLayer::getWirePos(G4double w, G4int a) {
    std::vector<G4double> val;
    G4int closestWire = -1;
    G4double minimumSoFar = 20*CLHEP::m;
    std::vector<G4double> wires = (det->getWires())[a];
    if(wires.size() < 2) return val;
    G4double distance = wires[1]-wires[0];
    G4int approxi = (int)((w - wires[0])/distance);
    G4int imin = (int)(fmax(0, approxi-2));
    G4int imax = (int)(fmin(wires.size()-1, approxi+2));
    for(G4int i = imin; i <= imax; i++) {
        G4double dist = fabs(wires[i] - w);
        if(dist < minimumSoFar) {
            minimumSoFar = dist;
            closestWire = i;
        }
    }
    if(closestWire == 0 || closestWire == wires.size()-1) {
	    val.push_back(-999);
	    val.push_back(-999);
    }
    val.push_back(wires[closestWire]);
    if(closestWire == wires.size()-1) {
        val.push_back(fabs(wires[wires.size()-1]-wires[wires.size()-2]));
    } else {
        val.push_back(fabs(wires[closestWire+1]-wires[closestWire]));
    }
//    if(imin<1){        
//        std::cout << "index "<< closestWire;        
//        std::cout <<" w " << w;
//        std::cout <<" approxi " << approxi;
//        std::cout <<" imin " << imin;
//        std::cout <<" imax " << imax;
//        std::cout <<" value " << wires[closestWire];
//        std::cout << std::endl;
//    }
    return val;

}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Setting the resolution of the hit                                    //
//----------------------------------------------------------------------//
void DriftChamberLayer::SetStructure(Detector *a) {
    det = a;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Needed by the parent class                                           //
//----------------------------------------------------------------------//
void DriftChamberLayer::EndOfEvent(G4HCofThisEvent* /*HCE*/) {
    ;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//

