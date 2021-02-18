#include "EventAction.hh"
#include "EventActionMessenger.hh"

#ifdef G4ANALYSIS_USE
#include "AnalysisManager.hh"
#endif

#include "Analysis.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

#include "DriftChamberLayerHit.hh"

#include "ConfigurationGeometry.hh"

#include <math.h>

#include <iostream>

#include <chrono>
#include <ctime> 

#define PI 3.141592653589793238462643383279502884



//----------------------------------------------------------------------//
// Constructor                                                          //
//----------------------------------------------------------------------//
EventAction::EventAction(ConfigurationGeometry *geom_) {

    geom = geom_;
    G4String colName;
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    G4int numberOfDetectors = SDman->GetCollectionCapacity();
    SDman->ListTree();
    G4int N_layers_det_1 = geom->getDetector1()->getNLayers();
    for( int c_layer = 0; c_layer < N_layers_det_1; c_layer++){
        DHCID_1.push_back(SDman->GetCollectionID("Detector1_Chamber_" + std::to_string(c_layer)+"/HitsCollection_Detector1_"+std::to_string(c_layer)));
    }
      
    G4int N_layers_det_2 = geom->getDetector2()->getNLayers();
    for( int c_layer = 0; c_layer < N_layers_det_2; c_layer++){
        DHCID_2.push_back(SDman->GetCollectionID("Detector2_Chamber_" + std::to_string(c_layer)+"/HitsCollection_Detector2_"+std::to_string(c_layer)));
    }
    
//    DHCID.push_back(SDman->GetCollectionID("Chamber1/HitsCollection1"));
//    DHCID.push_back(SDman->GetCollectionID("Chamber2/HitsCollection2"));

    messenger = new EventActionMessenger(this);
    
    numberOfEvents = 0;
    event_progress = 0;
    mean_deltaTime = 0;
    elapsedTime = 0;

}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Destructor                                                           //
//----------------------------------------------------------------------//
EventAction::~EventAction() {
    delete messenger;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//



//----------------------------------------------------------------------//
// Begin of the event                                                   //
//----------------------------------------------------------------------//
void EventAction::BeginOfEventAction(const G4Event*) {}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// End of the event -> the action goes here                             //
//----------------------------------------------------------------------//
void EventAction::EndOfEventAction(const G4Event* evt) {


    
    //Actualizamos la barra de progreso con cada evento.
    UpdateProgressBar();   
    
    
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
    std::vector<DriftChamberLayerHitsCollection* > DHC1;
    std::vector<DriftChamberLayerHitsCollection* > DHC2;
    //G4SDManager* SDman = G4SDManager::GetSDMpointer();
    
    G4int N_layers_det_1 = geom->getDetector1()->getNLayers();
    G4int N_layers_det_2 = geom->getDetector2()->getNLayers();

    if(HCE) {
        DriftChamberLayerHitsCollection *a = 0;     
        for( int c_layer = 0; c_layer < N_layers_det_1; c_layer++){
            a = (DriftChamberLayerHitsCollection*)(HCE->GetHC(DHCID_1.at(c_layer)));
            DHC1.push_back(a);
        }
        DriftChamberLayerHitsCollection *b = 0;             

        for( int c_layer = 0; c_layer < N_layers_det_2; c_layer++){
            b = (DriftChamberLayerHitsCollection*)(HCE->GetHC(DHCID_2.at(c_layer)));
            DHC2.push_back(b);
        }
    }
//    G4double E1=0;G4double E2=0;
    
    for(G4int iParticle = 0; iParticle < evt->GetPrimaryVertex(0)->GetNumberOfParticle(); iParticle++) {

        G4PrimaryParticle*   primary = evt->GetPrimaryVertex(0)->GetPrimary(iParticle);
        if(primary->GetParticleDefinition()->GetParticleName() != "microblackhole") continue;
        //if(primary->GetParticleDefinition()->GetParticleName() != "mu-" && primary->GetParticleDefinition()->GetParticleName() != "mu+") continue;
	G4int hits_1 = 0;
        for (int c_layer = 0; c_layer < N_layers_det_1; c_layer++){
            if(DHC1.at(c_layer)) {
                G4int n_hit = DHC1.at(c_layer)->entries();
                //if(n_hit < geom->getDetector1()->getNLayers()/2.0) return;
                if(n_hit == 0) return;
                for(G4int hit = 0; hit < n_hit; hit++) {
                    DriftChamberLayerHit* aHit = (*(DHC1.at(c_layer)))[hit];
                    G4ThreeVector pos = aHit->GetLocalPos();
                    G4ThreeVector meas = aHit->GetLocalMeas();
                    G4ThreeVector errormeas = aHit->GetLocalMeasError();
                    G4double e = aHit->GetEnergy();
                    G4int layerID = aHit->GetLayerID();
                    man->FillNtupleDColumn(8*layerID + 0, pos.x()/CLHEP::cm);
                    man->FillNtupleDColumn(8*layerID + 1, pos.y()/CLHEP::cm);
                    man->FillNtupleDColumn(8*layerID + 2, pos.z()/CLHEP::cm);
                    man->FillNtupleDColumn(8*layerID + 3, meas.x()/CLHEP::cm);
                    man->FillNtupleDColumn(8*layerID + 4, meas.y()/CLHEP::cm);
                    man->FillNtupleDColumn(8*layerID + 5, meas.z()/CLHEP::cm);
                    man->FillNtupleIColumn(8*layerID + 6, layerID);
                    man->FillNtupleDColumn(8*layerID + 7, e/CLHEP::MeV);
		    hits_1++;
//                    std::cout << "Added hit from det. 1"<<std::endl;
                }
            }
        }
        G4int hits_2 = 0;
        for (int c_layer = 0; c_layer < N_layers_det_2; c_layer++){
            if(DHC2.at(c_layer)) {
                G4int n_hit = DHC2.at(c_layer)->entries();
                for(G4int hit = 0; hit < n_hit; hit++) {
                    DriftChamberLayerHit* aHit = (*(DHC2.at(c_layer)))[hit];
                    G4ThreeVector pos = aHit->GetLocalPos();
                    G4ThreeVector meas = aHit->GetLocalMeas();
                    G4ThreeVector errormeas = aHit->GetLocalMeasError();
                    G4double e = aHit->GetEnergy();
//                    E1=e;
                    G4int layerID = aHit->GetLayerID();
                    man->FillNtupleDColumn(geom->getDetector1()->getNLayers() * 8 + 8*layerID + 0, pos.x()/CLHEP::cm);
                    man->FillNtupleDColumn(geom->getDetector1()->getNLayers() * 8 + 8*layerID + 1, pos.y()/CLHEP::cm);
                    man->FillNtupleDColumn(geom->getDetector1()->getNLayers() * 8 + 8*layerID + 2, pos.z()/CLHEP::cm);
                    man->FillNtupleDColumn(geom->getDetector1()->getNLayers() * 8 + 8*layerID + 3, meas.x()/CLHEP::cm);
                    man->FillNtupleDColumn(geom->getDetector1()->getNLayers() * 8 + 8*layerID + 4, meas.y()/CLHEP::cm);
                    man->FillNtupleDColumn(geom->getDetector1()->getNLayers() * 8 + 8*layerID + 5, meas.z()/CLHEP::cm);
                    man->FillNtupleIColumn(geom->getDetector1()->getNLayers() * 8 + 8*layerID + 6, layerID);
                    man->FillNtupleDColumn(geom->getDetector1()->getNLayers() * 8 + 8*layerID + 7, e/CLHEP::MeV);
		    hits_2++;
//                    std::cout << "Added hit from det. 2"<<std::endl;
                }
            }
        }
        // Solo imponemos que las trazas deban pasar por detector 1 y 2
        if(hits_1 == N_layers_det_1 && hits_2 == N_layers_det_2) man->AddNtupleRow(); 
    }


    //G4TrajectoryContainer *container = evt->GetTrajectoryContainer();
    //for(size_t i = 0; i < container->size(); i++) {
    //    if (((*(container))[i])->GetParticleName() == "mu-" || ((*(container))[i])->GetParticleName() == "mu+") {
    //        ((*(container))[i])->ShowTrajectory(); 
    //    }
//
//    }
//    evt->Print();

}

//----------------------------------------------------------------------//
//----------------------------------------------------------------------//

//----------------------------------------------------------------------//
// End of the event -> the action goes here                             //
//----------------------------------------------------------------------//
void EventAction::SetNumberOfEvents(G4int n) {
    numberOfEvents = n;

}

//----------------------------------------------------------------------//
//----------------------------------------------------------------------//




//----------------------------------------------------------------------//
// Update Progress Bar -> show progress events processed / total events //
//----------------------------------------------------------------------//
void EventAction::UpdateProgressBar(){
    
    event_progress++;
    thisEvent_time = std::chrono::system_clock::now();
    
    float progress = (float)event_progress/(float)numberOfEvents;
    int barWidth = 10;

    std::cout <<std::setprecision(2)<< "[" << std::setw(2);
//    int pos = barWidth * progress;
//    for (int i = 0; i < barWidth; ++i) {
//        if (i < pos) std::cout << "=";
//        else if (i == pos) std::cout << ">";
//        else std::cout << " ";
//    }
    if(event_progress == 1) std::cout << std::fixed << int(progress * 100.0) << "%] \r";
    else{
        std::chrono::duration<double> elapsed_seconds = thisEvent_time-lastEvent_time;
        elapsedTime += elapsed_seconds.count();
        float remainingTime = (numberOfEvents - event_progress) * elapsed_seconds.count()/60.;
        mean_deltaTime = elapsed_seconds.count()*1/event_progress + mean_deltaTime * (event_progress-1)/event_progress;
        
        std::cout << std::fixed 
                << int(progress * 100.0) << "%] Remaining: "  
                << (numberOfEvents - event_progress) * mean_deltaTime /60. 
                << " min Elapsed: " <<  elapsedTime/60. << " min\r";
        
        
    }
    
    std::cout.flush();
    if(int(progress) == 100) 
    {
        std::cout.unsetf(std::ios::fixed | std::ios::scientific);//Devuelvo el formato a su estado por defecto
        std::cout<< std::setprecision(0) << std::endl; //Cierro la linea
    }
    lastEvent_time = thisEvent_time;
        
}

//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// End of the event -> the action goes here                             //
//----------------------------------------------------------------------//
void EventAction::EndOfEventAction_old(const G4Event* evt) {

    G4AnalysisManager* man = G4AnalysisManager::Instance();
    G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
    std::vector<DriftChamberLayerHitsCollection* > DHC1;
    std::vector<DriftChamberLayerHitsCollection* > DHC2;
    //G4SDManager* SDman = G4SDManager::GetSDMpointer();

    if(HCE) {
        DriftChamberLayerHitsCollection *a = 0;
        a = (DriftChamberLayerHitsCollection*)(HCE->GetHC(DHCID.at(0)));
        DHC1.push_back(a);
        DriftChamberLayerHitsCollection *b = 0;
        b = (DriftChamberLayerHitsCollection*)(HCE->GetHC(DHCID.at(1)));
        DHC2.push_back(b);
    }

    for(G4int iParticle = 0; iParticle < evt->GetPrimaryVertex(0)->GetNumberOfParticle(); iParticle++) {

        G4PrimaryParticle*   primary = evt->GetPrimaryVertex(0)->GetPrimary(iParticle);
        if(primary->GetParticleDefinition()->GetParticleName() != "mu-" && primary->GetParticleDefinition()->GetParticleName() != "mu+") continue;

        if(DHC1.at(0)) {
            G4int n_hit = DHC1.at(0)->entries();
            //if(n_hit < geom->getDetector1()->getNLayers()/2.0) return;
            if(n_hit == 0) return;
            for(G4int hit = 0; hit < n_hit; hit++) {
                DriftChamberLayerHit* aHit = (*(DHC1.at(0)))[hit];
                G4ThreeVector pos = aHit->GetLocalPos();
                G4ThreeVector meas = aHit->GetLocalMeas();
                G4ThreeVector errormeas = aHit->GetLocalMeasError();
                G4double e = aHit->GetEnergy();
                G4int layerID = aHit->GetLayerID();
                man->FillNtupleDColumn(8*layerID + 0, pos.x()/CLHEP::cm);
                man->FillNtupleDColumn(8*layerID + 1, pos.y()/CLHEP::cm);
                man->FillNtupleDColumn(8*layerID + 2, pos.z()/CLHEP::cm);
                man->FillNtupleDColumn(8*layerID + 3, meas.x()/CLHEP::cm);
                man->FillNtupleDColumn(8*layerID + 4, meas.y()/CLHEP::cm);
                man->FillNtupleDColumn(8*layerID + 5, meas.z()/CLHEP::cm);
                man->FillNtupleIColumn(8*layerID + 6, layerID);
                man->FillNtupleDColumn(8*layerID + 7, e/CLHEP::MeV);
            }
        }

        if(DHC2.at(0)) {
            G4int n_hit = DHC2.at(0)->entries();
            for(G4int hit = 0; hit < n_hit; hit++) {
                DriftChamberLayerHit* aHit = (*(DHC2.at(0)))[hit];
                G4ThreeVector pos = aHit->GetLocalPos();
                G4ThreeVector meas = aHit->GetLocalMeas();
                G4ThreeVector errormeas = aHit->GetLocalMeasError();
                G4double e = aHit->GetEnergy();
                G4int layerID = aHit->GetLayerID();
                man->FillNtupleDColumn(geom->getDetector1()->getNLayers() * 8 + 8*layerID + 0, pos.x()/CLHEP::cm);
                man->FillNtupleDColumn(geom->getDetector1()->getNLayers() * 8 + 8*layerID + 1, pos.y()/CLHEP::cm);
                man->FillNtupleDColumn(geom->getDetector1()->getNLayers() * 8 + 8*layerID + 2, pos.z()/CLHEP::cm);
                man->FillNtupleDColumn(geom->getDetector1()->getNLayers() * 8 + 8*layerID + 3, meas.x()/CLHEP::cm);
                man->FillNtupleDColumn(geom->getDetector1()->getNLayers() * 8 + 8*layerID + 4, meas.y()/CLHEP::cm);
                man->FillNtupleDColumn(geom->getDetector1()->getNLayers() * 8 + 8*layerID + 5, meas.z()/CLHEP::cm);
                man->FillNtupleIColumn(geom->getDetector1()->getNLayers() * 8 + 8*layerID + 6, layerID);
                man->FillNtupleDColumn(geom->getDetector1()->getNLayers() * 8 + 8*layerID + 7, e/CLHEP::MeV);
            }
        }
        if(DHC1.at(0)) man->AddNtupleRow();
    }

}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//
