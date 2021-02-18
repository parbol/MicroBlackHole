#include "RunAction.hh"

#include "Analysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"



//----------------------------------------------------------------------//
// Constructor                                                          //
//----------------------------------------------------------------------//
RunAction::RunAction(std::string name, ConfigurationGeometry *gem_) {
    nameOfOutputFile = name;
    gem = gem_;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Destructor                                                           //
//----------------------------------------------------------------------//
RunAction::~RunAction() {}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// We set up the columns when the run starts                            //
//----------------------------------------------------------------------//
void RunAction::BeginOfRunAction(const G4Run* run) {

    G4AnalysisManager* man = G4AnalysisManager::Instance();
    man->OpenFile(nameOfOutputFile.c_str());

    // Create ntuple
    man->CreateNtuple("hits", "The hits");
    for(G4int a = 0; a < gem->getDetector1()->getNLayers(); a++) {
        char name_posx[20];
        sprintf(name_posx, "posX_1%d", a);
        char name_posy[20];
        sprintf(name_posy, "posY_1%d", a);
        char name_posz[20];
        sprintf(name_posz, "posZ_1%d", a);
        char name_measx[20];
        sprintf(name_measx, "measX_1%d", a);
        char name_measy[20];
        sprintf(name_measy, "measY_1%d", a);
        char name_measz[20];
        sprintf(name_measz, "measZ_1%d", a);
        char name_layerID[20];
        sprintf(name_layerID, "layerID_1%d", a);
        char name_e[20];
        sprintf(name_e, "e_1%d", a);
        man->CreateNtupleDColumn(name_posx);
        man->CreateNtupleDColumn(name_posy);
        man->CreateNtupleDColumn(name_posz);
        man->CreateNtupleDColumn(name_measx);
        man->CreateNtupleDColumn(name_measy);
        man->CreateNtupleDColumn(name_measz);
        man->CreateNtupleIColumn(name_layerID);
        man->CreateNtupleDColumn(name_e);
    }

    for(G4int a = 0; a < gem->getDetector2()->getNLayers(); a++) {
        char name_posx[20];
        sprintf(name_posx, "posX_2%d", a);
        char name_posy[20];
        sprintf(name_posy, "posY_2%d", a);
        char name_posz[20];
        sprintf(name_posz, "posZ_2%d", a);
        char name_measx[20];
        sprintf(name_measx, "measX_2%d", a);
        char name_measy[20];
        sprintf(name_measy, "measY_2%d", a);
        char name_measz[20];
        sprintf(name_measz, "measZ_2%d", a);
        char name_layerID[20];
        sprintf(name_layerID, "layerID_2%d", a);
        char name_e[20];
        sprintf(name_e, "e_2%d", a);
        man->CreateNtupleDColumn(name_posx);
        man->CreateNtupleDColumn(name_posy);
        man->CreateNtupleDColumn(name_posz);
        man->CreateNtupleDColumn(name_measx);
        man->CreateNtupleDColumn(name_measy);
        man->CreateNtupleDColumn(name_measz);
        man->CreateNtupleIColumn(name_layerID);
        man->CreateNtupleDColumn(name_e);
    }
    man->FinishNtuple();
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// At the end we store                                                  //
//----------------------------------------------------------------------//
void RunAction::EndOfRunAction(const G4Run* aRun)
{
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();
    delete G4AnalysisManager::Instance();
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//

