//----------------------------------------------//
//-------- Microblackhole simulation -----------//
//----------------------------------------------//
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "RunAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "MicroBlackHole.hh"
#include "MicroBlackHolePhysics.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <sys/ioctl.h>
#include <fcntl.h>
#include <linux/kd.h>


//----------------------------------------------------------------------//
// Methods defined in this file                                         //
//----------------------------------------------------------------------//
//Get options from the command line
bool getOptions(int , char **, G4String &, G4String &, G4int &, G4long &, G4String &);
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Start of the program                                                 //
//----------------------------------------------------------------------//
int main(int argc,char** argv) {

    G4String nameOfInputFile;
    G4String nameOfOutputFile;
    G4int    numberOfEvents;
    G4long   randomSeed = 0;
    G4String visual="";
        
    if(!getOptions(argc, argv, nameOfInputFile, nameOfOutputFile, numberOfEvents, randomSeed, visual)) {
        G4cerr << "\033[1;31m" << "Usage: ./MicroBlackHole --input nameOfInputFile --output outputfile --number numberOfEvents --seed seed --vis visualManager"  << "\033[0m" << G4endl;
        return -1;
    }

    //Number of events all right?
    if(numberOfEvents < 1) {
        G4cerr << "\033[1;31m" << "The number of events must be a positive integer greater than 0" << "\033[0m" << G4endl;
        return -1;
    }

    ConfigurationGeometry *geomConf = new ConfigurationGeometry(nameOfInputFile);
    if(!geomConf->isGood()) {
        G4cerr << "\033[1;31m" << "Problems in the configuration geometry file" << "\033[0m" << G4endl;
        return -1;
    }

    //Initializing runManager
    G4RunManager* runManager = new G4RunManager;

    runManager->SetVerboseLevel(4);   

    DetectorConstruction *myDetector = new DetectorConstruction(geomConf);
    if(myDetector == NULL) {
        G4cerr << "\033[1;31m" << "Problems in the construction of the detector" << "\033[0m" << G4endl;
        return -1;
    }

    runManager->SetUserInitialization(myDetector);

    // Using LHC standard list of physics processes
    G4PhysListFactory factory;
    G4String plName = "FTFP_BERT";
    G4VModularPhysicsList* phys = factory.GetReferencePhysList(plName);

    // Adding physics for the MicroBlackHole
    MicroBlackHolePhysics * theBH = new MicroBlackHolePhysics();
    phys->RegisterPhysics(theBH);

    runManager->SetUserInitialization(phys);

    PrimaryGeneratorAction *myPrimaryGeneratorAction = new PrimaryGeneratorAction(geomConf, "", randomSeed);
    if(myPrimaryGeneratorAction == NULL) {
        G4cerr << "\033[1;31m" << "Problems in PrimaryGeneratorAction" << "\033[0m" << G4endl;
        return -1;
    }

    runManager->Initialize();

    runManager->SetUserAction(myPrimaryGeneratorAction);

    RunAction *myRunAction = new RunAction(nameOfOutputFile, geomConf);
    if(myRunAction == NULL) {
        G4cerr << "\033[1;31m" << "Problems in RunAction" << "\033[0m" << G4endl;
        return -1;
    }

    runManager->SetUserAction(myRunAction);

    EventAction *myEventAction = new EventAction(geomConf);
    myEventAction->SetNumberOfEvents(numberOfEvents);
    if(myEventAction == NULL) {
        G4cerr << "\033[1;31m" << "Problems in EventAction" << "\033[0m" << G4endl;
        return -1;
    }

    runManager->SetUserAction(myEventAction);


    #ifdef G4VIS_USE
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
    //visManager->RegisterGraphicsSystem (new G4HepRep);
    #endif

     if(visual == "heprep"){
        G4UImanager* UI = G4UImanager::GetUIpointer();
        UI->ApplyCommand("/run/verbose 2");
        UI->ApplyCommand("/vis/scene/create A01Output.heprep");
        UI->ApplyCommand("/vis/open HepRepXML");
        UI->ApplyCommand("/vis/drawVolume world");
        UI->ApplyCommand("/vis/scene/add/volume");
        //UI->ApplyCommand("");
        //UI->ApplyCommand("");
        UI->ApplyCommand("/vis/viewer/flush");
        //UI->ApplyCommand("/vis/verbose");
        //UI->ApplyCommand("/vis/heprep/setFileName archivo.hep");
    
        //   runManager->BeamOn(numberOfEvents);
    
    }else if(visual == "opengl"){
        G4UImanager* UI = G4UImanager::GetUIpointer();
        G4UIExecutive* ui_e = new G4UIExecutive(argc, argv);
    
        UI->ApplyCommand("/run/verbose 2");
        UI->ApplyCommand("/vis/scene/create A01Output.heprep");


        UI->ApplyCommand("/vis/open OGL 600x600-0+0");

        UI->ApplyCommand("/vis/drawVolume world");

        UI->ApplyCommand("/vis/scene/add/volume");
        UI->ApplyCommand("/vis/scene/add/trajectories");
        UI->ApplyCommand("/vis/scene/endOfEventAction accumulate");
        //UI->ApplyCommand("");
        UI->ApplyCommand("/vis/viewer/flush");
        UI->ApplyCommand("/vis/viewer/set/style surface");

        
        ui_e->SessionStart();
        delete ui_e;
    }

    
    runManager->BeamOn(numberOfEvents);


    #ifdef G4VIS_USE
    delete visManager;
    #endif
    delete runManager;
    delete geomConf;
    G4cout << "The program finished successfully" << std::endl;
    return 0;

}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//



//----------------------------------------------------------------------//
// This method will put the command line input in the variables.        //
//----------------------------------------------------------------------//
bool getOptions(int argc, char **argv, G4String &nameOfInputFile, G4String &nameOfOutputFile, G4int &numberOfEvents, G4long &randomSeed, G4String &visual) {

    int option_iterator;
    int option_counter = 0;
    bool moreoptions = true;

    while (moreoptions) {
        static struct option long_options[] = {
            /* These options set a flag. */
            {"input",     required_argument, 0, 'i'},
            {"output",    required_argument, 0, 'o'},
            {"seed",      required_argument, 0, 's'},
            {"number",    required_argument, 0, 'n'},
            {"vis",       required_argument, 0, 'v'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        option_iterator = getopt_long(argc, argv, "d:", long_options, &option_index);
        if (option_iterator == -1) {
            moreoptions = false;
        } else {
            option_counter++;
            switch (option_iterator) {
            case 0:
                if (long_options[option_index].flag != 0)
                    break;
                if (optarg)
                    break;
            case 'i':
                nameOfInputFile = (G4String) optarg;
                break;
            case 'o':
                nameOfOutputFile = (G4String) optarg;
                break;
            case 'n':
                numberOfEvents = (G4int) atoi(optarg);
                break;
            case 's':
                randomSeed = (G4long) atoi(optarg);
                break;
            case 'v':
                visual = (G4String) optarg;
                break;
            case '?':
                return false;
                break;
            default:
                return false;
            }
        }
    }

    if (option_counter == 0) {
        return false;
    }
    return true;

}



