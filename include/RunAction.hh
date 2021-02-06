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
// RunAction.hh                                                         //
//----------------------------------------------------------------------//
// Class that handles the storage of variables to be stored.            //
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "ConfigurationGeometry.hh"


class G4Run;

class RunAction : public G4UserRunAction
{
public:
    RunAction(std::string, ConfigurationGeometry *);
    virtual ~RunAction();

    void BeginOfRunAction(const G4Run*);
    void   EndOfRunAction(const G4Run*);

    std::string nameOfOutputFile, mode;
    G4int nDet;
    ConfigurationGeometry *gem;

};


#endif

