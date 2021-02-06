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
// EventActionMessanger.hh                                              //
//----------------------------------------------------------------------//
// Handles the output from the Event Action.                            //
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


#ifndef EventActionMessenger_h
#define EventActionMessenger_h 1

class EventAction;
class G4UIcmdWithAnInteger;

#include "G4UImessenger.hh"
#include "globals.hh"

class EventActionMessenger: public G4UImessenger
{
public:
    EventActionMessenger(EventAction* mpga);
    ~EventActionMessenger();

public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

private:
    EventAction* target;

private: //commands
    G4UIcmdWithAnInteger*  verboseCmd;

};

#endif


