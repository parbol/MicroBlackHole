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
// EventActionMessanger.cc                                              //
//----------------------------------------------------------------------//
// Methods for EventActionMessanger.                                    //
// To do: Possibly we can skip this method.                             //
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//

#include "EventActionMessenger.hh"
#include "EventAction.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4ios.hh"

//----------------------------------------------------------------------//
// Constructor                                                          //
//----------------------------------------------------------------------//
EventActionMessenger::EventActionMessenger(EventAction * mpga) : target(mpga) {
    verboseCmd = new G4UIcmdWithAnInteger("/mydet/verbose",this);
    verboseCmd->SetGuidance("Verbose level for each event.");
    verboseCmd->SetGuidance(" Event summary will be displayed for every 'level' events.");
    verboseCmd->SetParameterName("level",true);
    verboseCmd->SetRange("level>=0");
    verboseCmd->SetDefaultValue(0);
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Destructor                                                           //
//----------------------------------------------------------------------//
EventActionMessenger::~EventActionMessenger() {
    delete verboseCmd;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// SetVerbosityLevel                                                    //
//----------------------------------------------------------------------//
void EventActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue) {
    if( command==verboseCmd )
    {
        target->SetVerbose(verboseCmd->GetNewIntValue(newValue));
    }
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Get Verbosity Level                                                  //
//----------------------------------------------------------------------//
G4String EventActionMessenger::GetCurrentValue(G4UIcommand * command) {
    G4String cv;
    if( command==verboseCmd )
    {
        cv = verboseCmd->ConvertToString(target->GetVerbose());
    }

    return cv;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//