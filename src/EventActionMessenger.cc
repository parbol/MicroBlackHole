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
