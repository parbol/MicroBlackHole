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
// EventAction.h                                                        //
//----------------------------------------------------------------------//
// Class for filling the hits in the proper placer.                     //
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//
#ifndef EventAction_h
#define EventAction_h 1


#include "G4UserEventAction.hh"
#include "globals.hh"
#include <vector>

#include <chrono>
#include <ctime> 

#include "ConfigurationGeometry.hh"

#ifdef G4ANALYSIS_USE
#include "AIDA/AIDA.h"
using namespace AIDA;
#endif


class EventActionMessenger;


class EventAction : public G4UserEventAction {
public:
    EventAction(ConfigurationGeometry *);
    virtual ~EventAction();

public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);
    void EndOfEventAction_old(const G4Event*);
    void SetNumberOfEvents(G4int);
private:
    void UpdateProgressBar();

private:
    std::vector<G4int> DHCID;
    std::vector<G4int> DHCID_1;
    std::vector<G4int> DHCID_2;

    EventActionMessenger *messenger;
    ConfigurationGeometry *geom;
    G4int verboseLevel;
    
    G4int numberOfEvents, event_progress;
    
private:
//    typedef std::chrono::steady_clock Clock;
//    typedef Clock::time_point time_point;
//    typedef Clock::period period;
//    typedef std::chrono::duration<float, period> duration;
//    std::uint64_t _total_ticks;
//    std::uint64_t _ticks_occurred;
//    time_point _begin;
    std::chrono::high_resolution_clock::time_point lastEvent_time;
    std::chrono::high_resolution_clock::time_point thisEvent_time;
    float mean_deltaTime, elapsedTime;
    
    

public:
    inline void SetVerbose(G4int val) {
        verboseLevel = val;
    }
    inline G4int GetVerbose() const {
        return verboseLevel;
    }
};


#endif
