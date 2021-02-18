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

