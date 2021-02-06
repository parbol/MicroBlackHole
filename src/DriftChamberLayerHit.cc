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
// DriftChamberLayerHit.cc                                              //
//----------------------------------------------------------------------//
// Methods for the DriftChamberLayerHit class.                          //
// This means the sizes of the universe, the blast furnace and their    //
// materials.                                                           //
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


#include "DriftChamberLayerHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"

//----------------------------------------------------------------------//
// Declaration of the allocator                                         //
//----------------------------------------------------------------------//
G4Allocator<DriftChamberLayerHit> DriftChamberLayerHitAllocator;
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Constructor without indicating layer ID                              //
//----------------------------------------------------------------------//
DriftChamberLayerHit::DriftChamberLayerHit() {
    layerID = -1;
    time = 0.;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Constructor indicating layer ID                                      //
//----------------------------------------------------------------------//
DriftChamberLayerHit::DriftChamberLayerHit(G4int z) {
    layerID = z;
    time = 0.;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Destructor                                                           //
//----------------------------------------------------------------------//
DriftChamberLayerHit::~DriftChamberLayerHit() {
    ;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Constructor from another hit                                         //
//----------------------------------------------------------------------//
DriftChamberLayerHit::DriftChamberLayerHit(const DriftChamberLayerHit &right): G4VHit() {
    layerID = right.layerID;
    localPos = right.localPos;
    localMeas = right.localMeas;
    time = right.time;
    energy = right.energy;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Overloading = operator                                               //
//----------------------------------------------------------------------//
const DriftChamberLayerHit& DriftChamberLayerHit::operator=(const DriftChamberLayerHit &right) {
    layerID = right.layerID;
    localPos = right.localPos;
    localMeas = right.localMeas;
    time = right.time;
    energy = right.energy;
    return *this;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Overloading == operator                                               //
//----------------------------------------------------------------------//
int DriftChamberLayerHit::operator==(const DriftChamberLayerHit &/*right*/) const {
    return 0;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Associate measurements and units                                     //
//----------------------------------------------------------------------//
const std::map<G4String,G4AttDef>* DriftChamberLayerHit::GetAttDefs() const {
    G4bool isNew;
    std::map<G4String,G4AttDef>* store = G4AttDefStore::GetInstance("DriftChamberLayerHit",isNew);
    if (isNew) {
        G4String HitType("HitType");
        (*store)[HitType] = G4AttDef(HitType,"Hit Type","Physics","","G4String");

        G4String ID("ID");
        (*store)[ID] = G4AttDef(ID,"ID","Physics","","G4int");

        G4String Time("Time");
        (*store)[Time] = G4AttDef(Time,"Time","Physics","G4BestUnit","G4double");

        G4String Energy("Energy");
        (*store)[Energy] = G4AttDef(Energy,"Energy","Physics","G4BestUnit","G4double");

        G4String Pos("Pos");
        (*store)[Pos] = G4AttDef(Pos, "Position", "Physics","G4BestUnit","G4ThreeVector");

    }
    return store;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Create specials fuields inside the hit                               //
//----------------------------------------------------------------------//
std::vector<G4AttValue>* DriftChamberLayerHit::CreateAttValues() const {
    std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

    values->push_back(G4AttValue("HitType","DriftChamberLayerHit",""));

    values->push_back
    (G4AttValue("ID",G4UIcommand::ConvertToString(layerID),""));

    values->push_back
    (G4AttValue("Time",G4BestUnit(time,"Time"),""));

    values->push_back
    (G4AttValue("Energy",G4BestUnit(energy,"Energy"),""));

    values->push_back
    (G4AttValue("Pos",G4BestUnit(worldPos,"Length"),""));

    return values;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Printing the information of the hit                                  //
//----------------------------------------------------------------------//
void DriftChamberLayerHit::Print()
{
    G4cout << "\033[1;34m"
           << "*******************************Hit**********************************" << G4endl
           << "Layer: " << "\033[1;33m" << layerID << "\033[1;34m" << G4endl
           << "Time: " << "\033[1;33m" << time << "\033[1;34m" << G4endl
           << "Energy: " << "\033[1;33m" << energy/CLHEP::GeV << "\033[1;34m" << G4endl
           << "Local pos: " << "\033[1;33m" << localPos.x()/CLHEP::cm << " " << localPos.y()/CLHEP::cm << " " << localPos.z()/CLHEP::cm << "\033[1;34m" << G4endl
           << "Local pos: " << "\033[1;33m" << localMeas.x()/CLHEP::cm << " " << localMeas.y()/CLHEP::cm << " " << localMeas.z()/CLHEP::cm << "\033[1;34m" << G4endl
           << "********************************************************************" << "\033[0m" << G4endl;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


