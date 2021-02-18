#ifndef DriftChamberLayerHit_h
#define DriftChamberLayerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

class DriftChamberLayerHit : public G4VHit
{
public:

    DriftChamberLayerHit();
    DriftChamberLayerHit(G4int z);
    virtual ~DriftChamberLayerHit();
    DriftChamberLayerHit(const DriftChamberLayerHit &right);
    const DriftChamberLayerHit& operator=(const DriftChamberLayerHit &right);
    int operator==(const DriftChamberLayerHit &right) const;

    inline void *operator new(size_t);
    inline void operator delete(void *aHit);

    inline float x();
    inline float y();

    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
    virtual std::vector<G4AttValue>* CreateAttValues() const;
    virtual void Print();

private:
    G4int layerID;
    G4double time;
    G4double energy;
    G4ThreeVector localPos;
    G4ThreeVector worldPos;
    G4ThreeVector localMeas;
    G4ThreeVector localMeasError;
    G4double efficiency;

public:
    inline void SetLayerID(G4int z) {
        layerID = z;
    }
    inline G4int GetLayerID() const {
        return layerID;
    }
    inline void SetTime(G4double t) {
        time = t;
    }
    inline G4double GetTime() const {
        return time;
    }
    inline void SetEnergy(G4double t) {
        energy = t;
    }
    inline G4double GetEnergy() const {
        return energy;
    }
    inline void SetLocalPos(G4ThreeVector xyz) {
        localPos = xyz;
    }
    inline G4ThreeVector GetLocalPos() const {
        return localPos;
    }
    inline void SetLocalMeas(G4ThreeVector xyz) {
        localMeas = xyz;
    }
    inline G4ThreeVector GetLocalMeas() const {
        return localMeas;
    }
    inline void SetLocalMeasError(G4ThreeVector xyz) {
        localMeasError = xyz;
    }
    inline G4ThreeVector GetLocalMeasError() const {
        return localMeasError;
    }

};


//----------------------------------------------------------------------//
// Collections of this                                                  //
//----------------------------------------------------------------------//
typedef G4THitsCollection<DriftChamberLayerHit> DriftChamberLayerHitsCollection;

extern G4Allocator<DriftChamberLayerHit> DriftChamberLayerHitAllocator;
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Redefinition of new operator                                         //
//----------------------------------------------------------------------//
inline void* DriftChamberLayerHit::operator new(size_t) {
    void* aHit;
    aHit = (void*)DriftChamberLayerHitAllocator.MallocSingle();
    return aHit;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Redefinition of delete operator                                      //
//----------------------------------------------------------------------//
inline void DriftChamberLayerHit::operator delete(void* aHit) {
    DriftChamberLayerHitAllocator.FreeSingle((DriftChamberLayerHit*) aHit);
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//

#endif


