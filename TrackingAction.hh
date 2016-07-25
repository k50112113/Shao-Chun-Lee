#ifndef TrackingAction_h
#define TrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class TrackingAction : public G4UserTrackingAction
{
  public:
	TrackingAction();
	virtual ~TrackingAction();

	virtual void PreUserTrackingAction(const G4Track* preTrack);
	virtual void PostUserTrackingAction(const G4Track* posTrack);
	  

};

#endif