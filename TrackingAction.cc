#include "TrackingAction.hh"

TrackingAction::TrackingAction() : G4UserTrackingAction()

{
}

TrackingAction::~TrackingAction()
{
}

void TrackingAction::PreUserTrackingAction(const G4Track* preTrack){
	/*
	if(preTrack->GetCreatorModelID() != -1){
			//G4cout << "Particle : " << preTrack->GetDefinition()->GetParticleName() << G4endl;
			//G4cout << "Pre GetCreatorModelName : " << preTrack->GetCreatorModelName() << G4endl;
			//G4cout << "Pre GetCreatorModelID : " << preTrack->GetCreatorModelID() << G4endl;
			if(preTrack->GetDefinition()->GetParticleName() == "e-" && preTrack->GetCreatorModelID() == 2){
			  G4cout << "in the track" << G4endl;
			  system("pause");
			}
		}
	*/
}

void TrackingAction::PostUserTrackingAction(const G4Track* posTrack){
	
}