#ifndef StackingAction_h
#define StackingAction_h 1

#include "G4UserStackingAction.hh"
#include "globals.hh"
#include "G4Track.hh"

class StackingAction : public G4UserStackingAction
{
	public:
		StackingAction(){
		  G4cout << "StackingAcition started" << G4endl;
		}
		virtual ~StackingAction(){
		  G4cout << "StackingAcition ended" << G4endl;
		}

	public:
		virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack){
			/*
			if(aTrack->GetDefinition()->GetParticleName() == "e-"){
			  G4cout << "electron is found : " << aTrack->GetCreatorModelName() << G4endl ;
			  system("pause");
			  return fKill;
			}*/
			return fUrgent;
		}

};

#endif