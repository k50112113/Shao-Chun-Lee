//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: EventAction.cc 75117 2013-10-28 09:38:37Z gcosmo $
//
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"
#include "Run.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

EventAction::EventAction(): G4UserEventAction()
{
} 

EventAction::~EventAction()
{
  delete pga;
}

void EventAction::BeginOfEventAction(const G4Event*)
{    
  //G4cout << "BeginOfEventAction" << G4endl;
	for(int i = 0 ; i < Slices ; i++)weightflag[i]=false;
}

void EventAction::EndOfEventAction(const G4Event*)
{   
  //G4cout << "EndOfEventAction" << G4endl;
  // accumulate statistics in Run
	
  Run* run 
    = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  //G4cout << "Adding dose with weight = " << pga->GetNowWeighting() << G4endl;
  for(int index = 0 ; index < Slices ; index++){
	if(weightflag[index]==true){
	  run->Addtotwei(pga->GetNowWeighting(),index);
	}
  }
}

void EventAction::Setweightflag(int index,bool fg){
	weightflag[index]=fg;
}
