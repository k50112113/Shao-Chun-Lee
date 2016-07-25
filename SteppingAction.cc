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
// $Id: SteppingAction.cc 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "Run.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include <string>

SteppingAction::SteppingAction(EventAction* EVE)
: G4UserSteppingAction(),eve(EVE)
{}

SteppingAction::~SteppingAction()
{}

void SteppingAction::SetRun(Run *arun){
	CurrentRun=arun;
}

void SteppingAction::SetNowWeighting(double nowweight){
  NowWeight=nowweight;
}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
 // G4cout << "SteppingAction" << G4endl;
  // get volume of the current step
  Volume = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();
  //G4cout << "getlogicalvolume" << G4endl;
  //G4cout << "************ " << Volume->GetInstanceID() << G4endl;
  /*
  char* str=NULL;
  int i=0 ;
  bool pass = false;
  strcpy(str,Volume->GetName().data());
  printf("%s\n",str);
  tmp=0;
  while(str[i]>='0' || str[i]<='9'){
     tmp*=10;
	 tmp+=(str[i]-'0');
	 i++;
	 pass=true;
  }
  if(!pass)return;
  G4cout << tmp << G4endl;
  // check if we are in scoring volume
  */
  sliindex=-1;
  if(Volume->GetName()=="5cm"){
	sliindex=0;
  }else if(Volume->GetName()=="11cm"){
    sliindex=1;
  }else if(Volume->GetName()=="13cm"){
    sliindex=2;
  }else if(Volume->GetName()=="15cm"){
    sliindex=3;
  }else if(Volume->GetName()=="17cm"){
    sliindex=4;
  }else if(Volume->GetName()=="19cm"){
    sliindex=5;
  }
  if(sliindex > -1){
	//G4cout<<sliindex<<"  "<<(int)(step->GetPreStepPoint()->GetPosition().getX()+120)<<"  "<<(int)(step->GetPreStepPoint()->GetPosition().getY()+120)<<G4endl;
    CurrentRun->AddEdep(step->GetTotalEnergyDeposit()*NowWeight,sliindex,(int)(step->GetPreStepPoint()->GetPosition().getX()+120),(int)(step->GetPreStepPoint()->GetPosition().getY()+120));
	eve->Setweightflag(sliindex,true);
  }
}

