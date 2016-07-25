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
/// \file electromagnetic/TestEm7/src/PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//
// $Id: PhysicsList.cc 82280 2014-06-13 14:45:31Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"

#include "G4StepLimiterPhysics.hh"
#include "G4EmStandardPhysics_option3.hh"

#include "G4DecayPhysics.hh"

#include "G4HadronElasticPhysics.hh"
#include "G4HadronInelasticQBBC.hh"

#include "G4VPhysicsConstructor.hh"

#include "G4SystemOfUnits.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList() : G4VUserPhysicsList()
{
  SetVerboseLevel(1);
  
  SetParticleCuts(34*cm,"gamma");
  SetParticleCuts(143.1*um,"e-");
  SetParticleCuts(143.1*um,"e+");
  SetParticleCuts(1*mm,"proton");

  fEmPhysicsList = new G4EmStandardPhysics_option3();   
  fDecPhysicsList = new G4DecayPhysics();
  fHadronPhys.push_back(new G4HadronElasticPhysics());
  fHadronPhys.push_back(new G4HadronInelasticQBBC());
  fUserSteplimitPhysicsList = new G4StepLimiterPhysics();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
  delete fEmPhysicsList;
  delete fDecPhysicsList;
  delete fUserSteplimitPhysicsList;
  for(size_t i=0; i<fHadronPhys.size(); i++) {delete fHadronPhys[i];}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  fDecPhysicsList->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  AddTransportation();
  fEmPhysicsList->ConstructProcess();
  fDecPhysicsList->ConstructProcess();
  fUserSteplimitPhysicsList->ConstructProcess();
  for(size_t i=0; i<fHadronPhys.size(); i++) {
    fHadronPhys[i]->ConstructProcess();
  }
}