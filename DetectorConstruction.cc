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
// $Id: DetectorConstruction.cc 75117 2013-10-28 09:38:37Z gcosmo $
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include <stdio.h>
#include <stdlib.h>
#include "G4UserLimits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(NULL),fStepLimit(NULL)
{
	fScoringVolume = new G4LogicalVolume*[Slices];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete [] fScoringVolume;
  delete fStepLimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Envelope parameters
  //

  G4Element* elH  = new G4Element("Hydrogen","H" , 1., 1.01*g/mole);
  G4Element* elO  = new G4Element("Oxygen"  ,"O" , 8., 16.00*g/mole);
  G4Material* H2O = new G4Material("Water",1.000*g/cm3,2);
  H2O->AddElement(elH, 2);
  H2O->AddElement(elO, 1);

  G4Material* water = nist->FindOrBuildMaterial("G4_WATER") ;

  double ExEng=75. ;
  
  //printf("%lf eV\n",H2O->GetIonisation()->GetMeanExcitationEnergy()* 1000000) ;
  //printf("%lf eV\n",water->GetIonisation()->GetMeanExcitationEnergy()* 1000000) ;
  //printf("Enter Mean Excitation Energy(eV) : ") ;
  //scanf("%lf",&ExEng) ;
  if(ExEng > 0){
    H2O->GetIonisation()->SetMeanExcitationEnergy(ExEng*eV) ;
    water->GetIonisation()->SetMeanExcitationEnergy(ExEng*eV) ;
  }
  

  //H2O->GetIonisation()->SetMeanExcitationEnergy(75.0*eV) ;
  
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4Material* vac = nist->FindOrBuildMaterial("G4_Galactic");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       600*cm,600*cm,600*cm);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        vac,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
  

  //     
  // Env Watertank
  // 
  
  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        500*cm, 500*cm, 20*cm); //its size

  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        H2O,             //its materialsi
                        "Envelope");         //its name
   
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,0),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  fStepLimit = new G4UserLimits() ; 
  G4double MINKINE=0.1*MeV;
  //printf("SetUserMinEkine(MeV) : ") ;
  //scanf("%lf",&MINKINE) ;
  fStepLimit->SetMaxAllowedStep(1*mm);
  fStepLimit->SetUserMinEkine(MINKINE);

  G4Box* abox;
  G4double SliceWidth=3*mm;

  abox = new G4Box("box",12*cm,12*cm,SliceWidth/2.);	  
  fScoringVolume[0] = new G4LogicalVolume(abox,H2O,"5cm",0,0,0,true);
  new G4PVPlacement(0,G4ThreeVector(0*cm,0*cm,-15*cm),fScoringVolume[0],"5cm",logicEnv,false,0,true);
  
  abox = new G4Box("box",12*cm,12*cm,SliceWidth/2.);	  
  fScoringVolume[1] = new G4LogicalVolume(abox,H2O,"11cm",0,0,0,true);
  new G4PVPlacement(0,G4ThreeVector(0*cm,0*cm,-9*cm),fScoringVolume[1],"11cm",logicEnv,false,0,true);

  abox = new G4Box("box",12*cm,12*cm,SliceWidth/2.);	  
  fScoringVolume[2] = new G4LogicalVolume(abox,H2O,"13cm",0,0,0,true);
  new G4PVPlacement(0,G4ThreeVector(0*cm,0*cm,-7*cm),fScoringVolume[2],"13cm",logicEnv,false,0,true);

  abox = new G4Box("box",12*cm,12*cm,SliceWidth/2.);	  
  fScoringVolume[3] = new G4LogicalVolume(abox,H2O,"15cm",0,0,0,true);
  new G4PVPlacement(0,G4ThreeVector(0*cm,0*cm,-5*cm),fScoringVolume[3],"15cm",logicEnv,false,0,true);

  abox = new G4Box("box",12*cm,12*cm,SliceWidth/2.);	  
  fScoringVolume[4] = new G4LogicalVolume(abox,H2O,"17cm",0,0,0,true);
  new G4PVPlacement(0,G4ThreeVector(0*cm,0*cm,-3*cm),fScoringVolume[4],"17cm",logicEnv,false,0,true);

  abox = new G4Box("box",12*cm,12*cm,SliceWidth/2.);	  
  fScoringVolume[5] = new G4LogicalVolume(abox,H2O,"19cm",0,0,0,true);
  new G4PVPlacement(0,G4ThreeVector(0*cm,0*cm,-1*cm),fScoringVolume[5],"19cm",logicEnv,false,0,true);
  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
