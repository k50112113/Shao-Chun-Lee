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
// $Id: B1PrimaryGeneratorAction.cc 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file B1PrimaryGeneratorAction.cc
/// \brief Implementation of the B1PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"
#include "SteppingAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(SteppingAction *STP)
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0),stp(STP)
{
  //**************SOBP**************
  FILE *fp = fopen("SOBPdata.txt","r") ;
  int counting=0 ; 
  double totalw = 0;
  fscanf(fp,"%d",&totalenergygroup);
  for(index = 0 ; index < totalenergygroup ; index++){
    fscanf(fp,"%lf",&energygroup[index]) ;
  }
  index=0;indexi=0;
  fscanf(fp,"%lf",&xpos[index]);
  while(~fscanf(fp,"%lf",&xpos[index])){
	  if(xpos[index] > -2147483646){
	    fscanf(fp,"%lf",&ypos[index]);
		energy[index]=energygroup[indexi];
	    index++;
	  }else{
	    indexi++;
		if(indexi >= 28)break;
	  }
  }
  index=0;
  while(~fscanf(fp,"%lf",&weighting[index])){
	  if(weighting[index] > -2147483646){
		totalw+=weighting[index];
	    index++;
	  }
  }
  counting=index;
  /*
  for(index=0;index < counting;index++){
    printf("%d\t%lf\t%lf\t%lf\t%lf\n",index,xpos[index],ypos[index],weighting[index],energy[index]);
  }
  */
  printf("energy group = %d | total = %d | total weight = %lf\n",totalenergygroup,counting,totalw) ;
  total=counting;
  fclose(fp);
  
  FILE *fpp = fopen("stddevdata.txt","r");
  index = 0 ;
  while(~fscanf(fpp,"%lf %lf %lf %lf %lf",&calE[index],&delE[index],&xstdd[index],&ystdd[index],&stdd[index])){
    index++ ;
  }
  caltot = index ;
  printf("caltotal = %d\n",caltot) ;
  fclose(fpp) ;
  for(index = 0 ; index < counting ; index++){
	for(indexi = 0 ; indexi < caltot - 1 ; indexi++){
		if(calE[indexi] <= energy[index] && energy[index] <= calE[indexi + 1]){
			energy[index] = energy[index]+delE[indexi] ;
			groupindex[index] = indexi;
			break ;
		}
	}
  }
  index = 0;
  weighting[index]=1;
  newcount=0;
  system("pause");
  //**************SOBP**************

  G4int n_particle = 1 ;
  fParticleGun  = new G4ParticleGun(n_particle);
  
  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  fproton = particleTable->FindParticle(particleName="proton");
    fParticleGun->SetParticleDefinition(fproton);
  //printf("Enter beam mean energy(MeV) : ") ;
  //double MeanEnergy;
  //scanf("%lf",&MeanEnergy);
  fParticleGun->SetParticleEnergy(100.5*MeV);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,-220*cm));

  srand((unsigned)time(NULL)) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	
  if(newcount % 5000 == 0){ 
//	  index++;
	  printf("No. %d ~ %d\n",newcount,newcount+4999) ;
//	  printf("Now shooting index = %d | weighting = %lf | Energy = %lf MeV | Standard Deviation = %lf | Spot Size (%.2lfmm,%.2lfmm)\n"
//	  ,index,weighting[index],energy[index],stdd[groupindex[index]],xstdd[groupindex[index]],ystdd[groupindex[index]]);
  }
  newcount++;
  index = (int)(10641*G4RandFlat::shoot()) ;
  stp->SetNowWeighting(weighting[index]);
  Energy = G4RandGauss::shoot(energy[index],stdd[groupindex[index]]) ;
  fParticleGun->SetParticleEnergy(Energy*MeV);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(xpos[index]/10.,ypos[index]/10.,215));
  px=G4RandGauss::shoot(0,xstdd[groupindex[index]]);
  py=G4RandGauss::shoot(0,ystdd[groupindex[index]]);
  fParticleGun->SetParticlePosition(G4ThreeVector(px*mm,py*mm,(-20 - 200)*cm));
  fParticleGun->SetParticlePosition(G4ThreeVector(0*mm,0*mm,(-20 - 200)*cm));
  
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

