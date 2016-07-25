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
// $Id: Run.cc 66536 2012-12-19 14:32:36Z ihrivnac $
//
/// \file Run.cc
/// \brief Implementation of the Run class

#include "Run.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run()
: G4Run()
{
  int isl,ix,iy;
  for(isl = 0 ; isl < Slices ; isl++){
	totWei[isl]=0;
    for(ix = 0 ; ix < MAXS ; ix++){
		for(iy = 0 ; iy < MAXS ; iy++){
	      fEdep[isl][ix][iy]=0;
		}
	}
  }
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{} 
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::AddEdep (G4double edep,int sli,int indexx,int indexy)
{
    fEdep[sli][indexx][indexy]  += edep;
	//fEdep2[sli][indexx][indexy] += edep*edep;
	//G4cout << sli << " " << indexx << "," << indexy << " GET edep " << edep << G4endl;	
}

void Run::Addtotwei (G4double wei,int sli){
    totWei[sli]+=wei;
	//G4cout << sli << " GET wei " << wei  << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


