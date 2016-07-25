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
// $Id: Run.hh 66536 2012-12-19 14:32:36Z ihrivnac $
//
/// \file Run.hh
/// \brief Definition of the Run class

#ifndef Run_h
#define Run_h 1
#include "DetectorConstruction.hh"
#include "G4Run.hh"
#include "globals.hh"

class G4Event;

/// Run class
///

class Run : public G4Run
{
  public:
    Run();
    virtual ~Run();

    // method from the base class
    virtual void Merge(const G4Run*);
    
    void AddEdep (G4double edep,int sli,int indexx,int indexy); 
	void Addtotwei(G4double wei,int sli);

    // get methods
    G4double GetEdep(int sli,int indexx,int indexy)  const { return fEdep[sli][indexx][indexy]; }
    G4double GetEdep2(int sli,int indexx,int indexy) const { return fEdep2[sli][indexx][indexy]; }
	G4double Gettotwei(int sli) const { return totWei[sli];}
  private:
    G4double  fEdep[Slices][MAXS][MAXS];
    G4double  fEdep2[Slices][MAXS][MAXS];
	G4double  totWei[Slices];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

