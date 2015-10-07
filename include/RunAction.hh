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
/// \file analysis/shared/include/RunAction.hh
/// \brief Definition of the RunAction class
//
//
// $Id: RunAction.hh 68015 2013-03-13 13:27:27Z gcosmo $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1
#include <fstream>
#include "G4UserRunAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;
class HistoManager; 

//class OutfileManager;

class RunAction : public G4UserRunAction
{
public:
  RunAction(HistoManager*);//,OutfileManager*);
  virtual ~RunAction();

  virtual void BeginOfRunAction(const G4Run*);
  virtual void   EndOfRunAction(const G4Run*);
    
  void fillPerEvent(G4double, G4double, G4double, G4double, G4double, 
		    G4double, G4double, G4double, G4double, G4double, 
		    G4double, G4double, G4double, G4double, G4double, 
		    G4double, G4double, G4double, G4double, G4double, 
		    G4double, G4double, G4double, G4double, G4double);
 
private:
  HistoManager* fHistoManager;
  //OutfileManager* fOutfileManager;

  G4double fSumEAbs1;
  G4double fSumEAbs2;
  G4double fSumEAbs3;
  G4double fSumEAbs4;
  G4double fSumEAbs5;
  G4double fSumEAbs6;
  G4double fSumEAbs7;
  G4double fSumEAl1;
  G4double fSumEAl2;
  G4double fSumEAl3;
  G4double fSumEAl4;
  G4double fSumEAl5;
  G4double fSumEAl6;
  G4double fSumEPCB1;
  G4double fSumEPCB2;
  G4double fSumEPCB3;
  G4double fSumEPCB4;
  G4double fSumEPCB5;
  G4double fSumEPCB6;
  G4double fSumEAr1;
  G4double fSumEAr2;
  G4double fSumEAr3;
  G4double fSumEAr4;
  G4double fSumEAr5;
  G4double fSumEAr6;
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

