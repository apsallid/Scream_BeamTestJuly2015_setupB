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
/// \file analysis/shared/include/EventAction.hh
/// \brief Definition of the EventAction class
//
//
// $Id: EventAction.hh 68015 2013-03-13 13:27:27Z gcosmo $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <fstream>

#include "TTree.h"
#include "TFile.h"

class RunAction;
class HistoManager;
class PrimaryGeneratorAction;
//class OutfileManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4int chambers =6;

class EventAction : public G4UserEventAction
{
public:
  //G4double fEnergyGap[10][10];
  EventAction(RunAction*, HistoManager*, PrimaryGeneratorAction*);
  virtual ~EventAction();
 
  virtual void  BeginOfEventAction(const G4Event*);
  virtual void    EndOfEventAction(const G4Event*);
  

  
  void AddAbs1(G4double de) {fEnergyAbs1 += de;};
  void AddAbs2(G4double de) {fEnergyAbs2 += de;};
  void AddAbs3(G4double de) {fEnergyAbs3 += de;};
  void AddAbs4(G4double de) {fEnergyAbs4 += de;};
  void AddAbs5(G4double de) {fEnergyAbs5 += de;};
  void AddAbs6(G4double de) {fEnergyAbs6 += de;};
  void AddAbs7(G4double de) {fEnergyAbs7 += de;};
  void AddAl1(G4double de) {fEnergyAl1 += de;};
  void AddAl2(G4double de) {fEnergyAl2 += de;};
  void AddAl3(G4double de) {fEnergyAl3 += de;};
  void AddAl4(G4double de) {fEnergyAl4 += de;};
  void AddAl5(G4double de) {fEnergyAl5 += de;};
  void AddAl6(G4double de) {fEnergyAl6 += de;};
  void AddAr1(G4double de, int i , int j)  {fEnergyAr1[i][j] += de;};
  void AddAr2(G4double de, int i , int j)  {fEnergyAr2[i][j] += de;};
  void AddAr3(G4double de, int i , int j)  {fEnergyAr3[i][j] += de;};
  void AddAr4(G4double de, int i , int j)  {fEnergyAr4[i][j] += de;};
  void AddAr5(G4double de, int i , int j)  {fEnergyAr5[i][j] += de;};
  void AddAr6(G4double de, int i , int j)  {fEnergyAr6[i][j] += de;};
  void AddPCB1(G4double de) {fEnergyPCB1 += de;};
  void AddPCB2(G4double de) {fEnergyPCB2 += de;};
  void AddPCB3(G4double de) {fEnergyPCB3 += de;};
  void AddPCB4(G4double de) {fEnergyPCB4 += de;};
  void AddPCB5(G4double de) {fEnergyPCB5 += de;};
  void AddPCB6(G4double de) {fEnergyPCB6 += de;};



  G4int ovfcounter;
  G4int counter;
  G4bool ovf;
  G4int  AdcValsAr1[10][10], AdcValsAr2[10][10], AdcValsAr3[10][10], AdcValsAr4[10][10], AdcValsAr5[10][10], AdcValsAr6[10][10];
  std::ofstream outfile;
 
    
private:
  RunAction*    fRunAct;
  HistoManager* fHistoManager;
  PrimaryGeneratorAction* fKin;
   
  G4double  fEnergyAbs1, fEnergyAbs2, fEnergyAbs3, fEnergyAbs4, fEnergyAbs5, fEnergyAbs6, fEnergyAbs7;
  G4double  fEnergyAl1, fEnergyAl2, fEnergyAl3, fEnergyAl4, fEnergyAl5, fEnergyAl6;
  G4double  fEnergyAr1[10][10], fEnergyAr2[10][10], fEnergyAr3[10][10], fEnergyAr4[10][10], fEnergyAr5[10][10], fEnergyAr6[10][10];
  G4double  fEnergyPCB1, fEnergyPCB2, fEnergyPCB3, fEnergyPCB4, fEnergyPCB5, fEnergyPCB6;

  //Energies
  G4double eneAr1, eneAr2, eneAr3, eneAr4, eneAr5, eneAr6, eneinallcha; 
  G4double eneAr1sm, eneAr2sm, eneAr3sm, eneAr4sm, eneAr5sm, eneAr6sm, eneinallchasm; 
  
  TFile *outF_;
  TTree *tree_;

  
                    
   G4int     fPrintModulo;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
