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
/// \file analysis/AnaEx01/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
//
// $Id: HistoManager.cc 74272 2013-10-02 14:48:50Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
//#include "OutfileManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
{
  fFileName[0] = "Scream";
  fFactoryOn = false;
  
  // histograms
  for (G4int k=0; k<MaxHisto; k++) {
    fHistId[k] = 0;
    fHistPt[k] = 0;    
 /*}
  // ntuple
  for (G4int k=0; k<MaxNtCol; k++) {
    fNtColId[k] = 0;
  }*/  
}
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(2);
  G4String extension = analysisManager->GetFileType();
  fFileName[1] = fFileName[0] + "." + extension;
      
  // Create directories 
  analysisManager->SetHistoDirectoryName("histo");
  //analysisManager->SetNtupleDirectoryName("ntuple");
    
  // Open an output file
  //
  G4bool fileOpen = analysisManager->OpenFile(fFileName[0]);
  if (!fileOpen) {
    G4cout << "\n---> HistoManager::book(): cannot open " << fFileName[1] 
           << G4endl;
    return;
  }
  
  // create selected histograms
  //
  analysisManager->SetFirstHistoId(1);
  //analysisManager->SetFirstNtupleId(1);

  // fHistId[1] = analysisManager->CreateH1("1","Edep in absorber (MeV)",
  //                                             1000, 0., 50*GeV);
  // fHistPt[1] = analysisManager->GetH1(fHistId[1]);
                                           
  // fHistId[2] = analysisManager->CreateH1("2","Edep in Detector1 (MeV)",
  //                                             1000, 0., 150*keV);
  // fHistPt[2] = analysisManager->GetH1(fHistId[2]);
  
 // fHistId[3] = analysisManager->CreateH1("3","Edep in Detector1 (MeV)",
  //                                            1000, 0., 150*keV);
 // fHistPt[3] = analysisManager->GetH1(fHistId[3]);
                                           
  fHistId[1] = analysisManager->CreateH1("1","Edep in Detector1 pad (1,5) ",
                                              1000, 0., 500*keV);
  fHistPt[1] = analysisManager->GetH1(fHistId[1]);

  fHistId[2] = analysisManager->CreateH1("2","Edep in Detector1 pad (2,5)",
                                              1000, 0., 500*keV);
  fHistPt[2] = analysisManager->GetH1(fHistId[2]);

  fHistId[3] = analysisManager->CreateH1("3","Edep in Detector1 pad (3,5) ",
                                              1000, 0., 500*keV);
  fHistPt[3] = analysisManager->GetH1(fHistId[3]);

  fHistId[4] = analysisManager->CreateH1("4","Edep in Detector1 pad (4,5) ",
                                              1000, 0., 500*keV);
  fHistPt[4] = analysisManager->GetH1(fHistId[4]);

  fHistId[5] = analysisManager->CreateH1("5","Edep in Detector1 pad (5,5) ",
                                              1000, 0., 500*keV);
  fHistPt[5] = analysisManager->GetH1(fHistId[5]);

  fHistId[6] = analysisManager->CreateH1("6","Edep in chamber1",
                                              100000, 0., 50*MeV);
  fHistPt[6] = analysisManager->GetH1(fHistId[6]);

  fHistId[7] = analysisManager->CreateH1("7","Edep in chamber2",
                                              100000, 0., 50*MeV);
  fHistPt[7] = analysisManager->GetH1(fHistId[7]);

  fHistId[8] = analysisManager->CreateH1("8","Edep in chamber3",
                                              100000, 0., 50*MeV);
  fHistPt[8] = analysisManager->GetH1(fHistId[8]);

  fHistId[9] = analysisManager->CreateH1("9","Edep in chamber4",
                                              100000, 0., 50*MeV);
  fHistPt[9] = analysisManager->GetH1(fHistId[9]);

  fHistId[10] = analysisManager->CreateH1("10","Edep in chamber5",
                                              100000, 0., 50*MeV);
  fHistPt[10] = analysisManager->GetH1(fHistId[10]);

  fHistId[11] = analysisManager->CreateH1("11","Edep in chamber6",
                                              100000, 0., 50*MeV);
  fHistPt[11] = analysisManager->GetH1(fHistId[11]);

  fHistId[12] = analysisManager->CreateH1("12","Edep in all chambers",
                                              100000, 0., 50*MeV);
  fHistPt[12] = analysisManager->GetH1(fHistId[12]);

  fHistId[13] = analysisManager->CreateH1("13","Smeared Edep in Detector1 pad (1,5) ",
                                              1000, 0., 500*keV);
  fHistPt[13] = analysisManager->GetH1(fHistId[13]);

  fHistId[14] = analysisManager->CreateH1("14","Smeared Edep in Detector1 pad (2,5)",
                                              1000, 0., 500*keV);
  fHistPt[14] = analysisManager->GetH1(fHistId[14]);

  fHistId[15] = analysisManager->CreateH1("15","Smeared Edep in Detector1 pad (3,5) ",
                                              1000, 0., 500*keV);
  fHistPt[15] = analysisManager->GetH1(fHistId[15]);

  fHistId[16] = analysisManager->CreateH1("16","Smeared Edep in Detector1 pad (4,5) ",
                                              1000, 0., 500*keV);
  fHistPt[16] = analysisManager->GetH1(fHistId[16]);

  fHistId[17] = analysisManager->CreateH1("17","Smeared Edep in Detector1 pad (5,5) ",
                                              1000, 0., 500*keV);
  fHistPt[17] = analysisManager->GetH1(fHistId[17]);

  fHistId[18] = analysisManager->CreateH1("18","Smeared Edep in chamber1",
                                              100000, 0., 50*MeV);
  fHistPt[18] = analysisManager->GetH1(fHistId[18]);

  fHistId[19] = analysisManager->CreateH1("19","Smeared Edep in chamber2",
                                              100000, 0., 50*MeV);
  fHistPt[19] = analysisManager->GetH1(fHistId[19]);

  fHistId[20] = analysisManager->CreateH1("20","Smeared Edep in chamber3",
                                              100000, 0., 50*MeV);
  fHistPt[20] = analysisManager->GetH1(fHistId[20]);

  fHistId[21] = analysisManager->CreateH1("21","Smeared Edep in chamber4",
                                              100000, 0., 50*MeV);
  fHistPt[21] = analysisManager->GetH1(fHistId[21]);

  fHistId[22] = analysisManager->CreateH1("22","Smeared Edep in chamber5",
                                              100000, 0., 50*MeV);
  fHistPt[22] = analysisManager->GetH1(fHistId[22]);

  fHistId[23] = analysisManager->CreateH1("23","Smeared Edep in chamber6",
                                              100000, 0., 50*MeV);
  fHistPt[23] = analysisManager->GetH1(fHistId[23]);

  fHistId[24] = analysisManager->CreateH1("24","Smeared Edep in all chambers",
                                              100000, 0., 50*MeV);
  fHistPt[24] = analysisManager->GetH1(fHistId[24]);

  fHistId[25] = analysisManager->CreateH1("25","The smear factor. It should be gaussian.",
                                              1000, -10., 10.);
  fHistPt[25] = analysisManager->GetH1(fHistId[25]);

  fHistId[26] = analysisManager->CreateH1("26","Edep in chamber1 in adc",
                                              100000, 0., 50*keV * 1000. * 30.);
  fHistPt[26] = analysisManager->GetH1(fHistId[26]);

  fHistId[27] = analysisManager->CreateH1("27","Edep in chamber2 in adc",
                                              100000, 0., 50*keV * 1000. * 30.);
  fHistPt[27] = analysisManager->GetH1(fHistId[27]);

  fHistId[28] = analysisManager->CreateH1("28","Edep in chamber3 in adc",
                                              100000, 0., 50*keV * 1000. * 30.);
  fHistPt[28] = analysisManager->GetH1(fHistId[28]);

  fHistId[29] = analysisManager->CreateH1("29","Edep in chamber4 in adc",
                                              100000, 0., 50*keV * 1000. * 30.);
  fHistPt[29] = analysisManager->GetH1(fHistId[29]);

  fHistId[30] = analysisManager->CreateH1("30","Edep in chamber5 in adc",
                                              100000, 0., 50*keV * 1000. * 30.);
  fHistPt[30] = analysisManager->GetH1(fHistId[30]);

  fHistId[31] = analysisManager->CreateH1("31","Edep in chamber6 in adc",
                                              100000, 0., 50*keV * 1000. * 30.);
  fHistPt[31] = analysisManager->GetH1(fHistId[31]);

  fHistId[32] = analysisManager->CreateH1("32","Edep in all chambers in adc",
					      100000, 0., 50*keV * 1000. * 30.);
  fHistPt[32] = analysisManager->GetH1(fHistId[32]);

  fHistId[33] = analysisManager->CreateH1("33","Energy deposited in absorber1",
  					      100000, 0., 5000*MeV);
  fHistPt[33] = analysisManager->GetH1(fHistId[33]);

  fHistId[34] = analysisManager->CreateH1("34","Energy deposited in absorber2",
  					      100000, 0., 5000*MeV);
  fHistPt[34] = analysisManager->GetH1(fHistId[34]);

  fHistId[35] = analysisManager->CreateH1("35","Energy deposited in absorber3",
  					      100000, 0., 5000*MeV);
  fHistPt[35] = analysisManager->GetH1(fHistId[35]);

  fHistId[36] = analysisManager->CreateH1("36","Energy deposited in absorber4",
  					      100000, 0., 5000*MeV);
  fHistPt[36] = analysisManager->GetH1(fHistId[36]);

  fHistId[37] = analysisManager->CreateH1("37","Energy deposited in absorber5",
  					      100000, 0., 5000*MeV);
  fHistPt[37] = analysisManager->GetH1(fHistId[37]);

  fHistId[38] = analysisManager->CreateH1("38","Energy deposited in absorber6",
  					      100000, 0., 5000*MeV);
  fHistPt[38] = analysisManager->GetH1(fHistId[38]);

  fHistId[39] = analysisManager->CreateH1("39","Energy deposited in absorber7",
  					      100000, 0., 5000*MeV);
  fHistPt[39] = analysisManager->GetH1(fHistId[39]);

  fHistId[40] = analysisManager->CreateH1("40","Energy deposited in absorber1 (percent of Eincident)",
  					      1100,0.,110.);
  fHistPt[40] = analysisManager->GetH1(fHistId[40]);

  fHistId[41] = analysisManager->CreateH1("41","Energy deposited in absorber2 (percent of Eincident)",
  					      1100,0.,110.);
  fHistPt[41] = analysisManager->GetH1(fHistId[41]);

  fHistId[42] = analysisManager->CreateH1("42","Energy deposited in absorber3 (percent of Eincident)",
  					      1100,0.,110.);
  fHistPt[42] = analysisManager->GetH1(fHistId[42]);

  fHistId[43] = analysisManager->CreateH1("43","Energy deposited in absorber4 (percent of Eincident)",
  					      1100,0.,110.);
  fHistPt[43] = analysisManager->GetH1(fHistId[43]);

  fHistId[44] = analysisManager->CreateH1("44","Energy deposited in absorber5 (percent of Eincident)",
  					      1100,0.,110.);
  fHistPt[44] = analysisManager->GetH1(fHistId[44]);

  fHistId[45] = analysisManager->CreateH1("45","Energy deposited in absorber6 (percent of Eincident)",
  					      1100,0.,110.);
  fHistPt[45] = analysisManager->GetH1(fHistId[45]);

  fHistId[46] = analysisManager->CreateH1("46","Energy deposited in absorber7 (percent of Eincident)",
  					      1100,0.,110.);
  fHistPt[46] = analysisManager->GetH1(fHistId[46]);

  fHistId[47] = analysisManager->CreateH1("47","Edep in chamber1 (percent of Eincident)",
                                              1100,0.,110.);
  fHistPt[47] = analysisManager->GetH1(fHistId[47]);

  fHistId[48] = analysisManager->CreateH1("48","Edep in chamber2 (percent of Eincident)",
                                              1100,0.,110.);
  fHistPt[48] = analysisManager->GetH1(fHistId[48]);

  fHistId[49] = analysisManager->CreateH1("49","Edep in chamber3 (percent of Eincident)",
                                              1100,0.,110.);
  fHistPt[49] = analysisManager->GetH1(fHistId[49]);

  fHistId[50] = analysisManager->CreateH1("50","Edep in chamber4 (percent of Eincident)",
                                              1100,0.,110.);
  fHistPt[50] = analysisManager->GetH1(fHistId[50]);

  fHistId[51] = analysisManager->CreateH1("51","Edep in chamber5 (percent of Eincident)",
                                              1100,0.,110.);
  fHistPt[51] = analysisManager->GetH1(fHistId[51]);

  fHistId[52] = analysisManager->CreateH1("52","Edep in chamber6 (percent of Eincident)",
                                              1100,0.,110.);
  fHistPt[52] = analysisManager->GetH1(fHistId[52]);

  fHistId[53] = analysisManager->CreateH1("53","Energy loss in absorber1",
  					      100, 0. , 100.*mm);
  fHistPt[53] = analysisManager->GetH1(fHistId[53]);

  fHistId[54] = analysisManager->CreateH1("54","Energy loss in absorber2",
  					      100, 0. , 100.*mm);
  fHistPt[54] = analysisManager->GetH1(fHistId[54]);

  fHistId[55] = analysisManager->CreateH1("55","Energy loss in absorber3",
  					      100, 0. , 100.*mm);
  fHistPt[55] = analysisManager->GetH1(fHistId[55]);

  fHistId[56] = analysisManager->CreateH1("56","Energy loss in absorber4",
  					      100, 0. , 100.*mm);
  fHistPt[56] = analysisManager->GetH1(fHistId[56]);

  fHistId[57] = analysisManager->CreateH1("57","Energy loss in absorber5",
  					      100, 0. , 100.*mm);
  fHistPt[57] = analysisManager->GetH1(fHistId[57]);

  fHistId[58] = analysisManager->CreateH1("58","Energy loss in absorber6",
  					      100, 0. , 100.*mm);
  fHistPt[58] = analysisManager->GetH1(fHistId[58]);

  fHistId[59] = analysisManager->CreateH1("59","Energy loss in absorber7",
  					      100, 0. , 100.*mm);
  fHistPt[59] = analysisManager->GetH1(fHistId[59]);


  // Create 1st ntuple (id = 1)
  //    
 /* analysisManager->CreateNtuple("101", "Edep");
  fNtColId[0] = analysisManager->CreateNtupleDColumn("Eabs");
  fNtColId[1] = analysisManager->CreateNtupleDColumn("Egap");
  analysisManager->FinishNtuple();

  // Create 2nd ntuple (id = 2)
  //    
  analysisManager->CreateNtuple("102", "TrackL");
  fNtColId[2] = analysisManager->CreateNtupleDColumn("Labs");
  fNtColId[3] = analysisManager->CreateNtupleDColumn("Lgap");
  analysisManager->FinishNtuple();
  
  */fFactoryOn = true; /*      
  G4cout << "\n----> Histogram Tree is opened in " << fFileName[1] << G4endl;*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::save()
{
  if (fFactoryOn) {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();    
    analysisManager->Write();
    analysisManager->CloseFile();  
   // G4cout << "\n----> Histogram Tree is saved in " << fFileName[1] << G4endl;
      
    delete G4AnalysisManager::Instance();
    fFactoryOn = false;
  }                    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillHisto(G4int ih, G4double xbin)//, G4double weight)
{
  if (ih > MaxHisto) {
    G4cout << "---> warning from HistoManager::FillHisto() : histo " << ih
           << "does note xist; xbin= " << xbin /*<< " w= " << weight */<< G4endl;
    return;
  }

  if (fHistPt[ih]) fHistPt[ih]->fill(xbin);//, weight);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*void HistoManager::Normalize(G4int ih, G4double fac)
{
  if (ih >= MaxHisto) {
    G4cout << "---> warning from HistoManager::Normalize() : histo " << ih
           << "  fac= " << fac << G4endl;
    return;
  }

  if (fHistPt[ih]) fHistPt[ih]->scale(fac);
}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*void HistoManager::FillNtuple(G4double energyAbs, G4double energyGap,
                              G4double trackLAbs, G4double trackLGap)
{                
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // Fill 1st ntuple ( id = 1)
  analysisManager->FillNtupleDColumn(1,fNtColId[0], energyAbs);
  analysisManager->FillNtupleDColumn(1,fNtColId[1], energyGap);
  analysisManager->AddNtupleRow(1);  
  // Fill 2nd ntuple ( id = 2)
  analysisManager->FillNtupleDColumn(2,fNtColId[2], trackLAbs);
  analysisManager->FillNtupleDColumn(2,fNtColId[3], trackLGap);
  analysisManager->AddNtupleRow(2);  
}  */

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*void HistoManager::PrintStatistic()
{
  if(fFactoryOn) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
       << " EAbs : mean = " << G4BestUnit(fHistPt[1]->mean(), "Energy") 
               << " rms = " << G4BestUnit(fHistPt[1]->rms(),  "Energy") 
               << G4endl;
    G4cout                
       << " EGap : mean = " << G4BestUnit(fHistPt[2]->mean(), "Energy") 
               << " rms = " << G4BestUnit(fHistPt[2]->rms(),  "Energy") 
               << G4endl;
    G4cout 
       << " LAbs : mean = " << G4BestUnit(fHistPt[3]->mean(), "Length") 
               << " rms = " << G4BestUnit(fHistPt[3]->rms(),  "Length") 
               << G4endl;
    G4cout 
       << " LGap : mean = " << G4BestUnit(fHistPt[4]->mean(), "Length") 
               << " rms = " << G4BestUnit(fHistPt[4]->rms(),  "Length") 
               << G4endl;
  }
}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


