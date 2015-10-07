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
/// \file analysis/shared/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
// $Id: SteppingAction.cc 68015 2013-03-13 13:27:27Z gcosmo $
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"
//#include "OutfileManager.hh"
#include "G4Step.hh"
#include "Randomize.hh"
#include "HistoManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det,
			       EventAction* evt)
: G4UserSteppingAction(), 
  fDetector(det), fEventAction(evt)                                   
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  // get volume of the current step
  G4VPhysicalVolume* volume 
  = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  
  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit();
 // G4double Earray[10][10]; 
  /*G4double stepl = 0.;
  if (aStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0.)
    stepl = aStep->GetStepLength();*/

  //Bragg curve
  G4StepPoint* prePoint  = aStep->GetPreStepPoint();
  G4StepPoint* postPoint = aStep->GetPostStepPoint();
   
  G4double x1 = prePoint->GetPosition().x();
  G4double x2 = postPoint->GetPosition().x();  
  //The absorber1 starts at x=0 so we do not need the third term here 
  G4double x  = x1 + fabs(G4UniformRand()*(x2-x1)); // + 0.5*(fDetector->GetAbsorber1Thickness());
  //The absorber2 starts at x= fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness 
  G4double startpos_abs2 = fDetector->GetAbsorber1Thickness() + fDetector->GetAl1Thickness() + fDetector->GetAr1Thickness() + fDetector->GetPCB1Thickness() + fDetector->GetAirThickness(); 
  G4double xx  = x1 + fabs(G4UniformRand()*(x2-x1)) - startpos_abs2; 

  //The absorber3 starts at 
  //x= fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + fAr2Thickness + fPCB2Thickness + fAirThickness
  G4double startpos_abs3 = fDetector->GetAbsorber1Thickness() + fDetector->GetAl1Thickness() + fDetector->GetAr1Thickness() + fDetector->GetPCB1Thickness() + fDetector->GetAirThickness() + fDetector->GetAbsorber2Thickness() + fDetector->GetAl2Thickness() + fDetector->GetAr2Thickness() + fDetector->GetPCB2Thickness() + fDetector->GetAirThickness(); 

  G4double xx3  = x1 + fabs(G4UniformRand()*(x2-x1)) - startpos_abs3; 
 
  //The absorber4 starts at 
  //x= fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + fAr2Thickness + fPCB2Thickness + + fAirThickness + fAbsorber3Thickness + fAl3Thickness + fAr3Thickness + fPCB3Thickness + fAirThickness
  G4double startpos_abs4 = fDetector->GetAbsorber1Thickness() + fDetector->GetAl1Thickness() + fDetector->GetAr1Thickness() + fDetector->GetPCB1Thickness() + fDetector->GetAirThickness() + fDetector->GetAbsorber2Thickness() + fDetector->GetAl2Thickness ()+ fDetector->GetAr2Thickness() + fDetector->GetPCB2Thickness() + fDetector->GetAirThickness() + fDetector->GetAbsorber3Thickness() + fDetector->GetAl3Thickness() + fDetector->GetAr3Thickness() + fDetector->GetPCB3Thickness() + fDetector->GetAirThickness(); 

  G4double xx4  = x1 + fabs(G4UniformRand()*(x2-x1)) - startpos_abs4; 
 //The absorber5 starts at 
  //x= fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + fAr2Thickness + fPCB2Thickness + fAirThickness + fAbsorber3Thickness + fAl3Thickness + fAr3Thickness + fPCB3Thickness + fAirThickness + fAbsorber4Thickness + fAl4Thickness + fAr4Thickness + fPCB4Thickness + fAirThickness 
  G4double startpos_abs5 = fDetector->GetAbsorber1Thickness() + fDetector->GetAl1Thickness() + fDetector->GetAr1Thickness() + fDetector->GetPCB1Thickness() + fDetector->GetAirThickness() + fDetector->GetAbsorber2Thickness() + fDetector->GetAl2Thickness()+ fDetector->GetAr2Thickness() + fDetector->GetPCB2Thickness() + fDetector->GetAirThickness() + fDetector->GetAbsorber3Thickness() + fDetector->GetAl3Thickness() + fDetector->GetAr3Thickness() + fDetector->GetPCB3Thickness() + fDetector->GetAirThickness() + fDetector->GetAbsorber4Thickness() + fDetector->GetAl4Thickness() + fDetector->GetAr4Thickness() + fDetector->GetPCB4Thickness() + fDetector->GetAirThickness(); 

  G4double xx5  = x1 + fabs(G4UniformRand()*(x2-x1)) - startpos_abs5; 
 //The absorber6 starts at 
  //x=fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + fAr2Thickness + fPCB2Thickness + fAirThickness + fAbsorber3Thickness + fAl3Thickness + fAr3Thickness + fPCB3Thickness + fAirThickness + fAbsorber4Thickness + fAl4Thickness + fAr4Thickness + fPCB4Thickness + fAirThickness + fAbsorber5Thickness + fAl5Thickness + Ar5Thickness + fPCB5Thickness + fAirThickness 
  G4double startpos_abs6 = fDetector->GetAbsorber1Thickness() + fDetector->GetAl1Thickness() + fDetector->GetAr1Thickness() + fDetector->GetPCB1Thickness() + fDetector->GetAirThickness() + fDetector->GetAbsorber2Thickness ()+ fDetector->GetAl2Thickness() + fDetector->GetAr2Thickness() + fDetector->GetPCB2Thickness() + fDetector->GetAirThickness() + fDetector->GetAbsorber3Thickness() + fDetector->GetAl3Thickness() + fDetector->GetAr3Thickness() + fDetector->GetPCB3Thickness() + fDetector->GetAirThickness() + fDetector->GetAbsorber4Thickness() + fDetector->GetAl4Thickness() + fDetector->GetAr4Thickness() + fDetector->GetPCB4Thickness() + fDetector->GetAirThickness()+ fDetector->GetAbsorber5Thickness() + fDetector->GetAl5Thickness() + fDetector->GetAr5Thickness() + fDetector->GetPCB5Thickness() + fDetector->GetAirThickness(); 

  G4double xx6  = x1 + fabs(G4UniformRand()*(x2-x1)) - startpos_abs6; 
 //The absorber7 starts at 
  //x=fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + fAr2Thickness + fPCB2Thickness + fAirThickness + fAbsorber3Thickness + fAl3Thickness + fAr3Thickness + fPCB3Thickness + fAirThickness + fAbsorber4Thickness + fAl4Thickness + fAr4Thickness + fPCB4Thickness + fAirThickness + fAbsorber5Thickness + fAl5Thickness + Ar5Thickness + fPCB5Thickness + fAirThickness + fAbsorber6Thickness + fAl6Thickness + fAr6Thickness + fPCB6Thickness + fAirThickness 
  G4double startpos_abs7 = fDetector->GetAbsorber1Thickness() + fDetector->GetAl1Thickness() + fDetector->GetAr1Thickness() + fDetector->GetPCB1Thickness() + fDetector->GetAirThickness() + fDetector->GetAbsorber2Thickness() + fDetector->GetAl2Thickness() + fDetector->GetAr2Thickness() + fDetector->GetPCB2Thickness() + fDetector->GetAirThickness() + fDetector->GetAbsorber3Thickness() + fDetector->GetAl3Thickness() + fDetector->GetAr3Thickness() + fDetector->GetPCB3Thickness() + fDetector->GetAirThickness() + fDetector->GetAbsorber4Thickness() + fDetector->GetAl4Thickness() + fDetector->GetAr4Thickness() + fDetector->GetPCB4Thickness() + fDetector->GetAirThickness() + fDetector->GetAbsorber5Thickness() + fDetector->GetAl5Thickness() + fDetector->GetAr5Thickness() + fDetector->GetPCB5Thickness() + fDetector->GetAirThickness() + fDetector->GetAbsorber6Thickness() + fDetector->GetAl6Thickness() + fDetector->GetAr6Thickness() + fDetector->GetPCB6Thickness() + fDetector->GetAirThickness(); 

  G4double xx7  = x1 + fabs(G4UniformRand()*(x2-x1)) - startpos_abs7; 
      // std::cout << "The absorber1 size is " << fDetector->GetAbsorber1Thickness() << std::endl;

  // G4double y1 = prePoint->GetPosition().y();
  // G4double y2 = postPoint->GetPosition().y();  
  // G4double y  = y1 + fabs(G4UniformRand()*(y2-y1));

  // G4double z1 = prePoint->GetPosition().z();
  // G4double z2 = postPoint->GetPosition().z();  
  // G4double z  = z1 + fabs(G4UniformRand()*(z2-z1));

  // G4double radius = sqrt( pow(y,2.) + pow(z,2.) ); // in mm
  
  // std::cout << "y position is  " <<  y << std::endl;
  // std::cout << "z position is  " <<  z << std::endl;
  // std::cout << "The radius is  " <<  sqrt( pow(y,2.) + pow(z,2.) ) << std::endl;

  //HistoManager* histo;
  //fHistoManager->FillHisto(33, x, edep);//Energy loss 
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
        
  if (volume == fDetector->GetAbsorber1()) {
    analysisManager->FillH1(53, x, edep);
  }
  if (volume == fDetector->GetAbsorber2()) {
    analysisManager->FillH1(54, xx, edep);
  }
  if (volume == fDetector->GetAbsorber3()) {
    analysisManager->FillH1(55, xx3, edep);
  }
  if (volume == fDetector->GetAbsorber4()) {
    analysisManager->FillH1(56, xx4, edep);
  }
  if (volume == fDetector->GetAbsorber5()) {
    analysisManager->FillH1(57, xx5, edep);
  }
  if (volume == fDetector->GetAbsorber6()) {
    analysisManager->FillH1(58, xx6, edep);
  }
  if (volume == fDetector->GetAbsorber7()) {
    analysisManager->FillH1(59, xx7, edep);
  }

  if (volume == fDetector->GetAbsorber1()) {fEventAction->AddAbs1(edep);}

  if (volume == fDetector->GetAbsorber2()) {fEventAction->AddAbs2(edep);}

  if (volume == fDetector->GetAbsorber3()) {fEventAction->AddAbs3(edep);}
				       				 
  if (volume == fDetector->GetAbsorber4()) {fEventAction->AddAbs4(edep);}
				       				 
  if (volume == fDetector->GetAbsorber5()) {fEventAction->AddAbs5(edep);}
				       				 
  if (volume == fDetector->GetAbsorber6()) {fEventAction->AddAbs6(edep);}
				       				 
  if (volume == fDetector->GetAbsorber7()) {fEventAction->AddAbs7(edep);}

  if (volume == fDetector->GetAl1()) {fEventAction->AddAl1(edep);}

  if (volume == fDetector->GetAl2()) {fEventAction->AddAl2(edep);}

  if (volume == fDetector->GetAl3()) {fEventAction->AddAl3(edep);}
				 			  
  if (volume == fDetector->GetAl4()) {fEventAction->AddAl4(edep);}
				 			  
  if (volume == fDetector->GetAl5()) {fEventAction->AddAl5(edep);}
				 			  
  if (volume == fDetector->GetAl6()) {fEventAction->AddAl6(edep);}

  if (volume == fDetector->GetPCB1()) {fEventAction->AddPCB1(edep);}

  if (volume == fDetector->GetPCB2()) {fEventAction->AddPCB2(edep);}

  if (volume == fDetector->GetPCB3()) {fEventAction->AddPCB3(edep);}
				  			    
  if (volume == fDetector->GetPCB4()) {fEventAction->AddPCB4(edep);}
				  			    
  if (volume == fDetector->GetPCB5()) {fEventAction->AddPCB5(edep);}
				  			    
  if (volume == fDetector->GetPCB6()) {fEventAction->AddPCB6(edep);}

  for(G4int j=0;j<10;j++){
    for(G4int i=0;i<10;i++){

      //Disk of radius 2 mm
      //bool radius = ( sqrt( pow(y,2.) + pow(z,2.) ) < 2. ) ; // in mm

      // if ((volume == fDetector->GetGap(i,j))  && (k==0) && radius ) fEventAction->AddGap(edep,k,i,j);

      // if ((volume == fDetector->GetGap1(i,j)) && (k==1) && radius ) fEventAction->AddGap(edep,k,i,j);

      // if ((volume == fDetector->GetGap2(i,j)) && (k==2) && radius ) fEventAction->AddGap(edep,k,i,j);
 
      // if ((volume == fDetector->GetGap3(i,j)) && (k==3) && radius ) fEventAction->AddGap(edep,k,i,j);

      // if ((volume == fDetector->GetGap4(i,j)) && (k==4) && radius ) fEventAction->AddGap(edep,k,i,j);

      // if ((volume == fDetector->GetGap5(i,j)) && (k==5) && radius ) fEventAction->AddGap(edep,k,i,j);
 
      if ( (volume == fDetector->GetAr1(i,j)) ) fEventAction->AddAr1(edep,i,j);

      if ( (volume == fDetector->GetAr2(i,j)) ) fEventAction->AddAr2(edep,i,j);

      if ( (volume == fDetector->GetAr3(i,j)) ) fEventAction->AddAr3(edep,i,j);
				       				    
      if ( (volume == fDetector->GetAr4(i,j)) ) fEventAction->AddAr4(edep,i,j);
				       				    
      if ( (volume == fDetector->GetAr5(i,j)) ) fEventAction->AddAr5(edep,i,j);
				       				    
      if ( (volume == fDetector->GetAr6(i,j)) ) fEventAction->AddAr6(edep,i,j);

    }//for i
  }//for j

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
