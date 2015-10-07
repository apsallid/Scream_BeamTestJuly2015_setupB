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
/// \file analysis/shared/src/EventAction.cc
/// \brief Implementation of the EventAction class
//
//
// $Id: EventAction.cc 68015 2013-03-13 13:27:27Z gcosmo $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
//#include "OutfileManager.hh"
#include "RunAction.hh"
#include "HistoManager.hh"
//#include <fstream>
#include "G4Event.hh"
#include "TRandom3.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



EventAction::EventAction(RunAction* run, HistoManager* histo, PrimaryGeneratorAction* kin)//, OutfileManager* out)
 :G4UserEventAction(),
  fRunAct(run),fHistoManager(histo),fKin(kin),//,fOutfileManager(out),
  fEnergyAbs1(0.), 
  fEnergyAbs2(0.), 
  fEnergyAbs3(0.), 
  fEnergyAbs4(0.), 
  fEnergyAbs5(0.), 
  fEnergyAbs6(0.), 
  fEnergyAbs7(0.), 
  fEnergyAl1(0.), 
  fEnergyAl2(0.), 
  fEnergyAl3(0.), 
  fEnergyAl4(0.), 
  fEnergyAl5(0.), 
  fEnergyAl6(0.), 
  fEnergyPCB1(0.), 
  fEnergyPCB2(0.),
  fEnergyPCB3(0.),
  fEnergyPCB4(0.),
  fEnergyPCB5(0.),
  fEnergyPCB6(0.),
  //, fEnergyGap1(0.), //fEnergyGap(0.),
  // fTrackLAbs(0.), fTrackLGap(0.),
  fPrintModulo(0)                             
{
  fPrintModulo = 100; 
  outF_=TFile::Open("ScreamTree.root","RECREATE");
  outF_->cd();
  tree_=new TTree("SCREAMTree","SCREAM simulation tree");
  tree_->Branch("fEnergyAbs1",&fEnergyAbs1,"fEnergyAbs1/D");
  tree_->Branch("fEnergyAbs2",&fEnergyAbs2,"fEnergyAbs2/D");
  tree_->Branch("fEnergyAbs3",&fEnergyAbs3,"fEnergyAbs3/D");
  tree_->Branch("fEnergyAbs4",&fEnergyAbs4,"fEnergyAbs4/D");
  tree_->Branch("fEnergyAbs5",&fEnergyAbs5,"fEnergyAbs5/D");
  tree_->Branch("fEnergyAbs6",&fEnergyAbs6,"fEnergyAbs6/D");
  tree_->Branch("fEnergyAbs7",&fEnergyAbs7,"fEnergyAbs7/D");
  tree_->Branch("eneAr1",&eneAr1,"eneAr1/D");
  tree_->Branch("eneAr2",&eneAr2,"eneAr2/D");
  tree_->Branch("eneAr3",&eneAr3,"eneAr3/D");
  tree_->Branch("eneAr4",&eneAr4,"eneAr4/D");
  tree_->Branch("eneAr5",&eneAr5,"eneAr5/D");
  tree_->Branch("eneAr6",&eneAr6,"eneAr6/D");
  tree_->Branch("eneinallcha",&eneinallcha,"eneinallcha/D");
  tree_->Branch("eneAr1sm",&eneAr1sm,"eneAr1sm/D");
  tree_->Branch("eneAr2sm",&eneAr2sm,"eneAr2sm/D");
  tree_->Branch("eneAr3sm",&eneAr3sm,"eneAr3sm/D");
  tree_->Branch("eneAr4sm",&eneAr4sm,"eneAr4sm/D");
  tree_->Branch("eneAr5sm",&eneAr5sm,"eneAr5sm/D");
  tree_->Branch("eneAr6sm",&eneAr6sm,"eneAr6sm/D");
  tree_->Branch("eneinallchasm",&eneinallcha,"eneinallchasm/D");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{ 
  outF_->cd();
  tree_->Write();
  outF_->Close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void EventAction::BeginOfEventAction(const G4Event* evt)
{  
  G4int evtNb = evt->GetEventID();
  if (evtNb%fPrintModulo == 0){ 
    G4cout << "\n---> Begin of event   : " << evtNb << G4endl;
  }
 counter=evtNb+1;
 // initialisation per event
 eneAr1 = 0.; eneAr2 = 0.; eneAr3 = 0.; eneAr4 = 0.; eneAr5 = 0.; eneAr6 = 0.; eneinallcha = 0.;
 eneAr1sm = 0.; eneAr2sm = 0.; eneAr3sm = 0.; eneAr4sm = 0.; eneAr5sm = 0.; eneAr6sm = 0.;eneinallchasm = 0.;
 //if (evtNb==1){

  
 for(G4int j=0;j<10;j++){
   for(G4int i=0;i<10;i++){

     fEnergyAbs1 = 0.; fEnergyAbs2 = 0.; fEnergyAbs3 = 0.; fEnergyAbs4 = 0.; fEnergyAbs5 = 0.; fEnergyAbs6 = 0.;
     fEnergyAl1 = 0.; fEnergyAl2 = 0.; fEnergyAl3 = 0.; fEnergyAl4 = 0.; fEnergyAl5 = 0.; fEnergyAl6 = 0.; 
     fEnergyPCB1 = 0.; fEnergyPCB2 = 0.; fEnergyPCB3 = 0.; fEnergyPCB4 = 0.; fEnergyPCB5 = 0.; fEnergyPCB6 = 0.;
     fEnergyAr1[i][j]=0.; fEnergyAr2[i][j]=0.; fEnergyAr3[i][j]=0.; fEnergyAr4[i][j]=0.; fEnergyAr5[i][j]=0.; fEnergyAr6[i][j]=0.; 
     AdcValsAr1[i][j]=0; AdcValsAr2[i][j]=0; AdcValsAr3[i][j]=0; AdcValsAr4[i][j]=0; AdcValsAr5[i][j]=0; AdcValsAr6[i][j]=0; 

   }
 }


}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//std::ofstream outfile("/home/Kalamaris/geant4/Test-build");

void EventAction::EndOfEventAction(const G4Event* evt)
{
  //accumulates statistic
  G4int evtNb = evt->GetEventID();
  G4int  Mapping[96]={89,91,93,95,49,51,53,55,
                      70,72,74,76,78,80,81,83,85,87,
                      77,84,79,82,42,44,46,48,66,68,
                      67,94,69,92,71,90,73,88,75,86,
                      41,56,43,54,45,52,47,50,65,96,
                      18,15,64,33,62,35,60,37,58,39,
                      28,5,26,07,24,9,22,11,20,13,
                      14,16,34,36,38,40,32,1,30,3,
                      25,27,29,31,2,4,6,8,10,12,
                      57,59,61,63,17,19,21,23};
  

  fRunAct->fillPerEvent(fEnergyAbs1, fEnergyAbs2, fEnergyAbs3, fEnergyAbs4, fEnergyAbs5, fEnergyAbs6, fEnergyAbs7, fEnergyAl1, fEnergyAl2, fEnergyAl3, fEnergyAl4, fEnergyAl5, fEnergyAl6, fEnergyPCB1, fEnergyPCB2, fEnergyPCB3, fEnergyPCB4, fEnergyPCB5, fEnergyPCB6, fEnergyAr1[0][0], fEnergyAr2[0][0], fEnergyAr3[0][0], fEnergyAr4[0][0], fEnergyAr5[0][0], fEnergyAr6[0][0]);
  
  G4double Ekin=fKin->GetParticleGun()->GetParticleEnergy();
  G4double mass=fKin->GetParticleGun()->GetParticleDefinition()->GetPDGMass();

  G4int padcounter=0;

  //if (evtNb==0) {outfile("Run.data",std::ios::out);}
  // outfile.open("Run.data", std::ofstream::out | std::ofstream::app);
  outfile.open("Run.data", std::ios::out | std::ios::app);
  outfile<<"START_RUN"<<G4endl;
  outfile<<evtNb<<G4endl;
  outfile<<evtNb<<G4endl;
  outfile<<"START_EVENT"<<G4endl;
  outfile<<evtNb+1<<G4endl;
  outfile<<evtNb+1<<G4endl;
  
  outfile.close();

  //Smearing
  TRandom3 lrndm;
  //If the seed is zero the seed is set to a random value
  lrndm.SetSeed(0);
  
  //Loop over the first Argon chamber and the 10x10 pads
  for(G4int j=0;j<10;j++){
    for(G4int i=0;i<10;i++){

      // Chamber1
      if ( fEnergyAr1[i][j] != 0 ){
	//std::cout<< lrndm.Gaus( 1., (0.237/sqrt(fEnergyGap[k][i][j] * 1000.) ) + 0.08  ) << std::endl;
	eneAr1 = eneAr1 + fEnergyAr1[i][j]; 
	eneAr1sm = eneAr1sm + (fEnergyAr1[i][j] *  lrndm.Gaus( 1., (0.237/sqrt(fEnergyAr1[i][j] * 1000.) ) + 0.08  ) );
      }

      outfile.open("Run.data", std::ofstream::out | std::ofstream::app);

      if ((i==0) && (j==0)){outfile<<96<<G4endl;}

      // AdcVals[k][i][j]=G4int((fEnergyGap[k][i][j]/0.0012)*150.+200.);
      AdcValsAr1[i][j]=G4int( (fEnergyAr1[i][j] * 30.45158) + 200.   ); //Conversion to from keV to adc hardcoded. Fix Me! This is for star. 

      ovf=false;
      if (AdcValsAr1[i][j]>1023) {ovf=true;AdcValsAr1[i][j]=1023;ovfcounter++;}
      else {ovf=false;}
  
      if (  ! (((i==0) && (j==0)) || ((i==0) && (j==9)) || ((i==9) && (j==0)) || ((i==9) && (j==9)) )   ){ 
	outfile/*<<"Adc Value ["<<k<<","<<i<<","<<j<<"]="*/<<(10000000*1 + 1000000*G4int(ovf) + 10000*(Mapping[padcounter]-1) + AdcValsAr1[i][j])<<G4endl;
	padcounter++;   
 
	if (padcounter>95) padcounter=0;
 
      }
      
      outfile.close();

      //fill histograms
      // fHistoManager->FillHisto(1, fEnergyAbs);
      // fHistoManager->FillHisto(2, fEnergyGap);

    } // end of loop over pad i
  } // end of loop over pad j
  


  //======================================================
  //Loop over the second Argon chamber and the 10x10 pads
  for(G4int j=0;j<10;j++){
    for(G4int i=0;i<10;i++){

      // Chamber2
      if ( fEnergyAr2[i][j] != 0 ){
	  eneAr2 = eneAr2 + fEnergyAr2[i][j]; 
	  eneAr2sm = eneAr2sm + (fEnergyAr2[i][j] *  lrndm.Gaus( 1., (0.237/sqrt(fEnergyAr2[i][j] * 1000.) ) + 0.08   ) );
      }

      outfile.open("Run.data", std::ofstream::out | std::ofstream::app);

      if ((i==0) && (j==0)){outfile<<96<<G4endl;}

      // AdcVals[k][i][j]=G4int((fEnergyGap[k][i][j]/0.0012)*150.+200.);
      AdcValsAr2[i][j]=G4int( (fEnergyAr2[i][j] * 30.45158) + 200.   ); //Conversion to from keV to adc hardcoded. Fix Me! This is for star. 

      ovf=false;
      //if ((i==0) && (j==0)){outfile<<96<<G4endl;}
      if (AdcValsAr2[i][j]>1023) {ovf=true;AdcValsAr2[i][j]=1023;ovfcounter++;}
      else {ovf=false;}
  
      if (  ! (((i==0) && (j==0)) || ((i==0) && (j==9)) || ((i==9) && (j==0)) || ((i==9) && (j==9)) )   ){ 
	outfile/*<<"Adc Value ["<<k<<","<<i<<","<<j<<"]="*/<<(10000000*1 + 1000000*G4int(ovf) + 10000*(Mapping[padcounter]-1) + AdcValsAr2[i][j])<<G4endl;
	padcounter++;   
 
	if (padcounter>95) padcounter=0;
 
      }
      
      outfile.close();

      //fill histograms
      // fHistoManager->FillHisto(1, fEnergyAbs);
      // fHistoManager->FillHisto(2, fEnergyGap);

    } // end of loop over pad i
  } // end of loop over pad j

  //======================================================
  //Loop over the third Argon chamber and the 10x10 pads
  for(G4int j=0;j<10;j++){
    for(G4int i=0;i<10;i++){

      // Chamber3
      if ( fEnergyAr3[i][j] != 0 ){
	  eneAr3 = eneAr3 + fEnergyAr3[i][j]; 
	  eneAr3sm = eneAr3sm + (fEnergyAr3[i][j] *  lrndm.Gaus( 1., (0.237/sqrt(fEnergyAr3[i][j] * 1000.) ) + 0.08   ) );
      }

      outfile.open("Run.data", std::ofstream::out | std::ofstream::app);

      if ((i==0) && (j==0)){outfile<<96<<G4endl;}

      // AdcVals[k][i][j]=G4int((fEnergyGap[k][i][j]/0.0012)*150.+200.);
      AdcValsAr3[i][j]=G4int( (fEnergyAr3[i][j] * 30.45158) + 200.   ); //Conversion to from keV to adc hardcoded. Fix Me! This is for star. 

      ovf=false;
      //if ((i==0) && (j==0)){outfile<<96<<G4endl;}
      if (AdcValsAr3[i][j]>1023) {ovf=true;AdcValsAr3[i][j]=1023;ovfcounter++;}
      else {ovf=false;}
  
      if (  ! (((i==0) && (j==0)) || ((i==0) && (j==9)) || ((i==9) && (j==0)) || ((i==9) && (j==9)) )   ){ 
	outfile/*<<"Adc Value ["<<k<<","<<i<<","<<j<<"]="*/<<(10000000*1 + 1000000*G4int(ovf) + 10000*(Mapping[padcounter]-1) + AdcValsAr3[i][j])<<G4endl;
	padcounter++;   
 
	if (padcounter>95) padcounter=0;
 
      }
      
      outfile.close();

      //fill histograms
      // fHistoManager->FillHisto(1, fEnergyAbs);
      // fHistoManager->FillHisto(2, fEnergyGap);

    } // end of loop over pad i
  } // end of loop over pad j
  
  //======================================================
  //Loop over the fourth Argon chamber and the 10x10 pads
  for(G4int j=0;j<10;j++){
    for(G4int i=0;i<10;i++){

      // Chamber4
      if ( fEnergyAr4[i][j] != 0 ){
	  eneAr4 = eneAr4 + fEnergyAr4[i][j]; 
	  eneAr4sm = eneAr4sm + (fEnergyAr4[i][j] *  lrndm.Gaus( 1., (0.237/sqrt(fEnergyAr4[i][j] * 1000.) ) + 0.08   ) );
      }

      outfile.open("Run.data", std::ofstream::out | std::ofstream::app);

      if ((i==0) && (j==0)){outfile<<96<<G4endl;}

      // AdcVals[k][i][j]=G4int((fEnergyGap[k][i][j]/0.0012)*150.+200.);
      AdcValsAr4[i][j]=G4int( (fEnergyAr4[i][j] * 30.45158) + 200.   ); //Conversion to from keV to adc hardcoded. Fix Me! This is for star. 

      ovf=false;
      //if ((i==0) && (j==0)){outfile<<96<<G4endl;}
      if (AdcValsAr4[i][j]>1023) {ovf=true;AdcValsAr4[i][j]=1023;ovfcounter++;}
      else {ovf=false;}
  
      if (  ! (((i==0) && (j==0)) || ((i==0) && (j==9)) || ((i==9) && (j==0)) || ((i==9) && (j==9)) )   ){ 
	outfile/*<<"Adc Value ["<<k<<","<<i<<","<<j<<"]="*/<<(10000000*1 + 1000000*G4int(ovf) + 10000*(Mapping[padcounter]-1) + AdcValsAr4[i][j])<<G4endl;
	padcounter++;   
 
	if (padcounter>95) padcounter=0;
 
      }
      
      outfile.close();

      //fill histograms
      // fHistoManager->FillHisto(1, fEnergyAbs);
      // fHistoManager->FillHisto(2, fEnergyGap);

    } // end of loop over pad i
  } // end of loop over pad j
  
  //======================================================
  //Loop over the fifth Argon chamber and the 10x10 pads
  for(G4int j=0;j<10;j++){
    for(G4int i=0;i<10;i++){

      // Chamber5
      if ( fEnergyAr5[i][j] != 0 ){
	  eneAr5 = eneAr5 + fEnergyAr5[i][j]; 
	  eneAr5sm = eneAr5sm + (fEnergyAr5[i][j] *  lrndm.Gaus( 1., (0.237/sqrt(fEnergyAr5[i][j] * 1000.) ) + 0.08   ) );
      }

      outfile.open("Run.data", std::ofstream::out | std::ofstream::app);

      if ((i==0) && (j==0)){outfile<<96<<G4endl;}

      // AdcVals[k][i][j]=G4int((fEnergyGap[k][i][j]/0.0012)*150.+200.);
      AdcValsAr5[i][j]=G4int( (fEnergyAr5[i][j] * 30.45158) + 200.   ); //Conversion to from keV to adc hardcoded. Fix Me! This is for star. 

      ovf=false;
      //if ((i==0) && (j==0)){outfile<<96<<G4endl;}
      if (AdcValsAr5[i][j]>1023) {ovf=true;AdcValsAr5[i][j]=1023;ovfcounter++;}
      else {ovf=false;}
  
      if (  ! (((i==0) && (j==0)) || ((i==0) && (j==9)) || ((i==9) && (j==0)) || ((i==9) && (j==9)) )   ){ 
	outfile/*<<"Adc Value ["<<k<<","<<i<<","<<j<<"]="*/<<(10000000*1 + 1000000*G4int(ovf) + 10000*(Mapping[padcounter]-1) + AdcValsAr5[i][j])<<G4endl;
	padcounter++;   
 
	if (padcounter>95) padcounter=0;
 
      }
      
      outfile.close();

      //fill histograms
      // fHistoManager->FillHisto(1, fEnergyAbs);
      // fHistoManager->FillHisto(2, fEnergyGap);

    } // end of loop over pad i
  } // end of loop over pad j
   
  //======================================================
  //Loop over the sixth Argon chamber and the 10x10 pads
  for(G4int j=0;j<10;j++){
    for(G4int i=0;i<10;i++){

      // Chamber6
      if ( fEnergyAr6[i][j] != 0 ){
	  eneAr6 = eneAr6 + fEnergyAr6[i][j]; 
	  eneAr6sm = eneAr6sm + (fEnergyAr6[i][j] *  lrndm.Gaus( 1., (0.237/sqrt(fEnergyAr6[i][j] * 1000.) ) + 0.08   ) );
      }

      outfile.open("Run.data", std::ofstream::out | std::ofstream::app);

      if ((i==0) && (j==0)){outfile<<96<<G4endl;}

      // AdcVals[k][i][j]=G4int((fEnergyGap[k][i][j]/0.0012)*150.+200.);
      AdcValsAr6[i][j]=G4int( (fEnergyAr6[i][j] * 30.45158) + 200.   ); //Conversion to from keV to adc hardcoded. Fix Me! This is for star. 

      ovf=false;
      //if ((i==0) && (j==0)){outfile<<96<<G4endl;}
      if (AdcValsAr6[i][j]>1023) {ovf=true;AdcValsAr6[i][j]=1023;ovfcounter++;}
      else {ovf=false;}
  
      if (  ! (((i==0) && (j==0)) || ((i==0) && (j==9)) || ((i==9) && (j==0)) || ((i==9) && (j==9)) )   ){ 
	outfile/*<<"Adc Value ["<<k<<","<<i<<","<<j<<"]="*/<<(10000000*1 + 1000000*G4int(ovf) + 10000*(Mapping[padcounter]-1) + AdcValsAr6[i][j])<<G4endl;
	padcounter++;   
 
	if (padcounter>95) padcounter=0;
 
      }
      
      outfile.close();

      //fill histograms
      // fHistoManager->FillHisto(1, fEnergyAbs);
      // fHistoManager->FillHisto(2, fEnergyGap);

    } // end of loop over pad i
  } // end of loop over pad j
 
  
  
  outfile.open("Run.data", std::ios::out | std::ios::app);
  outfile<<0<<G4endl;
  outfile<<0<<G4endl;
  outfile<<0<<G4endl;
  outfile<<0<<G4endl;
  
  outfile.close();

  
  fHistoManager->FillHisto(1, fEnergyAr1[0][4]);//First Detector pad (1,5) 
  fHistoManager->FillHisto(2, fEnergyAr1[1][4]);//First Detector pad (2,5)
  fHistoManager->FillHisto(3, fEnergyAr1[2][4]);//First Detector pad (3,5)
  fHistoManager->FillHisto(4, fEnergyAr1[3][4]);//First Detector pad (4,5)
  fHistoManager->FillHisto(5, fEnergyAr1[4][4]);//First Detector pad (5,5)
  fHistoManager->FillHisto(6, eneAr1);//Energy in chamber1
  fHistoManager->FillHisto(7, eneAr2);//Energy in chamber2
  fHistoManager->FillHisto(8, eneAr3);//Energy in chamber3
  fHistoManager->FillHisto(9, eneAr4);//Energy in chamber4
  fHistoManager->FillHisto(10, eneAr5);//Energy in chamber5
  fHistoManager->FillHisto(11, eneAr6);//Energy in chamber6

  eneinallcha = eneAr1 + eneAr2+ eneAr3+ eneAr4+ eneAr5+ eneAr6;//Energy in all chamber
  fHistoManager->FillHisto(12, eneinallcha);//Energy in all chambers

  //Smeared distributions
  if ( fEnergyAr1[0][4] != 0 ){
    fHistoManager->FillHisto(13, fEnergyAr1[0][4] * lrndm.Gaus( 1., (0.237/sqrt(fEnergyAr1[0][4] * 1000.) )  + 0.08  )  );//First Detector pad (1,5) 
  }
  if ( fEnergyAr1[1][4] != 0 ){
    fHistoManager->FillHisto(14, fEnergyAr1[1][4] * lrndm.Gaus( 1., (0.237/sqrt(fEnergyAr1[1][4] * 1000.) )  + 0.08  )  );//First Detector pad (2,5)
  }
  if ( fEnergyAr1[2][4] != 0 ){
    fHistoManager->FillHisto(15, fEnergyAr1[2][4] * lrndm.Gaus( 1., (0.237/sqrt(fEnergyAr1[2][4] * 1000.) )  + 0.08 )  );//First Detector pad (3,5)
  }
  if ( fEnergyAr1[3][4] != 0 ){
    fHistoManager->FillHisto(16, fEnergyAr1[3][4] * lrndm.Gaus( 1., (0.237/sqrt(fEnergyAr1[3][4] * 1000.) )  + 0.08 )  );//First Detector pad (4,5)
  }
  if ( fEnergyAr1[4][4] != 0 ){
    fHistoManager->FillHisto(17, fEnergyAr1[4][4] * lrndm.Gaus( 1., (0.237/sqrt(fEnergyAr1[4][4] * 1000.) )  + 0.08 )  );//First Detector pad (5,5)
  }
  fHistoManager->FillHisto(18, eneAr1sm );//Smear Energy in chamber1
  fHistoManager->FillHisto(19, eneAr2sm );//Smear Energy in chamber2
  fHistoManager->FillHisto(20, eneAr3sm);//Smear Energy in chamber3
  fHistoManager->FillHisto(21, eneAr4sm);//Smear Energy in chamber4
  fHistoManager->FillHisto(22, eneAr5sm);//Smear Energy in chamber5
  fHistoManager->FillHisto(23, eneAr6sm);//Smear Energy in chamber6

  eneinallchasm = eneAr1sm + eneAr2sm+ eneAr3sm+ eneAr4sm+ eneAr5sm+ eneAr6sm;//Energy in all chamber

  fHistoManager->FillHisto(24, eneinallchasm);//Energy in all chamber
  if ( fEnergyAr1[0][4] != 0 ){
    fHistoManager->FillHisto(25, lrndm.Gaus( 1., (0.237/sqrt(fEnergyAr1[0][4] * 1000.) ) + 0.08  ) );//Testing the smear factor only for the 1st chamber pad. This should be gaussian.
  }
  fHistoManager->FillHisto(26, eneAr1sm * 30.95062);//Energy in chamber1 in adc
  fHistoManager->FillHisto(27, eneAr2sm * 30.95062);//Energy in chamber2 in adc
  fHistoManager->FillHisto(28, eneAr3sm * 30.95062);//Energy in chamber3 in adc
  fHistoManager->FillHisto(29, eneAr4sm * 30.95062);//Energy in chamber4 in adc
  fHistoManager->FillHisto(30, eneAr5sm * 30.95062);//Energy in chamber5 in adc
  fHistoManager->FillHisto(31, eneAr6sm * 30.95062);//Energy in chamber6 in adc

  fHistoManager->FillHisto(32, eneinallcha * 30.95062);//Energy in all chamber in adc

  fHistoManager->FillHisto(33, fEnergyAbs1);//Energy deposited in absorber1
  fHistoManager->FillHisto(34, fEnergyAbs2);//Energy deposited in absorber2
  fHistoManager->FillHisto(35, fEnergyAbs3);//Energy deposited in absorber3
  fHistoManager->FillHisto(36, fEnergyAbs4);//Energy deposited in absorber4
  fHistoManager->FillHisto(37, fEnergyAbs5);//Energy deposited in absorber5
  fHistoManager->FillHisto(38, fEnergyAbs6);//Energy deposited in absorber6
  fHistoManager->FillHisto(39, fEnergyAbs7);//Energy deposited in absorber7

  fHistoManager->FillHisto(40, 100.*fEnergyAbs1/(Ekin+mass));//Energy deposited in absorber1 (percent of Eincident)
  fHistoManager->FillHisto(41, 100.*fEnergyAbs2/(Ekin+mass));//Energy deposited in absorber2 (percent of Eincident)
  fHistoManager->FillHisto(42, 100.*fEnergyAbs3/(Ekin+mass));//Energy deposited in absorber3 (percent of Eincident)
  fHistoManager->FillHisto(43, 100.*fEnergyAbs4/(Ekin+mass));//Energy deposited in absorber4 (percent of Eincident)
  fHistoManager->FillHisto(44, 100.*fEnergyAbs5/(Ekin+mass));//Energy deposited in absorber5 (percent of Eincident)
  fHistoManager->FillHisto(45, 100.*fEnergyAbs6/(Ekin+mass));//Energy deposited in absorber6 (percent of Eincident)
  fHistoManager->FillHisto(46, 100.*fEnergyAbs7/(Ekin+mass));//Energy deposited in absorber7 (percent of Eincident)

  fHistoManager->FillHisto(47, 100.*eneAr1/(Ekin+mass));//Energy in chamber1 (percent of Eincident)
  fHistoManager->FillHisto(48, 100.*eneAr2/(Ekin+mass));//Energy in chamber2 (percent of Eincident)
  fHistoManager->FillHisto(49, 100.*eneAr3/(Ekin+mass));//Energy in chamber3 (percent of Eincident)
  fHistoManager->FillHisto(50, 100.*eneAr4/(Ekin+mass));//Energy in chamber4 (percent of Eincident)
  fHistoManager->FillHisto(51, 100.*eneAr5/(Ekin+mass));//Energy in chamber5 (percent of Eincident)
  fHistoManager->FillHisto(52, 100.*eneAr6/(Ekin+mass));//Energy in chamber6 (percent of Eincident)

  tree_->Fill();
  //G4cout<< "First Detector pad (1,5) " << fEnergyAr1[0][4] <<G4endl;
  //G4cout<< "First Detector pad (2,5) " << fEnergyAr1[1][4] <<G4endl;

  //G4cout<<"Padcounter  "<<padcounter<<G4endl;
  //G4cout<<"End of event # "<<counter<<G4endl;
}

  /* for(G4int j=0;j<10;j++){
  for(G4int i=0;i<10;i++){

  if (AdcVals[i][j]=1023) {ovfcounter++;}

}


}
}
  G4cout<<"Overflow # per event "<<ovfcounter<<G4endl;
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
