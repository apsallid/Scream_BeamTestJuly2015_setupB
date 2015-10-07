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
/// \file analysis/shared/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
// $Id: DetectorConstruction.cc 77256 2013-11-22 10:10:23Z gcosmo $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"

#include "G4SDManager.hh"
//#include "G4VSensitiveDetector.hh"
//#include "SDGap.hh"
#include "G4PVReplica.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 :G4VUserDetectorConstruction(),
  fAbsorber1Material(0),
  fAbsorber2Material(0),
  fAbsorber3Material(0),
  fAbsorber4Material(0),
  fAbsorber5Material(0),
  fAbsorber6Material(0),
  fAbsorber7Material(0),
  fAl1Material(0),
  fAl2Material(0),
  fAl3Material(0),
  fAl4Material(0),
  fAl5Material(0),
  fAl6Material(0),
  fAr1Material(0),
  fAr2Material(0),
  fAr3Material(0),
  fAr4Material(0),
  fAr5Material(0),
  fAr6Material(0),
  fPCB1Material(0),
  fPCB2Material(0),
  fPCB3Material(0),
  fPCB4Material(0),
  fPCB5Material(0),
  fPCB6Material(0),
  fDefaultMaterial(0),
  fSolidWorld(0),fLogicWorld(0),fPhysiWorld(0),
  fSolidCalor(0),fLogicCalor(0),fPhysiCalor(0),
  fSolidAbsorber1(0),fLogicAbsorber1(0),fPhysiAbsorber1(0),
  fSolidAbsorber2(0),fLogicAbsorber2(0),fPhysiAbsorber2(0),
  fSolidAbsorber3(0),fLogicAbsorber3(0),fPhysiAbsorber3(0),
  fSolidAbsorber4(0),fLogicAbsorber4(0),fPhysiAbsorber4(0),
  fSolidAbsorber5(0),fLogicAbsorber5(0),fPhysiAbsorber5(0),
  fSolidAbsorber6(0),fLogicAbsorber6(0),fPhysiAbsorber6(0),
  fSolidAbsorber7(0),fLogicAbsorber7(0),fPhysiAbsorber7(0),
  fSolidAl1(0),fLogicAl1(0),fPhysiAl1(0),
  fSolidAl2(0),fLogicAl2(0),fPhysiAl2(0),
  fSolidAl3(0),fLogicAl3(0),fPhysiAl3(0),
  fSolidAl4(0),fLogicAl4(0),fPhysiAl4(0),
  fSolidAl5(0),fLogicAl5(0),fPhysiAl5(0),
  fSolidAl6(0),fLogicAl6(0),fPhysiAl6(0),
  fLogicAr1(0),fLogicAr2(0),fLogicAr3(0),
  fLogicAr4(0),fLogicAr5(0),fLogicAr6(0),
  //fSolidAr1(0),fLogicAr1(0),fPhysiAr1(0),
  //fSolidAr2(0),fLogicAr2(0),fPhysiAr2(0),
  fSolidPCB1(0),fLogicPCB1(0),fPhysiPCB1(0),
  fSolidPCB2(0),fLogicPCB2(0),fPhysiPCB2(0),
  fSolidPCB3(0),fLogicPCB3(0),fPhysiPCB3(0),
  fSolidPCB4(0),fLogicPCB4(0),fPhysiPCB4(0),
  fSolidPCB5(0),fLogicPCB5(0),fPhysiPCB5(0),
  fSolidPCB6(0),fLogicPCB6(0),fPhysiPCB6(0),
  fDetectorMessenger(0)
{
  // default parameter values of the calorimeter
  fAbsorber1Thickness = 4.0*cm; 
  fAbsorber2Thickness = 6.0*cm; 
  fAbsorber3Thickness = 6.0*cm; 
  fAbsorber4Thickness = 6.0*cm; 
  fAbsorber5Thickness = 6.0*cm; 
  fAbsorber6Thickness = 6.0*cm; 
  fAbsorber7Thickness = 6.0*cm; 
  fAl1Thickness      = 2.*mm;
  fAl2Thickness      = 2.*mm;
  fAl3Thickness      = 2.*mm;
  fAl4Thickness      = 2.*mm;
  fAl5Thickness      = 2.*mm;
  fAl6Thickness      = 2.*mm;
  fAr1Thickness      = 3.*mm;
  fAr2Thickness      = 3.*mm;
  fAr3Thickness      = 3.*mm;
  fAr4Thickness      = 3.*mm;
  fAr5Thickness      = 3.*mm;
  fAr6Thickness      = 3.*mm;
  fPCB1Thickness     = 2.*mm;
  fPCB2Thickness     = 2.*mm;
  fPCB3Thickness     = 2.*mm;
  fPCB4Thickness     = 2.*mm;
  fPCB5Thickness     = 2.*mm;
  fPCB6Thickness     = 2.*mm;
  fAirThickness      = 3.*mm;

  fEnvThickness      = 5.*cm;
  fCalorSizeYZ       = 10.*cm;
   
  // materials
  DefineMaterials();
  SetAbsorber1Material("G4_Fe");
  SetAbsorber2Material("G4_Fe");
  SetAbsorber3Material("G4_Fe");
  SetAbsorber4Material("G4_Fe");
  SetAbsorber5Material("G4_Fe");
  SetAbsorber6Material("G4_Fe");
  SetAbsorber7Material("G4_Fe");
  SetAl1Material("G4_Al");
  SetAl2Material("G4_Al");
  SetAl3Material("G4_Al");
  SetAl4Material("G4_Al");
  SetAl5Material("G4_Al");
  SetAl6Material("G4_Al");
  SetAr1Material("G4_Ar");
  SetAr2Material("G4_Ar");
  SetAr3Material("G4_Ar");
  SetAr4Material("G4_Ar");
  SetAr5Material("G4_Ar");
  SetAr6Material("G4_Ar");
  SetPCB1Material();
  SetPCB2Material();
  SetPCB3Material();
  SetPCB4Material();
  SetPCB5Material();
  SetPCB6Material();
  
  // create commands for interactive definition of the calorimeter
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructCalorimeter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 
// use G4-NIST materials data base
//
G4NistManager* man = G4NistManager::Instance();
fDefaultMaterial = man->FindOrBuildMaterial("G4_AIR");
man->FindOrBuildMaterial("G4_Fe");
man->FindOrBuildMaterial("G4_Ar");

// print table
//
G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructCalorimeter()
{
  G4double posx,posy,posz;
  // Clean old geometry, if any
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // complete the Calor parameters definition
  ComputeCalorParameters();
  
  // World
  fSolidWorld = new G4Box("World",                                //its name
                   fWorldSizeX/2,fWorldSizeYZ/2,fWorldSizeYZ/2);  //its size
                         
  fLogicWorld = new G4LogicalVolume(fSolidWorld,            //its solid
                                   fDefaultMaterial,        //its material
                                   "World");                //its name
                                   
  fPhysiWorld = new G4PVPlacement(0,                        //no rotation
                                 G4ThreeVector(),           //at (0,0,0)
                                 fLogicWorld,             //its logical volume
                                 "World",                   //its name
                                 0,                         //its mother  volume
                                 false,                  //no boolean operation
                                 0);                        //copy number
  
  // Calorimeter
  fSolidCalor=0; fLogicCalor=0; fPhysiCalor=0;
  // fSolidLayer=0; fLogicLayer=0; fPhysiLayer=0;
  
  if (fCalorThickness > 0.)  
    { fSolidCalor = new G4Box("Calorimeter",                //its name
			      4*fCalorThickness/2,4*fCalorSizeYZ,4*fCalorSizeYZ);//size
                                 
      fLogicCalor = new G4LogicalVolume(fSolidCalor,        //its solid
                                        fDefaultMaterial,   //its material
                                        "Calorimeter");     //its name
                                           
      fPhysiCalor = new G4PVPlacement(0,                    //no rotation
				      G4ThreeVector(),       //at (0,0,0)
				      fLogicCalor,           //its logical volume
				      "Calorimeter",         //its name
				      fLogicWorld,           //its mother  volume
				      false,              //no boolean operation
				      0);                    //copy number
  
    }
  // Absorber1
  fSolidAbsorber1=0; fLogicAbsorber1=0; fPhysiAbsorber1=0;  
  
  if (fAbsorber1Thickness > 0.) 
    { fSolidAbsorber1 = new G4Box("Absorber1",                //its name
				  fAbsorber1Thickness/2,fCalorSizeYZ/2,fCalorSizeYZ/2); 
                          
      fLogicAbsorber1 = new G4LogicalVolume(fSolidAbsorber1,    //its solid
                                            fAbsorber1Material, //its material
                                            fAbsorber1Material->GetName());//name
                                                
      fPhysiAbsorber1 = new G4PVPlacement(0,                   //no rotation
					  G4ThreeVector(fAbsorber1Thickness/2,0.,0.),  //its position
					  fLogicAbsorber1,     //its logical volume
					  fAbsorber1Material->GetName(), //its name
					  fLogicCalor,          //its mother
					  false,               //no boulean operat
					  0);                   //copy number
                                        
    }
  
  // Al1
  fSolidAl1 = 0; fLogicAl1 = 0; fPhysiAl1 = 0;
  if (fAl1Thickness > 0.){
    fSolidAl1 = new G4Box("Al1",
			  fAl1Thickness/2.,fCalorSizeYZ/2,fCalorSizeYZ/2);
                               
    fLogicAl1 = new G4LogicalVolume(fSolidAl1,
				    fAl1Material,
				    fAl1Material->GetName());
   
     
    posx = fAbsorber1Thickness + (fAl1Thickness/2.);

    fPhysiAl1= new G4PVPlacement(0,  //no rotation
				 G4ThreeVector(posx,0.,0.),       //at (0,0,0)
				 fLogicAl1,           //its logical volume
				 "Al1",               //its name
				 fLogicCalor,           //its mother  volume
				 false,               //no boolean operation
				 0);                    //copy number 
    
    
    
    
    PrintCalorParameters();     
  
    //                                        
    // Visualization attributes
    //
    fLogicWorld->SetVisAttributes (G4VisAttributes::Invisible);
    fLogicCalor->SetVisAttributes(G4VisAttributes::Invisible);
   
    // G4VisAttributes* magenta=new G4VisAttributes(true,G4Colour(1.,0.,1.));
    G4VisAttributes* blue=new G4VisAttributes(true,G4Colour(0.,0.,1.));
    G4VisAttributes* yellow=new G4VisAttributes(true,G4Colour(1.,1.,0.));
    fLogicAl1->SetVisAttributes(blue);
 
    fLogicAbsorber1->SetVisAttributes(yellow);



    G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    simpleBoxVisAtt->SetVisibility(false);
    fLogicCalor->SetVisAttributes(simpleBoxVisAtt);

  }

   //---------------------------------------------Ar1--------------------------------------------------------
  for(G4int j=0;j<10;j++){
    for(G4int i=0;i<10;i++)
      {
       

	fSolidAr1[i][j] = new G4Box("Ar1",
				    fAr1Thickness/2,fCalorSizeYZ/20,fCalorSizeYZ/20);
                               
	fLogicAr1 = new G4LogicalVolume(fSolidAr1[i][j],
					fAr1Material,
					fAr1Material->GetName());
   
     
	posx = fAbsorber1Thickness + fAl1Thickness + (fAr1Thickness/2.);
	 
	posy=-45+j*10;
	posz=-45+i*10;

	fPhysiAr1[i][j] = new G4PVPlacement(0,                  //no rotation
					    G4ThreeVector(posx,posy,posz),       //at (0,0,0)
					    fLogicAr1,           //its logical volume
					    "Ar1",               //its name
					    fLogicCalor,           //its mother  volume
					    false,               //no boolean operation
					    i);                    //copy number 



    
	PrintCalorParameters();     
  
	//                                        
	// Visualization attributes
	//
	fLogicWorld->SetVisAttributes (G4VisAttributes::Invisible);
	fLogicCalor->SetVisAttributes(G4VisAttributes::Invisible);
   
	G4VisAttributes* magenta=new G4VisAttributes(true,G4Colour(1.,0.,1.));
	G4VisAttributes* yellow=new G4VisAttributes(true,G4Colour(1.,1.,0.));
	fLogicAr1->SetVisAttributes(magenta);
	fLogicAbsorber1->SetVisAttributes(yellow);

	G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
	simpleBoxVisAtt->SetVisibility(false);
	fLogicCalor->SetVisAttributes(simpleBoxVisAtt);
      }
  }


   // PCB1
  fSolidPCB1 = 0; fLogicPCB1 = 0; fPhysiPCB1 = 0;
  if (fPCB1Thickness > 0.){
    fSolidPCB1 = new G4Box("PCB1",
			   fPCB1Thickness/2.,fCalorSizeYZ/2,fCalorSizeYZ/2);
                               
    fLogicPCB1 = new G4LogicalVolume(fSolidPCB1,
				     fPCB1Material,
				     fPCB1Material->GetName());
   
     
    posx = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + (fPCB1Thickness/2.);
    
    fPhysiPCB1= new G4PVPlacement(0,  //no rotation
				  G4ThreeVector(posx,0.,0.),       //at (0,0,0)
				  fLogicPCB1,           //its logical volume
				  "PCB1",               //its name
				  fLogicCalor,           //its mother  volume
				  false,               //no boolean operation
				  0);                    //copy number 
    
    
    
    
    PrintCalorParameters();     
  
    //                                        
    // Visualization attributes
    //
    fLogicWorld->SetVisAttributes (G4VisAttributes::Invisible);
    fLogicCalor->SetVisAttributes(G4VisAttributes::Invisible);
   
    // G4VisAttributes* magenta=new G4VisAttributes(true,G4Colour(1.,0.,1.));
    G4VisAttributes* greenn=new G4VisAttributes(true,G4Colour(0.0, 1.0, 0.0));
    G4VisAttributes* yellow=new G4VisAttributes(true,G4Colour(1.,1.,0.));
    fLogicPCB1->SetVisAttributes(greenn);
    fLogicAbsorber1->SetVisAttributes(yellow);

    G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    simpleBoxVisAtt->SetVisibility(false);
    fLogicCalor->SetVisAttributes(simpleBoxVisAtt);

  }
 
  // Absorber2
  fSolidAbsorber2=0; fLogicAbsorber2=0; fPhysiAbsorber2=0;  
  
  if (fAbsorber2Thickness > 0.) 
    { fSolidAbsorber2 = new G4Box("Absorber2",                //its name
				  fAbsorber2Thickness/2,fCalorSizeYZ/2,fCalorSizeYZ/2); 
                          
      fLogicAbsorber2 = new G4LogicalVolume(fSolidAbsorber2,    //its solid
                                            fAbsorber2Material, //its material
                                            fAbsorber2Material->GetName());//name

      posx = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + (fAbsorber2Thickness/2.);


                                                
      fPhysiAbsorber2 = new G4PVPlacement(0,                   //no rotation
					  G4ThreeVector(posx,0.,0.),  //its position
					  fLogicAbsorber2,     //its logical volume
					  fAbsorber2Material->GetName(), //its name
					  fLogicCalor,          //its mother
					  false,               //no boulean operat
					  0);                   //copy number
                                        
    }
  

  // Al2
  fSolidAl2 = 0; fLogicAl2 = 0; fPhysiAl2 = 0;
  if (fAl2Thickness > 0.){
    fSolidAl2 = new G4Box("Al2",
			  fAl2Thickness/2.,fCalorSizeYZ/2,fCalorSizeYZ/2);
                               
    fLogicAl2 = new G4LogicalVolume(fSolidAl2,
				    fAl2Material,
				    fAl2Material->GetName());
   
     
    posx = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + (fAl2Thickness/2.);
    
    fPhysiAl2= new G4PVPlacement(0,  //no rotation
				 G4ThreeVector(posx,0.,0.),       //at (0,0,0)
				 fLogicAl2,           //its logical volume
				 "Al2",               //its name
				 fLogicCalor,           //its mother  volume
				 false,               //no boolean operation
				 0);                    //copy number 
    
    
    
    
    PrintCalorParameters();     
  
    //                                        
    // Visualization attributes
    //
    fLogicWorld->SetVisAttributes (G4VisAttributes::Invisible);
    fLogicCalor->SetVisAttributes(G4VisAttributes::Invisible);
   
    G4VisAttributes* blue=new G4VisAttributes(true,G4Colour(0.,0.,1.));
    G4VisAttributes* yellow=new G4VisAttributes(true,G4Colour(1.,1.,0.));
    fLogicAl2->SetVisAttributes(blue);
 
    fLogicAbsorber1->SetVisAttributes(yellow);
    fLogicAbsorber2->SetVisAttributes(yellow);



    G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    simpleBoxVisAtt->SetVisibility(false);
    fLogicCalor->SetVisAttributes(simpleBoxVisAtt);

  }

   //---------------------------------------------Ar2--------------------------------------------------------
  for(G4int j=0;j<10;j++){
    for(G4int i=0;i<10;i++)
      {
       

	fSolidAr2[i][j] = new G4Box("Ar2",
				    fAr2Thickness/2.,fCalorSizeYZ/20,fCalorSizeYZ/20);
                               
	fLogicAr2 = new G4LogicalVolume(fSolidAr2[i][j],
					fAr2Material,
					fAr2Material->GetName());
   
     
	posx = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + (fAr2Thickness/2.);
	 	 
	posy=-45+j*10;
	posz=-45+i*10;

	fPhysiAr2[i][j] = new G4PVPlacement(0,                  //no rotation
					    G4ThreeVector(posx,posy,posz),       //at (0,0,0)
					    fLogicAr2,           //its logical volume
					    "Ar2",               //its name
					    fLogicCalor,           //its mother  volume
					    false,               //no boolean operation
					    i);                    //copy number 



    
	PrintCalorParameters();     
  
	//                                        
	// Visualization attributes
	//
	fLogicWorld->SetVisAttributes (G4VisAttributes::Invisible);
	fLogicCalor->SetVisAttributes(G4VisAttributes::Invisible);
   
	G4VisAttributes* magenta=new G4VisAttributes(true,G4Colour(1.,0.,1.));
	G4VisAttributes* yellow=new G4VisAttributes(true,G4Colour(1.,1.,0.));
	fLogicAr2->SetVisAttributes(magenta);
	fLogicAbsorber1->SetVisAttributes(yellow);
	fLogicAbsorber2->SetVisAttributes(yellow);

	G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
	simpleBoxVisAtt->SetVisibility(false);
	fLogicCalor->SetVisAttributes(simpleBoxVisAtt);
      }
  }

   // PCB2
  fSolidPCB2 = 0; fLogicPCB2 = 0; fPhysiPCB2 = 0;
  if (fPCB2Thickness > 0.){
    fSolidPCB2 = new G4Box("PCB2",
			   fPCB2Thickness/2.,fCalorSizeYZ/2,fCalorSizeYZ/2);
                               
    fLogicPCB2 = new G4LogicalVolume(fSolidPCB2,
				     fPCB2Material,
				     fPCB2Material->GetName());
   
     
    posx = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + fAr2Thickness + (fPCB2Thickness/2.);
    
    fPhysiPCB2= new G4PVPlacement(0,  //no rotation
				  G4ThreeVector(posx,0.,0.),       //at (0,0,0)
				  fLogicPCB2,           //its logical volume
				  "PCB2",               //its name
				  fLogicCalor,           //its mother  volume
				  false,               //no boolean operation
				  0);                    //copy number 
    
    PrintCalorParameters();     
  
    //                                        
    // Visualization attributes
    //
    fLogicWorld->SetVisAttributes (G4VisAttributes::Invisible);
    fLogicCalor->SetVisAttributes(G4VisAttributes::Invisible);
   
    G4VisAttributes* greenn=new G4VisAttributes(true,G4Colour(0.0, 1.0, 0.0));
    G4VisAttributes* yellow=new G4VisAttributes(true,G4Colour(1.,1.,0.));
    fLogicPCB2->SetVisAttributes(greenn);
    fLogicAbsorber1->SetVisAttributes(yellow);
    fLogicAbsorber2->SetVisAttributes(yellow);
 
    G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    simpleBoxVisAtt->SetVisibility(false);
    fLogicCalor->SetVisAttributes(simpleBoxVisAtt);

  }
  
  // Absorber3
  fSolidAbsorber3=0; fLogicAbsorber3=0; fPhysiAbsorber3=0;  
  
  if (fAbsorber3Thickness > 0.) 
    { fSolidAbsorber3 = new G4Box("Absorber3",                //its name
                          fAbsorber3Thickness/2,fCalorSizeYZ/2,fCalorSizeYZ/2); 
                          
      fLogicAbsorber3 = new G4LogicalVolume(fSolidAbsorber3,    //its solid
                                            fAbsorber3Material, //its material
                                            fAbsorber3Material->GetName());//name

      posx = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + fAr2Thickness + fPCB2Thickness + fAirThickness + (fAbsorber3Thickness/2.);

                                                
      fPhysiAbsorber3 = new G4PVPlacement(0,                   //no rotation
					  G4ThreeVector(posx,0.,0.),  //its position
					  fLogicAbsorber3,     //its logical volume
					  fAbsorber3Material->GetName(), //its name
					  fLogicCalor,          //its mother
					  false,               //no boulean operat
					  0);                   //copy number
                                        
    }


  // Al3
  fSolidAl3 = 0; fLogicAl3 = 0; fPhysiAl3 = 0;
  if (fAl3Thickness > 0.){
    fSolidAl3 = new G4Box("Al3",
			 fAl3Thickness/2.,fCalorSizeYZ/2,fCalorSizeYZ/2);
                               
    fLogicAl3 = new G4LogicalVolume(fSolidAl3,
				   fAl3Material,
				   fAl3Material->GetName());
   
     
    posx = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + fAr2Thickness + fPCB2Thickness + + fAirThickness + fAbsorber3Thickness + (fAl3Thickness/2.); 
   
    fPhysiAl3= new G4PVPlacement(0,  //no rotation
				 G4ThreeVector(posx,0.,0.),       //at (0,0,0)
				 fLogicAl3,           //its logical volume
				 "Al3",               //its name
				 fLogicCalor,           //its mother  volume
				 false,               //no boolean operation
				 0);                    //copy number 
    
    
    
    
    PrintCalorParameters();     
  
    //                                        
    // Visualization attributes
    //
    fLogicWorld->SetVisAttributes (G4VisAttributes::Invisible);
    fLogicCalor->SetVisAttributes(G4VisAttributes::Invisible);
   
    G4VisAttributes* blue=new G4VisAttributes(true,G4Colour(0.,0.,1.));
    G4VisAttributes* yellow=new G4VisAttributes(true,G4Colour(1.,1.,0.));
    fLogicAl3->SetVisAttributes(blue);
 
    fLogicAbsorber1->SetVisAttributes(yellow);
    fLogicAbsorber2->SetVisAttributes(yellow);



    G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    simpleBoxVisAtt->SetVisibility(false);
    fLogicCalor->SetVisAttributes(simpleBoxVisAtt);

  }

   //---------------------------------------------Ar3--------------------------------------------------------
   for(G4int j=0;j<10;j++){
     for(G4int i=0;i<10;i++)
       {
       

	 fSolidAr3[i][j] = new G4Box("Ar3",
				      fAr3Thickness/2.,fCalorSizeYZ/20,fCalorSizeYZ/20);
                               
	 fLogicAr3 = new G4LogicalVolume(fSolidAr3[i][j],
					  fAr3Material,
					  fAr3Material->GetName());
   
     
	 posx = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + fAr2Thickness + fPCB2Thickness + + fAirThickness + fAbsorber3Thickness + fAl3Thickness + (fAr3Thickness/2.); 

	 posy=-45+j*10;
	 posz=-45+i*10;

	 fPhysiAr3[i][j] = new G4PVPlacement(0,                  //no rotation
					    G4ThreeVector(posx,posy,posz),       //at (0,0,0)
					    fLogicAr3,           //its logical volume
					    "Ar3",               //its name
					    fLogicCalor,           //its mother  volume
					    false,               //no boolean operation
					    i);                    //copy number 



    
	 PrintCalorParameters();     
  
	 //                                        
	 // Visualization attributes
	 //
	 fLogicWorld->SetVisAttributes (G4VisAttributes::Invisible);
	 fLogicCalor->SetVisAttributes(G4VisAttributes::Invisible);
   
	 G4VisAttributes* magenta=new G4VisAttributes(true,G4Colour(1.,0.,1.));
	 G4VisAttributes* yellow=new G4VisAttributes(true,G4Colour(1.,1.,0.));
	 fLogicAr3->SetVisAttributes(magenta);
	 fLogicAbsorber1->SetVisAttributes(yellow);
	 fLogicAbsorber2->SetVisAttributes(yellow);

	 G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
	 simpleBoxVisAtt->SetVisibility(false);
	 fLogicCalor->SetVisAttributes(simpleBoxVisAtt);
       }
   }

   // PCB3
  fSolidPCB3 = 0; fLogicPCB3 = 0; fPhysiPCB3 = 0;
  if (fPCB3Thickness > 0.){
    fSolidPCB3 = new G4Box("PCB3",
			 fPCB3Thickness/2.,fCalorSizeYZ/2,fCalorSizeYZ/2);
                               
    fLogicPCB3 = new G4LogicalVolume(fSolidPCB3,
				     fPCB3Material,
				     fPCB3Material->GetName());
   
     
    posx = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + fAr2Thickness + fPCB2Thickness + + fAirThickness + fAbsorber3Thickness + fAl3Thickness + fAr3Thickness + (fPCB3Thickness/2.); 
  
    fPhysiPCB3= new G4PVPlacement(0,  //no rotation
				 G4ThreeVector(posx,0.,0.),       //at (0,0,0)
				 fLogicPCB3,           //its logical volume
				 "PCB3",               //its name
				 fLogicCalor,           //its mother  volume
				 false,               //no boolean operation
				 0);                    //copy number 
   
  
    PrintCalorParameters();     
  
    //                                        
    // Visualization attributes
    //
    fLogicWorld->SetVisAttributes (G4VisAttributes::Invisible);
    fLogicCalor->SetVisAttributes(G4VisAttributes::Invisible);
   
    G4VisAttributes* greenn=new G4VisAttributes(true,G4Colour(0.0, 1.0, 0.0));
    G4VisAttributes* yellow=new G4VisAttributes(true,G4Colour(1.,1.,0.));
    fLogicPCB3->SetVisAttributes(greenn);
    fLogicAbsorber1->SetVisAttributes(yellow);
    fLogicAbsorber2->SetVisAttributes(yellow);
    fLogicAbsorber3->SetVisAttributes(yellow);

    G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    simpleBoxVisAtt->SetVisibility(false);
    fLogicCalor->SetVisAttributes(simpleBoxVisAtt);

  }


   // Absorber4
    fSolidAbsorber4=0; fLogicAbsorber4=0; fPhysiAbsorber4=0;  
  
    if (fAbsorber4Thickness > 0.) 
      { fSolidAbsorber4 = new G4Box("Absorber4",                //its name
				    fAbsorber4Thickness/2,fCalorSizeYZ/2,fCalorSizeYZ/2); 
                          
	fLogicAbsorber4 = new G4LogicalVolume(fSolidAbsorber4,    //its solid
					      fAbsorber4Material, //its material
					      fAbsorber4Material->GetName());//name

	posx = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + fAr2Thickness + fPCB2Thickness + fAirThickness + fAbsorber3Thickness + fAl3Thickness + fAr3Thickness + fPCB3Thickness + fAirThickness + (fAbsorber4Thickness/2.); 
                               
	fPhysiAbsorber4 = new G4PVPlacement(0,                   //no rotation
					    G4ThreeVector(posx,0.,0.),  //its position
					    fLogicAbsorber4,     //its logical volume
					    fAbsorber4Material->GetName(), //its name
					    fLogicCalor,          //its mother
					    false,               //no boulean operat
					    0);                   //copy number
                                        
      }


 // Al4
  fSolidAl4 = 0; fLogicAl4 = 0; fPhysiAl4 = 0;
  if (fAl4Thickness > 0.){
    fSolidAl4 = new G4Box("Al4",
			 fAl4Thickness/2.,fCalorSizeYZ/2,fCalorSizeYZ/2);
                               
    fLogicAl4 = new G4LogicalVolume(fSolidAl4,
				   fAl4Material,
				   fAl4Material->GetName());
   
     
    posx = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + fAr2Thickness + fPCB2Thickness + + fAirThickness + fAbsorber3Thickness + fAl3Thickness + fAr3Thickness + fPCB3Thickness + fAirThickness + fAbsorber4Thickness + (fAl4Thickness/2.); 
    
    fPhysiAl4= new G4PVPlacement(0,  //no rotation
				 G4ThreeVector(posx,0.,0.),       //at (0,0,0)
				 fLogicAl4,           //its logical volume
				 "Al4",               //its name
				 fLogicCalor,           //its mother  volume
				 false,               //no boolean operation
				 0);                    //copy number 
    
    
    
    
    PrintCalorParameters();     
  
    //                                        
    // Visualization attributes
    //
    fLogicWorld->SetVisAttributes (G4VisAttributes::Invisible);
    fLogicCalor->SetVisAttributes(G4VisAttributes::Invisible);
   
    G4VisAttributes* blue=new G4VisAttributes(true,G4Colour(0.,0.,1.));
    G4VisAttributes* yellow=new G4VisAttributes(true,G4Colour(1.,1.,0.));
    fLogicAl4->SetVisAttributes(blue);
 
    fLogicAbsorber1->SetVisAttributes(yellow);
    fLogicAbsorber2->SetVisAttributes(yellow);



    G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    simpleBoxVisAtt->SetVisibility(false);
    fLogicCalor->SetVisAttributes(simpleBoxVisAtt);

  }

   //---------------------------------------------Ar4--------------------------------------------------------
   for(G4int j=0;j<10;j++){
     for(G4int i=0;i<10;i++)
       {
       

	 fSolidAr4[i][j] = new G4Box("Ar4",
				      fAr4Thickness/2.,fCalorSizeYZ/20,fCalorSizeYZ/20);
                               
	 fLogicAr4 = new G4LogicalVolume(fSolidAr4[i][j],
					  fAr4Material,
					  fAr4Material->GetName());
   
     
	 posx = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + fAr2Thickness + fPCB2Thickness + + fAirThickness + fAbsorber3Thickness + fAl3Thickness + fAr3Thickness + fPCB3Thickness + fAirThickness + fAbsorber4Thickness + fAl4Thickness + (fAr4Thickness/2.); 
	 
	 posy=-45+j*10;
	 posz=-45+i*10;

	 fPhysiAr4[i][j] = new G4PVPlacement(0,                  //no rotation
					    G4ThreeVector(posx,posy,posz),       //at (0,0,0)
					    fLogicAr4,           //its logical volume
					    "Ar4",               //its name
					    fLogicCalor,           //its mother  volume
					    false,               //no boolean operation
					    i);                    //copy number 



    
	 PrintCalorParameters();     
  
	 //                                        
	 // Visualization attributes
	 //
	 fLogicWorld->SetVisAttributes (G4VisAttributes::Invisible);
	 fLogicCalor->SetVisAttributes(G4VisAttributes::Invisible);
   
	 G4VisAttributes* magenta=new G4VisAttributes(true,G4Colour(1.,0.,1.));
	 G4VisAttributes* yellow=new G4VisAttributes(true,G4Colour(1.,1.,0.));
	 fLogicAr4->SetVisAttributes(magenta);
	 fLogicAbsorber1->SetVisAttributes(yellow);
	 fLogicAbsorber2->SetVisAttributes(yellow);

	 G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
	 simpleBoxVisAtt->SetVisibility(false);
	 fLogicCalor->SetVisAttributes(simpleBoxVisAtt);
       }
   }

   // PCB4
  fSolidPCB4 = 0; fLogicPCB4 = 0; fPhysiPCB4 = 0;
  if (fPCB4Thickness > 0.){
    fSolidPCB4 = new G4Box("PCB4",
			 fPCB4Thickness/2.,fCalorSizeYZ/2,fCalorSizeYZ/2);
                               
    fLogicPCB4 = new G4LogicalVolume(fSolidPCB4,
				     fPCB4Material,
				     fPCB4Material->GetName());
   
     
    posx = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + fAr2Thickness + fPCB2Thickness + + fAirThickness + fAbsorber3Thickness + fAl3Thickness + fAr3Thickness + fPCB3Thickness + fAirThickness + fAbsorber4Thickness + fAl4Thickness + fAr4Thickness + (fPCB4Thickness/2.); 
    
    fPhysiPCB4= new G4PVPlacement(0,  //no rotation
				 G4ThreeVector(posx,0.,0.),       //at (0,0,0)
				 fLogicPCB4,           //its logical volume
				 "PCB4",               //its name
				 fLogicCalor,           //its mother  volume
				 false,               //no boolean operation
				 0);                    //copy number 
   
  
    PrintCalorParameters();     
  
    //                                        
    // Visualization attributes
    //
    fLogicWorld->SetVisAttributes (G4VisAttributes::Invisible);
    fLogicCalor->SetVisAttributes(G4VisAttributes::Invisible);
   
    G4VisAttributes* greenn=new G4VisAttributes(true,G4Colour(0.0, 1.0, 0.0));
    G4VisAttributes* yellow=new G4VisAttributes(true,G4Colour(1.,1.,0.));
    fLogicPCB4->SetVisAttributes(greenn);
    fLogicAbsorber1->SetVisAttributes(yellow);
    fLogicAbsorber2->SetVisAttributes(yellow);
    fLogicAbsorber3->SetVisAttributes(yellow);

    G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    simpleBoxVisAtt->SetVisibility(false);
    fLogicCalor->SetVisAttributes(simpleBoxVisAtt);

  }


   // Absorber5
    fSolidAbsorber5=0; fLogicAbsorber5=0; fPhysiAbsorber5=0;  
  
    if (fAbsorber5Thickness > 0.) 
      { fSolidAbsorber5 = new G4Box("Absorber5",                //its name
				    fAbsorber5Thickness/2,fCalorSizeYZ/2,fCalorSizeYZ/2); 
                          
	fLogicAbsorber5 = new G4LogicalVolume(fSolidAbsorber5,    //its solid
					      fAbsorber5Material, //its material
					      fAbsorber5Material->GetName());//name

	posx = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + fAr2Thickness + fPCB2Thickness + fAirThickness + fAbsorber3Thickness + fAl3Thickness + fAr3Thickness + fPCB3Thickness + fAirThickness + fAbsorber4Thickness + fAl4Thickness + fAr4Thickness + fPCB4Thickness + fAirThickness + (fAbsorber5Thickness/2.); 

                                                
	fPhysiAbsorber5 = new G4PVPlacement(0,                   //no rotation
					    G4ThreeVector(posx,0.,0.),  //its position
					    fLogicAbsorber5,     //its logical volume
					    fAbsorber5Material->GetName(), //its name
					    fLogicCalor,          //its mother
					    false,               //no boulean operat
					    0);                   //copy number
                                        
      }



 // Al5
  fSolidAl5 = 0; fLogicAl5 = 0; fPhysiAl5 = 0;
  if (fAl5Thickness > 0.){
    fSolidAl5 = new G4Box("Al5",
			 fAl5Thickness/2.,fCalorSizeYZ/2,fCalorSizeYZ/2);
                               
    fLogicAl5 = new G4LogicalVolume(fSolidAl5,
				   fAl5Material,
				   fAl5Material->GetName());
   
     
    posx = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + fAr2Thickness + fPCB2Thickness + fAirThickness + fAbsorber3Thickness + fAl3Thickness + fAr3Thickness + fPCB3Thickness + fAirThickness + fAbsorber4Thickness + fAl4Thickness + fAr4Thickness + fPCB4Thickness + fAirThickness + fAbsorber5Thickness + (fAl5Thickness/2.); 

    fPhysiAl5= new G4PVPlacement(0,  //no rotation
				 G4ThreeVector(posx,0.,0.),       //at (0,0,0)
				 fLogicAl5,           //its logical volume
				 "Al5",               //its name
				 fLogicCalor,           //its mother  volume
				 false,               //no boolean operation
				 0);                    //copy number 
    
    
    
    
    PrintCalorParameters();     
  
    //                                        
    // Visualization attributes
    //
    fLogicWorld->SetVisAttributes (G4VisAttributes::Invisible);
    fLogicCalor->SetVisAttributes(G4VisAttributes::Invisible);
   
    G4VisAttributes* blue=new G4VisAttributes(true,G4Colour(0.,0.,1.));
    G4VisAttributes* yellow=new G4VisAttributes(true,G4Colour(1.,1.,0.));
    fLogicAl5->SetVisAttributes(blue);
 
    fLogicAbsorber1->SetVisAttributes(yellow);
    fLogicAbsorber2->SetVisAttributes(yellow);



    G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    simpleBoxVisAtt->SetVisibility(false);
    fLogicCalor->SetVisAttributes(simpleBoxVisAtt);

  }

   //---------------------------------------------Ar5--------------------------------------------------------
   for(G4int j=0;j<10;j++){
     for(G4int i=0;i<10;i++)
       {
       

	 fSolidAr5[i][j] = new G4Box("Ar5",
				      fAr5Thickness/2.,fCalorSizeYZ/20,fCalorSizeYZ/20);
                               
	 fLogicAr5 = new G4LogicalVolume(fSolidAr5[i][j],
					  fAr5Material,
					  fAr5Material->GetName());
   
     
	 posx = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + fAr2Thickness + fPCB2Thickness + fAirThickness + fAbsorber3Thickness + fAl3Thickness + fAr3Thickness + fPCB3Thickness + fAirThickness + fAbsorber4Thickness + fAl4Thickness + fAr4Thickness + fPCB4Thickness + fAirThickness + fAbsorber5Thickness + fAl5Thickness + (fAr5Thickness/2.); 
	 	 
	 posy=-45+j*10;
	 posz=-45+i*10;

	 fPhysiAr5[i][j] = new G4PVPlacement(0,                  //no rotation
					    G4ThreeVector(posx,posy,posz),       //at (0,0,0)
					    fLogicAr5,           //its logical volume
					    "Ar5",               //its name
					    fLogicCalor,           //its mother  volume
					    false,               //no boolean operation
					    i);                    //copy number 



    
	 PrintCalorParameters();     
  
	 //                                        
	 // Visualization attributes
	 //
	 fLogicWorld->SetVisAttributes (G4VisAttributes::Invisible);
	 fLogicCalor->SetVisAttributes(G4VisAttributes::Invisible);
   
	 G4VisAttributes* magenta=new G4VisAttributes(true,G4Colour(1.,0.,1.));
	 G4VisAttributes* yellow=new G4VisAttributes(true,G4Colour(1.,1.,0.));
	 fLogicAr5->SetVisAttributes(magenta);
	 fLogicAbsorber1->SetVisAttributes(yellow);
	 fLogicAbsorber2->SetVisAttributes(yellow);

	 G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
	 simpleBoxVisAtt->SetVisibility(false);
	 fLogicCalor->SetVisAttributes(simpleBoxVisAtt);
       }
   }

   // PCB5
  fSolidPCB5 = 0; fLogicPCB5 = 0; fPhysiPCB5 = 0;
  if (fPCB5Thickness > 0.){
    fSolidPCB5 = new G4Box("PCB5",
			   fPCB5Thickness/2.,fCalorSizeYZ/2,fCalorSizeYZ/2);
                               
    fLogicPCB5 = new G4LogicalVolume(fSolidPCB5,
				     fPCB5Material,
				     fPCB5Material->GetName());
   
     
    posx = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + fAr2Thickness + fPCB2Thickness + fAirThickness + fAbsorber3Thickness + fAl3Thickness + fAr3Thickness + fPCB3Thickness + fAirThickness + fAbsorber4Thickness + fAl4Thickness + fAr4Thickness + fPCB4Thickness + fAirThickness + fAbsorber5Thickness + fAl5Thickness + fAr5Thickness + (fPCB5Thickness/2.); 

    fPhysiPCB5= new G4PVPlacement(0,  //no rotation
				  G4ThreeVector(posx,0.,0.),       //at (0,0,0)
				  fLogicPCB5,           //its logical volume
				  "PCB5",               //its name
				  fLogicCalor,           //its mother  volume
				  false,               //no boolean operation
				  0);                    //copy number 
   
  
    PrintCalorParameters();     
  
    //                                        
    // Visualization attributes
    //
    fLogicWorld->SetVisAttributes (G4VisAttributes::Invisible);
    fLogicCalor->SetVisAttributes(G4VisAttributes::Invisible);
   
    G4VisAttributes* greenn=new G4VisAttributes(true,G4Colour(0.0, 1.0, 0.0));
    G4VisAttributes* yellow=new G4VisAttributes(true,G4Colour(1.,1.,0.));
    fLogicPCB5->SetVisAttributes(greenn);
    fLogicAbsorber1->SetVisAttributes(yellow);
    fLogicAbsorber2->SetVisAttributes(yellow);
    fLogicAbsorber3->SetVisAttributes(yellow);

    G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    simpleBoxVisAtt->SetVisibility(false);
    fLogicCalor->SetVisAttributes(simpleBoxVisAtt);

  }
  
  // Absorber6
  fSolidAbsorber6=0; fLogicAbsorber6=0; fPhysiAbsorber6=0;  
  
    if (fAbsorber6Thickness > 0.) 
      { fSolidAbsorber6 = new G4Box("Absorber6",                //its name
				    fAbsorber6Thickness/2,fCalorSizeYZ/2,fCalorSizeYZ/2); 
                          
	fLogicAbsorber6 = new G4LogicalVolume(fSolidAbsorber6,    //its solid
					      fAbsorber6Material, //its material
					      fAbsorber6Material->GetName());//name


	posx = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + fAr2Thickness + fPCB2Thickness + fAirThickness + fAbsorber3Thickness + fAl3Thickness + fAr3Thickness + fPCB3Thickness + fAirThickness + fAbsorber4Thickness + fAl4Thickness + fAr4Thickness + fPCB4Thickness + fAirThickness + fAbsorber5Thickness + fAl5Thickness + fAr5Thickness + fPCB5Thickness + fAirThickness + (fAbsorber6Thickness/2.); 


	fPhysiAbsorber6 = new G4PVPlacement(0,                   //no rotation
					    G4ThreeVector(posx,0.,0.),  //its position
					    fLogicAbsorber6,     //its logical volume
					    fAbsorber6Material->GetName(), //its name
					    fLogicCalor,          //its mother
					    false,               //no boulean operat
					    0);                   //copy number
                                        
      }





 // Al6
  fSolidAl6 = 0; fLogicAl6 = 0; fPhysiAl6 = 0;
  if (fAl6Thickness > 0.){
    fSolidAl6 = new G4Box("Al6",
			 fAl6Thickness/2.,fCalorSizeYZ/2,fCalorSizeYZ/2);
                               
    fLogicAl6 = new G4LogicalVolume(fSolidAl6,
				   fAl6Material,
				   fAl6Material->GetName());
   
     
    posx = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + fAr2Thickness + fPCB2Thickness + fAirThickness + fAbsorber3Thickness + fAl3Thickness + fAr3Thickness + fPCB3Thickness + fAirThickness + fAbsorber4Thickness + fAl4Thickness + fAr4Thickness + fPCB4Thickness + fAirThickness + fAbsorber5Thickness + fAl5Thickness + fAr5Thickness + fPCB5Thickness + fAirThickness + fAbsorber6Thickness + (fAl6Thickness/2.); 

    fPhysiAl6= new G4PVPlacement(0,  //no rotation
				 G4ThreeVector(posx,0.,0.),       //at (0,0,0)
				 fLogicAl6,           //its logical volume
				 "Al6",               //its name
				 fLogicCalor,           //its mother  volume
				 false,               //no boolean operation
				 0);                    //copy number 
    
    
    
    
    PrintCalorParameters();     
  
    //                                        
    // Visualization attributes
    //
    fLogicWorld->SetVisAttributes (G4VisAttributes::Invisible);
    fLogicCalor->SetVisAttributes(G4VisAttributes::Invisible);
   
    G4VisAttributes* blue=new G4VisAttributes(true,G4Colour(0.,0.,1.));
    G4VisAttributes* yellow=new G4VisAttributes(true,G4Colour(1.,1.,0.));
    fLogicAl6->SetVisAttributes(blue);
 
    fLogicAbsorber1->SetVisAttributes(yellow);
    fLogicAbsorber2->SetVisAttributes(yellow);



    G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    simpleBoxVisAtt->SetVisibility(false);
    fLogicCalor->SetVisAttributes(simpleBoxVisAtt);

  }

   //---------------------------------------------Ar6--------------------------------------------------------
   for(G4int j=0;j<10;j++){
     for(G4int i=0;i<10;i++)
       {
       

	 fSolidAr6[i][j] = new G4Box("Ar6",
				      fAr6Thickness/2.,fCalorSizeYZ/20,fCalorSizeYZ/20);
                               
	 fLogicAr6 = new G4LogicalVolume(fSolidAr6[i][j],
					  fAr6Material,
					  fAr6Material->GetName());
   
     
	 posx = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + fAr2Thickness + fPCB2Thickness + fAirThickness + fAbsorber3Thickness + fAl3Thickness + fAr3Thickness + fPCB3Thickness + fAirThickness + fAbsorber4Thickness + fAl4Thickness + fAr4Thickness + fPCB4Thickness + fAirThickness + fAbsorber5Thickness + fAl5Thickness + fAr5Thickness + fPCB5Thickness + fAirThickness + fAbsorber6Thickness + fAl6Thickness + (fAr6Thickness/2.); 

	 posy=-45+j*10;
	 posz=-45+i*10;

	 fPhysiAr6[i][j] = new G4PVPlacement(0,                  //no rotation
					    G4ThreeVector(posx,posy,posz),       //at (0,0,0)
					    fLogicAr6,           //its logical volume
					    "Ar6",               //its name
					    fLogicCalor,           //its mother  volume
					    false,               //no boolean operation
					    i);                    //copy number 



    
	 PrintCalorParameters();     
  
	 //                                        
	 // Visualization attributes
	 //
	 fLogicWorld->SetVisAttributes (G4VisAttributes::Invisible);
	 fLogicCalor->SetVisAttributes(G4VisAttributes::Invisible);
   
	 G4VisAttributes* magenta=new G4VisAttributes(true,G4Colour(1.,0.,1.));
	 G4VisAttributes* yellow=new G4VisAttributes(true,G4Colour(1.,1.,0.));
	 fLogicAr6->SetVisAttributes(magenta);
	 fLogicAbsorber1->SetVisAttributes(yellow);
	 fLogicAbsorber2->SetVisAttributes(yellow);

	 G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
	 simpleBoxVisAtt->SetVisibility(false);
	 fLogicCalor->SetVisAttributes(simpleBoxVisAtt);
       }
   }

   // PCB6
  fSolidPCB6 = 0; fLogicPCB6 = 0; fPhysiPCB6 = 0;
  if (fPCB6Thickness > 0.){
    fSolidPCB6 = new G4Box("PCB6",
			 fPCB6Thickness/2.,fCalorSizeYZ/2,fCalorSizeYZ/2);
                               
    fLogicPCB6 = new G4LogicalVolume(fSolidPCB6,
				     fPCB6Material,
				     fPCB6Material->GetName());
   
     
    posx = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + fAr2Thickness + fPCB2Thickness + fAirThickness + fAbsorber3Thickness + fAl3Thickness + fAr3Thickness + fPCB3Thickness + fAirThickness + fAbsorber4Thickness + fAl4Thickness + fAr4Thickness + fPCB4Thickness + fAirThickness + fAbsorber5Thickness + fAl5Thickness + fAr5Thickness + fPCB5Thickness + fAirThickness + fAbsorber6Thickness + fAl6Thickness + fAr6Thickness + (fPCB6Thickness/2.); 

    fPhysiPCB6= new G4PVPlacement(0,  //no rotation
				 G4ThreeVector(posx,0.,0.),       //at (0,0,0)
				 fLogicPCB6,           //its logical volume
				 "PCB6",               //its name
				 fLogicCalor,           //its mother  volume
				 false,               //no boolean operation
				 0);                    //copy number 
   
  
    PrintCalorParameters();     
  
    //                                        
    // Visualization attributes
    //
    fLogicWorld->SetVisAttributes (G4VisAttributes::Invisible);
    fLogicCalor->SetVisAttributes(G4VisAttributes::Invisible);
   
    G4VisAttributes* greenn=new G4VisAttributes(true,G4Colour(0.0, 1.0, 0.0));
    G4VisAttributes* yellow=new G4VisAttributes(true,G4Colour(1.,1.,0.));
    fLogicPCB6->SetVisAttributes(greenn);
    fLogicAbsorber1->SetVisAttributes(yellow);
    fLogicAbsorber2->SetVisAttributes(yellow);
    fLogicAbsorber3->SetVisAttributes(yellow);

    G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    simpleBoxVisAtt->SetVisibility(false);
    fLogicCalor->SetVisAttributes(simpleBoxVisAtt);

  }

   // Absorber7
    fSolidAbsorber7=0; fLogicAbsorber7=0; fPhysiAbsorber7=0;  
  
    if (fAbsorber7Thickness > 0.) 
      { fSolidAbsorber7 = new G4Box("Absorber7",                //its name
				    fAbsorber7Thickness/2,fCalorSizeYZ/2,fCalorSizeYZ/2); 
                          
	fLogicAbsorber7 = new G4LogicalVolume(fSolidAbsorber7,    //its solid
					      fAbsorber7Material, //its material
					      fAbsorber7Material->GetName());//name

	posx = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + fAr2Thickness + fPCB2Thickness + fAirThickness + fAbsorber3Thickness + fAl3Thickness + fAr3Thickness + fPCB3Thickness + fAirThickness + fAbsorber4Thickness + fAl4Thickness + fAr4Thickness + fPCB4Thickness + fAirThickness + fAbsorber5Thickness + fAl5Thickness + fAr5Thickness + fPCB5Thickness + fAirThickness + fAbsorber6Thickness + fAl6Thickness + fAr6Thickness + fPCB6Thickness + fAirThickness + (fAbsorber7Thickness/2.); 

                                                
	fPhysiAbsorber7 = new G4PVPlacement(0,                   //no rotation
					    G4ThreeVector(posx,0.,0.),  //its position
					    fLogicAbsorber7,     //its logical volume
					    fAbsorber7Material->GetName(), //its name
					    fLogicCalor,          //its mother
					    false,               //no boulean operat
					    0);                   //copy number
                                        

	G4VisAttributes* yellow=new G4VisAttributes(true,G4Colour(1.,1.,0.));
	fLogicAbsorber2->SetVisAttributes(yellow);
	fLogicAbsorber3->SetVisAttributes(yellow);
	fLogicAbsorber4->SetVisAttributes(yellow);
	fLogicAbsorber5->SetVisAttributes(yellow);
	fLogicAbsorber6->SetVisAttributes(yellow);
	fLogicAbsorber7->SetVisAttributes(yellow);


      }



  //always return the physical World
  //
  return fPhysiWorld;

  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintCalorParameters()
{
 /* G4cout << "\n------------------------------------------------------------"
         << "\n---> The calorimeter is " << fEnvThickness << " layers of: [ "
         << fAbsorber1Thickness/mm << "mm of " << fAbsorber1Material->GetName() 
         << " + "
         << fGapThickness/mm << "mm of " << fGapMaterial->GetName() << " ] " 
         << "\n------------------------------------------------------------\n";
*/}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorber1Material(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial =
  G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);  
  if (pttoMaterial)
  {
      fAbsorber1Material = pttoMaterial;
      if ( fLogicAbsorber1 )
      {
          fLogicAbsorber1->SetMaterial(fAbsorber1Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();
      }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorber2Material(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial =
  G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);      
  if (pttoMaterial)
  {
      fAbsorber2Material = pttoMaterial;
      if ( fLogicAbsorber2 )
      {
          fLogicAbsorber2->SetMaterial(fAbsorber2Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();
      }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorber3Material(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial =
  G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);      
  if (pttoMaterial)
  {
      fAbsorber3Material = pttoMaterial;
      if ( fLogicAbsorber3 )
      {
          fLogicAbsorber3->SetMaterial(fAbsorber3Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();
      }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAbsorber4Material(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial =
  G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);      
  if (pttoMaterial)
  {
      fAbsorber4Material = pttoMaterial;
      if ( fLogicAbsorber4 )
      {
          fLogicAbsorber4->SetMaterial(fAbsorber4Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();
      }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAbsorber5Material(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial =
  G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);      
  if (pttoMaterial)
  {
      fAbsorber5Material = pttoMaterial;
      if ( fLogicAbsorber5 )
      {
          fLogicAbsorber5->SetMaterial(fAbsorber5Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();
      }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAbsorber6Material(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial =
  G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);      
  if (pttoMaterial)
  {
      fAbsorber6Material = pttoMaterial;
      if ( fLogicAbsorber6 )
      {
          fLogicAbsorber6->SetMaterial(fAbsorber6Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();
      }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAbsorber7Material(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial =
  G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);      
  if (pttoMaterial)
  {
      fAbsorber7Material = pttoMaterial;
      if ( fLogicAbsorber7 )
      {
          fLogicAbsorber7->SetMaterial(fAbsorber7Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();
      }
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAl1Material(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =  
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  if (pttoMaterial)
    {
      fAl1Material = pttoMaterial;
      if ( fLogicAl1 )
	{
          fLogicAl1->SetMaterial(fAl1Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();

	}
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAl2Material(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =  
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  if (pttoMaterial)
    {
      fAl2Material = pttoMaterial;
      if ( fLogicAl2 )
	{
          fLogicAl2->SetMaterial(fAl2Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();

	}
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAl3Material(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =  
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  if (pttoMaterial)
    {
      fAl3Material = pttoMaterial;
      if ( fLogicAl3 )
	{
          fLogicAl3->SetMaterial(fAl3Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();

	}
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAl4Material(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =  
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  if (pttoMaterial)
    {
      fAl4Material = pttoMaterial;
      if ( fLogicAl4 )
	{
          fLogicAl4->SetMaterial(fAl4Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();

	}
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAl5Material(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =  
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  if (pttoMaterial)
    {
      fAl5Material = pttoMaterial;
      if ( fLogicAl5 )
	{
          fLogicAl5->SetMaterial(fAl5Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();

	}
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAl6Material(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =  
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  if (pttoMaterial)
    {
      fAl6Material = pttoMaterial;
      if ( fLogicAl6 )
	{
          fLogicAl6->SetMaterial(fAl6Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();

	}
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAr1Material(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =  
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);  
  if (pttoMaterial)
    {
      fAr1Material = pttoMaterial;
      if ( fLogicAr1 )
	{
	  fLogicAr1->SetMaterial(fAr1Material);
	  G4RunManager::GetRunManager()->PhysicsHasBeenModified();

	}
    }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAr2Material(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =  
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  if (pttoMaterial)
    {
      fAr2Material = pttoMaterial;
      if ( fLogicAr2 )
	{
          fLogicAr2->SetMaterial(fAr2Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();

	}
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAr3Material(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =  
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  if (pttoMaterial)
    {
      fAr3Material = pttoMaterial;
      if ( fLogicAr3 )
	{
          fLogicAr3->SetMaterial(fAr3Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();

	}
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAr4Material(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =  
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  if (pttoMaterial)
    {
      fAr4Material = pttoMaterial;
      if ( fLogicAr4 )
	{
          fLogicAr4->SetMaterial(fAr4Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();

	}
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAr5Material(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =  
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  if (pttoMaterial)
    {
      fAr5Material = pttoMaterial;
      if ( fLogicAr5 )
	{
          fLogicAr5->SetMaterial(fAr5Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();

	}
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAr6Material(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =  
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  if (pttoMaterial)
    {
      fAr6Material = pttoMaterial;
      if ( fLogicAr6 )
	{
          fLogicAr6->SetMaterial(fAr6Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();

	}
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetPCB1Material()
{
  // build the material 
 G4Material* pttoMaterial = new G4Material("G10",1.700*g/cm3,4);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(14), 1);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(8) , 2);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(6) , 3);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(1) , 3); 
 
 if (pttoMaterial)
    {
      fPCB1Material = pttoMaterial;
      if ( fLogicPCB1 )
	{
          fLogicPCB1->SetMaterial(fPCB1Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();

	}
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetPCB2Material()
{
  // build the material 
 G4Material* pttoMaterial = new G4Material("G10_0",1.700*g/cm3,4);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(14), 1);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(8) , 2);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(6) , 3);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(1) , 3); 
 
 if (pttoMaterial)
    {
      fPCB2Material = pttoMaterial;
      if ( fLogicPCB2 )
	{
          fLogicPCB2->SetMaterial(fPCB2Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();

	}
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetPCB3Material()
{
  // build the material 
 G4Material* pttoMaterial = new G4Material("G10_1",1.700*g/cm3,4);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(14), 1);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(8) , 2);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(6) , 3);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(1) , 3); 
 
 if (pttoMaterial)
    {
      fPCB3Material = pttoMaterial;
      if ( fLogicPCB3 )
	{
          fLogicPCB3->SetMaterial(fPCB3Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();

	}
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetPCB4Material()
{
  // build the material 
 G4Material* pttoMaterial = new G4Material("G10_2",1.700*g/cm3,4);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(14), 1);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(8) , 2);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(6) , 3);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(1) , 3); 
 
 if (pttoMaterial)
    {
      fPCB4Material = pttoMaterial;
      if ( fLogicPCB4 )
	{
          fLogicPCB4->SetMaterial(fPCB4Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();

	}
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetPCB5Material()
{
  // build the material 
 G4Material* pttoMaterial = new G4Material("G10_3",1.700*g/cm3,4);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(14), 1);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(8) , 2);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(6) , 3);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(1) , 3); 
 
 if (pttoMaterial)
    {
      fPCB5Material = pttoMaterial;
      if ( fLogicPCB5 )
	{
          fLogicPCB5->SetMaterial(fPCB5Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();

	}
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetPCB6Material()
{
  // build the material 
 G4Material* pttoMaterial = new G4Material("G10_4",1.700*g/cm3,4);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(14), 1);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(8) , 2);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(6) , 3);
 pttoMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement(1) , 3); 
 
 if (pttoMaterial)
    {
      fPCB6Material = pttoMaterial;
      if ( fLogicPCB6 )
	{
          fLogicPCB6->SetMaterial(fPCB6Material);
          G4RunManager::GetRunManager()->PhysicsHasBeenModified();

	}
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAbsorber1Thickness(G4double val)
{
  // change Absorber1 thickness and recompute the calorimeter parameters
  fAbsorber1Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAbsorber2Thickness(G4double val)
{
  // change Absorber2 thickness and recompute the calorimeter parameters
  fAbsorber2Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAbsorber3Thickness(G4double val)
{
  // change Absorber3 thickness and recompute the calorimeter parameters
  fAbsorber3Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAbsorber4Thickness(G4double val)
{
  // change Absorber4 thickness and recompute the calorimeter parameters
  fAbsorber4Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAbsorber5Thickness(G4double val)
{
  // change Absorber5 thickness and recompute the calorimeter parameters
  fAbsorber5Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAbsorber6Thickness(G4double val)
{
  // change Absorber6 thickness and recompute the calorimeter parameters
  fAbsorber6Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAbsorber7Thickness(G4double val)
{
  // change Absorber7 thickness and recompute the calorimeter parameters
  fAbsorber7Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAl1Thickness(G4double val)
{
  // change Al thickness and recompute the calorimeter parameters
  fAl1Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAl2Thickness(G4double val)
{
  // change Al thickness and recompute the calorimeter parameters
  fAl2Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAl3Thickness(G4double val)
{
  // change Al thickness and recompute the calorimeter parameters
  fAl3Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAl4Thickness(G4double val)
{
  // change Al thickness and recompute the calorimeter parameters
  fAl4Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAl5Thickness(G4double val)
{
  // change Al thickness and recompute the calorimeter parameters
  fAl5Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAl6Thickness(G4double val)
{
  // change Al thickness and recompute the calorimeter parameters
  fAl6Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAr1Thickness(G4double val)
{
  // change Ar thickness and recompute the calorimeter parameters
  fAr1Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAr2Thickness(G4double val)
{
  // change Ar thickness and recompute the calorimeter parameters
  fAr2Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAr3Thickness(G4double val)
{
  // change Ar thickness and recompute the calorimeter parameters
  fAr3Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAr4Thickness(G4double val)
{
  // change Ar thickness and recompute the calorimeter parameters
  fAr4Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAr5Thickness(G4double val)
{
  // change Ar thickness and recompute the calorimeter parameters
  fAr5Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAr6Thickness(G4double val)
{
  // change Ar thickness and recompute the calorimeter parameters
  fAr6Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetPCB1Thickness(G4double val)
{
  // change PCB thickness and recompute the calorimeter parameters
  fPCB1Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetPCB2Thickness(G4double val)
{
  // change PCB thickness and recompute the calorimeter parameters
  fPCB2Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetPCB3Thickness(G4double val)
{
  // change PCB thickness and recompute the calorimeter parameters
  fPCB3Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetPCB4Thickness(G4double val)
{
  // change PCB thickness and recompute the calorimeter parameters
  fPCB4Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetPCB5Thickness(G4double val)
{
  // change PCB thickness and recompute the calorimeter parameters
  fPCB5Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetPCB6Thickness(G4double val)
{
  // change PCB thickness and recompute the calorimeter parameters
  fPCB6Thickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetCalorSizeYZ(G4double val)
{
  // change the transverse size and recompute the calorimeter parameters
  fCalorSizeYZ = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetEnvThickness(G4int val)
{
  fEnvThickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}











//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
