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
/// \file analysis/shared/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
//
// $Id: DetectorConstruction.hh 77256 2013-11-22 10:10:23Z gcosmo $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4SDManager.hh"
#include "globals.hh"
//#include "G4VSensitiveDetector.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
    virtual ~DetectorConstruction();

  public:
     
     void SetAbsorber1Material (G4String);     
     void SetAbsorber1Thickness(G4double);     

     void SetAbsorber2Material (G4String);     
     void SetAbsorber2Thickness(G4double);     

     void SetAbsorber3Material (G4String);     
     void SetAbsorber3Thickness(G4double);     

     void SetAbsorber4Material (G4String);     
     void SetAbsorber4Thickness(G4double);     

     void SetAbsorber5Material (G4String);     
     void SetAbsorber5Thickness(G4double);     

     void SetAbsorber6Material (G4String);     
     void SetAbsorber6Thickness(G4double);     

     void SetAbsorber7Material (G4String);     
     void SetAbsorber7Thickness(G4double);     

     void SetAl1Material (G4String);     
     void SetAl1Thickness(G4double);
     
     void SetAl2Material (G4String);     
     void SetAl2Thickness(G4double);

     void SetAl3Material (G4String);     
     void SetAl3Thickness(G4double);

     void SetAl4Material (G4String);     
     void SetAl4Thickness(G4double);

     void SetAl5Material (G4String);     
     void SetAl5Thickness(G4double);

     void SetAl6Material (G4String);     
     void SetAl6Thickness(G4double);

     void SetAr1Material (G4String);     
     void SetAr1Thickness(G4double);
     
     void SetAr2Material (G4String);     
     void SetAr2Thickness(G4double);

     void SetAr3Material (G4String);     
     void SetAr3Thickness(G4double);

     void SetAr4Material (G4String);     
     void SetAr4Thickness(G4double);

     void SetAr5Material (G4String);     
     void SetAr5Thickness(G4double);

     void SetAr6Material (G4String);     
     void SetAr6Thickness(G4double);

     void SetPCB1Material ();     
     void SetPCB1Thickness(G4double);
     
     void SetPCB2Material ();     
     void SetPCB2Thickness(G4double);

     void SetPCB3Material ();     
     void SetPCB3Thickness(G4double);

     void SetPCB4Material ();     
     void SetPCB4Thickness(G4double);

     void SetPCB5Material ();     
     void SetPCB5Thickness(G4double);

     void SetPCB6Material ();     
     void SetPCB6Thickness(G4double);

     void SetCalorSizeYZ(G4double);          
     void SetEnvThickness (G4int);   
     
     virtual G4VPhysicalVolume* Construct();

  public:
  
     void PrintCalorParameters(); 
                    
     G4double GetWorldSizeX()           {return fWorldSizeX;}; 
     G4double GetWorldSizeYZ()          {return fWorldSizeYZ;};
     
     G4double GetCalorThickness()       {return fCalorThickness;}; 
     G4double GetCalorSizeYZ()          {return fCalorSizeYZ;};
      
     G4int GetEnvThickness()              {return fEnvThickness;}; 
     
     G4Material* GetAbsorber1Material()  {return fAbsorber1Material;};
     G4double    GetAbsorber1Thickness() {return fAbsorber1Thickness;};      
     
     G4Material* GetAbsorber2Material()  {return fAbsorber2Material;};
     G4double    GetAbsorber2Thickness() {return fAbsorber2Thickness;};      

     G4Material* GetAbsorber3Material()  {return fAbsorber3Material;};
     G4double    GetAbsorber3Thickness() {return fAbsorber3Thickness;};      

     G4Material* GetAbsorber4Material()  {return fAbsorber4Material;};
     G4double    GetAbsorber4Thickness() {return fAbsorber4Thickness;};      

     G4Material* GetAbsorber5Material()  {return fAbsorber5Material;};
     G4double    GetAbsorber5Thickness() {return fAbsorber5Thickness;};      

     G4Material* GetAbsorber6Material()  {return fAbsorber6Material;};
     G4double    GetAbsorber6Thickness() {return fAbsorber6Thickness;};      

     G4Material* GetAbsorber7Material()  {return fAbsorber7Material;};
     G4double    GetAbsorber7Thickness() {return fAbsorber7Thickness;};      

     G4Material* GetAl1Material()       {return fAl1Material;};
     G4double    GetAl1Thickness()      {return fAl1Thickness;};

     G4Material* GetAl2Material()       {return fAl2Material;};
     G4double    GetAl2Thickness()      {return fAl2Thickness;};
         
     G4Material* GetAl3Material()       {return fAl3Material;};
     G4double    GetAl3Thickness()      {return fAl3Thickness;};
         
     G4Material* GetAl4Material()       {return fAl4Material;};
     G4double    GetAl4Thickness()      {return fAl4Thickness;};
         
     G4Material* GetAl5Material()       {return fAl5Material;};
     G4double    GetAl5Thickness()      {return fAl5Thickness;};
         
     G4Material* GetAl6Material()       {return fAl6Material;};
     G4double    GetAl6Thickness()      {return fAl6Thickness;};
         
     G4Material* GetAr1Material()       {return fAr1Material;};
     G4double    GetAr1Thickness()      {return fAr1Thickness;};

     G4Material* GetAr2Material()       {return fAr2Material;};
     G4double    GetAr2Thickness()      {return fAr2Thickness;};
     
     G4Material* GetAr3Material()       {return fAr3Material;};
     G4double    GetAr3Thickness()      {return fAr3Thickness;};
     
     G4Material* GetAr4Material()       {return fAr4Material;};
     G4double    GetAr4Thickness()      {return fAr4Thickness;};
     
     G4Material* GetAr5Material()       {return fAr5Material;};
     G4double    GetAr5Thickness()      {return fAr5Thickness;};
     
     G4Material* GetAr6Material()       {return fAr6Material;};
     G4double    GetAr6Thickness()      {return fAr6Thickness;};
     
     G4Material* GetPCB1Material()       {return fPCB1Material;};
     G4double    GetPCB1Thickness()      {return fPCB1Thickness;};

     G4Material* GetPCB2Material()       {return fPCB2Material;};
     G4double    GetPCB2Thickness()      {return fPCB2Thickness;};
 
     G4Material* GetPCB3Material()       {return fPCB3Material;};
     G4double    GetPCB3Thickness()      {return fPCB3Thickness;};
 
     G4Material* GetPCB4Material()       {return fPCB4Material;};
     G4double    GetPCB4Thickness()      {return fPCB4Thickness;};
 
     G4Material* GetPCB5Material()       {return fPCB5Material;};
     G4double    GetPCB5Thickness()      {return fPCB5Thickness;};
 
     G4Material* GetPCB6Material()       {return fPCB6Material;};
     G4double    GetPCB6Thickness()      {return fPCB6Thickness;};
 
     G4double    GetAirThickness()      {return fAirThickness;};

     const G4VPhysicalVolume* GetphysiWorld()              {return fPhysiWorld;};           
     const G4VPhysicalVolume* GetAbsorber1()                {return fPhysiAbsorber1;};
     const G4VPhysicalVolume* GetAbsorber2()                {return fPhysiAbsorber2;};
     const G4VPhysicalVolume* GetAbsorber3()                {return fPhysiAbsorber3;};
     const G4VPhysicalVolume* GetAbsorber4()                {return fPhysiAbsorber4;};
     const G4VPhysicalVolume* GetAbsorber5()                {return fPhysiAbsorber5;};
     const G4VPhysicalVolume* GetAbsorber6()                {return fPhysiAbsorber6;};
     const G4VPhysicalVolume* GetAbsorber7()                {return fPhysiAbsorber7;};
     const G4VPhysicalVolume* GetAl1()                       {return fPhysiAl1;};
     const G4VPhysicalVolume* GetAl2()                       {return fPhysiAl2;};
     const G4VPhysicalVolume* GetAl3()                       {return fPhysiAl3;};
     const G4VPhysicalVolume* GetAl4()                       {return fPhysiAl4;};
     const G4VPhysicalVolume* GetAl5()                       {return fPhysiAl5;};
     const G4VPhysicalVolume* GetAl6()                       {return fPhysiAl6;};
     const G4VPhysicalVolume* GetAr1(int i, int j)        {return fPhysiAr1[i][j];};
     const G4VPhysicalVolume* GetAr2(int i, int j)        {return fPhysiAr2[i][j];};
     const G4VPhysicalVolume* GetAr3(int i, int j)        {return fPhysiAr3[i][j];};
     const G4VPhysicalVolume* GetAr4(int i, int j)        {return fPhysiAr4[i][j];};
     const G4VPhysicalVolume* GetAr5(int i, int j)        {return fPhysiAr5[i][j];};
     const G4VPhysicalVolume* GetAr6(int i, int j)        {return fPhysiAr6[i][j];};
     const G4VPhysicalVolume* GetPCB1()                       {return fPhysiPCB1;};
     const G4VPhysicalVolume* GetPCB2()                       {return fPhysiPCB2;};
     const G4VPhysicalVolume* GetPCB3()                       {return fPhysiPCB3;};
     const G4VPhysicalVolume* GetPCB4()                       {return fPhysiPCB4;};
     const G4VPhysicalVolume* GetPCB5()                       {return fPhysiPCB5;};
     const G4VPhysicalVolume* GetPCB6()                       {return fPhysiPCB6;};

            
  private:
     
     G4Material*        fAbsorber1Material;
     G4double           fAbsorber1Thickness;
     
     G4Material*        fAbsorber2Material;
     G4double           fAbsorber2Thickness;

     G4Material*        fAbsorber3Material;
     G4double           fAbsorber3Thickness;

     G4Material*        fAbsorber4Material;
     G4double           fAbsorber4Thickness;

     G4Material*        fAbsorber5Material;
     G4double           fAbsorber5Thickness;

     G4Material*        fAbsorber6Material;
     G4double           fAbsorber6Thickness;

     G4Material*        fAbsorber7Material;
     G4double           fAbsorber7Thickness;

     G4Material*        fAl1Material;
     G4double           fAl1Thickness;
     
     G4Material*        fAl2Material;
     G4double           fAl2Thickness;

     G4Material*        fAl3Material;
     G4double           fAl3Thickness;

     G4Material*        fAl4Material;
     G4double           fAl4Thickness;

     G4Material*        fAl5Material;
     G4double           fAl5Thickness;

     G4Material*        fAl6Material;
     G4double           fAl6Thickness;

     G4Material*        fAr1Material;
     G4double           fAr1Thickness;
     
     G4Material*        fAr2Material;
     G4double           fAr2Thickness;

     G4Material*        fAr3Material;
     G4double           fAr3Thickness;

     G4Material*        fAr4Material;
     G4double           fAr4Thickness;

     G4Material*        fAr5Material;
     G4double           fAr5Thickness;

     G4Material*        fAr6Material;
     G4double           fAr6Thickness;

     G4Material*        fPCB1Material;
     G4double           fPCB1Thickness;
     
     G4Material*        fPCB2Material;
     G4double           fPCB2Thickness;

     G4Material*        fPCB3Material;
     G4double           fPCB3Thickness;

     G4Material*        fPCB4Material;
     G4double           fPCB4Thickness;

     G4Material*        fPCB5Material;
     G4double           fPCB5Thickness;

     G4Material*        fPCB6Material;
     G4double           fPCB6Thickness;

     G4double           fAirThickness;

     G4int              fEnvThickness;
     G4double           fLayerThickness;
          
     G4double           fCalorSizeYZ;
     G4double           fCalorThickness;
     
     G4Material*        fDefaultMaterial;
     G4double           fWorldSizeYZ;
     G4double           fWorldSizeX;
            
     G4Box*             fSolidWorld;    //pointer to the solid World 
     G4LogicalVolume*   fLogicWorld;    //pointer to the logical World
     G4VPhysicalVolume* fPhysiWorld;    //pointer to the physical World

     G4Box*             fSolidCalor;    //pointer to the solid Calor 
     G4LogicalVolume*   fLogicCalor;    //pointer to the logical Calor
     G4VPhysicalVolume* fPhysiCalor;    //pointer to the physical Calor
     
   //  G4Box*             fSolidLayer;    //pointer to the solid Layer 
   //  G4LogicalVolume*   fLogicLayer;    //pointer to the logical Layer
   //  G4VPhysicalVolume* fPhysiLayer;    //pointer to the physical Layer
         
     G4Box*             fSolidAbsorber1; //pointer to the solid Absorber1
     G4LogicalVolume*   fLogicAbsorber1; //pointer to the logical Absorber1
     G4VPhysicalVolume* fPhysiAbsorber1; //pointer to the physical Absorber1
     
     G4Box*             fSolidAbsorber2; //pointer to the solid Absorber2
     G4LogicalVolume*   fLogicAbsorber2; //pointer to the logical Absorber2
     G4VPhysicalVolume* fPhysiAbsorber2; //pointer to the physical Absorber2

     G4Box*             fSolidAbsorber3; //pointer to the solid Absorber3
     G4LogicalVolume*   fLogicAbsorber3; //pointer to the logical Absorber3
     G4VPhysicalVolume* fPhysiAbsorber3; //pointer to the physical Absorber3

     G4Box*             fSolidAbsorber4; //pointer to the solid Absorber4
     G4LogicalVolume*   fLogicAbsorber4; //pointer to the logical Absorber4
     G4VPhysicalVolume* fPhysiAbsorber4; //pointer to the physical Absorber4

     G4Box*             fSolidAbsorber5; //pointer to the solid Absorber5
     G4LogicalVolume*   fLogicAbsorber5; //pointer to the logical Absorber5
     G4VPhysicalVolume* fPhysiAbsorber5; //pointer to the physical Absorber5

     G4Box*             fSolidAbsorber6; //pointer to the solid Absorber6
     G4LogicalVolume*   fLogicAbsorber6; //pointer to the logical Absorber6
     G4VPhysicalVolume* fPhysiAbsorber6; //pointer to the physical Absorber6

     G4Box*             fSolidAbsorber7; //pointer to the solid Absorber7
     G4LogicalVolume*   fLogicAbsorber7; //pointer to the logical Absorber7
     G4VPhysicalVolume* fPhysiAbsorber7; //pointer to the physical Absorber7

     G4Box*             fSolidAl1;      //pointer to the solid Al1
     G4LogicalVolume*   fLogicAl1;      //pointer to the logical Al1
     G4VPhysicalVolume* fPhysiAl1;      //pointer to the physical Al1

     G4Box*             fSolidAl2;      //pointer to the solid Al2
     G4LogicalVolume*   fLogicAl2;      //pointer to the logical Al2
     G4VPhysicalVolume* fPhysiAl2;      //pointer to the physical Al2

     G4Box*             fSolidAl3;      //pointer to the solid Al3
     G4LogicalVolume*   fLogicAl3;      //pointer to the logical Al3
     G4VPhysicalVolume* fPhysiAl3;      //pointer to the physical Al3

     G4Box*             fSolidAl4;      //pointer to the solid Al4
     G4LogicalVolume*   fLogicAl4;      //pointer to the logical Al4
     G4VPhysicalVolume* fPhysiAl4;      //pointer to the physical Al4

     G4Box*             fSolidAl5;      //pointer to the solid Al5
     G4LogicalVolume*   fLogicAl5;      //pointer to the logical Al5
     G4VPhysicalVolume* fPhysiAl5;      //pointer to the physical Al5

     G4Box*             fSolidAl6;      //pointer to the solid Al6
     G4LogicalVolume*   fLogicAl6;      //pointer to the logical Al6
     G4VPhysicalVolume* fPhysiAl6;      //pointer to the physical Al6

     G4Box*             fSolidAr1[10][10];      //pointer to the solid Ar1
     G4LogicalVolume*   fLogicAr1;      //pointer to the logical Ar1
     G4VPhysicalVolume* fPhysiAr1[10][10];      //pointer to the physical Ar1

     G4Box*             fSolidAr2[10][10];      //pointer to the solid Ar2
     G4LogicalVolume*   fLogicAr2;      //pointer to the logical Ar2
     G4VPhysicalVolume* fPhysiAr2[10][10];      //pointer to the physical Ar2

     G4Box*             fSolidAr3[10][10];      //pointer to the solid Ar3
     G4LogicalVolume*   fLogicAr3;      //pointer to the logical Ar3
     G4VPhysicalVolume* fPhysiAr3[10][10];      //pointer to the physical Ar3

     G4Box*             fSolidAr4[10][10];      //pointer to the solid Ar4
     G4LogicalVolume*   fLogicAr4;      //pointer to the logical Ar4
     G4VPhysicalVolume* fPhysiAr4[10][10];      //pointer to the physical Ar4

     G4Box*             fSolidAr5[10][10];      //pointer to the solid Ar5
     G4LogicalVolume*   fLogicAr5;      //pointer to the logical Ar5
     G4VPhysicalVolume* fPhysiAr5[10][10];      //pointer to the physical Ar5

     G4Box*             fSolidAr6[10][10];      //pointer to the solid Ar6
     G4LogicalVolume*   fLogicAr6;      //pointer to the logical Ar6
     G4VPhysicalVolume* fPhysiAr6[10][10];      //pointer to the physical Ar6

     G4Box*             fSolidPCB1;      //pointer to the solid PCB1
     G4LogicalVolume*   fLogicPCB1;      //pointer to the logical PCB1
     G4VPhysicalVolume* fPhysiPCB1;      //pointer to the physical PCB1
   
     G4Box*             fSolidPCB2;      //pointer to the solid PCB2
     G4LogicalVolume*   fLogicPCB2;      //pointer to the logical PCB2
     G4VPhysicalVolume* fPhysiPCB2;      //pointer to the physical PCB2
     
     G4Box*             fSolidPCB3;      //pointer to the solid PCB3
     G4LogicalVolume*   fLogicPCB3;      //pointer to the logical PCB3
     G4VPhysicalVolume* fPhysiPCB3;      //pointer to the physical PCB3
     
     G4Box*             fSolidPCB4;      //pointer to the solid PCB4
     G4LogicalVolume*   fLogicPCB4;      //pointer to the logical PCB4
     G4VPhysicalVolume* fPhysiPCB4;      //pointer to the physical PCB4
     
     G4Box*             fSolidPCB5;      //pointer to the solid PCB5
     G4LogicalVolume*   fLogicPCB5;      //pointer to the logical PCB5
     G4VPhysicalVolume* fPhysiPCB5;      //pointer to the physical PCB5
     
     G4Box*             fSolidPCB6;      //pointer to the solid PCB6
     G4LogicalVolume*   fLogicPCB6;      //pointer to the logical PCB6
     G4VPhysicalVolume* fPhysiPCB6;      //pointer to the physical PCB6
     
     DetectorMessenger* fDetectorMessenger;  //pointer to the Messenger
      
  private:
    
     void DefineMaterials();
     void ComputeCalorParameters();
     G4VPhysicalVolume* ConstructCalorimeter();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void DetectorConstruction::ComputeCalorParameters()
{
  // Compute derived parameters of the calorimeter
  fCalorThickness = fAbsorber1Thickness + fAl1Thickness + fAr1Thickness + fPCB1Thickness + fAirThickness + fAbsorber2Thickness + fAl2Thickness + fAr2Thickness + fPCB2Thickness + fAirThickness + fAbsorber3Thickness + fAl3Thickness + fAr3Thickness + fPCB3Thickness + fAirThickness + fAbsorber4Thickness + fAl4Thickness + fAr4Thickness + fPCB4Thickness + fAirThickness + fAbsorber5Thickness + fAl5Thickness + fAr5Thickness + fPCB5Thickness + fAirThickness + fAbsorber6Thickness + fAl6Thickness + fAr6Thickness + fPCB6Thickness + fAirThickness + fAbsorber7Thickness;

  fWorldSizeX = 1.2*fCalorThickness;
  fWorldSizeYZ = 1.2*fCalorSizeYZ;
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

