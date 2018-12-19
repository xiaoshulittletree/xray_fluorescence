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
// $Id: B5MagneticField.cc 101036 2016-11-04 09:00:23Z gcosmo $
//
/// \file B5MagneticField.cc
/// \brief Implementation of the B5MagneticField class

#include "XrayFluoMagneticField.hh"

#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include "G4AutoDelete.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

XrayFluoMagneticField::XrayFluoMagneticField()
: G4MagneticField(),
  fMessenger(nullptr), fB0(1.0*tesla)
{
  R = 125.0*um;//beam Radius
  // define commands for this class
  DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

XrayFluoMagneticField::~XrayFluoMagneticField()
{
  delete fMessenger;
}

void XrayFluoMagneticField::GetFieldValue(const G4double point[4],double *bField) const
{
  G4double r = sqrt(point[0]*point[0]+point[1]*point[1]);
  bField[0] = point[1]/(r+1.0e-10*um)*fB0*(r/R)*exp(-pow(r,2)/(R*R))/(sqrt(0.5)*exp(-0.5));
  bField[1] = -point[0]/(r+1.0e-10*um)*fB0*(r/R)*exp(-r*r/(R*R))/(sqrt(0.5)*exp(-0.5));
  bField[2] = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void XrayFluoMagneticField::DefineCommands()
{
  // Define /XrayFluo/field command directory using generic messenger class
  fMessenger = new G4GenericMessenger(this,
                                      "/XrayFluo/field/",
                                      "Field control");

  // fieldValue command
  auto& valueCmd
    = fMessenger->DeclareMethodWithUnit("value","tesla",
                                &XrayFluoMagneticField::SetField,
                                "Set field strength.");
  valueCmd.SetParameterName("field", true);
  valueCmd.SetDefaultValue("1.");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
