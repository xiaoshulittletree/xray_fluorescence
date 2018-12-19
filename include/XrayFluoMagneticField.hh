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
// $Id: B5MagneticField.hh 101036 2016-11-04 09:00:23Z gcosmo $
//
/// \file B5MagneticField.hh
/// \brief Definition of the B5MagneticField class

#ifndef XrayFluoMagneticField_H
#define XrayFluoMagneticField_H 1

#include "globals.hh"
#include "G4MagneticField.hh"

class G4GenericMessenger;

/// Magnetic field

class XrayFluoMagneticField : public G4MagneticField
{
  public:
    XrayFluoMagneticField();
    virtual ~XrayFluoMagneticField();

    virtual void GetFieldValue(const G4double point[4],double* bField ) const;

    void SetField(G4double val) { fB0 = val; }
    G4double GetField() const { return fB0; }

  private:
    void DefineCommands();

    G4GenericMessenger* fMessenger;
    G4double fB0;
    G4double R;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
