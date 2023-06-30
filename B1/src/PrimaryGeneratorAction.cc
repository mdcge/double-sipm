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
//
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the B1::PrimaryGeneratorAction class

#include "nain4.hh"

#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <memory>

namespace B1
{

G4Box* getEnvelope();

PrimaryGeneratorAction::PrimaryGeneratorAction()
  : fParticleGun{std::make_unique<G4ParticleGun>(1)}
{
  fParticleGun -> SetParticleDefinition(n4::find_particle("gamma"));
  fParticleGun -> SetParticleMomentumDirection({0, 0, 1});
  fParticleGun -> SetParticleEnergy(6 * MeV);
}

// This function is called at the beginning of each event
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

  if ( !fEnvelopeBox ) { fEnvelopeBox = getEnvelope(); }

  G4double envSizeXY = fEnvelopeBox -> GetXHalfLength() * 2;
  G4double envSizeZ  = fEnvelopeBox -> GetZHalfLength() * 2;

  G4double size = 0.8;
  G4double x0 = size * envSizeXY * (G4UniformRand() - 0.5);
  G4double y0 = size * envSizeXY * (G4UniformRand() - 0.5);
  G4double z0 = -0.5 * envSizeZ;

  fParticleGun -> SetParticlePosition({x0, y0, z0});
  fParticleGun -> GeneratePrimaryVertex(anEvent);
}

// In order to avoid dependence of PrimaryGeneratorAction
// on DetectorConstruction class we get Envelope volume
// from G4LogicalVolumeStore.
G4Box* getEnvelope() {
  auto logical_volume = n4::find_logical("Envelope");
  auto solid_box = dynamic_cast<G4Box*>(logical_volume -> GetSolid());
  if (!solid_box) {
    G4ExceptionDescription msg;
    msg << "Envelope volume of box shape not found.\n";
    msg << "Perhaps you have changed geometry.\n";
    msg << "The gun will be placed at the centre.";
    G4Exception("PrimaryGeneratorAction::GeneratePrimaries()",
     "MyCode0002",JustWarning,msg);
  }
  return solid_box;
}

}


