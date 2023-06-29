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
/// \file DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "nain4.hh"

namespace B1 {

DetectorConstruction:: DetectorConstruction(){}
DetectorConstruction::~DetectorConstruction(){}

G4VPhysicalVolume* DetectorConstruction::Construct() {

  // Materials
  auto water  = n4::material("G4_WATER");
  auto air    = n4::material("G4_AIR");
  auto tissue = n4::material("G4_A-150_TISSUE");
  auto bone   = n4::material("G4_BONE_COMPACT_ICRU");

  // Dimensions

  // Envelope parameters
  G4double env_sizeXY = 20*cm, env_sizeZ = 30*cm;
  // World
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  // Conical section shape
  G4double cone_rmina =  0.*cm, cone_rmaxa = 2.*cm;
  G4double cone_rminb =  0.*cm, cone_rmaxb = 4.*cm;
  G4double cone_hz = 3.*cm;
  G4double cone_phimin = 0.*deg, cone_phimax = 360.*deg;
  // Trapezoid shape
  G4double trapezoid_dxa = 12*cm, trapezoid_dxb = 12*cm;
  G4double trapezoid_dya = 10*cm, trapezoid_dyb = 16*cm;
  G4double trapezoid_dz  = 6*cm;

  // Volumes
  auto world     = n4::volume<G4Box>("World", air, 0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);
  auto envelope  = n4::volume<G4Box>("Envelope", water, 0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ);
  auto cone      = n4::volume<G4Cons>("Tissue Cone", tissue, cone_rmina, cone_rmaxa, cone_rminb, cone_rmaxb, cone_hz, cone_phimin, cone_phimax);
  auto trapezoid = n4::volume<G4Trd>("Bone Trapezoid", bone, 0.5*trapezoid_dxa, 0.5*trapezoid_dxb, 0.5*trapezoid_dya, 0.5*trapezoid_dyb, 0.5*trapezoid_dz);

  // Set Trapezoid as scoring volume
  this->fScoringVolume = trapezoid;

  // Placement
  n4::       place(envelope) .in(world)                      .check(true).now();
  n4::       place(cone)     .in(envelope).at(0, 2*cm, -7*cm).check(true).now();
  n4::       place(trapezoid).in(envelope).at(0, -1*cm, 7*cm).check(true).now();
  return n4::place(world)                                    .check(true).now();

}

} // namespace::B1
