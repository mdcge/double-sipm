#include "geometry.hh"

#include "nain4.hh"
#include "g4-mandatory.hh"

#include <CLHEP/Vector/ThreeVector.h>
#include <FTFP_BERT.hh>
#include <G4EmStandardPhysics_option4.hh>
#include <G4LogicalVolume.hh>
#include <G4MaterialPropertiesTable.hh>
#include <G4OpticalPhysics.hh>
#include <G4RotationMatrix.hh>
#include <G4SubtractionSolid.hh>
#include <G4RunManagerFactory.hh>
#include <G4SystemOfUnits.hh>
#include <G4RandomDirection.hh>
#include <G4ThreeVector.hh>
#include <G4Tubs.hh>

G4PVPlacement* geometry() {
    auto fLXe = new G4Material("LXe", 54., 131.29 * g / mole, 3.020 * g / cm3);

    std::vector<G4double> lxe_Energy = {7.0 * eV, 7.07 * eV, 7.14 * eV};

    std::vector<G4double> lxe_SCINT = {0.1, 1.0, 0.1};
    std::vector<G4double> lxe_RIND = {1.59, 1.57, 1.54};
    std::vector<G4double> lxe_ABSL = {35. * cm, 35. * cm, 35. * cm};
    G4MaterialPropertiesTable *fLXe_mt = n4::material_properties()
        .add("SCINTILLATIONCOMPONENT1", lxe_Energy, lxe_SCINT)
        .add("SCINTILLATIONCOMPONENT2", lxe_Energy, lxe_SCINT)
        .add("RINDEX", lxe_Energy, lxe_RIND)
        .add("ABSLENGTH", lxe_Energy, lxe_ABSL)
        .add("SCINTILLATIONYIELD", 12000. / MeV)
        .add("RESOLUTIONSCALE", 1.0)
        .add("SCINTILLATIONTIMECONSTANT1", 20. * ns)
        .add("SCINTILLATIONTIMECONSTANT2", 45. * ns)
        .add("SCINTILLATIONYIELD1", 1.0)
        .add("SCINTILLATIONYIELD2", 0.0)
        .done();
    fLXe->SetMaterialPropertiesTable(fLXe_mt);

    std::vector<G4double> energy = {1.239841939*eV/0.25, 1.239841939*eV/0.9}; // denominator is wavelength in micrometres

    auto csi = n4::material("G4_CESIUM_IODIDE");
    std::vector<G4double> rindex_csi = {2.2094, 1.7611};
    G4MaterialPropertiesTable *mpt_csi = n4::material_properties()
        .add("RINDEX", energy, rindex_csi)
        .add("SCINTILLATIONYIELD", 100000./eV)
        .add("SCINTILLATIONYIELD1", 0.03)
        .done();
    csi -> SetMaterialPropertiesTable(mpt_csi);

    auto air = n4::material("G4_AIR");
    std::vector<G4double> rindex_air = {1.0, 1.0};
    G4MaterialPropertiesTable *mpt_air = n4::material_properties()
        .add("RINDEX", energy, rindex_air)
        .done();
    air -> SetMaterialPropertiesTable(mpt_air);

    auto teflon = n4::material("G4_TEFLON");
    std::vector<G4double> rindex_teflon = {1.305, 1.28}; // taken from "Optical properties of Teflon AF amorphous fluoropolymers" by Yang, French & Tokarsky (using AF2400, Fig.6)
    G4MaterialPropertiesTable *mpt_teflon = n4::material_properties()
        .add("RINDEX", energy, rindex_teflon)
        .done();
    teflon -> SetMaterialPropertiesTable(mpt_teflon);

    // G4double rmin = 0, rmax = 10*cm, half_z = 0.5*cm, min_phi = 0*deg, max_phi = 360*deg;
    G4double half_scint_x = 1.5*mm, half_scint_y = 1.5*mm, half_scint_z = 10*mm;
    G4double half_world_size = 50*mm;
    G4double coating_thck = 0.5*mm;
    // auto cylinder = n4::volume<G4Tubs>("Cylinder", copper, rmin, rmax, half_z, min_phi, max_phi);
    auto scintillator_r = n4::volume<G4Box>("ScintillatorR", fLXe, half_scint_x, half_scint_y, half_scint_z);
    auto scintillator_l = n4::volume<G4Box>("ScintillatorL", csi, half_scint_x, half_scint_y, half_scint_z);

    G4VSolid* coating_ext = new G4Box("CoatingExterior", half_scint_x+coating_thck, half_scint_y+coating_thck, half_scint_z+(coating_thck)/2);
    G4VSolid* coating_int = new G4Box("CoatingInt", half_scint_x, half_scint_y, half_scint_z);
    G4VSolid* coating_solid = new G4SubtractionSolid("Coating", coating_ext, coating_int, 0, G4ThreeVector(0, 0, -coating_thck/2));
    G4LogicalVolume* coating_logical = new G4LogicalVolume(coating_solid, teflon, "Coating", 0, 0, 0);
    auto rot180 = new G4RotationMatrix();
    rot180 -> rotateY(180*deg);

    auto world = n4::volume<G4Box>("World", air, half_world_size, half_world_size, half_world_size);

    // auto cylinder_offset = 1.5*cm;
    auto scintillator_offset = 22.5*mm;
    // n4::place(cylinder).in(world).at({0, 0, cylinder_offset}).now();
    n4::place(scintillator_r).in(world).at({0, 0, scintillator_offset}).check(true).now();
    n4::place(scintillator_l).in(world).at({0, 0, -scintillator_offset}).check(true).now();
    n4::place(coating_logical).in(world).rotate(*rot180).at({0, 0, scintillator_offset-(coating_thck/2)}).check(true).now();
    n4::place(coating_logical).in(world).at({0, 0, -scintillator_offset+(coating_thck/2)}).copy_no(1).check(true).now();
    return n4::place(world).now();
}
