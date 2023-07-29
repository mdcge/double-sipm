#include "geometry.hh"

#include "nain4.hh"
#include "n4-volumes.hh"
#include "g4-mandatory.hh"

#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Vector/ThreeVector.h>
#include <FTFP_BERT.hh>
#include <G4EmStandardPhysics_option4.hh>
#include <G4LogicalVolume.hh>
#include <G4MaterialPropertiesTable.hh>
#include <G4OpticalPhysics.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4RotationMatrix.hh>
#include <G4SubtractionSolid.hh>
#include <G4RunManagerFactory.hh>
#include <G4SystemOfUnits.hh>
#include <G4RandomDirection.hh>
#include <G4ThreeVector.hh>
#include <G4Tubs.hh>
#include <G4Types.hh>

using vec_double = std::vector<G4double>;

const vec_double OPTPHOT_ENERGY_RANGE{1*eV, 8.21*eV};

G4PVPlacement* make_geometry() {
    auto LXe = new G4Material("LXe", 54., 131.29 * g / mole, 3.020 * g / cm3);

    auto       lxe_energy    = n4::scale_by(eV, { 7.0 ,  7.07,  7.14});
    vec_double lxe_rindex    =                  { 1.59,  1.57,  1.54};
    vec_double lxe_scint     =                  { 0.1 ,  1.0 ,  0.1 };
    auto       lxe_abslength = n4::scale_by(cm, {35   , 35   , 35   });
    G4MaterialPropertiesTable *LXe_properties = n4::material_properties()
        .add("RINDEX"                 , lxe_energy, lxe_rindex)
        .add("SCINTILLATIONCOMPONENT1", lxe_energy, lxe_scint)
        .add("SCINTILLATIONCOMPONENT2", lxe_energy, lxe_scint)
        .add("ABSLENGTH"              , lxe_energy, lxe_abslength)
        .add("SCINTILLATIONTIMECONSTANT1",    20 * ns )
        .add("SCINTILLATIONTIMECONSTANT2",    45 * ns )
        .add("SCINTILLATIONYIELD"        , 12000 / MeV)
        .add("SCINTILLATIONYIELD1"       ,     1.0    )
        .add("SCINTILLATIONYIELD2"       ,     0.0    )
        .add("RESOLUTIONSCALE"           ,     1.0    )
        .done();
    LXe -> SetMaterialPropertiesTable(LXe_properties);

    auto csi = n4::material("G4_CESIUM_IODIDE");

    // csi_rindex: values taken from "Optimization of Parameters for a CsI(Tl) Scintillator Detector Using GEANT4-Based Monte Carlo..." by Mitra et al (mainly page 3)
    //  csi_scint: Fig. 2 in the paper
    auto hc = 1.239841939;
    auto csi_energy       = n4::scale_by(hc*eV, {1/0.35, 1/0.54, 1/0.7, 1/0.9}) ; // denominator is wavelength in micrometres
    vec_double csi_rindex =                     {1.79  , 1.79  , 1.79 , 1.79 };   //vec_double csi_rindex = {2.2094, 1.7611};
    vec_double  csi_scint =                     {0.0   , 1.0   , 0.1  , 0.0  };
    auto    csi_abslength = n4::scale_by(m    , {5     , 5     , 5    , 5    });
    G4MaterialPropertiesTable *csi_mpt = n4::material_properties()
        .add("RINDEX"                 , csi_energy, csi_rindex)
        .add("SCINTILLATIONCOMPONENT1", csi_energy,  csi_scint)
        .add("SCINTILLATIONCOMPONENT2", csi_energy,  csi_scint)
        .add("ABSLENGTH"              , csi_energy, csi_abslength)
        .add("SCINTILLATIONTIMECONSTANT1",   700 * ns )
        .add("SCINTILLATIONTIMECONSTANT2",  3500 * ns )
        .add("SCINTILLATIONYIELD"        , 65000 / MeV)
        .add("SCINTILLATIONYIELD1"       ,     0.57   )
        .add("SCINTILLATIONYIELD2"       ,     0.43   )
        .add("RESOLUTIONSCALE"           ,     1.0    )
        .done();
    csi -> GetIonisation() -> SetBirksConstant(0.00152 * mm/MeV);
    csi -> SetMaterialPropertiesTable(csi_mpt);

    auto air = n4::material("G4_AIR");
    G4MaterialPropertiesTable *mpt_air = n4::material_properties()
        .add("RINDEX", OPTPHOT_ENERGY_RANGE, {1, 1})
        .done();
    air -> SetMaterialPropertiesTable(mpt_air);

    auto teflon = n4::material("G4_TEFLON");
    // Values could be taken from "Optical properties of Teflon AF amorphous fluoropolymers" by Yang, French & Tokarsky (using AF2400, Fig.6)
    // but are also stated in the same paper as above
    G4MaterialPropertiesTable *mpt_teflon = n4::material_properties()
        .add("RINDEX", OPTPHOT_ENERGY_RANGE, {1.35, 1.35})
        .done();
    teflon -> SetMaterialPropertiesTable(mpt_teflon);

    G4double scint_xy = 3*mm, scint_z = 2*cm;
    G4double world_size = 10*cm;
    G4double coating_thck = 0.5*mm;
    // auto cylinder = n4::tubs("Cylinder").r(10*cm).z(1*cm).volume(copper);
    auto coating_interior = n4::box("CoatingInterior").cube(scint_xy                 ).z(scint_z                );
    auto coating_logical  = n4::box("CoatingExterior").cube(scint_xy + coating_thck*2).z(scint_z  + coating_thck)
        .subtract(coating_interior)
        .at(0, 0, -coating_thck/2)
        .name("Coating")
        .volume(teflon);

    auto rot180 = new G4RotationMatrix(); rot180 -> rotateY(180*deg);

    auto world = n4::box{"World"}.cube(world_size).volume(air);

    // auto cylinder_offset = 1.5*cm;
    // n4::place(cylinder).in(world).at({0, 0, cylinder_offset}).now();
    auto scintillator_offset = 22.5*mm;
    auto scintillator = coating_interior.name("Scintillator").volume(csi);
    n4::place(scintillator)   .in(world)                .at({0, 0,  scintillator_offset                   }).copy_no(0).now();
    n4::place(scintillator)   .in(world)                .at({0, 0, -scintillator_offset                   }).copy_no(1).now();
    n4::place(coating_logical).in(world).rotate(*rot180).at({0, 0,  scintillator_offset - (coating_thck/2)}).copy_no(0).now();
    n4::place(coating_logical).in(world)                .at({0, 0, -scintillator_offset + (coating_thck/2)}).copy_no(1).now();

    // Check which world's daughter is which object
    //G4cout << "XXXXXXXXXXXXXXXX " << world->GetDaughter(2)->GetName() << G4endl;

    // TO BE CHANGED!!!!!!
    G4OpticalSurface* OpSurface = new G4OpticalSurface("name");

    // world's 2nd daughter is the right teflon coating, world's 0th daughter is the right scintillator
    G4LogicalBorderSurface* Surface = new G4LogicalBorderSurface("name", world->GetDaughter(2), world->GetDaughter(0), OpSurface);

    OpSurface->SetType(dielectric_dielectric);
    OpSurface->SetModel(unified);
    OpSurface->SetFinish(groundbackpainted);
    OpSurface->SetSigmaAlpha(0.1);

    vec_double pp = {2.038*eV, 4.144*eV};
    OpSurface -> SetMaterialPropertiesTable(
        n4::material_properties{}
           .add("RINDEX"                , pp, {0.3 , 0.3 })
           .add("SPECULARLOBECONSTANT"  , pp, {0.2 , 0.2 })
           .add("SPECULARSPIKECONSTANT" , pp, {0.1 , 0.1 })
           .add("BACKSCATTERCONSTANT"   , pp, {1.35, 1.40})
           .add("REFLECTIVITY"          , pp, {0.3 , 0.5 })
           .add("EFFICIENCY"            , pp, {0.8 , 0.1 })
           .done()
    );

    return n4::place(world).now();
}
