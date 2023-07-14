#include "geometry.hh"

#include "nain4.hh"
#include "n4-volumes.hh"
#include "g4-mandatory.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Vector/ThreeVector.h>
#include <FTFP_BERT.hh>
#include <G4EmStandardPhysics_option4.hh>
#include <G4HCofThisEvent.hh>
#include <G4LogicalVolume.hh>
#include <G4MaterialPropertiesTable.hh>
#include <G4OpticalPhysics.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4SubtractionSolid.hh>
#include <G4RunManagerFactory.hh>
#include <G4SystemOfUnits.hh>
#include <G4RandomDirection.hh>
#include <G4ThreeVector.hh>
#include <G4Tubs.hh>
#include <G4Types.hh>

using vec_double = std::vector<G4double>;
using vec_int    = std::vector<G4int>;

const vec_double OPTPHOT_ENERGY_RANGE{1*eV, 8.21*eV};

G4PVPlacement* make_geometry(vec_int& photon_count) {
    auto csi = n4::material("G4_CESIUM_IODIDE");

    // csi_rindex: values taken from "Optimization of Parameters for a CsI(Tl) Scintillator Detector Using GEANT4-Based Monte Carlo..." by Mitra et al (mainly page 3)
    //  csi_scint: Fig. 2 in the paper
    auto hc = CLHEP::h_Planck * CLHEP::c_light;
    auto csi_energy       = n4::scale_by(hc*eV, {1/0.35, 1/0.54, 1/0.7, 1/0.9}) ; // denominator is wavelength in micrometres
    vec_double csi_rindex =                     {1.79  , 1.79  , 1.79 , 1.79 };   //vec_double csi_rindex = {2.2094, 1.7611};
    vec_double  csi_scint =                     {0.0   , 1.0   , 0.1  , 0.0  };
    auto    csi_abslength = n4::scale_by(m    , {5     , 5     , 5    , 5    });
    // Values from "Temperature dependence of pure CsI: scintillation light yield and decay time" by Amsler et al
    // "cold" refers to ~77K, i.e. liquid nitrogen temperature
    G4double csi_scint_yield      =  3200 / MeV;
    G4double csi_scint_yield_cold = 50000 / MeV;
    G4double csi_time_fast        =     6 * ns;
    G4double csi_time_slow        =    28 * ns;
    G4double csi_time_fast_cold   =  1015 * ns; // only one component at cold temps!
    G4double csi_time_slow_cold   =  1015 * ns;
    G4MaterialPropertiesTable *csi_mpt = n4::material_properties()
        .add("RINDEX"                 , csi_energy, csi_rindex)
        .add("SCINTILLATIONCOMPONENT1", csi_energy,  csi_scint)
        .add("SCINTILLATIONCOMPONENT2", csi_energy,  csi_scint)
        .add("ABSLENGTH"              , csi_energy, csi_abslength)
        .add("SCINTILLATIONTIMECONSTANT1", csi_time_fast)
        .add("SCINTILLATIONTIMECONSTANT2", csi_time_slow)
        .add("SCINTILLATIONYIELD"        , csi_scint_yield)
      //.add("SCINTILLATIONYIELD"        ,   100 / MeV) // for testing
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

    auto plastic = n4::material("G4_POLYCARBONATE"); // probably wrong

    G4double scint_xy = 3*mm, scint_z = 2*cm;
    G4double world_size = 10*cm;
    G4double coating_thck = 0.25*mm;
    // auto cylinder = n4::tubs("Cylinder").r(10*cm).z(1*cm).volume(copper);
    auto coating_interior = n4::box("CoatingInterior").cube(scint_xy                 ).z(scint_z                );
    auto coating_logical  = n4::box("CoatingExterior").cube(scint_xy + coating_thck*2).z(scint_z  + coating_thck)
        .subtract(coating_interior)
        .at(0, 0, -coating_thck/2)
        .name("Coating")
        .volume(teflon);

    G4int nb_detectors_per_side = 3;
    G4double half_detector_width = scint_xy / 2 / nb_detectors_per_side; // assumes the detectors are square
    G4double half_detector_depth = half_detector_width; // this will make the detectors cubes
    auto detector = n4::volume<G4Box>("Detector", air, half_detector_width, half_detector_width, half_detector_depth); // material doesn't matter

    auto world = n4::box{"World"}.cube(world_size).volume(air);

    auto scintillator_offset = 23*mm;
    auto scintillator = coating_interior.name("Scintillator").volume(csi);
    n4::place(scintillator)   .in(world)                  .at({0, 0,  scintillator_offset                   }).copy_no(0).now();
    n4::place(scintillator)   .in(world)                  .at({0, 0, -scintillator_offset                   }).copy_no(1).now();
    n4::place(coating_logical).in(world).rotate_y(180*deg).at({0, 0,  scintillator_offset - (coating_thck/2)}).copy_no(0).now();
    n4::place(coating_logical).in(world)                  .at({0, 0, -scintillator_offset + (coating_thck/2)}).copy_no(1).now();

    G4double source_ring_rmax = 12.5*mm; G4double source_ring_rmin = 9.5*mm; G4double source_ring_thck = 3*mm;
    auto source_ring = n4::volume<G4Tubs>("SourceRing", plastic, source_ring_rmin, source_ring_rmax, source_ring_thck, 0*deg, 360*deg);

    n4::place(source_ring).in(world).rotate_y(90*deg).at({0, 0, 0}).now();

    for (G4int side=0; side<2; side++) {
        for (G4int i=0; i<nb_detectors_per_side; i++) {
            for (G4int j=0; j<nb_detectors_per_side; j++) {
                G4double xpos = (i - ((float) nb_detectors_per_side/2 - 0.5)) * 2*half_detector_width;
                G4double ypos = (j - ((float) nb_detectors_per_side/2 - 0.5)) * 2*half_detector_width;
                G4double zpos;
                if (side == 0) { zpos =   scintillator_offset + scint_z/2 + half_detector_depth;  }
                else           { zpos = -(scintillator_offset + scint_z/2 + half_detector_depth); }
                n4::place(detector)
                    .in(world)
                    .at({xpos, ypos, zpos})
                    .copy_no(side*pow(nb_detectors_per_side, 2) + i*nb_detectors_per_side + j)
                    .now();
            }
        }
    }

    auto process_hits = [nb_detectors_per_side, &photon_count](G4Step* step) {
        G4Track* track = step -> GetTrack();
        track -> SetTrackStatus(fStopAndKill);

        G4int copy_nb = step -> GetPreStepPoint() -> GetTouchable() -> GetCopyNumber();
        if (copy_nb < pow(nb_detectors_per_side, 2)) { photon_count[0]++; }
        else                                         { photon_count[1]++; }

        G4int time = step -> GetPreStepPoint() -> GetGlobalTime();
        // G4cout << "XXXXXXXXXXXXXXXXXXXX " << time << G4endl;
        return true;
    };
    auto end_of_event = [](G4HCofThisEvent* what) {};
    auto sens_detector = (new n4::sensitive_detector{"Detector", process_hits}) -> end_of_event(end_of_event);
    detector -> SetSensitiveDetector(sens_detector);

    // Border surface ------------------------------------------
    // Check which world's daughter is which object
    //G4cout << "XXXXXXXXXXXXXXXX " << world->GetDaughter(2)->GetName() << G4endl;

    G4OpticalSurface* surface = new G4OpticalSurface("OpticalSurface");

    // Values from same paper as above ("Optimization of Parameters...")
    // "groundfrontpainted" (I think) only considers whether the photon is reflected or absorbed, so there will be no shine through visible in the simulation
    surface->SetType(dielectric_dielectric);
    surface->SetModel(unified);
    surface->SetFinish(groundfrontpainted);
    surface->SetSigmaAlpha(0.0);

    // world's 2nd daughter is the right teflon coating, world's 0th daughter is the right scintillator
    // this seems to apply the surface to the two physical objects without needing assignment
    G4LogicalBorderSurface* border0 = new G4LogicalBorderSurface("OpticalSurface", world->GetDaughter(0), world->GetDaughter(2), surface);
    G4LogicalBorderSurface* border1 = new G4LogicalBorderSurface("OpticalSurface", world->GetDaughter(1), world->GetDaughter(3), surface);

    vec_double pp = {2.038*eV, 4.144*eV};
    // According to the docs, for UNIFIED, dielectric_dielectric surfaces only the Lambertian reflection is turned on
    surface -> SetMaterialPropertiesTable(
        n4::material_properties{}
           .add("REFLECTIVITY", pp, {0.98 , 0.98 })
           .done()
    );

    return n4::place(world).now();
}
