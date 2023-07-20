#include "geometry.hh"

#include "nain4.hh"
#include "g4-mandatory.hh"

#include <CLHEP/Vector/ThreeVector.h>
#include <FTFP_BERT.hh>
#include <G4EmStandardPhysics_option4.hh>
#include <G4HCofThisEvent.hh>
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
#include <Randomize.hh>


G4int photon_count_0 = 0;
G4int photon_count_1 = 0;
std::vector<G4double> times_of_arrival_0;
std::vector<G4double> times_of_arrival_1;


G4double detection_probability(G4double energy, std::vector<G4double> energies, std::vector<G4double> scintillation);

G4PVPlacement* make_geometry() {

    // must be in increasing ENERGY order (decreasing wavelength) for scintillation to work properly
    std::vector<G4double> energies = {1.239841939*eV/0.9, 1.239841939*eV/0.7, 1.239841939*eV/0.54, 1.239841939*eV/0.35}; // denominator is wavelength in micrometres

    auto csi = n4::material("G4_CESIUM_IODIDE");
    //std::vector<G4double> rindex_csi = {2.2094, 1.7611};
    std::vector<G4double> rindex_csi = {1.79, 1.79, 1.79, 1.79}; // values taken from "Optimization of Parameters for a CsI(Tl) Scintillator Detector Using GEANT4-Based Monte Carlo..." by Mitra et al (mainly page 3)
    std::vector<G4double> scint_csi = {0.0, 0.1, 1.0, 0.0}; // Fig. 2 in the paper
    // Values from "Temperature dependence of pure CsI: scintillation light yield and decay time" by Amsler et al
    // "cold" refers to ~77K, i.e. liquid nitrogen temperature
    G4double csi_scint_yield = 3200. / MeV;
    G4double csi_scint_yield_cold = 50000. / MeV;
    G4double csi_time_fast = 6 * ns;
    G4double csi_time_slow = 28 * ns;
    G4double csi_time_fast_cold = 1015 * ns; // only one component at cold temps!
    G4double csi_time_slow_cold = 1015 * ns;
    G4MaterialPropertiesTable *mpt_csi = n4::material_properties()
        .add("RINDEX", energies, rindex_csi)
        .add("SCINTILLATIONYIELD", csi_scint_yield)
        // .add("SCINTILLATIONYIELD", 100. / MeV) // for testing
        .add("SCINTILLATIONTIMECONSTANT1", csi_time_fast)
        .add("SCINTILLATIONTIMECONSTANT2", csi_time_slow)
        .add("RESOLUTIONSCALE", 1.0)
        .add("SCINTILLATIONYIELD1", 0.57)
        .add("SCINTILLATIONYIELD2", 0.43)
        .add("SCINTILLATIONCOMPONENT1", energies, scint_csi)
        .add("SCINTILLATIONCOMPONENT2", energies, scint_csi)
        .add("ABSLENGTH", energies, {5.*m, 5.*m, 5.*m, 5.*m})
        .done();
    csi -> GetIonisation() -> SetBirksConstant(0.00152 * mm/MeV);
    csi -> SetMaterialPropertiesTable(mpt_csi);

    auto air = n4::material("G4_AIR");
    std::vector<G4double> rindex_air = {1.0, 1.0, 1.0, 1.0};
    G4MaterialPropertiesTable *mpt_air = n4::material_properties()
        .add("RINDEX", energies, rindex_air)
        .done();
    air -> SetMaterialPropertiesTable(mpt_air);

    auto teflon = n4::material("G4_TEFLON");
    // Values could be taken from "Optical properties of Teflon AF amorphous fluoropolymers" by Yang, French & Tokarsky (using AF2400, Fig.6)
    // but are also stated in the same paper as above
    std::vector<G4double> rindex_teflon = {1.35, 1.35, 1.35, 1.35};
    G4MaterialPropertiesTable *mpt_teflon = n4::material_properties()
        .add("RINDEX", energies, rindex_teflon)
        .done();
    teflon -> SetMaterialPropertiesTable(mpt_teflon);

    auto plastic = n4::material("G4_POLYCARBONATE"); // probably wrong

    G4double half_scint_x = 1.5*mm, half_scint_y = 1.5*mm, half_scint_z = 10*mm;
    auto scintillator = n4::volume<G4Box>("Scintillator", csi, half_scint_x, half_scint_y, half_scint_z);

    G4double coating_thck = 0.25*mm;
    G4VSolid* coating_ext = new G4Box("CoatingExterior", half_scint_x+coating_thck, half_scint_y+coating_thck, half_scint_z+(coating_thck)/2);
    G4VSolid* coating_int = new G4Box("CoatingInterior", half_scint_x, half_scint_y, half_scint_z);
    G4VSolid* coating_solid = new G4SubtractionSolid("Coating", coating_ext, coating_int, 0, G4ThreeVector(0, 0, -coating_thck/2));
    G4LogicalVolume* coating_logical = new G4LogicalVolume(coating_solid, teflon, "Coating", 0, 0, 0);
    auto rotY180 = new G4RotationMatrix();
    rotY180 -> rotateY(180*deg);

    G4int nb_detectors_per_side = 3;
    G4double half_detector_width = half_scint_x / nb_detectors_per_side; // assumes the detectors are square
    G4double half_detector_depth = half_detector_width; // this will make the detectors cubes
    auto detector = n4::volume<G4Box>("Detector", air, half_detector_width, half_detector_width, half_detector_depth); // material doesn't matter

    G4double half_world_size = 50*mm;
    auto world = n4::volume<G4Box>("World", air, half_world_size, half_world_size, half_world_size);

    G4double source_ring_rmax = 12.5*mm; G4double source_ring_rmin = 9.5*mm; G4double source_ring_thck = 3*mm;
    auto source_ring = n4::volume<G4Tubs>("SourceRing", plastic, source_ring_rmin, source_ring_rmax, source_ring_thck, 0*deg, 360*deg);
    auto rotY90 = new G4RotationMatrix();
    rotY90 -> rotateY(90*deg);

    auto scintillator_offset = 23*mm;

    n4::place(scintillator).in(world).at({0, 0,  scintillator_offset}).copy_no(0).now();
    n4::place(scintillator).in(world).at({0, 0, -scintillator_offset}).copy_no(1).now();
    n4::place(coating_logical).in(world).rotate(*rotY180).at({0, 0, scintillator_offset-(coating_thck/2)}).copy_no(0).now();
    n4::place(coating_logical).in(world).at({0, 0, -scintillator_offset+(coating_thck/2)}).copy_no(1).now();
    n4::place(source_ring).in(world).rotate(*rotY90).at({0, 0, 0}).now();

    for (G4int side=0; side<2; side++) {
        for (G4int i=0; i<nb_detectors_per_side; i++) {
            for (G4int j=0; j<nb_detectors_per_side; j++) {
                G4double xpos = (i - ((float) nb_detectors_per_side/2 - 0.5)) * 2*half_detector_width;
                G4double ypos = (j - ((float) nb_detectors_per_side/2 - 0.5)) * 2*half_detector_width;
                G4double zpos;
                if (side == 0) { zpos = scintillator_offset+half_scint_z+half_detector_depth; }
                else { zpos = -(scintillator_offset+half_scint_z+half_detector_depth); }
                n4::place(detector)
                    .in(world)
                    .at({xpos, ypos, zpos})
                    .copy_no(side*pow(nb_detectors_per_side, 2) + i*nb_detectors_per_side + j)
                    .now();
            }
        }
    }

    // Detector physics -------------------------------------

    auto process_hits = [nb_detectors_per_side](G4Step* step) {
        G4Track* track = step -> GetTrack();
        track -> SetTrackStatus(fStopAndKill);

        G4int copy_nb = step -> GetPreStepPoint() -> GetTouchable() -> GetCopyNumber();
        G4int time = step -> GetPreStepPoint() -> GetGlobalTime();
        G4ThreeVector photon_momentum = step -> GetPreStepPoint() -> GetMomentum();
        G4double photon_energy = photon_momentum.mag();

        // Pixel pitch 25 um
        std::vector<G4double> sipm_energies = {1.239841939*eV/0.9, 1.239841939*eV/0.7, 1.239841939*eV/0.5, 1.239841939*eV/0.46, 1.239841939*eV/0.4, 1.239841939*eV/0.32};
        std::vector<G4double> sipm_pdes = {0.03, 0.1, 0.245, 0.255, 0.23, 0.02};
        if (G4UniformRand() < detection_probability(photon_energy, sipm_energies, sipm_pdes)) {
            if (copy_nb < pow(nb_detectors_per_side, 2)) {
                times_of_arrival_0.push_back(time);
                photon_count_0 += 1;
            }
            else {
                times_of_arrival_1.push_back(time);
                photon_count_1 += 1;
            }
        }

        return true;
    };

    auto end_of_event = [](G4HCofThisEvent* what) {};
    auto sens_detector = new n4::sensitive_detector("Detector", process_hits, end_of_event);
    detector -> SetSensitiveDetector(sens_detector);

    // Border surface ------------------------------------------
    // Check which world's daughter is which object
    //G4cout << "XXXXXXXXXXXXXXXX " << world->GetDaughter(2)->GetName() << G4endl;

    G4OpticalSurface* csi_teflon_surface = new G4OpticalSurface("CsI-TeflonSurface");

    // Values from same paper as above ("Optimization of Parameters...")
    // "groundfrontpainted" (I think) only considers whether the photon is reflected or absorbed, so there will be no shine through visible in the simulation
    csi_teflon_surface->SetType(dielectric_dielectric);
    csi_teflon_surface->SetModel(unified);
    csi_teflon_surface->SetFinish(groundfrontpainted);
    csi_teflon_surface->SetSigmaAlpha(0.0);

    // world's 2nd daughter is the right teflon coating, world's 0th daughter is the right scintillator
    // this seems to apply the surface to the two physical objects without needing assignment
    G4LogicalBorderSurface* border0 = new G4LogicalBorderSurface("CsI-TeflonSurface", world->GetDaughter(0), world->GetDaughter(2), csi_teflon_surface);
    G4LogicalBorderSurface* border1 = new G4LogicalBorderSurface("CsI-TeflonSurface", world->GetDaughter(1), world->GetDaughter(3), csi_teflon_surface);

    std::vector<G4double> reflectivity = {0.98, 0.98, 0.98, 0.98};

    // According to the docs, for UNIFIED, dielectric_dielectric surfaces only the Lambertian reflection is turned on
    G4MaterialPropertiesTable* csi_teflon_surface_mpt = n4::material_properties()
        .add("REFLECTIVITY", energies, reflectivity)
        .done();
    csi_teflon_surface -> SetMaterialPropertiesTable(csi_teflon_surface_mpt);

    return n4::place(world).now();
}


G4double detection_probability(G4double energy, std::vector<G4double> energies, std::vector<G4double> scintillation) {
    G4int index;
    G4int length = (G4int) energies.size() - 1;
    if (energy < energies[0] || energies[length] < energy) {
        return 0;
    }
    for (G4int i=0; i<length; i++) {
        if (energies[i] <= energy && energy <= energies[i+1]) {
            index = i;
        }
    }
    G4double y0 = scintillation[index]; G4double y1 = scintillation[index+1];
    G4double x0 = energies[index]; G4double x1 = energies[index+1];
    return y0 + ((y1 - y0) / (x1 - x0)) * (energy - x0);
}
