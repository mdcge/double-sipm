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

#include <G4UIExecutive.hh>
#include <G4UImanager.hh>
#include <G4VisExecutive.hh>
#include <G4VisManager.hh>
#include <iostream>
#include <memory>

void generate_back_to_back_511_keV_gammas(G4Event* event, G4ThreeVector position, G4double time) {
    auto gamma = nain4::find_particle("gamma");
    auto direction = G4ThreeVector(3*G4UniformRand()-1.5, 3*G4UniformRand()-1.5, 12.5).unit(); // random unit vector which hits scintillator
    auto p = 0.511*MeV * direction;
    auto vertex = new G4PrimaryVertex(position, time);
    vertex -> SetPrimary(new G4PrimaryParticle(gamma,  p.x(),  p.y(),  p.z()));
    vertex -> SetPrimary(new G4PrimaryParticle(gamma, -p.x(), -p.y(), -p.z()));
    event -> AddPrimaryVertex(vertex);
}

int main(int argc, char *argv[]) {

    auto geometry = [] () {
        std::vector<G4double> energy = {1.239841939*eV/0.25, 1.239841939*eV/0.9};

        auto csi = n4::material("G4_CESIUM_IODIDE");
        std::vector<G4double> rindex_csi = {2.2094, 1.7611};
        G4MaterialPropertiesTable *mpt_csi = n4::material_properties()
            .add("RINDEX", energy, rindex_csi)
            .add("SCINTILLATIONYIELD", 100000./eV).done();
        csi -> SetMaterialPropertiesTable(mpt_csi);

        auto air = n4::material("G4_AIR");
        std::vector<G4double> rindex_air = {1.0, 1.0};
        G4MaterialPropertiesTable *mpt_air = n4::material_properties()
            .add("RINDEX", energy, rindex_air)
            .done();
        air -> SetMaterialPropertiesTable(mpt_air);

        auto teflon = n4::material("G4_TEFLON");

        // G4double rmin = 0, rmax = 10*cm, half_z = 0.5*cm, min_phi = 0*deg, max_phi = 360*deg;
        G4double half_scint_x = 1.5*mm, half_scint_y = 1.5*mm, half_scint_z = 10*mm;
        G4double half_world_size = 50*mm;
        G4double coating_thck = 0.5*mm;
        // auto cylinder = n4::volume<G4Tubs>("Cylinder", copper, rmin, rmax, half_z, min_phi, max_phi);
        auto scintillator_r = n4::volume<G4Box>("ScintillatorR", csi, half_scint_x, half_scint_y, half_scint_z);
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
    };

    auto two_gammas = [](auto event){ generate_back_to_back_511_keV_gammas(event, {}, 0); };

    G4int verbosity = 0;
    auto physics_list = new FTFP_BERT{verbosity};
    physics_list -> ReplacePhysics(new G4EmStandardPhysics_option4());
    physics_list -> RegisterPhysics(new G4OpticalPhysics{});

    auto run_manager = std::unique_ptr<G4RunManager>
        {G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial)};

    // Physics list must be attached to run manager before instantiating other user action classes
    run_manager -> SetUserInitialization(physics_list);
    run_manager -> SetUserInitialization(new n4::actions{two_gammas});
    run_manager -> SetUserInitialization(new n4::geometry{geometry});


    // Initialize visualization
    std::unique_ptr<G4VisManager> visManager = std::make_unique<G4VisExecutive>();
    // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
    // G4VisManager* visManager = new G4VisExecutive("Quiet");
    visManager -> Initialize();

    // Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    std::unique_ptr<G4UIExecutive> ui;
    if ( argc == 1 ) { ui = std::make_unique<G4UIExecutive>(argc, argv); }
    if ( ! ui ) {
        // batch mode
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command+fileName);
    } else {
        // interactive mode
        UImanager->ApplyCommand("/control/execute init_vis.mac");
        ui->SessionStart();
    }

}
