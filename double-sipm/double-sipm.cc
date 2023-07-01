#include "nain4.hh"
#include "g4-mandatory.hh"

#include <FTFP_BERT.hh>
#include <G4EmStandardPhysics_option4.hh>
#include <G4OpticalPhysics.hh>
#include <G4RunManagerFactory.hh>
#include <G4SystemOfUnits.hh>
#include <G4RandomDirection.hh>
#include <G4Tubs.hh>

#include <G4UIExecutive.hh>
#include <G4UImanager.hh>
#include <G4VisExecutive.hh>
#include <G4VisManager.hh>
#include <iostream>
#include <memory>

void generate_back_to_back_511_keV_gammas(G4Event* event, G4ThreeVector position, G4double time) {
    auto gamma = nain4::find_particle("gamma");
    auto p = 511*keV * G4RandomDirection();
    auto vertex =      new G4PrimaryVertex(position, time);
    vertex -> SetPrimary(new G4PrimaryParticle(gamma,  p.x(),  p.y(),  p.z()));
    vertex -> SetPrimary(new G4PrimaryParticle(gamma, -p.x(), -p.y(), -p.z()));
    event -> AddPrimaryVertex(vertex);
}

int main(int argc, char *argv[]) {
    std::cout << "Hello World!" << std::endl;

    auto geometry = [] () {
        auto copper = n4::material("G4_Cu");
        auto csi = n4::material("G4_CESIUM_IODIDE");
        auto air = n4::material("G4_AIR");

        G4double rmin = 0, rmax = 10*cm, half_z = 0.5*cm, min_phi = 0*deg, max_phi = 360*deg;
        G4double scint_x = 2*cm, scint_y = 1*cm, scint_z = 1*cm;
        G4double world_size = 20*cm;
        auto cylinder = n4::volume<G4Tubs>("Cylinder", copper, rmin, rmax, half_z, min_phi, max_phi);
        auto scintillator_r = n4::volume<G4Box>("ScintillatorR", csi, scint_x, scint_y, scint_z);
        auto scintillator_l = n4::volume<G4Box>("ScintillatorL", csi, scint_x, scint_y, scint_z);
        auto world = n4::volume<G4Box>("World", air, world_size, world_size, world_size);

        auto cylinder_offset = 1.5*cm;
        auto scintillator_offset = 9*cm;
        n4::place(cylinder).in(world).at({0, 0, cylinder_offset}).now();
        n4::place(scintillator_r).in(world).at({scintillator_offset, 0, 0}).now();
        n4::place(scintillator_l).in(world).at({-scintillator_offset, 0, 0}).now();
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
