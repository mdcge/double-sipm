#include "nain4.hh"
#include "g4-mandatory.hh"
#include "geometry.hh"
#include "generator.hh"

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
#include <G4Run.hh>
#include <G4Step.hh>

#include <G4UIExecutive.hh>
#include <G4UImanager.hh>
#include <G4VisExecutive.hh>
#include <G4VisManager.hh>
#include <iostream>
#include <memory>


void add_step_edep (G4double& total_edep, G4Step const* step) {
    G4double step_edep;
    G4String step_mat_name = step -> GetPreStepPoint() -> GetMaterial() -> GetName();
    if (step_mat_name == "LXe") {
        G4double pre_energy = step -> GetPreStepPoint() -> GetTotalEnergy();
        G4double post_energy = step -> GetPostStepPoint() -> GetTotalEnergy();
        step_edep = pre_energy - post_energy;
    }
    else {
        step_edep = 0;
    }
    total_edep += step_edep;
}

int main(int argc, char *argv[]) {

    auto two_gammas = [](auto event){ generate_back_to_back_511_keV_gammas(event, {}, 0); };

    G4double total_edep = 0;

    auto reset_total_edep = [&total_edep] (G4Run const* run) { total_edep = 0; };
    auto print_total_edep = [&total_edep] (G4Run const* run) { G4cout << G4endl << "Total deposited energy: " << total_edep << G4endl << G4endl; };
    // This can be used to move the closure "out of scope" (add_step_edep is outside of main)
    auto accumulate_energy = [&total_edep] (G4Step const* step) { add_step_edep(total_edep, step); };
    auto kill_secondaries = [] (G4Track const* track) {
        G4int parent_ID = track -> GetParentID();
        if (parent_ID > 0) {
            return G4ClassificationOfNewTrack::fKill;
        }
        else {
            return G4ClassificationOfNewTrack::fUrgent;
        }
    };

    G4int verbosity = 0;
    auto physics_list = new FTFP_BERT{verbosity};
    physics_list -> ReplacePhysics(new G4EmStandardPhysics_option4());
    physics_list -> RegisterPhysics(new G4OpticalPhysics{});

    auto run_manager = std::unique_ptr<G4RunManager>
        {G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial)};

    // Physics list must be attached to run manager before instantiating other user action classes
    run_manager -> SetUserInitialization(physics_list);
    run_manager -> SetUserInitialization((new n4::actions{two_gammas})
                                         -> set((new n4::run_action())
                                                -> begin(reset_total_edep)
                                                -> end(print_total_edep))
                                         -> set((new n4::stepping_action{accumulate_energy}))
                                         -> set((new n4::stacking_action())
                                                -> classify(kill_secondaries)));
    run_manager -> SetUserInitialization(new n4::geometry{make_geometry});


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
