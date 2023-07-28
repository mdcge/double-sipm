#include "n4_run_manager.hh"
#include "nain4.hh"
#include "n4_ui.hh"
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

#include <G4Types.hh>
#include <G4UIExecutive.hh>
#include <G4UImanager.hh>
#include <G4VisExecutive.hh>
#include <G4VisManager.hh>
#include <iostream>
#include <memory>

auto physics_list() {
    G4int verbosity;
    auto physics_list = new FTFP_BERT{verbosity = 0};
    physics_list ->  ReplacePhysics(new G4EmStandardPhysics_option4());
    physics_list -> RegisterPhysics(new G4OpticalPhysics{});
    return physics_list;
}

G4double delta_total_energy(G4Step const * step) {
    auto  pre_energy = step ->  GetPreStepPoint() -> GetTotalEnergy();
    auto post_energy = step -> GetPostStepPoint() -> GetTotalEnergy();
    return pre_energy - post_energy;
}

// Main part of run action
void add_step_edep(G4double& total_edep_0, G4double& total_edep_1, G4Step const* step) {
    G4double step_edep_0 = 0;
    G4double step_edep_1 = 0;

    auto step_solid_name = step -> GetPreStepPoint() -> GetTouchable() -> GetVolume() -> GetName();

    if      (step_solid_name == "Scintillator-0") { step_edep_0 = delta_total_energy(step); }
    else if (step_solid_name == "Scintillator-1") { step_edep_1 = delta_total_energy(step); }
    else {                                          step_edep_0 = step_edep_1 = 0;          }

    total_edep_0 += step_edep_0;
    total_edep_1 += step_edep_1;
}

int main(int argc, char *argv[]) {

    G4double total_edep_0 = 0;
    G4double total_edep_1 = 0;
    std::ofstream data_file;

    // User action functions ------------------------
    // Generator
    auto two_gammas = [](auto event){ generate_back_to_back_511_keV_gammas(event, {}, 0); };

    // Run actions
    auto  open_file = [&data_file] (G4Run const*) { data_file.open("G4_data_test.csv"); };
    auto close_file = [&data_file] (G4Run const*) { data_file.close(); };

    // Event actions
    auto reset_total_edep = [&total_edep_0, &total_edep_1] (G4Event const*) {
        total_edep_0 = 0;
        total_edep_1 = 0;
    };
    auto print_total_edep = [&data_file, &total_edep_0, &total_edep_1] (G4Event const*) {
        G4cout << "\nTotal deposited energy in scintillator 0: " << total_edep_0
               << "\nTotal deposited energy in scintillator 1: " << total_edep_1 << G4endl << G4endl;

        data_file << total_edep_0 << "," << total_edep_1 << std::endl;
    };

    // Stepping action
    // This can be used to move the closure "out of scope" (add_step_edep is outside of main)
    auto accumulate_energy = [&total_edep_0, &total_edep_1] (G4Step const* step) { add_step_edep(total_edep_0, total_edep_1, step); };

    // Stacking action
    auto kill_secondaries = [] (G4Track const* track) {
        G4int parent_ID = track -> GetParentID();
        if (parent_ID > 0) { return G4ClassificationOfNewTrack::fKill;   }
        else               { return G4ClassificationOfNewTrack::fUrgent; }
    };

    // Setting mandatory G4 objects --------------------------
    auto run_manager = n4::run_manager::create()
        .physics(physics_list)
        .geometry(make_geometry)
        .actions([&] {return (new n4::actions{two_gammas})
            -> set((new n4::run_action())
                   -> begin(open_file)
                   -> end(close_file))
            -> set((new n4::event_action())
                   -> begin(reset_total_edep)
                   -> end(print_total_edep))
            -> set((new n4::stepping_action{accumulate_energy}))
            -> set((new n4::stacking_action())
                   -> classify(kill_secondaries));});

    n4::ui(argc, argv);

}
