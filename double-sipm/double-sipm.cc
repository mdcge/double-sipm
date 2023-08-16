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
    return step -> GetTotalEnergyDeposit();
}

// Add energies deposited in this step to running totals of deposited energies in whole event
void add_step_edep(std::vector<G4double>& total_edep, G4Step const* step) {
    auto step_solid_name = step -> GetPreStepPoint() -> GetTouchable() -> GetVolume() -> GetName();
    if      (step_solid_name == "Scintillator-0") { total_edep[0] += delta_total_energy(step); }
    else if (step_solid_name == "Scintillator-1") { total_edep[1] += delta_total_energy(step); }
}

int main(int argc, char *argv[]) {

    // Each event produces a pair of back-to-back 511 keV gammas
    auto two_gammas = [](auto event){ generate_back_to_back_511_keV_gammas(event, {}, 0); };

    // Open output file at start of run, close it at the end of the run

    std::ofstream data_file_0; std::ofstream data_file_1;
    std::vector<std::ofstream*> data_file{&data_file_0, &data_file_1};
    auto open_file = [&data_file] (G4Run const*) {
        data_file[0] -> open("G4_photon_times_0.csv");
        data_file[1] -> open("G4_photon_times_1.csv");
    };
    auto close_file = [&data_file] (G4Run const*) { data_file[0] -> close(); data_file[1] -> close(); };

    // Accumulators for energy and photons observed in each scintillator during a single event
    std::vector<G4double> total_edep{0, 0};
    std::vector<std::vector<G4double>> times_of_arrival{{}, {}};
    // At the start of each event: reset the accumulators to zero
    auto reset_photon_count = [&total_edep, &times_of_arrival] (G4Event const*) {
        total_edep  [0] = total_edep  [1] = 0;
        times_of_arrival[0].clear();
        times_of_arrival[1].clear();
    };

    G4int double_hits = 0;
    auto write_photon_count = [&data_file, &double_hits, &total_edep, &times_of_arrival] (G4Event const* event) {
        //G4cout << "Event number: " << event -> GetEventID() << G4endl;
        // G4cout << G4endl << "Total deposited energy in scintillator 0: " << total_edep_0 << G4endl;
        // G4cout << "Total deposited energy in scintillator 1: " << total_edep_1 << G4endl;
        // G4cout << "Photon count 0: " << photon_count[0] << G4endl;
        // G4cout << "Photon count 1: " << photon_count[1] << G4endl << G4endl;

        G4cout << "Number of double events: " << double_hits << "/" << event -> GetEventID() << " events" << G4endl;
        // G4cout << "Number of times: " << times_of_arrival[0].size() << " and " << times_of_arrival[1].size() << G4endl;
        std::vector<size_t> photon_count{times_of_arrival[0].size(), times_of_arrival[0].size()};
        std::vector<size_t> photon_count{times_of_arrival[0].size(), times_of_arrival[1].size()};
        if (photon_count[0] > 0 && photon_count[1] > 0) { double_hits += 1; }
        for (int side=0; side < 2; ++side) {
            for (size_t i=0; i<times_of_arrival[side].size(); i++) {
                if (photon_count[side] != 0) { *data_file[side] << times_of_arrival[side][i] << ","; }
                else                         { *data_file[side] << 0; }
            }
            *data_file[side] << std::endl;
        }
    };

    // At every step: increment running total of deposited energy during the event
    auto accumulate_energy = [&total_edep] (G4Step const* step) { add_step_edep(total_edep, step); };

    // If needed, use this as stacking_action -> classify, to disable tracking
    // of scintillation products
    auto kill_secondaries = [] (G4Track const* track) {
        G4int parent_ID = track -> GetParentID();
        if (parent_ID > 0) { return G4ClassificationOfNewTrack::fKill;   }
        else               { return G4ClassificationOfNewTrack::fUrgent; }
    };

    // Setting mandatory G4 objects --------------------------
    auto run_manager = n4::run_manager::create()
        .physics(physics_list)
        .geometry([&]{ return make_geometry(times_of_arrival); })
        .actions( [&]{ return (new n4::actions{two_gammas})
            -> set((new n4::run_action())
                   -> begin(open_file)
                   -> end (close_file))
            -> set((new n4::event_action())
                   -> begin(reset_photon_count)
                   -> end(write_photon_count))
            -> set((new n4::stepping_action{accumulate_energy}));});
            // -> set((new n4::stacking_action())
            //        -> classify(kill_secondaries));});

    // Run the simulation

    // + No CLI arguments: open GUI (using the settings in macs/vis.mac)

    // + 1 CLI argument (a macro file such as `macs/run.mac`): run in
    //     batch mode, with the specified macro

    // Batch mode will run the simulation much more quickly, because it will not
    // spend resources on drawing trajectories.
    n4::ui(argc, argv);
}
