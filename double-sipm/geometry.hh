#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <G4PVPlacement.hh>

extern G4int photon_count_0;
extern G4int photon_count_1;
extern std::vector<G4double> times_of_arrival_0;
extern std::vector<G4double> times_of_arrival_1;

G4PVPlacement* make_geometry(std::vector<int>& photon_count);

#endif // GEOMETRY_H_
