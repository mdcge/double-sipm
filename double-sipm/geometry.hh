#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <G4PVPlacement.hh>

G4PVPlacement* make_geometry(std::vector<int>& photon_count, std::vector<std::vector<G4double>>&);

#endif // GEOMETRY_H_
