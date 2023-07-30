#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <G4PVPlacement.hh>

G4PVPlacement* make_geometry(std::vector<std::vector<G4double>>& times_of_arrival);

#endif // GEOMETRY_H_
