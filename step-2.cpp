#include <iomanip>

#include "NBodySimulation.h"
#include <algorithm>
#include <vector>

class NBodySimulationMolecularForces : public NBodySimulation {

  double r_c;
  // double* neighbour_delta_coordinate = new double[2];

  // Override the force calculation method
  double NBodySimulation::force_calculation (int i, int j, int direction, double r_c) override {
    
    // Euclidean distance
    const double distance = sqrt(
                                (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
                                (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
                                (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
                                );
    minDx = std::min( minDx,distance );

    // Calculate the molecular force instead
    if (distance <= r_c) {
      return 10 * (pow((0.1 / distance), 13) - pow((0.1 / distance), 9)) * distance;
    } else {
      return 0;
    }
  }

  void NBodySimulation::linked_cell (double lx, double ly, double r_c, int cell_index, double* neighbour_delta_coordinate, int a, double domain) {
    // Take in i and see what particles are close-by
    // Create a hash table

    double* list_particles = new double[NumberOfBodies];

    for (int j=0; j<NumberOfBodies; j++){
      double particle = list_particles[j]
      for (int i=0; i<NumberOfBodies; i++) {
          if (i != j) {
            //
          }

    }

  };
  

void NBodySimulation::updateBody () {

  timeStepCounter++;

  // force0 = force along x direction
  // force1 = force along y direction
  // force2 = force along z direction
  double* force0 = new double[NumberOfBodies];
  double* force1 = new double[NumberOfBodies];
  double* force2 = new double[NumberOfBodies];
  // Instead of resizing an array, add the merged particles to a vector list
  std::vector<int> otherIndexList;

  for (int j=0; j<NumberOfBodies; j++) {
    // Calculate only for the existing particles
    if (!(std::count(otherIndexList.begin(), otherIndexList.end(), j))) {
      int newIndex = 0;
      int otherIndex = 0;
      maxV   = 0.0;
      minDx  = std::numeric_limits<double>::max();

      force0[j] = 0.0;
      force1[j] = 0.0;
      force2[j] = 0.0;

      // Calculate the forces only on the nearest neigbours
      // Based on the classical linked-cell algorithm
      for (int i=0; i<NumberOfBodies; i++) {
        if (i != j) {
          // x,y,z forces acting on particle 0.
          force0[j] += force_calculation(i,j,0);
          force1[j] += force_calculation(i,j,1);
          force2[j] += force_calculation(i,j,2);
        }

        if (abs(x[i] - x[j]) <= 10*(-2) * mass[i] + mass[j]) {
          // Find the lower index
          if (i < j) {
            newIndex = i;
            otherIndexList.push_back(j);
          } else {
            newIndex = j;
            otherIndex = i;
            otherIndexList.push_back(i);
          }

          // Redefine the mass
          mass[newIndex] = mass[i] + mass[j];
          // Instead of resizing the array, set all the variables for the second particle to zero
          mass[otherIndex] = 0;

          // Redefine the position
          x[newIndex][0] = (mass[i] * x[i][0] + mass[j] * x[j][0]) / (mass[i] + mass[j]);
          x[newIndex][1] = (mass[i] * x[i][1] + mass[j] * x[j][1]) / (mass[i] + mass[j]);
          x[newIndex][2] = (mass[i] * x[i][2] + mass[j] * x[j][2]) / (mass[i] + mass[j]);
          x[otherIndex] = 0;

          // Redefine the velocity
          v[newIndex][0] = (mass[i] * v[i][0] + mass[j] * v[j][0]) / (mass[i] + mass[j]);
          v[newIndex][1] = (mass[i] * v[i][1] + mass[j] * v[j][1]) / (mass[i] + mass[j]);
          v[newIndex][2] = (mass[i] * v[i][2] + mass[j] * v[j][2]) / (mass[i] + mass[j]);
          v[otherIndex] = 0;
        }
      }
    }

    x[j][0] = x[j][0] + timeStepSize * v[j][0];
    x[j][1] = x[j][1] + timeStepSize * v[j][1];
    x[j][2] = x[j][2] + timeStepSize * v[j][2];

    v[j][0] = v[j][0] + timeStepSize * force0[j] / mass[j];
    v[j][1] = v[j][1] + timeStepSize * force1[j] / mass[j];
    v[j][2] = v[j][2] + timeStepSize * force2[j] / mass[j];

    maxV = std::sqrt( v[j][0]*v[j][0] + v[j][1]*v[j][1] + v[j][2]*v[j][2] );
  }

  // Allow a variable, stable timestep
  t += timeStepSize;

  // Variables defined in the makefile
  // Calculate the error
  f = (1 + lambda * timeStepSize) * f;
  e = f - f0 * exp(lambda * t);

  // if (abs(1 + lambda * timeStepSize) > 1) {
  // }

  delete[] force0;
  delete[] force1;
  delete[] force2;
}

int main (int argc, char** argv) {

  std::cout << std::setprecision(15);

  // Code that initialises and runs the simulation.
  NBodySimulationMolecularForces nbs;
  nbs.setUp(argc,argv);
  nbs.openParaviewVideoFile();
  nbs.takeSnapshot();

  while (!nbs.hasReachedEnd()) {
    nbs.updateBody();
    nbs.takeSnapshot();
  }

  nbs.printSummary();
  nbs.closeParaviewVideoFile();

  return 0;
}
