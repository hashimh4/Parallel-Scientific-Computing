#include "NBodySimulation.h"
#include <algorithm>
#include <vector>

class NBodySimulationVectorised : public NBodySimulation {
    // Perform automatic vectorisation 
    // -O3 -fopt-info
    // Unroll the right number of iterations of the loop

    //Unroll safelen()

    double NBodySimulationVectorised::force_calculation (int i, int j, int direction){
    // Euclidean distance
    const double distance = sqrt(
                                (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
                                (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
                                (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
                                );
    const double distance3 = distance * distance * distance;
    minDx = std::min( minDx,distance );

    return (x[i][direction]-x[j][direction]) * mass[i]*mass[j] / distance3;
    }

    // Assume there is no aliassing 
    #pragma omp simd
    void NBodySimulationVectorised::updateBody () {

    timeStepCounter++;

    // force0 = force along x direction
    // force1 = force along y direction
    // force2 = force along z direction
    double* force0 = new double[NumberOfBodies];
    double* force1 = new double[NumberOfBodies];
    double* force2 = new double[NumberOfBodies];
    // Instead of resizing an array, add the merged particles to a vector list
    std::vector<int> otherIndexList;

    // Vectorise both loops
    #pragma omp simd collapse(2)
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

        // For the loop containing the calculation function
        #pragma omp declare simd
        for (int i=0; i<NumberOfBodies; i++) {
            if (i != j) {
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


};
