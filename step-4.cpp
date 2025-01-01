#include <iomanip>

#include "NBodySimulationVectorised.cpp"

class NBodySimulationParallelised : public NBodySimulationVectorised {
  // Copy paste and then ...
  
  #pragma omp parallel

  #pragma omp for



};

int main (int argc, char** argv) {

  std::cout << std::setprecision(15);

  // Code that initialises and runs the simulation.
  NBodySimulationVectorised nbs;
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
