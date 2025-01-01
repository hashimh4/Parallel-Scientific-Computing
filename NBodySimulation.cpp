#include "NBodySimulation.h"
#include <algorithm>
#include <vector>

NBodySimulation::NBodySimulation () :
  t(0), tFinal(0), tPlot(0), tPlotDelta(0), NumberOfBodies(0),
  x(nullptr), v(nullptr), mass(nullptr),
  timeStepSize(0), maxV(0), minDx(0), videoFile(nullptr),
  snapshotCounter(0), timeStepCounter(0) {};

NBodySimulation::~NBodySimulation () {
  if (x != nullptr) {
    for (int i=0; i<NumberOfBodies; i++)
      delete [] x[i];
    delete [] x;
  }
  if (v != nullptr) {
    for (int i=0; i<NumberOfBodies; i++)
      delete [] v[i];
    delete [] v;
  }
  if (mass != nullptr) {
    delete [] mass;
  }
}

void NBodySimulation::checkInput(int argc, char** argv) {
    if (argc==1) {
    std::cerr << "usage: " << std::string(argv[0])
              << " plot-time final-time dt objects" << std::endl
              << " Details:" << std::endl
              << " ----------------------------------" << std::endl
              << "  plot-time:       interval after how many time units to plot."
                 " Use 0 to switch off plotting" << std::endl
              << "  final-time:      simulated time (greater 0)" << std::endl
              << "  dt:              time step size (greater 0)" << std::endl
              << "  objects:         any number of bodies, specified by position, velocity, mass" << std::endl
              << std::endl
              << "Examples of arguments:" << std::endl
              << "+ One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
              << "    0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0" << std::endl
              << "+ One body spiralling around the other" << std::endl
              << "    0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0     0.0 1.0 0.0  1.0 0.0 0.0  1.0" << std::endl
              << "+ Three-body setup from first lecture" << std::endl
              << "    0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0" << std::endl
              << "+ Five-body setup" << std::endl
              << "    0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0     2.0 1.0 0.0  0.0 0.0 0.0  1.0     2.0 0.0 1.0  0.0 0.0 0.0  1.0" << std::endl
              << std::endl;

    throw -1;
  }
  else if ( (argc-4)%7!=0 ) {
    std::cerr << "error in arguments: each body is given by seven entries"
                 " (position, velocity, mass)" << std::endl;
    std::cerr << "got " << argc << " arguments"
                 " (three of them are reserved)" << std::endl;
    std::cerr << "run without arguments for usage instruction" << std::endl;
    throw -2;
  }
}

void NBodySimulation::setUp (int argc, char** argv) {

  checkInput(argc, argv);

  NumberOfBodies = (argc-4) / 7;

  x    = new double*[NumberOfBodies];
  v    = new double*[NumberOfBodies];
  mass = new double [NumberOfBodies];
  
  int readArgument = 1;

  tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
  tFinal       = std::stof(argv[readArgument]); readArgument++;
  timeStepSize = std::stof(argv[readArgument]); readArgument++;

  for (int i=0; i<NumberOfBodies; i++) {
    x[i] = new double[3];
    v[i] = new double[3];

    x[i][0] = std::stof(argv[readArgument]); readArgument++;
    x[i][1] = std::stof(argv[readArgument]); readArgument++;
    x[i][2] = std::stof(argv[readArgument]); readArgument++;

    v[i][0] = std::stof(argv[readArgument]); readArgument++;
    v[i][1] = std::stof(argv[readArgument]); readArgument++;
    v[i][2] = std::stof(argv[readArgument]); readArgument++;

    mass[i] = std::stof(argv[readArgument]); readArgument++;

    if (mass[i]<=0.0 ) {
      std::cerr << "invalid mass for body " << i << std::endl;
      exit(-2);
    }
  }

  std::cout << "created setup with " << NumberOfBodies << " bodies"
            << std::endl;

  if (tPlotDelta<=0.0) {
    std::cout << "plotting switched off" << std::endl;
    tPlot = tFinal + 1.0;
  }
  else {
    std::cout << "plot initial setup plus every " << tPlotDelta
              << " time units" << std::endl;
    tPlot = 0.0;
  }
}


double NBodySimulation::force_calculation (int i, int j, int direction){
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


bool NBodySimulation::hasReachedEnd () {
  return t > tFinal;
}

void NBodySimulation::takeSnapshot () {
  if (t >= tPlot) {
    printParaviewSnapshot();
    printSnapshotSummary();
    tPlot += tPlotDelta;
  }
}


void NBodySimulation::openParaviewVideoFile () {
  videoFile.open("paraview-output/result.pvd");
  videoFile << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\""
    " version=\"0.1\""
    " byte_order=\"LittleEndian\""
    " compressor=\"vtkZLibDataCompressor\">" << std::endl
            << "<Collection>";
}

void NBodySimulation::closeParaviewVideoFile () {
  videoFile << "</Collection>"
            << "</VTKFile>" << std::endl;
  videoFile.close();
}

void NBodySimulation::printParaviewSnapshot () {
  static int counter = -1;
  counter++;
  std::stringstream filename, filename_nofolder;
  filename << "paraview-output/result-" << counter <<  ".vtp";
  filename_nofolder << "result-" << counter <<  ".vtp";
  std::ofstream out( filename.str().c_str() );
  out << "<VTKFile type=\"PolyData\" >" << std::endl
      << "<PolyData>" << std::endl
      << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
      << "  <Points>" << std::endl
      << "   <DataArray type=\"Float64\""
    " NumberOfComponents=\"3\""
    " format=\"ascii\">";

  for (int i=0; i<NumberOfBodies; i++) {
    out << x[i][0]
        << " "
        << x[i][1]
        << " "
        << x[i][2]
        << " ";
  }

  out << "   </DataArray>" << std::endl
      << "  </Points>" << std::endl
      << " </Piece>" << std::endl
      << "</PolyData>" << std::endl
      << "</VTKFile>"  << std::endl;

  out.close();

  videoFile << "<DataSet timestep=\"" << counter
            << "\" group=\"\" part=\"0\" file=\"" << filename_nofolder.str()
            << "\"/>" << std::endl;
}

void NBodySimulation::printSnapshotSummary () {
  std::cout << "plot next snapshot"
            << ",\t time step=" << timeStepCounter
            << ",\t t="         << t
            << ",\t dt="        << timeStepSize
            << ",\t v_max="     << maxV
            << ",\t dx_min="    << minDx
            << std::endl;
}

void NBodySimulation::printSummary () {
  std::cout << "Number of remaining objects: " << NumberOfBodies << std::endl;
  std::cout << "Position of first remaining object: "
            << x[0][0] << ", " << x[0][1] << ", " << x[0][2] << std::endl;
}
