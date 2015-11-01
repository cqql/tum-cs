#include "VTKStencil.h"

#define GHOST_OFFSET 2

#include <sstream>
#include <fstream>
#include <algorithm>

VTKStencil::VTKStencil(const Parameters& parameters)
  : FieldStencil(parameters) {
}

void VTKStencil::apply(FlowField& flowField, int i, int j) {

}

void VTKStencil::apply(FlowField& flowField, int i, int j, int k) {

}

void VTKStencil::write(FlowField& flowField, int timeStep) {
  GeometricParameters geom = this->_parameters.geometry;
  Meshsize* mesh = this->_parameters.meshsize;
  IntScalarField& flags = flowField.getFlags();

  // +3 is for the ghost layer
  int cellsX = geom.sizeX;
  int cellsY = geom.sizeY;
  int cellsZ = geom.sizeZ;
  int cells = cellsX * cellsY * cellsZ;
  int pointsX = cellsX + 1;
  int pointsY = cellsY + 1;
  int pointsZ = 1;
  float pressures[cellsZ][cellsY][cellsX];
  float velocitiesX[cellsZ][cellsY][cellsX];
  float velocitiesY[cellsZ][cellsY][cellsX];
  float velocitiesZ[cellsZ][cellsY][cellsX];

  if (geom.dim == 3) {
    pointsZ = cellsZ + 1;
  }

  int points = pointsX * pointsY * pointsZ;

  for (int k = 0; k < cellsZ; k++) {
    for (int j = 0; j < cellsY; j++) {
      for (int i = 0; i < cellsX; i++) {
        FLOAT p;
        FLOAT v[3];
        int x = i + GHOST_OFFSET;
        int y = j + GHOST_OFFSET;
        int z = k + GHOST_OFFSET;

        if ((flags.getValue(x, y, k) & OBSTACLE_SELF) == 0) {
          if (geom.dim == 2) {
            flowField.getPressureAndVelocity(p, v, x, y);
          } else {
            flowField.getPressureAndVelocity(p, v, x, y, z);
          }
        } else {
          p = 0.0;
        }

        pressures[k][j][i] = p;
      }
    }
  }

  for (int k = 0; k < cellsZ; k++) {
    for (int j = 0; j < cellsY; j++) {
      for (int i = 0; i < cellsX; i++) {
        FLOAT p;
        FLOAT v[3] = {0, 0, 0};
        int x = i + GHOST_OFFSET;
        int y = j + GHOST_OFFSET;
        int z = k + GHOST_OFFSET;

        if ((flags.getValue(x, y, k) & OBSTACLE_SELF) == 0) {
          if (geom.dim == 2) {
            flowField.getPressureAndVelocity(p, v, x, y);
          } else {
            flowField.getPressureAndVelocity(p, v, x, y, z);
          }
        }

        velocitiesX[k][j][i] = v[0];
        velocitiesY[k][j][i] = v[1];
        velocitiesZ[k][j][i] = v[2];
      }
    }
  }

  std::stringstream fstream;
  fstream << _parameters.vtk.prefix << "." << timeStep << ".vtk";
  std::string filename = fstream.str();
  std::ofstream file;
  file.open(filename.c_str(), std::ios::out);

  // Print floats with fixed precision
  file << std::fixed;

  file << "# vtk DataFile Version 2.0\n";
  file << "NS-EOF\n";
  file << "ASCII\n";
  file << "DATASET STRUCTURED_GRID\n";
  file << "DIMENSIONS " << pointsX << " " << pointsY << " " << pointsZ << "\n";
  file << "POINTS " << points << " float\n";

  for (int k = GHOST_OFFSET; k < pointsZ + GHOST_OFFSET; k++) {
    for (int j = GHOST_OFFSET; j < pointsY + GHOST_OFFSET; j++) {
      for (int i = GHOST_OFFSET; i < pointsX + GHOST_OFFSET; i++) {
        FLOAT posX = mesh->getPosX(i, j, k);
        FLOAT posY = mesh->getPosY(i, j, k);
        FLOAT posZ = mesh->getPosZ(i, j, k);

        file << posX << " " << posY << " " << posZ << "\n";
      }
    }
  }

  file << "CELL_DATA " << cells << "\n";

  file << "SCALARS pressure float 1\n";
  file << "LOOKUP_TABLE default\n";

  for (int i = 0; i < cells; i++) {
    file << ((float*)pressures)[i] << "\n";
  }

  file << "VECTORS velocity float\n";

  for (int i = 0; i < cells; i++) {
    file << ((float*)velocitiesX)[i] << " "
         << ((float*)velocitiesY)[i] << " "
         << ((float*)velocitiesZ)[i] << "\n";
  }

  file.close();
}
