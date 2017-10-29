//
//  mol_water.cpp
//  gromacs-water-analysis
//
//  Created by Yiming Tang on 08/01/2017.
//  Copyright © 2017 Yiming Tang. All rights reserved.
//

#include "mol_water.hpp"
#include <cmath>
#include <iostream>

using namespace std;

mol_water::mol_water()
{
    atom_O.X = atom_O.Y = atom_O.Z = 0;
    atom_H1.X = atom_H1.Y = atom_H1.Z = 0;
    atom_H2.X = atom_H2.Z = atom_H2.Z = 0;
}

mol_water::mol_water(coor input_O, coor input_H1, coor input_H2)
{
    atom_O = input_O;
    atom_H1 = input_H1;
    atom_H2 = input_H2;
}

void mol_water::update_O(coor input_O)
{
    atom_O = input_O;
}

void mol_water::update_H1(coor input_H1)
{
    atom_H1 = input_H1;
}

void mol_water::update_H2(coor input_H2)
{
    atom_H2 = input_H2;
}


void mol_water::update_O(double x, double y, double z)
{
    atom_O.X = x;
    atom_O.Y = y;
    atom_O.Z = z;
}

void mol_water::update_H1(double x, double y, double z)
{
    atom_H1.X = x;
    atom_H1.Y = y;
    atom_H1.Z = z;
}

void mol_water::update_H2(double x, double y, double z)
{
    atom_H2.X = x;
    atom_H2.Y = y;
    atom_H2.Z = z;
}



coor mol_water::getCoordinate()
{
    if (atom_O.X == 0 && atom_O.Y == 0 && atom_O.Z == 0)
    {
        std::cerr << "Atom not initialized!" << std::endl;
        exit(0);
    }
    return atom_O;
}

coor mol_water::getDirection()
{
    if (atom_O.X == 0 && atom_O.Y == 0 && atom_O.Z == 0)
    {
        std::cerr << "Atom not initialized!" << std::endl;
        exit(0);
    }
    if (atom_H1.X == 0 && atom_H1.Y == 0 && atom_H1.Z == 0)
    {
        std::cerr << "Atom not initialized!" << std::endl;
        exit(0);
    }
    if (atom_H2.X == 0 && atom_H2.Y == 0 && atom_H2.Z == 0)
    {
        std::cerr << "Atom not initialized!" << std::endl;
        exit(0);
    }
    struct coor outputDirection = { \
        (atom_H1.X + atom_H2.X) / 2 - atom_O.X , \
        (atom_H1.Y + atom_H2.Y) / 2 - atom_O.Y , \
        (atom_H1.Z + atom_H2.Z) / 2 - atom_O.Z };
    return outputDirection;
}

float mol_water::getDirectionAngle()
{
    coor tempCoor = getDirection();
    return acos(tempCoor.Z / sqrt(tempCoor.X * tempCoor.X + tempCoor.Y * tempCoor.Y + tempCoor.Z * tempCoor.Z)) * 180.0 / M_PI;
}










