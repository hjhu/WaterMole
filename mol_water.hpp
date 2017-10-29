//
//  mol_water.hpp
//  gromacs-water-analysis
//
//  Created by Yiming Tang on 08/01/2017.
//  Copyright © 2017 Yiming Tang. All rights reserved.
//

#ifndef mol_water_hpp
#define mol_water_hpp

#include <stdio.h>
#include "xdrfile_xtc.h"

struct coor
{
    double X;
    double Y;
    double Z;
};

class mol_water
{
private:
    coor atom_O;
    coor atom_H1;
    coor atom_H2;
    
public:
    mol_water();
    mol_water(coor input_O, coor input_H1, coor input_H2);
    
    void update_O(coor input_O);
    void update_H1(coor input_H1);
    void update_H2(coor input_H2);
    
    void update_O(double x, double y, double z);
    void update_H1(double x, double y, double z);
    void update_H2(double x, double y, double z);
    
    coor getCoordinate();
    coor getDirection();
    float getDirectionAngle();
};



#endif /* mol_water_hpp */
