//
//  main.cpp
//  gromacs-water-analysis
//
//  Created by Yiming Tang on 08/01/2017.
//  Copyright © 2017 Yiming Tang. All rights reserved.
//

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <fstream>
#include <cmath>
#include <vector>
#include "xdrfile_xtc.h"
#include "mol_water.hpp"
#include "mol_water.cpp"
#define partBreak  "***************************************************" << endl;
#define MOMENT_SINGLE 1.0

#define Direction_Step 1
#define Z_Step 0.05
#define STEP_STEP 200


using namespace std;

int main(int argc, const char * argv[]) {

    // Initializing and Testing Arguments from Command Line
    
    cout << "Welcome to GROMACS-WATER-ANALYSIS by Yiming Tang @ Fudan" << endl;
    cout << "WARNING: MUST CONVERT XTC FILE INTO WATER MOLECULE ONLY!" << endl;
    
    if (argc != 5)
    {
        
        cout << "usage: tym_gmx_water waterModel inputFile output_for_Coordinate output_for_DipoleMoment" << endl;
        cout << "Thanks~" << endl;
        exit(0);
    }
    
    if ( (strcmp(argv[1],"spce")))
    {
        cerr << "Cannot process water model " << argv[1];
        exit(0);
    }
    
    // Getting File Name From Command Line
    
    char inputFileName[100], outputCoorName[100], outputDirecName[100];
    strcpy(inputFileName, argv[2]);
    strcpy(outputCoorName, argv[3]);
    strcpy(outputDirecName, argv[4]);
    
    ofstream outFile_Coor_Up;
    ofstream outFile_Coor_Down;
    ofstream outFile_Moment_Up;
    ofstream outFile_Moment_Down;

    char outputCoorName_up[100], outputCoorName_down[100], outputMomentName_up[100], outputMomentName_down[100];


    strcpy(outputCoorName_up,outputCoorName);
    strcpy(outputCoorName_down,outputCoorName);
    strcpy(outputMomentName_up,outputDirecName);
    strcpy(outputMomentName_down,outputDirecName);

    strcat(outputCoorName_up,"_up.txt");
    strcat(outputCoorName_down, "_down.txt");
    strcat(outputMomentName_up, "_up.txt");
    strcat(outputMomentName_down, "_down.txt");
    cout << "Will Read File From       " << inputFileName << endl;
    cout << "Will Write Coordinate To  " << outputCoorName_up << endl;
    cout << "Will Write Direction To   " << outputMomentName_up << endl;
    cout << partBreak;
    
    // opening File
    
    int atom_count;
    int step_current;
    float time_current;
    float p;
    matrix box;
    rvec *coor_current;
    XDRFILE *xtc;
    float Z_Max;
    int mystep=1;
    
    mol_water tempWater;
    int *molCount_up;
    int *molCount_down;
    double *directionAverage_up;
    double *directionAverage_down;
    

    xtc = xdrfile_open(inputFileName,"r");
    if ( xtc == NULL )
    {
        cerr << "Error Opening Input File!" << endl;
        exit(0);
    }
    cout << "Succeed in Opening Input File!" << endl;


    outFile_Coor_Up.open(outputCoorName_up);
    outFile_Coor_Down.open(outputCoorName_down);
    outFile_Moment_Up.open(outputMomentName_up);
    outFile_Moment_Down.open(outputMomentName_down);

    if(!(outFile_Coor_Up.is_open()) || !(outFile_Coor_Down.is_open()))
    {
        cerr << "Error Opening Output File for Coordination!" << endl;
        exit(0);
    }
    cout << "Succeed in Opening Output File for Coordination!" << endl;
    
    if(!(outFile_Moment_Up.is_open()) || !(outFile_Moment_Down.is_open()))
    {
        cerr << "Error Opening Output File for Direction!" << endl;
        exit(0);
    }
    cout << "Succeed in Opening Output File for Direction!" << endl;
    
    int read_xtc_return = read_xtc_natoms(inputFileName, &atom_count);
    
    if (atom_count % 3 != 0)
    {
        cerr << "Processing total number of atoms not a multiple of 3!" << endl;
        cerr << "Check if your xtc file contains molecules other than water!" << endl;
        exit(0);
    }
    
        cout << partBreak;
    
    coor_current = (rvec *)calloc(atom_count, sizeof(coor_current[0]));
    
    // Reading first frame
    
    read_xtc_return = read_xtc(xtc,atom_count,&step_current,&time_current,box,coor_current,&p);
    
    if(read_xtc_return != 0)
    {
        cerr << "Cannot Read the first frame!" << endl;
        exit(0);
    }
    cout << "Succeed in Reading the first frame. Now processing!" << endl;
    cout << "Will process a total number of water molecules as " << atom_count/4 << endl;
    Z_Max = box[2][2];
    cout << "Processing Box With Maximum Z of " << Z_Max << endl;

    molCount_up = new int[(int)ceil(Z_Max / Z_Step)];
    molCount_down = new int[(int)ceil(Z_Max / Z_Step)];
    directionAverage_up = new double[(int)ceil(Z_Max / Z_Step)];
    directionAverage_down = new double[(int)ceil(Z_Max / Z_Step)];
    
    for(int i = 0; i < (int)ceil(Z_Max/Z_Step); i++)
    {
        molCount_up[i] = 0;
        molCount_down[i] = 0;
        directionAverage_up[i] = 0;
        directionAverage_down[i] = 0;
    }

    for(int atom_number = 1; atom_number <= atom_count ; atom_number++)
    {
        switch (atom_number % 3)
        {
            case 1:
                tempWater.update_O(coor_current[atom_number][0], coor_current[atom_number][1], coor_current[atom_number][2]);
                break;
            case 2:
                tempWater.update_H1(coor_current[atom_number][0], coor_current[atom_number][1], coor_current[atom_number][2]);
                break;
            case 3:
                tempWater.update_H2(coor_current[atom_number][0], coor_current[atom_number][1], coor_current[atom_number][2]);
                if(tempWater.getCoordinate().Z>0)
                {
                    molCount_up[(int)floor(tempWater.getCoordinate().Z / Z_Step)]++;
                    directionAverage_up[(int)floor(tempWater.getCoordinate().Z / Z_Step)] += MOMENT_SINGLE * cos(tempWater.getDirectionAngle() / 180 * M_PI);
                }
                else
                {
                    molCount_down[(int)floor(tempWater.getCoordinate().Z / Z_Step)]++;
                    directionAverage_down[(int)floor(tempWater.getCoordinate().Z / Z_Step)] += MOMENT_SINGLE * cos(tempWater.getDirectionAngle() / 180 * M_PI);
                }
                break;
        }
    }
    
    cout << "Succeed in Processing the first frame. Going to futher frames!" << endl;
    
    outFile_Coor_Up << "STEP" << '\t';
    outFile_Coor_Down << "STEP" << '\t';
    outFile_Moment_Up << "STEP" << '\t';
    outFile_Moment_Down << "STEP" << '\t';

    for(int i = 0; i < (int)ceil(Z_Max/Z_Step); i++)
    {
        outFile_Coor_Up << i*Z_Step << "~" << (i+1)*Z_Step << '\t' ;
        outFile_Coor_Down << i*Z_Step << "~" << (i+1)*Z_Step << '\t' ;
        outFile_Moment_Up << i*Z_Step << "~" << (i+1)*Z_Step << '\t' ;
        outFile_Moment_Down << i*Z_Step << "~" << (i+1)*Z_Step << '\t' ;
    }
    outFile_Coor_Up << endl;
    outFile_Coor_Down << endl;
    outFile_Moment_Up << endl;
    outFile_Moment_Up << endl;
    
    while(1)
    {
        read_xtc_return = read_xtc(xtc,atom_count,&step_current,&time_current,box,coor_current,&p);
        
        if(read_xtc_return != 0)
        {
            break;
        }
        
        mystep = mystep + 1;
        
        for(int atom_number = 1; atom_number <=  atom_count ; atom_number++)
        {
            switch (atom_number % 3)
            {
                case 1:
                    tempWater.update_O(coor_current[atom_number][0], coor_current[atom_number][1], coor_current[atom_number][2]);
                    break;
                case 2:
                    tempWater.update_H1(coor_current[atom_number][0], coor_current[atom_number][1], coor_current[atom_number][2]);
                    break;
                case 3:
                    tempWater.update_H2(coor_current[atom_number][0], coor_current[atom_number][1], coor_current[atom_number][2]);
                    if(tempWater.getCoordinate().Z>0)
                    {
                        molCount_up[(int)floor(tempWater.getCoordinate().Z / Z_Step)]++;
                        directionAverage_up[(int)floor(tempWater.getCoordinate().Z / Z_Step)] += MOMENT_SINGLE * cos(tempWater.getDirectionAngle() / 180 * M_PI);
                    }
                    else
                    {
                        molCount_down[(int)floor(tempWater.getCoordinate().Z / Z_Step)]++;
                        directionAverage_down[(int)floor(tempWater.getCoordinate().Z / Z_Step)] += MOMENT_SINGLE * cos(tempWater.getDirectionAngle() / 180 * M_PI);
                    }
                    break;
            }
        }

        
        if (mystep % STEP_STEP == 0)
        {
            outFile_Coor_Up << mystep - STEP_STEP + 1 << "~" << mystep << '\t';
            outFile_Coor_Down << mystep - STEP_STEP + 1 << "~" << mystep << '\t';
            outFile_Moment_Up << mystep - STEP_STEP + 1 << "~" << mystep << '\t';
            outFile_Moment_Down << mystep - STEP_STEP + 1 << "~" << mystep << '\t';
            for(int i = 0; i < (int)ceil(Z_Max/Z_Step) ; i++)
            {
                outFile_Coor_Up << molCount_up[i]  / (float)STEP_STEP  << '\t' ;
                outFile_Coor_Down << molCount_down[i] / (float)STEP_STEP << '\t' ;
                outFile_Moment_Up << directionAverage_up[i] / molCount_up[i] << '\t';
                outFile_Moment_Down << directionAverage_down[i] / molCount_down[i] << '\t';
            }
            outFile_Coor_Up << endl;
            outFile_Coor_Down << endl;
            outFile_Moment_Up << endl;
            outFile_Moment_Down << endl;
            
            for(int i = 0; i < (int)ceil(Z_Max/Z_Step); i++)
            {
                molCount_up[i] = 0;
                molCount_down[i] = 0;
                directionAverage_up[i] = 0;
                directionAverage_down[i] = 0;
            }
        }
    }
    
    cout << "Succeed in Processing all frames. Now write file!" << endl;
    
    xdrfile_close(xtc);
    outFile_Coor_Up.close();
    outFile_Coor_Down.close();
    outFile_Moment_Up.close();
    outFile_Moment_Down.close();
    
    std::cout << "Hello, World!\n";
    return 0;
}
