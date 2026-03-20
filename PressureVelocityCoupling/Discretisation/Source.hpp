#pragma once

#include <iostream>
#include <map>
#include <array>
#include "FaceAddressedMesh2D.hpp"
#include "FVMatrix.hpp"

using namespace std;

/* *******************************************************************
    Source.hpp - This class contains the source assembly 
    function. 
******************************************************************* */

class Source
{

private:
    FVMatrix &fvMatrix_;
    const FaceAddressedMesh2D &mesh_;
    const map<int, array<double, 2>> &source_;

public:
    
    // Constructs the source term assembler using the target matrix, mesh, and per-cell source coefficients.
    Source(FVMatrix &fvMatrix, FaceAddressedMesh2D &mesh,
           map<int, array<double, 2>> &source) : fvMatrix_(fvMatrix),
                                                 mesh_(mesh),
                                                 source_(source) {}

    // This function assembles. Source is defined as Su - Sp*PhiP
    void assemble()
    {
        for (int i = 0; i < mesh_.numCells(); ++i)
        {
            double volume = mesh_.cellArea(i) ;
            fvMatrix_.addDiagonal(i, source_.at(i)[1]*volume) ;
            fvMatrix_.addToSource(i, source_.at(i)[0]*volume) ;
        }
    }
};