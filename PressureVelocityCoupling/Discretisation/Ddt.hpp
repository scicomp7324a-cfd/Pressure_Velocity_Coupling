# pragma once

# include <iostream>
# include <map>
# include "FaceAddressedMesh2D.hpp"
# include "FVMatrix.hpp"

using namespace std ;

/* *******************************************************************
    Ddt.hpp - This class contains the temporal assembly 
    function. 
******************************************************************* */

class Ddt{

    private:

        FVMatrix& fvMatrix_ ;
        const FaceAddressedMesh2D& mesh_ ;
        const map<int, double>& phiOld_ ;
        double deltaT_ ;

    public: 

        // Constructs the temporal term assembler using the mesh, old field values and time steps.
        Ddt(FVMatrix& fvMatrix, FaceAddressedMesh2D& mesh, map<int, 
                                        double>& phiOld, double deltaT):fvMatrix_(fvMatrix), 
                                                          mesh_(mesh),
                                                          phiOld_(phiOld),
                                                          deltaT_(deltaT){}

        // This functions assembles the diagonal and source contributions of the temporal term.
        void assemble(){
            for(int i = 0; i < mesh_.numCells(); ++i){
                double volume = mesh_.cellArea(i) ;
                fvMatrix_.addDiagonal(i, volume/deltaT_) ;
                fvMatrix_.addToSource(i, volume*phiOld_.at(i)/deltaT_) ;
            }
        }

} ;