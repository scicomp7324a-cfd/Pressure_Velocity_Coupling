# pragma once

# include <iostream>
# include <vector>
# include <array>
# include <math.h>
# include <map>
# include "FaceAddressedMesh2D.hpp"
# include "SparseLinearVector.hpp"
# include "FaceFluxField.hpp"
# include "FVMatrix.hpp"
# include "Field.hpp"
# include "Patch.hpp"
# include "ScalarProperty.hpp"
# include "Patch.hpp"

using namespace std ;

/* *******************************************************************
    PressureLaplacian.hpp -  This class assembles the pressure-correction 
    Laplacian equation from momentum diagonal coefficients, face fluxes, 
    and pressure boundary conditions.
******************************************************************* */

class PressureLaplacian{

    private:
    FVMatrix& fvMatrix_ ;
    FVMatrix& Ux_ ;
    FVMatrix& Uy_ ;
    const FaceAddressedMesh2D& mesh_ ;
    FaceFluxField& flux_ ;
    Field& fieldInitBC_ ;


    public:

    // Constructs the pressure Laplacian assembler using pressure, momentum, mesh, flux, and boundary field data.
    PressureLaplacian(FVMatrix& fvMatrix, FVMatrix& Ux, FVMatrix& Uy, const FaceAddressedMesh2D& mesh, FaceFluxField& flux, Field& fieldInitBC):fvMatrix_(fvMatrix),
                 Ux_(Ux),
                 Uy_(Uy),
                 mesh_(mesh),
                 flux_(flux),
                 fieldInitBC_(fieldInitBC){}
    
    // Builds the pressure equation by adding internal-face diffusion terms, boundary contributions, and a reference pressure constraint.
    void assemble(){

        for(int faceIndex: mesh_.internalFaces()){

            array<int,2> cellOwnerNeighbour = mesh_.faceOwnerNeighbour(faceIndex) ;
            array<double,2> faceAreaVector = mesh_.faceVector(faceIndex) ;
            double df = distanceCellToCell(faceIndex) ;
            double sf = sqrt(faceAreaVector[0]*faceAreaVector[0] + faceAreaVector[1]*faceAreaVector[1]) ;

            double invAxf = (1/Ux_.sys().A()(cellOwnerNeighbour[0],cellOwnerNeighbour[0]) + 1/Ux_.sys().A()(cellOwnerNeighbour[1], cellOwnerNeighbour[1]))*0.5 ;
            double invAyf = (1/Uy_.sys().A()(cellOwnerNeighbour[0],cellOwnerNeighbour[0]) + 1/Uy_.sys().A()(cellOwnerNeighbour[1], cellOwnerNeighbour[1]))*0.5 ;

            double sfx2 = faceAreaVector[0]*faceAreaVector[0] ;
            double sfy2 = faceAreaVector[1]*faceAreaVector[1] ;
            
            double Df = (sfx2*invAxf + sfy2*invAyf)/(df*sf) ;

            fvMatrix_.addDiagonal(cellOwnerNeighbour[0], Df) ;
            fvMatrix_.addCoefficient(cellOwnerNeighbour[0], cellOwnerNeighbour[1], -Df) ;
            fvMatrix_.addToSource(cellOwnerNeighbour[0], -flux_(faceIndex)) ;

            fvMatrix_.addDiagonal(cellOwnerNeighbour[1], Df) ;
            fvMatrix_.addCoefficient(cellOwnerNeighbour[1], cellOwnerNeighbour[0], -Df) ;
            fvMatrix_.addToSource(cellOwnerNeighbour[1], flux_(faceIndex)) ;
    
        }

        for(Patch p: mesh_.boundaryPatches()){
            map<int, double> boundaryValues = fieldInitBC_.boundary()[p.index()]->boundaryValues() ;
            for(const auto& [faceIndex, val]: boundaryValues){
                int cellOwner = mesh_.faceOwnerNeighbour(faceIndex)[0] ;
                array<double,2> faceAreaVector = mesh_.faceVector(faceIndex) ;
                double Sf = sqrt(faceAreaVector[0]*faceAreaVector[0] + faceAreaVector[1]*faceAreaVector[1]) ;
                array<double,2> normal = {faceAreaVector[0]/Sf, faceAreaVector[1]/Sf} ;
                double gamma = normal[0]*normal[0]/Ux_.sys().A()(cellOwner, cellOwner) + normal[1]*normal[1]/Uy_.sys().A()(cellOwner, cellOwner) ;
                fvMatrix_.addToSource(cellOwner, -gamma*Sf*val) ;
                fvMatrix_.addToSource(cellOwner, -flux_(faceIndex)) ;
            }

        }
    
        
        // Setting pressure reference by setting the pressure at cell 0 to be 0
        vector<int> neighbours = mesh_.cellNeighbours(0) ;
        fvMatrix_.setSource(0, 0.0) ;
        fvMatrix_.setDiagonal(0, 1.0) ;

        for(int neighbour: neighbours){
            fvMatrix_.setCoefficient(0, neighbour, 0.0) ;
        }
        

    }

    // Utility function to calculate distance between the centre of two cells
    double distanceCellToCell(const int& faceIndex){

        array<double,2> faceAreaVec = mesh_.faceVector(faceIndex) ;
        double sfMag = sqrt(faceAreaVec[0]*faceAreaVec[0] + faceAreaVec[1]*faceAreaVec[1]) ;
        faceAreaVec = {faceAreaVec[0]/sfMag, faceAreaVec[1]/sfMag};
        array<double,2> p1 = mesh_.cellCentre(mesh_.faceOwnerNeighbour(faceIndex)[0]) ;
        array<double,2> p2 = mesh_.cellCentre(mesh_.faceOwnerNeighbour(faceIndex)[1]) ;
        double d = (p1[0]-p2[0])*faceAreaVec[0] + (p1[1]-p2[1])*faceAreaVec[1] ;
        return abs(d);
            
    }

} ;
