# pragma once

# include "PathConfig.hpp"
# include<iostream>
# include<vector>
# include "SparseLinearSystem.hpp"
# include "GaussSeidelSmoother.hpp"

using namespace std ;

/* *******************************************************************
    GaussSeidelSolver.hpp - This class contains the solve function that 
    runs the Gauss Seidel smoother. It also controls the number of 
    sweeps per check and the maximum number of smoother iterations.
******************************************************************* */

class GaussSeidelSolver: public GaussSeidelSmoother{

    private:
        
        SparseLinearSystem& sys_ ;
        int maxIter = 200;
        int sweepsPerCheck = 10 ;
        double resInitial = 0.0 ;
        double resFinal = 0.0 ;

    public:
        // Constructs the Gauss seidel solver object using the sparse linear system.
        GaussSeidelSolver(SparseLinearSystem& sys) :GaussSeidelSmoother(sys), sys_(sys){}
        
        // Function that runs the smoother. It also logs the residual values to the output file.
        void solve(double relTol, bool verbose, ostream& outInner, ostream& outFinal, int loop_count){

            resInitial = 0.0 ;
            resFinal = 0.0 ;

            int numSweepsCompleted = 0 ;
            sys_.calculateResidual() ;
            resInitial = sys_.twoNorm() ;

            while(numSweepsCompleted != maxIter){

                smooth(sweepsPerCheck) ;
                sys_.calculateResidual() ;
                resFinal = sys_.twoNorm() ;
                double ratio = resFinal/resInitial ;
                if(relTol > ratio){
                    numSweepsCompleted += sweepsPerCheck ;
                    break ;
                }

                numSweepsCompleted += sweepsPerCheck ;
                outInner << numSweepsCompleted << "," << resInitial << "," << resFinal << "," << ratio << endl;
            }
            outInner << endl ;

            if(verbose){
                cout << "Initial Residual: " << resInitial << ", Final Residual: " << resFinal << endl ;
            }     
            outFinal << loop_count << "," << resInitial << "," << resFinal << "," << resFinal/resInitial << endl ;

        }

} ;