# pragma once

# include "PathConfig.hpp"

# include<iostream>
# include<vector>
# include<array>
# include<map>

# include "MeshFileReader2D.hpp"
# include "FaceAddressedMesh2D.hpp"

# include "SparseAddress.hpp"
# include "SparseMatrix.hpp"
# include "SparseLinearVector.hpp"
# include "SparseLinearSystem.hpp"
# include "GaussSeidelSmoother.hpp"
# include "GaussSeidelSolver.hpp"

# include "Convection.hpp"
# include "Diffusion.hpp"
# include "Ddt.hpp"
# include "Source.hpp"
# include "FVMatrix.hpp"
# include "FaceFluxField.hpp"
# include "ScalarProperty.hpp"
# include "BoundaryCondition.hpp"
# include "FixedValueBC.hpp"
# include "GradientBC.hpp"
# include "Patch.hpp"
# include "FieldFileReader.hpp"
# include "Field.hpp"
# include "PressureLaplacian.hpp"
# include "GaussGradient.hpp"
# include "UnderRelaxation.hpp"

using namespace std ;

/* *******************************************************************
    SIMPLE.hpp - This class solves the 2D incompressible lid-driven 
    cavity problem using a SIMPLE-based pressure-velocity coupling 
    algorithm.
******************************************************************* */

class SIMPLE{

        private:

        const SparseAddress addr_ ;
        const FaceAddressedMesh2D mesh_ ;

        FVMatrix Ux_ ;
        FVMatrix Uy_ ;
        FVMatrix P_ ;

        int n_ ;
        double re_ ;
        double l = 0.1 ;
        double nu ;
        double velocityRelativeTolerance = 0.0000001 ;
        double pressureRelativeTolerance = 0.05 ;
        double uxVelocity = 1 ;

        ScalarProperty gamma_ ;

        FaceFluxField flux_ ;

        FieldFileReader UxReader ;
        FieldFileReader UyReader ;
        FieldFileReader PReader ;

        Field UxInitF_ ;
        Field UyInitF_ ;
        Field PInitF_ ; 

        public:

        // Constructs the SIMPLE solver from mesh files, field files, grid size, and Reynolds number.
        SIMPLE(vector<string> meshFiles, vector<string> fieldFiles, int dimension, int re):addr_(meshFiles),
                                                                    mesh_(addr_.mesh()),
                                                                    Ux_(addr_),
                                                                    Uy_(addr_),
                                                                    P_(addr_),
                                                                    n_(dimension),
                                                                    re_(re),
                                                                    gamma_(mesh_, 0.0),
                                                                    flux_(addr_),
                                                                    UxReader(mesh_, Ux_, flux_, gamma_, fieldFiles[0]),
                                                                    UyReader(mesh_, Uy_, flux_, gamma_, fieldFiles[1]),
                                                                    PReader(mesh_, P_, flux_, gamma_, fieldFiles[2]),
                                                                    UxInitF_(UxReader.internal(), UxReader.boundary()),
                                                                    UyInitF_(UyReader.internal(), UyReader.boundary()),
                                                                    PInitF_(PReader.internal(), PReader.boundary()){}

        // Executes the full SIMPLE iteration loop, writes convergence diagnostics, and exports the final solution fields.
        void runSIMPLE(){

           /*
                Creating output streams to report necessary data during the SIMPLE loop run
           */

            PathConfig paths ;
            paths.createOutputDirectory(n_,re_) ;

            // Writes the basic case information such as mesh size, Reynolds number, and cavity length.
            ofstream outCase(paths.getOutputDirectory(n_, re_) / "case.txt") ;
            outCase << "Mesh Size: " << sqrt(mesh_.numCells()) << "x" << sqrt(mesh_.numCells()) << endl ;
            outCase << "Reynolds Number: " << re_ << endl ;
            outCase << "length: " << l << endl ;

            // Writes the pre-correction face flux at every face for each SIMPLE outer iteration.
            ofstream outFprePerFace(paths.getOutputDirectory(n_, re_) / "outputFprePerFace.csv") ;
            outFprePerFace << "Outer_Iter," ;
            for(int i=0; i<mesh_.numFaces(); ++i){ outFprePerFace << "Face_Index_" << i << "," ;}
            outFprePerFace << endl ;

            // Writes the pre-correction flux on boundary faces and their total sum for each outer iteration.
            ofstream outFpreBoundaryFace(paths.getOutputDirectory(n_, re_) / "outputFpreBoundaryFace.csv") ;
            outFpreBoundaryFace << "Outer_Iter," ;
            for(int i: mesh_.boundaryFaces()){ outFpreBoundaryFace << "Face_Index_" << i << "," ; }
            outFpreBoundaryFace << "sum," << endl ;

            // Writes the divergence of the pre-correction flux in every cell for each outer iteration.
            ofstream outDivFpre(paths.getOutputDirectory(n_, re_) / "outputDivFprePerCell.csv") ;
            outDivFpre << "Outer_Iter," ;
            for(int i=0; i<mesh_.numCells(); ++i){ outDivFpre << "Cell_" << i << "," ;}  
            outDivFpre << endl ;          

            // Writes the corrected face flux at every face for each SIMPLE outer iteration.
            ofstream outFcorrPerFace(paths.getOutputDirectory(n_, re_) / "outputFcorrPerFace.csv") ;
            outFcorrPerFace << "Outer_Iter," ;
            for(int i=0; i<mesh_.numFaces(); ++i){ outFcorrPerFace << "Face_Index_" << i << "," ;}
            outFcorrPerFace << endl ;

            // Writes the divergence of the corrected flux in every cell for each outer iteration.
            ofstream outDivFcorr(paths.getOutputDirectory(n_, re_) / "outputDivFcorrPerCell.csv") ;
            outDivFcorr << "Outer_Iter," ;
            for(int i=0; i<mesh_.numCells(); ++i){ outDivFcorr << "Cell_" << i << "," ; }  
            outDivFcorr << endl ;  

            // Writes global flux-balance diagnostics including divergence sums, maximum imbalances, and pressure-correction magnitude.
            ofstream outGlobalFluxBalance(paths.getOutputDirectory(n_, re_) / "outputGlobalFluxBalance.csv") ;
            outGlobalFluxBalance << "Outer_Iter,sumDivFpre,sumDivFcorr,maxAbsDivFpre,maxAbsDivFcorr,maxAbsPPrime,sumFPreAllFace,sumFCorrAllFace" << endl ;

            // Writes the final converged velocity and pressure values for every cell.
            ofstream outFinalFields(paths.getOutputDirectory(n_, re_) / "outputFinalFields.csv") ;
            outFinalFields << "Cell Number,UxFinal,UyFinal,PFinal" << endl ;

            // Writes the maximum, minimum, and absolute range of all fields at each outer iteration.
            ofstream outFieldsRange(paths.getOutputDirectory(n_, re_) / "outputFieldsRange.csv") ;
            outFieldsRange << "Outer_Iter,MaxUx,MinUx,MaxAbsUx,MinAbsUx,MaxUy,MinUy,MaxAbsUy,MinAbsUy,MaxP,MinP,MaxAbsP,MinAbsP" << endl ;

            // Writes the maximum change in Ux, Uy, and pressure between successive outer iterations.
            ofstream outMaxFieldChanges(paths.getOutputDirectory(n_, re_) / "outputMaxFieldChanges.csv") ;
            outMaxFieldChanges << "OuterIter,Ux,Uy,P" << endl ;

            // Writes the final assembled matrix coefficients for the Ux, Uy, and pressure systems.
            ofstream outCoefficients(paths.getOutputDirectory(n_, re_) / "outputFinalCoeffcients.csv") ;
            outCoefficients << "i,j,Ux_A,Uy_A,P_A" << endl ;

            // Writes the vertical centreline Ux values for each outer iteration.
            ofstream outMidLineUx(paths.getOutputDirectory(n_, re_) / "outputMidLineUx.csv") ;
            outMidLineUx << "Outer_Iter," ;
            int cellsPerRow = sqrt(mesh_.numCells()) ;
            int midLine = cellsPerRow/2 ;
            for(int i=midLine; i<mesh_.numCells(); i+=cellsPerRow){ outMidLineUx << "Cell_" << i << "," ;}
            outMidLineUx << endl ;

            // Writes the vertical centreline Uy values for each outer iteration.
            ofstream outMidLineUy(paths.getOutputDirectory(n_, re_) / "outputMidLineUy.csv") ;
            outMidLineUy << "Outer_Iter," ;
            cellsPerRow = sqrt(mesh_.numCells()) ;
            midLine = cellsPerRow/2 ;
            for(int i=midLine; i<mesh_.numCells(); i+=cellsPerRow){ outMidLineUy << "Cell_" << i << "," ;}
            outMidLineUy << endl ;

            // Writes the vertical centreline pressure values for each outer iteration.
            ofstream outMidLineP(paths.getOutputDirectory(n_, re_) / "outputmidLineP.csv") ;
            outMidLineP << "Outer_Iter," ;
            cellsPerRow = sqrt(mesh_.numCells()) ;
            midLine = cellsPerRow/2 ;
            for(int i=midLine; i<mesh_.numCells(); i+=cellsPerRow){ outMidLineP << "Cell_" << i << "," ;}
            outMidLineP << endl ;

            // Writes the inner Gauss-Seidel residual history for the Ux solve.
            ofstream outInnerResidualUx(paths.getOutputDirectory(n_, re_) / "innerResidualUx.csv") ;

            // Writes the inner Gauss-Seidel residual history for the Uy solve.
            ofstream outInnerResidualUy(paths.getOutputDirectory(n_, re_) / "innerResidualUy.csv") ;

            // Writes the inner Gauss-Seidel residual history for the pressure solve.
            ofstream outInnerResidualP(paths.getOutputDirectory(n_, re_) / "innerResidualP.csv") ;

            outInnerResidualUx << "Gauss Seidel Sweeps,Initial Residual,Final Residual,Ratio" << endl ;
            outInnerResidualUy << "Gauss Seidel Sweeps,Initial Residual,Final Residual,Ratio" << endl ;
            outInnerResidualP << "Gauss Seidel Sweeps,Initial Residual,Final Residual,Ratio" << endl ;

            // Writes the final residual summary for the Ux equation at each SIMPLE loop iteration.
            ofstream outFinalResidualUx(paths.getOutputDirectory(n_, re_) / "finalResidualUx.csv") ;

            // Writes the final residual summary for the Uy equation at each SIMPLE loop iteration.
            ofstream outFinalResidualUy(paths.getOutputDirectory(n_, re_) / "finalResidualUy.csv") ;

            // Writes the final residual summary for the pressure equation at each SIMPLE loop iteration.
            ofstream outFinalResidualP(paths.getOutputDirectory(n_, re_) / "finalResidualP.csv") ;
            
            outFinalResidualUx << "Loop Count,Initial Residual,Final Residual,Ratio" << endl ;
            outFinalResidualUy << "Loop Count,Initial Residual,Final Residual,Ratio" << endl ;
            outFinalResidualP << "Loop Count,Initial Residual,Final Residual,Ratio" << endl ;

           /*
                Initializing the kinematic viscocity object. Assuming only x component of lid velocity is initially non zero. 
           */ 

            nu = uxVelocity*l/re_ ;
            for(int i=0; i<mesh_.numCells(); ++i){ gamma_(i) = nu ;}

            /*
                Clearing all the matrices and initializing the solution variables as per the field
            */

            Ux_.clear() ;
            Uy_.clear() ;
            P_.clear() ;

            Ux_.initializeSolution(UxInitF_.internal()) ;
            Uy_.initializeSolution(UyInitF_.internal()) ;
            P_.initializeSolution(PInitF_.internal()) ;

            /*
                Initial flux calculation and under-relaxation objects
            */

            flux_.updateInitialFlux(UxInitF_, UyInitF_) ;

            Convection UxConvection(Ux_, flux_, UxInitF_) ;
            Diffusion UxDiffusion(Ux_, mesh_, UxInitF_, gamma_) ;

            Convection UyConvection(Uy_, flux_, UyInitF_) ;
            Diffusion UyDiffusion(Uy_, mesh_, UyInitF_, gamma_) ;

            PressureLaplacian pressureLaplacian(P_, Ux_, Uy_, mesh_, flux_, PInitF_) ;
            GaussGradient pGrad(P_, PInitF_, mesh_) ;

            UnderRelaxation underRelaxation(0.8, 0.2) ;

            int count = 0 ;
            SparseLinearVector pPrime(mesh_.numCells()) ;

            double maxFieldChangeUx = -1e8 ;
            double maxFieldChangeUy = -1e8 ;
            double maxFieldChangeP = -1e8 ;
            double maxAbsDivFlux = 0.0 ;
            double tolUpdateU = 1e-7 ;
            double tolUpdateP = 5e-5 ;
            double tolCont = 1e-16 ;
            bool fieldsAreChanging = true ;

            double maxAbsDivFluxCorr = 0.0 ;
            double maxAbsDivFluxPre = 0.0 ;
            double sumDivFCorr = 0.0 ;
            double sumDivFPre = 0.0 ;
            double maxAbsPPrime = 0.0 ;
            double sumFPreAllFace = 0.0 ;
            double sumFCorrAllFace = 0.0 ;
            int maxIter = 1500 ;

        
            /*
                SIMPLE LOOP START
            */

            while(fieldsAreChanging){

                cout << "Loop Count: " << count << endl ; 
                /*
                    Clearing the matrices before each iteration
                */
                Ux_.clear() ;
                Uy_.clear() ;
                P_.clear() ;

                SparseLinearVector UxFieldChange_ = Ux_.sys().X() ;
                SparseLinearVector UyFieldChange_ = Uy_.sys().X() ;
                SparseLinearVector PFieldChange_ = P_.sys().X() ;

                /*
                    Ux and Uy momentum predictor step with solution to find velocity field
                */

                UxConvection.assemble() ;
                UxDiffusion.assemble() ;
                underRelaxation.underRelaxMomentum(Ux_) ;
                
                UyConvection.assemble() ;
                UyDiffusion.assemble() ;
                underRelaxation.underRelaxMomentum(Uy_) ;
                
                Ux_.solve(velocityRelativeTolerance, true, outInnerResidualUx, outFinalResidualUx, count) ;
                Uy_.solve(velocityRelativeTolerance, true, outInnerResidualUy, outFinalResidualUy, count) ;

                /*
                    Updating the pre flux (Fpre) before assembling pressure equation assembly and correction
                */

                flux_.updateFlux(Ux_.sys().X(), Uy_.sys().X()) ;
                flux_.updateBoundaryFlux(UxInitF_, UyInitF_) ;

                outFprePerFace << count+1 << "," ;
                for(int i=0; i<mesh_.numFaces(); ++i){
                    outFprePerFace << flux_(i) << "," ;
                    sumFPreAllFace += flux_(i) ;
                }
                outFprePerFace << endl ;

                outFpreBoundaryFace << count+1 << "," ;
                double sum = 0 ;
                for(int i: mesh_.boundaryFaces()){
                    outFpreBoundaryFace << flux_(i) << "," ;
                    sum += flux_(i) ;
                }
                outFpreBoundaryFace << sum << "," ;
                outFpreBoundaryFace << endl ;

                outDivFpre << count+1 << "," ;
                for(int i=0; i<mesh_.numCells(); ++i){
                    vector<int> faces = mesh_.cellFaces(i) ;
                    double sum = 0 ;
                    for(int faceIndex: faces){
                        int cellOwner = mesh_.faceOwnerNeighbour(faceIndex)[0] ;
                        if(cellOwner == i){
                            sum += flux_(faceIndex) ;
                        }
                        else{
                            sum -= flux_(faceIndex) ;
                        }
                    }
                    outDivFpre << sum << "," ;
                    sumDivFPre += sum ;
                    if(maxAbsDivFluxPre < abs(sum)){
                        maxAbsDivFluxPre = abs(sum) ;
                    }
                }
                outDivFpre << endl ;
                maxAbsDivFluxPre = 0.0 ;

                /*
                    assembling and solving pressure laplacian
                */

                pressureLaplacian.assemble() ;
                SparseLinearVector pOld_ = P_.sys().X() ;
                
                for(int i=0; i<P_.sys().size(); ++i){
                    P_.sys().X()(i) = pPrime(i) ;
                }

                P_.solve(pressureRelativeTolerance, true, outInnerResidualP, outFinalResidualP, count) ;
                maxAbsPPrime = P_.sys().X().maxAbsVal() ;

                for(int i=0; i<P_.sys().size(); ++i){
                    pPrime(i) = P_.sys().X()(i) ;
                }                
                
                /*
                    Correcting flux and velocity.
                */
                
                pGrad.computeGradient() ;
                for(int i=0; i<mesh_.numCells(); ++i){
                    Ux_.sys().X()(i) = Ux_.sys().X()(i) - (1/Ux_.sys().A()(i,i))*pGrad.X()(i) ;
                    Uy_.sys().X()(i) = Uy_.sys().X()(i) - (1/Uy_.sys().A()(i,i))*pGrad.Y()(i) ;
                }
                flux_.correctFlux(Ux_, Uy_, P_) ;
                underRelaxation.underRelaxPressure(P_, pOld_) ;

                UxFieldChange_ = UxFieldChange_ - Ux_.sys().X() ;
                UyFieldChange_ = UyFieldChange_ - Uy_.sys().X() ;
                PFieldChange_ = PFieldChange_ - P_.sys().X() ;

                maxFieldChangeUx = UxFieldChange_.maxAbsVal() ;
                maxFieldChangeUy = UyFieldChange_.maxAbsVal() ;
                maxFieldChangeP = PFieldChange_.maxAbsVal() ;

                outFcorrPerFace << count+1 << "," ;
                for(int i=0; i<mesh_.numFaces(); ++i){
                    outFcorrPerFace << flux_(i) << "," ;
                    sumFCorrAllFace += flux_(i) ;
                }
                outFcorrPerFace << endl ;

                outDivFcorr << count+1 << "," ;
                for(int i=0; i<mesh_.numCells(); ++i){
                    vector<int> faces = mesh_.cellFaces(i) ;
                    sum = 0.0 ;
                    for(int faceIndex: faces){
                        int cellOwner = mesh_.faceOwnerNeighbour(faceIndex)[0] ;
                        if(cellOwner == i){
                            sum += flux_(faceIndex) ;
                        }
                        else{
                            sum -= flux_(faceIndex) ;
                        }
                    }
                    outDivFcorr << sum << "," ;
                    sumDivFCorr += sum ;
                    if(maxAbsDivFlux < abs(sum)){
                        maxAbsDivFlux = abs(sum) ;
                    }

                }
                outDivFcorr << endl ;

                cout << "maxAbsDivFlux: " << maxAbsDivFlux << " tolCont: " << tolCont << endl ;
                cout << "maxFieldChangeUx: " << maxFieldChangeUx << endl ;
                cout << "maxFieldChangeUy: " << maxFieldChangeUy << endl ;
                cout << "maxFieldChangeP: " << maxFieldChangeP << endl ;
                outMaxFieldChanges << count + 1 << "," << maxFieldChangeUx << "," << maxFieldChangeUy << "," << maxFieldChangeP << endl ;

                outMidLineUx << count+1 << "," ;
                int cellsPerRow = sqrt(mesh_.numCells()) ;
                int midLine = cellsPerRow/2 ;
                for(int i=midLine; i<mesh_.numCells(); i+=cellsPerRow){outMidLineUx << Ux_.sys().X()(i) << "," ;}
                outMidLineUx << endl ;

                outMidLineUy << count+1 << "," ;
                cellsPerRow = sqrt(mesh_.numCells()) ;
                midLine = cellsPerRow/2 ;
                for(int i=midLine; i<mesh_.numCells(); i+=cellsPerRow){outMidLineUy << Uy_.sys().X()(i) << "," ;}
                outMidLineUy << endl ;

                outMidLineP << count+1 << "," ;
                cellsPerRow = sqrt(mesh_.numCells()) ;
                midLine = cellsPerRow/2 ;
                for(int i=midLine; i<mesh_.numCells(); i+=cellsPerRow){ outMidLineP << P_.sys().X()(i) << "," ;}
                outMidLineP << endl ;

                outFieldsRange << count+1 << "," 
                               << Ux_.sys().X().maxVal() << "," << Ux_.sys().X().minVal() << "," << Ux_.sys().X().maxAbsVal() << "," << Ux_.sys().X().minAbsVal() << ","
                               << Uy_.sys().X().maxVal() << "," << Uy_.sys().X().minVal() << "," << Uy_.sys().X().maxAbsVal() << "," << Uy_.sys().X().minAbsVal() << ","
                               << P_.sys().X().maxVal() << "," << P_.sys().X().minVal() << "," << P_.sys().X().maxAbsVal() << "," << P_.sys().X().minAbsVal() << ","
                               << endl ;

                fieldsAreChanging = !((maxFieldChangeUx < tolUpdateU) && (maxFieldChangeUy < tolUpdateU) && (maxFieldChangeP < tolUpdateP) && (maxAbsDivFlux < tolCont)) && (count < maxIter);                
                maxAbsDivFluxCorr = maxAbsDivFlux ;
                ++count ;


                outGlobalFluxBalance << count << "," << sumDivFPre << "," << sumDivFCorr << "," << maxAbsDivFluxPre << "," << maxAbsDivFluxCorr << "," << maxAbsPPrime << ","
                                     << sumFPreAllFace << "," << sumFCorrAllFace << endl ;

                maxAbsDivFlux = 0.0 ;
                maxAbsDivFluxCorr = 0.0 ;
                maxAbsDivFluxPre = 0.0 ;
                sumDivFCorr = 0.0 ;
                sumDivFPre = 0.0 ;
                maxAbsPPrime = 0.0 ;
                sumFPreAllFace = 0.0 ;
                sumFCorrAllFace = 0.0 ;

                cout << "***************************************************************************************" << endl ;
            }

            for(int i=0; i< mesh_.numCells(); ++i){ 
                outCoefficients << i << "," << i << "," << Ux_.sys().A()(i,i) << "," << Uy_.sys().A()(i,i) << "," << P_.sys().A()(i,i) << endl ;
                outFinalFields << i << "," << Ux_.sys().X()(i) << "," << Uy_.sys().X()(i) << "," << P_.sys().X()(i) << endl ;
                for(int j = addr_.rowArray()[i]; j < addr_.rowArray()[i+1]; ++j){
                    outCoefficients << i << "," << addr_.columnIndexArray()[j] << "," 
                    << Ux_.sys().A()(i,addr_.columnIndexArray()[j]) << ","
                    << Uy_.sys().A()(i,addr_.columnIndexArray()[j]) << ","
                    << P_.sys().A()(i,addr_.columnIndexArray()[j]) << endl ;
                }
            }
            ofstream outVTKFile(paths.getOutputDirectory(n_, re_) / "SimulationFile.vtk") ;
            writeVTKRectilinearFile(outVTKFile) ;
        }
    
        // Writes the final pressure and velocity fields to a rectilinear VTK file for visualization.
        void writeVTKRectilinearFile(ostream& out){

            if(!out){throw runtime_error("Could not open VTK file for writing") ;}

            out << "# vtk DataFile Version 3.0" << endl ;
            out << "Lid-driven cavity result" << endl ;
            out << "ASCII" << endl ;
            out << "DATASET RECTILINEAR_GRID" << endl ;
            out << "DIMENSIONS " << (static_cast<int>(sqrt(mesh_.numCells())) + 1) << " " << (static_cast<int>(sqrt(mesh_.numCells())) + 1) << " 1" << endl << endl ;

            out << "X_COORDINATES " << (static_cast<int>(sqrt(mesh_.numCells())) + 1) << " float" << endl ;
            for (int i = 0; i <= static_cast<int>(sqrt(mesh_.numCells())); ++i) {
                double x = l * static_cast<double>(i) / static_cast<double>(static_cast<int>(sqrt(mesh_.numCells())));
                out << setprecision(10) << x << (i < static_cast<int>(sqrt(mesh_.numCells())) ? " " : "\n");
            }
            out << endl ;

            out << "Y_COORDINATES " << (static_cast<int>(sqrt(mesh_.numCells())) + 1) << " float" << endl ;
            for (int j = 0; j <= static_cast<int>(sqrt(mesh_.numCells())); ++j) {
                double y = l * static_cast<double>(j) / static_cast<double>(static_cast<int>(sqrt(mesh_.numCells())));
                out << setprecision(10) << y << (j < static_cast<int>(sqrt(mesh_.numCells())) ? " " : "\n");
            }
            out << endl ;

            out << "Z_COORDINATES 1 float" << endl ;
            out << "0.0" << endl << endl ;

            out << "CELL_DATA " << mesh_.numCells() << endl ;

            out << "SCALARS pressure float 1\n";
            out << "LOOKUP_TABLE default\n";
            for (int c = 0; c < mesh_.numCells(); ++c) {
                out << setprecision(10) << P_.sys().X()(c) << endl ;
            }
            out << endl ;

            out << "VECTORS velocity float\n";
            for (int c = 0; c < mesh_.numCells(); ++c) {
                out << std::setprecision(10)
                    << Ux_.sys().X()(c) << " " << Uy_.sys().X()(c) << " 0.0" << endl ;
            }

        }
   
    } ;