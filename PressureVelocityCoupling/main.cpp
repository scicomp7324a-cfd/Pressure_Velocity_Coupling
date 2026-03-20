# include "PathConfig.hpp"

# include<iostream>
# include<vector>
# include<array>
# include<map>
# include<filesystem>

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
# include "SIMPLE.hpp"


using namespace std ;

/* *******************************************************************
    Main driver for running the pressure-velocity coupling solver 
    over multiple mesh sizes and Reynolds numbers.
******************************************************************* */

int main(){

    // Creates and validates all required input/output paths before starting the simulations.
    PathConfig paths ;

    // Stops execution if any required path or output directory setup fails.
    try{
        paths.validate() ;
        paths.ensureOutputDir() ;
    }
    catch(const exception& e){
        cerr << "Path error: " << e.what() << endl ;
        return 1 ;
    }

    // Stores the mesh and field file lists used to initialise each simulation case.
    vector<string> meshFiles, fieldFiles ;

    // Defines the field files for the initial Ux, Uy, and pressure data.
    fieldFiles = {paths.UxFieldFile(), paths.UyFieldFile(), paths.PFieldFile()} ;
    
    // Lists the mesh sizes used in the study.
    vector<int> nArr = {20,40,60,80} ;

    // Lists the Reynolds numbers used in the study.
    vector<int> reArr = {100, 400, 1000} ;

    // Loops over all mesh-size and Reynolds-number combinations and runs the SIMPLE solver for each case.
    for(int n: nArr){
        for(int re: reArr){
            // Prints the output directory for the current case.
            cout << paths.getOutputDirectory(n,re) << endl ;

            // Retrieves the mesh files for the current mesh configuration.
            meshFiles = paths.getMeshFilenames(n) ;

            // Constructs the SIMPLE solver for the current case and executes the full simulation.
            SIMPLE(meshFiles, fieldFiles, n, re).runSIMPLE() ;
        }
    }
}








