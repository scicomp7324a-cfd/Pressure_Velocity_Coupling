# pragma once

# include<iostream>
# include<vector>
# include<array>
# include "Patch.hpp"
# include "BoundaryCondition.hpp"
# include "SparseLinearVector.hpp"

using namespace std ;

/* *******************************************************************
    Field.hpp -  This class contains access functions to 
    manage initial internal and boundary values of the velocity
    and pressure field. 
******************************************************************* */

class Field{

    private:
    SparseLinearVector internal_ ;
    vector<unique_ptr<BoundaryCondition>> boundary_ ;

    public:
    // Constructs the field using the internal and boundary field values
    Field(SparseLinearVector& internal, vector<unique_ptr<BoundaryCondition>>& boundary):internal_(internal),
                                                                          boundary_(std::move(boundary)){}

    // Non-constant access to the initial internal field
    SparseLinearVector& internal() {
        return internal_ ;
    }

    // Non-constant access to the boundary faces using unique pointer to the specific boundary condition based on patch id
    vector<unique_ptr<BoundaryCondition>>& boundary() {
        return boundary_ ;
    }

} ;
