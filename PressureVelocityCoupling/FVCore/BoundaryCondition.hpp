# pragma once

/* *******************************************************************
    BoundaryCondition.hpp -  This is an interface class for boundary 
    condition required to construct polymorphic functions for 
    applying various types of conditions at the boundary face for 
    diffusion and convection. 
******************************************************************* */

class BoundaryCondition{

    public:
    virtual void applyDiffusionBC() = 0 ;
    virtual void applyConvectionBC() = 0 ;
    virtual map<int, double>& boundaryValues() = 0 ;

    virtual ~BoundaryCondition() = default ;

} ;


