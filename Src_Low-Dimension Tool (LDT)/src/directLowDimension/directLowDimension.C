/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Namespace
    Foam::incompressible::lowDimensionMethod::directLowDimension

Description
    Namespace for low-dimensional tools(LDT), which is a post-process procedure.

Class
    Foam::incompressible::lowDimensionMethod::directLowDimension

Description
    direct low-dimension process of Low-Dimension Tool (LDT)

SourceFiles
    directLowDimension.C

\*---------------------------------------------------------------------------*/

#include "directLowDimension.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace lowDimensionMethod
{
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
directLowDimension::directLowDimension
(  
   const volScalarField& C,
   const word& type
)
:
    lowDimensionMethod(C,type)
//though same data construction as ABC, but to updata it according to dict set!
{

    dimensionedVector ld_Dimension_
    (
        dimensioned<vector>::lookupOrAddToDict
        (
            "lowDimensionSet",
            coeffDict_,
            pTraits<vector>::one
        )
    );
      this->ld_Dimension_=ld_Dimension_;
      fieldNums();
      setPoints();
      Info<<"Derived class: directLowDimension construction finished!"<<endl<<endl;
};


directLowDimension::directLowDimension
(  
   const lowDimensionMethod& LDM
)
:   lowDimensionMethod(LDM) 
{
Info<<"Calling directLowDimension lowDimensionMethod& LDM:"<<nl;

};



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

 //- decide whether to update the ld-settings, then do it and return a flag

     bool directLowDimension::update()

     { 
 	
        Info<<"Direct method need no update!"<<endl<<endl;
        return false;   //direct method need no update
     } 


 //- show the process info and results
     void directLowDimension::resultShow()
     {
        Info<<"Low-Dimension Tool (LDT) works successfully!"
            <<endl <<endl 
            <<"Results as below:"
            <<endl <<endl  
            <<"The computational domain has been divided into " 
            <<ld_fieldNums_
            <<" zones, with " 
            <<ld_Dimension_.value()
            <<" pieces in each (X,Y,Z) coordinate!"
            <<endl<<endl 
            <<"The total L-D ERROR equals "
            <<ld_Error_.value()*100<<"%"
            <<endl<<endl 
            <<"Congrats! Direct Low-Dimension finallized!"<<endl;
      }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace lowDimensionMethod
} // End namespace incompressible
} // End namespace Foam

