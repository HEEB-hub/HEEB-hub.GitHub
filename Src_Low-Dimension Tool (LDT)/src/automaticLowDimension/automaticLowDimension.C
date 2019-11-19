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
    Foam::incompressible::automaticLowDimension

Description
    Namespace for low-dimensional tools(LDT), which is a post-process procedure.

Class
    Foam::incompressible::automaticLowDimension

Description
     automatic low-dimension process of Low-Dimension Tool (LDT)

SourceFiles
    automaticLowDimension.C

\*---------------------------------------------------------------------------*/

#include "automaticLowDimension.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace lowDimensionMethod
{
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
automaticLowDimension::automaticLowDimension
(  
   const volScalarField& C,
   const word& type
)
:
    directLowDimension(C,type),

    ld_Initial_set_
    (
        dimensioned<vector>::lookupOrAddToDict
        (
            "lowDimensionSet",
            coeffDict_,
            vector::one
        )
    ),

    ld_LimitError_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "errorLimit",
            coeffDict_,
            0.1
        )
    ),

    ld_updateRatio_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "updateRatio",
            coeffDict_,
            1
        )
    ),

    ld_maxStep_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "maxStep",
            coeffDict_,
            1000
        )
    ),
    Step(1),
    Valid_step(1)
{
  Info<<"Derived class: automaticLowDimension construction finished!"<<endl<<endl;
};


automaticLowDimension::automaticLowDimension
(
    const directLowDimension& LDM
):
    directLowDimension(LDM),

    ld_Initial_set_(pTraits<vector>::one),

    ld_LimitError_(0.1),

    ld_updateRatio_(1),

    ld_maxStep_(1000),
 
    Step(0),

    Valid_step(0)
{
Info<<"Calling automaticLowDimension directLowDimension& LDM:"<<nl;
};


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

 //- decide whether to update the ld-settings, then do it and return a flag

bool automaticLowDimension::update()   //automatic method need update

{ 
     if (ld_Error_<ld_LimitError_)
        { 
          Info<<"Finish step "<<Step<<nl; 
          Info<<"Relative error is now less than the creterion! : "
              <<ld_LimitError_.value()*100<<"%"<<endl
              <<"Iteration finished!"<<endl<<endl<<endl<<endl<<endl;
        return false;
        }

     else if (Step>=ld_maxStep_.value())
        { 
          Info<<"Finish step "<<Step<<nl;
          Info<<"Reach the max step limit! "
              <<"Iteration finished!"<<endl<<endl<<endl<<endl<<endl;
        return false;
        }

     else
     {  
        int max=0;
        for (int i=1;i<3;i++)
        if  (fieldVector.value()[max]<fieldVector.value()[i])
              max=i;
       const vector normalsize = fieldVector.value()/fieldVector.value()[max];
       const vector oldsize=ld_Dimension_.value();   
      do
      {
        Info<<"Finish step "<<Step<<", getting a new lowDimensionSet... "<<nl
            <<"Checking if it is available..."<<endl<<nl;
       vector temp(normalsize*Step*ld_updateRatio_.value());         
       ld_Dimension_.value() = vector(                  
                               int (0.5+temp.x()+ld_Initial_set_.value().x()),
                               int (0.5+temp.y()+ld_Initial_set_.value().y()),
                               int (0.5+temp.z()+ld_Initial_set_.value().z())
                               );
        Step++;
        if (Step>ld_maxStep_.value())
        { Info<<"No new lowDimensionSet available"<<nl
              <<"Reach the max step limit! "
              <<"Iteration finished!"<<endl<<endl<<endl<<endl<<endl;
        return false;
        }
      }
      while (ld_Dimension_.value()[0]==oldsize[0] &&
             ld_Dimension_.value()[1]==oldsize[1] &&
             ld_Dimension_.value()[2]==oldsize[2] );
      Valid_step++;
       Info<<"The lowDimensionSet updated successfully!"<<endl
              <<"Going to start next iteration!"<<endl<<endl<<endl<<endl<<endl;
       fieldNums();
       setPoints();
        return true;
      }
}

 //- show the process info and results

     void automaticLowDimension::resultShow()
     {
        Info<<"Low-Dimension Tool (LDT) works successfully!"
            <<endl<<endl
            <<"Results as below:"
            <<endl<<endl
            <<"The computational domain has been divided into " 
            <<ld_fieldNums_
            <<" zones, with " 
            <<ld_Dimension_.value()
            <<"pieces in each (X,Y,Z) coordinate!"
            <<endl<<endl
            <<"The total L-D ERROR equals "
            <<ld_Error_.value()*100<<"%"
            <<endl<<endl
            <<"Total iteration: "<<Step
            <<" steps, "
            <<"with "<<Valid_step<<" Steps valid."<<endl<<endl
            <<"Congrats! Automatic Low-Dimension finallized!"<<endl<<endl;
      }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace lowDimensionMethod
} // End namespace incompressible
} // End namespace Foam

