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
    the Free Software Foundation, either version 6 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Namespace
    Foam::incompressible::lowDimensionMethod

Description
    Namespace for low-dimensional tools(LDT), which is a post-process procedure.

Class
    Foam::incompressible::lowDimensionMethod::directLowDimension

Description
    direct low-dimension process of Low-Dimension Tool (LDT)

SourceFiles
    directLowDimension.C

\*---------------------------------------------------------------------------*/

/**

@class directLowDimension

@brief Derived class from lowDimensionMethod, for direct low-dimension process.
    
@author Jie Ren

@date 14 October 2019

@file directLowDimension.H

*/

/*---------------------------------------------------------------------------*/

#ifndef directLowDimension_H
#define directLowDimension_H

#include "lowDimensionMethod.H"

namespace Foam
{
namespace incompressible
{
namespace lowDimensionMethod
{
/*---------------------------------------------------------------------------*\
                           Class directLowDimension Declaration
\*---------------------------------------------------------------------------*/
class directLowDimension 
:
    public lowDimensionMethod
{

private:
 // Private data
    
 // none data, totally based on ABC

public:

 // Constructors

        /**
        @brief the normal constructor
        the normal constructor, construct from a given ScalarField, model type default as direct
        */
        directLowDimension
        (
           const volScalarField& C,
           const word& type="direct"
        );
        /**
        @brief the type promotion constructor
        the type promotion constructor, construct from its base class-lowDimensionMethod, make the conversion valid
        */
 	directLowDimension
        (
          const lowDimensionMethod&
        );

 //- Destructor
  virtual ~directLowDimension()   {}

  // Member Functions

        /**
        @brief decide whether to update the settings
        decide whether to update the low-dimension settings, default as false for the class
        reloaded version of the abstract base class 
        */
        virtual bool update();

        /**
        @brief show the releated result information
        show the process information and results 
        reloaded version of the abstract base class 
        */
        virtual void resultShow();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace lowDimensionMethod
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
