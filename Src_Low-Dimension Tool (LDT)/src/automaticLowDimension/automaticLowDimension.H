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
    Foam::incompressible::lowDimensionMethod::automaticLowDimension

Description
    automatic low-dimension process of Low-Dimension Tool (LDT)

SourceFiles
    automaticLowDimension.C

\*---------------------------------------------------------------------------*/

/**

@class automaticLowDimension

@brief Derived class from directLowDimension, for auto-updated low-dimension process.
    
@author Jie Ren

@date 20 October 2019

@file automaticLowDimension.H

*/

/*---------------------------------------------------------------------------*/

#ifndef automaticLowDimension_H
#define automaticLowDimension_H

#include "directLowDimension.H"

namespace Foam
{
namespace incompressible
{
namespace lowDimensionMethod
{
/*---------------------------------------------------------------------------*\
                           Class automaticLowDimension Declaration
\*---------------------------------------------------------------------------*/
class automaticLowDimension 
:
    public directLowDimension 
{

protected:

 // Protected data
    
    // Model coefficients setting
    /**initial low-dimension setting of AUTO mode*/
    const dimensionedVector   ld_Initial_set_; 
    /**Error limitation of the low-dimension process*/
    const dimensionedScalar   ld_LimitError_;  
    /**update ratio for update process*/
    const dimensionedScalar   ld_updateRatio_;
    /**max step limitation for low-dimension update process*/ 
    const dimensionedScalar   ld_maxStep_;
    /**number to record the number of update steps*/    
    int                       Step;  
    int 		      Valid_step;         
  
public:

 // Constructors

         /**
        @brief the normal constructor
        the normal constructor, construct from a given ScalarField, model type default as automatic
        */
        automaticLowDimension
        (
           const volScalarField& C,
           const word& type="automatic"
        );
        /**
        @brief the type promotion constructor
        the type promotion constructor, construct from its base class-directLowDimension, make the conversion valid
        */
 	automaticLowDimension
        (
           const directLowDimension&
        );


 //- Destructor
  virtual ~automaticLowDimension()   {}

  // Member Functions

        /**
        @brief decide whether to update the settings
        decide whether to update the low-dimension settings, by using the parameters read from the dictionary, we can control the update speed and set threshold for it.
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
