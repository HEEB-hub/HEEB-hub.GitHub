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
    Foam::incompressible::lowDimensionMethod::ioLowDimension

Description
    Creating low-dimension mesh and field process of Low-Dimension Tool (LDT)

SourceFiles
    ioLowDimension.C

\*---------------------------------------------------------------------------*/

/**

@class ioLowDimension

@brief Derived class from lowDimensionMethod, realize the data Output function.
    
@author Jie Ren

@date 21 October 2019

@file ioLowDimension.H

*/

/*---------------------------------------------------------------------------*/

#ifndef ioLowDimension_H
#define ioLowDimension_H

#include "lowDimensionMethod.H"

namespace Foam
{
namespace incompressible
{
namespace lowDimensionMethod
{
/*---------------------------------------------------------------------------*\
                           Class ioLowDimension Declaration
\*---------------------------------------------------------------------------*/
class ioLowDimension 
:
    public lowDimensionMethod
{

protected:

 // Protected data
    
    // coefficients setting
    /** Switch to execute the IO process, read from the dict*/
    Switch       ldDataOutput_;
    /** Time name for folder construction, from the input of constructor*/ 
    word         timeName;
    /** generated according the low-dimension mesh properties, protected class member to share the ownership within the class */     
    dictionary   meshSet;        

public:
   
 // Constructors

        /**
        @brief the normal & promotion constructor
        the normal & promotion constructor, construct from the base class and time info
        */
        ioLowDimension 
        (
           const lowDimensionMethod&, word
        );

 //- Destructor
  virtual ~ioLowDimension() {} 

  // Member Functions

        /**
        @brief return bool to decide whether to execute the l-d mesh and field construction
        */
        bool dataOutput()     { return ldDataOutput_;}
     
        /**
        @brief to generate the Mesh dict for mesh construction 
        generate the blockMesh dictionary for mesh construction, update it in the protected data, but only construct the dictionary, no connection with the IOobject within it
        */
        void ldMeshDict();

        /**
        @brief to generate the mesh by calling the 'blockMesh' utility
        use the updated dictionary to generate an IOdictionary, then write it out in the file, call the "blockMesh" utility in system to execute the low-dimension mesh construction.
        */
        void ldMeshCreate();

        /**
        @brief to generate low-dimension scalar field based the existing results and the newly built low-dimension mesh
        use the low-dimension mesh generated in "ldMeshCreate()", generate the low dimension field for visualization and further use, including match the data in 'ld_Scalar_[]' of base class with the newly built low-dimension mesh.
        */
        void ldFieldCreate();

        /**
        @brief show the releated result information
        show the information releated to the output process.
        reloaded version of the abstract base class 
        */
        virtual void resultShow();

        /**
        @brief decide whether to update the settings
        default as false for the class, unnecessary function, never used.
        reloaded version of the abstract base class 
        */
        virtual bool update() {return false;}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace lowDimensionMethod
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif
