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
    Foam::incompressible::lowDimensionMethod::lowDimensionMethod

Description
    Abstract base class of low-dimensional tools(LDT).

SourceFiles
    lowDimensionMethod.C

\*---------------------------------------------------------------------------*/

/**

@class lowDimensionMethod

@brief Abstract base class for low-dimensional tools (LDT) methods.
    
@author Jie Ren

@date 31 October 2019

@file lowDimensionMethod.H

*/

/*---------------------------------------------------------------------------*/
#ifndef lowDimensionMethod_H
#define lowDimensionMethod_H


#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "surfaceFieldsFwd.H"
#include "surfaceFields.H"
#include "Switch.H"
#include "Vector.H"
#include <vector>
#include <algorithm>
//#include "autoPtr.H"

namespace Foam
{
// Forward declarations
class fvMesh;
namespace incompressible
{
namespace lowDimensionMethod
{
/*---------------------------------------------------------------------------*\
                           Class lowDimensionMethod Declaration
\*---------------------------------------------------------------------------*/
class lowDimensionMethod 
:
    public IOdictionary
{

protected:

 // protected data
   

     /** Mesh read from original case*/
     const fvMesh& mesh_;
 
     /** Scalar Field read from original case*/
     const volScalarField& C_;

     /** Model coefficients dictionary*/
     dictionary coeffDict_;

     /** Switch to print the model coeffs*/ 
     Switch printCoeffs_;

     /** field informaion: the start point offset from the original point*/
      dimensionedVector      minVector;

     /** field informaion: the field range of the domain presented as 3 components*/
      dimensionedVector      fieldVector;

     /** Low-Dimension data: low-Dimension quantities at present*/
      dimensionedVector      ld_Dimension_;

     /** Low-Dimension data: relative error result*/
      dimensionedScalar      ld_Error_;

     /** Low-Dimension data: array of the new low-dimension scalar data*/
      dimensionedScalar*     ld_Scalar_;

     /** Low-Dimension data: array of the new low-dimension sub-zone error*/
      dimensionedScalar*     ld_SubError_;

     /** Low-Dimension data: multiple vector to save the cut point vector of each direction*/
      std::vector<std::vector<double>> 
			     ld_cutpoint_;

     /** Low-Dimension data: total low-dimension field/zone numbers*/
      int                    ld_fieldNums_;

  // Protected Member Functions

    /**
      @brief calculate the total fieldNums
      function to calculate the total fieldNums under the given low-dimension settings
    */  
   inline void fieldNums()
    	 { ld_fieldNums_=
                      int(
                      ld_Dimension_.value().x()*
                      ld_Dimension_.value().y()*
                      ld_Dimension_.value().z());
      delete [] ld_Scalar_; 
      delete [] ld_SubError_;
      ld_Scalar_  =new dimensionedScalar[ld_fieldNums_];//update the L-D size
      ld_SubError_=new dimensionedScalar[ld_fieldNums_];//update the L-D size
   	  };

    /**
      @brief calculate the field range
      calculate the fieldRange vector(minVector & fieldVector) of the mesh domain
    */
    void fieldRange();
 
    /**      
      @brief set the uniform cut point
      set the cut point for the low-dimension conversion, only used for the uniform cut strategy 
    */
    void setPoints();
	
  
private:

 // Private Member Functions

     /**
        Disallow default bitwise assignment
     */
     void operator=(const lowDimensionMethod&);

public:
     
 // Constructors

        /**
        @brief the normal constructor
        the normal constructor, construct from a given ScalarField, abstract base class, no specific object
        */
        lowDimensionMethod
        (
            const volScalarField& C,
            const word& type
        );

        /**
        @brief the copy constructor
        the copy constructor, construct from same base class, deep copy of the array.
        */
        lowDimensionMethod(const lowDimensionMethod&);


 //- Destructor
  /**
     delete the dynamic array pointer. 
  */
  virtual ~lowDimensionMethod()   { delete [] ld_Scalar_; delete [] ld_SubError_; }

  // Member Functions

        /**
        @brief show setting coeffs
        print Coeffs of subDict based on dictionary settings
       */
        void printCoeffs() 
        { if (printCoeffs_)
          Info<< "Low-Dimension Coeffs:" << coeffDict_ << endl;}
       
        /**
        @brief execute the low-dimension process
         low-dimension process of scalarfield, which is the main process of LDT 
         the improved process is both feasible for uniform and non-uniform pattern 
        */
        void ldProcess();

        /**
        @brief calculate the low-dimension error
         low-dimension error calculation, also the main process of LDT
         the improved process is both feasible for uniform and non-uniform pattern 
        */
        void ldError();
        /**
         @brief totally virtual function of ABC
         show the process info and results 
         totally virtual here, will be reloaded in derived class 
        */ 
        virtual void resultShow()=0;

        /**
         @brief totally virtual function of ABC 
         decide whether to update the ld-settings, then update the l-d settings
         vital important stage for automatic method
         totally virtual here, will be reloaded in derived class 
        */ 
        virtual bool update()=0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace lowDimensionMethod
} // End namespace incompressible
} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif
// ************************************************************************* //
