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
    Foam::incompressible::lowDimensionMethod::nonUniformMethod

Description
    non-Uniform divide method for low-dimension process of Low-Dimension Tool (LDT).

SourceFiles
    nonUniformMethod.C

\*---------------------------------------------------------------------------*/

/**

@class nonUniformMethod

@brief Derived class from automaticLowDimension, for the construction of non-uniform cut points.
    
@author Jie Ren

@date 31 Oct 2019

@file nonUniformMethod.H

*/

/*---------------------------------------------------------------------------*/

#ifndef nonUniformMethod_H
#define nonUniformMethod_H

#include "automaticLowDimension.H"

namespace Foam
{
namespace incompressible
{
namespace lowDimensionMethod
{
/*---------------------------------------------------------------------------*\
                           Class nonUniformMethod Declaration
\*---------------------------------------------------------------------------*/
class nonUniformMethod
:
    public automaticLowDimension
{

protected:

 // Protected data
    
    /** Dict read from user's input */
    dictionary                nUnicoeffDict_;
    /** user's setting: the selection of non-uniform model */
    word                      priority_;
    /** depend on the input object in constructor, to classify the low-dimension model */
    bool                      update_;
    /** user's setting: the normalized minimum zone thickness */
    const vector              nDimMinZoneThick_; 
    /** user's setting: the sublayer thickness, normalized to the max mesh size */
    const vector              nDimSubLayerThick_;
    /** calculated result: the max mesh size, normalized to each domain length */
          vector              nDimMaxMeshSize_;
    /** calculated result: the thickness of the sublayers in each direction */
          vector              sublayer_Thick;
    /** calculated result: the minimum number of sublayers of one zone in each direction */
    std::vector<int>          least_LayerNum_;
    /** calculated result: the number of total sublayer in each direction */
    std::vector<int>          sublayer_Num_;
    /** calculated result: the total volume of the mesh, used for normalization */
         scalar               total_volume_;
    /** calculated result: current best reletive error under exist cut point settings */
         vector               best_error;
    /** calculated result: the volume of each sublayer in each direction */
    std::vector<std::vector<double>>    	      
                              sublayer_V_;
    /** calculated result: the value of the scalar field of each sublayer in each direction, a volume-averaged value  */
    std::vector<std::vector<double>> 
			      sublayer_Scalar_;

 // Protected member function

    /**
      @brief calculate the max and min mesh size
      function to calculate the max and min mesh as reference for settings, all mesh size will be normalized to the domain length
    */ 
    void genMeshSizeL();

    /**
      @brief calculate the sublayer, get the number, vloume, scalar field of them.
      function to divide the sublayer based on the user's layer thickness settings, use the mesh size value calculated from genMeshSizeL(), generate the sublayer information. 
    */
    void genSubLayerC();

    /** 
      @brief the checkSet() function to avoid invalid settings and warning.
      to check if the number of sublayer, the low-dimension quantities, the layer thickness are valid. 
    */
    void checkSet();

    /** 
      @brief a subfunction of accuracyPoint(), to compare the current input point set
      input the cut point set and the direction, judge if it is better than current best setting, to update and record the best cut point strategy.
    */
    inline void accura_isbetter(std::vector<int>, const int);

    /** 
      @brief one of the two main method to generate the cut point for non-uniform divide, find the most suitable one but may a cost long time.
      by using a 'smart' combination algorithm, producing all the possible point settings, and then by calling subfunction accura_isbetter() to find the best solution
    */
    bool accuracyPoint();

    /** 
      @brief the first of two subfunctions of fastPoint(), to find the cut section
      input the current cutpoint set--std::vector<int> and the direction of it--const int, then calculate the error, find the highest one section and return the number of it.
    */
    int fast_spaceFind(std::vector<int>, const int);

    /** 
      @brief the second of two subdunctions of fastPoint(), to find the cut point under appointed section
       int put the start and end of the section--int & int, together with the direction--const int. then find the best cutpoint between the start and end.
    */
    int fast_pointFind(int, int, const int);

    /** 
      @brief one of the two main method to generate the cutpoint for non-uniform divide, the fast but not accurate enough one.
      find the point by using bisection thoughts, first find the worst section--which means the largest reletive error, then in that section, find a point to cut the section into two part producing the least error.
    */
    bool fastPoint();

public:

 // Constructors

        /**
        @brief the normal constructor
        the normal constructor, construct from a given base class--automaticLowDimension, this base class can be a class promotion one, depending on the input object.
        */
        nonUniformMethod(const automaticLowDimension&);

 //- Destructor
  virtual ~nonUniformMethod()   {}

  // Member Functions

 	/**
        @brief decide whether to update the settings
        decide whether to update the low-dimension settings.Realized by optionally calling the former update() function, depending on the input object type when initialization
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
