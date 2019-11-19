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
    Foam::incompressible::lowDimensionMethod
 
Description
    Namespace for low-dimensional tools(LDT), which is a post-process procedure.

Class
    Foam::incompressible::lowDimensionMethod

Description
    Realization for low-dimensional tools(LDT).

SourceFiles
    lowDimensionMethod.C

\*---------------------------------------------------------------------------*/

#include "lowDimensionMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace lowDimensionMethod
{
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

 //- construct from assigned field and mesh
lowDimensionMethod::lowDimensionMethod
(  const volScalarField& C,
   const word& type
):

    IOdictionary
    (
        IOobject
        (
            "lowDimensionProperties",
            C.time().constant(),
            C.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mesh_(C.mesh()),
    C_(C),
    coeffDict_(subOrEmptyDict(type + "Coeffs")),
    printCoeffs_(lookupOrDefault<Switch>("printCoeffs", false)),
    minVector("minVector",mesh_.C().dimensions(),vector::zero),
    fieldVector("fieldVector",mesh_.C().dimensions(),vector::zero),
    ld_Dimension_("ld_Dimension_",dimless,vector::one),
    ld_cutpoint_(3),
    ld_fieldNums_(1)
{   
    ld_Scalar_  = new dimensionedScalar[ld_fieldNums_];
    ld_SubError_= new dimensionedScalar[ld_fieldNums_];
    ld_Error_=0;
    fieldRange();
    Info<<"Base class construction finished!"<<endl<<endl;
}

 //- Allow default bitwise copy construct
lowDimensionMethod::lowDimensionMethod (const lowDimensionMethod& LDM)
:
 IOdictionary (LDM),

        mesh_ (LDM.mesh_),

           C_ (LDM.C_),

   coeffDict_ (LDM.coeffDict_),

 printCoeffs_ (LDM.printCoeffs_),

    minVector (LDM.minVector),
    
  fieldVector (LDM.fieldVector),

ld_Dimension_ (LDM.ld_Dimension_),

    ld_Error_ (LDM.ld_Error_),

ld_cutpoint_ (LDM.ld_cutpoint_),

ld_fieldNums_ (LDM.ld_fieldNums_)

{
ld_Scalar_ =  new dimensionedScalar [ld_fieldNums_]; 
ld_SubError_= new dimensionedScalar [ld_fieldNums_];
for (int i=0; i< ld_fieldNums_; i++)
    {  
   ld_Scalar_[i] = LDM.ld_Scalar_[i];
   ld_SubError_[i] = LDM.ld_SubError_[i];
    }
 Info<<"Base class copy constructor runs successfully!"<< nl;

}
// * * * * * * * * * * * * * * *  Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Protected Member Functions  * * * * * *  * * * * //


/**
     calculate the fieldRange vector(min & whole) of the mesh domain	      */

void lowDimensionMethod::fieldRange() 
{
   const pointField & pointlist = mesh_.points();

     scalar xMax(max(pointlist & vector(1,0,0)));
     scalar yMax(max(pointlist & vector(0,1,0)));
     scalar zMax(max(pointlist & vector(0,0,1)));

     dimensionedVector    maxvector("maxvector",mesh_.C().dimensions(),vector(xMax,yMax,zMax));
      
     scalar xMin(min(pointlist & vector(1,0,0)));
     scalar yMin(min(pointlist & vector(0,1,0)));
     scalar zMin(min(pointlist & vector(0,0,1)));

     dimensionedVector minvector("minvector",mesh_.C().dimensions(),vector(xMin,yMin,zMin));
     dimensionedVector fieldvector("fieldvector",mesh_.C().dimensions(),maxvector.value()-minvector.value());
  
     minVector   = minvector;
     fieldVector = fieldvector;
}

/**
      calculate the cut point for the low-dimension conversion, average points   */


void lowDimensionMethod::setPoints()
{
   const vector subpoint_num(
       ld_Dimension_.value().x()+1,
       ld_Dimension_.value().y()+1,
       ld_Dimension_.value().z()+1);

   const vector subzone_len(
       fieldVector.value().x()/ld_Dimension_.value().x(),
       fieldVector.value().y()/ld_Dimension_.value().y(),
       fieldVector.value().z()/ld_Dimension_.value().z());

     for(int i=0;i<3;i++)
     {
       ld_cutpoint_[i].resize(subpoint_num[i]);
       for(unsigned int j=0;j<ld_cutpoint_[i].size();j++)
	{
            ld_cutpoint_[i][j]=j*subzone_len[i];
        }
      }
}


// * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * //


 //- low-dimension process of scalarfield
void lowDimensionMethod::ldProcess()
     {
       Info<<"Starting the Low-Dimension Dividing ...please wait"<<endl<<endl;
      
       dimensionedScalar* totalv =new dimensionedScalar[ld_fieldNums_];
       dimensionedScalar* totalvc=new dimensionedScalar[ld_fieldNums_];

       forAll (C_,celli)
       {
        const vector abs_coord=mesh_.C()[celli]-minVector.value();
              vector temp_pos(0,0,0);
        for(int i=0;i<3;i++)
          {
           for(unsigned int j=0;j<ld_Dimension_.value()[i];j++)
	     {
            if(abs_coord[i]  >  ld_cutpoint_[i][j] &&
	       abs_coord[i]  <= ld_cutpoint_[i][j+1] )
               { 
		temp_pos[i] = j;
		break;
	       }
             }
           }
          label ZoneLabel=temp_pos[0]+temp_pos[1]*ld_Dimension_.value().x()+temp_pos[2]*ld_Dimension_.value().x()*ld_Dimension_.value().y();
          totalv[ZoneLabel]  += mesh_.V()[celli];
          totalvc[ZoneLabel] += mesh_.V()[celli]*C_[celli];
       }

        for (int i=0;i<ld_fieldNums_;i++)
        {  
          ld_Scalar_[i]=totalvc[i]/totalv[i];
        }

       Info<<"Dividing into "<<ld_fieldNums_<<" zones successfully, with " 
            <<ld_Dimension_.value()<<" pieces in each (X,Y,Z) coordinate!"<<endl<<endl;

        delete [] totalv;
        delete [] totalvc;
     }


 //- low-dimension error calculation
void lowDimensionMethod::ldError()
     {

       Info<<"Calculating the Low-Dimension Error ...please wait"<<endl<<endl;
       dimensionedScalar totalError2(0);
       dimensionedScalar totalZone2(0);

       dimensionedScalar* totalv    = new dimensionedScalar[ld_fieldNums_];
       dimensionedScalar* subError2 = new dimensionedScalar[ld_fieldNums_];
       dimensionedScalar* subZone2  = new dimensionedScalar[ld_fieldNums_];

       dimensionedScalar errorVar(0);
       dimensionedScalar totalV(0);

        forAll (C_,celli)
         {
          const vector abs_coord=mesh_.C()[celli]-minVector.value();
          vector temp_pos(0,0,0);
          for(int i=0;i<3;i++)
           {
             for(unsigned int j=0;j<(ld_cutpoint_[i].size()-1);j++)
	      {
              if(abs_coord[i]  >  ld_cutpoint_[i][j] &&
	         abs_coord[i]  <= ld_cutpoint_[i][j+1] )
                { 
	       	temp_pos[i] = j;
	  	break;
	        }
              }
           }
          label ZoneLabel=temp_pos[0]+temp_pos[1]*ld_Dimension_.value().x()+temp_pos[2]*ld_Dimension_.value().x()*ld_Dimension_.value().y();
      
          totalError2 += mesh_.V()[celli]*pow(ld_Scalar_[ZoneLabel]-C_[celli],2);
          totalZone2  += mesh_.V()[celli]*pow(C_[celli],2);

          totalv[ZoneLabel]    += mesh_.V()[celli];
          subError2[ZoneLabel] += mesh_.V()[celli]*pow(ld_Scalar_[ZoneLabel]-C_[celli],2);
          subZone2[ZoneLabel]  += mesh_.V()[celli]*pow(C_[celli],2);

      //    totalError2 += pow(ld_Scalar_[ZoneLabel]-C_[celli],2);
      //    totalZone2  += pow(C_[celli],2);
        }

         ld_Error_=sqrt(totalError2/totalZone2);

        for (int i=0;i<ld_fieldNums_;i++)
        {  
          ld_SubError_[i]=sqrt(subError2[i]/subZone2[i]);
          errorVar += totalv[i]*pow(ld_SubError_[i]-ld_Error_,2);
          totalV += totalv[i];
        }

          errorVar=sqrt(errorVar/totalV); 
        
      Info<<"The Low-Dimension Discretization Error of "<<ld_Dimension_.value()<<" set equals "<<ld_Error_.value()*100<<"%"<<endl<<endl;
      Info<<"The Low-Dimension Discretization Error variance of "<<ld_Dimension_.value()<<" set equals "<<errorVar.value()*100<<"%"<<endl<<endl;

        delete [] totalv;
        delete [] subError2;
        delete [] subZone2;
     }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace namespace lowDimensionMethod
} // End namespace incompressible
} // End namespace Foam

