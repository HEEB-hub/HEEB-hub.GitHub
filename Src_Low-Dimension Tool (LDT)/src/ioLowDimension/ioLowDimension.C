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
    Foam::incompressible::lowDimensionMethod::ioLowDimension

Description
    Namespace for low-dimensional tools(LDT), which is a post-process procedure.

Class
    Foam::incompressible::lowDimensionMethod::ioLowDimension

Description
    Creating low-dimension mesh and field process of Low-Dimension Tool (LDT)

SourceFiles
    ioLowDimension.C

\*---------------------------------------------------------------------------*/

#include "ioLowDimension.H"
#include "OSspecific.H" 
#include <stdlib.h>
#include <sstream>

using std::stringstream;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
namespace incompressible
{
namespace lowDimensionMethod
{
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
ioLowDimension::ioLowDimension
(  const lowDimensionMethod& LDM, word time_
)
:
  lowDimensionMethod(LDM),
  ldDataOutput_(lookupOrDefault<Switch>("ldDataOutput",false)),
  timeName(time_)
{ 
  meshSet=dictionary(word("meshSet"));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

  //- show the process info and results
void ioLowDimension::resultShow()
    { 
       if(ldDataOutput_)
       {
       Info<<"IO utility works successfully! \n\nLow-Dimension Mesh and Field now can be found in ./LD folder!\n\nInput 'cd LD' and 'paraFoam' command to check the lowdimension field in ParaView!\n\n"<<endl;
       }
    }

  //generate the blockMesh dict for mesh construction
void ioLowDimension::ldMeshDict()
    {
   
    meshSet.add("convertToMeters",1); //keyword 1 complete

    stringstream vert;
    vector* verts=new vector[8];
    vert << "(";

    for (int k=0;k<2;k++)
    {  for (int j=0;j<2;j++)
       {  for (int i=0;i<2;i++)
         {
          int PN=k*4+j*2+i;
          verts[PN]=vector(
                           fieldVector.value().x()*i+minVector.value().x(),
                           fieldVector.value().y()*j+minVector.value().y(), 
                           fieldVector.value().z()*k+minVector.value().z()
                                     );
          vert << " ( "<<verts[PN].x()<<" "
                      <<verts[PN].y()<<" "
                      <<verts[PN].z()
                      <<" )";
         }
       }
    }

    vert << " )";
    string vertices = vert.str();
    delete [] verts;

    meshSet.add("vertices",word(vertices),1);  //keyword 2 complete

    word nam[3]={"Xpos","Ypos","Zpos"};
    for (int i=0;i<3;i++)
    {
      stringstream pointset;
      pointset << " (";
      for (int j=0;j<ld_Dimension_.value()[i];j++)
      pointset << " ( " << ld_cutpoint_[i][j+1]-ld_cutpoint_[i][j] << " 1 1 ) ";
      pointset << " )";
      string ptset = pointset.str();
      meshSet.add(nam[i],word(ptset),1);       //keyword 3 complete
    }					       

    stringstream blockset;
    blockset << "( hex ( 0 1 3 2 4 5 7 6 ) ";
    blockset << "( "<<ld_Dimension_.value().x()<<" "
                    <<ld_Dimension_.value().y()<<" "
                    <<ld_Dimension_.value().z()
                    <<" )";
    blockset << " simpleGrading ( $Xpos $Ypos $Zpos ) )";
    string bkset = blockset.str();
    meshSet.add("blocks",word(bkset),1);       //keyword 4 complete


    stringstream bdset;

    bdset << "( )";

    string bddset = bdset.str();

    meshSet.add("boundary",word(bddset),1);    //keyword 5 complete

    Info<<"ldMesh dictionary cerate successfully:"<<endl<<meshSet<<endl; }


  //generate the blockMesh for ld field, from IOdictionary
void ioLowDimension::ldMeshCreate()
    {
  
    mkDir("LD/system"); // create folders

    IOdictionary mymeshDict //create IOdictionary
    (
        IOobject
        (
            "blockMeshDict",
            "constant",
            C_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        ),
            meshSet);

    mymeshDict.Foam::regIOobject::write();

    int oa=system(" mv  ./constant/blockMeshDict  ./LD/system"); 
    int ob=system(" cp  ./system/controlDict      ./LD/system"); 
    int oc=system(" cp  ./system/fvSchemes        ./LD/system"); 
    int od=system(" cp  ./system/fvSolution       ./LD/system");     
    int oe=system(" cd LD ; blockMesh"); 

    if(oa*4==(ob+oc+od+oe))
    Info<<"ldMesh create successfully!\nCreating low-dimension field..."<<endl; }



  //generate the low dimension field for visualization
void ioLowDimension::ldFieldCreate()
    {

  int oA=system(" mv  ./constant/polyMesh     ./system"); 
  int oB=system(" cp -r ./LD/constant/polyMesh  ./constant"); 

  mkDir(word("LD/"+timeName));

    fvMesh ldmesh
        (
        IOobject
            (
                fvMesh::defaultRegion,
                "constant",
                C_.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

    word fieldName(lookupOrDefault<word>("lowDimensionField","C"));

    word fieldName1=fieldName+"_ld";

    volScalarField c_ld
    (
        IOobject
        (
                fieldName1,
                "constant",
                ldmesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
        ),
         ldmesh,
         dimensionedScalar(fieldName1, C_.dimensions() , 0)
    );        

   word fieldName2=fieldName+"_ld_error";

   volScalarField c_ld_error
    (
        IOobject
        (
                fieldName2,
                "constant",
                ldmesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
        ),
         ldmesh,
         dimensionedScalar(fieldName2, dimless , 0)
    );    

    forAll (c_ld,celli)
         {
          const vector abs_coord=ldmesh.C()[celli]-minVector.value();
                vector temp_pos(0,0,0);
        for(int i=0;i<3;i++)
          {
           for(int j=0;j<ld_Dimension_.value()[i];j++)
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
          c_ld[celli]=ld_Scalar_[ZoneLabel].value();
          c_ld_error[celli]=ld_SubError_[ZoneLabel].value();
         };

    c_ld.Foam::regIOobject::write();
    c_ld_error.Foam::regIOobject::write();
    word com1(" mv ./constant/"+fieldName1+"  ./LD/"+timeName);
    word com2(" mv ./constant/"+fieldName2+"  ./LD/"+timeName);

    int oC=system(com1);
    int oD=system(com2);
    int oE=system(" rm -r ./constant/polyMesh"); 
    int oF=system(" mv ./system/polyMesh     ./constant"); 

    if(oA*5==(oB+oC+oD+oE+oF))
    Info<<"ldField create successfully!\n-------------------------------------------------\n"<<endl; }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace namespace lowDimensionMethod
} // End namespace incompressible
} // End namespace Foam

