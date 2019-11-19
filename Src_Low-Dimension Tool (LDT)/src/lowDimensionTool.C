/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    
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

Application
    lowDimensionTool

Description
    Create lowDimension field for scalar, including the construction of low-dimension mesh, low-dimension field and the discretization error.

\*---------------------------------------------------------------------------*/

/**
@mainpage
This is the <B>Low-Dimension Tool (LDT)</B>, which is a self-developed dimension reduction utility based on the open-source CFD platform OpenFOAM.

LDT performs a dimension reduction procedure based on the proposed <B>Self-adaptive Low-dimension Ventilation Model</B>.

The proposed <B>Self-adaptive Low-dimension Model</B>, including <B>direct</B> and <B>automatic</B> modes, together with an advanced <B>non-uniform</B> dividing module, is able to execute a series of low-dimension conversions of given high-resolution scalar field.

At the end of the processing, a new low-dimension <B>mesh</B> and corresponding<B> scalar field</B> will be constructed in the <B> ./LD folder </B>, which can be easily accessed by OpenFOAM  (i.e.  ParaView) for visualization. 

<B>Methods implemented</B>
   - lowDimensionMethod
        - import (CFD scalar field & mesh) 
        - execution of low-dimension process
        - low-dimension error calculation
        - advanced non-uniform dividing method
        - low-dimension magnitude update
        - export (Low-Dimension mesh & low-dimension scalar field)


@author Jie Ren,  Shijie Cao

@date November  2019

@file lowDimensionTool.C

*/

#include "argList.H"
#include "nonUniformMethod.H"
#include "ioLowDimension.H"

using namespace Foam;
using namespace Foam::incompressible::lowDimensionMethod;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
   #include "setRootCase.H"
   #include "createTime.H"
   #include "createMesh.H"
   #include "readPreSettings.H"

   // Check the settings and create the corresponding class for further step.  
   
   Info<<"Starting Low-Dimension process"<<endl<<endl;

   word judge=lowDimensionProperties.lookup("lowDimensionMethod");      

   directLowDimension* LDT;

   if("direct"==judge)         
   LDT = new directLowDimension(C);

   else if ("automatic"==judge) 
   LDT = new automaticLowDimension(C);

   else
    {
     Info<<"Invalid Settings of lowDimensionMethod, program terminated!"<<endl;      
     return 0;
    }
   
   LDT->printCoeffs();

   Switch nonUniform = lowDimensionProperties.lookup("nonUniformDivide"); 
   
   //Check and decide whether to execute the non-Uniform process, take advantage of RTTI feature of C++, realizing the dynamic type check & construction.
 
   if (nonUniform)
   {
   nonUniformMethod* temp;
     if (auto LDA = dynamic_cast <automaticLowDimension *> (LDT))
       temp = new nonUniformMethod(*LDA);
     else
       temp = new nonUniformMethod(*LDT);
   delete LDT;
   LDT = temp;
   }
   
   //Starting the low-dimension conversion and update loop. The class has been designed by virtual function feature, which simplified the code and process!
   
   do
   {
     LDT->ldProcess();
     LDT->ldError();
   }
   while (LDT->update());
   
   //Check and decide whether to output the low-dimension mesh and scalar field.
    
   ioLowDimension ioLDT(*LDT,runTime.timeName());
   if(ioLDT.dataOutput())
   {
     ioLDT.ldMeshDict();
     ioLDT.ldMeshCreate();
     ioLDT.ldFieldCreate();
   }
     LDT->resultShow();

     ioLDT.resultShow();

   delete LDT;

   Info<<"\nEnd."<<nl;
   return 0;
}


// ************************************************************************* //
