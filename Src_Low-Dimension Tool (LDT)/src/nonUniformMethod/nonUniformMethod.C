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
    Foam::incompressible::nonUniformMethod

Description
    non-Uniform divide method for low-dimension process of Low-Dimension Tool (LDT)

SourceFiles
    nonUniformMethod.C

\*---------------------------------------------------------------------------*/

#include "nonUniformMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace lowDimensionMethod
{
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


nonUniformMethod::nonUniformMethod(const automaticLowDimension& LDM)
:
    automaticLowDimension(LDM),
  
    nUnicoeffDict_(subOrEmptyDict("nonUniformCoeffs")),

    priority_ (nUnicoeffDict_.lookupOrDefault<word>
		  (  "dividePriority", "accuracy" )),
    update_(true),
  
    nDimMinZoneThick_ (nUnicoeffDict_.lookupOrDefault<vector>
		          (   "minZoneThick", vector::zero )),

    nDimSubLayerThick_ (nUnicoeffDict_.lookupOrDefault<vector>
		          (   "subLayerThick", vector::one )),

    nDimMaxMeshSize_(vector::zero),

    sublayer_Thick(vector::zero),

    least_LayerNum_(3),

    sublayer_Num_(3),

    total_volume_(0),

    best_error(vector::zero),
    
    sublayer_V_(3),
    
    sublayer_Scalar_(3)
{
    if (Step==0)  update_=false;
    genMeshSizeL();
    genSubLayerC();
    if (accuracyPoint() || fastPoint() ) 
    Info<<"Non-Uniform finished under "<< priority_ <<"method!"<<nl<<nl;
    Info<<"Derived class: nonUniformMethod construction finished!"<<nl<<nl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void nonUniformMethod::genMeshSizeL()
{
  const faceList & ff = mesh_.faces();//mesh faces presented as 4 points;
  const pointField & pp = mesh_.points(); //mesh points' specific (x,y,z);
  std::vector < std::vector<double> > mlen(3);
  forAll (mesh_.C(), celli)
   {
    const cell & cc = mesh_.cells()[celli];//single mesh presented as 6 faces
    labelList pLabels(cc.labels(ff));//single mesh presented as points' label, important
    pointField pLocal(pLabels.size(), vector::zero);

     forAll (pLabels, pointi) //find the (x,y,z) of every point of mesh
      {
          pLocal[pointi] = pp[pLabels[pointi]];
      }
     scalar xDim=Foam::max(pLocal & vector(1,0,0)) - Foam::min(pLocal & vector(1,0,0));
     scalar yDim=Foam::max(pLocal & vector(0,1,0)) - Foam::min(pLocal & vector(0,1,0));
     scalar zDim=Foam::max(pLocal & vector(0,0,1)) - Foam::min(pLocal & vector(0,0,1));
     vector Dim(xDim,yDim,zDim);
     for(int i=0;i<3;i++)
       mlen[i].push_back(Dim[i]);
   }
   vector nDimMinMeshSize_(vector::zero);
  for(int i=0;i<3;i++)
   {   
         auto maxPosition = max_element(mlen[i].begin(), mlen[i].end());
	 auto minPosition = min_element(mlen[i].begin(), mlen[i].end());
           nDimMaxMeshSize_[i] = *maxPosition/fieldVector.value()[i];
 	   nDimMinMeshSize_[i] = *minPosition/fieldVector.value()[i];
  Info<<"nonDimension Max Mesh Size (percent of domain length) in " 
      << char('x'+i)<<" direction: "<< nDimMaxMeshSize_[i] << nl;
  Info<<"nonDimension Min Mesh Size (percent of domain length) in " 
      << char('x'+i)<<" direction: "<< nDimMinMeshSize_[i] << nl <<nl;
   }
  Info<<"nonDimension Mesh Size above for setting reference!\n" 
      <<"Press 'Enter' to continue!"<<nl;
  cin.get();
}


void nonUniformMethod::genSubLayerC()
{
  for(int i=0;i<3;i++)
  {
    double nDimLay_Thick_= nDimMaxMeshSize_[i] * nDimSubLayerThick_[i];
    sublayer_Num_[i] = int (pow(nDimLay_Thick_ , -1));
    double minLayNum_temp = nDimMinZoneThick_[i] / (1/double(sublayer_Num_[i]));
    least_LayerNum_[i] = int (minLayNum_temp);
    if (minLayNum_temp - least_LayerNum_[i] > 0 || least_LayerNum_[i]==0) 
    least_LayerNum_[i]++;
    sublayer_Thick[i] = fieldVector.value()[i] / sublayer_Num_[i];
    dimensionedScalar* laytotalv  = new dimensionedScalar[sublayer_Num_[i]];
    dimensionedScalar* laytotalvc = new dimensionedScalar[sublayer_Num_[i]];
    forAll (C_,celli)
    {
      const double abs_pos = mesh_.C()[celli][i]-minVector.value()[i];
      int temp_pos = int (abs_pos / sublayer_Thick[i]) ;
      laytotalv  [temp_pos] += mesh_.V()[celli];
      laytotalvc [temp_pos] += mesh_.V()[celli] * C_[celli];
    }
    for (int j=0;j<sublayer_Num_[i];j++)
    {
       if (!i) total_volume_ += laytotalv[j].value();
       if (laytotalv[j].value() != 0 )
       {
         sublayer_V_[i].push_back ( laytotalv[j].value() );
         sublayer_Scalar_[i].push_back ( (laytotalvc[j] / laytotalv[j]).value() );
       }
       else
       {
         Info<<"Too small 'subLayerThick' set for present mesh!"<<nl
             <<"Please check and enlarge!"<<nl;
         std::abort();
       }
    }
    Info<<"Dividing into "<<sublayer_Num_[i]<<" sub-layers in "
        << char('x'+i) <<" direction  successfully! "<<nl
        <<"The minimun sublayer number is "<< least_LayerNum_[i] <<nl<<nl;
    delete [] laytotalv;
    delete [] laytotalvc;
  }
  checkSet();
}


void nonUniformMethod::checkSet()
{
for (int i=0; i<3; i++)
  {
    if(
      nDimMinZoneThick_[i] * ld_Dimension_.value()[i] > 1 || 
      sublayer_Num_[i] < ld_Dimension_.value()[i] ||
      least_LayerNum_[i] * ld_Dimension_.value()[i] > sublayer_Num_[i] 
      )
    {
      Info<< "sublayer_Num_ "<< char('x'+i)<<" :"<< sublayer_Num_[i]<<nl;
      Info<< "nDimMinZoneThick_ "<< char('x'+i)<<" :"<< nDimMinZoneThick_[i]<<nl;
      Info<< "least_LayerNum_ "<< char('x'+i)<<" :"<< least_LayerNum_[i]<<nl;
      Info<<"LDT program terminated in "<<ld_Dimension_.value()
          <<" low-dimension divide pattern!"<<nl
          <<"Reasons may be below:"<<nl
          <<"1. Too big 'minZoneThick' set for L-D dimension above."<<nl
          <<"2. Too big 'subLayerThick' set for L-D dimension above."<<nl
          <<"3. Too strict error limit set."<<nl
          <<"Luck to you, bye~"<<nl;
      std::abort();
    }
  }
}

inline void nonUniformMethod::accura_isbetter( std::vector<int> combine_set ,const int i)
{
  int k = least_LayerNum_[i];
  const int m = ld_Dimension_.value()[i]-1;
  scalar temp_error=0;
  scalar temp_v[m+1]{};
  scalar temp_c[m+1]{};
  scalar temp_e[m+1]{};
  combine_set.push_back(sublayer_Num_[i]);
  for (int kk=1; kk<m+1; kk++)
  combine_set[kk] += kk*(k-1);
  for (int lay=0; lay < sublayer_Num_[i]; lay++)
  {
    for (int p=0; p < m+1; p++)
    {
      if (lay >= combine_set[p] && lay < combine_set[p+1])
      {     
      temp_v[p] += sublayer_V_[i][lay];
      temp_c[p] += sublayer_V_[i][lay]*sublayer_Scalar_[i][lay];
      break;
      }
    }
  }
  for (int p=0; p < m+1; p++)
     temp_c[p] = temp_c[p] / temp_v[p];
  for (int lay=0; lay < sublayer_Num_[i]; lay++)
  {
    for (int p=0; p < m+1; p++)
    {
      if (lay >= combine_set[p] && lay < combine_set[p+1])
      {  
        temp_e[p] += sublayer_V_[i][lay]*pow(sublayer_Scalar_[i][lay]-temp_c[p],2);
        break;
      }
    }
  }
  for (int p=0; p < m+1; p++)
  { 
    temp_e[p] = Foam::sqrt(temp_e[p]/temp_v[p]);
    temp_error += temp_e[p]*temp_v[p]/total_volume_;
  }
  if (best_error[i]==0) best_error[i] = temp_error; 
  if (temp_error < best_error[i])
  {
    best_error[i]=temp_error;
    for(unsigned int j=0;j<ld_cutpoint_[i].size();j++)
      ld_cutpoint_[i][j] = combine_set[j]*sublayer_Thick[i];
  }
}


bool nonUniformMethod::accuracyPoint() //realization of the 'precise sub-layer recombination' algorithm to find the accuracy cut point
{
checkSet();
best_error=vector::zero;
  if( priority_ == "accuracy" )
  {
    for(int i=0;i<3;i++)
    {
      const int k = least_LayerNum_[i];
      const int m = ld_Dimension_.value()[i]-1;
      const int n = sublayer_Num_[i]-ld_Dimension_.value()[i]*(k-1)-1;
      std::vector<int>  temp_comb(m+1); 
      temp_comb[0]=0;
      int level(2); 
      const int maxlevel(m); 
      const int volumn(n-m);
      if (m==0) continue;
      if (m==1) 
      {
        for(int j=1; j < n+1; j++) 
        {
          temp_comb[1]=j;
          accura_isbetter(temp_comb,i);
        }
      }
      else
      {
        temp_comb[1]=1;
        while(1)
        {
          temp_comb[level] = temp_comb[level-1] + 1;
          if(level==maxlevel) 
          {
            accura_isbetter(temp_comb,i);
            for(int j=temp_comb[level]+1; j<n+1; j++)
            {
              temp_comb[level] = j;
              accura_isbetter(temp_comb,i);
            }
            level--;
            while(temp_comb[level]==level+volumn && level>=1) 
            level--;
            if(level<1)  break;
            temp_comb[level]++;
          }
          level++;
        }
      }
      for(unsigned int j=0;j<ld_cutpoint_[i].size();j++)
      {
        Info <<"Optimal cutpoint in "<< char('x'+i) << " direction: " 
             << ld_cutpoint_[i][j] << nl;
      }
      Info <<nl<<nl;
    }
  return true;
  }
  else
  return false;
}


int nonUniformMethod::fast_spaceFind(std::vector<int> input_cut, const int i)
{
  int  m = input_cut.size()-1;
  int  space_find=0;
  scalar temp_v[m]{};
  scalar temp_c[m]{};
  scalar temp_e[m]{};
  for (int lay=0; lay < sublayer_Num_[i]; lay++)
  {
    for (int p=0; p < m; p++)
    {
      if (lay >= input_cut[p] && lay < input_cut[p+1])
      {     
        temp_v[p] += sublayer_V_[i][lay];
        temp_c[p] += sublayer_V_[i][lay]*sublayer_Scalar_[i][lay];
        //Hx=();
        break;
      }
    }
  }
  for (int p=0; p < m; p++)
    temp_c[p] = temp_c[p] / temp_v[p];
  for (int lay=0; lay < sublayer_Num_[i]; lay++)
  {
    for (int p=0; p < m; p++)
    {
      if (lay >= input_cut[p] && lay < input_cut[p+1])
      {  
        temp_e[p] += sublayer_V_[i][lay]*pow(sublayer_Scalar_[i][lay]-temp_c[p],2);
        break;
      }
    }
  }
  for (int p=0; p < m; p++)
  {
    temp_e[p] = Foam::sqrt(temp_e[p]/temp_v[p]);
  if (p>0)
  if (temp_e[p]>temp_e[space_find] && 
      input_cut[p+1]-input_cut[p]>=2*least_LayerNum_[i])
    space_find=p;
  }
  return space_find;
}


int nonUniformMethod::fast_pointFind(int start, int end, const int i)
{
  int  point_find=start+least_LayerNum_[i];
  double min_error=0;
  for (int j=start+least_LayerNum_[i];j<=end-least_LayerNum_[i];j++)
  {
    scalar temp_error=0;
    scalar temp_v[2]{};
    scalar temp_c[2]{};
    scalar temp_e[2]{};
    for (int lay=start; lay < end; lay++)
    {
      if (lay<j)
      {
        temp_v[0] += sublayer_V_[i][lay];
        temp_c[0] += sublayer_V_[i][lay]*sublayer_Scalar_[i][lay];
      }
      else 
      {
        temp_v[1] += sublayer_V_[i][lay];
        temp_c[1] += sublayer_V_[i][lay]*sublayer_Scalar_[i][lay];
      }
    }

    for (int p=0; p < 2; p++)
      temp_c[p] = temp_c[p] / temp_v[p];

    for (int lay=start; lay < end; lay++)
    {
      if (lay<j)          
      temp_e[0] += sublayer_V_[i][lay]*pow(sublayer_Scalar_[i][lay]-temp_c[0],2);
      else
      temp_e[1] += sublayer_V_[i][lay]*pow(sublayer_Scalar_[i][lay]-temp_c[1],2);
    }

    for (int p=0; p < 2; p++)
    { 
      temp_e[p] = Foam::sqrt(temp_e[p]/temp_v[p]);
      temp_error += temp_e[p]*temp_v[p]/(temp_v[0]+temp_v[1]);
    }
    if (min_error==0)
      min_error = temp_error;
    if (temp_error < min_error)
    {
      min_error = temp_error;
      point_find=j;
    }
  }
  return point_find;
}


bool nonUniformMethod::fastPoint()  //a dividing method based on the information entropy evaluation.
{
checkSet();
  if( priority_ == "fast" )
  {
    for (int i=0;i<3;i++)
    {
      std::vector<int> temp_point{0,sublayer_Num_[i]};
      while (temp_point.size()<ld_Dimension_.value()[i]+1)
      {
      int cutSpace=fast_spaceFind(temp_point, i);
      int cutPoint=fast_pointFind(temp_point[cutSpace], temp_point[cutSpace+1], i);
      temp_point.push_back(cutPoint);
      sort(temp_point.begin(), temp_point.end());
      }
      for(unsigned int j=0;j<ld_cutpoint_[i].size();j++)
      {
      Info << "temp_point[" << j <<"]: "<< temp_point[j]<<nl;
      ld_cutpoint_[i][j] = temp_point[j]*sublayer_Thick[i];
      }
    }
  return true;
  }
  else
  return false;
}

// * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * //

//- decide whether to update the ld-settings, then do it and return a flag

bool nonUniformMethod::update()   //update the field number and cut point
{
if (update_)
  {    
    if (automaticLowDimension::update())
    return (accuracyPoint() || fastPoint());
    else
    return false;
  }
else
return directLowDimension::update();
}

 //- show the process info and results
void nonUniformMethod::resultShow()
{
  Info<<"Low-Dimension Tool (LDT) works successfully!"
      <<endl<<endl
      <<"Results as below:"
      <<endl<<endl
      <<"The computational domain has been divided into " 
      <<ld_fieldNums_
      <<" different size zones, with " 
      <<ld_Dimension_.value()
      <<"pieces in each (X,Y,Z) coordinate!"
      <<endl<<endl
      <<"The total L-D ERROR equals "
      <<ld_Error_.value()*100<<"%"
      <<endl<<endl;
  if (update_)   
  {
    Info<<"Total iteration: "<<Step
        <<" steps."<<endl<<endl
        <<"Congrats! Non-Uniform Automatic Low-Dimension finallized!"<<endl<<endl;
  }
  else 
  Info<<"Congrats! Non-Uniform direct Low-Dimension finallized!"<<endl<<endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace lowDimensionMethod
} // End namespace incompressible
} // End namespace Foam

