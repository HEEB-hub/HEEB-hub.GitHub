Description
    Low-Dimension Tool - a self-adaptive low-dimension package/utility for dimension reduction 
	                 of high-resolution CFD field that contains the distribution of a given 
			 scalar (e.g. concentration, temperature). 						 
    P.S. The open-source CFD platform OpenFOAM is used for the package development of this 
	self-adaptive low-dimension tool (LDT).

Author
    Jie Ren
    Research Assistant, Guangzhou University, China
	
Date: Nov. 2019

License
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


Installation instructions (for Linux):
--------------------------

1. Ensure that OpenFOAM and doxygen are installed and properly configured.

   Doxygen can be not necessary if one already knows the usage of LDT.

2. Open a terminal, navigate to folder as following, then type: wmake. 

   We strongly suggest you putting this file in the OpenFOAM's user 
   directory called "applications", and check if the folder name is 
   "lowDimensionTool". 
   If error happens when compiling, pls check the settings in "options".

3. If all goes well, an application called "lowDimensionTool" is now created.

4. You can optionally create the code documentation by running:
   doxygen doxygen_log.txt


Usage and notice:
------
Navigate to folder of your CFD case for dimension reduction.
Then type the command: lowDimensionTool

Examples are provided to illustrate some processing results of the lowDimensionTool 
in our corresponding research paper, which will be available in the author's RG page:
https://www.researchgate.net/profile/Jie_Ren91

Notice:

 1. A scalar field must be given in advance. (concentration, temperature etc.).
 
 2. The given file "lowDimensionProperties" is necessary as the pre-settings of the tool.
    It should be located in the "constant" folder in your case directory of OpenFOAM.
	
 3. More details of dimension reduction settings in "lowDimensionProperties" file. 
    Please check these settings before running the application.
