*****************************************************************************
LocoMouse tracking algorithm from “A quantitative framework for whole-body 
coordination reveals specific deficits in freely walking ataxic mice"
Machado, Darmohray, Fayad, Marques, and Carey.
eLife (2015) http://dx.doi.org/10.7554/eLife.07892.

Code written by Joao Fayad (joaofayad@gmail.com).

------------------ PLEASE NOTE -----------------------------------------------
This version of LocoMouse retains the analysis procedures as used in the 
cited paper. An updated version (Matlab R2014B and higher) can be found under 
http://github.com/careylab/LocoMouse_Dev (careylab's development version)
------------------ Dennis Eckmeier, 2015 ------------------------------------

If you find this code useful, please reference our paper.
*****************************************************************************

--|Intro|-- 
LocoMouse_Tracker is a software developed in MATLAB R2013 for tracking
features of locomoting mice when observing them from the side and
bottom view simultaneously. It uses Support Vector Machines to learn
image filters for each of the features and a multi-target tracking
framework to resolve the most likely tracks from the image observations.
This software was developed for a specific setup, and so replicating such
conditions is essential for it to work as intended. The general framework
is, however, flexible enough to be modified for a different setup.

Please read the following files before using this code:

- "README" (this file) - All other license and readme files contained in
the "additional_packages" folder.

--|License & Disclaimer|--
Copyright (C) 2015  Joao Fayad (joaofayad@gmail.com)
 
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


--|Requirements|-- 

Matlab toolboxes: 
- Image processing toolbox 
- Parallel toolbox (not strictly required, but very slow otherwise). 
- Signal processing toolbox - Statistics toolbox

Additional packages (already included): 
- A modified version of the following algorithm:

Efficient Second Order Multi-Target Tracking with Exclusion Constraints
Russel C., Setti F., Agapito L. British Machine Vision Conference (BMVC
2011), Dundee, UK, 2011.

- sc (https://github.com/ojwoodford/sc)

- combinator
(http://www.mathworks.com/matlabcentral/fileexchange/24325-combinator-combinations-and-permutations)

These packages are already included in this release. Please read the
respective copyright and license files before using this software.

--|Usage|-- 
To use this software please edit the MTF_main.m file so the different
paths point to the appropriate folders. An example under a typical
installation of MATLAB is provided. Make sure the path to all the
functions are added to your MATLAB path, including the external packages
folder.

To run the example you will need a video example. Please go to
https://www.dropbox.com/sh/hfu0sayfzwoqwp1/AADcMwYPl0UHw8jTDIRIHBtoa?dl=0
 and download 'G6AE1_98_28_1_control_S1T1.avi' and
 'G6AE1_98_28_1_control_S1T.png' into the
 LocoMouse/movies/3_11_2013_S1/G6AE1 folder.

Optionally, you can also download the contents of
'G6AE1_98_28_1_control_S1T1_corrected' into LocoMouse/output/Distortion
images/G6AE1_98_28_1_control_S1T1_corrected. This will save computational
time.

The list of specific files to track function follows the specific naming
convention of our system and therefore should be used carefully.

Binaries for 64bit Linux and Windows are included. If using a different
system, please compile the 'combinator' and 'tracking' packages under
'LocoMouse/external packages'.

In recent versions of MATLAB the behaviour of VideoReader has changed.
Therefore, if problems arise, move (and rename accordingly) the functions
inside 'LocoMouse/tracking code/matlab2014' to the 'LocoMouse/tracking
code' folder.

--|Practical Use and Limitations|-- 
The SVM models are dependent on the image conditions used to train them.
To use the provided SVM models, please make sure the following conditions
apply:

* Images must be grayscale.

* Images must be resized such that features have the expected size in
pixels as no multiresolution analysis is performed (for reference, the
size of the bottom view detector for the paw is 30x30 pixels). The system
is capable of rescaling, but needs user input on the scaling factor.

* Background must be subtracted.

* Mice must be black. For other colours, consider training a new SVM
model. However, for mice that look white on the image it might not be
possible to distinguish paws from body on the bottom view, which would
break the system.

--|Acknowledgements|-- 
This software was developed by Joao Fayad in the Neural Circuits and
Behavior Lab at the Champalimaud Neuroscience Programme. It is inspired
by previously existing tracking code (unreleased) developed by Ana S.
Machado, Dana Darmohray and Megan R. Carey.

Ana S. Machado, Dana Darmohray, Hugo G. Marques and Megan R. Carey
contributed with discussions, suggestions and testing of the software.
