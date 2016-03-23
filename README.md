# TrackNTrace
--------------

TrackNTrace is an open source MATLAB framework for single molecule localization, tracking, and super-resolution applications written by Simon Christoph Stein and Jan Thiart from the University of Goettingen.

The software facilitates development, distribution, and comparison of methods in the community by providing a readily extendable, plugin-based system and combining it with an easy-to-use graphical user interface (GUI). This GUI incorporates possibilities for quick inspection of localization and tracking results, giving direct feedback of the quality achieved with the chosen algorithms and parameter values, as well as possible errors, a feature neglected in most software packages available. The plugin system simplifies adapting and tailoring methods towards any research problem's individual requirements. We provide a set of plugins implementing state-of-the-art methods together with the basic program, alongside tools for common post-processing steps such as (d)STORM image generation, or drift correction.

To get started with using TrackNTrace, have a look at the PDF file inside the manual folder.
A publication introducing the TrackNTrace framework is currently in review.


## Licensing
--------------
(See License.txt for more information)

If not stated otherwise, all files part of TrackNTrace are under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
You can find a copy of the license at <http://www.gnu.org/licenses/>.

If not stated otherwise, the following copyright applies:
 Copyright (C) 2016  Simon Christoph Stein, scstein@phys.uni-goettingen.de
 Copyright (C) 2016  Jan Thiart, jthiart@phys.uni-goettingen.de

Plugin files/(sub)functions (mostly in the plugins subfolder) might come with their own license.

TrackNTrace uses variety of external programs / modules:

* distinguishable_colors.m, 
	Copyright 2010-2011 Timothy E. Holy, 
	from the MATLAB file exchange, 
	BSD license
* u-track, 
	Copyright (C) 2014 LCCB, 
	URL: http://lccb.hms.harvard.edu/software.html, 
	GPLv3 license
* GPUGauss_MLEv2, 
	Copyright (C) 2011 Peter Relich
	URL: http://omictools.com/gaussmlev2-tool, 
	GPLv3 license
* nanoflann library, 
	Copyright 2008-2009  Marius Muja (mariusm@cs.ubc.ca). All rights reserved.
	Copyright 2008-2009  David G. Lowe (lowe@cs.ubc.ca). All rights reserved.
	Copyright 2011-2013  Jose Luis Blanco (joseluisblancoc@gmail.com).
	License see 2) below
* ceres solver library	
	Copyright 2016 Google Inc. All rights reserved.
	License see 3) below
* "FastPsfFitting" and "NearestNeighborTracker",
    Copyright (c) 2016 Simon Christoph Stein and Jan Thiart,
    Free BSD license (see 1) below))

