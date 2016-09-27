cs2cs for Matlab
================

This is a wrapper-function to call 
[cs2cs](http://proj.maptools.org/man_cs2cs.html) from 
[GNU Octave](https://www.gnu.org/software/octave/) or 
[Mathworks Matlab](http://www.mathworks.com/products/matlab/).


Summary
-------
Cs2cs performs transformation between the source and destination cartographic 
coordinate system on a set of input points. The coordinate system 
transformation can include translation between projected and geographic 
coordinates as well as the application of datum shifts.


Example Usage
-------------
<pre>
prj4_params = ['-f  "%.12f" +proj=tmerc +lat_0=0 +lon_0=34 +k=1 +x_0=0 +y_0=0 +ellps=bessel ' ...
               '+towgs84=577.326,90.129,463.919,5.137,1.474,5.297,2.4232 +units=m +no_defs'];
x = [-1300043.7691, -1300143.76919];
y = [5484682.0535, 5484382.0535];
[lon, lat] = cs2cs(x, y, prj4_params);
</pre>


Background
----------
Arguments are passed to the binary of cs2cs using temporary files. This may be 
slow for too many function calls. 


Dependencies
------------
Linux:
Please install the package proj-bin which contains cs2cs:
<pre>
aptitude install proj-bin
</pre>

Windows:
The binaries of cs2cs are shipped with this release.


License
-------
Copyright (c) 2013, Erwin Nindl<br/>
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
- Neither the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Contact
-------
Please use the contact-form provided by github.
