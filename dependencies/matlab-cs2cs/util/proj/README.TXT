	PROJ 4.4.6 WIN32 Binary Release
	===============================

This is a package of pre-built binaries for PROJ 4.4.6 for Win32 systems.
After the package is unpacked the PROJ\BIN directory should be added to 
your path.

eg. 

C:\> set PATH=%PATH%;C:\PROJ\BIN

If you unpack the distribution somewhere other than C:\, you will also
have to set the PROJ_LIB environment variable so that PROJ can find 
initialization files. 

eg. 

C:\> set PROJ_LIB=C:\Software\PROJ\NAD

Once these steps are complete, you should be able to use the various PROJ
utilities from the commandline. 

eg. 

C:\> proj.exe  -I +proj=utm +zone=11 +ellps=WGS84
332000 4000000                                           (input)
118d52'1.182"W  36d7'48.593"N                            (output)
<cntl-Z>
C:\>

See the PROJ.4 web site for more details on PROJ.

 http://www.remotesensing.org/proj

Those wanting to use PROJ.4 DLLs for other applications like MapServer
can just copy the PROJ\BIN\PROJ.DLL file, and ignore the .EXE files. 
However, the PROJ_LIB environment variable will still need to be set to
enable use of initialization and datum grid shift files. 

License
-------

 Copyright (c) 2000, Frank Warmerdam

 Permission is hereby granted, free of charge, to any person obtaining a
 copy of this software and associated documentation files (the "Software"),
 to deal in the Software without restriction, including without limitation
 the rights to use, copy, modify, merge, publish, distribute, sublicense,
 and/or sell copies of the Software, and to permit persons to whom the
 Software is furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included
 in all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 DEALINGS IN THE SOFTWARE.
