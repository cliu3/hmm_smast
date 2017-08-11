HMM geolocation toolbox


Installation

Unzip hmmgeoloc.zip file into a directory of your choice. You should be aware that working with the toolbox might create large files (>100 MB for long tracks or large databases) so be sure to have plenty of available disk space. Once the file is unzipped start Matlab and go to the file menu and pick "Set path...". Choose "Add with subfolders.." and find the root directory of the geolocation toolbox (the one you unzipped the files into). This makes Matlab able to "see" the toolbox even if you are working in another directory. Now click "Save" and "Close" in the "Set path..." window and you should be ready to start geolocating some tags.


Running the example

The toolbox comes with an example script file which shows its basic functionality. The example script file is named codexample.m and is stored in the tbworkdir directory. The example is executed in Matlab by typing "codexample" in the prompt. The example reads the database files included in the toolbox and the raw tag data file which is also included. It performs the steps that are always required for geolocating a tag. The script file can be viewed, edited and copied from by typing "edit codexample" in the prompt. The example will typically run very fast so it might be a good idea to run one line at a time to get a feel for what happens in each step. Remember that there is many more functions than the ones used in this basic example.

This example might be the easiest approach to the toolbox, both for setting up a script file but also for getting comfortable with the format of the database and tag files, which are viewable as well, see the section on dataformats in the reference manual.