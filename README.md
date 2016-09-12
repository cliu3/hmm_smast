# hmm_smast
`hmm_smast` is a MATLAB package for geolocating New England groundfish using Data Storage Tag (DST) data. This package is a fork of Martin W Pedersen's MATLAB geolocation toolbox, available at http://mwpedersen.dk/tracking.html.

## Dependencies
This code requires MATLAB Mapping Toolbox.

This code requires [t_tide](https://www.eoas.ubc.ca/~rich/) to perform tidal elevation prediction and the [Matlab Google Earth Toolbox](http://www.mathworks.com/matlabcentral/fileexchange/12954-google-earth-toolbox) to plot the most probable track in .kml format. Please [follow the instruction here](dependencies/README.md) to obtain these two packages and put them in the `dependencies` directory.

## Preprocessing
An R code package for preprocessing raw DST data for the use of hmm_smast is available at https://github.com/cliu3/R_HMM_Preprocessing. 

## Publications
* Liu, C., Cowles, G., Zemeckis, D.R., Cadrin, S.X, and Dean, M.J. *In prep*. **Developing and validating geolocation methods based on hidden Markov Models for groundfish species off New England**.
* Pedersen, M.W., Righton, D., Thygesen, U.H., Andersen, K.H., Madsen, H., 2008. **Geolocation of North Sea cod (*Gadus morhua*) using hidden Markov models and behavioural switching**. Canadian Journal of Fisheries and Aquatic Sciences 65, 2367–2377.
