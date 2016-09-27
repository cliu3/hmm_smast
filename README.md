# hmm_smast
`hmm_smast` is a MATLAB package for geolocating New England groundfish using Data Storage Tag (DST) data. This package is a fork of Martin W Pedersen's MATLAB geolocation toolbox, available at http://mwpedersen.dk/tracking.html.

## Dependencies
This code requires [t_tide](https://www.eoas.ubc.ca/~rich/) to perform tidal elevation prediction and the [Matlab Google Earth Toolbox](http://www.mathworks.com/matlabcentral/fileexchange/12954-google-earth-toolbox) to plot the most probable track in .kml format. These packages are included in the `dependencies` directory.

## Documentation

View the documentation generated by m2html [here](http://htmlpreview.github.io/?https://github.com/cliu3/hmm_smast/blob/dev/doc/index.html).

## Preprocessing
As an alternative to the preprocessing code in MATLAB provided in [test/preprocessing], an R code package for preprocessing raw DST data for the use of hmm_smast is available at https://github.com/cliu3/R_HMM_Preprocessing. 

## Publications
* Liu, C., Cowles, G., Zemeckis, D.R., Cadrin, S.X, and Dean, M.J. *In review*. **Developing and validating geolocation methods based on hidden Markov Models for groundfish species off New England**.
* Pedersen, M.W., Righton, D., Thygesen, U.H., Andersen, K.H., Madsen, H., 2008. **Geolocation of North Sea cod (*Gadus morhua*) using hidden Markov models and behavioural switching**. Canadian Journal of Fisheries and Aquatic Sciences 65, 2367–2377.
