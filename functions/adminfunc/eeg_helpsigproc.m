%EEGLAB signal processing functions (sigprocfunc folder):
%  <a href="matlab:helpwin acsobiro">acsobiro</a>             - A. Chickocki's robust Second-Order Blind Identification (SOBI)...
%  <a href="matlab:helpwin axcopy">axcopy</a>               - Copy a Matlab figure axis and its graphic objects to a new pop-up window...
%  <a href="matlab:helpwin binica">binica</a>               - Run stand-alone binary version of RUNICA from the...
%  <a href="matlab:helpwin biosig2eeglab">biosig2eeglab</a>        - Convert BIOSIG structure to EEGLAB structure...
%  <a href="matlab:helpwin biosig2eeglabevent">biosig2eeglabevent</a>   - Convert biosig events to EEGLAB event structure...
%  <a href="matlab:helpwin blockave">blockave</a>             - Make block average of concatenated data sets of same size...
%  <a href="matlab:helpwin cart2topo">cart2topo</a>            - Convert xyz-cartesian channel coordinates...
%  <a href="matlab:helpwin cbar">cbar</a>                 - Display full or partial color bar...
%  <a href="matlab:helpwin celltomat">celltomat</a>            - Convert cell array to matrix...
%  <a href="matlab:helpwin chancenter">chancenter</a>           - Recenter cartesian X,Y,Z channel coordinates...
%  <a href="matlab:helpwin changeunits">changeunits</a>          - Takes one or more points in one axes and gives its position...
%  <a href="matlab:helpwin compvar">compvar</a>              - Project selected components and compute the variance of...
%  <a href="matlab:helpwin condstat">condstat</a>             - Accumulate surrogate data for comparing two data conditions...
%  <a href="matlab:helpwin convertlocs">convertlocs</a>          - Convert electrode locations between coordinate systems...
%  <a href="matlab:helpwin copyaxis">copyaxis</a>             - Helper function for AXCOPY...
%  <a href="matlab:helpwin coregister">coregister</a>           - Co-register measured or template electrode locations with a...
%  <a href="matlab:helpwin dipoledensity">dipoledensity</a>        - Compute and optionally plot a measure of the 3-D spatial...
%  <a href="matlab:helpwin eegfilt">eegfilt</a>              - (high|low|band)-pass filter data using two-way least-squares...
%  <a href="matlab:helpwin eegfiltfft">eegfiltfft</a>           - (high|low|band)-pass filter data using inverse fft...
%  <a href="matlab:helpwin eegplot">eegplot</a>              - Scroll (horizontally and/or vertically) through multichannel data.
%  <a href="matlab:helpwin eegplot2event">eegplot2event</a>        - Convert EEGPLOT rejections into events...
%  <a href="matlab:helpwin eegplot2trial">eegplot2trial</a>        - Convert EEGPLOT rejections into trial and electrode...
%  <a href="matlab:helpwin eegplot_readkey">eegplot_readkey</a>      - Eegplot helper function to read key strokes...
%  <a href="matlab:helpwin eegrej">eegrej</a>               - Reject/excise arbitrary periods from continuous EEG data...
%  <a href="matlab:helpwin eegthresh">eegthresh</a>            - Reject trials with out-of-bounds channel values within a...
%  <a href="matlab:helpwin entropy_rej">entropy_rej</a>          - Calculation of entropy of a 1D, 2D or 3D array and...
%  <a href="matlab:helpwin env">env</a>                  - Return envelope of rows of a data matrix, or optionally...
%  <a href="matlab:helpwin envtopo">envtopo</a>              - No help information...
%  <a href="matlab:helpwin epoch">epoch</a>                - Extract epochs time locked to specified events from continuous EEG data.
%  <a href="matlab:helpwin erpimage">erpimage</a>             - Plot a colored image of a collection of single-trial data epochs, optionally...
%  <a href="matlab:helpwin eventalign">eventalign</a>           - Function called by POP_IMPORTEVENT to find the best...
%  <a href="matlab:helpwin eventlock">eventlock</a>            - DEPRECATED: Please use EEGALIGN instead.
%  <a href="matlab:helpwin eyelike">eyelike</a>              - Calculate a permutation matrix P and a scaling (diagonal) maxtrix S...
%  <a href="matlab:helpwin fastif">fastif</a>               - Fast if function.
%  <a href="matlab:helpwin floatread">floatread</a>            - Read matrix from float file ssuming four byte floating point number...
%  <a href="matlab:helpwin floatwrite">floatwrite</a>           - Write data matrix to float file.
%  <a href="matlab:helpwin forcelocs">forcelocs</a>            - Rotate location in 3-D so specified electrodes...
%  <a href="matlab:helpwin gettempfolder">gettempfolder</a>        - Return the temporary folder...
%  <a href="matlab:helpwin headplot">headplot</a>             - Plot a spherically-splined EEG field map on a semi-realistic...
%  <a href="matlab:helpwin icaact">icaact</a>               - Compute ICA activation waveforms = weights*sphere*(data-meandata)...
%  <a href="matlab:helpwin icadefs">icadefs</a>              - Function to read in a set of EEGLAB system-wide (i.e. lab-wide)...
%  <a href="matlab:helpwin icaproj">icaproj</a>              - Project ICA component activations through the...
%  <a href="matlab:helpwin icavar">icavar</a>               - Project ICA component activations through the ICA weight matrices...
%  <a href="matlab:helpwin imagesctc">imagesctc</a>            - DEPRECATED. never completed or documented.
%  <a href="matlab:helpwin isscript">isscript</a>             - Function checking if a specific file is a script or a function...
%  <a href="matlab:helpwin jader">jader</a>                - Blind separation of real signals using JADE (v1.5, Dec. 1997).
%  <a href="matlab:helpwin jointprob">jointprob</a>            - Rejection of odd columns of a data array  using...
%  <a href="matlab:helpwin kmeanscluster">kmeanscluster</a>        - Simple k means clustering algorithm...
%  <a href="matlab:helpwin kurt">kurt</a>                 - Return kurtosis of input data distribution...
%  <a href="matlab:helpwin loadtxt">loadtxt</a>              - Load ascii text file into numeric or cell arrays...
%  <a href="matlab:helpwin lookupchantemplate">lookupchantemplate</a>   - Look up channel template.
%  <a href="matlab:helpwin matsel">matsel</a>               - Select rows, columns, and epochs from given multi-epoch data matrix...
%  <a href="matlab:helpwin mattocell">mattocell</a>            - Convert matrix to cell array...
%  <a href="matlab:helpwin metaplottopo">metaplottopo</a>         - Plot concatenated multichannel data epochs in a topographic or...
%  <a href="matlab:helpwin movav">movav</a>                - Perform a moving average of data indexed by xvals.
%  <a href="matlab:helpwin moveaxes">moveaxes</a>             - Move, resize, or copy Matlab axes using the mouse...
%  <a href="matlab:helpwin mri3dplot">mri3dplot</a>            - Plot 3-D density image translucently on top of the mean MR...
%  <a href="matlab:helpwin nan_mean">nan_mean</a>             - Average, not considering NaN values...
%  <a href="matlab:helpwin parsetxt">parsetxt</a>             - Parse text input into cell array...
%  <a href="matlab:helpwin phasecoher">phasecoher</a>           - Implements inter-trial amp/coherence using Gaussian wavelets.
%  <a href="matlab:helpwin plotchans3d">plotchans3d</a>          - Plots the 3-D configuration from a Polhemus ELP...
%  <a href="matlab:helpwin plotcurve">plotcurve</a>            - Plot curve(s) with optional significance highlighting.
%  <a href="matlab:helpwin plotdata">plotdata</a>             - Plot concatenated multichannel data epochs in two-column format...
%  <a href="matlab:helpwin ploterp">ploterp</a>              - Plot a selected multichannel data epoch on a common timebase...
%  <a href="matlab:helpwin plotmesh">plotmesh</a>             - Plot mesh defined by faces and vertex...
%  <a href="matlab:helpwin plotsphere">plotsphere</a>           - This function is used to plot a sphere and...
%  <a href="matlab:helpwin plottopo">plottopo</a>             - Plot concatenated multichannel data epochs in a topographic...
%  <a href="matlab:helpwin posact">posact</a>               - Make RUNICA activations all RMS-positive.
%  <a href="matlab:helpwin projtopo">projtopo</a>             - Plot projections of one or more ICA components along with...
%  <a href="matlab:helpwin qqdiagram">qqdiagram</a>            - Empirical quantile-quantile diagram.
%  <a href="matlab:helpwin quantile">quantile</a>             - Computes the quantiles of the data sample from a distribution X...
%  <a href="matlab:helpwin readedf">readedf</a>              - Read eeg data in EDF format.
%  <a href="matlab:helpwin readeetraklocs">readeetraklocs</a>       - Read 3-D location files saved using the EETrak...
%  <a href="matlab:helpwin readegi">readegi</a>              - Read EGI Simple Binary datafile (versions 2,3,4,5,6,7).
%  <a href="matlab:helpwin readegihdr">readegihdr</a>           - Read header information from EGI (versions 2,3,4,5,6,7) data file.
%  <a href="matlab:helpwin readegilocs">readegilocs</a>          - Look up locations for EGI EEG dataset.
%  <a href="matlab:helpwin readelp">readelp</a>              - Read electrode locations from an .elp (electrode positions)...
%  <a href="matlab:helpwin readlocs">readlocs</a>             - Read electrode location coordinates and other information from a file.
%  <a href="matlab:helpwin readtxtfile">readtxtfile</a>          - Read text file into a Matlab variable...
%  <a href="matlab:helpwin realproba">realproba</a>            - Compute the effective probability of the value...
%  <a href="matlab:helpwin rejkurt">rejkurt</a>              - Calculation of kutosis of a 1D, 2D or 3D array and...
%  <a href="matlab:helpwin rejstatepoch">rejstatepoch</a>         - Reject bad eeg trials based a statistical measure. Can...
%  <a href="matlab:helpwin rejtrend">rejtrend</a>             - Detect linear trends in EEG activity and reject the...
%  <a href="matlab:helpwin reref">reref</a>                - Convert common reference EEG data to some other common reference...
%  <a href="matlab:helpwin rmbase">rmbase</a>               - Subtract basevector channel means from multi-epoch data matrix...
%  <a href="matlab:helpwin runica">runica</a>               - Perform Independent Component Analysis (ICA) decomposition...
%  <a href="matlab:helpwin runica_ml">runica_ml</a>            - Perform Independent Component Analysis (ICA) decomposition...
%  <a href="matlab:helpwin runica_ml2">runica_ml2</a>           - Perform Independent Component Analysis (ICA) decomposition...
%  <a href="matlab:helpwin runica_mlb">runica_mlb</a>           - Perform Independent Component Analysis (ICA) decomposition...
%  <a href="matlab:helpwin sbplot">sbplot</a>               - Create axes in arbitrary subplot grid positions and sizes...
%  <a href="matlab:helpwin shuffle">shuffle</a>              - Shuffle a given dimension in an array...
%  <a href="matlab:helpwin signalstat">signalstat</a>           - Computes and plots statistical characteristics of a signal,...
%  <a href="matlab:helpwin slider">slider</a>               - Add slider to a figure...
%  <a href="matlab:helpwin snapread">snapread</a>             - Read data in Snap-Master Standard Binary Data File Format...
%  <a href="matlab:helpwin sobi">sobi</a>                 - Second Order Blind Identification (SOBI) by joint diagonalization of...
%  <a href="matlab:helpwin spec">spec</a>                 - Power spectrum. This function replaces PSD function if the signal...
%  <a href="matlab:helpwin spectopo">spectopo</a>             - Plot the power spectral density (PSD) of winsize length segments of data...
%  <a href="matlab:helpwin sph2topo">sph2topo</a>             - Convert from a 3-column headplot file in spherical coordinates...
%  <a href="matlab:helpwin spher">spher</a>                - Return the sphering matrix for given input data...
%  <a href="matlab:helpwin spherror">spherror</a>             - CHANCENTER sub function to compute minimum distance...
%  <a href="matlab:helpwin strmultiline">strmultiline</a>         - Format a long string as a multi-line string.
%  <a href="matlab:helpwin textsc">textsc</a>               - Places text in screen coordinates and places...
%  <a href="matlab:helpwin timefdetails">timefdetails</a>         - Details of the TIMEF function for time/frequency analysis...
%  <a href="matlab:helpwin timtopo">timtopo</a>              - Plot all channels of a data epoch on the same axis...
%  <a href="matlab:helpwin topo2sph">topo2sph</a>             - Convert a TOPOPLOT style 2-D polar-coordinate...
%  <a href="matlab:helpwin topoplot">topoplot</a>             - Plot a topographic map of a scalp data field in a 2-D circular view...
%  <a href="matlab:helpwin transformcoords">transformcoords</a>      - Select nazion and inion in anatomical MRI images.
%  <a href="matlab:helpwin trial2eegplot">trial2eegplot</a>        - Convert eeglab format to eeplot format of rejection window...
%  <a href="matlab:helpwin uigetfile2">uigetfile2</a>           - Same as uigetfile but remember folder location.
%  <a href="matlab:helpwin uiputfile2">uiputfile2</a>           - Same as uigputfile but remember folder location.
%  <a href="matlab:helpwin writelocs">writelocs</a>            - Write a file containing channel location, type and gain information...
