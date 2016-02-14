# coronagraph_noise
Coronagraph noise modeling routines

The enclosed IDL routines are coronagraph noise simulators based on the models 
described in: 
	Robinson, T.D., Stapelfeldt, K.R., and Marley, M.S. (2016). "Characterizing Rocky 
	and Gaseous Exoplanets with 2m Class Space-based Coronagraphs." Proceedings of 
	the Astronomical Society of the Pacific, 128:025003

If you use these tools, please cite this paper.  Any questions can be sent to:
	Tyler Robinson (robinson.tyler.d@gmail.com)

Code is extensively commented, and symbols names generally follow those in the 
aforementioned manuscript.

Note that the phase function parameter ("Phi") is a common source of misunderstanding.  
This paramater is meant to be a gray adjustment to the albedo spectrum (supplied in 
"Ahr").  If the albedo spectrum already includes phase effects, then simply set Phi 
equal to unity.  However, if Ahr is, say, a geometric albedo spectrum, then one should 
set Phi to, for example, the Lambert phase function value at the desired phase angle 
("alpha").  The simulators apply no internal phase correction, and the phase angle is 
only used to compute the planet-star separation.

The simulators are:
	WFIRST_SPCS - WFIRST Shaped Pupil Coronagraph in "spectroscopy" mode
	WFRIST_SPCI - WFIRST Shaped Pupil Coronagraph in "imaging" mode
	WFRIST_HLC  - WFIRST Hybrid Lyor Coronagraph
    CORONAGRAPH - a generalized coronagraph routine

Also included is an example that calls the different simulators.  The example case is 
a 2-AU Jupiter model (Cahoy et al. 2010, ApJ 724:189) around a Sun-like star at 6 pc.  
Spectra with faux noise are also created, assuming a 20 hr integration per coronagraph 
bandpass.  A figure is generated showing the noise-free and noisy results for all the 
different simulators.