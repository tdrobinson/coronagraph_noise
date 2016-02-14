PRO EXAMPLE
  
  ;;;;;;;;;;;;;;;;;;;;;
  ;;;  reflectance  ;;;
  ;;;;;;;;;;;;;;;;;;;;;
  
  ; read in 2 AU jupiter reflectance spectrum from albedo code
  fn = 'jupiter_2.0AU_3x_90deg.alb'
  READCOL, fn, k, lamhr, a, Ahr, a, a, a
  
  ; convolve with R=70 gaussian lineshape function
  res = 70.
  Ahr = GAUSSFOLD(ALOG(lamhr),Ahr,1/res)
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; planet/orbit properties ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  alpha = 90.    ; phase angle of model (deg)
  Phi   = 1.     ; set phase function adjustment to unity, as phase effects already in model
  Rp    = 11.2   ; planet radius, equals Jupiter radius in Earth radii
  r     = 2.0    ; orbital distance (au)
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; star/system properties ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  Teff  = 5778.   ; effective temperature (K)
  Rs    = 1.0     ; radius in solar radii
  d     = 6.0     ; distance to system (pc)
  Nez   = 1.      ; number of exo-zodi
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; call WFIRST noise models ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  WFIRST_SPCS, Ahr, lamhr, alpha, Phi, Rp, Teff, Rs, r, d, Nez, $
               lam_spcs, dlam_spcs, A_spcs, Cratio_spcs, cp_spcs, $
               csp_spcs, cz_spcs, cez_spcs, cD_spcs, cR_spcs, $
               ctot_spcs, DtSNR_spcs
  WFIRST_SPCI, Ahr, lamhr, alpha, Phi, Rp, Teff, Rs, r, d, Nez, $
               lam_spci, dlam_spci, A_spci, Cratio_spci, cp_spci, $
               csp_spci, cz_spci, cez_spci, cD_spci, cR_spci, $
               ctot_spci, DtSNR_spci
  WFIRST_HLC,  Ahr, lamhr, alpha, Phi, Rp, Teff, Rs, r, d, Nez, $
               lam_hlc, dlam_hlc, A_hlc, Cratio_hlc, cp_hlc, $
               csp_hlc, cz_hlc, cez_hlc, cD_hlc, cR_hlc, $
               ctot_hlc, DtSNR_hlc

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; call generalized noise model ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  lammin = 0.4    ; minimum wavelength (um)
  lammax = 1.0    ; maximum wavelength (um)
  Res    = 70.    ; spectral resolution
  X      = 1.5    ; photometric aperture size (lambda/D)
  diam   = 2.4    ; telescope diameter (m)
  Tput   = 0.05   ; system throughput
  cont   = 2.e-10 ; design raw contrast
  IWA    = 2.8    ; inner working angle (lambda/D)
  OWA    = 10.    ; outer working angle (lambda/D)
  De     = 5.d-4  ; dark current (s**-1)
  DNhpix = 3      ; horizontal pixel spread of IFS spectrum
  Re     = 0.1    ; read noise per pixel
  Dtmax  = 1.     ; maximum exposure time (hr)
  CORONAGRAPH, Ahr, lamhr, alpha, Phi, Rp, Teff, Rs, r, d, Nez, lammin, lammax, $
               Res, X, diam, Tput, cont, IWA, OWA, De, DNhpix, Re, Dtmax, $
               lam, dlam, A, q, Cratio, cp, csp, cz, cez, cD, cR, $
               ctot, DtSNR, /COMPUTE_LAM, /COMPUTE_DLAM

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;; example spectra w/noise  ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  Dt    = 20.       ; integration time (hr)
  
  ; WFIRST instruments
  SNR_spcs   = cp_spcs*Dt*3600./SQRT((cp_spcs + $
               2*(csp_spcs+cz_spcs+cez_spcs+cD_spcs+cR_spcs))*Dt*3600.)    ; signal-to-noise
  sigma_spcs = Cratio_spcs/SNR_spcs ; standard deviation
  cont_spcs  = Cratio_spcs + RANDOMN(seed,N_ELEMENTS(lam_spcs))*sigma_spcs ; contrast spectrum w/random noise
  i    = WHERE(cont_spcs LT 0) ; any points with contrast < 0 are set to zero
  IF(i[0] NE -1) THEN cont_spcs[i] = 0.
  
  SNR_spci   = cp_spci*Dt*3600./SQRT((cp_spci + $
               2*(csp_spci+cz_spci+cez_spci+cD_spci+cR_spci))*Dt*3600.)    ; signal-to-noise
  sigma_spci = Cratio_spci/SNR_spci ; standard deviation
  cont_spci  = Cratio_spci + RANDOMN(seed,N_ELEMENTS(lam_spci))*sigma_spci ; contrast spectrum w/random noise
  i    = WHERE(cont_spci LT 0) ; any points with contrast < 0 are set to zero
  IF(i[0] NE -1) THEN cont_spci[i] = 0.
  
  SNR_hlc    = cp_hlc*Dt*3600./SQRT((cp_hlc + $
               2*(csp_hlc+cz_hlc+cez_hlc+cD_hlc+cR_hlc))*Dt*3600.)    ; signal-to-noise
  sigma_hlc  = Cratio_hlc/SNR_hlc ; standard deviation
  cont_hlc   = Cratio_hlc + RANDOMN(seed,N_ELEMENTS(lam_hlc))*sigma_hlc ; contrast spectrum w/random noise
  i    = WHERE(cont_hlc LT 0) ; any points with contrast < 0 are set to zero
  IF(i[0] NE -1) THEN cont_hlc[i] = 0.
  
  ; generalized model
  SNR   = cp*Dt*3600./SQRT((cp + 2*(csp+cz+cez+cD+cR))*Dt*3600.) ; signal-to-noise
  sigma = Cratio/SNR ; standard deviation
  cont  = Cratio + RANDOMN(seed,N_ELEMENTS(lam))*sigma ; contrast spectrum w/random noise
  i    = WHERE(cont LT 0) ; any points with contrast < 0 are set to zero
  IF(i[0] NE -1) THEN cont[i] = 0.  

  ;;;;;;;;;;;;
  ;;; plot ;;;
  ;;;;;;;;;;;;
  
  ; plot parameters
  cth = 5.0   ;character thickness
  csi = 1.0   ;character size
  th  = 5.0   ;line thickness
  v = FINDGEN(17) * (!PI*2/16.)
  USERSYM, COS(v), SIN(v), /FILL
  
  ; make plot
  SET_PLOT, 'ps'
  DEVICE, FILENAME='wfirst_example_cool_jupiter_20hr.eps', /ENCAPSULATED
  DEVICE, /COLOR
  LOADCT, 0
  PLOT, lam, Cratio*1.e9, XRANGE=[0.4,1.0], XSTYLE=1, YSTYLE=1, $
        YRANGE=[0,15.], XTHICK=th, YTHICK=th, CHARSIZE=csi, $
        XTITLE=TeXtoIDL('Wavelength (\mum)'), CHARTHICK=cth, $
        YTITLE=TeXtoIDL('Planet-Star Flux Ratio \times 10^{9}'), $
        POSITION=[0.13,0.14,0.97,0.98], /NODATA
  LOADCT, 15
  OPLOT, lam, Cratio*1.e9, THICK=th, COLOR=40, PSYM=10
  LOADCT, 0
  OPLOTERROR, lam, cont*1.e9, sigma*1.e9, THICK=0.75*th, ERRCOLOR=125, PSYM=3
  OPLOTERROR, lam_spcs, cont_spcs*1.e9, sigma_spcs*1.e9, THICK=0.75*th, ERRCOLOR=125, PSYM=3
  OPLOTERROR, lam_spci, cont_spci*1.e9, sigma_spci*1.e9, THICK=0.75*th, ERRCOLOR=125, PSYM=3
  LOADCT, 15
  OPLOTERROR, lam_spci, cont_spci*1.e9, dlam_spci/2., sigma_spci*1.e9, /NOHAT, THICK=0.75*th, ERRCOLOR=75, PSYM=3
  LOADCT, 0
  OPLOTERROR, lam_hlc, cont_hlc*1.e9, sigma_hlc*1.e9, THICK=0.75*th, ERRCOLOR=125, PSYM=3
  LOADCT, 15
  OPLOTERROR, lam_hlc, cont_hlc*1.e9, dlam_hlc/2., sigma_hlc*1.e9, /NOHAT, THICK=0.75*th, ERRCOLOR=125, PSYM=3
  LOADCT, 15
  OPLOT, lam, cont*1.e9, PSYM=8, SYMSIZE=0.1*th, COLOR=175
  OPLOT, lam_spcs, cont_spcs*1.e9, PSYM=8, SYMSIZE=0.1*th, COLOR=75
  OPLOT, lam_spci, cont_spci*1.e9, PSYM=8, SYMSIZE=0.1*th, COLOR=75
  OPLOT, lam_hlc, cont_hlc*1.e9, PSYM=8, SYMSIZE=0.1*th, COLOR=125
  LEGEND, ['SPC','HLC','General','Noise-Free'], LINESTYLE=[0,0,0,0], COLOR=[75,125,175,40], THICK=th, CHARSIZE=1.5*csi, $
          CHARTHICK=cth, /BOTTOM, /LEFT, BOX=0
  DEVICE, /CLOSE
  SET_PLOT, 'X'
  LOADCT, 0
  
END