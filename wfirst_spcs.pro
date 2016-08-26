;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;      WFIRST-AFTA SPC (spectr.) noise model - Tyler D. Robinson
;  inputs:
;    Ahr  - hi-res planetary albedo spectrum (ideally gaussian convolved at R=70)
;  lamhr  - wavelength grid for Ahr (um)
;  alpha  - phase angle (deg)
;    Phi  - phase function evaluated at alpha
;     Rp  - planetary radius (Rearth)
;   Teff  - stellar effective temperature (K)
;     Rs  - stellar radius (Rsun)
;      r  - orbital separation (au)
;      d  - distance to system (pc)
;    Nez  - number of exozodis
;
;  outputs:
;    lam  - low-res wavelength grid centers (um)
;   dlam  - spectral bin widths (um)
;      A  - planetary albedo spectrum at low-res
; Cratio  - planet-star contrast (flux) ratio
;     cp  - planet count rate (s**-1)
;    csp  - speckle count rate (s**-1)
;     cz  - zodiacal light count rate (s**-1) 
;    cez  - exozodiacal light count rate (s**-1)
;     cD  - dark current count rate (s**-1)
;     cR  - read noise count rate (s**-1)
;   ctot  - total count rate (s**-1)
;  DtSNR  - integration time to SNR=1 (hr)
;
;  options:
;EMCCD_EM - use EMCCD in electron-multiplying mode
;EMCCD_PC - use EMCCD in photon counting model
;  SILENT - do not print warnings
;
;  updates:
;    Jul 29, 2016 - increased CCD read noise to more realistic value
;    Aug 24, 2016 - updated coronagraph design contrasts to newest simulations
;    Aug 24, 2016 - adjusted photometric aperture to 1 lam/D (Traub et al., 2016)
;    Aug 24, 2016 - updated QE using Figure 4 from Traub et al. 2016
;    Aug 24, 2016 - updated throughput to be non-grey, using Traub et al. (2016)
;    Aug 24, 2016 - updated to use diffuse and point throughputs as functions of lam/D
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO WFIRST_SPCS, Ahr, lamhr, alpha, Phi, Rp, Teff, Rs, r, d, Nez, $
                 lam, dlam, A, Cratio, cp, csp, cz, cez, cD, cR, cC, $
                 ctot, csig, cnoise, DtSNR, EMCCD_EM = emccd_em, $
                 EMCCD_PC = emccd_pc, SILENT = silent

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;   set system parameters  ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  De     = 5.d-4     ; dark current (s**-1)
  diam   = 2.4       ; telescope diameter (m)
  Re     = 3.0       ; read noise per pixel
  Tput   = 0.037     ; system throughput
  theta  = 0.026     ; angular size of lenslet (arcsec)
  IWA    = 2.7       ; inner working angle (lambda/D)
  OWA    = 10.       ; outer working angle (550 nm / D) -- set by SPC, not lenslets --
  Dtmax  = 1.        ; maximum exposure time (hr)
  X      = 1.0       ; size of photometric aperture (lambda/D)
  Cfloor = 1.d-10    ; contrast floor
  lammin = 0.60      ; minimum wavelength (um)
  lammax = 0.97      ; maximum wavelength (um)
  Res    = 70.       ; spectral resolution
  DNhpix = 2         ; horizontal pixel spread of IFS spectrum
  fpa    = 1.        ; fraction of planetary signal in Airy pattern; set to unity
                     ; as the pointsource throughput now includes this
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; adjust parameters if EMCCD mode ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  IF (KEYWORD_SET(emccd_em) OR KEYWORD_SET(emccd_pc)) THEN BEGIN
    De     = 4.d-4     ; dark current (s**-1)
    Re     = 0.2       ; read noise per pixel
    G      = 5.e3      ; gain
    Af     = SQRT(2.)  ; amplification noise factor
    Rc     = 3.e-3     ; clock induced charge (e/pixel/read)
  ENDIF
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; set astrophys parameters ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  MzV  = 23.0 ; zodiacal light surface brightness (mag/arcsec**2)
  MezV = 22.0 ; exozodiacal light surface brightness (mag/arcsec**2)
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;   set wavelength grid    ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  lam  = lammin ;in [um]
  Nlam = 1
  WHILE lam LT lammax DO BEGIN
    lam  = lam + lam/res
    Nlam = Nlam +1
  ENDWHILE 
  lam    = FLTARR(Nlam)
  lam[0] = lammin
  FOR j=1,Nlam-1 DO BEGIN
    lam[j] = lam[j-1] + lam[j-1]/res
  ENDFOR
  dlam = FLTARR(Nlam) ;grid widths (um)
  FOR j=1,Nlam-2 DO BEGIN
    dlam[j] = 0.5*(lam[j+1]+lam[j]) - 0.5*(lam[j-1]+lam[j])
  ENDFOR
  ;widths at edges are same as neighbor
  dlam[0] = dlam[1]
  dlam[Nlam-1] = dlam[Nlam-2]
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; angular separation ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  sep  = r/d*SIN(alpha*!DPI/180.)/(lam*1.e-6/diam*180/!DPI*3600) ;separation in lambda/D
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;   set design contrasts   ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  fn   = './input_data/wfirst_spc_contrast.dat'
  READCOL, fn, lamD, contrast00,  contrast04,  contrast08,  contast16, /SILENT
  C    = DBLARR(Nlam)
  C    = EXP(INTERPOL(ALOG(contrast04),lamD,sep)) ;interpolate data to given separation
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;    check IWA and OWA     ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  i    = WHERE(sep GT OWA) ;points outside OWA
  IF (i[0] NE -1) THEN BEGIN
    C[i] = 1.d0 ;set contrast to cruddy outside OWA
    IF ~KEYWORD_SET(silent) THEN PRINT, 'WARNING: portion of spectrum outside OWA'
  ENDIF
  i    = WHERE(sep LT IWA) ;points inside IWA
  IF (i[0] NE -1) THEN BEGIN
    C[i] = 1.d0 ;set contrast to cruddy inside IWA
    IF ~KEYWORD_SET(silent) THEN PRINT, 'WARNING: portion of spectrum inside IWA'
  ENDIF
 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;   set "all" throughput   ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  fn   = './input_data/wfirst_throughput_all.dat'
  READCOL, fn, lamT, Tall, /SILENT
  T    = DBLARR(Nlam)
  T    = INTERPOL(Tall,lamT,lam) ;interpolate throughput to wavelength grid
  i    = WHERE( T LT 0)
  IF (i[0] NE -1) THEN T[i] = 0.d
 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;  set quantum efficiency  ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  IF (KEYWORD_SET(emccd_em) OR KEYWORD_SET(emccd_pc)) THEN BEGIN
    fn = './input_data/qe_emccd.dat'
    READCOL, fn, lamqe, qe, /SILENT
    q = INTERPOL(qe,lamqe,lam)
    i = WHERE(q LT 0)
    IF( i[0] NE -1) THEN q[i] = 0.
  ENDIF ELSE BEGIN
    fn = './input_data/qe_ccd.dat'
    READCOL, fn, lamqe, qe, /SILENT
    q = INTERPOL(qe,lamqe,lam)
    i = WHERE(q LT 0)
    IF( i[0] NE -1) THEN q[i] = 0.
  ENDELSE

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;    set SPC throughput    ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  fn     = './input_data/wfirst_spc_throughput.dat'
  READCOL, fn, lamD, Td, Tp, /SILENT
  Tdiff  = DBLARR(Nlam)
  Tpoint = DBLARR(Nlam)
  Tdiff  = INTERPOL(Td,lamD,sep) ;interpolate diffuse throughput to given separation
  Tpoint = INTERPOL(Tp,lamD,sep) ;interpolate diffuse throughput to given separation
  i      = WHERE( Tdiff LT 0)
  IF (i[0] NE -1) THEN Tdiff[i] = 0.d
  i      = WHERE( Tpoint LT 0)
  IF (i[0] NE -1) THEN Tpoint[i] = 0.d
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; degrade albedo spectrum  ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  A = DEGRADE_SPEC(Ahr,lamhr,lam,DLAM=dlam)
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;      compute fluxes      ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  Fs = Fstar(lam, Teff, Rs, r, /AU) ;stellar flux on planet (W/m**2/um)
  Fp = Fplan(A, Phi, Fs, Rp, d)     ;planet flux at telescope (W/m**2/um)
  Cratio = FpFs(A, Phi, Rp, r)      ;planet-to-star flux ratio
  IF( MIN(Cratio) LE Cfloor ) THEN $
    IF ~KEYWORD_SET(silent) THEN PRINT, 'WARNING: portions of spectrum are below contrast floor'
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;    compute count rates   ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  cp     =  cplan(q, fpa, T*Tpoint, lam, dlam, Fp, diam)
  cz     =  czodi(q, X, T*Tdiff, lam, dlam, diam, MzV)
  cez    = cezodi(q, X, T*Tdiff, lam, dlam, diam, r, Fstar(lam,Teff,Rs,1.,/AU), Nez, MezV)
  csp    = cspeck(q, T*Tpoint, C, lam, dlam, Fstar(lam,Teff,Rs,d), diam)
  cD     =  cdark(De, X, lam, diam, theta, DNhpix)
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; if not using EMCCD ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  IF ~(KEYWORD_SET(emccd_em) OR KEYWORD_SET(emccd_pc)) THEN BEGIN
    cR     =  cread(Re, X, lam, diam, theta, DNHpix, Dtmax)
    ctot   = cp + cz + cez + csp + cD + cR
    csig   = cp
    cnoise = cp + 2.*(cz + cez + csp + cD + cR) ;assumes roll for background subtraction
  ENDIF
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; if using photon counting mode ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  IF KEYWORD_SET(emccd_pc) THEN BEGIN
    max_c  = MAX([cp,cz,cez,csp]) ;maximum count rate (1/s)
    max_c  = max_c*3600.          ;convert to 1/hr
    IF (max_c GT 0) THEN Dtmax  = 1.d/max_c/10.      ;read at 10x maximum count rate
    cR     =  cread(Re, X, lam, diam, theta, DNHpix, Dtmax)
    cC     = G*ccic(Rc, X, lam, diam, theta, DNHpix, Dtmax)
    cp     = G*cp
    cz     = G*cz
    cez    = G*cez
    csp    = G*csp
    cD     = G*cD
    ctot   = cp + cz + cez + csp + cD + cR + cC
    csig   = cp
    cnoise = (G*cp + 2.*(cR + G*(cz + cez + csp + cD + cC)))
  ENDIF
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; if using electron multiplying mode ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  IF KEYWORD_SET(emccd_em) THEN BEGIN
    cR     =  cread(Re, X, lam, diam, theta, DNHpix, Dtmax)
    cC     = G*ccic(Rc, X, lam, diam, theta, DNHpix, Dtmax)
    cp     = G*cp
    cz     = G*cz
    cez    = G*cez
    csp    = G*csp
    cD     = G*cD
    ctot   = cp + cz + cez + csp + cD + cR + cC
    csig   = cp
    cnoise = (G*Af^2.*cp + 2.*(cR + G*Af^2.*(cz + cez + csp + cD + cC)))
  ENDIF
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;  exposure time to SNR=1  ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  DtSNR = DBLARR(Nlam)
  DtSNR[*] = 0.d
  i     = WHERE(csig GT 0)
  IF (i[0] NE -1) THEN DtSNR[i] = cnoise[i]/csig[i]^2./3600. ; (hr)
  
END
