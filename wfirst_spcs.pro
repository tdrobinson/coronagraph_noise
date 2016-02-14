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
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO WFIRST_SPCS, Ahr, lamhr, alpha, Phi, Rp, Teff, Rs, r, d, Nez, $
                 lam, dlam, A, Cratio, cp, csp, cz, cez, cD, cR, ctot, $
                 DtSNR

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;   set system parameters  ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  De     = 5.d-4     ; dark current (s**-1)
  diam   = 2.4       ; telescope diameter (m)
  Re     = 0.2       ; read noise per pixel
  Tput   = 0.037     ; system throughput
  theta  = 0.017     ; angular size of lenslet (arcsec)
  IWA    = 2.7       ; inner working angle (lambda/D)
  OWA    = 10.       ; outer working angle (550 nm / D) -- set by SPC, not lenslets --
  Dtmax  = 1.        ; maximum exposure time (hr)
  X      = 1.5       ; size of photometric aperture (lambda/D)
  Cfloor = 1.d-10    ; contrast floor
  lammin = 0.60      ; minimum wavelength (um)
  lammax = 0.97      ; maximum wavelength (um)
  Res    = 70.       ; spectral resolution
  DNhpix = 3         ; horizontal pixel spread of IFS spectrum
  fpa    = f_airy(X) ; fraction of planetary signal in Airy pattern
  
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
  ;;;   set mask throughput    ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  T    = DBLARR(Nlam)
  T[*] = Tput
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;    check IWA and OWA     ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  i    = WHERE(sep GT OWA) ;points outside OWA
  IF (i[0] NE -1) THEN BEGIN
    C[i] = 1.d0 ;set contrast to cruddy outside OWA
    PRINT, 'WARNING: portion of spectrum outside OWA'
  ENDIF
  i    = WHERE(sep LT IWA) ;points inside IWA
  IF (i[0] NE -1) THEN BEGIN
    C[i] = 1.d0 ;set contrast to cruddy inside IWA
    PRINT, 'WARNING: portion of spectrum inside IWA'
  ENDIF
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;  set quantum efficiency  ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  q = FLTARR(Nlam)
  FOR j=0,Nlam-1 DO BEGIN
    IF (lam[j] LE 0.7) THEN BEGIN
      q[j] = 0.9
    ENDIF ELSE BEGIN
      q[j] = 0.9*(1.0 - (lam[j]-0.7)/(1.0-0.7))
    ENDELSE
    IF q[j] LT 0 THEN q[j] = 0.
  ENDFOR
  
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
    PRINT, 'WARNING: portions of spectrum are below contrast floor'
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;    compute count rates   ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  cp     =  cplan(q, fpa, T, lam, dlam, Fp, diam)
  cz     =  czodi(q, X, T, lam, dlam, diam, MzV)
  cez    = cezodi(q, X, T, lam, dlam, diam, r, Fstar(lam,Teff,Rs,1.,/AU), Nez, MezV)
  csp    = cspeck(q, T, C, lam, dlam, Fstar(lam,Teff,Rs,d), diam)
  cD     =  cdark(De, X, lam, diam, theta, DNhpix)
  cR     =  cread(Re, X, lam, diam, theta, DNHpix, Dtmax)
  ctot   = cp + cz + cez + csp + cD + cR
  cnoise = cp + 2*(cz + cez + csp + cD + cR) ;assumes roll for background subtraction
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;  exposure time to SNR=1  ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  DtSNR = DBLARR(Nlam)
  DtSNR[*] = 0.d
  i     = WHERE(cp GT 0)
  IF (i[0] NE -1) THEN DtSNR[i] = cnoise[i]/cp[i]^2./3600. ; (hr)
  
END
