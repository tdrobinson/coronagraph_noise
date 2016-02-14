;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        General coronagraph noise model - Tyler D. Robinson
;  inputs:
;    Ahr  - hi-res planetary albedo spectrum
;  lamhr  - wavelength grid for Ahr (um)
;  alpha  - phase angle (deg)
;    Phi  - phase function evaluated at alpha
;     Rp  - planetary radius (Rearth)
;   Teff  - stellar effective temperature (K)
;     Rs  - stellar radius (Rsun)
;      r  - orbital separation (au)
;      d  - distance to system (pc)
;    Nez  - number of exozodis
;  lammin - minimum wavelength (um)
;  lammax - maximum wavelength (um)
;     Res - spectral resolution (lambda/Dlambda)
;       X - size of photometric aperture (lambda/D)
;    diam - telescope diameter (m)
;    Tput - system throughput
;    cont - raw contrast
;     IWA - inner working angle (lambda/D)
;     OWA - outer working angle (lambda/D; unless /FIXOWA)
;      De - dark current (s**-1)
;  DNhpix - horizontal pixel spread of IFS spectrum
;      Re - read noise per pixel
;   Dtmax - maximum exposure time (hr)
;
;  outputs:
;    lam  - low-res wavelength grid (um)
;   dlam  - spectral bin widths (um)
;      A  - planetary albedo spectrum at low-res
;      q  - quantum efficiency
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
;    COMPUTE_LAM - set to compute lo-res wavelength grid, otherwise 
;                  the grid input as variable 'lam' is used
;   COMPUTE_DLAM - set to compute lo-res wavelength grid widths,  
;                  otherwise the grid input as variable 'dlam' is used
;            NIR - re-adjusts pixel size in NIR, as would occur if a 
;                  second instrument was designed to handle the NIR
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO CORONAGRAPH, Ahr, lamhr, alpha, Phi, Rp, Teff, Rs, r, d, Nez, lammin, lammax, $
                 Res, X, diam, Tput, cont, IWA, OWA, De, DNhpix, Re, Dtmax, $
                 lam, dlam, A, q, Cratio, cp, csp, cz, cez, cD, cR, $
                 ctot, DtSNR, COMPUTE_LAM = compute_lam, $
                 COMPUTE_DLAM = compute_dlam, NIR = nir

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; set key system parameters ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  fpa    = f_airy(X) ; fraction of planetary signal in Airy pattern
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; set astrophys parameters ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  MzV  = 23.0 ; zodiacal light surface brightness (mag/arcsec**2)
  MezV = 22.0 ; exozodiacal light surface brightness (mag/arcsec**2)
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;   set wavelength grid    ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  IF KEYWORD_SET(compute_lam) THEN BEGIN
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
  ENDIF
  Nlam = N_ELEMENTS(lam)
  IF KEYWORD_SET(compute_dlam) THEN BEGIN
    dlam = FLTARR(Nlam) ;grid widths (um)
    FOR j=1,Nlam-2 DO BEGIN
      dlam[j] = 0.5*(lam[j+1]+lam[j]) - 0.5*(lam[j-1]+lam[j])
    ENDFOR
    ;widths at edges are same as neighbor
    dlam[0] = dlam[1]
    dlam[Nlam-1] = dlam[Nlam-2]
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
  ;;;  angular size of lenslet ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  theta = lammin/1.d6/diam/2.*(180/!DPI*3600.) ;diameter (assumes sampled at ~lambda/2D) (arcsec)
  IF KEYWORD_SET(nir) THEN BEGIN ;adjust NIR parameters, if keyword set
    theta = DBLARR(Nlam)
    iVIS  = WHERE(lam LE 1.0)
    iNIR  = WHERE(lam GT 1.0)
    theta[iVIS] = lammin/1.d6/diam/2.*(180/!DPI*3600.)
    theta[iNIR] = 1.0/1.d6/diam/2.*(180/!DPI*3600.)
  ENDIF
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;      set throughput      ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  T    = DBLARR(Nlam)
  T[*] = Tput

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;    check IWA and OWA     ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  sep  = r/d*SIN(alpha*!DPI/180.)/(lam*1.e-6/diam*180/!DPI*3600) ;separation in lambda/D
  C    = DBLARR(Nlam)
  C[*] = cont
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
  ;;; degrade albedo spectrum  ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  IF KEYWORD_SET(compute_lam)  THEN A = DEGRADE_SPEC(Ahr,lamhr,lam,DLAM=dlam)
  IF ~KEYWORD_SET(compute_lam) THEN A = Ahr
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;      compute fluxes      ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  Fs = Fstar(lam, Teff, Rs, r, /AU) ;stellar flux on planet (W/m**2/um)
  Fp = Fplan(A, Phi, Fs, Rp, d)     ;planet flux at telescope (W/m**2/um)
  Cratio = FpFs(A, Phi, Rp, r)      ;planet-to-star flux ratio
  
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
