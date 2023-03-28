#' WOPP
#'
#' @param Temp Temperature \[degC\]
#' @param Sal Salinity \[PSU\]
#' @param wavelen Wavelength \[nm\]
#'
#' @note
#'
#' The WOPP code is wrapper from Schoenfeld's codes
#'
#' This code is copied from https://github.com/bishun945/RrsTrans/blob/master/R/WOPP.R
#'
#' @return A list of IOPs of pure water
#' @export
#'
#' @references
#'
#' Röttgers, R, R Doerffer, D McKee, and W Schönfeld. “The Water Optical
#' Properties Processor (WOPP): Pure Water Spectral Absorption, Scattering and
#' Real Part of Refractive Index Model.” Technical Report No WOPP-ATBD/WRD6,
#' 2016. https://calvalportal.ceos.org/tools.
#'
#' @examples
#' WOPP(Temp = 25, Sal = 10)
WOPP <- function(Temp = 25, Sal = 10, wavelen = seq(350, 900, 5)) {

  # WOPP_path <- "/WOPP/"

  aw = min(wavelen)
  ew = max(wavelen)

  S = Sal
  Tc = Temp

  #functions to compute the IOPs
  read_refri_std <- function(){
    #read the standard spectrum of refraction index
    #from file 'fn'
    #from wavelength aw to wavelength ew
    #return vector of wavelength, vector of refractive index

    # if(FALSE) {
    #   fn = file.path("./WOPP", "computed_refri_T27_S0_180_4000nm.dat")
    #   WOPP_computed_refri_T27_S0_180_4000nm <-
    #     read.table(fn, comment.char = "%", col.names = c("wl", "refri"))
    #   save(WOPP_computed_refri_T27_S0_180_4000nm,
    #        file = "./data/WOPP_computed_refri_T27_S0_180_4000nm.rda")
    # }

    d <- WOPP_computed_refri_T27_S0_180_4000nm
    idx <- which(d$wl >= aw & d$wl <= ew)
    refri_err <- numeric(nrow(d))
    refri_err[d$wl < 700] <- 0.00003
    refri_err[d$wl >= 700] <- 0.0005

    refri_err <- refri_err[idx]
    wl <- d$wl[idx]
    refri <- d$refri[idx]

    return(list(wl = wl, refri = refri, refri_err = refri_err))

  }

  absorption <- function(fn = NULL, version = 3) {

    # absorption 300 - 4000 nm
    # Temp coeff 300 - 4000 nm
    # Sal coeff 300 - 4000 nm

    if(!(version %in% c(1, 2, 3))) stop("Version should be 1, 2, or 3")

    # if(FALSE) {
    #
    #   version = 3
    #   fn = file.path(WOPP_path, sprintf("purewater_abs_coefficients_v%s.dat", version))
    #
    #   WOPP_purewater_abs_coefficients_v3 <-
    #     read.table(
    #       fn,
    #       comment.char = "%",
    #       col.names = c("wl", "a0", "T_coeff", "S_coeff", "a0_err", "T_err", "S_err")
    #     )
    #
    #   save(WOPP_purewater_abs_coefficients_v3, file = "./data/WOPP_purewater_abs_coefficients_v3.rda")
    #
    # }

    d <- WOPP_purewater_abs_coefficients_v3

    idx <- which(d$wl >= aw & d$wl <= ew)

    d <- d[idx,]

    wl <- d$wl
    abso_ref <- d$a0
    S_coeff <- d$S_coeff
    T_coeff <- d$T_coeff
    a0_err <- d$a0_err
    S_err <- d$S_err
    T_err <- d$T_err

    T_ref= 20.0
    S_ref= 0.0

    abso = abso_ref + (Tc - T_ref) * T_coeff  + (S - S_ref) * S_coeff
    a_err = a0_err + (Tc - T_ref) * T_err  + (S - S_ref) * S_err

    return (list(wl = wl,  abso = abso,  a_err = a_err))

  }

  scattering <- function(wl, theta = 90, std_RI) {

    # Xiaodong Zhang, Lianbo Hu, and Ming-Xia He (2009), Scattering by pure
    # seawater: Effect of salinity, Optics Express, Vol. 17, No. 7, 5698-5710
    #
    # wl (nm): wavelength
    # Tc: temperauter in degree Celsius, must be a scalar
    # S: salinity, must be scalar
    # delta: depolarization ratio, if not provided, default = 0.039 will be
    # used.
    # betasw: volume scattering at angles defined by theta. Its size is [x y],
    # where x is the number of angles (x = length(theta)) and y is the number
    # of wavelengths in wl (y = length(wl))
    # beta90sw: volume scattering at 90 degree. Its size is [1 y]
    # bw: total scattering coefficient. Its size is [1 y]
    # for backscattering coefficients, divide total scattering by 2
    #
    # Xiaodong Zhang, March 10, 2009
    delta = 0.039
    # values of the constants
    Na = 6.0221417930e23     #  Avogadro's constant
    Kbz = 1.3806503e-23      #  Boltzmann constant
    Tk = Tc + 273.15           #  Absolute tempearture
    M0 = 18e-3               #  Molecular weigth of water in kg/mol

    rad = theta * pi / 180  # angle in radian as a colum variable

    # nsw: absolute refractive index of seawater
    # dnds: partial derivative of seawater refractive index w.r.t. salinity
    ##nsw, dnds = RInw(wl,Tc,S,  alt_RI)
    res_refri = refractive_index(std_RI)
    nsw <- res_refri$nsw
    dnds <- res_refri$dnswds

    #    print 'nsw',  nsw, wl
    # isothermal compressibility is from Lepple & Millero (1971,Deep
    # Sea-Research), pages 10-11
    # The error ~ +/-0.004e-6 bar^-1
    IsoComp = BetaT()

    # density of water and seawater,unit is Kg/m^3, from UNESCO,38,1981
    density_sw = rhou_sw()
    ##print 'density_sw', density_sw
    # water activity data of seawater is from Millero and Leung (1976,American
    # Journal of Science,276,1035-1077). Table 19 was reproduced using
    # Eq.(14,22,23,88,107) then were fitted to polynominal equation.
    # dlnawds is partial derivative of natural logarithm of water activity
    # w.r.t.salinity
    dlnawds = dlnasw_ds()
    ##print 'dlnawds', dlnawds
    # density derivative of refractive index from PMH model
    DFRI = PMH(nsw)   ## PMH model
    #print 'DFRI', DFRI.shape, wl.shape
    # volume scattering at 90 degree due to the density fluctuation
    beta_df = pi * pi / 2 * ((wl * 1e-9) ^ (-4)) *
      Kbz * Tk * IsoComp * DFRI ^ 2 * (6 + 6 * delta) / (6 - 7 * delta)
    ##print 'beta_df',  beta_df
    # volume scattering at 90 degree due to the concentration fluctuation
    flu_con = S * M0 * dnds ^ 2 / density_sw / (-dlnawds) / Na
    ##print 'flu_con',  flu_con
    beta_cf = 2 * pi * pi * ((wl * 1e-9) ^ (-4)) * nsw ^ 2 *
      (flu_con) * (6 + 6 * delta) / (6 - 7 * delta)
    ##print 'beta_cf',  beta_cf
    # total volume scattering at 90 degree
    beta90sw = beta_df + beta_cf
    ##print 'beta90sw',  beta90sw
    bsw = 8 * pi / 3 * beta90sw * (2 + delta) / (1 + delta)
    ##print 'bsw', bsw
    betasw = matrix(nrow = length(rad), ncol = length(wl))
    rownames(betasw) <- theta
    colnames(betasw) <- wl
    for(i in 1:length(wl)) {
      betasw[,i] <- beta90sw[i]*(1+((cos(rad))^2)*(1-delta)/(1+delta))
    }
    betasw=t(betasw)
    err_betasw=betasw*0.02
    err_beta90sw=beta90sw*0.02
    err_bsw=bsw*0.02

    return (
      list(
        betasw = (betasw),
        beta90sw = beta90sw,
        bsw = bsw,
        err_betasw = (err_betasw),
        err_beta90sw = (err_beta90sw),
        err_bsw = err_bsw,
        nsw = nsw
      )
    )

  }

  refractive_index <- function(std_RI) {

    # std_RI from `read_refri_std` refri_err

    # real part of R_efractive In_dex of Sea_w_ater
    # Refractive Index of air is from Ciddor (1996,Applied Optics)
    n_air = 1.0+(5792105.0/(238.0185-1/(wl/1e3)^2)+167917.0/(57.362-1/(wl/1e3)^2))/1e8
    #    print 'n_air',  n_air
    # refractive index of seawater is from Quan and Fry (1994, Applied Optics)
    n0 = 1.31405
    n1 = 1.779e-4
    n2 = -1.05e-6
    n3 = 1.6e-8
    n4 = -2.02e-6
    n5 = 15.868
    n6 = 0.01155
    n7 = -0.00423
    n8 = -4382.
    n9 = 1.1455e6
    #    print Tc,  S,  wl
    nsw = n0+(n1+n2*Tc+n3*Tc^2)*S+n4*Tc^2+(n5+n6*S+n7*Tc)/wl+n8/wl^2+n9/wl^3  # pure seawater
    nsw = nsw*n_air
    #    print 'nsw',  nsw
    dnswds = (n1+n2*Tc+n3*Tc^2+n6/wl)*n_air

    idx <- which(wl == 800)
    offs <- nsw[idx] - std_RI[idx]
    nsw[wl >= 800] <- std_RI[wl >= 800] + offs

    return (list(nsw = nsw, dnswds = dnswds))

  }

  BetaT <- function() {

    # pure water secant bulk Millero (1980, Deep-sea Research)
    kw = 19652.21 + 148.4206 * Tc - 2.327105 * Tc ^ 2 +
      1.360477e-2 * Tc ^ 3 - 5.155288e-5 * Tc ^ 4
    Btw_cal = 1 / kw
    # isothermal compressibility from Kell sound measurement in pure water
    # Btw = (50.88630+0.717582*Tc+0.7819867e-3*Tc^2+31.62214e-6*Tc^3-0.1323594e-6*Tc^4+0.634575e-9*Tc^5)/(1+21.65928e-3*Tc)*1e-6
    # seawater secant bulk
    a0 = 54.6746 - 0.603459 * Tc + 1.09987e-2 * Tc ^ 2 - 6.167e-5 * Tc ^ 3
    b0 = 7.944e-2 + 1.6483e-2 * Tc - 5.3009e-4 * Tc ^ 2
    Ks = kw + a0 * S + b0 * S ^ 1.5
    # calculate seawater isothermal compressibility from the secant bulk
    IsoComp = 1 / Ks * 1e-5  # unit is pa
    return (IsoComp)
  }

  rhou_sw <- function() {

    # density of water and seawater,unit is Kg/m^3, from UNESCO,38,1981
    a0 = 8.24493e-1
    a1 = -4.0899e-3
    a2 = 7.6438e-5
    a3 = -8.2467e-7
    a4 = 5.3875e-9
    a5 = -5.72466e-3
    a6 = 1.0227e-4
    a7 = -1.6546e-6
    a8 = 4.8314e-4
    b0 = 999.842594
    b1 = 6.793952e-2
    b2 = -9.09529e-3
    b3 = 1.001685e-4
    b4 = -1.120083e-6
    b5 = 6.536332e-9
    # density for pure water
    density_w = b0+b1*Tc+b2*Tc^2+b3*Tc^3+b4*Tc^4+b5*Tc^5
    # density for pure seawater
    density_sw = density_w +((a0+a1*Tc+a2*Tc^2+a3*Tc^3+a4*Tc^4)*S+(a5+a6*Tc+a7*Tc^2)*S^1.5+a8*S^2)
    return (density_sw)

  }

  dlnasw_ds <- function() {

    # water activity data of seawater is from Millero and Leung (1976,American
    # Journal of Science,276,1035-1077). Table 19 was reproduced using
    # Eqs.(14,22,23,88,107) then were fitted to polynominal equation.
    # dlnawds is partial derivative of natural logarithm of water activity
    # w.r.t.salinity
    # lnaw = (-1.64555e-6-1.34779e-7*Tc+1.85392e-9*Tc^2-1.40702e-11*Tc^3)+......
    #            (-5.58651e-4+2.40452e-7*Tc-3.12165e-9*Tc^2+2.40808e-11*Tc^3)*S+......
    #            (1.79613e-5-9.9422e-8*Tc+2.08919e-9*Tc^2-1.39872e-11*Tc^3)*S^1.5+......
    #            (-2.31065e-6-1.37674e-9*Tc-1.93316e-11*Tc^2)*S^2

    dlnawds = (-5.58651e-4+2.40452e-7*Tc-3.12165e-9*Tc^2+2.40808e-11*Tc^3)+1.5*
      (1.79613e-5-9.9422e-8*Tc+2.08919e-9*Tc^2-1.39872e-11*Tc^3)*S^0.5+
      2*(-2.31065e-6-1.37674e-9*Tc-1.93316e-11*Tc^2)*S

    return(dlnawds)

  }

  # density derivative of refractive index from PMH model
  PMH <- function(n_wat) {

    n_wat2 = n_wat ^ 2
    ##print 'n_wat2', n_wat2
    n_density_derivative = (n_wat2 - 1) *
      (1 + 2.0 / 3.0 * (n_wat2 + 2.0) *
         (n_wat / 3.0 - 1.0 / 3.0 / n_wat) ^ 2)
    ##print 'n_density_derivative',  n_density_derivative
    return (n_density_derivative)

  }


  ############################################################################## MAIN

  res_refri <- read_refri_std()

  a = absorption()

  wl = a$wl

  b = scattering(wl, std_RI = res_refri$refri)

  # output
  #   - wavelen
  #   - absorption
  #   - scattering

  # do approx?
  water_abs = a$abso
  water_sca = b$bsw
  nw = b$nsw

  water_abs_approx <- approx(wl, water_abs, wavelen)$y
  water_sca_approx <- approx(wl, water_sca, wavelen)$y
  nw_approx <- approx(wl, nw, wavelen)$y

  return(list(wavelen = wavelen, a = water_abs_approx, b = water_sca_approx,
              Temp = Temp, Sal = Sal, nw = nw_approx))

}


