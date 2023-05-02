#' Four-term and Two-term IOP models by Bi et al. (202x)
#'
#' @param Chl Chlorophyll a concentration, as a proxy of biomass concentration
#'   of phytoplankton, in the unit \[mg/m^3\]
#' @param ag440 CDOM absorption coefficient at 440 nm, as a proxy of
#'   concentration of CDOM (or called gelbstoff), in the unit \[m^-1\]
#' @param ISM Inorganic suspended matter concentration in the unit \[g/m3\]
#' @param Temp Water temperature in the unit of \[degC\]
#' @param Sal Water salinity in the unit of \[PSU\]
#' @param qt_bd Quantile for biogenic detritus absorption, default as 0.5
#' @param qt_md Quantile for minerogenic detritus absorption, default as 0.5
#' @param frac_phyto Fraction of phytoplankton groups. Vector with a length of 7.
#' @param S_cdom Spectral slope for CDOM absorption coefficient. Not used if
#'   \code{lib_cdom} is \code{TRUE}
#' @param bbdtilde Backscattering probability of detritus
#' @param wavelen Wavelength in the unit of \[nm\]
#' @param A_d Scaled albedo parameters of detritus
#' @param G_d Power law exponent of spectral attenuation of detritus
#' @param lib_cdom Option to use the library of CDOM absorption. Default as
#'   \code{TRUE}
#' @param aw_version Pure water absorption version: 1, 2, and 3
#' @param ag_seed Random seed for ag model. Default as 1234. Fixed to reproduce
#'   the same outputs.
#' @param ... Parameters passed to \code{IOP_ph_B22_C2}
#'
#' @return A list including spectral of IOPs and key parameters
#'
#' @export
#' @rdname IOP_model
#'
#' @examples
#'
#' str(IOP_four_comp(Chl = 1, ag440 = 0.18, ISM = 1))
IOP_four_comp <- function(
    Chl = 1,
    ag440 = 0.18,
    ISM = 1,
    Temp = 20, Sal = 15,
    qt_bd = 0.5, qt_md = 0.5,
    frac_phyto = rand_frac(),
    S_cdom = 0.017422,
    bbdtilde = 0.0216,
    wavelen = wavelen_IOP,
    A_d = NULL, G_d = NULL,
    lib_cdom = TRUE,
    aw_version = 3,
    ag_seed = NULL,
    ...
) {

  # spectra
  list_d     <- IOP_d_B22_C2(ISM, Chl, qt_bd = qt_bd, qt_md = qt_md,
                             wavelen = wavelen,
                             A = A_d, G = G_d)
  list_ph    <- IOP_ph_B22_C2(Chl, frac_phyto, wavelen = wavelen, ...)

  if("attr_aph676" %in% names(attributes(list_ph))) {
    attr_aph676 <- attr(list_ph, "attr_aph676")
  } else {
    attr_aph676 <- NULL
  }

  if(lib_cdom) {
    list_cdom <- IOP_cdom_B22_C2_lib(ag440, wavelen = wavelen, ag_seed = ag_seed)
  } else {
    list_cdom <- IOP_cdom_B22_C2(ag440, S_cdom, wavelen = wavelen)
  }

  list_WOPP  <- WOPP(Temp, Sal, wavelen = wavelen, aw_version = aw_version)

  list_output <- list(
    wavelen = list_ph$wavelen,
    aph     = list_ph$aph,
    ad      = list_d$ad,
    acdom   = list_cdom$acdom,
    bph     = list_ph$bph,
    bd      = list_d$bd,
    aw      = list_WOPP$a,
    bw      = list_WOPP$b
  )

  list_output$ap  = with(list_output, aph + ad)
  list_output$agp = with(list_output, ap + acdom)
  list_output$bp  = with(list_output, bph + bd)

  list_output$at  = with(list_output, agp + aw)
  list_output$bt  = with(list_output, bp + bw)

  # list_output

  # bbtilde blend
  bbphtilde <- sum(frac_phyto * phytodive_iop$bbtilde_phyto)

  # bb spectra
  list_output$bbph = with(list_output, bph * bbphtilde)
  list_output$bbd  = with(list_output, bd * bbdtilde)
  list_output$bbp  = with(list_output, bbph + bbd)
  list_output$bbw  = with(list_output, bw / 2)
  list_output$bbt  = with(list_output, bbp + bbw)

  # Rrs by Lee et al 2011
  list_output$Rrs_L11  = with(list_output, IOP2AOP_Lee_11(bbw, bbp, at, bbt))
  list_output$Rrs_AM03 = with(list_output, IOP2AOP_Albert_Mobley_03(a = at, bb = bbt))

  r <- list(
    spec = list_output,
    parm = list(Chl   = Chl,
                ag440 = ag440,
                ISM   = ISM,
                Temp  = Temp,
                Sal   = Sal,
                qt_bd = attr(list_d, "parm")$qt_bd,
                qt_md = attr(list_d, "parm")$qt_md,
                A_d   = attr(list_d, "parm")$A_cd,
                G_d   = attr(list_d, "parm")$G_cd,
                Albedo550_d = attr(list_d, "parm")$Albedo_550,
                frac_phyto = frac_phyto,
                S_cdom    = S_cdom,
                bbdtilde  = bbdtilde,
                bbphtilde = bbphtilde,
                aw_version = aw_version,
                ag_seed = ag_seed)
  )

  attr(r, "attr_aph676") <- attr_aph676

  r

}


#' @rdname IOP_model
#' @param seed Seed number of the two-term model
#' @export
#' @examples
#' str(IOP_two_comp(Chl = 1))
IOP_two_comp <- function(
    Chl,
    frac_phyto = NULL,
    Temp = 20, Sal = 15,
    wavelen = wavelen_IOP,
    aw_version = 3,
    seed = NULL,
    ...
) {

  if(is.null(frac_phyto)) frac_phyto <- rand_frac(case = 1, Chl)

  set.seed(seed)

  ## water ##

  list_WOPP  <- WOPP(Temp, Sal, wavelen = wavelen, aw_version = aw_version)

  aw <- list_WOPP$a
  bw <- list_WOPP$b
  bbw <- list_WOPP$b * 0.5

  ## ph ##

  # Hereon blending model
  res_aphs <- generate_aphs(Chl, ...)
  list_ph <-
    1:length(res_aphs) %>%
    lapply(function(i) {
      tmp <-
        as.data.table(res_aphs[[i]]) %>%
        melt.data.table(id.vars = 1)
      tmp$comp = names(res_aphs)[i]
      tmp
    }) %>%
    data.table::rbindlist() %>%
    .[, {
      approx(.SD$wv, .SD$value, xout = wavelen, rule = 2) %>%
        as.data.table() %>%
        setNames(., c("wavelen", "value"))
    }, by = .(variable, comp)] %>%
    dcast.data.table(wavelen+comp~variable, value.var = "value") %>%
    .[,c("wavelen", "comp", names(phytodive_iop$name_phyto)), with = FALSE] %>%
    split(.$comp) %>%
    lapply(function(x) {
      r <- data.table(
        wavelen = x$wavelen,
        value = rowSums(as.matrix(x[,-c(1:2)]) *
                          vec_to_mat(frac_phyto, n = nrow(x), by = 1))
      )
      r
    })

  aph = list_ph$aph$value * Chl
  bph = list_ph$bph$value * Chl
  cph = list_ph$cph$value * Chl
  bbphtilde <- sum(frac_phyto * phytodive_iop$bbtilde_phyto)
  bbph = bph * bbphtilde
  aph440 = approx(wavelen, aph, 440)$y

  ## d ##

  p1 = round(0.1 + 0.5 * runif(1) * aph440 / (0.05 + aph440), 2) # 0.1~0.6
  ad440 = p1 * aph440
  Sad = round(runif(1, 0.007, 0.015), 4)
  ad = ad440 * exp(-Sad * (wavelen - 440))

  gamma_d <- -0.5 + (2.0 + 1.2 * runif(1)) / (1 + Chl ^ 0.5)
  p4 <- round(runif(1, 0.06, 0.6), 2)
  bd550 <- p4 * Chl ^ 0.766
  bd = bd550 * (550/wavelen) ^ gamma_d

  cd = ad + bd
  bbdtilde <- 0.0183
  # bbdtilde <- 0.01
  bbd <- bbdtilde * bd

  ## cdom ##

  p2 = round(0.3 + 5.7 * runif(1) * aph440 / (0.02 + aph440), 1) # 0.3~6.0
  ag440 = p2 * aph440
  Sag = round(runif(1, 0.01, 0.02), 4)
  ag = ag440 * exp(-Sag * (wavelen - 440))

  spec <- list(
    wavelen = wavelen,
    aph = aph,
    bph = bph,
    cph = cph,
    bbph = bbph,
    ad = ad,
    bd = bd,
    cd = cd,
    bbd = bbd,
    ap = aph + ad,
    bp = bph + bd,
    cp = cph + cd,
    bbp = bbph + bbd,
    acdom = ag,
    agp = aph + ad + ag,
    aw = aw,
    bw = bw,
    bbw = bbw,
    at = aph + ad + ag + aw,
    bt = bph + bd + bw,
    bbt = bbph + bbd + bbw
  )

  # Rrs by Lee et al 2011
  spec$Rrs_L11 = with(spec, IOP2AOP_Lee_11(bbw, bbp, at, bbt))
  spec$Rrs_AM03 = with(spec, IOP2AOP_Albert_Mobley_03(a = at, bb = bbt))

  parm <- list(
    Chl = Chl,
    aph440 = aph440,
    frac_phyto = frac_phyto,
    bbphtilde = bbphtilde,
    Temp = Temp,
    Sal = Sal,
    p1 = p1,
    p2 = p2,
    p4 = p4,
    acdom440 = ag440,
    Sag = Sag,
    ad440 = ad440,
    Sad = Sad,
    gamma_d = gamma_d,
    p4 = p4,
    bd550 = bd550,
    bbdtilde = bbdtilde,
    aw_version = aw_version,
    seed = seed
  )

  result <- list(
    spec = spec,
    parm = parm
  )

  return(result)

}










