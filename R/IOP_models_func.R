#### Detritus ####

#' @noRd
IOP_d_B22_C2 <- function(ISM,
                         Chl,
                         A = NULL,
                         G = NULL,
                         A_bd = NULL,
                         S_bd = NULL,
                         C_bd = NULL,
                         A_md = NULL,
                         S_md = NULL,
                         C_md = NULL,
                         qt_bd = 0.5,
                         qt_md = 0.5,
                         wavelen = wavelen_IOP) {

  use_qt <- FALSE

  if (is.null(A_bd) | is.null(S_bd) | is.null(C_bd) |
      is.null(A_md) | is.null(S_md) | is.null(C_bd)) {

    use_qt = TRUE

  } else {

    coef_bd <- list(A_bd = A_bd,
                    S_bd = S_bd,
                    C_bd = C_bd)
    coef_md <- list(A_md = A_md,
                    S_md = S_md,
                    C_md = C_md)

  }

  if (use_qt) {

    if (is.null(qt_bd))
      qt_bd <- round(runif(1), 2)
    if (is.null(qt_md))
      qt_md <- round(runif(1), 2)

    coef_bd <- coef_exp_ads$abds %>% melt.data.table(id.vars = 1) %>%
      .[, .(value = with(.SD, approx(quantile, value, qt_bd, rule = 2)$y)), by = .(variable)]
    coef_bd <-
      coef_bd$value %>% setNames(., paste0(coef_bd$variable, "_bd")) %>% as.list()

    coef_md <- coef_exp_ads$amds %>% melt.data.table(id.vars = 1) %>%
      .[, .(value = with(.SD, approx(quantile, value, qt_md, rule = 2)$y)), by = .(variable)]
    coef_md <-
      coef_md$value %>% setNames(., paste0(coef_md$variable, "_md")) %>% as.list()

  }

  # exponential functions for absorption
  abds <- with(coef_bd, A_bd * exp(-S_bd * (wavelen - 550)) + C_bd)
  amds <- with(coef_md, A_md * exp(-S_md * (wavelen - 550)) + C_md)
  ad <- Chl * abds + ISM * amds
  ad550 <- approx(wavelen, ad, xout = 550)$y

  A_mean <- -1.3390
  A_sd <- 0.0618
  G_mean <-  0.3835
  G_sd <- 0.1277

  if (is.null(A))
    A <- rnorm(1, mean = A_mean, sd = A_sd)
  if (is.na(A))
    A <- rnorm(1, mean = A_mean, sd = A_sd)
  if (is.null(G))
    G <- rnorm_bound(1,
                     mean = G_mean,
                     sd = G_sd,
                     lo = 0,
                     up = NA)
  if (is.na(G))
    G <- rnorm_bound(1,
                     mean = G_mean,
                     sd = G_sd,
                     lo = 0,
                     up = NA)

  if (A < 0) {
    Albedo_550 <- 1 - 10 ^ A
  } else {
    Albedo_550 <- A
  }

  # powerlaw function for attenuation
  cd550 <- ad550 / (1 - Albedo_550)

  # Albedo reference
  # Stramski, Dariusz, Marcel Babin, and Slawomir B. Woźniak. “Variations in the
  # Optical Properties of Terrigenous Mineral-Rich Particulate Matter Suspended
  # in Seawater.” Limnology and Oceanography 52, no. 6 (November 2007): 2418–33.
  # https://doi.org/10.4319/lo.2007.52.6.2418.

  cd <- cd550 * (wavelen / 550) ^ -G
  bd <- cd - ad

  if (any(bd < 0)) {
    if (use_qt) {
      cat("Negative bd! Rational random values will be used!\n")
    }

    while (any(bd < 0)) {
      qt_bd <- round(runif(1), 2)
      qt_md <- round(runif(1), 2)

      coef_bd <- coef_exp_ads$abds %>% melt.data.table(id.vars = 1) %>%
        .[, .(value = with(.SD, approx(quantile, value, qt_bd)$y)), by = .(variable)]
      coef_bd <-
        coef_bd$value %>% setNames(., paste0(coef_bd$variable, "_bd")) %>% as.list()

      coef_md <- coef_exp_ads$amds %>% melt.data.table(id.vars = 1) %>%
        .[, .(value = with(.SD, approx(quantile, value, qt_md)$y)), by = .(variable)]
      coef_md <-
        coef_md$value %>% setNames(., paste0(coef_md$variable, "_md")) %>% as.list()

      abds <-
        with(coef_bd, A_bd * exp(-S_bd * (wavelen - 550)) + C_bd)
      amds <-
        with(coef_md, A_md * exp(-S_md * (wavelen - 550)) + C_md)
      ad <- Chl * abds + ISM * amds
      ad550 <- approx(wavelen, ad, xout = 550)$y

      A <- rnorm(1, mean = A_mean, sd = A_sd * 2)
      G <- rnorm_bound(
        1,
        mean = G_mean,
        sd = G_sd * 2,
        lo = 0,
        up = NA
      )

      Albedo_550 <- 1 - 10 ^ A
      cd550 <- ad550 / (1 - Albedo_550)

      cd <- cd550 * (wavelen / 550) ^ -G
      bd <- cd - ad

    }

  }

  r <- list(
    wavelen = wavelen,
    ad = ad,
    bd = bd,
    cd = cd
  )

  parm <- list(
    Chl = Chl,
    ISM = ISM,
    A_cd = A,
    G_cd = G,
    qt_bd = qt_bd,
    qt_md = qt_md,
    Albedo_550 = Albedo_550,
    A_bd = coef_bd$A_bd,
    S_bd = coef_bd$S_bd,
    C_bd = coef_bd$C_bd,
    A_md = coef_md$A_md,
    S_md = coef_md$S_md,
    C_md = coef_md$C_md
  )

  attr(r, "parm") <- parm

  return(r)

}

#### Phytoplankton ####

#' @noRd
IOP_ph_B22_C2 <- function(Chl,
                          frac_phyto = rand_frac(),
                          wavelen = wavelen_IOP,
                          aphs = NULL,
                          # if they are NULL, call `generate_aphs`
                          bphs = NULL,
                          aphs_func = calc_aph676,
                          ...) {
  if (is.null(aphs) | is.null(bphs)) {
    siop_ph <- generate_aphs(Chl, aphs_fun = aphs_func, ...)
    aphs <- siop_ph$aphs
    bphs <- siop_ph$bphs
    attr_aph676 <- attr(siop_ph, "attr_aph676")
  } else {
    attr_aph676 <- NULL
  }

  # aphs <- phytodive_iop$aphs
  # aphs <- as.matrix(aphs[aphs$wv %in% wavelen, -1])
  aphs <- as.data.table(aphs) %>%
    melt.data.table(id.vars = 1) %>%
    .[, {
      approx(.SD$wv, .SD$value, xout = wavelen, rule = 2) %>% as.data.table()
    }, by = .(variable)] %>%
    dcast.data.table(x ~ variable, value.var = "y") %>% .[, -1] %>% as.matrix()

  # bphs <- phytodive_iop$bphs
  # bphs <- as.matrix(bphs[bphs$wv %in% wavelen, -1])
  bphs <- as.data.table(bphs) %>%
    melt.data.table(id.vars = 1) %>%
    .[, {
      approx(.SD$wv, .SD$value, xout = wavelen, rule = 2) %>% as.data.table()
    }, by = .(variable)] %>%
    dcast.data.table(x ~ variable, value.var = "y") %>% .[, -1] %>% as.matrix()


  frac_mat <- vec_to_mat(frac_phyto, n = nrow(aphs), by = 1)

  aphs_sum <- rowSums(frac_mat * aphs)
  bphs_sum <- rowSums(frac_mat * bphs)
  cphs_sum <- aphs_sum + bphs_sum

  aph <- Chl * aphs_sum
  bph <- Chl * bphs_sum
  cph <- Chl * cphs_sum

  r <- list(
    wavelen = wavelen,
    aph = aph,
    bph = bph,
    cph = cph
  )

  parm <- list(
    Chl = Chl,
    aphs = aphs,
    bphs = bphs,
    cphs = aphs + bphs,
    frac_phyto = frac_phyto
  )

  attr(r, "parm") <- parm

  attr(r, "attr_aph676") <- attr_aph676

  return(r)

}

#' @noRd
calc_aph676 <- function(Chl, varyA = TRUE) {

  A    = 0.0239
  A_sd = 0.1
  A_up = 0.0501
  A_lo = 0.0112
  E    = 0.8938

  if (varyA) {
    A = 10 ^ rnorm_bound(
      1,
      mean = log10(A),
      sd = A_sd,
      up = log10(A_up),
      lo = log10(A_lo)
    )
  }

  # split them in two Chl ranges
  r1 <- A * Chl ^ 1 # linear Chl <  1
  r2 <- A * Chl ^ E # powerl Chl >= 1
  r <- r1
  w <- which(Chl >= 1)
  r[w] <- r2[w]

  attr(r, "A_aphs676") <- A
  attr(r, "E_aphs676") <- E
  attr(r, "varyA") <- varyA

  class(r) <- "Hereon_aphs_func"

  r

}

#' @noRd
generate_aphs <- function(Chl,
                          albedo_ph676 = c(0.8952, 0.8952, 0.8800, 0.9000, 0.9000, 0.9562, 0.8485),
                          gamma_ph676  = c(NA,   0.00, 0.00, 0.00, 0.00, NA,   NA),
                          vary_cph = FALSE,
                          aphs_fun = calc_aph676,
                          ...) {

  # Some albedo reference
  # Figure 8 (w532) of Nardelli, Schuyler C., and Michael S. Twardowski. “Assessing the
  # Link between Chlorophyll Concentration and Absorption Line Height at 676 Nm
  # over a Broad Range of Water Types.” Optics Express 24, no. 22 (October 31,
  # 2016): A1374. https://doi.org/10.1364/OE.24.0A1374.


  # 1. [Chl] -> aphs676
  # 2. [albedo] -> cphs676
  # 3. [Normalized aphs and cphs] -> bphs

  if (any(albedo_ph676 > 1) | any(albedo_ph676 < 0)) {
    stop("albedo of phytoplankton at 676 nm should between 0 and 1")
  }
  names(albedo_ph676) <- names(phytodive_iop$name_phyto)

  if (vary_cph) {
    albedo_ph676_new <- albedo_ph676
    for (i in 1:length(albedo_ph676)) {
      albedo_ph676_new[i] <-
        rnorm_bound(1, mean = albedo_ph676[i], sd = 0.03, lo = NA, up = 0.99)
    }
    albedo_ph676 <- albedo_ph676_new
  }

  if (any(gamma_ph676 < -1, na.rm = TRUE)) {
    stop("Gamma of phytoplankton attenuation normalized at 676 nm should > -1")
  }
  names(gamma_ph676) <- names(phytodive_iop$name_phyto)

  if (vary_cph) {
    gamma_ph676_new <- gamma_ph676
    for (i in 1:length(gamma_ph676)) {
      if (is.na(gamma_ph676_new[i])) {
        next # skip for gamma is NA
      } else {
        gamma_ph676_new[i] <-
          rnorm_bound(1, mean = gamma_ph676[i], sd = 0.01, lo = 0.0, up = 0.8)
      }
    }
    gamma_ph676 <- albedo_ph676_new
  }

  # determine aph676 from the input function

  aph676 <- aphs_fun(Chl, ...)

  if ("Hereon_aphs_func" %in% class(aph676)) {
    attr_aph676 <- attributes(aph676)
  } else {
    attr_aph676 <- list(aphs_fun)
  }

  aphs_676_new <- aph676 / Chl

  aphs_mat <- as.matrix(phytodive_iop$aphs[, -1])
  aphs_676 <-
    unlist(phytodive_iop$aphs[phytodive_iop$aphs$wv == 676, -1])
  aphs_mat_norm <-
    aphs_mat / vec_to_mat(aphs_676, n = nrow(aphs_mat))
  aphs_mat_new  <-
    aphs_mat_norm * aphs_676_new # using new aphs from Chl

  cphs_mat <- as.matrix(phytodive_iop$cphs[, -1])
  cphs_676 <-
    unlist(phytodive_iop$cphs[phytodive_iop$cphs$wv == 676, -1])
  cphs_mat_norm <-
    cphs_mat / vec_to_mat(cphs_676, n = nrow(cphs_mat), by = 1) # normalized at cph676 - shape
  # Find group to be changed for the shape
  PG_TBC <- which(!is.na(gamma_ph676))
  for (i in PG_TBC) {
    cphs_mat_norm[, i] <- (676 / phytodive_iop$cphs$wv) ^ gamma_ph676[i]
  }
  coef_lin <-
    aphs_676_new / (1 - albedo_ph676) # the coefficients in Eq 16 - magnitude
  cphs_mat_new <-
    cphs_mat_norm * vec_to_mat(coef_lin, n = nrow(cphs_mat), by = 1)
  bphs_mat_new <- cphs_mat_new - aphs_mat_new

  aphs <- data.table(wv = phytodive_iop$aphs$wv, aphs_mat_new)
  bphs <- data.table(wv = phytodive_iop$bphs$wv, bphs_mat_new)
  cphs <- data.table(wv = phytodive_iop$cphs$wv, cphs_mat_new)

  r <- list(aphs = aphs,
            bphs = bphs,
            cphs = cphs)

  parm <- list(
    aph676 = aph676,
    albedo_ph676 = albedo_ph676,
    gamma_ph676 = gamma_ph676,
    vary_cph = vary_cph
  )

  attr(r, "parm") <- parm

  attr(r, "attr_aph676") <- attr_aph676

  return(r)

}


#' @noRd
rand_frac <-
  function(case = 2,
           Chl = 0.5,
           r_up = NULL,
           r_lo = NULL) {

    # Separate Case 1 and Case 2
    nr = nrow(frac_phyto_lib)
    nc = ncol(frac_phyto_lib)

    if (!is.null(r_up) & !is.null(r_lo)) {
      if (length(r_up) != nc | length(r_lo) != nc) {
        stop("The defined lower and upper bounds should be same to the number of phytoplankton!")
      }

      if (any(r_lo > r_up)) {
        stop("The lower limit should be lower than the upper one!")
      }

      prob_col <- rep(1, nc)
      prob_row <- rep(1, nr)

    } else {
      if (case == 1) {
        r_up <- c(0.5, 0.5, 0.3, 0.5, 0.5, 0.1, 1.0)
        r_lo <- rep(0, nc)

        if (Chl < 0.05) {
          r <- c(rep(0, nc - 1), 1)
          names(r) <- names(phytodive_iop$name_phyto)
          return(r)
        }
        if (Chl >= 0.05 & Chl < 0.2) {
          r_lo <- c(rep(0, nc - 1), 0.5)
        }
        if (Chl >= 0.2 & Chl < 1) {
          r_lo <- c(rep(0, nc - 1), 0.2)
        }
        if (Chl > 1) {
          r_lo <- c(0.2, rep(0, nc - 2), 0.2)
        }

        prob_col <- rep(1, nc)
        prob_row <- rep(1, nr)

      } else if (case == 2) {
        r_up <- rep(1, nc)
        r_lo <- rep(0, nc)

        prob_col <- rep(1, nc)
        prob_row <- rep(1, nr)

      }

    }

    r <-
      frac_phyto_lib[sample(1:nr, 1, prob = prob_row),
                     sample(nc, nc, prob = prob_col)]
    names(r) <- names(phytodive_iop$name_phyto)

    cond = TRUE
    iter = 1
    from_lib = TRUE
    while (cond) {
      if (from_lib) {
        r <-
          frac_phyto_lib[sample(1:nr, 1, prob = prob_row),
                         sample(nc, nc, prob = prob_col)]
      } else {
        r <-
          apply(cbind(r_lo, r_up), 1, function(x)
            runif(1, x[1], x[2])) %>% {
              round(. / sum(.), 4)
            }
      }
      names(r) <- names(phytodive_iop$name_phyto)
      iter <- iter + 1
      if (all((r - r_up) <= 0) & all((r - r_lo) >= 0)) {
        cond = FALSE
      }
      if (iter == 50 & from_lib == FALSE) {
        cond = FALSE
      }
      if (iter == 50) {
        from_lib = FALSE
        iter <- 1
      }
    }

    r

  }


#### CDOM ####

#' @noRd
IOP_cdom_B22_C2 <- function(
    ag440 = NULL,
    S = NULL,
    wavelen = wavelen_IOP
) {

  if(is.null(ag440)) ag440 = 10^rnorm(1, mean = -0.9425989, sd = 0.4291600/2)
  if(is.na(ag440))   ag440 = 10^rnorm(1, mean = -0.9425989, sd = 0.4291600/2)
  if(is.null(S)) S = round(rnorm(1, mean = 0.017400895, sd = 0.001395171/2), 4)
  if(is.na(S))   S = round(rnorm(1, mean = 0.017400895, sd = 0.001395171/2), 4)

  r <-
    list(
      wavelen = wavelen,
      acdom = ag440 * exp(-S * (wavelen - 440))
    )

  parm <- list(ag440 = ag440, S_cdom = S)

  attr(r, "parm") <- parm

  return(r)

}

#' @noRd
IOP_cdom_B22_C2_lib <- function(
    ag440 = NULL,
    wavelen = wavelen_IOP,
    ag_seed = 1234
) {

  set.seed(ag_seed)

  if(min(wavelen) < 354) {
    ag_norm <-
      ag_lib[!(SampleID %in% c("I080919"))] %>%
      .[SampleID == sample(SampleID, 1), .(wavelen, value)]
  } else {
    ag_norm <- ag_lib[SampleID == sample(SampleID, 1), .(wavelen, value)]
  }

  # TBD
  # if wavelen < 250, the CDOM returns NA since the lib does not include that range

  ag_norm <- approx(ag_norm$wavelen, ag_norm$value, xout = wavelen)$y

  w_nan <- which(is.na(ag_norm) & wavelen >= 600)

  ag_norm[w_nan] <- 0

  r <-
    list(
      wavelen = wavelen,
      acdom = ag440 * ag_norm
    )

  parm <- list(ag440 = ag440)

  attr(r, "parm") <- parm

  return(r)

}

#### Alternative bph functions ####

#' @noRd
specific_bph_GM83 <-
  function(Chl,
           wavelen = wavelen_IOP,
           b0 = 0.3,
           n = 0.62,
           m = 1,
           wavelen0 = 550) {
    b0 * Chl ^ n * (wavelen0 / wavelen) ^ m / Chl
  }

#' @noRd
specific_bph_LM88 <-
  function(Chl,
           wavelen = wavelen_IOP,
           b0 = 0.416,
           n = 0.766,
           wavelen0 = 550) {
    m = ifelse(Chl >= 2, 0, 0.5 * (log10(Chl) - 0.3))
    b0 * Chl ^ n * (wavelen / wavelen0) ^ m / Chl
  }

#' @noRd
specific_bph_G99 <-
  function(Chl,
           wavelen = wavelen_IOP,
           b0 = 0.5,
           n = 0.62,
           m = -0.00113,
           i = 1.62517,
           wavelen0 = 550) {
    b0 * ((m * wavelen + i) / (m * wavelen0 + i)) * Chl ^ n / Chl
  }

#### Alternative aph functions ####

#' @noRd
specific_aph_B04 <- function(
    Chl = 1,
    a = 0.0654,
    b = 0.7280,
    # frac_phyto = rand_frac(),
    frac_phyto = rep(1/7, 7),
    wavelen = wavelen_IOP
) {

  aphs <- phytodive_iop$phs_n440$aphs
  aphs <- as.matrix(aphs[aphs$wv %in% wavelen, -1])

  frac_mat <- vec_to_mat(frac_phyto, n = nrow(aphs), by = 1)

  aphs_440 <- a * Chl ^ b / Chl

  aphs_sum <- rowSums(frac_mat * aphs) * aphs_440

  aphs_sum

}

#' @noRd
specific_aph_B98 <- function(
    Chl = 1,
    wavelen = wavelen_IOP
) {

  AE <- B98_AE_midUVabs

  AE_approx <- list(
    A = approx(AE$wavelen, AE$A, wavelen, rule = 2) %>% as.data.frame(),
    B = approx(AE$wavelen, AE$E, wavelen, rule = 2) %>% as.data.frame()
  )
  AE_approx <-
    data.table::merge.data.table(AE_approx$A, AE_approx$B, by = "x") %>%
    setNames(., c("wavelen", "A", "E"))

  aphs <- AE_approx$A * Chl ^ AE_approx$E / Chl

  aphs

}








