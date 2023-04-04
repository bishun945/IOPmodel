
require(ggplot2)
require(magrittr)
require(data.table)

rm(list = ls())

# for WOPP
load("~/Documents/GitHub/HEShun/data/WOPP_computed_refri_T27_S0_180_4000nm.rda")

load("~/Documents/GitHub/HEShun/data/WOPP_purewater_abs_coefficients.rda")

WOPP_purewater_abs_coefficients <-
  WOPP_purewater_abs_coefficients %>%
  .[, {
    # tmp <- WOPP_purewater_abs_coefficients[version == 3]
    tmp <- copy(.SD)
    smooth_a0 <-
      loess(data = tmp, a0 ~ wl, span = 0.009) %>%
      predict()
    tmp[, a0_raw := a0]
    tmp[, a0 := smooth_a0]
  }, by = .(version)]

# smooth Phi_Sal
WOPP_purewater_abs_coefficients <-
  WOPP_purewater_abs_coefficients %>%
  as.data.table() %>%
  .[, {
    tmp <- copy(.SD)
    smoothed_S_coeff_pred <-
      loess(data = tmp[wl <= 506], S_coeff ~ wl, span = 0.3) %>%
      predict()
    wl_used <- tmp[wl <= 506, wl]
    tmp[, S_coeff_raw := S_coeff]
    offset <-
      tmp[wl == 500, S_coeff_raw] -
      smoothed_S_coeff_pred[wl_used == 500]
    tmp[between(wl, 380, 500), S_coeff := (smoothed_S_coeff_pred[between(wl_used, 380, 500)] + offset)]
    tmp[between(wl, NA, 394), S_coeff := 0]
    tmp$S_coeff <-
      loess(data = tmp, S_coeff ~ wl, span = 0.003) %>%
      predict()
    tmp
  }, by = .(version)]

# WOPP_purewater_abs_coefficients %>%
#     # .[version == 3 & between(wl, NA, 900)] %>%
#     .[version == 3 & between(wl, NA, 600)] %>%
#     ggplot(aes(x = wl)) +
#     geom_path(aes(y = S_coeff_raw, col = "S_coeff_raw")) +
#     geom_path(aes(y = S_coeff, col = "S_coeff")) +
#     # coord_cartesian(ylim = c(-5e-5, 5e-5))
#   coord_cartesian(ylim = c(-8e-5, 8e-5)) +
#   # geom_vline(xintercept = c(378, 380, 386, 388, 390, 500), linetype = 2)
#   geom_vline(xintercept = c(394, 500), linetype = 2)


# list_WOPP$a <- as.data.table(list_WOPP) %>%
#   loess(data = ., a ~ wavelen,  span = 0.04) %>% predict()

# for detrtius absorption table
load("~/Documents/GitHub/HEShun/data/coef_exp_ads.rda")

# import phytodive
load("~/Documents/GitHub/HEShun/data/phytodive_iop.rda")

phytodive_iop <-
  phytodive_iop[c("aphs",
                  "bphs",
                  "cphs",
                  "name_phyto",
                  "col_phyto",
                  "lty_phyto",
                  "bbtilde_phyto",
                  "phs_n440")]

# frac_phyto_lib
load("~/Documents/GitHub/HEShun/data/frac_phyto_lib.rda")

# ag_lib
load("~/Documents/GitHub/HEShun/data/ag_lib.rda")

# B98_AE_midUVabs
load("~/Documents/GitHub/HEShun/data/B98_AE_midUVabs.rda")

# DavidMcKee
load("~/Documents/GitHub/HEShun/data/DavidMcKee.rda")


save(list = ls(), file = "./R/sysdata.rda", compression_level = 9)

rm(list = ls())

devtools::document()


###### test ######

# sub modules
str(WOPP())
str(IOP_d_B22_C2(1, 1))
str(IOP_ph_B22_C2(1))
str(rand_frac())
str(IOP_cdom_B22_C2_lib(1))


str(specific_bph_G99(1))
str(specific_bph_GM83(1))
str(specific_bph_LM88(1))

str(specific_aph_B04(1))
str(specific_aph_B98(1))

# main functions
str(IOP_four_comp(1, 1, 1))
str(IOP_two_comp(1))

# test benchmark models
str(L23_IOP(aph = phytodive_iop$aphs$Brown, wavelen = phytodive_iop$aphs$wv))
str(R18_IOP(1, 1, 1))

IOP_four_comp(0, 0, 0, Sal = 30) %$% spec %>% as.data.table() %>%
  ggplot(aes(x = wavelen, y = Rrs_L11)) + geom_path()

IOP_two_comp(0, Sal = 30, aw_version = 3) %$% spec %>% as.data.table() %>%
  ggplot(aes(x = wavelen, y = Rrs_L11)) + geom_path()






