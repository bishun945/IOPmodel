rm(list = ls())

# for WOPP
load("~/Documents/GitHub/HEShun/data/WOPP_purewater_abs_coefficients_v3.rda")
load("~/Documents/GitHub/HEShun/data/WOPP_computed_refri_T27_S0_180_4000nm.rda")

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


# test

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








