rm(list = ls())

devtools::document()

load("./R/sysdata.rda")


library(data.table)

output_dir = "~/Documents/GitHub/pyIOPmodel"

fwrite(ag_lib, file.path(output_dir, "data", "ag_hat_lib.csv"))

fwrite(coef_exp_ads, file.path(output_dir, "data", "coef_exp_ads.csv"))

fwrite(frac_phyto_lib, file.path(output_dir, "data", "frac_phyto_lib.csv"))

fwrite(WOPP_computed_refri_T27_S0_180_4000nm,
       file.path(output_dir, "WOPP", "data", "WOPP_computed_refri_T27_S0_180_4000nm.csv"))

fwrite(WOPP_purewater_abs_coefficients,
       file.path(output_dir, "WOPP", "data", "WOPP_purewater_abs_coefficients.csv"))

for(var in c("aphs", "bphs", "cphs")) {
  fwrite(phytodive_iop[[var]], file.path(output_dir, "data", sprintf("Phyto_lib_%s.csv", var)))
}

colhex = NA
for(i in 1:length(phytodive_iop$col_phyto)) {
  tmp <- col2rgb(phytodive_iop$col_phyto[i])
  colhex[i] <- rgb(tmp[1], tmp[2], tmp[3], maxColorValue = 255)
}

Phyto_lib_info = data.table(
  ShortName = names(phytodive_iop$name_phyto),
  LongName = phytodive_iop$name_phyto,
  ColorName = phytodive_iop$col_phyto,
  ColorHex = colhex,
  bbtilde = phytodive_iop$bbtilde_phyto
)

fwrite(Phyto_lib_info, file.path(output_dir, "data", "Phyto_lib_info.csv"))





