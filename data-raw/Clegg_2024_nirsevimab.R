# Load the Clegg_2024_nirsevimab model from the paper supplement.

# The code here was used to initially prepare the model. Additional
# simplifications and updates were manually performed following the code here
# for loading.

# The control stream below is Supplementary Material 1 of the source paper and
# is not redistributed here; download it from the article to re-run this
# script.  Clegg L, Freshwater E, Leach A, Villafana T, Wahlby Hamren U.
# Population Pharmacokinetics of Nirsevimab in Preterm and Term Infants.
# J Clin Pharmacol. 2024;64(5):555-567. doi:10.1002/jcph.2401

library(nonmem2rx)

nonmem_model <-
  nonmem2rx(
    "Clegg_2024_nirsevimab_jcph2401-sup-0001-suppmat.ctl",
    thetaNames =
      c(
        TH01_CL = "cl",
        TH02_V2 = "vc",
        TH03_Q = "q",
        TH04_V3 = "vp",
        TH05_KA = "ka",
        TH06_F1 = "f_m",
        TH07_PROP = "propSd",
        TH08_ADD = "addSd",
        TH09_CLBETA_CL = "emax_cl_age",
        TH10_CLT50 = "et50_cl_age",
        TH11_CLAL = "e_cl_wt",
        TH12_V2AL = "e_vc_wt",
        TH13_CLADACAT1 = "e_cl_ada",
        TH14_CLRACEN1 = "e_cl_black",
        TH15_CLRACEN2 = "e_cl_asian_amind_mult",
        TH16_V2RACEN3 = "e_vc_asian_amind_mult",
        TH17_CLSEASON2 = "e_cl_season2"
      )
  )
