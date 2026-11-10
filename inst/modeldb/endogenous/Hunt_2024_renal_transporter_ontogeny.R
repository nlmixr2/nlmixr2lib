Hunt_2024_renal_transporter_ontogeny <- function() {
  description <- paste(
    "PBPK system-parameter (renal transporter ontogeny) model. Maximum",
    "likelihood ontogeny profiles of ten renal transporters (OAT1, OAT3,",
    "OCT2, P-gp, URAT1, BCRP, MATE1, MATE-2K, MRP2 and MRP4) from human",
    "post-mortem kidney cortex expression data, preterm neonate to adult",
    "(Hunt 2024). Each profile gives the ontogeny scaling factor (OSF),",
    "the fraction of mature transporter activity, as a lognormal",
    "distribution whose median is a Hill function of postmenstrual age,",
    "OSF = [1 + (a / PMA)^b]^-1 with PMA in years (Equation 1), and whose",
    "log-scale standard deviation shrinks as the transporter matures,",
    "sigma = s * exp(-OSF) (Equation 3). In a PBPK model the OSF",
    "multiplies the adult transporter concentration (Equation 4). This is",
    "a SYSTEM layer intended to be coupled to a drug model: it has no",
    "drug, no dosing, no compartments and no ODEs, and every output is an",
    "algebraic function of the PAGE covariate. Each subject draws one",
    "standard-normal deviate per transporter, so a stochastic solve",
    "samples the age-appropriate OSF distribution and a zeroRe() solve",
    "returns the median profiles.",
    sep = " "
  )
  reference <- paste(
    "Hunt JP, Dubinsky S, McKnite AM, Cheung KWK, van Groen BD,",
    "Giacomini KM, de Wildt SN, Edginton A, Watt KM. Maximum likelihood",
    "estimation of renal transporter ontogeny profiles for pediatric PBPK",
    "modeling. CPT Pharmacometrics Syst Pharmacol. 2024;13:576-588.",
    "doi:10.1002/psp4.13102.",
    "Equations 1-4 are from the main article; the profile parameters a,",
    "b and s are from Supplementary Table S1 and the tabulated profiles",
    "used to check them are Supplementary Tables S2-S11.",
    sep = " "
  )
  vignette <- "Hunt_2024_renal_transporter_ontogeny"
  units <- list(
    time = paste(
      "n/a -- the model has no dynamics. Each record gives the ontogeny",
      "scaling factors implied by its PAGE value, so the time column is",
      "ignored."
    ),
    dosing = "n/a (transporter ontogeny system model; no dosing)",
    concentration = paste(
      "n/a (no drug concentration). Outputs are the unitless ontogeny",
      "scaling factors ont_<transporter> (fraction of mature activity),",
      "their medians ont_<transporter>_med, and the geometric standard",
      "deviations gsd_<transporter> = exp(sigma)."
    )
  )

  covariateData <- list(
    PAGE = list(
      description = "Postmenstrual age (gestational age at birth plus postnatal age)",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Equation 1 is written with x = PMA in YEARS, so this model takes",
        "PAGE in years rather than the register's default of months.",
        "Supply PAGE > 0. The expression data span 25 weeks PMA (about",
        "0.48 years) to over 30 years, except URAT1 and BCRP, which had",
        "no preterm neonate data."
      ),
      source_name = "PMA"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 71L,
    n_studies = 1L,
    age_range = paste(
      "25 weeks postmenstrual age to over 30 years (URAT1 and BCRP: no",
      "preterm neonate samples)"
    ),
    weight_range = "Not reported",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported",
    disease_state = paste(
      "Post-mortem renal cortical tissue samples from biobanks in the",
      "United States and the Netherlands (Cheung 2019,",
      "doi:10.1002/cpt.1516). Transporter protein abundance (and, for",
      "preterm neonates and for MRP2 and MRP4, mRNA expression) was",
      "normalised to the mean of the mature samples, defined as PMA over",
      "5.5 years (288 weeks; n = 23 per transporter)."
    ),
    dose_range = "n/a (system model; no drug)",
    regions = "United States, the Netherlands",
    n_data_points = paste(
      "Per transporter (Supplementary Table S1, protein / preterm mRNA):",
      "OAT3, OAT1, P-gp, MATE1 and MATE-2K 62 / 9; BCRP 61 / 0; MRP2 and",
      "MRP4 0 / 62 (mRNA only); URAT1 62 / 0; OCT2 61 / 9."
    ),
    notes = paste(
      "n_subjects is the largest per-transporter sample count (n = 71",
      "data points, one per tissue donor). The paper also builds",
      "PK-Sim v11 pediatric PBPK models of acyclovir, ciprofloxacin,",
      "furosemide and meropenem that apply the OAT3, MATE1 and P-gp",
      "profiles; those platform models are not part of this file (see",
      "the vignette)."
    )
  )

  # The etas are standard-normal deviates scaled in model() by the
  # PMA-dependent log-scale SD of Equation 3, not variances on a typical
  # value, so no `lont_<transporter>` parameter exists to pair with them.
  paper_specific_etas <- c(
    "etalont_oat3",
    "etalont_oat1",
    "etalont_pgp",
    "etalont_mate1",
    "etalont_mate2k",
    "etalont_bcrp",
    "etalont_mrp2",
    "etalont_mrp4",
    "etalont_urat1",
    "etalont_oct2"
  )

  ini({
    # Supplementary Table S1, first block of a / b / s rows: the
    # lognormally distributed profiles from combined protein and preterm
    # neonate mRNA data -- the final profiles that the paper tabulates in
    # Tables S2-S11 and applies in its PBPK models. a is the PMA (years)
    # at which the median OSF is one half, b is the Hill exponent, and s
    # is the log-scale standard deviation at OSF = 0 (Equation 3).
    # Values marked * in Table S1 were set by manual, iterative
    # adjustment ('Man.') rather than by a converged optimisation.

    # ---- OAT3 (SLC22A8); Table S1 result 'Math.' ----
    oat3_t50 <- 1.081
    label("OAT3: PMA at half of mature activity, a (years)") # Table S1, OAT3, lognormal a = 1.081
    oat3_hill <- 5.948
    label("OAT3: Hill exponent on PMA, b (unitless)") # Table S1, OAT3, lognormal b = 5.948
    oat3_sd0 <- 1.279
    label("OAT3: log-scale SD of the OSF at zero maturity, s (unitless)") # Table S1, OAT3, lognormal s = 1.279

    # ---- OAT1 (SLC22A6); 'Math.' ----
    oat1_t50 <- 1.000
    label("OAT1: PMA at half of mature activity, a (years)") # Table S1, OAT1, lognormal a = 1.000
    oat1_hill <- 7.602
    label("OAT1: Hill exponent on PMA, b (unitless)") # Table S1, OAT1, lognormal b = 7.602
    oat1_sd0 <- 1.31
    label("OAT1: log-scale SD of the OSF at zero maturity, s (unitless)") # Table S1, OAT1, lognormal s = 1.31

    # ---- P-gp (ABCB1); 'Math.' ----
    pgp_t50 <- 1.001
    label("P-gp: PMA at half of mature activity, a (years)") # Table S1, P-gp, lognormal a = 1.001
    pgp_hill <- 4.565
    label("P-gp: Hill exponent on PMA, b (unitless)") # Table S1, P-gp, lognormal b = 4.565
    pgp_sd0 <- 1.066
    label("P-gp: log-scale SD of the OSF at zero maturity, s (unitless)") # Table S1, P-gp, lognormal s = 1.066

    # ---- MATE1 (SLC47A1); 'Math.' ----
    mate1_t50 <- 0.824
    label("MATE1: PMA at half of mature activity, a (years)") # Table S1, MATE1, lognormal a = 0.824
    mate1_hill <- 5.488
    label("MATE1: Hill exponent on PMA, b (unitless)") # Table S1, MATE1, lognormal b = 5.488
    mate1_sd0 <- 1.116
    label("MATE1: log-scale SD of the OSF at zero maturity, s (unitless)") # Table S1, MATE1, lognormal s = 1.116

    # ---- MATE-2K (SLC47A2; 'MATE2' in Table S1); 'Math.' ----
    mate2k_t50 <- 0.661
    label("MATE-2K: PMA at half of mature activity, a (years)") # Table S1, MATE2, lognormal a = 0.661
    mate2k_hill <- 8.221
    label("MATE-2K: Hill exponent on PMA, b (unitless)") # Table S1, MATE2, lognormal b = 8.221
    mate2k_sd0 <- 1.256
    label("MATE-2K: log-scale SD of the OSF at zero maturity, s (unitless)") # Table S1, MATE2, lognormal s = 1.256

    # ---- BCRP (ABCG2); 'Man.'; protein only, no preterm neonate data ----
    bcrp_t50 <- 0.7
    label("BCRP: PMA at half of mature activity, a (years); manually adjusted") # Table S1, BCRP, lognormal a = 0.7*
    bcrp_hill <- 2.7
    label("BCRP: Hill exponent on PMA, b (unitless); manually adjusted") # Table S1, BCRP, lognormal b = 2.7*
    bcrp_sd0 <- 1.1
    label("BCRP: log-scale SD of the OSF at zero maturity, s (unitless); manually adjusted") # Table S1, BCRP, lognormal s = 1.1*

    # ---- MRP2 (ABCC2); 'Man.'; mRNA only ----
    mrp2_t50 <- 1.0
    label("MRP2: PMA at half of mature activity, a (years); manually adjusted") # Table S1, MRP2, lognormal a = 1.0*
    mrp2_hill <- 3.9
    label("MRP2: Hill exponent on PMA, b (unitless); manually adjusted") # Table S1, MRP2, lognormal b = 3.9*
    mrp2_sd0 <- 1.7
    label("MRP2: log-scale SD of the OSF at zero maturity, s (unitless); manually adjusted") # Table S1, MRP2, lognormal s = 1.7*

    # ---- MRP4 (ABCC4); 'Man.'; mRNA only ----
    mrp4_t50 <- 0.6
    label("MRP4: PMA at half of mature activity, a (years); manually adjusted") # Table S1, MRP4, lognormal a = 0.6*
    mrp4_hill <- 5.0
    label("MRP4: Hill exponent on PMA, b (unitless); manually adjusted") # Table S1, MRP4, lognormal b = 5.0*
    mrp4_sd0 <- 1.3
    label("MRP4: log-scale SD of the OSF at zero maturity, s (unitless); manually adjusted") # Table S1, MRP4, lognormal s = 1.3*; Table S9 GeoSD 3.706 at PMA 0.1 implies 1.31 (vignette, Errata)

    # ---- URAT1 (SLC22A12); 'Math.'; protein only, no preterm neonate data ----
    urat1_t50 <- 0.982
    label("URAT1: PMA at half of mature activity, a (years)") # Table S1, URAT1, lognormal a = 0.982
    urat1_hill <- 8.732
    label("URAT1: Hill exponent on PMA, b (unitless)") # Table S1, URAT1, lognormal b = 8.732
    urat1_sd0 <- 1.294
    label("URAT1: log-scale SD of the OSF at zero maturity, s (unitless)") # Table S1, URAT1, lognormal s = 1.294

    # ---- OCT2 (SLC22A2); 'Math.' ----
    oct2_t50 <- 0.902
    label("OCT2: PMA at half of mature activity, a (years)") # Table S1, OCT2, lognormal a = 0.902
    oct2_hill <- 8.670
    label("OCT2: Hill exponent on PMA, b (unitless)") # Table S1, OCT2, lognormal b = 8.670
    oct2_sd0 <- 1.496
    label("OCT2: log-scale SD of the OSF at zero maturity, s (unitless)") # Table S1, OCT2, lognormal s = 1.496

    # ---- Between-subject variability ----
    # Standard-normal deviates, one per transporter. The OSF is lognormal
    # with a log-scale SD that depends on PMA (Equation 3), so the
    # deviate is scaled by that SD in model() rather than carrying a
    # fixed variance here. The paper samples each transporter's OSF
    # independently in PK-Sim, so the deviates are uncorrelated.
    etalont_oat3 ~ fixed(1) # standard normal; scaled by Equation 3 sigma in model()
    etalont_oat1 ~ fixed(1) # standard normal; scaled by Equation 3 sigma in model()
    etalont_pgp ~ fixed(1) # standard normal; scaled by Equation 3 sigma in model()
    etalont_mate1 ~ fixed(1) # standard normal; scaled by Equation 3 sigma in model()
    etalont_mate2k ~ fixed(1) # standard normal; scaled by Equation 3 sigma in model()
    etalont_bcrp ~ fixed(1) # standard normal; scaled by Equation 3 sigma in model()
    etalont_mrp2 ~ fixed(1) # standard normal; scaled by Equation 3 sigma in model()
    etalont_mrp4 ~ fixed(1) # standard normal; scaled by Equation 3 sigma in model()
    etalont_urat1 ~ fixed(1) # standard normal; scaled by Equation 3 sigma in model()
    etalont_oct2 ~ fixed(1) # standard normal; scaled by Equation 3 sigma in model()
  })

  model({
    pma <- PAGE

    # ---- Median OSF, Equation 1: f = [1 + (a / PMA)^b]^-1 ----
    ont_oat3_med <- 1 / (1 + (oat3_t50 / pma)^oat3_hill)
    ont_oat1_med <- 1 / (1 + (oat1_t50 / pma)^oat1_hill)
    ont_pgp_med <- 1 / (1 + (pgp_t50 / pma)^pgp_hill)
    ont_mate1_med <- 1 / (1 + (mate1_t50 / pma)^mate1_hill)
    ont_mate2k_med <- 1 / (1 + (mate2k_t50 / pma)^mate2k_hill)
    ont_bcrp_med <- 1 / (1 + (bcrp_t50 / pma)^bcrp_hill)
    ont_mrp2_med <- 1 / (1 + (mrp2_t50 / pma)^mrp2_hill)
    ont_mrp4_med <- 1 / (1 + (mrp4_t50 / pma)^mrp4_hill)
    ont_urat1_med <- 1 / (1 + (urat1_t50 / pma)^urat1_hill)
    ont_oct2_med <- 1 / (1 + (oct2_t50 / pma)^oct2_hill)

    # ---- Log-scale SD, Equation 3: sigma = s * exp(-f) ----
    sd_oat3 <- oat3_sd0 * exp(-ont_oat3_med)
    sd_oat1 <- oat1_sd0 * exp(-ont_oat1_med)
    sd_pgp <- pgp_sd0 * exp(-ont_pgp_med)
    sd_mate1 <- mate1_sd0 * exp(-ont_mate1_med)
    sd_mate2k <- mate2k_sd0 * exp(-ont_mate2k_med)
    sd_bcrp <- bcrp_sd0 * exp(-ont_bcrp_med)
    sd_mrp2 <- mrp2_sd0 * exp(-ont_mrp2_med)
    sd_mrp4 <- mrp4_sd0 * exp(-ont_mrp4_med)
    sd_urat1 <- urat1_sd0 * exp(-ont_urat1_med)
    sd_oct2 <- oct2_sd0 * exp(-ont_oct2_med)

    # Geometric SD, the 'GeoSD' column of Supplementary Tables S2-S11.
    gsd_oat3 <- exp(sd_oat3)
    gsd_oat1 <- exp(sd_oat1)
    gsd_pgp <- exp(sd_pgp)
    gsd_mate1 <- exp(sd_mate1)
    gsd_mate2k <- exp(sd_mate2k)
    gsd_bcrp <- exp(sd_bcrp)
    gsd_mrp2 <- exp(sd_mrp2)
    gsd_mrp4 <- exp(sd_mrp4)
    gsd_urat1 <- exp(sd_urat1)
    gsd_oct2 <- exp(sd_oct2)

    # ---- Individual OSF: lognormal about the median (Equation 2) ----
    ont_oat3 <- ont_oat3_med * exp(sd_oat3 * etalont_oat3)
    ont_oat1 <- ont_oat1_med * exp(sd_oat1 * etalont_oat1)
    ont_pgp <- ont_pgp_med * exp(sd_pgp * etalont_pgp)
    ont_mate1 <- ont_mate1_med * exp(sd_mate1 * etalont_mate1)
    ont_mate2k <- ont_mate2k_med * exp(sd_mate2k * etalont_mate2k)
    ont_bcrp <- ont_bcrp_med * exp(sd_bcrp * etalont_bcrp)
    ont_mrp2 <- ont_mrp2_med * exp(sd_mrp2 * etalont_mrp2)
    ont_mrp4 <- ont_mrp4_med * exp(sd_mrp4 * etalont_mrp4)
    ont_urat1 <- ont_urat1_med * exp(sd_urat1 * etalont_urat1)
    ont_oct2 <- ont_oct2_med * exp(sd_oct2 * etalont_oct2)
  })
}
