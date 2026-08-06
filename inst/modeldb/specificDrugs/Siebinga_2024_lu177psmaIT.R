Siebinga_2024_lu177psmaIT <- function() {
  description <- paste(
    "Five-compartment population PK model with a sequentially fitted PSA",
    "pharmacodynamic layer for the PSMA-targeted radioligand",
    "[177Lu]Lu-PSMA-I&T in men with metastatic castration-resistant prostate",
    "cancer treated with ~7.4 GBq per cycle. PK states are central (1),",
    "salivary glands (2), kidneys (3), tumor lesions (4) and a lumped",
    "remaining-tissue compartment (5); the tissue states carry radioactivity",
    "amounts (MBq) and the central state is observed as a concentration",
    "(MBq/L) through a fixed volume V1 = 10.3 L. Renal excretion leaves the",
    "central compartment (k10 = 0.253 1/h, i.e. 2.61 L/h). Salivary-gland",
    "uptake is saturable (Bmax 134 MBq); all other exchange is first order.",
    "Tumor uptake carries a power effect of segmented tumor volume (exponent",
    "1.08), a structural decline over treatment cycles (73%, 50% and 44% of",
    "the cycle-1 rate in cycles 2, 3 and 4-7), interindividual variability",
    "(63.4% CV) and interoccasion variability across cycles (37.8% CV), and",
    "salivary-gland uptake carries a weak tumor-volume power effect. All PK",
    "parameters are allometrically scaled to the 79 kg median body weight.",
    "The PD layer is a single PSA state growing exponentially (kG 0.000408",
    "1/h) and eliminated by a direct effect linear in the physical tumor",
    "activity concentration plus a delayed effect linear in an",
    "effect-compartment concentration (ke0 0.00128 1/h). Baseline PSA is",
    "fixed at 140 ug/L with a linear tumor-volume effect. Every state",
    "carries a decay-corrected radioactivity amount, matching the source",
    "SPECT/CT data and the predecessor [177Lu]Lu-PSMA-617 model.",
    sep = " "
  )
  reference <- paste(
    "Siebinga H, de Wit-van der Veen BJ, de Vries-Huizing DMV, Vogel WV,",
    "Hendrikx JJMA, Huitema ADR. Quantification of biochemical PSA dynamics",
    "after radioligand therapy with [177Lu]Lu-PSMA-I&T using a population",
    "pharmacokinetic/pharmacodynamic model. EJNMMI Phys. 2024;11(1):39.",
    "doi:10.1186/s40658-024-00642-2",
    sep = " "
  )
  vignette <- "Siebinga_2024_lu177psmaIT"
  units <- list(time = "h", dosing = "MBq", concentration = "MBq/L")

  covariateData <- list(
    WT = list(
      description = "Total body weight at baseline (Table 1 median 79 kg, range 61-116).",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling was added to all PK parameters (Results,",
        "'Population pharmacokinetic model'; dOFV -11.3). Methods state that",
        "body weight was scaled 'in relation to the median body weight, where",
        "exponents for volume of distribution (V) and rate constants (k) were",
        "set to 1 and -0.25 respectively', so both exponents are fixed a",
        "priori rather than estimated and the reference is the Table 1 median",
        "of 79 kg. Bmax is neither a volume nor a rate constant and the paper",
        "names only V and k, so Bmax is left unscaled. Baseline value, held",
        "constant per subject.",
        sep = " "
      ),
      source_name = "Weight"
    ),
    TUM_VOL = list(
      description = paste(
        "Tumor volume of the segmented (target) lesions, determined on",
        "pre-treatment diagnostic PSMA-PET/CT with a 45% SULmax threshold.",
        "Reported in L in Table 1 (median 0.0443 L, range 0.000122-0.546);",
        "carried here in the canonical mm^3 (1 L = 1e6 mm^3), so the",
        "reference is 44300 mm^3.",
        sep = " "
      ),
      units = "mm^3",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters three places. (1) Tumor uptake rate as a power function",
        "(Table 2 'Tumor volume on k14' = 1.08); the paper's own check is",
        "that a two-fold volume gives a 2.11-fold uptake rate, and",
        "2^1.08 = 2.114. (2) Salivary-gland uptake rate as a power function",
        "(Table 2 'Tumor volume on k12' = 0.0910). (3) Baseline PSA as a",
        "linear function (Table 3, 57.5 ug/L). The same volume also converts",
        "the tumor radioactivity amount to the concentration that drives the",
        "PD layer. Only target-lesion volume was available, not total tumor",
        "volume (Methods, 'Pharmacokinetic model development'). Baseline",
        "value, held constant per subject.",
        sep = " "
      ),
      source_name = "Tumor volume of segmented tumors"
    ),
    CYCLE = list(
      description = "Treatment-cycle number (1 = first ~7.4 GBq administration). Patients received up to eight cycles (Table 1).",
      units = "(count)",
      type = "count",
      reference_category = NULL,
      notes = paste(
        "Carries two distinct effects on the tumor uptake rate constant and",
        "is decomposed inside model() into binary indicators rather than",
        "entering as an integer. (1) The structural cycle effect of Equation",
        "8: uptake in each cycle is a fraction of the cycle-1 uptake, with",
        "cycles 4 and later lumped into a single level because few patients",
        "received more than four cycles (Table 2: 0.731, 0.498, 0.436 for",
        "cycles 2, 3 and 4-7). (2) The occasion index for interoccasion",
        "variability, which the paper added only to the tumor uptake rate",
        "(Table 2 IOV k14 = 37.8% CV). Occasions 1-8 are provided to cover",
        "the observed range of cycles.",
        sep = " "
      ),
      source_name = "Number of cycles received"
    )
  )

  # Screened during covariate analysis but NOT retained in the final model, so
  # they are documented here rather than in covariateData (they are never
  # referenced in model()).
  covariatesDataExcluded <- list(
    CRCL = list(
      description = "Glomerular filtration rate estimated with the Cockcroft-Gault equation (Table 1 median 80.9 mL/min, range 25.9-181).",
      units = "mL/min",
      type = "continuous",
      notes = paste(
        "Tested as a linear covariate on the renal excretion rate constant",
        "k10, which the paper considered pharmacologically the most suitable",
        "form because glomerular filtration is the only elimination route.",
        "It was not retained: 'Serum creatinine clearance was not identified",
        "as a covariate on k10, since the OFV and model fit did not improve",
        "(dOFV -0.207)' (Results). No coefficient is reported, so the effect",
        "cannot be encoded. Note that the predecessor [177Lu]Lu-PSMA-617",
        "model (Siebinga_2023_lu177psma617.R) did scale k10 by CRCL a priori.",
        sep = " "
      ),
      source_name = "GFR (calculated using Cockcroft Gault)"
    ),
    AGE = list(
      description = "Age at baseline (Table 1 median 73 years, range 48-91).",
      units = "years",
      type = "continuous",
      notes = paste(
        "Tested as a power and a linear covariate on baseline PSA on the",
        "hypothesis that age tracks disease severity; not retained ('Age did",
        "not explain IIV on the baseline PSA (dOFV -1.01)', Results). No",
        "coefficient is reported.",
        sep = " "
      ),
      source_name = "Age"
    ),
    HCT = list(
      description = "Hematocrit at baseline (Table 1 median 0.35, range 0.22-0.44).",
      units = "(fraction)",
      type = "continuous",
      notes = "Reported among the baseline characteristics in Table 1 but never tested as a covariate and not part of any model.",
      source_name = "Hematocrit"
    )
  )

  compartmentData <- list(
    central         = list(analyte = "[177Lu]Lu-PSMA-I&T radioactivity", units = "MBq", specimen = "whole blood", verified = TRUE),
    salivary_gland  = list(analyte = "[177Lu]Lu-PSMA-I&T radioactivity", units = "MBq", specimen = "tissue", verified = TRUE),
    kidney          = list(analyte = "[177Lu]Lu-PSMA-I&T radioactivity", units = "MBq", specimen = "tissue", verified = TRUE),
    tumor           = list(analyte = "[177Lu]Lu-PSMA-I&T radioactivity", units = "MBq", specimen = "tumor", verified = TRUE),
    other           = list(analyte = "[177Lu]Lu-PSMA-I&T radioactivity", units = "MBq", specimen = "tissue", verified = FALSE),
    effect          = list(analyte = "[177Lu]Lu-PSMA-I&T effect-compartment activity concentration", units = "MBq/L", specimen = "not applicable", verified = FALSE),
    PSA             = list(analyte = "prostate-specific antigen", units = "ug/L", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 76,
    n_studies = 1,
    age_range = "48-91 years",
    age_median = "73 years",
    weight_range = "61-116 kg",
    weight_median = "79 kg",
    sex_female_pct = 0,
    disease_state = paste(
      "Metastatic castration-resistant prostate cancer with PSMA-positive",
      "lesions on pre-treatment diagnostic [68Ga]Ga-PSMA PET/CT and no",
      "PSMA-negative lesions on contrast-enhanced CT (VISION criteria).",
      "Baseline PSA median 260 ug/L (0.12-4896); segmented target-tumor",
      "volume median 0.0443 L (0.000122-0.546); hematocrit median 0.35.",
      sep = " "
    ),
    renal_function = "Cockcroft-Gault GFR median 80.9 mL/min (range 25.9-181)",
    dose_range = paste(
      "Intravenous [177Lu]Lu-PSMA-I&T, median injected activity 7378 MBq",
      "(5605-7763) per cycle, up to eight cycles. Two schedules were used:",
      "four cycles six weeks apart ('4 x 6', n = 10) and two cycles two weeks",
      "apart repeated after twelve weeks ('2 x 2 - repeated after twelve",
      "weeks', n = 66).",
      sep = " "
    ),
    regions = "Netherlands (Netherlands Cancer Institute / Antoni van Leeuwenhoek; IRBd21288)",
    notes = paste(
      "Retrospective cohort treated September 2019 to January 2023. The PK",
      "model was built on 409 quantitative SPECT/CT scans from 76 patients",
      "(3 of 79 screened patients had only lesions < 20 mm and were",
      "excluded); scans were acquired at 4.6 +/- 0.95 h, 23.8 +/- 4.5 h and",
      "6.8 +/- 0.46 days post injection. Kidney, salivary-gland and tumor",
      "readings were lumped into one state each; whole-body data were not",
      "available. All radioactivity data were corrected for decay back to the",
      "time of injection, and SPECT-derived blood concentrations were",
      "recalibrated with the [177Lu]Lu-PSMA-617 regression (intercept 6.27",
      "MBq/L, slope 0.828) taken from the predecessor model. The PD layer was",
      "fitted sequentially to 566 PSA observations.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural PK parameters -- all from Table 2 ("Model estimates of the
    # final pharmacokinetic model"). The paper's numeric compartment indices
    # are 1 = central, 2 = salivary gland, 3 = kidney, 4 = tumor and
    # 5 = remaining tissue; the kin_<tissue> / kout_<tissue> names below
    # carry the same meaning without relying on the index.
    # ------------------------------------------------------------------
    lkel <- log(0.253); label("Renal excretion rate constant from central, k10 (1/h)")                      # Table 2: k10 = 0.253 (RSE 5.2%); k10 * V1 = 2.61 L/h, the renal clearance quoted in the Discussion
    lvc <- fixed(log(10.3)); label("Central volume of distribution, V1 (L)")                                # Table 2: V1 = 10.3, footnote c "Fixed parameter" -- blood data from scans were too unreliable to estimate it

    lkin_salivary_gland <- log(0.0105); label("Central-to-salivary-gland uptake rate constant, k12 (1/h)")  # Table 2: k12 = 0.0105 (RSE 9.2%)
    lkout_salivary_gland <- log(0.0629); label("Salivary-gland-to-central rate constant, k21 (1/h)")        # Table 2: k21 = 0.0629 (RSE 6.3%); log(2)/0.0629 = 11.0 h, the "organs-at-risk ~11 h" half-life in the Discussion
    lbmax <- log(134); label("Maximum salivary-gland binding capacity, Bmax (MBq)")                         # Table 2: Bmax compartment 2 = 134 (RSE 15.5%)

    lkin_kidney <- log(0.0321); label("Central-to-kidney uptake rate constant, k13 (1/h)")                  # Table 2: k13 = 0.0321 (RSE 6.9%)
    lkout_kidney <- log(0.0625); label("Kidney-to-central rate constant, k31 (1/h)")                        # Table 2: k31 = 0.0625 (RSE 5.4%)

    lkin_tumor <- log(0.00967); label("Central-to-tumor uptake rate constant, k14 (1/h)")                   # Table 2: k14 = 0.00967 (RSE 9.9%)
    lkout_tumor <- log(0.0150); label("Tumor-to-central rate constant, k41 (1/h)")                          # Table 2: k41 = 0.0150 (RSE 3.4%); log(2)/0.0150 = 46.2 h, the "tumors ~46 h" half-life in the Discussion

    lkin_other <- log(0.275); label("Central-to-remaining-tissue rate constant, k15 (1/h)")                 # Table 2: k15 = 0.275 (RSE 6.8%)
    lkout_other <- log(0.0247); label("Remaining-tissue-to-central rate constant, k51 (1/h)")               # Table 2: k51 = 0.0247 (RSE 6.5%)

    # ------------------------------------------------------------------
    # PK covariate effects.
    # ------------------------------------------------------------------
    e_tum_vol_kin_tumor <- 1.08; label("Power exponent of segmented tumor volume on kin_tumor (unitless)")            # Table 2 "Tumor volume on k14" = 1.08 (RSE 7.5%), footnote a power function; 2^1.08 = 2.114 reproduces the paper's "two-times higher volume -> 2.11-fold increased tumor uptake rate"
    e_tum_vol_kin_salivary_gland <- 0.0910; label("Power exponent of segmented tumor volume on kin_salivary_gland (unitless)") # Table 2 "Tumor volume on k12" = 0.0910 (RSE 30.0%), footnote a power function (see the vignette Errata for the sign conflict with the Results narrative)

    # Equation 8 fractional cycle effect: uptake in cycle n is a fraction of
    # the cycle-1 uptake. Cycles 4 and later share one level.
    e_cycle2_kin_tumor <- 0.731; label("Fraction of cycle-1 tumor uptake rate retained in cycle 2 (unitless)")        # Table 2 "Cycle 2 on k14" = 0.731 (RSE 8.1%); Abstract and Conclusions round this to 73%
    e_cycle3_kin_tumor <- 0.498; label("Fraction of cycle-1 tumor uptake rate retained in cycle 3 (unitless)")        # Table 2 "Cycle 3 on k14" = 0.498 (RSE 11.3%); rounded to 50% in the Abstract
    e_cycle47_kin_tumor <- 0.436; label("Fraction of cycle-1 tumor uptake rate retained in cycles 4-7 (unitless)")    # Table 2 "Cycle 4-7 on k14" = 0.436 (RSE 11.6%); rounded to 44% in the Abstract

    # Allometric scaling on median body weight. Both exponents were set a
    # priori, not estimated, so both are fixed.
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on the central volume (unitless)")                 # Methods: "exponents for volume of distribution (V) and rate constants (k) were set to 1 and -0.25, respectively"
    e_wt_k <- fixed(-0.25); label("Allometric exponent of body weight on every rate constant (unitless)")             # Methods: as above; applied to all rate constants, and Results confirm "Allometric scaling was added to all PK parameters (dOFV -11.3)"

    # ------------------------------------------------------------------
    # Structural PD parameters -- all from Table 3 ("Pharmacodynamic model
    # estimates of the final model").
    # ------------------------------------------------------------------
    rbase <- fixed(140); label("Baseline PSA concentration (ug/L)")                                                   # Table 3: Baseline PSA = 140, footnote a "Fixed parameter" -- estimated at 140 then fixed to avoid overparameterization
    e_tum_vol_rbase <- 57.5; label("Linear effect of segmented tumor volume on baseline PSA (ug/L)")                  # Table 3 "Tumor volume on baseline PSA" = 57.5 (RSE 38.9%), footnote b linear function Pcov = Ppop + theta * (COV / COVmedian)
    lkgrow <- log(0.000408); label("First-order exponential PSA growth rate, kG (1/h)")                               # Table 3: PSA growth rate kG = 0.000408 (RSE 14.2%)
    lkd_direct <- log(0.00335); label("Direct drug-induced PSA elimination coefficient, kD,direct (L/day/GBq)")       # Table 3: kD,direct = 0.00335 (RSE 40.1%, 95%CI 0.000961-0.006147). The Results text prints 0.000335; Table 3's estimate, RSE and CI are mutually consistent and the text value is not -- see the vignette Errata
    lke0 <- log(0.00128); label("Effect-compartment equilibration rate constant, ke0 (1/h)")                          # Table 3: ke0 = 0.00128 (RSE 13.1%)
    lkd_delay <- log(0.0000328); label("Delayed drug-induced PSA elimination coefficient, kD,delay (L/day/MBq)")      # Table 3: kD,delay = 0.0000328 (RSE 17.4%)
    boxcox_lkgrow <- -0.822; label("Box-Cox shape parameter for the interindividual distribution of kG (unitless)")   # Table 3: Box-cox shape parameter = -0.822 (RSE 28.3%); applied to the kG eta only (Results: "For IIV on kG, a box-cox transformation was applied")

    # ------------------------------------------------------------------
    # Interindividual variability. Tables 2 and 3 report IIV as CV%; the
    # internal log-normal variances below are omega^2 = log(CV^2 + 1),
    # because Equation 2 is the exponential model P_i = P_pop * exp(eta_i).
    # ------------------------------------------------------------------
    etalkel ~ 0.1039644              # Table 2 IIV k10  = 33.1% CV -> log(0.331^2 + 1)
    etalkin_kidney ~ 0.1057608       # Table 2 IIV k13  = 33.4% CV -> log(0.334^2 + 1)
    etalkin_tumor ~ 0.3378684        # Table 2 IIV k14  = 63.4% CV -> log(0.634^2 + 1)
    etalkin_other ~ 0.0812859        # Table 2 IIV k15  = 29.1% CV -> log(0.291^2 + 1)
    etalbmax ~ 0.3708046             # Table 2 IIV Bmax compartment 2 = 67.0% CV -> log(0.670^2 + 1)

    etarbase ~ 1.4360602             # Table 3 IIV baseline PSA = 179% CV -> log(1.79^2 + 1)
    etalkgrow ~ 0.6022817            # Table 3 IIV kG = 90.9% CV -> log(0.909^2 + 1); this eta is Box-Cox transformed in model()
    etalkd_direct ~ 1.0851893        # Table 3 IIV kD,direct = 140% CV -> log(1.40^2 + 1)
    etalkd_delay ~ 0.5685047         # Table 3 IIV kD,delay = 87.5% CV -> log(0.875^2 + 1)

    # Interoccasion variability on the tumor uptake rate constant, one
    # occasion per treatment cycle (Table 2 IOV k14 = 37.8% CV). nlmixr2 has
    # no NONMEM $OMEGA BLOCK(1) SAME shortcut, so occasions 2-8 each get
    # their own eta with the variance fixed equal to the occasion-1 estimate
    # (Siebinga_2023_lu177psma617.R / Jonsson_2011_ethambutol.R precedent).
    # Eight occasions cover the maximum number of cycles received (Table 1).
    etaiov_kin_tumor_1 ~ 0.1335549
    etaiov_kin_tumor_2 ~ fix(0.1335549)
    etaiov_kin_tumor_3 ~ fix(0.1335549)
    etaiov_kin_tumor_4 ~ fix(0.1335549)
    etaiov_kin_tumor_5 ~ fix(0.1335549)
    etaiov_kin_tumor_6 ~ fix(0.1335549)
    etaiov_kin_tumor_7 ~ fix(0.1335549)
    etaiov_kin_tumor_8 ~ fix(0.1335549)

    # ------------------------------------------------------------------
    # Residual unexplained variability. Equation 3 (combined) applies to the
    # central compartment and Equation 4 (proportional) to the observed
    # tissue compartments and to PSA. Table 2 reports the proportional terms
    # as CV%, which for these linear-scale models is the proportional SD.
    # The remaining-tissue compartment was never observed and so has no
    # residual-error term.
    # ------------------------------------------------------------------
    propSd <- 0.555; label("Proportional residual error, central concentration (fraction)")                # Table 2 RUV proportional error compartment 1 = 55.5% CV (RSE 29.0%)
    addSd <- 9.57; label("Additive residual error, central concentration (MBq/L)")                         # Table 2 RUV additive error compartment 1 = 9.57 MBq/L (RSE 16.8%)
    propSd_Asalivary_gland <- 0.397; label("Proportional residual error, salivary-gland activity (fraction)") # Table 2 RUV proportional error compartment 2 = 39.7% CV (RSE 8.6%)
    propSd_Akidney <- 0.319; label("Proportional residual error, kidney activity (fraction)")               # Table 2 RUV proportional error compartment 3 = 31.9% CV (RSE 8.5%)
    propSd_Atumor <- 0.327; label("Proportional residual error, tumor activity (fraction)")                 # Table 2 RUV proportional error compartment 4 = 32.7% CV (RSE 9.7%)
    propSd_PSA <- 0.293; label("Proportional residual error, PSA (fraction)")                               # Table 3 RUV proportional error = 29.3% CV (RSE 9.2%)
  })

  model({
    # ------------------------------------------------------------------
    # Population parameters that are used outside a mu-referenced
    # exp(theta + eta) position are first aliased to locals, each on its own
    # simple line. rxode2's mu-referencing pass otherwise fails with
    # "subscript out of bounds" (Siebinga_2023_lu177psma617.R precedent).
    # Aliasing changes nothing numerically.
    aloV <- e_wt_vc
    aloK <- e_wt_k
    rbasePop <- rbase
    eTumVolRbase <- e_tum_vol_rbase
    bcLambda <- boxcox_lkgrow

    # ------------------------------------------------------------------
    # Cycle decomposition. CYCLE is carried as an integer but never enters
    # an equation as one: it is expanded into the Equation 8 structural
    # indicators and into the interoccasion-variability multiplexer.
    # ------------------------------------------------------------------
    cyc2 <- (CYCLE == 2)
    cyc3 <- (CYCLE == 3)
    cyc47 <- (CYCLE >= 4)

    fcyc2 <- e_cycle2_kin_tumor
    fcyc3 <- e_cycle3_kin_tumor
    fcyc47 <- e_cycle47_kin_tumor
    cycleFactor <- fcyc2^cyc2 * fcyc3^cyc3 * fcyc47^cyc47

    iov_kin_tumor <-
      (CYCLE == 1) * etaiov_kin_tumor_1 +
      (CYCLE == 2) * etaiov_kin_tumor_2 +
      (CYCLE == 3) * etaiov_kin_tumor_3 +
      (CYCLE == 4) * etaiov_kin_tumor_4 +
      (CYCLE == 5) * etaiov_kin_tumor_5 +
      (CYCLE == 6) * etaiov_kin_tumor_6 +
      (CYCLE == 7) * etaiov_kin_tumor_7 +
      (CYCLE >= 8) * etaiov_kin_tumor_8

    # ------------------------------------------------------------------
    # Allometric scaling on the 79 kg median body weight (Methods): the
    # volume scales with exponent 1 and every rate constant with -0.25.
    # ------------------------------------------------------------------
    alloV <- (WT / 79)^aloV
    alloK <- (WT / 79)^aloK

    # ------------------------------------------------------------------
    # Individual PK parameters. Equation 2 is the exponential IIV / IOV
    # model; Equation 5 is the power covariate model.
    # ------------------------------------------------------------------
    kel <- exp(lkel + etalkel) * alloK
    vc <- exp(lvc) * alloV

    kin_salivary_gland <- exp(lkin_salivary_gland) * alloK * (TUM_VOL / 44300)^e_tum_vol_kin_salivary_gland
    kout_salivary_gland <- exp(lkout_salivary_gland) * alloK
    bmax <- exp(lbmax + etalbmax)

    kin_kidney <- exp(lkin_kidney + etalkin_kidney) * alloK
    kout_kidney <- exp(lkout_kidney) * alloK

    kin_tumor <- exp(lkin_tumor + etalkin_tumor + iov_kin_tumor) * alloK *
      (TUM_VOL / 44300)^e_tum_vol_kin_tumor * cycleFactor
    kout_tumor <- exp(lkout_tumor) * alloK

    kin_other <- exp(lkin_other + etalkin_other) * alloK
    kout_other <- exp(lkout_other) * alloK

    # ------------------------------------------------------------------
    # Individual PD parameters. The tumor-volume effect on baseline PSA is
    # the Equation 6 linear form, applied to the population value before the
    # Equation 2 exponential IIV.
    # ------------------------------------------------------------------
    basePop <- rbasePop + eTumVolRbase * (TUM_VOL / 44300)
    psa0 <- basePop * exp(etarbase)

    # Box-Cox transformed eta on the PSA growth rate (Results: the random
    # effect deviated strongly from log-normal). Petersson 2009 form:
    # phi = (exp(eta)^lambda - 1) / lambda, then kG = kG_pop * exp(phi).
    phi_lkgrow <- (exp(etalkgrow)^bcLambda - 1) / bcLambda
    kgrow <- exp(lkgrow + phi_lkgrow)

    ke0 <- exp(lke0)
    kd_direct <- exp(lkd_direct + etalkd_direct)
    kd_delay <- exp(lkd_delay + etalkd_delay)

    # ------------------------------------------------------------------
    # Saturable salivary-gland binding, Equation 1: uptake into the target
    # is scaled by the unoccupied fraction of the maximum binding capacity
    # while release back to the central compartment stays first order. Only
    # the salivary gland was retained as saturable -- Bmax was
    # non-identifiable for the kidney and tumor compartments.
    # ------------------------------------------------------------------
    fluxSalivaryGland <- kin_salivary_gland * central * (1 - salivary_gland / bmax)

    # ------------------------------------------------------------------
    # Five-compartment radioligand disposition (Figure 1). Every state is a
    # decay-corrected radioactivity amount in MBq, matching the source data
    # ("All radioactivity data were corrected for decay to the time of
    # injection", Methods) and the predecessor model
    # Siebinga_2023_lu177psma617.R. Only the central state is converted to
    # a concentration.
    # ------------------------------------------------------------------
    d/dt(central) <-
      -kel * central -
      fluxSalivaryGland + kout_salivary_gland * salivary_gland -
      kin_kidney * central + kout_kidney * kidney -
      kin_tumor * central + kout_tumor * tumor -
      kin_other * central + kout_other * other
    d/dt(salivary_gland) <- fluxSalivaryGland - kout_salivary_gland * salivary_gland
    d/dt(kidney) <- kin_kidney * central - kout_kidney * kidney
    d/dt(tumor) <- kin_tumor * central - kout_tumor * tumor
    d/dt(other) <- kin_other * central - kout_other * other

    # ------------------------------------------------------------------
    # PD layer. Ctumor is the estimated activity concentration in the tumor
    # compartment (Methods, below Equation 11); the segmented tumor volume
    # converts the MBq amount to MBq/L (1 L = 1e6 mm^3). Equation 9 gives
    # the direct effect and the effect compartment gives the delayed one.
    # kD,direct is reported per GBq and kD,delay per MBq, and both are per
    # day, so each is converted to the model's 1/h time base.
    # ------------------------------------------------------------------
    vtumor <- TUM_VOL / 1e6
    Ctumor <- tumor / vtumor
    d/dt(effect) <- ke0 * (Ctumor - effect)

    edrug <- kd_direct * Ctumor / 1000 / 24
    edelayed <- kd_delay * effect / 24

    # Equation 12. No PSA decrease occurs before the first administration
    # because both effects are driven by tumor activity, which is zero.
    d/dt(PSA) <- kgrow * PSA - (edrug + edelayed) * PSA
    PSA(0) <- psa0

    # ------------------------------------------------------------------
    # Observation channels. Every radioactivity record in the source data
    # was corrected for decay back to the time of injection, so all states
    # and all four channels are on the decay-corrected scale. The
    # remaining-tissue compartment was never observed and has no channel.
    # ------------------------------------------------------------------
    Cc <- central / vc
    Asalivary_gland <- salivary_gland
    Akidney <- kidney
    Atumor <- tumor

    Cc ~ add(addSd) + prop(propSd)
    Asalivary_gland ~ prop(propSd_Asalivary_gland)
    Akidney ~ prop(propSd_Akidney)
    Atumor ~ prop(propSd_Atumor)
    PSA ~ prop(propSd_PSA)
  })
}
