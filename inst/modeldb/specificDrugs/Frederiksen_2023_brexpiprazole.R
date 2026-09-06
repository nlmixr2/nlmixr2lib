Frederiksen_2023_brexpiprazole <- function() {
  description <- paste(
    "Joint parent + two-metabolite population PK model for oral",
    "brexpiprazole and its major (DM-3411, CYP3A4/CYP2D6-formed) and",
    "minor (DM-3412, CYP2D6-formed) metabolites, pooled across 13 phase I",
    "and phase II studies in 826 healthy subjects and patients with",
    "schizophrenia, major depressive disorder, or attention deficit",
    "hyperactivity disorder. Brexpiprazole is two-compartment with",
    "first-order absorption and an absorption lag; DM-3412 is",
    "two-compartment and DM-3411 is one-compartment, each fed by its own",
    "formation clearance from the brexpiprazole central compartment, so",
    "the parent's total elimination is the sum of the two formation",
    "clearances. Covariates retained: food state on the absorption rate",
    "constant, and body-mass index on the DM-3411 central volume. The",
    "brexpiprazole-side parameters (ka, lag time, V2, V3, Q, and the IIV",
    "on ka and V3) were fixed to a prior brexpiprazole-only model to",
    "stabilise the joint fit. The paper's purpose was to use the",
    "DM-3412:brexpiprazole metabolic ratio, which reduces to",
    "CL_form_DM3412 / CL_DM3412, as an in vivo surrogate of CYP2D6",
    "activity for CYP2D6 allele activity-score estimation; the CYP2D6",
    "genotype itself is NOT a covariate in the population PK model.",
    "Fit in NONMEM 7.4."
  )
  reference <- paste(
    "Frederiksen T, Areberg J, Raoufinia A, Schmidt E, Stage TB,",
    "Brosen K. Estimating the In Vivo Function of CYP2D6 Alleles through",
    "Population Pharmacokinetic Modeling of Brexpiprazole.",
    "Clin Pharmacol Ther. 2023 Feb;113(2):360-369.",
    "doi:10.1002/cpt.2791. PMCID: PMC10099095.",
    sep = " "
  )
  vignette <- "Frederiksen_2023_brexpiprazole"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    FED = list(
      description        = "Fed state at the time of dosing (1 = fed, 0 = fasted)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = paste(
        "Frederiksen 2023 Table 2 reports two absorption rate constants,",
        "'Absorption rate constant fasted' = 1.02 /h (fixed) and",
        "'Absorption rate constant fed' = 0.535 /h (23.4 %RSE); the model",
        "here encodes the pair as the canonical linear-deviation form",
        "ka * (1 + e_fed_ka * FED), so e_fed_ka = 0.535/1.02 - 1.",
        "Food effect on ka was the single most significant covariate in",
        "the forward-inclusion step (OFV -281; Results, 'Population",
        "pharmacokinetic analysis'). The Discussion notes the effect",
        "appeared to be driven by subjects from one multiple-ascending-dose",
        "study and that a dedicated food-effect study found no impact on",
        "brexpiprazole PK, so the coefficient should be read as a",
        "study-confounded absorption-rate difference rather than an",
        "established prandial effect. Per dose record.",
        "One of the 13 pooled studies (Table S1) was an 'AME, food effect'",
        "study in 24 healthy men at a single 2 mg dose."
      ),
      source_name        = "food effect"
    ),
    BMI = list(
      description        = "Body mass index at baseline",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Applied to the DM-3411 central volume V6 as the linear-deviation",
        "form vc_dm3411 * (1 + e_bmi_vc_dm3411 * (BMI - 26)), the form",
        "documented for BMI in inst/references/covariate-columns.md.",
        "Frederiksen 2023 Table 2 reports the coefficient 'Body mass index",
        "on V6' = -0.0107 (95.3 %RSE, 95% CI [-0.0174; -0.0000303]) but",
        "does NOT print the covariate equation or the centring value.",
        "The reference 26 kg/m^2 is the population-PK-analysis-set median",
        "BMI (Table 1, n = 782; IQR 23-30, range 16-57). A power form",
        "(BMI/26)^-0.0107 is ruled out because it would change V6 by under",
        "1% across the whole observed BMI range and could not produce the",
        "reported 19-point OFV drop. BMI on V6 was the second and last",
        "covariate retained (forward inclusion OFV -19; retained in",
        "backward elimination). Time-fixed at baseline."
      ),
      source_name        = "BMI"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Table 1: median 76 kg, IQR 65-90, range 32-159 (n = 782).",
        "Screened in the stepwise covariate analysis (Methods,",
        "'Population pharmacokinetic modeling'); not retained."
      )
    ),
    HT = list(
      description = "Body height",
      units       = "cm",
      type        = "continuous",
      notes       = paste(
        "Table 1: median 170 cm, IQR 164-177, range 131-196 (n = 782).",
        "Screened; not retained."
      )
    ),
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Table 1: median 36 years, IQR 27-46, range 18-65 (n = 826).",
        "Screened; not retained."
      )
    ),
    SEXF = list(
      description        = "Female sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Source column was 'gender'. Table 1: 295 of 826 female (35.7%).",
        "Screened; not retained."
      )
    ),
    CRCL = list(
      description = "Creatinine clearance estimated by the Cockcroft-Gault formula",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Table 1: median 119 mL/min, IQR 99-141, range 42-333 (n = 782);",
        "Table 1 footnote a. Screened; not retained."
      )
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "IU/L",
      type        = "continuous",
      notes       = paste(
        "Reported as ALAT in Table 1: median 20 IU/L, IQR 15-28,",
        "range 5-107 (n = 791). Screened; not retained."
      )
    ),
    DOSE_BREX_MG = list(
      description = "Administered brexpiprazole dose",
      units       = "mg",
      type        = "continuous",
      notes       = paste(
        "Dose was screened as a covariate (Methods) to test for",
        "dose-dependent, i.e. non-linear, PK; not retained, so the final",
        "model is linear in dose across 0.15-12 mg."
      )
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "Brexpiprazole", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "Brexpiprazole", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "Brexpiprazole", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    central_dm3412 = list(
      analyte = "DM-3412", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1_dm3412 = list(
      analyte = "DM-3412", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    central_dm3411 = list(
      analyte = "DM-3411", units = "mg",
      specimen = "plasma", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 826,
    n_studies      = 13,
    age_range      = "18-65 years",
    age_median     = "36 years (IQR 27-46)",
    weight_range   = "32-159 kg",
    weight_median  = "76 kg (IQR 65-90)",
    height_range   = "131-196 cm (median 170, IQR 164-177)",
    bmi_range      = "16-57 kg/m^2 (median 26, IQR 23-30)",
    sex_female_pct = 35.7,
    race_ethnicity = c(White = 61.1, Black = 17.9, Asian = 19.5, Other = 1.5),
    disease_state  = paste(
      "Pooled healthy subjects (179, 21.7%) and patients with",
      "schizophrenia (362, 43.8%), major depressive disorder",
      "(127, 15.4%), or attention deficit hyperactivity disorder",
      "(158, 19.1%); Table 1. Patients on moderate or strong CYP2D6",
      "inhibitors (fluoxetine, paroxetine, sertraline) were excluded from",
      "the analysis and concomitant CYP3A4 inducers and inhibitors were",
      "disallowed in all studies."
    ),
    dose_range     = paste(
      "Oral brexpiprazole 0.15-12 mg, single dose (phase I ascending-dose",
      "and PET studies) and multiple dose (multiple-ascending-dose and",
      "phase II studies), as monotherapy (N = 541) or adjunctive treatment",
      "(N = 285); Table S1."
    ),
    regions        = "United States, Europe, Japan, Korea",
    genotype       = paste(
      "CYP2D6 genotyped in all studies (TaqMan, DNA sequencing, and",
      "gel-based assays). CYP2D6 predicted phenotype in the population PK",
      "set: ultra-rapid 15 (1.8%), normal 369 (44.7%), intermediate 215",
      "(26.0%), poor 23 (2.8%), missing 204 (24.7%). CYP2D6 genotype is",
      "NOT a covariate in the population PK model; it enters the paper",
      "only through the downstream metabolic-ratio regression."
    ),
    notes          = paste(
      "Nine phase I and four phase II studies (Table S1). Plasma sampling",
      "was rich (>20 samples/subject, N = 245), semi-sparse (8/subject,",
      "N = 314), or sparse (3-4/subject, N = 267). Brexpiprazole, DM-3411,",
      "and DM-3412 assayed by LC-MS/MS, linear 0.300-100 ng/mL for all",
      "three analytes; LLOQ 0.300 ng/mL. The final data set held 8,604",
      "brexpiprazole, 8,156 DM-3411, and 5,113 DM-3412 concentrations",
      "above the LLOQ. 160 of the 826 subjects had no DM-3412",
      "concentration above the LLOQ. The CYP2D6 genotype-phenotype",
      "sub-analysis used the 496 subjects who had both measurable DM-3412",
      "and genotype data (Table 1, right-hand columns)."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Brexpiprazole (parent) structural PK. Two-compartment with
    # first-order absorption and an absorption lag time.
    # Frederiksen 2023 Table 2. Every parameter in this block except the
    # fed absorption rate constant was FIXED to the final estimate of an
    # initial brexpiprazole-only population PK model, to stabilise the
    # joint parent + metabolite fit (Results: "To stabilize the model,
    # parameters related to the PK of brexpiprazole (ka, ALAG, V2, V3, Q,
    # IIV_ka, and IIV_V3) were fixed to the final estimates from the
    # initial PopPK model based on brexpiprazole data only").
    # All clearances and volumes are APPARENT (CL/F, V/F): every study
    # dosed orally and no absolute bioavailability was estimated.
    # ------------------------------------------------------------------

    lka <- fixed(log(1.02))
    label("Brexpiprazole absorption rate constant, fasted (1/h)")
    # Table 2, "Absorption rate constant fasted" = 1.02 (fixed).

    e_fed_ka <- -0.47549
    label("Fractional change in ka when dosed fed (unitless)")
    # Derived from Table 2: ka fed = 0.535 (23.4 %RSE) and ka fasted =
    # 1.02 (fixed), so the linear-deviation coefficient is
    # 0.535 / 1.02 - 1 = -0.4754902, i.e. a 47.5% lower absorption rate
    # constant in the fed state. Not itself printed in the paper; it is
    # an exact re-parameterisation of the two printed ka values. Left
    # estimated (not fixed) because the fed ka was estimated.

    ltlag <- fixed(log(0.411))
    label("Brexpiprazole absorption lag time (h)")
    # Table 2, "Lag-time" = 0.411 (fixed).

    lvc <- fixed(log(81.2))
    label("Brexpiprazole apparent central volume V2 (L)")
    # Table 2, "Volume of distribution, brexpiprazole central
    # compartment (V2)" = 81.2 (fixed).

    lvp <- fixed(log(40.1))
    label("Brexpiprazole apparent peripheral volume V3 (L)")
    # Table 2, "Volume of distribution, brexpiprazole peripheral
    # compartment (V3)" = 40.1 (fixed).

    lq <- fixed(log(0.714))
    label("Brexpiprazole apparent inter-compartmental clearance Q (L/h)")
    # Table 2, "Inter-compartmental clearance, brexpiprazole (Q)" =
    # 0.714 (fixed).

    # ------------------------------------------------------------------
    # Formation clearances out of the brexpiprazole central compartment.
    # These two ARE the parent's total elimination: the model carries no
    # separate non-formation clearance arm, so
    # CL_brexpiprazole_total = cl_form_dm3412 + cl_form_dm3411
    #                        = 0.0224 + 1.12 = 1.1424 L/h, and the
    # fraction of parent clearance forming DM-3412 is 0.0196.
    # ------------------------------------------------------------------

    lcl_form_dm3412 <- log(0.0224)
    label("Brexpiprazole-to-DM-3412 apparent formation clearance CLM1 (L/h)")
    # Table 2, "Clearance from brexpiprazole to DM-3412 (CLM1)" =
    # 0.0224 (22.7 %RSE). The Results note that the Table 2 bootstrap
    # 95% CI for this and the other DM-3412 parameters does NOT contain
    # the point estimate, driven by the <LLOQ DM-3412 observations; the
    # point estimate is used here, as it is the final model's estimate.

    lcl_form_dm3411 <- log(1.12)
    label("Brexpiprazole-to-DM-3411 apparent formation clearance CLM2 (L/h)")
    # Table 2, "Clearance from brexpiprazole to DM-3411 (CLM2)" =
    # 1.12 (2.7 %RSE).

    # ------------------------------------------------------------------
    # DM-3412 (minor, CYP2D6-formed metabolite) disposition:
    # two-compartment. Table 2.
    # ------------------------------------------------------------------

    lvc_dm3412 <- log(0.447)
    label("DM-3412 apparent central volume V4 (L)")
    # Table 2, "Volume of distribution, DM-3412 central compartment
    # (V4)" = 0.447 (41.4 %RSE).

    lvp_dm3412 <- log(20.4)
    label("DM-3412 apparent peripheral volume V5 (L)")
    # Table 2, "Volume of distribution, DM-3412 peripheral compartment
    # (V5)" = 20.4 (22.0 %RSE).

    lq_dm3412 <- log(2.16)
    label("DM-3412 apparent inter-compartmental clearance QMET (L/h)")
    # Table 2, "Inter-compartmental clearance, DM-3412 (QMET)" =
    # 2.16 (21.9 %RSE).

    lcl_dm3412 <- log(0.420)
    label("DM-3412 apparent elimination clearance CLMET1 (L/h)")
    # Table 2, "Clearance from DM-3412 central compartment (CLMET1)" =
    # 0.420 (21.8 %RSE).

    # ------------------------------------------------------------------
    # DM-3411 (major metabolite) disposition: one-compartment. Table 2.
    # ------------------------------------------------------------------

    lvc_dm3411 <- log(2.44)
    label("DM-3411 apparent central volume V6 (L)")
    # Table 2, "Volume of distribution, DM-3411 central compartment
    # (V6)" = 2.44 (8.0 %RSE). This is the typical value at the
    # reference BMI of 26 kg/m^2.

    lcl_dm3411 <- log(3.33)
    label("DM-3411 apparent elimination clearance CLMET2 (L/h)")
    # Table 2, "Clearance from DM-3411 central compartment (CLMET2)" =
    # 3.33 (3.5 %RSE).

    e_bmi_vc_dm3411 <- -0.0107
    label("Fractional change in DM-3411 V6 per 1 kg/m^2 above BMI 26 (1/(kg/m^2))")
    # Table 2, "Body mass index on V6" = -0.0107 (95.3 %RSE). The
    # covariate equation and centring value are not printed; see the
    # BMI entry in covariateData for the form and the choice of 26
    # kg/m^2 (the Table 1 population-PK-set median BMI).

    # ------------------------------------------------------------------
    # Interindividual variability. Exponential model on each parameter.
    # Table 2 reports IIV as a percent coefficient of variation with the
    # footnote %CV = sqrt(exp(omega^2) - 1) * 100, so each variance below
    # is omega^2 = log((CV/100)^2 + 1).
    #
    # The paper states three IIV covariance blocks were estimated
    # (V2-CLM1-CLM2, V4-CLMET1, and V6-CLMET2) but reports NO covariance
    # or correlation values anywhere in the article, Table 2, or the
    # Supplementary Material (Tables S1-S3). The variances are therefore
    # encoded as a diagonal Omega; see the vignette Errata. Simulated
    # between-subject correlation among these parameters will be absent.
    # ------------------------------------------------------------------

    etalka ~ fixed(1.72678)
    # Table 2, IIV on ka = 215% (fixed); log(2.15^2 + 1) = 1.72678.

    etalvc ~ 0.23688
    # Table 2, IIV on V2 = 51.7% (4.3 %RSE); log(0.517^2 + 1) = 0.23688.

    etalvp ~ fixed(0.44917)
    # Table 2, IIV on V3 = 75.3% (fixed); log(0.753^2 + 1) = 0.44917.

    etalcl_form_dm3412 ~ 0.52801
    # Table 2, IIV on CLM1 = 83.4% (9.3 %RSE);
    # log(0.834^2 + 1) = 0.52801.

    etalcl_form_dm3411 ~ 0.30396
    # Table 2, IIV on CLM2 = 59.6% (5.8 %RSE);
    # log(0.596^2 + 1) = 0.30396.

    etalvc_dm3412 ~ 1.70373
    # Table 2, IIV on V4 = 212% (33.7 %RSE); log(2.12^2 + 1) = 1.70373.

    etalcl_dm3412 ~ 0.21916
    # Table 2, IIV on CLMET1 = 49.5% (14.6 %RSE);
    # log(0.495^2 + 1) = 0.21916.

    etalvc_dm3411 ~ 0.95073
    # Table 2, IIV on V6 = 126% (7.6 %RSE); log(1.26^2 + 1) = 0.95073.

    etalcl_dm3411 ~ 0.29259
    # Table 2, IIV on CLMET2 = 58.3% (6.7 %RSE);
    # log(0.583^2 + 1) = 0.29259.

    # ------------------------------------------------------------------
    # Residual unexplained variability: a proportional error model for
    # each of the three analytes (Methods: "The residual unexplained
    # variability (RUV) was modeled separately for each of the three
    # compounds"; Results: "The RUV was modeled using proportional error
    # models for all three compounds").
    #
    # Table 2 prints the three residual-error rows as bare numbers with
    # no scale label. They are read here as NONMEM $SIGMA VARIANCES, so
    # the nlmixr2 proportional SD is their square root. The paper
    # converts the OMEGA values to %CV under an explicit footnote but
    # gives the SIGMA rows no such footnote, i.e. it leaves them on the
    # native NONMEM variance scale; and the resulting 20-31% residual CVs
    # are the expected magnitude for pooled rich + sparse phase I/II
    # data, whereas reading them as SDs would imply an implausible
    # 4-9% residual CV. See the vignette Errata.
    # ------------------------------------------------------------------

    propSd <- 0.30610
    label("Brexpiprazole proportional residual SD (fraction)")
    # Table 2, "Residual error brexpiprazole" = 0.0937 (variance);
    # sqrt(0.0937) = 0.30610.

    propSd_dm3412 <- 0.20543
    label("DM-3412 proportional residual SD (fraction)")
    # Table 2, "Residual error DM-3412" = 0.0422 (variance);
    # sqrt(0.0422) = 0.20543.

    propSd_dm3411 <- 0.29120
    label("DM-3411 proportional residual SD (fraction)")
    # Table 2, "Residual error DM-3411" = 0.0848 (variance);
    # sqrt(0.0848) = 0.29120.
  })

  model({
    # 1. Individual parameters. Exponential IIV; the food effect is a
    #    linear-deviation multiplier on ka and the BMI effect a
    #    linear-deviation multiplier on the DM-3411 central volume.
    ka   <- exp(lka + etalka) * (1 + e_fed_ka * FED)
    tlag <- exp(ltlag)
    vc   <- exp(lvc + etalvc)
    vp   <- exp(lvp + etalvp)
    q    <- exp(lq)

    cl_form_dm3412 <- exp(lcl_form_dm3412 + etalcl_form_dm3412)
    cl_form_dm3411 <- exp(lcl_form_dm3411 + etalcl_form_dm3411)

    vc_dm3412 <- exp(lvc_dm3412 + etalvc_dm3412)
    vp_dm3412 <- exp(lvp_dm3412)
    q_dm3412  <- exp(lq_dm3412)
    cl_dm3412 <- exp(lcl_dm3412 + etalcl_dm3412)

    vc_dm3411 <- exp(lvc_dm3411 + etalvc_dm3411) *
      (1 + e_bmi_vc_dm3411 * (BMI - 26))
    cl_dm3411 <- exp(lcl_dm3411 + etalcl_dm3411)

    # 2. ODE system (Frederiksen 2023 Figure 1b). Brexpiprazole leaves
    #    the central compartment only via the two formation clearances,
    #    so no separate parent elimination term appears.
    d/dt(depot) <- -ka * depot

    d/dt(central) <-
      ka * depot +
      q * peripheral1 / vp -
      (q + cl_form_dm3412 + cl_form_dm3411) * central / vc

    d/dt(peripheral1) <-
      q * central / vc - q * peripheral1 / vp

    # DM-3412: minor metabolite, formed almost exclusively by CYP2D6,
    # two-compartment disposition.
    d/dt(central_dm3412) <-
      cl_form_dm3412 * central / vc +
      q_dm3412 * peripheral1_dm3412 / vp_dm3412 -
      (q_dm3412 + cl_dm3412) * central_dm3412 / vc_dm3412

    d/dt(peripheral1_dm3412) <-
      q_dm3412 * central_dm3412 / vc_dm3412 -
      q_dm3412 * peripheral1_dm3412 / vp_dm3412

    # DM-3411: major metabolite, formed mainly by CYP3A4 with a CYP2D6
    # contribution, one-compartment disposition.
    d/dt(central_dm3411) <-
      cl_form_dm3411 * central / vc -
      cl_dm3411 * central_dm3411 / vc_dm3411

    # 3. Absorption lag.
    alag(depot) <- tlag

    # 4. Observations. Dose in mg over an apparent volume in L gives
    #    mg/L; the assay reports ng/mL, so multiply by 1000
    #    (1 mg/L = 1000 ng/mL). The metabolite states carry
    #    brexpiprazole-mass equivalents: the formation clearances are
    #    apparent, so the molar-mass ratio and the metabolite
    #    bioavailability are absorbed into V4, V5, and V6 exactly as they
    #    were in the NONMEM fit.
    Cc        <- 1000 * central        / vc
    Cc_dm3412 <- 1000 * central_dm3412 / vc_dm3412
    Cc_dm3411 <- 1000 * central_dm3411 / vc_dm3411

    Cc        ~ prop(propSd)
    Cc_dm3412 ~ prop(propSd_dm3412)
    Cc_dm3411 ~ prop(propSd_dm3411)
  })
}
