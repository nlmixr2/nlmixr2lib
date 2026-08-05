Gong_2023_pemigatinib <- function() {
  description <- paste(
    "Two-compartment population PK model for oral pemigatinib (a selective",
    "fibroblast growth factor receptor [FGFR] 1-3 inhibitor) in healthy",
    "participants and patients with advanced solid tumors including",
    "cholangiocarcinoma (CCA) and myeloid/lymphoid neoplasms with FGFR1",
    "rearrangement (Gong 2023; N = 467 across seven FIGHT studies,",
    "4552 concentration records). Disposition is a two-compartment model",
    "with sequential zero-order and first-order oral absorption (dose",
    "enters the depot over a zero-order window of duration D1 = 0.810 h,",
    "then is absorbed first-order at ka = 3.67 h^-1) and linear",
    "elimination (CL/F = 10.7 L/h, Vc/F = 118 L, Vp/F = 95.0 L,",
    "Q/F = 25.2 L/h). Covariate effects retained after stepwise forward",
    "addition and backward elimination: concomitant phosphate-binding",
    "agent use reduces CL/F by 15.5%; male sex increases CL/F by 26.2%",
    "and reduces ka by 58.3% (female is the reference); concomitant",
    "proton-pump inhibitor use reduces ka by 62.0%; concomitant",
    "histamine-2 receptor antagonist use reduces D1 by 33.4%; and",
    "baseline body weight scales Vc/F and Vp/F as power functions of",
    "(WT / 73.9 kg) with exponents 0.842 and 1.13 respectively. The",
    "paper reports all of these covariate effects as statistically",
    "significant but not clinically significant on steady-state",
    "exposure. Inter-individual variability is a correlated CL/F-Vc/F",
    "block (log-scale variances 0.201 and 0.191, covariance 0.129,",
    "correlation 0.659) plus independent terms on ka (1.69) and D1",
    "(0.371); IIV on Q/F and Vp/F was fixed to zero in the final model.",
    "Residual error is additive on the log scale (approximately",
    "proportional in linear space) at 32.3%. Concentrations are",
    "expressed in nM using the pemigatinib molecular weight of",
    "487.5 g/mol, matching the h*nM exposure units used throughout the",
    "paper's exposure-response analyses.",
    sep = " "
  )
  reference <- paste(
    "Gong X, Akil A, Ndi A, Ji T, Liu X, Lovern M, Chen X. (2023).",
    "Population pharmacokinetic and exposure-response analyses of",
    "pemigatinib in patients with advanced solid tumors including",
    "cholangiocarcinoma.",
    "CPT Pharmacometrics Syst Pharmacol 12(11):1784-1794.",
    "doi:10.1002/psp4.13064",
    sep = " "
  )
  vignette <- "Gong_2023_pemigatinib"
  units    <- list(time = "h", dosing = "mg", concentration = "nM")

  # States hold pemigatinib amounts in mg (the dosing unit); the observation Cc
  # converts the central amount to nM via the Appendix S1 scaling
  # S2 = V2/1000000*487.5, i.e. pemigatinib MW 487.5 g/mol.
  compartmentData <- list(
    depot       = list(analyte = "pemigatinib", units = "mg", specimen = "administration site", verified = FALSE),
    central     = list(analyte = "pemigatinib", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "pemigatinib", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed baseline body weight. Cohort mean 75.7 kg (SD 19.7),",
        "median 73.9 kg, range 39.8-156 kg across 467 subjects (Gong 2023",
        "Table 1). Enters as a power function normalized to the cohort",
        "median 73.9 kg on both apparent volumes, per the NONMEM control",
        "stream in Appendix S1:",
        "V2BWT = ((BWT/73.9)**THETA(12)) and V3BWT = ((BWT/73.9)**THETA(13)).",
        "The 73.9 kg reference is the population median from Table 1, not a",
        "rounded standard weight. Gong 2023 Table 2 lists the two exponents",
        "in rows labelled 'Body weight on Vc/F, %' and 'Body weight to",
        "Vp/F, %'; the '%' in those row labels is a table-formatting",
        "artifact inherited from the categorical covariate rows above them",
        "-- the values 0.842 and 1.13 are dimensionless power exponents,",
        "as the control-stream `**THETA` form confirms. No weight effect",
        "was retained on CL/F or Q/F."
      ),
      source_name        = "BWT"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female) -- the reference (most common) level in Gong 2023",
      notes              = paste(
        "1 = female, 0 = male. 219 of 467 subjects (46.9%) were men",
        "(Gong 2023 Table 1), so female is the more common level and is",
        "the reference category. The source column is SEXN, coded 1 = male",
        "and 2 = female; the NONMEM control stream in Appendix S1 sets",
        "CLSEXN = 1 and KASEXN = 1 for SEXN.EQ.2 ('Most common') and",
        "applies (1 + THETA) for SEXN.EQ.1. Canonical SEXF is therefore",
        "the inversion SEXF = 2 - SEXN, and the male deviation is applied",
        "in model() via the (1 - SEXF) indicator. Both retained sex effects",
        "are reported in Gong 2023 Table 2 as 'Male sex on CL/F, %' (+26.2%)",
        "and 'Male sex on Ka, %' (-58.3%), confirming that SEXN = 1 is male."
      ),
      source_name        = "SEXN"
    ),
    CONMED_PPI = list(
      description        = "Concomitant proton-pump inhibitor use",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant PPI use)",
      notes              = paste(
        "1 = patient on a concomitant proton-pump inhibitor, 0 = not.",
        "114 of 467 subjects (24.4%) were PPI users (Gong 2023 Table 1).",
        "Applied multiplicatively to ka as (1 + e_ppi_ka * CONMED_PPI) per",
        "the Appendix S1 control stream:",
        "IF(PPI.EQ.0) KAPPI = 1; IF(PPI.EQ.1) KAPPI = (1 + THETA(10)).",
        "Mechanistically consistent with pemigatinib's pH-dependent",
        "solubility (poorly soluble at pH > 2; Gong 2023 Discussion)."
      ),
      source_name        = "PPI"
    ),
    CONMED_H2RA = list(
      description        = "Concomitant histamine-2 receptor antagonist use",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant H2RA use)",
      notes              = paste(
        "1 = patient on a concomitant histamine-2 receptor antagonist",
        "(H2B in the source), 0 = not. 39 of 467 subjects (8.4%) were H2RA",
        "users (Gong 2023 Table 1). Applied multiplicatively to the",
        "zero-order input duration D1 as (1 + e_h2ra_d1 * CONMED_H2RA) per",
        "the Appendix S1 control stream:",
        "IF(H2B.EQ.0) D1H2B = 1; IF(H2B.EQ.1) D1H2B = (1 + THETA(9)).",
        "This is the least precisely estimated effect in the final model",
        "(%RSE 38.9; 95% CI -58.9%, -8.0%) and Gong 2023 Discussion notes",
        "the H2RA result is confounded by sparse sampling and the small",
        "number of H2RA users."
      ),
      source_name        = "H2B"
    ),
    CONMED_PHOSBINDER = list(
      description        = "Concomitant phosphate-binding agent use",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant phosphate-binder use)",
      notes              = paste(
        "1 = patient on a concomitant phosphate-binding agent, 0 = not.",
        "61 of 467 subjects (13.1%) were phosphate-binder users, all of",
        "them patients rather than healthy participants (Gong 2023",
        "Table 1). Phosphate binders are used to manage the",
        "hyperphosphatemia that is the on-target class effect of FGFR",
        "inhibition. Applied multiplicatively to CL/F as",
        "(1 + e_phosbinder_cl * CONMED_PHOSBINDER) per the Appendix S1",
        "control stream:",
        "IF(BINDER.EQ.0) CLBINDER = 1; IF(BINDER.EQ.1) CLBINDER = (1 + THETA(7)).",
        "Gong 2023 Discussion states the 15.5% CL/F reduction is 'not",
        "readily explained' and is not clinically significant."
      ),
      source_name        = "BINDER"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened on CL/F and Vc/F in the stepwise covariate search but",
        "not retained in the final model (Gong 2023 Methods 'Population PK",
        "analysis'). Cohort mean 54.3 years (SD 14.5), median 56.0,",
        "range 19-83 (Table 1)."
      )
    ),
    TUMTP_OTHER = list(
      description = "Tumor type / participant type indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Participant type (healthy participant; patient with CCA, other",
        "solid tumors, or MLN with FGFR1 rearrangement) was screened on",
        "CL/F and was significant at the forward-addition step, but was",
        "removed during final model refinement 'because the 95% CI for",
        "both parameters included 0' (Gong 2023 Results). Not in the",
        "final model."
      )
    ),
    CRCL = list(
      description = "Renal function (eGFR by MDRD), used for renal-impairment classification",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = paste(
        "Renal impairment classified by MDRD eGFR was screened on CL/F and",
        "was not significant (Gong 2023 Results). Renal clearance of",
        "pemigatinib is only 1-2% of total clearance (Discussion), so mild",
        "or moderate renal impairment minimally affects PK. Not in the",
        "final model."
      )
    ),
    HEPATIC_NCI = list(
      description = "Hepatic impairment category (NCI Organ Dysfunction Working Group)",
      units       = "(category)",
      type        = "categorical",
      notes       = paste(
        "Screened on CL/F and not significant (Gong 2023 Results). Only 12",
        "patients had moderate hepatic impairment, which prevented",
        "comparison of CL/F in that group. Not in the final model."
      )
    ),
    CONMED_CYP3A4_IND = list(
      description = "Concomitant CYP3A4 inducer use",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened on CL/F and not significant (Gong 2023 Results). Only 36",
        "subjects were on weak and 1 on a moderate CYP3A4 inducer",
        "(Table 1); the paper states there was insufficient data to support",
        "inducer-based dose adjustment. Not in the final model."
      )
    ),
    CONMED_CYP3A4_INH = list(
      description = "Concomitant CYP3A4 inhibitor use",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened on CL/F and not significant (Gong 2023 Results). 98",
        "subjects were on weak and 14 on a moderate CYP3A4 inhibitor",
        "(Table 1). The paper retains the label recommendation to adjust",
        "dose with moderate/strong CYP3A4 inhibitors on the basis of a",
        "dedicated DDI study and prior PBPK work, not this popPK analysis.",
        "Not in the final model."
      )
    ),
    FASTED = list(
      description = "Fasting status at the dose record",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened on ka (Gong 2023 Methods) and not retained in the final",
        "model. Source column FASTFL / FOOD in the Appendix S1 $INPUT",
        "record."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 467,
    n_studies      = 7,
    n_observations = 4552,
    age_range      = "19-83 years",
    age_median     = "56 years",
    age_mean_sd    = "54.3 (14.5) years",
    weight_range   = "39.8-156 kg",
    weight_median  = "73.9 kg",
    weight_mean_sd = "75.7 (19.7) kg",
    sex_female_pct = 53.1,
    disease_state  = paste(
      "Pooled healthy participants (78; 16.7%) and patients with advanced",
      "malignancies: cholangiocarcinoma (163; 34.9%), myeloid/lymphoid",
      "neoplasms with FGFR1 rearrangement (34; 7.3%), and other cancers",
      "(192; 41.1%)."
    ),
    dose_range     = paste(
      "1-20 mg orally once daily; continuous once-daily and intermittent",
      "(2 weeks on / 1 week off, 21-day cycle) regimens. The approved and",
      "recommended dose is 13.5 mg once daily on the intermittent schedule."
    ),
    studies        = paste(
      "FIGHT-101 (NCT02393248), FIGHT-102 (NCT03235570), FIGHT-104,",
      "FIGHT-105, FIGHT-106, FIGHT-202, FIGHT-203 (Gong 2023 Table S1)."
    ),
    regions        = "United States, Denmark, Japan",
    renal_function = paste(
      "Unimpaired 250 (53.5%), mild 169 (36.2%), moderate 48 (10.3%);",
      "classified by MDRD eGFR."
    ),
    hepatic_function = paste(
      "Unimpaired 339 (72.6%), mild 116 (24.8%), moderate 12 (2.6%);",
      "classified by the NCI Hepatic Dysfunction Working Group criteria."
    ),
    co_medication  = paste(
      "Phosphate binders 61 (13.1%), PPI 114 (24.4%), H2RA 39 (8.4%);",
      "89.3% received pemigatinib monotherapy."
    ),
    notes          = paste(
      "Baseline demographics and covariate summary are Gong 2023 Table 1.",
      "This model is a refinement of an earlier pemigatinib popPK model",
      "built on FIGHT-101 (n = 157), FIGHT-102 (n = 25), and FIGHT-202",
      "(n = 136); the update adds MLN-FGFR1 patients and replaces the",
      "first-order absorption model with a sequential zero-/first-order",
      "absorption model (objective function -872.76 -> -2262.23)."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural PK -- Gong 2023 Table 2 'Final parameter estimate /
    # Population Mean' column. Units are L, L/h, and h, matching the
    # paper. Reference subject for the typical values is a FEMALE with
    # baseline body weight 73.9 kg (the cohort median), not on a
    # phosphate binder, not on a PPI, and not on an H2RA -- i.e. every
    # covariate multiplier equals 1.
    #
    # The NONMEM control stream in Appendix S1 carries INITIAL estimates
    # (CL 10.41, V2 120.717, Q 23.7627, V3 97.9081, KA 3.49239,
    # D1 0.804979); the FINAL estimates below come from Table 2.
    # ---------------------------------------------------------------------
    lcl <- log(10.7)  ; label("Apparent clearance CL/F (L/h)")                        # Gong 2023 Table 2: CL/F 10.7 L/h, %RSE 2.7, 95% CI 10.1-11.3
    lvc <- log(118)   ; label("Apparent central volume of distribution Vc/F (L)")     # Gong 2023 Table 2: Vc/F 118 L, %RSE 3.4, 95% CI 110-125
    lvp <- log(95.0)  ; label("Apparent peripheral volume of distribution Vp/F (L)")  # Gong 2023 Table 2: Vp/F 95.0 L, %RSE 3.0, 95% CI 89.4-101
    lq  <- log(25.2)  ; label("Apparent inter-compartmental clearance Q/F (L/h)")     # Gong 2023 Table 2: Q/F 25.2 L/h, %RSE 7.7, 95% CI 21.4-29.0
    lka <- log(3.67)  ; label("First-order absorption rate constant ka (1/h)")        # Gong 2023 Table 2: Ka 3.67 1/h, %RSE 12.2, 95% CI 2.79-4.55
    ld1 <- log(0.810) ; label("Zero-order absorption input duration D1 (h)")          # Gong 2023 Table 2: D1 0.810 h, %RSE 3.0, 95% CI 0.763-0.857

    # ---------------------------------------------------------------------
    # Covariate effects. Gong 2023 Table 2 reports the five categorical
    # effects as percentage changes; the Appendix S1 control stream shows
    # each is a multiplicative (1 + THETA) term applied to the typical
    # value, so the percentages divide by 100 directly. The two body-weight
    # effects are power exponents on (WT / 73.9), NOT percentages -- see
    # the WT covariateData note above.
    # ---------------------------------------------------------------------
    e_phosbinder_cl <- -0.155 ; label("Fractional change in CL/F with concomitant phosphate-binder use (unitless)")   # Gong 2023 Table 2: Phosphate binder on CL/F -15.5%, %RSE 14.6, 95% CI -19.9 to -11.1
    e_male_cl       <-  0.262 ; label("Fractional change in CL/F for male vs female (unitless)")                      # Gong 2023 Table 2: Male sex on CL/F +26.2%, %RSE 18.5, 95% CI 16.7-35.8
    e_ppi_ka        <- -0.620 ; label("Fractional change in ka with concomitant PPI use (unitless)")                  # Gong 2023 Table 2: PPI on Ka -62.0%, %RSE 12.2, 95% CI -76.8 to -47.2
    e_male_ka       <- -0.583 ; label("Fractional change in ka for male vs female (unitless)")                        # Gong 2023 Table 2: Male sex on Ka -58.3%, %RSE 11.1, 95% CI -71.0 to -45.6
    e_h2ra_d1       <- -0.334 ; label("Fractional change in D1 with concomitant H2RA use (unitless)")                 # Gong 2023 Table 2: H2B antagonist on D1 -33.4%, %RSE 38.9, 95% CI -58.9 to -8.0
    e_wt_vc         <-  0.842 ; label("Body-weight power exponent on Vc/F (unitless, reference 73.9 kg)")             # Gong 2023 Table 2: Body weight on Vc/F 0.842, %RSE 14.3, 95% CI 0.605-1.08; Appendix S1 V2BWT = (BWT/73.9)**THETA(12)
    e_wt_vp         <-  1.13  ; label("Body-weight power exponent on Vp/F (unitless, reference 73.9 kg)")             # Gong 2023 Table 2: Body weight to Vp/F 1.13, %RSE 13.0, 95% CI 0.842-1.42; Appendix S1 V3BWT = (BWT/73.9)**THETA(13)

    # ---------------------------------------------------------------------
    # Inter-individual variability. Gong 2023 Table 2 reports IIV as
    # '%CV'. Those values are 100 * omega (the log-scale SD), NOT
    # 100 * sqrt(exp(omega^2) - 1). This is settled exactly by the
    # paper's own reported correlation coefficient: the Table 2 footnote
    # gives cov(CL, Vc) = 0.129 and the body of the table gives the
    # correlation as 0.659.
    #   omega^2 = CV^2         -> 0.129 / (0.448 * 0.437)          = 0.6589  (matches)
    #   omega^2 = log(1 + CV^2)-> 0.129 / sqrt(0.1829 * 0.1748)    = 0.7215  (does not)
    # The Appendix S1 $OMEGA initial estimates corroborate the same
    # reading (e.g. IIV KA initial 1.65573 -> sqrt = 1.287 = 128.7%,
    # against the 130% final; IIV D1 initial 0.429359 -> 0.655 = 65.5%,
    # against the 60.9% final).
    #
    # IIV on Q/F and Vp/F is 'Fixed to 0' in Table 2 (Appendix S1:
    # '$OMEGA 0 FIX ; IIV Q' and '0 FIX ; IIV V3'), so no etalq / etalvp
    # term is declared here -- q and vp are typical-value-only.
    # ---------------------------------------------------------------------
    etalcl + etalvc ~ c(0.200704,
                        0.129, 0.190969)                                                                             # Gong 2023 Table 2: IIV CL/F 44.8 %CV -> 0.448^2 = 0.200704; IIV Vc/F 43.7 %CV -> 0.437^2 = 0.190969; covariance 0.129 (Table 2 footnote a, %RSE 13.6) giving correlation 0.659
    etalka          ~ 1.69                                                                                           # Gong 2023 Table 2: IIV Ka 130 %CV -> 1.30^2 = 1.69 (95% CI 113-144)
    etald1          ~ 0.370881                                                                                       # Gong 2023 Table 2: IIV D1 60.9 %CV -> 0.609^2 = 0.370881 (95% CI 54.4-66.8)

    # ---------------------------------------------------------------------
    # Residual error. Appendix S1 $ERROR fits the log-transformed
    # prediction with a single additive term on the log scale
    # (IPRED = LOG(F); Y = IPRED + EPS(1)), which is equivalent to
    # proportional error in linear space. Table 2 reports RUV as 32.3%,
    # i.e. sqrt(sigma^2); the $SIGMA initial 0.102961 -> sqrt = 0.3209
    # corroborates the same SD-scale reading used for the IIV terms.
    # ---------------------------------------------------------------------
    propSd <- 0.323 ; label("Proportional residual SD (fraction; log-scale additive in Gong 2023)")                   # Gong 2023 Table 2: RUV 32.3%, %RSE 3.0, 95% CI 30.4-34.1; Appendix S1 $SIGMA initial 0.102961
  })

  model({
    # -------------------------------------------------------------------
    # 1. Covariate multipliers. Every effect is multiplicative on the
    #    typical value, exactly as written in the Appendix S1 $PK block:
    #      CLCOV = CLBINDER * CLSEXN
    #      KACOV = KAPPI    * KASEXN
    #      D1COV = D1H2B
    #      V2COV = V2BWT,  V3COV = V3BWT
    #    The sex terms fire for MALES. Canonical SEXF is 1 = female, and
    #    female is the reference level, so the male indicator is
    #    (1 - SEXF).
    # -------------------------------------------------------------------
    male   <- 1 - SEXF

    cl_cov <- (1 + e_phosbinder_cl * CONMED_PHOSBINDER) * (1 + e_male_cl * male)
    ka_cov <- (1 + e_ppi_ka        * CONMED_PPI)        * (1 + e_male_ka * male)
    d1_cov <- (1 + e_h2ra_d1       * CONMED_H2RA)

    # 2. Individual PK parameters. IIV is on CL/F, Vc/F, ka, and D1 only;
    #    Q/F and Vp/F are typical-value-only (IIV fixed to 0 in Table 2).
    cl <- exp(lcl + etalcl) * cl_cov
    vc <- exp(lvc + etalvc) * (WT / 73.9) ^ e_wt_vc
    vp <- exp(lvp)          * (WT / 73.9) ^ e_wt_vp
    q  <- exp(lq)
    ka <- exp(lka + etalka) * ka_cov
    d1 <- exp(ld1 + etald1) * d1_cov

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 4. Two-compartment ODE system with sequential zero- and first-order
    #    absorption (NONMEM ADVAN4 TRANS4 with a D1 zero-order input on
    #    the depot). The dose is released into `depot` at a constant rate
    #    over D1 hours, then absorbed first-order at ka into `central`.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. Zero-order input duration. Dose records must carry rate = -2 so
    #    that rxode2 applies the modeled duration rather than treating the
    #    dose as an instantaneous bolus.
    dur(depot) <- d1

    # 6. Observation. Doses are in mg and vc is in L, so central / vc is
    #    mg/L; the paper reports concentrations and exposures in nM and
    #    h*nM. The Appendix S1 control stream applies exactly this
    #    conversion through the NONMEM scaling parameter
    #      S2 = V2 / 1000000 * 487.5  ; nM
    #    where 487.5 g/mol is the pemigatinib molecular weight, so
    #      C[nM] = A[mg] / V2[L] * 1e6 / 487.5.
    Cc <- central / vc * 1e6 / 487.5

    # 7. Residual error: log-scale additive in NONMEM == proportional here.
    Cc ~ prop(propSd)
  })
}
