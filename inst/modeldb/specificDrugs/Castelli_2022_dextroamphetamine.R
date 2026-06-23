Castelli_2022_dextroamphetamine <- function() {
  description <- "One-compartment population PK model for dextroamphetamine transdermal system (d-ATS) in adults and children with ADHD (Castelli 2022 APNA poster), with sequential zero- and first-order absorption (zero-order release over duration D1 into depot, then first-order Ka into central), power-law body-weight scaling on CL/F (exponent 0.47), V/F (0.53), and Ka (-0.29) at an assumed 70 kg reference, independent IIV on CL/F, V/F, Ka, and D1, bioavailability anchored at F = 1, and residual error not reported in the conference poster (encoded fixed at 0; see vignette Errata)."
  reference <- paste(
    "Castelli M, Suzuki K, Starling B, Balakrishnan K, Meeves S, Komaroff M,",
    "Lennie J, Mondick JT, Faraone SV.",
    "Extrapolation of the Efficacy of a Dextroamphetamine Transdermal System",
    "Investigated in Pediatric Populations to Adults Using Pharmacokinetic Modeling.",
    "Poster, American Psychiatric Nurses Association (APNA), 14 June 2022.",
    "https://metrumrg.com/wp-content/uploads/2022/10/CastelliC_APNA-PK-Ped-Adult-Poster-2022-06-14-1.pdf",
    sep = " "
  )
  vignette <- "Castelli_2022_dextroamphetamine"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Castelli 2022 Methods: 'Body weight was predefined as a covariate on CL/F, V/F, and Ka' (full-covariate-modeling approach). Table 1 reports the weight-effect estimates 0.47 on CL/F, 0.53 on V/F, and -0.29 on Ka, interpreted here as power-law exponents of the form parameter_i = TVparameter * (WT_i / WT_ref)^exponent. The reference weight WT_ref is NOT stated in the poster; encoded here as 70 kg per the rounded-standard policy for an unreported adult-pediatric pooled reference (study-population median ~66 kg; adult median 73.2 kg, pediatric median 40.6 kg). See vignette Errata.",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 156L,
    n_studies      = 6L,
    n_observations = 6607L,
    age_range      = "6-62 years (adults 18-62, children 6-12)",
    age_median     = "adults 33 years; children 10 years",
    weight_range   = "23.1-101 kg (adults 43.8-101, children 23.1-63.5)",
    weight_median  = "adults 73.2 kg; children 40.6 kg",
    sex_female_pct = NA_real_,
    disease_state  = "Attention-deficit/hyperactivity disorder (ADHD); a pooled population of 122 adults (78%; median age 33 years, range 18-62) and 34 children (22%; median age 10 years, range 6-12) with ADHD across six PK studies of the dextroamphetamine transdermal system (d-ATS; Xelstrym). Sex distribution and race/ethnicity are not tabulated in the conference poster.",
    dose_range     = "5 mg, 10 mg, 15 mg, and 20 mg d-ATS transdermal patches; the 5/10/15/20 mg doses were the levels evaluated in the pivotal pediatric study and were the reference dose levels for the exposure-matching extrapolation to adults reported in the poster. Patches were applied to the hip for 9 hours per day (standard d-ATS wear time per the FDA label) in the pivotal pediatric study; the poster does not state the wear time used in each contributing PK study.",
    regions        = "Not stated in the poster.",
    blq_handling   = "Below-limit-of-quantification (BLQ) concentrations in the elimination phase were excluded from parameter estimation; BLQ records in the absorption/distribution phase were retained.",
    notes          = "Castelli 2022 Methods + Results: 156 subjects pooled from six PK studies of d-ATS, totalling 6607 amphetamine plasma concentration data points. NONMEM 7.4.3 with first-order conditional estimation. A full-covariate-modeling approach was taken with body weight pre-specified on CL/F, V/F, and Ka. Inter-occasion variability was estimated on D1 and bioavailability (F); the IOV magnitudes are NOT reported in the poster and are NOT encoded here (see vignette Errata). The residual-error model and magnitudes are NOT reported in the poster either."
  )

  ini({
    # Structural PK parameters -- Castelli 2022 Table 1 final-model estimates
    # (with 95% CIs). Reference subject: 70 kg adult (the WT reference is not
    # stated in the poster and is encoded here as the rounded-standard 70 kg
    # per the unreported-reference-value policy; see vignette Errata).
    lcl <- log(18.4)  ; label("Apparent clearance CL/F at WT 70 kg (L/h)")             # Castelli 2022 Table 1 final CL/F = 18.4 L/h (95% CI 17.6, 19.2)
    lvc <- log(51.9)  ; label("Apparent volume of distribution V/F at WT 70 kg (L)")    # Castelli 2022 Table 1 final V/F = 51.9 L (95% CI 48.1, 56.1)
    lka <- log(0.070) ; label("First-order absorption rate constant Ka at WT 70 kg (1/h)") # Castelli 2022 Table 1 final Ka = 0.070 1/h (95% CI 0.067, 0.073)
    ld1 <- log(1.9)   ; label("Zero-order absorption duration D1 (h)")                 # Castelli 2022 Table 1 final D1 = 1.9 h (95% CI 1.8, 2.1)

    # Bioavailability anchor. d-ATS is a transdermal patch; absolute F is not
    # identifiable from oral/transdermal-only data and the poster does not
    # report a typical-value F (only the existence of IOV on F). Anchor at 1
    # by the standard popPK convention for an apparent-parameter
    # parameterisation (CL/F, V/F).
    lfdepot <- fixed(log(1.0)) ; label("Bioavailability into depot (F, fixed at 1)")    # Castelli 2022: typical-value F not reported (Methods: 'Inter-occasion variability was estimated for D1 and bioavailability (F)'); anchored at 1 per the apparent-parameter convention -- see vignette Errata

    # Body-weight covariate effects. Castelli 2022 Methods describes a
    # full-covariate modeling approach with body weight pre-specified on
    # CL/F, V/F, and Ka; Table 1 reports the "weight effect estimate" with
    # 95% CIs and no functional form. Encoded as power-law exponents:
    #   CL/F_i = 18.4 * (WT/70)^0.47
    #   V/F_i  = 51.9 * (WT/70)^0.53
    #   Ka_i   = 0.070 * (WT/70)^-0.29
    # This is the most common interpretation of a "weight effect estimate"
    # in a full-covariate popPK model; see vignette Errata.
    e_wt_cl <-  0.47  ; label("Body-weight power exponent on CL/F (WT/70 kg; unitless)") # Castelli 2022 Table 1 weight effect on CL/F = 0.47 (95% CI 0.36, 0.58)
    e_wt_vc <-  0.53  ; label("Body-weight power exponent on V/F (WT/70 kg; unitless)")  # Castelli 2022 Table 1 weight effect on V/F = 0.53 (95% CI 0.31, 0.75)
    e_wt_ka <- -0.29  ; label("Body-weight power exponent on Ka (WT/70 kg; unitless)")    # Castelli 2022 Table 1 weight effect on Ka = -0.29 (95% CI -0.43, -0.15)

    # Inter-individual variability -- Castelli 2022 Table 1 reports IIV as
    # CV%. Convert to log-scale variance via omega^2 = log(CV^2 + 1):
    #   CL/F: 20.1% CV -> log(1 + 0.201^2) = 0.0396
    #   V/F : 37.4% CV -> log(1 + 0.374^2) = 0.1309
    #   Ka  : 20.3% CV -> log(1 + 0.203^2) = 0.0404
    #   D1  : 48.6% CV -> log(1 + 0.486^2) = 0.2121
    # The poster does not state whether the etas are correlated; encoded as
    # independent (the simplest assumption consistent with what is reported).
    etalcl ~ log(1 + 0.201^2) # Castelli 2022 Table 1: IIV CL/F 20.1% CV (95% CI 16.6, 23.1)
    etalvc ~ log(1 + 0.374^2) # Castelli 2022 Table 1: IIV V/F  37.4% CV (95% CI 30.3, 43.7)
    etalka ~ log(1 + 0.203^2) # Castelli 2022 Table 1: IIV Ka   20.3% CV (95% CI 15.7, 24.1)
    etald1 ~ log(1 + 0.486^2) # Castelli 2022 Table 1: IIV D1   48.6% CV (95% CI 40.6, 55.8)

    # Residual unexplained variability. The poster does NOT report the
    # residual-error model or its magnitude. Encoded as fixed(0) per the
    # standing "unreported RUV -> fixed(0) + Errata" policy so the model
    # remains usable for typical-value / IIV-only forward simulation; users
    # wanting a stochastic VPC must supply their own residual magnitude.
    # See vignette Errata.
    propSd <- fixed(0) ; label("Proportional residual SD (fraction; FIXED at 0 -- not reported in poster)") # Castelli 2022: residual error model not reported in the conference poster
  })

  model({
    # Individual PK parameters with body-weight power-law scaling on CL/F,
    # V/F, and Ka. Reference subject: 70 kg adult (see Errata for the
    # rationale behind 70 kg). D1 carries IIV but no weight effect per
    # Castelli 2022 Table 1.
    cl <- exp(lcl + etalcl) * (WT / 70) ^ e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70) ^ e_wt_vc
    ka <- exp(lka + etalka) * (WT / 70) ^ e_wt_ka
    d1 <- exp(ld1 + etald1)

    kel <- cl / vc

    # One-compartment PK with sequential zero- and first-order absorption:
    # the d-ATS dose enters the depot via a zero-order window of duration
    # D1 (representing the initial patch release into the skin reservoir),
    # then is absorbed first-order at rate Ka into central. Bioavailability
    # F is anchored at 1.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    dur(depot) <- d1
    f(depot)   <- exp(lfdepot)

    # Plasma amphetamine concentration. Dose units mg, vc L -> internal
    # mg/L; multiply by 1000 to report in ng/mL, the standard plasma unit
    # for amphetamine PK.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
