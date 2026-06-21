Kassir_2014_tacrolimus <- function() {
  description <- "Two-compartment population PK model with first-order absorption and an absorption lag time for twice-daily oral tacrolimus in paediatric liver transplant recipients (Kassir 2014). Apparent oral clearance CL/F and apparent inter-compartmental clearance Q2/F scale allometrically with body weight at a fixed exponent of 0.75 referenced to the cohort median weight of 20 kg; apparent central volume V1/F and apparent peripheral volume V2/F scale at a fixed exponent of 1.0 to the same 20 kg reference; the first-order absorption rate constant ka carries an allometric exponent of -0.25 per Anderson and Holford theory. Apparent peripheral volume V2/F was fixed to 290 L during estimation to stabilise the model (Kassir 2014 Table 4 footnote). Inter-individual variability is diagonal on CL/F, V1/F, and Q2/F (no IIV on ka, tlag, or V2/F). Residual error is a proportional model. No covariates beyond body weight were retained after stepwise covariate analysis -- age, sex, type of transplant, age of liver donor, time post-transplantation, liver function tests, albumin, renal function (serum creatinine and creatinine clearance), haematocrit, use of steroids, presence of clinically relevant CYP3A4 inhibitors, and drug formulation were all screened and dropped (Kassir 2014 Results 'Analysis of covariates and sources of variability')."
  reference   <- "Kassir N, Labbe L, Delaloye J-R, Mouksassi M-S, Lapeyraque A-L, Alvarez F, Lallier M, Beaunoyer M, Theoret Y, Litalien C. Population pharmacokinetics and Bayesian estimation of tacrolimus exposure in paediatric liver transplant recipients. Br J Clin Pharmacol. 2014;77(6):1051-1063. doi:10.1111/bcp.12276"
  vignette    <- "Kassir_2014_tacrolimus"
  units       <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Body weight at the time of the PK profile. Enters all five PK parameters via fixed allometric exponents referenced to the cohort median weight of 20 kg (Kassir 2014 Methods 'Population pharmacokinetic analysis' equation block citing references [20, 21], i.e., Anderson and Holford theory-based scaling): CL/F = 12.1 * (WT/20)^0.75; V1/F = 31.3 * (WT/20)^1; Q2/F = 30.7 * (WT/20)^0.75; V2/F = 290 * (WT/20)^1; ka = 0.342 * (WT/20)^(-0.25). The cohort weight range was 4.5-57.8 kg (median 20.4 kg, Table 3); model application outside this range is not directly supported by the data and should be considered extrapolation.",
      source_name        = "WT"
    )
  )

  population <- list(
    species               = "human",
    n_subjects            = 30L,
    n_observations        = 341L,
    n_pk_profiles         = 38L,
    n_studies             = 1L,
    age_range             = "0.4-18.4 years",
    age_median            = "7.3 years",
    weight_range          = "4.5-57.8 kg",
    weight_median         = "20.4 kg (also used as the 20 kg allometric reference)",
    sex_female_pct        = 43.3,
    race_ethnicity        = NULL,
    disease_state         = "Paediatric liver transplant recipients on tacrolimus-based immunosuppression. Underlying diagnoses included biliary atresia (n = 12), tyrosinaemia (n = 8), North American Indian childhood cirrhosis (n = 3), fulminant hepatitis (n = 2), Alagille syndrome (n = 2), histiocytosis X (n = 1), sclerosing cholangitis (n = 1), and autoimmune hepatitis (n = 1). Median time post-transplantation 2.5 months (range 0.5-188.2 months); no patients were studied during the first 2 weeks post-transplant.",
    dose_range            = "Oral tacrolimus twice daily as capsules (n = 15) or as a 0.5 mg/mL extemporaneously compounded oral suspension (n = 15). Dose was individualised by the transplant team to maintain a trough concentration target of 5-15 ng/mL.",
    regions               = "Canada (single-centre paediatric liver transplant cohort at Centre Hospitalier Universitaire Sainte-Justine, Montreal, Quebec).",
    transplant_type       = "Cut-down liver = 20 (66.7%); full liver = 10 (33.3%) (Kassir 2014 Table 2). Donor age range 0.58-66 years.",
    co_medications        = "CYP3A4 inhibitors (amlodipine, lansoprazole, or diltiazem) in n = 8 (26.7%); concomitant steroids (prednisone) in n = 20 (66.7%); 17% with hepatic impairment defined as total bilirubin > 68.4 mmol/L and/or ALT > 2x age-group upper limit (Kassir 2014 Table 2).",
    sampling_window       = "Twelve-hour intensive PK profiles obtained at steady state (>= 3 days on the same dose), with samples at 0 (trough), 0.5, 1, 1.5, 2, 3, 4, 8, and 12 hours post-dose. Median 9 concentrations per profile (range 8-10). Total 341 samples across 38 profiles in 30 patients.",
    assay                 = "Tacrolimus whole-blood concentrations determined by microparticle enzyme immunoassay (MEIA IMx, Abbott Laboratories, Abbott Park, IL, USA). Limits of detection 1.5-30 ng/mL; between-run CVs 14.10%, 11.15%, and 10.21% at 5, 11, and 22 ng/mL.",
    cyp3a5_distribution   = "CYP3A5 genotype not available for recipients or donors in this study; not tested as a covariate.",
    notes                 = "Retrospective analysis of intensive PK profiles obtained between July 2006 and May 2011. AUC(0-12) measurement was triggered by clinical indications (nephrotoxicity despite in-range troughs; high intra-individual trough variability; or initiation of mycophenolate mofetil). The published model is intended for paediatric liver transplant recipients beyond the first 2 weeks post-transplant, primarily within the first year. Application to young infants (< 1 year) is cautioned by the authors because only 4 such patients were included; CYP3A5 expresser status and donor age effects could not be evaluated. A typical 20 kg child has CL/F = 12.1 L/h, V1/F = 31.3 L, Q2/F = 30.7 L/h, V2/F = 290 L (sum Vd/F = 321.3 L), ka = 0.342 /h, tlag = 0.433 h."
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Subject age in years",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested in the stepwise covariate analysis (Kassir 2014 Methods 'Covariate analysis and sources of variability') but did not significantly explain variability in PK parameters after weight-based allometric scaling and was not retained in the final model. Cohort median 7.3 years, range 0.4-18.4 years (Table 3)."
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Tested in the stepwise covariate analysis but not retained in the final model (Kassir 2014 Results 'Analysis of covariates and sources of variability'). Cohort 43.3% female (Table 2)."
    ),
    TX_FULL = list(
      description        = "Whole-liver-graft indicator (vs cut-down liver graft)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (cut-down liver graft)",
      notes              = "Tested in the stepwise covariate analysis but not retained in the final model (Kassir 2014 Results 'Analysis of covariates and sources of variability' and Discussion: contrast with Staatz et al. 2001 paediatric liver transplant cohort, where cut-down liver from an adult donor was associated with ~7x lower CL/F vs whole liver from a child donor). Cohort 33.3% whole liver, 66.7% cut-down liver (Table 2)."
    ),
    DONOR_AGE = list(
      description        = "Liver donor age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested in the stepwise covariate analysis but not retained in the final model (Kassir 2014 Results and Discussion). Donor age range 0.58-66 years (Kassir 2014 Results)."
    ),
    POD = list(
      description        = "Time post-transplantation (months)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested in the stepwise covariate analysis and not retained (Kassir 2014 Results). Time post-transplantation has been reported as a CL/F covariate by several other paediatric liver transplant tacrolimus authors (Sam 2003, Garcia Sanchez 2001, Wallin 2011, Fukudo 2006); Kassir et al. note that no patients were studied during the first 2 weeks post-transplant and the median time post-transplant was 2.5 months, which may explain non-detection. Cohort range 0.5-188.2 months (Table 3)."
    ),
    ALT = list(
      description        = "Alanine aminotransferase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested in the stepwise covariate analysis (also evaluated as a discrete < 45 vs >= 45 indicator per Kassir 2014 Table 2 and reference [12]) and not retained in the final model (Kassir 2014 Results). Cohort median 44.5 U/L, range 14-140 (Table 3)."
    ),
    AST = list(
      description        = "Aspartate aminotransferase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested in the stepwise covariate analysis and not retained in the final model (Kassir 2014 Results). Cohort median 33 U/L (Table 3)."
    ),
    GGT = list(
      description        = "Gamma-glutamyl transferase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested in the stepwise covariate analysis and not retained in the final model (Kassir 2014 Results). Cohort median 54.5 U/L (Table 3)."
    ),
    TBIL = list(
      description        = "Total bilirubin",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested in the stepwise covariate analysis and not retained in the final model (Kassir 2014 Results). Also used as a binary indicator of hepatic impairment (TBIL > 68.4 umol/L and/or ALT > 2x age-group upper limit; 17% of cohort, Table 2). Cohort median 11.0 umol/L (Table 3)."
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested in the stepwise covariate analysis (also as a discrete <= 32 vs > 32 g/L indicator) and not retained in the final model (Kassir 2014 Results). Cohort median 36 g/L (Table 3)."
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested in the stepwise covariate analysis and not retained in the final model (Kassir 2014 Results). Cohort median 0.4 mg/dL (Table 3)."
    ),
    CRCL = list(
      description        = "Creatinine clearance",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Calculated by the Schwartz equation for patients <= 18 years and by Cockcroft-Gault for patients > 18 years (Kassir 2014 Table 3 footnote). Tested in the stepwise covariate analysis and not retained in the final model. Cohort median 155.6 mL/min (Table 3)."
    ),
    HCT = list(
      description        = "Haematocrit",
      units              = "percent",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested in the stepwise covariate analysis (also as a discrete < 33% vs >= 33% indicator) and not retained in the final model (Kassir 2014 Results), in contrast to several other published popPK models (Andrews 2017, Zhao 2009). Cohort median 34.5%, range 25-44.1 (Table 3)."
    ),
    STEROID = list(
      description        = "Concomitant steroid use (prednisone)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no steroid)",
      notes              = "Tested in the stepwise covariate analysis and not retained in the final model (Kassir 2014 Results and Discussion). Cohort 66.7% on prednisone (Table 2)."
    ),
    CYP3A4_INH = list(
      description        = "Concomitant clinically relevant CYP3A4 inhibitor use",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no CYP3A4 inhibitor)",
      notes              = "Inhibitors included amlodipine, lansoprazole, and diltiazem (Kassir 2014 Table 2 footnote). Tested in the stepwise covariate analysis and not retained in the final model. Cohort 26.7% on a CYP3A4 inhibitor (Table 2)."
    ),
    FORM = list(
      description        = "Drug formulation indicator (1 = oral suspension; 0 = capsule)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (capsule)",
      notes              = "Tested in the stepwise covariate analysis and not retained in the final model (Kassir 2014 Results). Cohort 50% capsule, 50% extemporaneously compounded 0.5 mg/mL oral suspension (Table 2). The compounded suspension was prepared by CHU Sainte-Justine hospital or community pharmacies (possibly with different reconditioning methods), so its bioavailability is unknown -- this could not be evaluated as a covariate effect in the present cohort (Discussion)."
    )
  )

  ini({
    # Structural PK -- Kassir 2014 Table 4 final-model point estimates plus the
    # paper-stated final-model equations on p. 1056. Reference subject: a
    # paediatric liver transplant recipient at the cohort median weight of
    # 20 kg, beyond the first 2 weeks post-transplant, on standard tacrolimus
    # therapy.
    lka   <- log(0.342) ; label("First-order absorption rate constant ka at WT = 20 kg (1/h)")          # Kassir 2014 Table 4 final ka = 0.342 1/h (RSE 33.3%)
    ltlag <- log(0.433) ; label("Absorption lag time tlag (h)")                                          # Kassir 2014 Table 4 final tlag = 0.433 h (RSE 4.2%)
    lcl   <- log(12.1)  ; label("Apparent oral clearance CL/F at WT = 20 kg (L/h)")                      # Kassir 2014 Table 4 final CL/F = 12.1 L/h (RSE 10.1%)
    lvc   <- log(31.3)  ; label("Apparent central volume V1/F at WT = 20 kg (L)")                        # Kassir 2014 Table 4 final V1/F = 31.3 L (RSE 42.8%)
    lq    <- log(30.7)  ; label("Apparent inter-compartmental clearance Q2/F at WT = 20 kg (L/h)")       # Kassir 2014 Table 4 final Q2/F = 30.7 L/h (RSE 29.3%)
    lvp   <- fixed(log(290)) ; label("Apparent peripheral volume V2/F at WT = 20 kg (L; fixed)")         # Kassir 2014 Table 4 final V2/F = 290 L (fixed; Table 4 footnote: "V2/F was fixed to a value estimated from a previous run in order to stabilize the model")

    # Allometric exponents -- Kassir 2014 Methods 'Population pharmacokinetic
    # analysis' equation block (citing references [20, 21], i.e., Anderson and
    # Holford theory-based allometric scaling). All exponents are fixed at the
    # theory values; the paper's results state "Bodyweight was included in all
    # pharmacokinetic parameters as an allometric fixed term" (Table 4 footnote)
    # to confirm the fixed status.
    e_wt_cl <- fixed(0.75)  ; label("Allometric exponent of (WT/20 kg) on CL/F (unitless; fixed)")        # Kassir 2014 Methods equation block (Anderson-Holford theory)
    e_wt_q  <- fixed(0.75)  ; label("Allometric exponent of (WT/20 kg) on Q2/F (unitless; fixed)")        # Kassir 2014 Methods equation block
    e_wt_vc <- fixed(1)     ; label("Allometric exponent of (WT/20 kg) on V1/F (unitless; fixed)")        # Kassir 2014 Methods equation block
    e_wt_vp <- fixed(1)     ; label("Allometric exponent of (WT/20 kg) on V2/F (unitless; fixed)")        # Kassir 2014 Methods equation block
    e_wt_ka <- fixed(-0.25) ; label("Allometric exponent of (WT/20 kg) on ka (unitless; fixed)")          # Kassir 2014 Methods equation block (ka = theta * (WT/WTmedian)^(-0.25))

    # Diagonal inter-individual variability on CL/F, V1/F, and Q2/F. Kassir 2014
    # Table 4 reports BSV as %CV from an exponential (log-normal) random-effect
    # model; no inter-eta covariance was reported in the final model (the
    # paper's Methods state "Covariance between parameters was also examined"
    # but Table 4 lists only diagonal BSV entries). Convert %CV to the
    # internal log-scale variance via omega^2 = log(1 + CV^2):
    #   CL/F  CV 55.6%  -> log(1 + 0.556^2) = 0.269367
    #   V1/F  CV 126.1% -> log(1 + 1.261^2) = 0.951705
    #   Q2/F  CV 84.0%  -> log(1 + 0.840^2) = 0.533917
    etalcl ~ 0.269367   # Kassir 2014 Table 4 BSV CL/F = 55.6% CV (RSE 9.6%)
    etalvc ~ 0.951705   # Kassir 2014 Table 4 BSV V1/F = 126.1% CV (RSE 18%)
    etalq  ~ 0.533917   # Kassir 2014 Table 4 BSV Q2/F = 84.0% CV (RSE 21.3%)

    # Residual unexplained variability -- Kassir 2014 Methods 'Population
    # pharmacokinetic analysis' final model used a proportional error structure
    # for tacrolimus whole-blood concentrations. Table 4 reports the
    # proportional residual error as 20.3% (RSE 12.1%).
    propSd <- 0.203 ; label("Proportional residual error (fraction)")                                     # Kassir 2014 Table 4 residual proportional error = 20.3% (RSE 12.1%)
  })

  model({
    # Body-weight scaling factor relative to the 20 kg cohort median
    # (Kassir 2014 final-model equations on p. 1056). Applied to all five
    # PK parameters at the paper's fixed allometric exponents.
    wt20 <- WT / 20

    # Individual PK parameters with the Kassir 2014 covariate equations.
    # ka, tlag, and V2/F carry no IIV (Kassir 2014 Table 4 reports BSV on
    # CL/F, V1/F, and Q2/F only).
    ka   <- exp(lka)         * wt20 ^ e_wt_ka
    cl   <- exp(lcl + etalcl) * wt20 ^ e_wt_cl
    vc   <- exp(lvc + etalvc) * wt20 ^ e_wt_vc
    q    <- exp(lq  + etalq)  * wt20 ^ e_wt_q
    vp   <- exp(lvp)         * wt20 ^ e_wt_vp
    tlag <- exp(ltlag)

    # Two-compartment oral disposition. Dose lands in `depot`; bioavailability
    # is implicit in the apparent CL/F and V/F parameterisation. Lag time is
    # applied via alag(depot).
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # Tacrolimus whole-blood concentrations reported in ng/mL. Dose in mg, vc
    # in L, so central/vc is in mg/L = ug/mL; multiply by 1000 to convert to
    # ng/mL.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
