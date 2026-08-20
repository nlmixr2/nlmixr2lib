Wojciechowski_2023_ritlecitinib_final <- function() {
  description <- "Final model (iteration 3 of 3) for oral ritlecitinib: two-compartment disposition with first-order absorption and a direct-response non-stationary (autoinhibitory) Imax effect of the peripheral-compartment concentration on both apparent clearance and bioavailability. Parameters were re-estimated with the NONMEM PRIOR (NWPRI) frequentist-prior subroutine, using the updated model as the prior, on an evaluation dataset dominated by sparsely sampled phase IIb/III alopecia areata patients plus healthy participants and severe renal impairment participants (601 individuals, 2944 concentrations). Carries an alopecia areata effect on CL/F, a severe renal impairment effect on F, and inflammatory-disease scaling of the IIV and proportional residual error magnitudes. Allometric weight scaling with fixed exponents 0.75 on clearances and 1 on volumes referenced to 70 kg. This is the iteration that underwrote the approved Litfulo (ritlecitinib) product label."
  reference <- paste(
    "Wojciechowski J, Purohit VS, Huh Y, Banfield C, Nicholas T.",
    "Evolution of Ritlecitinib Population Pharmacokinetic Models During",
    "Clinical Drug Development. Clin Pharmacokinet. 2023;62(12):1765-1779.",
    "doi:10.1007/s40262-023-01318-3. PMCID PMC10684409.",
    "Final-model column of Table 2; complete NONMEM control stream,",
    "including the $PRIOR NWPRI blocks, in the Electronic Supplementary",
    "Material.",
    sep = " "
  )
  vignette <- "Wojciechowski_2023_ritlecitinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot       = list(analyte = "ritlecitinib", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "ritlecitinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ritlecitinib", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline total body weight; allometrically scales CL/F, Q/F, Vc/F and Vp/F",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline weight (NONMEM column BWT). Reference 70 kg with fixed exponents 0.75 on CL/F and Q/F and 1.00 on Vc/F and Vp/F; Wojciechowski 2023 Table 2 footnote. Final-model analysis population median 68.5 kg (range 29.6-131), which includes adolescents from 12 years of age. Wojciechowski 2023 Sect. 4 concludes that allometric weight scaling accounted for adolescent exposure differences and no weight-based dose adjustment was warranted.",
      source_name        = "BWT"
    ),
    DIS_ALOPECIA_AREATA = list(
      description        = "Alopecia areata patient indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy participant or severe renal impairment participant)",
      notes              = "Derived from the NONMEM PTST patient-type column (PTST = 3). Carries a multiplicative effect on CL/F and is the only member of the inflammatory-disease group (PTST in {1, 2, 3, 5}) present in the final-model evaluation dataset, so it alone scales the IIV and residual-error magnitudes here. The control stream names healthy participants and severe renal impairment participants as the CL/F reference.",
      source_name        = "PTST"
    ),
    RENALIMP_SEV = list(
      description        = "Severe renal impairment indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal renal function; healthy participants and alopecia areata patients)",
      notes              = "Derived from the NONMEM PTST patient-type column (PTST = 7), from the dedicated phase I renal impairment study NCT04037865 (n = 8 in the final-model dataset). Carries a multiplicative effect on F only. Severe renal impairment participants are NOT part of the inflammatory-disease group that scales the IIV and residual-error magnitudes (the control stream tests PTST > 0 AND PTST < 6). Wojciechowski 2023 Sect. 4 attributes the higher F to reduced first-pass metabolism in chronic kidney disease and reports an average 47% increase in steady-state AUC(tau) and 41% increase in Cmax.",
      source_name        = "PTST"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 601L,
    n_studies      = 3L,
    n_observations = 2944L,
    age_range      = "12.0-72.0 years",
    age_median     = "32.0 years",
    weight_range   = "29.6-131 kg",
    weight_median  = "68.5 kg",
    sex_female_pct = 60.6,
    race_ethnicity = c(White = 66.4, Asian = 27.3, Black = 3.8, Other = 1.7, Missing = 0.8),
    disease_state  = "Alopecia areata patients (584; 97.2%), healthy participants (9; 1.5%) and severe renal impairment participants (8; 1.3%).",
    dose_range     = "10-200 mg/day oral ritlecitinib in the phase IIb/III alopecia areata studies, plus phase I healthy participant and renal impairment doses; overall programme range 5-800 mg/day",
    regions        = "Global (phase IIb/III alopecia areata programme plus phase I studies; see Table S1 of the Electronic Supplementary Material)",
    renal_function = "Normal, plus 8 participants with severe renal impairment (NCT04037865)",
    hepatic_function = "Normal",
    notes          = "Demographics from Wojciechowski 2023 Table 1, final-model column (the Table 1 total row reads 599 while the counts sum to 601 and the text states 601 individuals; the text value is used here). Median baseline creatinine clearance 117 mL/min (range 14.7-288). This is the first iteration to include adolescents (age range starts at 12 years), supporting the Litfulo label in adults and adolescents aged 12 years and older. Below-limit-of-quantification observations were excluded during estimation."
  )

  # Implementation notes (see the vignette 'Assumptions and deviations'
  # section for the full justification of each item):
  # * Structure. Identical to the base and updated models; Sect. 2.5.1
  #   states the structure of the updated model was not re-evaluated, and
  #   Sect. 4 confirms the structural elements were retained across all
  #   three iterations. The deposited final-model $DES block is
  #     AUTOINH = 1 + IMAXP*PERIC/(IC50P + PERIC)
  #     INHF    = 1 + (1 - AUTOINH)
  #     FA      = COVPTSTF
  #     DADT(1) = -KA*A(1)
  #     DADT(2) = KA*A(1)*INHF*FA - Q*CENTC + Q*PERIC - CL*AUTOINH*CENTC
  #     DADT(3) = Q*CENTC - Q*PERIC
  #   Note that the final-model control stream carries NO food, dose,
  #   or formulation effects on ka and no depot-loss term: the
  #   evaluation dataset had no high-fat-meal, 800 mg or capsule records
  #   to inform them, so those eight covariate effects from the updated
  #   model are absent here rather than fixed. Simulate updated-model
  #   scenarios (RA, UC, vitiligo, moderate hepatic impairment, high-fat
  #   meal, capsule) with Wojciechowski_2023_ritlecitinib_updated, as
  #   Wojciechowski 2023 Fig. 4 does.
  # * Frequentist prior. Sect. 2.5.2 and the ESM $PRIOR NWPRI blocks:
  #   the updated model's final estimates were used as normally
  #   distributed priors on the fixed effects (with the updated model's
  #   variance-covariance matrix as the informative prior variance) and
  #   as inverse-Wishart priors on the IIV. The residual-error SD and the
  #   new severe renal impairment covariate carried no prior information
  #   ($THETAPV 10000 for PTVPTSTF7). The values encoded below are the
  #   resulting POSTERIOR estimates reported in the Table 2 final-model
  #   column, not the priors.
  # * Inlined state expressions, units, IIV scale, residual error and the
  #   inflammatory-disease SD scaling all follow the base model; see the
  #   implementation notes in Wojciechowski_2023_ritlecitinib_base.R.
  # * Residual-error interval. Table 2 prints the final-model
  #   sigma^2_pro as 35.6 (95% CI 24.9-36.4). The lower bound does not
  #   bracket symmetrically and is almost certainly a typographical error
  #   for 34.9 (which gives midpoint 35.65). Only the point estimate
  #   enters the model. See the vignette Errata.
  ini({
    # ----- Structural parameters (Wojciechowski 2023 Table 2, final-model column) -----
    lcl   <- log(107);   label("Apparent clearance CL/F at 70 kg, uninhibited (L/h)")                          # Table 2 final: 107 (95% CI 98.6-116)
    lvc   <- log(151);   label("Apparent central volume Vc/F at 70 kg (L)")                                    # Table 2 final: 151 (95% CI 147-156)
    lq    <- log(0.297); label("Apparent inter-compartmental clearance Q/F at 70 kg (L/h)")                    # Table 2 final: 0.297 (95% CI 0.262-0.332)
    lvp   <- log(4.87);  label("Apparent peripheral volume Vp/F at 70 kg (L)")                                 # Table 2 final: 4.87 (95% CI 4.55-5.20)
    lka   <- log(7.91);  label("First-order absorption rate constant ka (1/h)")                                # Table 2 final: 7.91 (95% CI 6.58-9.25)

    # ----- Non-stationary (autoinhibition) parameters -----
    imax  <- -0.452;     label("Maximum non-stationary fractional effect Imax,P on CL/F and F (unitless)")     # Table 2 final: -0.452 (95% CI -0.485 to -0.419)
    lic50 <- log(16.5);  label("Peripheral concentration giving half the maximum non-stationary effect IC50,P (ng/mL)")  # Table 2 final: 16.5 (95% CI 13.3-19.7)

    # ----- Allometric exponents (fixed, Table 2 footnote) -----
    e_wt_cl_q  <- fixed(0.75); label("Allometric WT exponent shared by CL/F and Q/F (unitless)")               # Table 2 footnote: exponent 0.75 on clearance parameters, reference 70 kg
    e_wt_vc_vp <- fixed(1.00); label("Allometric WT exponent shared by Vc/F and Vp/F (unitless)")              # Table 2 footnote: exponent 1 on volume parameters, reference 70 kg

    # ----- Categorical covariate effects -----
    e_alopecia_areata_cl   <- -0.260; label("Fractional effect of alopecia areata on CL/F (unitless)")                       # Table 2 final: -0.260 (95% CI -0.313 to -0.207); ESM final-model $PK THETA(12) PTSTCL3
    e_renalimp_sev_fdepot  <-  0.353; label("Fractional effect of severe renal impairment on bioavailability F (unitless)")  # Table 2 final: 0.353 (95% CI 0.218-0.488); ESM final-model $PK THETA(13) PTSTF7

    # ----- Inflammatory-disease scaling of the variability magnitudes -----
    # These multiply the SD (eta and residual SD), not the variance. In
    # this dataset the inflammatory-disease group is the alopecia areata
    # cohort alone.
    e_inflam_iiv_cl <- 1.61;  label("Fractional effect of inflammatory disease on the SD of IIV on CL/F (unitless)")  # Table 2 final: 1.61 (95% CI 1.29-1.93); ESM final-model $PK THETA(10) VARCLPTST
    e_inflam_iiv_vc <- 1.43;  label("Fractional effect of inflammatory disease on the SD of IIV on Vc/F (unitless)")  # Table 2 final row 2 of the duplicated omega^2 label: 1.43 (95% CI 0.931-1.94); ESM final-model $PK THETA(11) VARV2PTST
    e_inflam_propsd <- 0.290; label("Fractional effect of inflammatory disease on the proportional residual SD (unitless)")  # Table 2 final: 0.290 (95% CI 0.255-0.325); ESM final-model $ERROR THETA(9) RUVPROPTST

    # ----- Between-subject variability (diagonal) -----
    etalcl ~ 0.035344  # Table 2 final: 18.8 percent-CV column (= 100*omega) per the Table 2 footnote, so omega^2 = 0.188^2
    etalvc ~ 0.013225  # Table 2 final: 11.5 percent-CV column (= 100*omega) per the Table 2 footnote, so omega^2 = 0.115^2

    # ----- Residual error -----
    propSd <- 0.356; label("Proportional residual SD for the reference (healthy participant) group (fraction)")  # Table 2 final: 35.6 "% CV"; ESM final-model $ERROR THETA(8) RUVPRO with $SIGMA 1 FIX
  })
  model({
    ref_wt <- 70  # Table 2 footnote: allometric reference weight

    # ----- Derived covariate terms -----
    # Inflammatory-disease group = PTST in {1, 2, 3, 5}. Only alopecia
    # areata is present in the final-model evaluation dataset; healthy
    # participants (PTST = 0) and severe renal impairment (PTST = 7) are
    # the reference.
    inflam <- DIS_ALOPECIA_AREATA

    covptstcl <- 1 + e_alopecia_areata_cl * DIS_ALOPECIA_AREATA
    covptstf  <- 1 + e_renalimp_sev_fdepot * RENALIMP_SEV
    covvarcl  <- 1 + e_inflam_iiv_cl * inflam
    covvarvc  <- 1 + e_inflam_iiv_vc * inflam

    # ----- Individual parameters -----
    cl <- exp(lcl + etalcl * covvarcl) * (WT / ref_wt)^e_wt_cl_q * covptstcl
    vc <- exp(lvc + etalvc * covvarvc) * (WT / ref_wt)^e_wt_vc_vp
    q  <- exp(lq)  * (WT / ref_wt)^e_wt_cl_q
    vp <- exp(lvp) * (WT / ref_wt)^e_wt_vc_vp

    ka <- exp(lka)

    # IC50 is reported in ng/mL; the peripheral concentration below is in
    # mg/L, so convert exactly as the control stream does (THETA(6)/1000).
    ic50 <- exp(lic50) / 1000

    # ----- ODE system -----
    # Concentrations are written inline; see implementation notes.
    d/dt(depot) <- -ka * depot
    d/dt(central) <-
      ka * depot *
        (2 - (1 + imax * (peripheral1 / vp) / (ic50 + peripheral1 / vp))) *
        covptstf -
      q * (central / vc) + q * (peripheral1 / vp) -
      cl * (1 + imax * (peripheral1 / vp) / (ic50 + peripheral1 / vp)) *
        (central / vc)
    d/dt(peripheral1) <- q * (central / vc) - q * (peripheral1 / vp)

    # ----- Observation and residual error -----
    Cc <- central / vc * 1000
    propSdEff <- propSd * (1 + e_inflam_propsd * inflam)
    Cc ~ prop(propSdEff)
  })
}
