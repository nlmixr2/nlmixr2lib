Jang_2021_oxcarbazepine <- function() {
  description <- "One-compartment population PK model with first-order absorption and elimination for the monohydroxy derivative (MHD, licarbazepine) of oral oxcarbazepine in Korean adults with epilepsy (Jang 2021). Oxcarbazepine doses drive MHD directly (parent not modelled); apparent clearance CL/F and apparent volume V/F scale with body weight by estimated power exponents centred at 66 kg. Inter-individual variability on CL/F and on the absorption rate constant ka (the control stream fixes the V/F variance to zero), with a proportional residual error."
  reference <- paste(
    "Jang Y, Yoon S, Kim TJ, Lee S, Yu KS, Jang IJ, Chu K, Lee SK.",
    "Population pharmacokinetic model development and its relationship with",
    "adverse events of oxcarbazepine in adult patients with epilepsy.",
    "Sci Rep. 2021;11:6370. doi:10.1038/s41598-021-85920-0",
    sep = " "
  )
  vignette <- "Jang_2021_oxcarbazepine"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects on CL/F (exponent 0.67) and V/F (exponent 0.96) centred at 66 kg (Jang 2021 Table 2; supplementary NONMEM code CL = TVCL*EXP(ETA(1))*(BW/66)**THETA(4), V = TVV*EXP(ETA(2))*(BW/66)**THETA(5)). Cohort mean 65.8 kg (range 39-116 kg; Jang 2021 Results, Patient characteristics).",
      source_name = "BW"
    )
  )

  covariatesDataExcluded <- list(
    CONMED_EIAED = list(
      description = "Concomitant enzyme-inducing antiseizure medication indicator (carbamazepine, phenytoin, phenobarbital or valproic acid as grouped by the authors)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no enzyme-inducing antiseizure medication)",
      notes = "Tested on CL/F (an estimated 7% increase) but not retained because it improved neither the OFV nor the goodness-of-fit plots (Jang 2021 Results, Population PK analysis). No point estimate beyond the 7% is reported.",
      source_name = "EIASMs"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in stepwise covariate selection; not retained (Jang 2021 Results, Population PK analysis).",
      source_name = "age"
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened in stepwise covariate selection; not retained (Jang 2021 Results, Population PK analysis).",
      source_name = "sex"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "oxcarbazepine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "monohydroxy derivative (MHD)", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 487L,
    n_studies = 2L,
    n_observations = 748L,
    age_range = "16-80 years (registry cohort)",
    age_mean = "39.2 +/- 14.7 years (registry cohort)",
    weight_range = "39-116 kg (registry cohort)",
    weight_mean = "65.8 +/- 12.5 kg (registry cohort)",
    sex_female_pct = 40.9,
    race_ethnicity = c(Korean = 100),
    disease_state = "Adults with epilepsy on chronic oxcarbazepine therapy (Study 1, Epilepsy Registry Cohort; 65.3% on at least one concomitant antiseizure medication, most often levetiracetam or topiramate; 5.1% on an enzyme-inducing antiseizure medication; 63.8% seizure-free), pooled with a single-dose oral-loading study in patients with epilepsy (Study 2).",
    dose_range = "Study 1: stable oral oxcarbazepine for at least one month, mean 999 mg/day (range 150-2100 mg/day), 438 of 447 patients twice daily, 5 once daily, 4 three times daily. Study 2: single 30 mg/kg oral loading dose with serial sampling at 2, 4, 6, 8, 12, 14, 16 and 24 h.",
    regions = "South Korea (Seoul National University Hospital)",
    notes = "The PK dataset held 748 MHD concentrations from 487 patients (Jang 2021 Methods): the 447 registry patients with complete dosing and sampling records plus the 40 single-dose-study patients (Kim DW et al. Epilepsia 2012;53:e9-12). Demographics (Jang 2021 Table 1 and Results) describe the 447 registry patients only. Plasma MHD was measured by validated LC-MS/MS."
  )

  ini({
    # Final estimates: Jang 2021 Table 2 ('Estimates (RSE)' column). The
    # supplementary control stream (ADVAN2 TRANS2) fixes the structure; its
    # $THETA/$OMEGA values are initial estimates and are not used.
    lcl <- log(1.65); label("Apparent clearance CL/F of MHD at 66 kg (L/h)") # Jang 2021 Table 2 theta1 = 1.65 (1.8% RSE)
    lvc <- log(59.0); label("Apparent volume of distribution V/F of MHD at 66 kg (L)") # Jang 2021 Table 2 theta3 = 59.0 (4.5% RSE)
    lka <- log(0.34); label("First-order absorption rate constant ka (1/h)") # Jang 2021 Table 2 ka = 0.34 (9.7% RSE)

    e_wt_cl <- 0.67; label("Power exponent of body weight on CL/F (unitless)") # Jang 2021 Table 2 theta2 = 0.67 (14.2% RSE)
    e_wt_vc <- 0.96; label("Power exponent of body weight on V/F (unitless)") # Jang 2021 Table 2 theta4 = 0.96 (18.1% RSE)

    # IIV reported as %CV of an exponential model; omega^2 = log(1 + CV^2).
    # The control stream carries ETA(1) on CL, ETA(2) on V with $OMEGA '0 FIX'
    # and ETA(3) on KA, so the second estimated IIV in Table 2 (labelled
    # 'IIV V2/F') is the ka variability; see the vignette for the check against
    # the Figure 1 VPC.
    etalcl ~ 0.08183 # Jang 2021 Table 2 'IIV CL/F' 29.2 %CV: log(1 + 0.292^2)
    etalka ~ 0.15952 # Jang 2021 Table 2 second IIV row 41.5 %CV (ETA(3) on KA in the supplementary code): log(1 + 0.415^2)

    # Supplementary code: W = SQRT(THETA(6)**2 + THETA(7)**2 * IPRED**2) with
    # THETA(6) = 0.0001 ng/mL FIX and $SIGMA 1 FIX, so THETA(7) is the
    # proportional SD; the negligible fixed additive term is omitted.
    propSd <- 0.13; label("Proportional residual error (fraction)") # Jang 2021 Table 2 'Proportional error (SD)' = 0.13 (10.8% RSE)
  })

  model({
    ref_wt <- 66

    cl <- exp(lcl + etalcl) * (WT / ref_wt)^e_wt_cl
    vc <- exp(lvc) * (WT / ref_wt)^e_wt_vc
    ka <- exp(lka + etalka)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Source DV in ng/mL (SC = V/1000 with dose in mg); here mg/L, numerically
    # 1/1000 of the source scale.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
