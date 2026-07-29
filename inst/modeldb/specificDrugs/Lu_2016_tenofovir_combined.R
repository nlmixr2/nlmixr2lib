Lu_2016_tenofovir_combined <- function() {
  description <- "Two-compartment population PK model with first-order absorption and an absorption lag time for tenofovir (300 mg oral TDF once daily) in HIV-1-uninfected African adults receiving once-daily preexposure prophylaxis (Lu 2016, Partners PrEP Study). Combined variant: parameters estimated using a combined data set in which patient-reported dosing records were replaced with MEMS electronic adherence monitoring records where available. Absorption rate constant Ka is fixed at 1.5 /h; absorption lag time ALAG1 = 0.41 h. Apparent oral clearance (CL/F) carries a power-form covariate effect on creatinine clearance (raw Cockcroft-Gault, mL/min) centred at the cohort mean 106 mL/min. Diagonal IIV on CL/F, V1/F, and Ka; combined additive + proportional residual error."
  reference <- "Lu Y, Goti V, Chaturvedula A, Haberer JE, Fossler MJ, Sale ME, Bangsberg D, Baeten JM, Celum CL, Hendrix CW. Population pharmacokinetics of tenofovir in HIV-1-uninfected members of serodiscordant couples and effect of dose reporting methods. Antimicrob Agents Chemother. 2016;60(9):5379-5386. doi:10.1128/AAC.00559-16"
  vignette <- "Lu_2016_tenofovir"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CLCR. Cockcroft-Gault creatinine clearance in raw mL/min (NOT BSA-normalized). Stored under the canonical CRCL column per inst/references/covariate-columns.md (CRCL accepts raw mL/min when the source paper does not apply BSA normalization, with the per-model description recording the assay form). Reference value 106 mL/min is the cohort mean (Lu 2016 Table 1; SD 31 mL/min). The covariate enters CL/F as a power form (CRCL/106)^e_crcl_cl, where e_crcl_cl = 0.376. The paper's covariate-form selection (Methods: 'Models with linear, power, and exponential functions were tested') is not explicitly stated for the retained final-model form; the power interpretation is adopted here because (i) the reported coefficient 0.376 is dimensionless and in the typical range for power exponents on renal-function covariates, (ii) the linear/exponential alternatives with the same numeric coefficient give physically implausible CL values at the extremes of the observed CrCl range, and (iii) NONMEM's standard divisive centring (CRCL/median) for a power model matches the paper's prose 'centered to the median values'.",
      source_name        = "CLCR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 404L,
    n_studies      = 1L,
    age_range      = "mean 35 years (SD 8)",
    age_median     = "35 years",
    weight_range   = "mean 61 kg (SD 11)",
    weight_median  = "61 kg",
    sex_female_pct = 45,
    race_ethnicity = "Sub-Saharan African (Kenyan and Ugandan cohort); race not stratified further in the paper",
    disease_state  = "HIV-1-seronegative; HIV-1-uninfected members of HIV-1-serodiscordant heterosexual couples enrolled in the Partners PrEP Study",
    dose_range     = "300 mg oral tenofovir disoproxil fumarate (TDF) once daily as preexposure prophylaxis (alone or in fixed-dose combination with emtricitabine 200 mg)",
    regions        = "Kenya and Uganda",
    crcl_range     = "mean 106 mL/min (SD 31)",
    samples        = "1278 tenofovir plasma concentrations from 404 participants in the combined data set; MEMS records substituted for patient-reported dosing where available (26% of samples paired with MEMS records). The ancillary adherence substudy contributed 211 participants with MEMS-monitored dosing.",
    adherence      = "97.2-99.1% by unannounced pill counts and MEMS in the ancillary substudy; 97% by pill counts in the main trial",
    notes          = "Cohort demographics per Lu 2016 Table 1 (Combined data set column). Phase 3 randomised, double-blind, placebo-controlled HIV-1 PrEP trial; the PK substudy enrolled HIV-1-seronegative partners on active TDF or TDF+FTC arms. In the Combined data set, patient-reported dosing times were replaced with MEMS electronic-adherence records for the 26% of samples associated with MEMS-recorded openings; PRDI with steady-state assumption was retained for the remaining 74%. Lu 2016 Discussion notes that Ka was fixed (via a local search) due to numerical instabilities and that estimation of additional IIV on V1/F and Ka was supported only in the Combined model."
  )

  ini({
    # Lu 2016 Table 2 'Combined data Final model' column (point estimates).
    # Reference subject: CRCL = 106 mL/min (cohort mean, Table 1).
    # All apparent oral clearances (CL/F, Q/F) in L/h; apparent volumes (Vc/F,
    # Vp/F) in L; ka in 1/h; tlag in h.
    # Ka is FIXED (Lu 2016 Discussion: 'The absorption rate constant was fixed
    # by a local search due to numerical instabilities posed during modeling of
    # the combined data set').
    lka  <- fixed(log(1.5)); label("First-order absorption rate constant ka (1/h, fixed)")          # Lu 2016 Table 2 Combined final Ka = 1.5 (fixed)
    ltlag <- log(0.41);      label("Absorption lag time ALAG1 (h)")                                 # Lu 2016 Table 2 Combined final ALAG1 = 0.41 h
    lcl  <- log(61.5);       label("Apparent oral clearance CL/F at CRCL = 106 mL/min (L/h)")       # Lu 2016 Table 2 Combined final CL = 61.5 L/h
    lvc  <- log(345);        label("Apparent central volume of distribution V1/F (L)")              # Lu 2016 Table 2 Combined final V1 = 345 L
    lq   <- log(231);        label("Apparent inter-compartmental clearance Q/F (L/h)")              # Lu 2016 Table 2 Combined final Q = 231 L/h
    lvp  <- log(830);        label("Apparent peripheral volume of distribution Vp/F (L)")           # Lu 2016 Table 2 Combined final Vp = 830 L

    # Power-form CrCl effect on CL/F. Lu 2016 reports theta_CL-CR = 0.376 as
    # the bootstrap median for the final-model covariate coefficient (Table 2,
    # Combined final model column; 95% CI 0.27-0.49). Applied as
    # (CRCL/106)^e_crcl_cl.
    e_crcl_cl <- 0.376; label("Power exponent of (CRCL/106) on CL/F (unitless)")                     # Lu 2016 Table 2 Combined final theta_CL-CR = 0.376 (median 0.38, 95% CI 0.27-0.49)

    # Diagonal IIV. Lu 2016 Table 2 Combined final-model row:
    #   IIV on CL  = 16% CV -> log(1 + 0.16^2) = log(1.0256) = 0.02528
    #   IIV on V1  = 25% CV -> log(1 + 0.25^2) = log(1.0625) = 0.06062
    #   IIV on Ka  = 61% CV -> log(1 + 0.61^2) = log(1.3721) = 0.31636
    # No inter-eta correlations reported in the source.
    etalcl ~ 0.02528  # Lu 2016 Table 2 Combined final IIV on CL = 16% CV
    etalvc ~ 0.06062  # Lu 2016 Table 2 Combined final IIV on V1 = 25% CV
    etalka ~ 0.31636  # Lu 2016 Table 2 Combined final IIV on Ka = 61% CV

    # Combined additive + proportional residual error (Lu 2016 Table 2
    # Combined final-model row). Concentrations are in ng/mL; additive SD is
    # in ng/mL. Lu 2016 also reports 'IIV on additive error' = 136% CV for the
    # Combined final model; this represents an inter-individual scaling of the
    # additive residual-error magnitude per the Karlsson reference cited in
    # Methods (reference 21). It is NOT encoded structurally here because
    # nlmixr2's canonical residual-error syntax does not naturally express
    # IIV on the error magnitude without breaking parameter-naming
    # conventions; the omission affects only the lower tail of the prediction
    # interval near the LLOQ and does not alter typical-value disposition.
    # See the vignette Assumptions and deviations section.
    addSd  <- 30.2; label("Additive residual error (ng/mL)")                                         # Lu 2016 Table 2 Combined final additive = 30.2 ng/mL
    propSd <- 0.20; label("Proportional residual error (fraction)")                                  # Lu 2016 Table 2 Combined final proportional = 20% CV
  })
  model({
    # Individual PK parameters. Power-form CrCl effect on CL/F centred at the
    # cohort mean 106 mL/min; no covariates on Ka, Vc, Q, or Vp in the
    # Combined final model (Lu 2016 Section 'Population PK model for
    # tenofovir'). Ka has IIV but no Ka covariate; V1 has IIV.
    ka   <- exp(lka + etalka)
    tlag <- exp(ltlag)
    cl   <- exp(lcl + etalcl) * (CRCL / 106) ^ e_crcl_cl
    vc   <- exp(lvc + etalvc)
    q    <- exp(lq)
    vp   <- exp(lvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # Dose in mg, vc in L -> central / vc has units mg/L; * 1000 -> ng/mL.
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
