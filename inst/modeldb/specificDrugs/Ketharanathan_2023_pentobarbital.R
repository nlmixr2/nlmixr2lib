Ketharanathan_2023_pentobarbital <- function() {
  description <- "One-compartment population PK model for intravenous pentobarbital in critically ill children admitted to a paediatric intensive care unit for refractory status epilepticus or severe traumatic brain injury. Body weight is allometrically scaled onto clearance (exponent 0.75) and central volume (exponent 1) against a 70 kg reference, with both exponents fixed. Serum creatinine enters clearance as a power function normalised to 36 umol/L, and C-reactive protein as a second power function normalised to 70 mg/L that is switched on only at or above that threshold. Together the two covariates removed 84% of the inter-individual variability in clearance present in the base model."
  reference <- paste(
    "Ketharanathan N, Lili A, Penning de Vries JM, Wildschut ED, de Hoog M,",
    "Koch BCP, de Winter BCM.",
    "A Population Pharmacokinetic Model of Pentobarbital for Children with",
    "Status Epilepticus and Severe Traumatic Brain Injury.",
    "Clin Pharmacokinet. 2023;62(7):1011-1022.",
    "doi:10.1007/s40262-023-01249-z.",
    sep = " "
  )
  vignette <- "Ketharanathan_2023_pentobarbital"

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric size descriptor, normalised to a 70 kg reference. Exponents fixed at 0.75 on CL and 1 on Vc per Ketharanathan 2023 Sect. 2.4.1 ('Pharmacokinetic parameters were allometrically scaled with fixed exponents (0.75 for CL and 1 for Vd)'). Sect. 3.3 records that the CL exponent was also estimated during covariate analysis and that the confidence interval included the fixed value of 0.75, so the fixed value was retained. Cohort weights span 3-87 kg (Table 1) with an overall median of 10 kg (Abstract); the model is therefore an extrapolation well above its data at the 70 kg reference itself.",
      source_name        = "WGHT"
    ),
    CREAT = list(
      description        = "Serum creatinine concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying power-function covariate on CL, normalised to 36 umol/L (the model-building cohort median). Source NONMEM code: COV=(CREAT/36)**THETA(3) (Supplementary Information Fig. 3). Cohort medians by diagnosis were 23 umol/L (status epilepticus) and 49 umol/L (severe traumatic brain injury); observed range 12-109 umol/L (Table 1). Ketharanathan 2023 Sect. 3.5.1 summarises the effect as 'a doubling of creatinine results roughly in halving of pentobarbital clearance' (2^-0.919 = 0.53). The dosing simulations in Figs. 2 and 4 used 26 umol/L as the scenario input value; this is NOT the normalisation constant, which is 36 umol/L.",
      source_name        = "CREAT"
    ),
    CRP = list(
      description        = "C-reactive protein concentration (standard assay; time-varying)",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Threshold power-function covariate on CL: no effect at or below 70 mg/L, and (CRP/70)^-0.883 above it. Source NONMEM code: IF(CRP.GE.70) COVCL=COV*(CRP/70)**THETA(4) (Supplementary Information Fig. 3). Ketharanathan 2023 Sect. 3.3 reports that several cut-off values were evaluated and that > 70 mg/L gave the best correlation and largest OFV drop. The main-text equations state the condition as CRP > 70 whereas the control stream uses CRP >= 70; the two agree numerically because (70/70)^-0.883 = 1 exactly. Observed range up to 312.8 mg/L (Table 1). The model is only defined for CRP > 0.",
      source_name        = "CRP"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 36L,
    n_studies      = 1L,
    age_range      = "0.05-16 years",
    age_median     = "1.3 years (overall); 0.31 years (status epilepticus), 12 years (severe traumatic brain injury)",
    weight_range   = "3-87 kg",
    weight_median  = "10 kg (overall); 6.3 kg (status epilepticus), 40 kg (severe traumatic brain injury)",
    sex_female_pct = 42, # Table 1: male 14/22 (64%) status epilepticus + 7/14 (50%) severe traumatic brain injury = 21/36; female 8 + 7 = 15/36 = 42%. Table 1 is typeset with its sub-row labels offset one line below their values - the "Gender, n (%)" and "Pentobarbital*" category headers each carry a value row while the last sub-row of each block ("Female", "Max Infusion") is blank. The offset direction is pinned by the sibling block: read with the same offset, the severe traumatic brain injury loading dose is 1 [0-20] mg/kg, matching the per-protocol 1 mg/kg of Sect. 2.1, and the infusion min/max ranges reproduce the 0.5-10 and 0.05-5 mg/kg/h ranges quoted in Sect. 3.1. Read without the offset they do not.

    race_ethnicity = "Not reported",
    disease_state  = "Critically ill children (< 18 years) admitted to a paediatric intensive care unit with either refractory status epilepticus (n = 22) or severe traumatic brain injury defined by a Glasgow Coma Scale <= 8 with refractory intracranial hypertension (n = 14). Diagnosis itself was tested as a covariate and was not retained in the final model.",
    dose_range     = "Continuous intravenous infusion. Per-protocol dosing: status epilepticus 5 mg/kg loading dose then 3 mg/kg/h; traumatic brain injury 1 mg/kg loading dose then 2 mg/kg/h. Observed infusion rates 0.5-10 mg/kg/h (status epilepticus, median 3) and 0.05-5 mg/kg/h (severe traumatic brain injury, median 2); infusion durations 1-15 days.",
    regions        = "The Netherlands (single centre: Erasmus MC-Sophia Children's Hospital, Rotterdam)",
    renal_function = "Serum creatinine 12-109 umol/L; estimated GFR (Schwartz, k = 0.365) 30-199 mL/min/1.73 m^2. Ketharanathan 2023 cautions that 67% of heights were imputed from P50 growth curves, so eGFR and KDIGO stage were unreliable in this cohort and neither was retained as a covariate.",
    notes          = "Retrospective single-centre cohort, January 2007 - September 2021; 178 pentobarbital serum concentrations from 36 children (Table 1). Median measured concentration 27.5 mg/L (range 0.1-106). Seven samples (3%) outside the 0.5-90 mg/L assay range were retained using extrapolated concentrations. An independent cohort of 9 children (6 severe traumatic brain injury, 3 status epilepticus; October 2019 - February 2023, 60 samples) was used for external validation by stratified VPC and is not part of the 36 fitted here. Patients on ECMO or concurrent haemodialysis were excluded."
  )

  ini({
    # Structural parameters - Ketharanathan 2023 Table 2, "Final model" column.
    # Both are reference values at WT = 70 kg, CREAT = 36 umol/L and CRP <= 70 mg/L.
    lcl <- log(3.59); label("Clearance at the 70 kg reference (L/h)")                    # Table 2 final model: CL = 3.59 L/h/70 kg (11% RSE); bootstrap median 3.69, 95% CI 2.80-4.74. Also THETA(1) in Supplementary Fig. 3.
    lvc <- log(142);  label("Central volume of distribution at the 70 kg reference (L)") # Table 2 final model: Vd = 142 L/70 kg (6% RSE); bootstrap median 143, 95% CI 120-158. Also THETA(2) in Supplementary Fig. 3.

    # Allometric exponents - fixed, not estimated. Sect. 2.4.1: "allometrically
    # scaled with fixed exponents (0.75 for CL and 1 for Vd)". Sect. 3.3 confirms
    # the CL exponent was re-estimated as a covariate and the fixed value of 0.75
    # fell inside its confidence interval, so the fixed value was kept. Both appear
    # as hardcoded literals (not THETAs) in the Supplementary Fig. 3 control stream.
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL (unitless)") # Sect. 2.4.1 / Sect. 3.3; control stream TVCL=THETA(1)*(WGHT/70)**0.75
    e_wt_vc <- fixed(1);    label("Allometric exponent of body weight on Vc (unitless)") # Sect. 2.4.1; control stream V=THETA(2)*(WGHT/70)

    # Covariate effects on CL, both power functions on a normalised covariate.
    e_creat_cl <- -0.919; label("Power exponent of serum creatinine on CL (unitless)")                      # Table 2 final model: creatinine effect on CL = -0.919 (21% RSE); bootstrap median -0.909, 95% CI -1.192 to -0.242. THETA(3), applied as (CREAT/36)^e_creat_cl.
    e_crp_cl   <- -0.883; label("Power exponent of C-reactive protein on CL above 70 mg/L (unitless)")      # Table 2 final model: CRP effect on CL = -0.883 (14% RSE), footnote "*If CRP>70 mg/L"; bootstrap median -0.923, 95% CI -1.25 to -0.613. THETA(4), applied as (CRP/70)^e_crp_cl only when CRP >= 70.

    # Inter-individual variability. Exponential IIV on CL only; Sect. 3.2 records
    # that "Including an IIV resulted in a model improvement for CL only", so there
    # is no IIV on Vc. The control stream $OMEGA is 0.121 on the log scale, and
    # sqrt(0.121) = 0.348 reproduces the 34.8% reported in Table 2 - i.e. the paper
    # tabulates the approximate CV, not the variance.
    etalcl ~ 0.121 # Supplementary Fig. 3 $OMEGA 0.121; Table 2 final model IIV CL = 34.8% (24% RSE) [23% shrinkage]; bootstrap median 32.1%, 95% CI 15.8-49.6%.

    # Residual error - additive only. Sect. 3.2: "We described the residual error
    # with a combined error model first. The residual error was described as an
    # additive error, as the proportional error was estimated close to zero."
    #
    # DISCREPANCY (see the vignette Assumptions and deviations section): the
    # Supplementary Fig. 3 control stream prints its $ERROR block as
    #   W = IPRED+THETA(5) ; Y = IPRED+W*ERR(1)  with $SIGMA 1 FIX
    # which would be a residual SD of IPRED + 7.22 mg/L, i.e. a 100% proportional
    # term with a hardcoded coefficient of exactly 1 on top of the additive term.
    # That reading is rejected here as a transcription artefact of the abandoned
    # combined-error model (deleting "*THETA(6)" from "W = IPRED*THETA(6)+THETA(5)"
    # leaves exactly the printed line). It contradicts the quoted Sect. 3.2 text, it
    # contradicts Table 2 - which reports a single residual row in concentration
    # units, "Additional (mg/L) 7.22", with 10% RSE, 6% shrinkage and a bootstrap
    # 95% CI of 5.72-8.54 mg/L - and there is no sixth THETA in the $THETA block for
    # a proportional coefficient to have been estimated into. A residual SD of about
    # 35 mg/L at the cohort median observation of 27.5 mg/L is also irreconcilable
    # with the VPCs in Fig. 1. The additive-only encoding below matches every
    # reported estimate.
    addSd <- 7.22; label("Additive residual error (mg/L)") # Table 2 final model: additive residual = 7.22 mg/L (10% RSE) [6% shrinkage]; bootstrap median 7.11, 95% CI 5.72-8.54. THETA(5) in Supplementary Fig. 3.
  })

  model({
    # 1. Derived covariate terms.

    # Allometric size scaling to a 70 kg reference (Sect. 2.4.1).
    sizeCl <- (WT / 70)^e_wt_cl
    sizeVc <- (WT / 70)^e_wt_vc

    # Serum creatinine power effect on CL, normalised to the 36 umol/L cohort
    # median. Control stream: COV=(CREAT/36)**THETA(3).
    creatEff <- (CREAT / 36)^e_creat_cl

    # C-reactive protein power effect on CL, normalised to the 70 mg/L cut-off and
    # applied only at or above it. Control stream:
    # COVCL=COV ; IF(CRP.GE.70) COVCL=COV*(CRP/70)**THETA(4).
    crpEff <- 1
    if (CRP >= 70) {
      crpEff <- (CRP / 70)^e_crp_cl
    }

    # 2. Individual parameters. Exponential IIV on CL only.
    cl <- exp(lcl + etalcl) * sizeCl * creatEff * crpEff
    vc <- exp(lvc) * sizeVc

    # 3. Micro-constants. NONMEM ADVAN1 TRANS2 (one compartment, CL/V) with S1 = V.
    kel <- cl / vc

    # 4. ODE system. Intravenous administration only - loading doses as bolus and
    # maintenance as a continuous infusion both enter the central compartment.
    d/dt(central) <- -kel * central

    # 5. Observation and error.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
