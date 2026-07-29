Takama_2007_darbepoetin <- function() {
  description <- "Two-compartment intravenous population PK model for darbepoetin alfa in Japanese adult haemodialysis (HD) and peritoneal dialysis (PD) patients with an additive endogenous erythropoietin baseline concentration (Takama 2007). Body weight enters as a linear-deviation effect (centred on 54 kg) on clearance and central volume; peritoneal-dialysis modality adds a +17% multiplicative increment to central volume relative to the HD reference."
  reference <- paste(
    "Takama H, Tanaka H, Nakashima D, Ogata H, Uchida E, Akizawa T, Koshikawa S.",
    "Population pharmacokinetics of darbepoetin alfa in haemodialysis and",
    "peritoneal dialysis patients after intravenous administration.",
    "Br J Clin Pharmacol. 2007;63(3):300-309.",
    "doi:10.1111/j.1365-2125.2006.02756.x",
    sep = " "
  )
  vignette <- "Takama_2007_darbepoetin"
  units <- list(
    time          = "h",
    dosing        = "ug",
    concentration = "ng/mL"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight at study entry",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Takama 2007 Table 1 patient characteristics: median 54.7 kg (range 35.5-132.0). Reference value 54 kg (population median; covariate-deviation centre used in Methods equation 'CL = theta_CL * [1 + theta_CL_WT * (WT - 54)]'). Linear-deviation effect on both clearance (theta_CL_WT = 0.0195 kg^-1) and central volume (theta_V1_WT = 0.0163 kg^-1).",
      source_name        = "WT"
    ),
    PERIT_DIAL = list(
      description        = "Peritoneal-dialysis modality indicator (1 = PD, 0 = HD)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (intermittent haemodialysis, HD)",
      notes              = "Takama 2007 Methods and Table 4: dichotomous covariate DIA with DIA = 0 for HD and DIA = 1 for PD (63 HD / 68 PD subjects). Enters the structural model as an additive deviation term on central volume: V1 = theta_V1 * [1 + theta_V1_WT * (WT - 54) + theta_V1_DIA * DIA], with theta_V1_DIA = 0.170 so V1 in PD subjects is 17% higher than in HD subjects at the reference 54 kg weight. The paper attributes the V1 difference to higher extracellular body fluid in PD than post-HD intracorporeal fluid (Discussion, citing reference [22] of Takama 2007).",
      source_name        = "DIA"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 131L,
    n_studies        = 4L,
    age_range        = "23-84 years",
    age_median       = "60 years",
    weight_range     = "35.5-132.0 kg",
    weight_median    = "54.7 kg",
    sex_female_pct   = 38.2,
    race_ethnicity   = c(Japanese = 100),
    disease_state    = "Japanese adult chronic kidney disease patients on dialysis: 63 haemodialysis (HD) and 68 peritoneal dialysis (PD). HD inclusion: chronic kidney disease undergoing HD for >=3 months with stable IV rHuEPO three times weekly for >=2 months before the study, mean baseline haemoglobin 9-12 g/dL (haematocrit 25-35%). PD inclusion: end-stage renal failure on stable PD with mean baseline haemoglobin <10 g/dL (no rHuEPO) or 9-12 g/dL (on rHuEPO). Both groups excluded major surgery within 3 months (excluding vascular access surgery) and blood transfusions within 1 month before enrolment.",
    dose_range       = "10-90 ug intravenous, single or multiple administration every 1 or 2 weeks. HD subjects: 10-60 ug single doses with serial sampling at predose and 0.5, 1, 2, 5, 8, 12, 24, 48, 96, 168 h (90 ug dose group additionally sampled at 336 h). PD subjects: dosed every 1 or 2 weeks with sparse sampling predose in weeks 1-4 plus optional 0.5-1 h post-dose in week 2.",
    regions          = "Japan",
    n_concentrations = 917L,
    renal_function   = "Anuric / dialysis-dependent (chronic dialysis cohort). Baseline serum creatinine median 10.60 mg/dL (range 3.68-16.50); serum albumin median 3.8 g/dL (range 2.1-4.8); red blood cell counts median 317 x 10^4 /uL (range 211-422); white blood cell counts median 5500 /uL (range 2400-12500); platelet cell counts median 17.7 x 10^4 /uL (range 5.1-50.1).",
    notes            = "Four clinical studies pooled (three HD studies + one PD study). Serum darbepoetin alfa concentrations were quantified with the Quantikine in-vitro-diagnostic rHuEPO ELISA kit (R&D Systems) using a darbepoetin alfa standard curve; the assay limit of quantification was 0.078 ng/mL and inter-assay precision was within 20% for clinical samples. Baseline endogenous EPO concentration in the cohort (theta_k0) was 0.167 ng/mL (Table 4) and is retained in the model as an additive offset on the observed serum concentration. Bootstrap validation: the final model was refitted to 1000 bootstrap replicates; mean parameter estimates were close to those of the original dataset (Table 4 right panel). Initial values and a one-compartment alternative were rejected by Akaike's information criterion (Table 2)."
  )

  ini({
    # Structural parameters (Takama 2007 Table 4, "Typical values" column from the
    # original dataset; reference subject is a 54-kg HD patient). All values are
    # NONMEM final estimates (not initial values) with the standard-error and 95%
    # confidence-interval columns reported in Table 4 alongside the bootstrap
    # check. The paper parameterises clearance and volumes on the linear scale
    # (CL_i = theta_CL * exp(eta_i)); we store on the log scale here to match the
    # nlmixr2lib convention.
    lcl <- log(0.0807);   label("Clearance CL at WT = 54 kg (L/h)")                              # Takama 2007 Table 4: theta_CL = 0.0807 L/h (SE 0.00195)
    lq  <- log(0.0616);   label("Intercompartmental clearance Q (L/h)")                           # Takama 2007 Table 4: theta_Q = 0.0616 L/h (SE 0.00755)
    lvc <- log(2.51);     label("Central volume V1 at WT = 54 kg, HD (L)")                        # Takama 2007 Table 4: theta_V1 = 2.51 L (SE 0.0629)
    lvp <- log(0.522);    label("Peripheral volume V2 (L)")                                       # Takama 2007 Table 4: theta_V2 = 0.522 L (SE 0.0407)
    lrbase <- log(0.167); label("Baseline endogenous erythropoietin concentration k0 (ng/mL)")    # Takama 2007 Table 4: theta_k0 = 0.167 ng/mL (SE 0.00629)

    # Covariate effects. Takama 2007 uses linear-deviation centring:
    #   CL = theta_CL * [1 + theta_CL_WT * (WT - 54)]
    #   V1 = theta_V1 * [1 + theta_V1_WT * (WT - 54) + theta_V1_DIA * DIA]
    # All three coefficients are estimated (Table 4 reports SE and 95% CI for
    # each); none are held fixed by the paper.
    e_wt_cl         <- 0.0195; label("Linear-deviation coefficient of (WT - 54) on CL (1/kg)")              # Takama 2007 Table 4: theta_CL_WT = 0.0195 kg^-1 (SE 0.00246; 95% CI 0.0147-0.0243)
    e_wt_vc         <- 0.0163; label("Linear-deviation coefficient of (WT - 54) on V1 (1/kg)")              # Takama 2007 Table 4: theta_V1_WT = 0.0163 kg^-1 (SE 0.00214; 95% CI 0.0121-0.0205)
    e_perit_dial_vc <- 0.170;  label("Additive PD-modality effect on V1 (fractional increment; PERIT_DIAL = 1 for PD)")  # Takama 2007 Table 4: theta_V1_DIA = 0.170 (SE 0.0501; 95% CI 0.0718-0.268)

    # Inter-individual variability. Methods describe an exponential model
    # (CL_i = theta_CL * exp(eta_i), eta_i ~ N(0, omega^2)) with all five
    # parameters (CL, Q, V1, V2, k0) carrying IIV. Table 4 reports each IIV
    # as percent CV; conversion to the log-normal variance scale uses
    # omega^2 = log(CV^2 + 1). The paper does not report a covariance matrix
    # between etas, so the IIVs are encoded as a diagonal block (independent
    # etas) matching the standard NONMEM diagonal $OMEGA convention.
    etalcl    ~ 0.04153  # log(0.206^2 + 1) = 0.04153; Takama 2007 Table 4: omega_CL = 20.6%
    etalq     ~ 0.56449  # log(0.871^2 + 1) = 0.56449; Takama 2007 Table 4: omega_Q  = 87.1%
    etalvc    ~ 0.04644  # log(0.218^2 + 1) = 0.04644; Takama 2007 Table 4: omega_V1 = 21.8%
    etalvp    ~ 0.21205  # log(0.486^2 + 1) = 0.21205; Takama 2007 Table 4: omega_V2 = 48.6%
    etalrbase ~ 0.18741  # log(0.454^2 + 1) = 0.18741; Takama 2007 Table 4: omega_k0 = 45.4%

    # Residual error: combination of additive and constant-coefficient-of-
    # variation (proportional) error (Takama 2007 Methods equation).
    # Cp_obs = Cp_pred * (1 + eps1) + eps2 with eps1, eps2 ~ N(0, sigma^2).
    # Table 4 reports sigma1 = 6.53% (proportional) and sigma2 = 0.0634 ng/mL
    # (additive). Both are estimated, not fixed.
    addSd  <- 0.0634; label("Additive residual error (ng/mL)")          # Takama 2007 Table 4: sigma_2 = 0.0634 ng/mL
    propSd <- 0.0653; label("Proportional residual error (fraction)")    # Takama 2007 Table 4: sigma_1 = 6.53%
  })

  model({
    # Reference covariate value used for the linear-deviation centring (54 kg
    # is the cohort weight median from Takama 2007 Table 1 / Methods).
    wt_ref <- 54

    # Individual PK parameters and the endogenous baseline. Linear-deviation
    # form preserved exactly as Takama 2007 wrote it; the form remains positive
    # across the cohort WT range (35.5-132 kg) but can become negative for WT
    # far below the cohort, so simulation weights should stay within the
    # study range. The authors did not reparameterise to a power form.
    cl    <- exp(lcl    + etalcl)    * (1 + e_wt_cl * (WT - wt_ref))
    vc    <- exp(lvc    + etalvc)    * (1 + e_wt_vc * (WT - wt_ref) + e_perit_dial_vc * PERIT_DIAL)
    q     <- exp(lq     + etalq)
    vp    <- exp(lvp    + etalvp)
    rbase <- exp(lrbase + etalrbase)

    # Micro-constants for the two-compartment ODE system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment open model with intravenous bolus dosing into central.
    # Endogenous EPO is treated as an additive constant offset on the observed
    # concentration (Methods: 'k0_j is the baseline erythropoietin level for
    # the j-th individual'; the observed serum concentration includes both the
    # drug-related signal and the baseline endogenous EPO).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # Total observed concentration: drug contribution plus endogenous baseline.
    Cc <- central / vc + rbase
    Cc ~ add(addSd) + prop(propSd)
  })
}
