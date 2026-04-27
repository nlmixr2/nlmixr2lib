Fau_2020_isatuximab <- function() {
  description <- "Two-compartment population PK model for intravenous isatuximab (anti-CD38 IgG1) in adults with relapsed/refractory multiple myeloma, with parallel time-varying linear and Michaelis-Menten eliminations from the central compartment (Fau 2020). The linear clearance follows a sigmoidal Emax decay from baseline to steady state; the magnitude of the decay differs by multiple-myeloma immunoglobulin type."
  reference <- "Fau JB, El-Cheikh R, Brillac C, et al. Drug-Disease Interaction and Time-Dependent Population Pharmacokinetics of Isatuximab in Relapsed/Refractory Multiple Myeloma Patients. CPT Pharmacometrics Syst Pharmacol. 2020;9(11):649-658. doi:10.1002/psp4.12561"
  vignette <- "Fau_2020_isatuximab"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on CLinf (steady-state linear clearance, exponent 0.621), Vc (0.472), Vp (0.719), and Q (0.477). Reference 75.6 kg (Fau 2020 Table S1 population median).",
      source_name        = "WT"
    ),
    B2M = list(
      description        = "Baseline serum beta-2-microglobulin",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on CLinf (exponent 0.343). Reference 3.90 mg/L (Fau 2020 Table S1 population median).",
      source_name        = "B2M"
    ),
    MM_NIGG = list(
      description        = "Multiple-myeloma immunoglobulin type indicator: 1 = non-IgG MM, 0 = IgG MM",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (IgG MM)",
      notes              = "Exponential effect on both CLinf (-0.751) and KCL (-0.931). Patients secreting non-IgG monoclonal protein have lower linear clearance and a faster transition to steady-state CL than IgG-secretors; the mechanism proposed by Fau 2020 is that endogenous IgG M-protein in IgG-MM patients competes with the therapeutic mAb for FcRn-mediated salvage, increasing isatuximab clearance. Source column 'Ig_type' takes the value 1 for non-IgG MM and 0 for IgG MM.",
      source_name        = "Ig_type"
    ),
    FORM_P2F2 = list(
      description        = "Drug-material indicator: 1 = P2F2 (phase III / commercial-bound material), 0 = P1F1 (early-phase material)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (P1F1)",
      notes              = "Exponential effect on Vc (-0.137). The P2F2 material was used in the phase III ICARIA-MM study (EFC14335) and is the intended commercial formulation; P1F1 was used in the earlier phase I/II/Ib studies. Set to 1 to simulate the marketed material.",
      source_name        = "Drug_mat"
    ),
    RACE_ASIAN = list(
      description        = "Race indicator: 1 = Asian, 0 = non-Asian",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian)",
      notes              = "Exponential effect on Vc (-0.275). Asians made up 5.3% (n=25/476) of the analysis population (Fau 2020 Table S2).",
      source_name        = "Asian"
    ),
    SEXF = list(
      description        = "Biological sex indicator: 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Exponential effect on Vc (-0.126). Females had ~12% lower Vc than males in the typical patient.",
      source_name        = "Sex"
    )
  )

  population <- list(
    n_subjects     = 476L,
    n_observations = 7697L,
    n_studies      = 4L,
    age_range      = "49.0-79.0 years (5th-95th percentile)",
    age_median     = "65.0 years",
    weight_range   = "51.3-110 kg (5th-95th percentile)",
    weight_median  = "75.6 kg",
    sex_female_pct = 43.5,
    race_ethnicity = c(Caucasian = 79.2, Black = 3.8, Asian = 5.3, Other_or_Missing = 11.8),
    disease_state  = "Relapsed / refractory multiple myeloma (RRMM); 55% IgG MM, 45% non-IgG MM (Fau 2020 Table S2). 76.7% with detectable serum M-protein at baseline. ECOG 0/1/>=2 = 31.9% / 58.2% / 9.9%. ISS stage I/II/III = 36.3% / 38.2% / 25.4%.",
    dose_range     = "Isatuximab 1-20 mg/kg as ~1-h IV infusion; 91% of subjects received the marketed 10 or 20 mg/kg doses. Phase III combination dose 10 mg/kg q4w (cycle 1, weekly x 4) then q2w with pomalidomide-dexamethasone. Some early-phase cohorts dosed q2w throughout.",
    regions        = "Multinational pooled across 4 trials: TED10893 (phase I/II monotherapy, n=258), TED14154 (phase I, n=26), TCD14079 (phase Ib combination, n=44), and EFC14335 / ICARIA-MM (phase III combination, n=148).",
    notes          = "Pooled analysis from 4 clinical trials supporting the FDA approval of isatuximab + pomalidomide + dexamethasone for RRMM (March 2020). Baseline demographics from Fau 2020 Tables S1-S2. Drug material: 40.1% P1F1, 59.9% P2F2 (commercial-bound). Median baseline beta-2-microglobulin 3.90 mg/L. Median baseline serum M-protein 18.0 g/L. eGFR-MDRD median 68.3 mL/min/1.73 m2. The model uses Monolix's combined additive + proportional residual-error parameterization; CLm is normally distributed while all other inter-individual variabilities are log-normal."
  )

  ini({
    # Structural typical-value parameters at the reference patient: 75.6 kg
    # body weight, 3.90 mg/L baseline B2M, IgG MM, P1F1 drug material,
    # non-Asian, male (Fau 2020 Table S3 final-model estimates).
    lclinf <- log(0.00955); label("Linear CLinf at infinite time at reference covariates (L/hour)") # Fau 2020 Table S3: CLinf
    lvc    <- log(5.13);    label("Central volume of distribution Vc at reference (L)")             # Fau 2020 Table S3: Vc
    lq     <- log(0.0432);  label("Intercompartmental clearance Q at reference (L/hour)")           # Fau 2020 Table S3: Q
    lvp    <- log(3.62);    label("Peripheral volume of distribution Vp at reference (L)")          # Fau 2020 Table S3: Vp

    # Parallel Michaelis-Menten (target-mediated approximation) elimination
    # from the central compartment. Vm has units of mg/(L*h) in the source
    # equation dL/dt = ... - Vm/(Km + L) * L. The amount-form ODE multiplies
    # by Vc, giving a mass-rate term Vc * Vm * Cc / (Km + Cc) in mg/h.
    lvm <- log(0.136); label("Michaelis-Menten Vm of nonlinear CL (mg/(L*h))") # Fau 2020 Table S3: Vm
    lkm <- log(0.300); label("Michaelis-Menten Km of nonlinear CL (mg/L)")     # Fau 2020 Table S3: Km

    # Time-varying linear CL via a sigmoidal Emax form (Fau 2020):
    #   CLlin(t) = CLinf * exp( CLm * (1 - t^gamma / (KCL^gamma + t^gamma)) )
    # CLm is the signed log-ratio log(CL_0 / CLinf) — for a typical patient
    # CLm = 0.664 implies a baseline CL exp(0.664) ~= 1.94x the steady-state
    # value, i.e., a 48% decrease toward steady state. CLm is normally
    # distributed (the only normal IIV in this model); KCL is the time at
    # which the sigmoid has progressed half-way; gamma is the stiffness
    # exponent.
    clm  <- 0.664;      label("Magnitude of CLlin change from baseline to steady state, log(CL_0/CLinf) (unitless)") # Fau 2020 Table S3: CLm
    lkcl <- log(1055);  label("Half-time scale KCL of the time-on-CLlin sigmoid (hour)")                              # Fau 2020 Table S3: KCL
    lgam <- log(3.91);  label("Sigmoidal stiffness exponent gamma of the time-on-CLlin function (unitless)")          # Fau 2020 Table S3: gamma

    # Allometric body-weight power exponents (reference 75.6 kg).
    e_wt_clinf <- 0.621; label("Power exponent of WT/75.6 on CLinf (unitless)") # Fau 2020 Table S3: CLinf ~ Wght
    e_wt_vc    <- 0.472; label("Power exponent of WT/75.6 on Vc (unitless)")    # Fau 2020 Table S3: Vc ~ Wght
    e_wt_vp    <- 0.719; label("Power exponent of WT/75.6 on Vp (unitless)")    # Fau 2020 Table S3: Vp ~ Wght
    e_wt_q     <- 0.477; label("Power exponent of WT/75.6 on Q (unitless)")     # Fau 2020 Table S3: Q ~ Wght (RSE 57.5%; least-precisely estimated coefficient)

    # B2M power covariate on CLinf.
    e_b2m_clinf <- 0.343; label("Power exponent of B2M/3.90 on CLinf (unitless)") # Fau 2020 Table S3: CLinf ~ B2M

    # Categorical covariate effects, applied as exp(coef * indicator).
    e_nigg_clinf <- -0.751; label("Exponential coefficient of non-IgG MM (vs IgG MM) on CLinf")     # Fau 2020 Table S3: CLinf ~ IgType=Not_IgG
    e_nigg_kcl   <- -0.931; label("Exponential coefficient of non-IgG MM (vs IgG MM) on KCL")       # Fau 2020 Table S3: KCL ~ IgType=Not_IgG (main text rounds to -0.930)
    e_p2f2_vc    <- -0.137; label("Exponential coefficient of P2F2 (vs P1F1) drug material on Vc")  # Fau 2020 Table S3: Vc ~ Formulation=P2F2
    e_asian_vc   <- -0.275; label("Exponential coefficient of Asian race (vs non-Asian) on Vc")     # Fau 2020 Table S3: Vc ~ Race=Asian
    e_sexf_vc    <- -0.126; label("Exponential coefficient of female sex (vs male) on Vc")          # Fau 2020 Table S3: Vc ~ Sex=Female

    # Inter-individual variability. Fau 2020 Table S3 footnote defines
    # "ω: Between-subject coefficient of variation"; convert to nlmixr2
    # variances using:
    #   log-normal:  omega^2 = log(CV^2 + 1)   (CLinf, KCL, gamma, Vc, Q, Vp, Vm, Km)
    #   normal CLm:  var     = (CLm * CV)^2    (additive eta on the linear scale)
    # CLinf 47.5% -> log(0.475^2 + 1) = 0.2035
    # KCL  115%   -> log(1.15^2 + 1)  = 0.8427
    # gam  118%   -> log(1.18^2 + 1)  = 0.8723
    # Vc    25.7% -> log(0.257^2 + 1) = 0.0640
    # Q     85.8% -> log(0.858^2 + 1) = 0.5519
    # Vp    45.6% -> log(0.456^2 + 1) = 0.1889
    # Vm    61.5% -> log(0.615^2 + 1) = 0.3206
    # Km    88.9% -> log(0.889^2 + 1) = 0.5823
    # CLm   97.2% -> (0.664 * 0.972)^2 = 0.4164
    # IIVs are reported independently (no off-block correlations in the source).
    etalclinf ~ 0.2035  # Fau 2020 Table S3: ω(CLinf) 47.5%
    etaclm    ~ 0.4164  # Fau 2020 Table S3: ω(CLm) 97.2% (additive normal eta)
    etalkcl   ~ 0.8427  # Fau 2020 Table S3: ω(KCL) 115%
    etalgam   ~ 0.8723  # Fau 2020 Table S3: ω(gamma) 118%
    etalvc    ~ 0.0640  # Fau 2020 Table S3: ω(Vc) 25.7%
    etalq     ~ 0.5519  # Fau 2020 Table S3: ω(Q) 85.8%
    etalvp    ~ 0.1889  # Fau 2020 Table S3: ω(Vp) 45.6%
    etalvm    ~ 0.3206  # Fau 2020 Table S3: ω(Vm) 61.5%
    etalkm    ~ 0.5823  # Fau 2020 Table S3: ω(Km) 88.9%

    # Combined additive + proportional residual error (Fau 2020 Table S3).
    propSd <- 0.225;   label("Proportional residual error (fraction)") # Fau 2020 Table S3: σprop 22.5%
    addSd  <- 0.00196; label("Additive residual error (mg/L)")         # Fau 2020 Table S3: σadd 0.00196 mg/L
  })
  model({
    # Individual time-zero structural parameters with covariate effects
    # (Fau 2020 final-model equations for CLinf, KCL, Vc, Vp, Q).
    clinf <- exp(lclinf + etalclinf) *
             (WT  / 75.6)^e_wt_clinf *
             (B2M / 3.90)^e_b2m_clinf *
             exp(e_nigg_clinf * MM_NIGG)

    kcl <- exp(lkcl + etalkcl) * exp(e_nigg_kcl * MM_NIGG)

    # CLm is normally distributed; etaclm is additive on the linear scale.
    clm_i <- clm + etaclm

    gam <- exp(lgam + etalgam)

    vc <- exp(lvc + etalvc) *
          (WT / 75.6)^e_wt_vc *
          exp(e_p2f2_vc  * FORM_P2F2) *
          exp(e_asian_vc * RACE_ASIAN) *
          exp(e_sexf_vc  * SEXF)

    vp <- exp(lvp + etalvp) * (WT / 75.6)^e_wt_vp
    q  <- exp(lq  + etalq)  * (WT / 75.6)^e_wt_q

    vm <- exp(lvm + etalvm)
    km <- exp(lkm + etalkm)

    # Time-varying linear CL: starts at CLinf*exp(CLm) and decays to CLinf
    # as t -> infinity. KCL is the time at which the sigmoid is half-way
    # through (when gamma is large enough that the inflection is sharp).
    cllin <- clinf * exp(clm_i * (1 - t^gam / (kcl^gam + t^gam)))

    # Two-compartment IV-input PK with parallel time-varying linear and
    # Michaelis-Menten eliminations from the central compartment.
    # The MM term in Fau 2020's concentration-form ODE (-Vm/(Km+L) * L)
    # converts to amount form as -Vc * Vm * Cc / (Km + Cc) [mg/h], where
    # Cc = central / Vc is the central concentration.
    Cc <- central / vc

    d/dt(central)     <- -(cllin / vc) * central -
                          vc * vm * Cc / (km + Cc) -
                          (q / vc) * central +
                          (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central -
                          (q / vp) * peripheral1

    Cc ~ add(addSd) + prop(propSd)
  })
}
