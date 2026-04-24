Kuchimanchi_2018_evolocumab <- function() {
  description <- "One-compartment population PK model for evolocumab with first-order SC absorption and parallel linear plus Michaelis-Menten (target-mediated) elimination from the central compartment, in healthy adults and patients with hypercholesterolemia (Kuchimanchi 2018)"
  reference <- "Kuchimanchi M, Monine M, Kandadi Muralidharan K, Woodhead JL, Horner TJ. Population pharmacokinetics and exposure-response modeling and simulation for evolocumab in healthy volunteers and patients with hypercholesterolemia. J Pharmacokinet Pharmacodyn. 2019;46(2):133-148. doi:10.1007/s10928-018-9592-y"
  vignette <- "Kuchimanchi_2018_evolocumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate (WT/84)^exponent on CL (0.276), V (1.04), and Vmax (0.145). Reference 84 kg = mean body weight of the pooled phase 1-3 analysis population (Kuchimanchi 2018 Table 2 and Methods, reference-patient definition).",
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Female-sex exponent on V (1.11) — multiplicative factor of 1.11 on V for female subjects (Kuchimanchi 2018 Table 3). Reference patient is male (Methods, exposure-response reference patient).",
      source_name        = "SEXF"
    ),
    STATIN_MONO = list(
      description        = "Concomitant statin monotherapy indicator, 1 = patient on a statin only (no other lipid-lowering comedication), 0 = otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not on statin monotherapy)",
      notes              = "Kuchimanchi 2018 defines the statin covariate narrowly as 'patients on a statin only and no other comedication' (Methods, PopPK analysis). Multiplicative exponent 1.13 on Vmax (Table 3). Reference patient is not on any lipid-lowering medication. Mutually compatible with EZE (a patient can be 1 on STATIN_MONO xor 1 on EZE).",
      source_name        = "STATIN_MONO"
    ),
    EZE = list(
      description        = "Concomitant ezetimibe indicator, 1 = patient taking ezetimibe (with or without other lipid-lowering comedication), 0 = not on ezetimibe",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not on ezetimibe)",
      notes              = "Kuchimanchi 2018 defines the ezetimibe covariate as 'all patients on ezetimibe, regardless of comedications' (Methods, PopPK analysis). In the popPK dataset ~79% of ezetimibe users were also on a statin, so the effect effectively captures statin+ezetimibe combination therapy (the paper notates the exponent as 'Statin + ezetimibe' in Table 3). Multiplicative exponent 1.20 on Vmax (Table 3). Reference patient is not on ezetimibe.",
      source_name        = "EZE"
    ),
    PCSK9 = list(
      description        = "Baseline unbound PCSK9 (proprotein convertase subtilisin/kexin type 9) serum concentration",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate (PCSK9/425)^0.194 on Vmax (Kuchimanchi 2018 Table 3). Reference 425 ng/mL (= 5.9 nM) is the population median used for the reference patient (Methods, exposure-response reference patient). PCSK9 is reported in Kuchimanchi 2018 Table 2 in ng/mL; the nM equivalent uses a PCSK9 molecular weight of ~72 kDa. Baseline (time-fixed) covariate; patients with missing baseline PCSK9 were excluded from analyses that included PCSK9 as a covariate.",
      source_name        = "PCSK9"
    )
  )

  population <- list(
    n_subjects       = 3414L,
    n_observations   = 16179L,
    n_studies        = 11L,
    age_range        = "18-80 years",
    age_median       = "57 years (mean; Table 2 reports SD 58 which appears to be a typographical error)",
    weight_range     = "41-175 kg",
    weight_median    = "84.2 kg (mean)",
    sex_female_pct   = 50,
    race_ethnicity   = c(White = 87, Black = 7, Asian = 4, Hispanic = 0, Other = 1, AmericanIndianAlaska = 0, NativeHawaiianPacific = 0, Multiple = 0),
    disease_state    = "Pooled adults: healthy volunteers (phase 1a) and patients with hypercholesterolemia (phase 1b, 2, and 3), including patients with heterozygous familial hypercholesterolemia (9%), diabetes (11%), and statin-intolerance cohorts. Most subjects received concomitant lipid-lowering therapy (statins 72%, ezetimibe 12%).",
    dose_range       = "Evolocumab 7-420 mg IV or SC, single- and multiple-dose across Q2W and QM regimens. Phase 3 studies used the commercial regimens 140 mg SC Q2W and 420 mg SC QM.",
    regions          = "Multi-regional (11 pooled clinical studies spanning phase 1, 2, and 3).",
    pcsk9_baseline   = "Mean 402 ng/mL (SD 375), range 15.5-1233 ng/mL; median used for reference patient = 425 ng/mL (= 5.9 nM).",
    notes            = "Baseline characteristics from Kuchimanchi 2018 Table 2 (phase 1, 2, and 3 pooled column; N = 3414). Of the 5474 patients contributing data, 3414 received evolocumab and were included in the final popPK analysis; 1312 from 4 phase 2 studies were included in the exposure-response analysis (a separate Emax model on LDL-C; not packaged here — nlmixr2lib focuses on the popPK model)."
  )

  ini({
    # ---- Structural PK parameters (Table 3, updated phase 3 popPK model) ----
    # Reference body weight = 84 kg (paper Methods: 'an 84 kg male patient ... was considered as the reference patient').
    lka     <- fixed(log(0.319)); label("Absorption rate constant ka (1/day)")                       # Table 3: ka = 0.319 day^-1 (FIXED from phase 1a dense-data modeling)
    lcl     <- log(0.105);        label("Linear clearance CL for an 84 kg reference patient (L/day)") # Table 3: CL = 0.105 L/day (updated phase 3 popPK model)
    lvc     <- log(5.18);         label("Central volume of distribution V for an 84 kg male reference patient (L)") # Table 3: V = 5.18 L
    lfdepot <- fixed(log(0.72));  label("Subcutaneous bioavailability F (fraction; IV reference)")    # Table 3: F = 0.72 (FIXED; SC absolute bioavailability relative to IV)

    # ---- Michaelis-Menten (target-mediated) elimination from central ----
    # Paper parameterization: dC/dt-like rate = Vmax * C / (km + C), with both
    # Vmax and km reported in target-concentration units (nM). In the rxode2
    # mass-based ODE we multiply by V to get mass/day (see model() block).
    lvmax   <- fixed(log(9.85));  label("Nonlinear clearance capacity Vmax (nM/day; concentration-per-time form)") # Table 3: Vmax = 9.85 nM/day (FIXED in updated phase 3 popPK model)
    lkm     <- fixed(log(27.3));  label("Michaelis-Menten constant km (nM)")                          # Table 3: km = 27.3 nM (FIXED)

    # ---- Covariate exponents (Table 3) ----
    e_wt_cl     <- 0.276;  label("Power exponent of WT/84 on CL (unitless)")                          # Table 3: Body weight exponent on CL = 0.276 (RSE 30.4%)
    e_wt_v      <- 1.04;   label("Power exponent of WT/84 on V (unitless)")                           # Table 3: Body weight exponent on V = 1.04 (RSE 4.05%)
    e_sexf_v    <- 1.11;   label("Female-vs-male multiplicative exponent on V (unitless)")            # Table 3: Female exponent on V = 1.11 (RSE 1.42%)
    e_wt_vmax   <- 0.145;  label("Power exponent of WT/84 on Vmax (unitless)")                        # Table 3: Body weight exponent on Vmax = 0.145 (RSE 33.0%)
    e_smono_vmax <- 1.13;  label("Statin-monotherapy multiplicative exponent on Vmax (unitless)")      # Table 3: Statin exponent on Vmax = 1.13 (RSE 1.02%)
    e_eze_vmax  <- 1.20;   label("Ezetimibe (combination-therapy) multiplicative exponent on Vmax (unitless)") # Table 3: Statin + ezetimibe exponent on Vmax = 1.20 (RSE 1.59%)
    e_pcsk9_vmax <- 0.194; label("Power exponent of PCSK9/425 (ng/mL) on Vmax (unitless)")             # Table 3: PCSK9 baseline exponent on Vmax = 0.194 (RSE 7.47%)

    # ---- IIV (Table 3) ----
    # Paper reports between-subject variability as CV% for each parameter and
    # states that CL, V, and Vmax share a full-block variance matrix, with
    # independent diagonal IIV on ka and no IIV on km. The text mentions 'high
    # correlation between individual random effects on CL and Vmax (Table 4)',
    # but Table 4 in the paper is the exposure-response model and does not list
    # the PK omega-block correlations; they are not published. Consequently,
    # this file approximates the block matrix with independent diagonal IIV on
    # CL, V, and Vmax. This is documented as a deviation in the vignette's
    # Assumptions and deviations section.
    #   omega^2 = log(CV^2 + 1) for log-normal parameters:
    #   CV 54.3% -> log(0.543^2 + 1) = 0.25839 (CL)
    #   CV 28.3% -> log(0.283^2 + 1) = 0.07704 (V)
    #   CV 31.1% -> log(0.311^2 + 1) = 0.09232 (Vmax)
    #   CV 74.6% -> log(0.746^2 + 1) = 0.44245 (ka; FIXED in the paper)
    etalcl   ~ 0.25839                         # Table 3 IIV CV 54.3% on CL
    etalvc   ~ 0.07704                         # Table 3 IIV CV 28.3% on V
    etalvmax ~ 0.09232                         # Table 3 IIV CV 31.1% on Vmax
    etalka   ~ fix(0.44245)                    # Table 3 IIV CV 74.6% FIXED on ka

    # ---- Residual error (Table 3) ----
    # Paper reports proportional error as a fraction (0.282 = 28.2% CV) and
    # additive error in nM (5.41 nM). The additive residual is converted to
    # ug/mL inside model() using the evolocumab molecular weight of 141,800
    # g/mol (MW_EVO) so the residual applies directly to Cc (ug/mL).
    propSd    <- 0.282;  label("Proportional residual error (fraction)")                              # Table 3: Residual proportional error = 0.282 (RSE 1.12%)
    addSd_nM  <- 5.41;   label("Additive residual error (nM; converted to ug/mL inside model)")       # Table 3: Residual additive error = 5.41 nM (RSE 2.50%)
  })

  model({
    # ---- Physical constants ----
    # Evolocumab is a fully human IgG2 monoclonal antibody; molecular weight
    # 141,800 g/mol per the FDA-approved Repatha prescribing information.
    # Used to convert nM (paper units for km, Vmax, additive residual) to
    # ug/mL (the rxode2 concentration scale with dose in mg and V in L).
    MW_EVO  <- 141800            # Evolocumab molecular weight, g/mol
    nM_per_ugmL <- 1e6 / MW_EVO  # 1 ug/mL = 1e6 ng/L / MW (g/mol) nM; = 7.052 nM/(ug/mL)

    # ---- Individual PK parameters ----
    # Reference body weight 84 kg; reference patient is male, not on lipid-
    # lowering medication, with baseline PCSK9 = 425 ng/mL (paper Methods).
    ka    <- exp(lka + etalka)
    cl    <- exp(lcl + etalcl)   * (WT / 84)^e_wt_cl
    vc    <- exp(lvc + etalvc)   * (WT / 84)^e_wt_v  * e_sexf_v^SEXF
    vmax  <- exp(lvmax + etalvmax) * (WT / 84)^e_wt_vmax *
             e_smono_vmax^STATIN_MONO * e_eze_vmax^EZE *
             (PCSK9 / 425)^e_pcsk9_vmax   # Vmax in nM/day
    km    <- exp(lkm)            # km in nM
    fdepot <- exp(lfdepot)

    # ---- Unit bridging between nM (paper) and ug/mL (model state) ----
    # central amount is in mg; vc in L; so Cc = central/vc is in mg/L = ug/mL.
    # The MM arithmetic uses the paper's nM scale: Cc_nM = Cc (ug/mL) * 1e6 / MW.
    Cc      <- central / vc                                 # ug/mL
    Cc_nM   <- Cc * nM_per_ugmL                              # nM
    mm_rate_nM <- vmax * Cc_nM / (km + Cc_nM)                # nM/day (concentration-per-time)
    # Convert the nonlinear elimination rate back to mass/time for the ODE:
    # V * (nM/day) * (MW / 1e6 (ug/nmol)) = V [L] * (nM/day) * (ug/mL / nM) * (1 mL/0.001 L)  ->  mg/day.
    mm_rate_mgday <- vc * mm_rate_nM / nM_per_ugmL           # ug/mL/day * L = mg/L/day * L = mg/day

    # ---- ODE system (Figure 1a of Kuchimanchi 2018) ----
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - cl * Cc - mm_rate_mgday

    # ---- Bioavailability on SC depot dosing (IV doses go directly to central) ----
    f(depot) <- fdepot

    # ---- Observation and residual error ----
    # Convert additive residual from nM (paper units) to ug/mL (observation scale).
    addSd <- addSd_nM / nM_per_ugmL                          # ug/mL
    Cc ~ add(addSd) + prop(propSd)
  })
}
