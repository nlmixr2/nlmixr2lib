Bukkems_2021_raltegravir <- function() {
  description <- paste(
    "Two-compartment population PK model for oral raltegravir in a pooled",
    "cohort of 221 adults (healthy volunteers, non-pregnant adults living",
    "with HIV, and pregnant women living with HIV) across the 400 mg BID,",
    "800 mg QD, and 1200 mg QD (two 600 mg tablets) regimens (Bukkems 2021).",
    "Absorption is modelled as a chain of four sequential first-order",
    "compartments: depot -> transit1 -> transit2 -> transit3 -> central,",
    "with the first three transitions governed by a shared rate constant",
    "ktr = 3 / mat (paper's mean transit time MAT parameter) and the final",
    "transit3 -> central transition governed by a separate first-order",
    "absorption rate constant ka. Disposition is a two-compartment linear",
    "model with apparent central and peripheral volumes and apparent",
    "clearance and inter-compartmental clearance. Body weight enters as",
    "allometric scaling with fixed exponents 0.75 on CL and Q and 1.0 on",
    "Vc and Vp, referenced to 70 kg. Six covariate-parameter relationships",
    "are retained in the final model: FED on MAT (+160% with any food),",
    "CONMED_ATAZANAVIR on CL (-17%), FORM_RAL_600 on F (+21%),",
    "FED_LOWFAT on F (-46%), PREG on F (-49%), and CONMED_EFV on F (-17%).",
    "Inter-individual variability is a diagonal eta on Vc plus a correlated",
    "3-parameter block on CL, Q, and Vp with the CL-Vp off-diagonal fixed",
    "to zero. Residual error is a proportional model with a time-varying",
    "magnitude that switches at 3 h after dose (43.5% CV before, 29.0% CV",
    "after, tracking the paper's empirical time-varying residual). The",
    "paper's inter-occasion variability on F and MAT and its additional",
    "eta on the residual magnitude are omitted from this packaged model;",
    "see the validation vignette Assumptions and deviations section.",
    sep = " "
  )
  reference <- paste(
    "Bukkems VE, Post TM, Colbers AP, Burger DM, Svensson EM.",
    "A population pharmacokinetics analysis assessing the exposure of",
    "raltegravir once-daily 1200 mg in pregnant women living with HIV.",
    "CPT Pharmacometrics Syst Pharmacol. 2021;10(2):161-172.",
    "doi:10.1002/psp4.12586.",
    sep = " "
  )
  vignette <- "Bukkems_2021_raltegravir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on flow parameters (CL, Q) with fixed exponent 0.75 and on volume parameters (Vc, Vp) with fixed exponent 1.0, referenced to 70 kg (Bukkems 2021 Methods 'Development population pharmacokinetic model' paragraph 4). For pregnant women, postpartum weight is used because allometric scaling in pregnant women has not been established and would confound the pregnancy effect; when postpartum weight is unavailable (n = 4 subjects in the source), the source paper imputes it as third-trimester weight * 0.93 (the mean 7% decrease from third trimester to postpartum).",
      source_name        = "WT"
    ),
    PREG = list(
      description        = "Pregnancy status indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-pregnant)",
      notes              = "Time-fixed per subject. Applied as a linear additive effect on typical bioavailability F: F *= (1 + e_preg_fdepot * PREG) with e_preg_fdepot = -0.487 (Bukkems 2021 Table 2 'Factor change in F pregnancy'). Pregnant subjects (all in the third trimester at the intensive PK sampling, gestational age ~33 weeks) are 22 women from the PANNA study using the 400 mg BID regimen; the paper's sensitivity analysis (Supporting Information S2) showed near-identical primary-endpoint conclusions whether the pregnancy effect was placed on F (dOFV -17.82) or on CL (dOFV -12.62), and the F parameterisation was retained because it was the more significant of the two.",
      source_name        = "PREG"
    ),
    CONMED_ATAZANAVIR = list(
      description        = "Concomitant atazanavir coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant atazanavir)",
      notes              = "Time-fixed per subject within the analysis window. Applied as a linear additive effect on apparent clearance: CL *= (1 + e_atazanavir_cl * CONMED_ATAZANAVIR) with e_atazanavir_cl = -0.17 (Bukkems 2021 Table 2 'Factor change in CL with atazanavir'). Atazanavir-associated raltegravir CL is 17% lower than non-atazanavir reference, consistent with atazanavir-mediated UGT1A1 inhibition of raltegravir glucuronidation. In the pooled dataset the covariate primarily flags the 18 HIV-infected non-pregnant adults on a 400 mg BID atazanavir-based regimen who are otherwise indistinguishable from the healthy 400 mg BID reference cohort.",
      source_name        = "ATV"
    ),
    CONMED_EFV = list(
      description        = "Concomitant efavirenz coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant efavirenz)",
      notes              = "Time-fixed per subject within the analysis window. Applied as a linear additive effect on typical bioavailability F: F *= (1 + e_efv_fdepot * CONMED_EFV) with e_efv_fdepot = -0.167 (Bukkems 2021 Table 2 'Factor change in F efavirenz co-administration'). Efavirenz-associated raltegravir F is 17% lower than the no-efavirenz reference, consistent with efavirenz-mediated UGT1A1 induction (opposite direction from the atazanavir CL effect). The efavirenz-cohort raltegravir samples come from the 1200 mg QD sub-study reported as Bukkems 2021 reference 20.",
      source_name        = "EFV"
    ),
    FORM_RAL_600 = list(
      description        = "Raltegravir 600 mg film-coated tablet formulation indicator (Isentress HD)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (400 mg film-coated tablet, Isentress)",
      notes              = "Per-dose-record indicator. Applied as TWO effects: (1) linear additive on typical bioavailability F: F *= (1 + e_ral_600_fdepot * FORM_RAL_600) with e_ral_600_fdepot = 0.209 (Bukkems 2021 Table 2 'Factor change in F 600 mg formulation'); (2) multiplicative on the inter-occasion variability magnitude on F, which reduces the IOV F standard deviation by 72% for 600 mg records vs 400 mg records (this file omits IOV per convention). The 600 mg tablet is the polymer-based reformulation that disintegrates and dissolves faster than the original 400 mg tablet at physiological gastric pH; the two-tablet 1200 mg QD regimen using the 600 mg tablet is marketed as Isentress HD.",
      source_name        = "NEW"
    ),
    FED = list(
      description        = "Fed (any food) vs fasted dose-record indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = "Per-dose-record indicator, irrespective of meal type (low-fat, moderate-fat, high-fat all coded as FED = 1). Applied as a linear additive effect on the mean transit time MAT: MAT *= (1 + e_fed_mat * FED) with e_fed_mat = 1.6 (Bukkems 2021 Table 2 'Factor change in MAT fed'); MAT approximately triples with any food (0.336 h fasted -> 0.874 h fed), delaying absorption. The paper found the moderate-fat-meal effect could not be tested independently (individual meal-type data unavailable) and assumed the same MAT delay for any-food conditions based on prior raltegravir food-effect data.",
      source_name        = "FED"
    ),
    FED_LOWFAT = list(
      description        = "Low-fat meal (389 kcal, 6.9% fat) dose-record indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted, moderate-fat, or high-fat meal)",
      notes              = "Per-dose-record indicator specific to the low-fat meal condition. Applied as a linear additive effect on typical bioavailability F: F *= (1 + e_lowfat_fdepot * FED_LOWFAT) with e_lowfat_fdepot = -0.459 (Bukkems 2021 Table 2 'Factor change in F low-fat meal'); low-fat-meal raltegravir F is 46% lower than the fasted / moderate-fat reference. The low-fat definition (389 kcal, 6.9% fat) is inherited from Rizk et al. 2012 (Bukkems 2021 reference 8) which studied the low-fat food effect for both the 400 mg and 600 mg raltegravir tablets under identical protocols; the paper assumes the same low-fat food effect applies to both formulations.",
      source_name        = "LOWFAT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 221L,
    n_studies      = 11L,
    n_observations = 4016L,
    age_range      = "18-75 years (pooled across the 11 constituent studies)",
    weight_range   = "43-111 kg (pooled)",
    sex_female_pct = 44,
    disease_state  = paste(
      "Pooled cohort of (i) healthy adult volunteers taking raltegravir 400",
      "mg BID or 1200 mg QD (two 600 mg tablets), (ii) HIV-infected adults",
      "taking 400 mg BID (with concomitant atazanavir) or 1200 mg QD, and",
      "(iii) 22 HIV-infected pregnant women in the third trimester (target",
      "gestational age 33 weeks) taking 400 mg BID. See Bukkems 2021 Table 1",
      "for the 11 constituent studies (references 8, 10, 14-21).",
      sep = " "
    ),
    dose_range     = paste(
      "Raltegravir 400 mg BID (six studies), 800 mg QD (one study),",
      "and 1200 mg QD (five studies; the 1200 mg dose is administered as",
      "two 600 mg tablets). Both single-dose and steady-state sampling",
      "protocols are represented; intensive sampling schedules with",
      "typically 10-14 samples per curve over a 12- or 24-h post-dose",
      "window (up to 263 h for single-dose studies with terminal-phase",
      "sampling).",
      sep = " "
    ),
    regions        = "Europe (Netherlands PANNA multicentre network, other European sites, and healthy-volunteer studies).",
    pregnancy_cohort = paste(
      "22 European HIV-infected pregnant women receiving raltegravir 400 mg",
      "BID as part of a boosted or unboosted combination ART regimen.",
      "Intensive PK sampling in the third trimester (preferably 33 weeks",
      "gestation) and 4-6 weeks postpartum. Two additional women taking",
      "1200 mg QD raltegravir were sampled as clinical cases for the",
      "simulation-vs-observed comparison but were not included in the",
      "model-building dataset.",
      sep = " "
    ),
    notes          = paste(
      "Baseline demographics summarised from Bukkems 2021 Table 1 across",
      "the 11 constituent studies. Data was pooled from 226 individuals",
      "and 5772 sampling points; after exclusions for interacting",
      "comedications (1164 samples), imputable and non-imputable BLOQ",
      "handling (99 BLOQ excluded, 70 BLOQ imputed as LLOQ/2), and",
      "non-evaluable samples (493), the final PopPK model was built with",
      "221 individuals and 4016 sampling points. Model fit with NONMEM 7.4",
      "using FOCE with eta-epsilon interaction (Pirana 2.9.7 interface).",
      sep = " "
    )
  )

  ini({
    # ---- Structural fixed effects (Bukkems 2021 Table 2 'Parameter estimate' column, at 70 kg reference) ----
    lka  <- log(0.741); label("First-order absorption rate constant ka from transit3 to central (1/h)")            # Table 2: Ka = 0.741 h^-1 (RSE 2%)
    lmat <- log(0.336); label("Mean transit time MAT through the three ktr-governed transit steps at fasted reference (h)") # Table 2: MAT = 0.336 h fasted (RSE 8%)
    lvc  <- log(44.3);  label("Apparent central volume of distribution Vc/F at 70 kg (L)")                          # Table 2: V_c/F = 44.3 L at 70 kg (RSE 7%)
    lcl  <- log(55.8);  label("Apparent oral clearance CL/F at 70 kg (L/h)")                                        # Table 2: CL/F = 55.8 L/h at 70 kg (RSE 5%)
    lq   <- log(5.68);  label("Apparent inter-compartmental clearance Q/F at 70 kg (L/h)")                          # Table 2: Q/F = 5.68 L/h at 70 kg (RSE 7%)
    lvp  <- log(92.8);  label("Apparent peripheral volume of distribution Vp/F at 70 kg (L)")                       # Table 2: V_p/F = 92.8 L at 70 kg (RSE 9%)
    lfdepot <- fixed(log(1)); label("Typical bioavailability F at reference covariates (unitless)")     # Table 2: F = 1 FIXED (no IV data; absolute F not identifiable)

    # ---- Allometric exponents (Bukkems 2021 Methods 'Development population pharmacokinetic model' paragraph 4, fixed) ----
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL and Q (unitless)")                              # Methods p4: fixed allometric exponent 0.75 for flow parameters
    e_wt_vc <- fixed(1.0);  label("Allometric exponent on Vc and Vp (unitless)")                             # Methods p4: fixed allometric exponent 1 for volume parameters

    # ---- Covariate effects on structural parameters (Bukkems 2021 Table 2 covariate rows) ----
    # All applied via the linear additive form: parameter *= (1 + coefficient * indicator).
    # Table 2 footnote (a) makes the parameterisation explicit for each of the six retained
    # covariate-parameter relationships.
    e_fed_mat            <-  1.6;    label("Factor change in MAT with any food (fraction; +160% with FED = 1)")                                 # Table 2: Factor change in MAT fed = 1.6 (RSE 19%)
    e_atazanavir_cl      <- -0.17;   label("Factor change in CL with concomitant atazanavir (fraction; -17% with CONMED_ATAZANAVIR = 1)")       # Table 2: Factor change in CL with atazanavir = -0.17 (RSE 25%)
    e_ral_600_fdepot     <-  0.209;  label("Factor change in F for the 600 mg tablet (fraction; +21% with FORM_RAL_600 = 1)")                   # Table 2: Factor change in F 600 mg formulation = 0.209 (RSE 26%)
    e_lowfat_fdepot      <- -0.459;  label("Factor change in F for a low-fat meal (fraction; -46% with FED_LOWFAT = 1)")                        # Table 2: Factor change in F low-fat meal = -0.459 (RSE 9%)
    e_preg_fdepot        <- -0.487;  label("Factor change in F for pregnancy (fraction; -49% with PREG = 1)")                                   # Table 2: Factor change in F pregnancy = -0.487 (RSE 14%)
    e_efv_fdepot         <- -0.167;  label("Factor change in F with concomitant efavirenz (fraction; -17% with CONMED_EFV = 1)")                # Table 2: Factor change in F efavirenz co-administration = -0.167 (RSE 37%)

    # ---- Inter-individual variability (Bukkems 2021 Table 2 IIV block) ----
    # Table 2 reports IIV as %CV with the log-normal-variance-to-%CV transform
    # sqrt(exp(variance) - 1); the internal omega^2 on the log scale is recovered as
    # omega^2 = log(1 + CV^2). See Table 2 footnote (d).
    #   CV(Vc/F) = 69.7%  -> omega^2 = log(1 + 0.697^2) = 0.39596
    #   CV(CL/F) = 28.6%  -> omega^2 = log(1 + 0.286^2) = 0.07862
    #   CV(Q/F)  = 71.5%  -> omega^2 = log(1 + 0.715^2) = 0.41292
    #   CV(Vp/F) = 115.2% -> omega^2 = log(1 + 1.152^2) = 0.84462
    # The BLOCK(3) covariance structure for (CL, Q, Vp) has corr(CL, Q) = 0.18 and
    # corr(Q, Vp) = 0.59 (Table 2); the corr(CL, Vp) off-diagonal is fixed to zero in the
    # NONMEM control stream (Supporting Information S4 $OMEGA BLOCK(3) with an explicit
    # 0 in the (Vp, CL) slot). Covariances:
    #   cov(logCL, logQ)  = 0.18 * sqrt(0.07862 * 0.41292) = 0.03243
    #   cov(logQ,  logVp) = 0.59 * sqrt(0.41292 * 0.84462) = 0.34843
    #   cov(logCL, logVp) = 0                              (fixed)
    etalvc                     ~ 0.39596                                                                    # Table 2: IIV Vc/F = 69.7% CV
    etalcl + etalq + etalvp   ~ c(0.07862,
                                   0.03243, 0.41292,
                                   fixed(0), 0.34843, 0.84462)                                              # Table 2: IIV CL/F = 28.6%, IIV Q/F = 71.5%, IIV Vp/F = 115.2%; corr(CL,Q) = 0.18, corr(Q,Vp) = 0.59, corr(CL,Vp) = 0 (fixed)

    # ---- Residual error (Bukkems 2021 Table 2 residual rows) ----
    # Time-varying proportional residual with two magnitudes: one for time-after-dose
    # PTAD <= 3 h (the absorption-dominated phase) and one for PTAD > 3 h (the disposition
    # phase). The switch time 3 h is the average observed Tmax across the pooled dataset
    # (Methods paragraph on residual-error model, and Table 2 footnote e mirroring
    # Karlsson 1998 time-varying residual). The paper's additional log-normal per-subject
    # scaling of the residual magnitude (Table 2 'IIV residual error' 25.6% CV) is omitted
    # from this packaged model; see vignette 'Assumptions and deviations'.
    propSd_early <- 0.435; label("Proportional residual SD for PTAD <= 3 h after dose (fraction)")          # Table 2: Prop residual <= 3 h = 43.5% (RSE 3%)
    propSd_late  <- 0.290; label("Proportional residual SD for PTAD > 3 h after dose (fraction)")           # Table 2: Prop residual >  3 h = 29.0% (RSE 2%)
  })

  model({
    # ---- Derived multiplicative covariate factors ----
    # Each evaluates to 1 at the reference covariate values (all indicators at 0),
    # preserving the paper's structural typical values in Table 2.
    mat_fed_factor  <- 1 + e_fed_mat            * FED
    cl_atv_factor   <- 1 + e_atazanavir_cl      * CONMED_ATAZANAVIR
    f_form_factor   <- 1 + e_ral_600_fdepot     * FORM_RAL_600
    f_lowfat_factor <- 1 + e_lowfat_fdepot      * FED_LOWFAT
    f_preg_factor   <- 1 + e_preg_fdepot        * PREG
    f_efv_factor    <- 1 + e_efv_fdepot         * CONMED_EFV

    # ---- Individual PK parameters with allometric weight scaling on 70 kg ----
    # Vc, Vp scale allometrically with exponent 1 on (WT / 70); CL, Q scale with
    # exponent 0.75. The atazanavir factor multiplies CL (linear additive form). For
    # pregnant women the paper uses postpartum weight (or an imputed 0.93 * third-
    # trimester weight when postpartum is missing) as the allometric size descriptor, so
    # pregnancy does NOT change the allometric term itself; the pregnancy effect enters
    # only through the F multiplier below.
    ka  <- exp(lka)
    mat <- exp(lmat) * mat_fed_factor
    vc  <- exp(lvc  + etalvc) * (WT / 70)^e_wt_vc
    cl  <- exp(lcl  + etalcl) * (WT / 70)^e_wt_cl * cl_atv_factor
    q   <- exp(lq   + etalq)  * (WT / 70)^e_wt_cl
    vp  <- exp(lvp  + etalvp) * (WT / 70)^e_wt_vc

    # Transit-chain rate constant. NONMEM code $PK: KTR = 3 / AMRT, where AMRT is the
    # paper's MAT parameter. Three sequential compartments (depot, transit1, transit2)
    # all lose their content at rate ktr into the next; transit3 loses at rate ka into
    # central. Total mean-absorption-time contribution from the four compartments is
    # 3 / ktr + 1 / ka = MAT + 1 / ka.
    ktr <- 3 / mat

    # Two-compartment disposition micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---- ODE system: depot -> transit1 -> transit2 -> transit3 -> central
    #                  central <-> peripheral1 (two-compartment) ----
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot     - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1  - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2  - ka  * transit3
    d/dt(central)     <-  ka  * transit3  - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                     k12 * central                 - k21 * peripheral1

    # ---- Bioavailability ----
    # F = 1 * (formulation factor) * (low-fat factor) * (pregnancy factor) * (efavirenz factor)
    # matches the Bukkems 2021 Table 2 footnote (a) construction:
    #   F = 1 * (1 + factor change in F 600 mg formulation) * (1 + factor change in F low-fat meal)
    #         * (1 + factor change in F pregnancy) * (1 + factor change in F efavirenz co-admin).
    # No IIV on F is estimated in the paper; the F multipliers act on the deterministic
    # typical value.
    f(depot) <- exp(lfdepot) * f_form_factor * f_lowfat_factor * f_preg_factor * f_efv_factor

    # ---- Observation ----
    Cc <- central / vc

    # ---- Time-varying proportional residual error ----
    # Switch at 3 h after the most recent dose. tad() returns the time since the last
    # dose event (rxode2 built-in). early = 1 for PTAD <= 3 h (absorption-dominated,
    # higher variability), 0 otherwise; the residual SD is a linear combination of the
    # two magnitudes, matching the Andrews_2017_tacrolimus and Cirincione_2017_exenatide
    # per-record residual-switch precedent.
    early_flag <- tad() <= 3
    propSd     <- propSd_early * early_flag + propSd_late * (1 - early_flag)
    Cc ~ prop(propSd)
  })
}
