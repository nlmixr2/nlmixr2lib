Friberg_2012_voriconazole <- function() {
  description <- "Integrated population pharmacokinetic model for voriconazole in children, adolescents, and adults (Friberg 2012). Two-compartment with first-order oral absorption and mixed linear plus nonlinear (Michaelis-Menten with time-dependent Vmax) elimination; allometric scaling on all clearance terms (exponent 0.75) and on volumes (exponent 1.0) with 70 kg reference; population-specific Vmax,inh, Q, ka, and Alag for children, adolescents, and adults; CYP2C19 heterozygous extensive or poor metabolizer adults have fully blocked nonlinear clearance (Vmax,inh = 100%)."
  reference <- "Friberg LE, Ravva P, Karlsson MO, Liu P. Integrated population pharmacokinetic analysis of voriconazole in children, adolescents, and adults. Antimicrobial Agents and Chemotherapy. 2012;56(6):3032-3042. doi:10.1128/AAC.05761-11"
  vignette <- "Friberg_2012_voriconazole"
  paper_specific_etas <- c("etalkm_vmax1", "etalvmax1_ped", "etalgtf1_other", "etalgtf1_adult", "etalka_nonadult", "eta_re_nonadult")
  paper_specific_residual_sds <- c("expSdStdy1", "expSdStdy2", "expSdStdy34", "expSdStdy5Iv", "expSdStdy5Oral")
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on all clearance terms (CL, Q, Vmax,1) with exponent 0.75 and on the central and peripheral volumes (V2, V3) with exponent 1.0; reference weight 70 kg adult. Methods 'Structural-model description' and Table 3 footnote c.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used to derive a logit-scale shift on Vmax,inh for subjects under 12 years old (the AGE < 12 covariate of Table 3). The age 12 cutoff also matches the pediatric / adolescent population boundary used by the paper. Adult (study 5) status is derived from STUDY_VORI rather than AGE because the paper uses study membership rather than an age boundary for adult-only effects.",
      source_name        = "AGE"
    ),
    CYP2C19_IM = list(
      description        = "CYP2C19 intermediate-metabolizer indicator (paper-era 'HEM' heterozygous extensive metabolizer)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP2C19 EM or UM)",
      notes              = "1 = subject is a CYP2C19 heterozygous extensive metabolizer in the paper's nomenclature (modern intermediate metabolizer). The paper grouped HEM with PM because PM exposures in children and adolescents were not outliers relative to other subjects in studies 2-4. In adults (STUDY_VORI == 5) the combined HEM-or-PM indicator forces Vmax,inh to 100% (fully blocked nonlinear clearance). Source 'CYP2C19 status' column with categories UM / EM / HEM / PM (Table 2); HEM maps to CYP2C19_IM = 1. Four children in study 3 and two adolescents in study 4 with missing CYP2C19 data were assumed to be EM (CYP2C19_IM = CYP2C19_PM = 0) per Table 2 footnote d.",
      source_name        = "CYP2C19 status (HEM level)"
    ),
    CYP2C19_PM = list(
      description        = "CYP2C19 poor-metabolizer indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP2C19 EM, UM, or IM)",
      notes              = "1 = subject is a CYP2C19 poor metabolizer (two loss-of-function alleles). Combined with CYP2C19_IM in adults (STUDY_VORI == 5) to force Vmax,inh = 100% (fully blocked nonlinear clearance). Source 'CYP2C19 status' column with categories UM / EM / HEM / PM (Table 2); PM maps to CYP2C19_PM = 1.",
      source_name        = "CYP2C19 status (PM level)"
    ),
    STUDY_VORI = list(
      description        = "Friberg 2012 voriconazole study indicator (integer 1-5)",
      units              = "(integer 1-5)",
      type               = "categorical",
      reference_category = "5 (healthy adult study)",
      notes              = "Subject-level integer identifying which of the five pooled PK studies a subject belongs to: 1 / 2 / 3 = immunocompromised children (2 to <12 y), 4 = immunocompromised adolescents (12 to <17 y), 5 = healthy adults (22-55 y). Drives (i) the Study 1 pediatric typical-value modifier on Km and Vmax,1 (-0.382), (ii) the population-specific ka and Alag, (iii) the non-adult uplift on Q (+0.637), (iv) the F1 IIV magnitude (paper estimates separate omegas for adult vs non-adult), (v) the residual-error magnitude per study (with studies 3 and 4 sharing one magnitude), and (vi) the CL IIV scaling multiplier in non-adult studies (1 + 1.70). Each subject is assigned exactly one study; the value is fixed per subject. Population strata (child / adolescent / adult) are aligned with study membership (no crossovers between populations). Renamed from canonical STDY_VORI to STUDY_VORI on 2026-06-19 per the canonical-register standardization audit (typo correction: STDY was a missing-vowel abbreviation of STUDY).",
      source_name        = "Study"
    ),
    ORAL_VORI = list(
      description        = "Friberg 2012 voriconazole observation-during-oral-dose-phase indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (observation during IV phase or no dosing yet)",
      notes              = "1 = observation collected when the most recent administered dose was oral (powder for oral suspension or tablet); 0 = observation collected when the most recent administered dose was IV. Per-observation (record-level) indicator. Used only in adults (STUDY_VORI == 5) to switch the residual-error magnitude between IV-only (sigma_iv = 0.0912) and oral (sqrt(sigma_iv^2 + sigma_oral_extra^2) = sqrt(0.0912^2 + 0.132^2) = 0.160). Crossover studies (e.g. study 3 IV-then-oral) have ORAL_VORI = 0 on records during the IV period and ORAL_VORI = 1 on records during the oral period for the same subject; in studies that combine both arms the indicator is set per record from the most recent dose's route.",
      source_name        = "Route (IV vs PO/POS)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 173L,
    n_studies      = 5L,
    age_range      = "2-55 years",
    age_median     = "Children 5 y (range 2-11); Adolescents 13 y (range 12-16); Adults 34 y (range 22-55) (Table 2)",
    weight_range   = "10.8-97.0 kg",
    weight_median  = "Children 20.1 kg (range 10.8-54.9); Adolescents 57.1 kg (range 30.4-92.2); Adults 76.0 kg (range 49.0-97.0) (Table 2)",
    sex_female_pct = 41.6,
    race_ethnicity = c(Caucasian = 71.7, Black = 12.1, Asian = 3.5, Other = 13.3),
    cyp2c19_pct    = c(UM = 2.3, EM = 56.6, HEM_or_IM = 38.7, PM = 2.9),
    disease_state  = "Immunocompromised pediatric (studies 1-3) and adolescent (study 4) patients (mostly hospitalized transplant recipients with concomitant medications); healthy adult volunteers under tight control (study 5). Voriconazole was administered for treatment / prophylaxis of invasive fungal infections in the immunocompromised cohorts; healthy adults received voriconazole under a controlled PK protocol.",
    dose_range     = "Children 3-8 mg/kg IV q12h and 4-6 mg/kg or 200 mg oral suspension (POS) q12h; adolescents 6 mg/kg IV q12h on day 1, 4 mg/kg IV q12h thereafter, and 300 mg oral tablet q12h; adults 6 mg/kg IV q12h on day 1, 4 mg/kg IV q12h thereafter, and 200 mg oral tablet q12h. Maximum IV infusion rate ~3 mg/kg/h. Oral doses given at least 1 h before or after a meal. Table 1.",
    regions        = "Multinational pediatric and adolescent studies; the adult study was conducted in healthy volunteers (region not specified in the paper).",
    notes          = "Pooled data from 5 PK studies. Study 1: 12 children Jul-Dec 2000. Study 2: 24 children Jun 2003-Jun 2004. Study 3: 40 children Dec 2008-Oct 2009. Study 4: 26 adolescents Jun 2008-Dec 2009. Study 5: 35 healthy adults Apr 2009-Jul 2009. Total 3336 voriconazole plasma concentrations (children 2022, adolescents 554, adults 760). Two outlier records with absolute CWRES > 6 were excluded from the final dataset (Methods 'Population pharmacokinetic modeling')."
  )

  ini({
    # Structural parameters - reference is a 70 kg adult in study 5 (no study or
    # age modifiers active). Time scale hours; concentration ug/mL; dose mg.

    # Km (Michaelis-Menten constant)
    lkm                       <- log(1.15)   ; label("Michaelis-Menten constant Km (ug/mL)")                              # Table 3: thetaKm = 1.15
    e_stdy1_ped_km            <- -0.382      ; label("Additive-fractional modifier on Km and Vmax,1 for Study 1 pediatric (unitless)")  # Table 3: thetaSTDY1,ped = -0.382 (shared between Km and Vmax,1)

    # Vmax at 1 h (peak nonlinear elimination rate per 70 kg)
    lvmax1                    <- log(114)    ; label("Maximum nonlinear elimination rate Vmax at 1 h (mg/h per 70 kg)")    # Table 3: thetaVmax,1 = 114
    e_wt_vmax1                <- fixed(0.75) ; label("Allometric exponent on Vmax,1 (unitless)")                            # Methods 'Structural-model description': fixed at 0.75

    # Vmax,inh (maximum fractional inhibition of Vmax over time; logit-scale)
    lgt_vmax_inh              <- 1.50        ; label("logit(Vmax,inh): logit of the maximum fractional inhibition of Vmax (unitless)")  # Table 3: thetaVmax,inh = 1.50; expit(1.50) = 0.818
    e_age_lt12_lgt_vmax_inh   <- -0.390      ; label("Additive shift on logit(Vmax,inh) for AGE < 12 y (unitless)")          # Table 3: thetaAGE<12 = -0.390; child Vmax,inh = expit(1.50 - 0.390) = expit(1.11) = 0.752

    # T50 (time at which half the maximum inhibition occurs)
    lt50                      <- log(2.41)   ; label("T50 of the time-dependent Vmax inhibition (h)")                       # Table 3: thetaT50 = 2.41

    # CL (linear clearance per 70 kg)
    lcl                       <- log(6.16)   ; label("Linear clearance CL (L/h per 70 kg)")                                 # Table 3: thetaCL = 6.16
    e_wt_cl                   <- fixed(0.75) ; label("Allometric exponent on CL (unitless)")                                 # Methods 'Structural-model description': fixed at 0.75

    # V2 (central volume per 70 kg)
    lvc                       <- log(79.0)   ; label("Central volume of distribution V2 (L per 70 kg)")                     # Table 3: thetaV2 = 79.0
    e_wt_vc                   <- fixed(1.0)  ; label("Allometric exponent on V2 (unitless)")                                 # Methods 'Structural-model description': fixed at 1.0

    # V3 (peripheral volume per 70 kg)
    lvp                       <- log(103)    ; label("Peripheral volume of distribution V3 (L per 70 kg)")                  # Table 3: thetaV3 = 103
    e_wt_vp                   <- fixed(1.0)  ; label("Allometric exponent on V3 (unitless)")                                 # Methods 'Structural-model description': fixed at 1.0

    # Q (intercompartmental clearance per 70 kg, adult reference)
    lq                        <- log(15.5)   ; label("Intercompartmental clearance Q in adults (L/h per 70 kg)")            # Table 3: thetaQ = 15.5
    e_wt_q                    <- fixed(0.75) ; label("Allometric exponent on Q (unitless)")                                  # Methods 'Structural-model description': fixed at 0.75
    e_notstdy5_q              <- 0.637       ; label("Additive-fractional uplift on Q for non-adult studies (unitless)")     # Table 3: thetaQ,notSTDY5,adult = 0.637 -> Q in non-adult studies = Q * 1.637

    # F1 (oral bioavailability, logit-scale)
    lgt_f1                    <- 0.585       ; label("logit(F1): logit of typical oral bioavailability (unitless)")          # Table 3: thetaF1 = 0.585; expit(0.585) = 0.642 (typical F1 ~ 0.64)

    # ka (absorption rate constant)
    lka_ped                   <- log(1.19)   ; label("Absorption rate ka in pediatric studies (1/h)")                        # Table 3: thetaka = 1.19 (children studies 1-3)
    e_stdy4_adol_ka           <- -0.615      ; label("Additive-fractional shift on ka for Study 4 (adolescent) (unitless)")  # Table 3: thetaSTDY4,adol = -0.615 -> ka in adolescents = 1.19 * (1 - 0.615) = 0.458
    ka_adult                  <- fixed(100)  ; label("Absorption rate ka in adults (1/h) -- structurally fixed at a high value because Alag absorbs the lag")  # Table 3: thetaSTDY5,adult ka = 100 (fixed; oral bioavailability and Alag are estimated against this high ka)

    # Alag (absorption lag time)
    lalag_adult               <- log(0.949)  ; label("Absorption lag time Alag in adults (h)")                               # Table 3: thetaAlag = 0.949
    e_notstdy5_alag           <- -0.874      ; label("Additive-fractional shift on Alag for non-adult studies (unitless)")   # Table 3: thetaAlag,notSTDY5,adult = -0.874 -> non-adult Alag = 0.949 * (1 - 0.874) = 0.120

    # ---- Vmax,1 IIV scaling factors (paper Table 3 'Vmax,1' SD column footnote e) ----
    # The paper parameterises Vmax,1 IIV by scaling the joint Km / Vmax,1 eta differently in
    # each population: adults use the joint eta scaled by 0.584, adolescents use the joint eta
    # scaled by 0.208, and children use a separate eta with its own variance (omega = 0.239).
    # The two scalars are independent of the omega values (they multiply the realised eta in
    # the typical-value equation for Vmax,1), so they enter ini() as estimated thetas.
    e_stdy5_adult_vmax1_scale <- 0.584       ; label("Adult scaling factor on the Km/Vmax,1 joint eta for Vmax,1 (unitless)")        # Table 3: thetaVmax,scale,adult = 0.584
    e_stdy4_adol_vmax1_scale  <- 0.208       ; label("Adolescent scaling factor on the Km/Vmax,1 joint eta for Vmax,1 (unitless)")  # Table 3: thetaVmax,scale,adol = 0.208

    # ---- IIVs ----
    # (1) Joint Km / Vmax,1 block: etalkm_vmax1 is the shared eta used by Km in every population
    # and by Vmax,1 in adolescents and adults (scaled). etalvmax1_ped is the pediatric-specific
    # eta used by Vmax,1 in children. The block carries the 81% Km / Vmax,1 children correlation
    # reported in Results 'Voriconazole parameter estimates' (covariance = 0.81 * 1.36 * 0.239).
    etalkm_vmax1 + etalvmax1_ped ~ c(
      1.36 * 1.36,                # var(etalkm_vmax1)         # Table 3: omega_Km-Vmax,1 = 1.36
      0.81 * 1.36 * 0.239,        # cov(etalkm_vmax1, etalvmax1_ped) (81% Km / Vmax,1 children correlation)  # Results: Km and Vmax,1 in children 81% correlation
      0.239 * 0.239               # var(etalvmax1_ped)         # Table 3: omega_Vmax1,ped = 0.239
    )

    # (2) CL log-normal IIV; the realised eta is scaled by (1 + 1.70 * notSTDY5) in non-adult studies.
    etalcl ~ 0.435 * 0.435                                                                                                    # Table 3: omega_CL = 0.435 (adults)
    e_notstdy5_etalcl         <- 1.70        ; label("Multiplicative factor on the CL IIV magnitude in non-adult studies (unitless)")  # Table 3: theta_eta_CL,notSTDY5,adult = 1.70

    # (3) V2, V3, Q simple log-normal IIVs.
    etalvc ~ 0.136 * 0.136                                                                                                    # Table 3: omega_V2 = 0.136
    etalvp ~ 0.769 * 0.769                                                                                                    # Table 3: omega_V3 = 0.769
    etalq  ~ 0.424 * 0.424                                                                                                    # Table 3: omega_Q = 0.424

    # (4) F1 IIVs: study-specific magnitudes with Box-Cox-transformed eta on logit(F1).
    etalgtf1_other ~ 1.67 * 1.67                                                                                              # Table 3: omega_F,notSTDY5,adult = 1.67
    etalgtf1_adult ~ 0.686 * 0.686                                                                                            # Table 3: omega_F,STDY5,adult = 0.686
    bc_f1                     <- 0.367       ; label("Box-Cox transformation lambda on the F1 eta (unitless)")               # Table 3: thetaBC-F = 0.367

    # (5) ka log-normal IIV, only applied to non-adult studies.
    etalka_nonadult ~ 0.898 * 0.898                                                                                           # Table 3: omega_ka = 0.898

    # (6) Residual error: study-specific log-scale SDs.
    expSdStdy1                <- 0.593       ; label("Log-scale residual SD in Study 1 (pediatric)")                          # Table 3: thetaSTDY1,ped = 0.593
    expSdStdy2                <- 0.425       ; label("Log-scale residual SD in Study 2 (pediatric)")                          # Table 3: thetaSTDY2,ped = 0.425
    expSdStdy34               <- 0.365       ; label("Log-scale residual SD in Studies 3 and 4 (pediatric / adolescent)")    # Table 3: thetaSTDY3,4,ped,adol = 0.365
    expSdStdy5Iv              <- 0.0912      ; label("Log-scale residual SD in Study 5 (adult) IV phase")                    # Table 3: thetaSTDY5,adult = 0.0912
    expSdStdy5Oral            <- 0.132       ; label("Additional log-scale residual SD term for Study 5 (adult) oral phase") # Table 3: thetaSTDY5,adult,oral = 0.132 (combined adult oral SD = sqrt(0.0912^2 + 0.132^2) = 0.160)

    # (7) Individual residual-error eta, active only in non-adult studies.
    eta_re_nonadult ~ 0.456 * 0.456                                                                                            # Table 3: omega_RE = 0.456
  })

  model({
    # ---- 1. Derived population and study indicators ----
    is_child           <- (AGE < 12)
    is_stdy1_ped       <- (STUDY_VORI == 1)
    is_stdy4_adol      <- (STUDY_VORI == 4)
    is_stdy5_adult     <- (STUDY_VORI == 5)
    is_not_stdy5       <- 1.0 - is_stdy5_adult
    is_pediatric       <- (STUDY_VORI <= 3) * 1.0
    is_hem_or_pm_adult <- is_stdy5_adult * (CYP2C19_IM + CYP2C19_PM)

    # ---- 2. Vmax,inh (typical and HEM/PM-adult override) ----
    lgt_vmax_inh_i    <- lgt_vmax_inh + e_age_lt12_lgt_vmax_inh * is_child
    vmax_inh_typical  <- 1.0 / (1.0 + exp(-lgt_vmax_inh_i))
    vmax_inh          <- vmax_inh_typical * (1.0 - is_hem_or_pm_adult) + 1.0 * is_hem_or_pm_adult

    # ---- 3. T50 and time-dependent Vmax driver ----
    t50  <- exp(lt50)
    teff <- t - 1.0
    if (teff < 0.0) teff <- 0.0

    # ---- 4. Individual parameters with allometric and population effects ----
    # Km: uses the joint Km / Vmax,1 eta in every population; Study 1 pediatric carries the
    # -0.382 typical-value shift.
    km   <- exp(lkm + etalkm_vmax1) * (1.0 + e_stdy1_ped_km * is_stdy1_ped)

    # Vmax,1: the random-effect contribution is population-specific (joint eta scaled in
    # adolescents and adults; pediatric-specific eta in children); allometric scaling and the
    # Study 1 pediatric typical-value shift are applied multiplicatively.
    eta_vmax1 <- etalkm_vmax1   * e_stdy5_adult_vmax1_scale * is_stdy5_adult +
                 etalkm_vmax1   * e_stdy4_adol_vmax1_scale  * is_stdy4_adol  +
                 etalvmax1_ped                              * is_pediatric
    vmax1 <- exp(lvmax1 + eta_vmax1) *
             (WT / 70.0)^e_wt_vmax1 *
             (1.0 + e_stdy1_ped_km * is_stdy1_ped)

    # CL: log-normal IIV scaled in non-adult studies.
    etalcl_eff <- etalcl * (1.0 + e_notstdy5_etalcl * is_not_stdy5)
    cl <- exp(lcl + etalcl_eff) * (WT / 70.0)^e_wt_cl

    # V2 and V3.
    vc <- exp(lvc + etalvc) * (WT / 70.0)^e_wt_vc
    vp <- exp(lvp + etalvp) * (WT / 70.0)^e_wt_vp

    # Q: non-adult studies receive a +63.7% multiplicative bump on the typical Q.
    q <- exp(lq + etalq) * (WT / 70.0)^e_wt_q * (1.0 + e_notstdy5_q * is_not_stdy5)

    # F1: study-specific eta variance, with Box-Cox transformation applied on exp(eta).
    eta_f1_raw <- etalgtf1_adult * is_stdy5_adult + etalgtf1_other * is_not_stdy5
    eta_f1_bc  <- (exp(bc_f1 * eta_f1_raw) - 1.0) / bc_f1
    lgt_f1_i   <- lgt_f1 + eta_f1_bc
    f1         <- 1.0 / (1.0 + exp(-lgt_f1_i))

    # ka: population-specific typical value plus a non-adult log-normal IIV.
    ka_typical_ped  <- exp(lka_ped)
    ka_typical_adol <- ka_typical_ped * (1.0 + e_stdy4_adol_ka)
    ka_typical      <- ka_typical_ped  * is_pediatric  +
                       ka_typical_adol * is_stdy4_adol +
                       ka_adult        * is_stdy5_adult
    ka <- ka_typical * exp(etalka_nonadult * is_not_stdy5)

    # Alag.
    alag_typical <- exp(lalag_adult) * (1.0 + e_notstdy5_alag * is_not_stdy5)

    # ---- 5. Time-dependent Vmax (mg/h) ----
    vmax_t <- vmax1 * (1.0 - vmax_inh * teff / (teff + (t50 - 1.0)))

    # ---- 6. ODE system ----
    # central in mg, vc in L => Cp = central / vc in mg/L = ug/mL.
    Cp            <- central / vc
    cl_nonlin_amt <- vmax_t / (Cp + km)

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - (cl + cl_nonlin_amt) * Cp - q * Cp + q * (peripheral1 / vp)
    d/dt(peripheral1) <-                                            q * Cp - q * (peripheral1 / vp)

    f(depot)    <- f1
    alag(depot) <- alag_typical

    # ---- 7. Observation and study-specific residual error ----
    Cc <- Cp

    expSdStdy5Combined <- sqrt(expSdStdy5Iv^2 + expSdStdy5Oral^2 * ORAL_VORI)
    w_typical <-
      expSdStdy1         * (STUDY_VORI == 1) +
      expSdStdy2         * (STUDY_VORI == 2) +
      expSdStdy34        * ((STUDY_VORI == 3) + (STUDY_VORI == 4)) +
      expSdStdy5Combined * is_stdy5_adult
    w_indiv <- w_typical * exp(eta_re_nonadult * is_not_stdy5)

    Cc ~ lnorm(w_indiv)
  })
}
