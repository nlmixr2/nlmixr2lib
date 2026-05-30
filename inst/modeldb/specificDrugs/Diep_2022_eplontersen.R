Diep_2022_eplontersen <- function() {
  description <- "Two-compartment population PK and indirect-response PD model for the GalNAc3-conjugated antisense oligonucleotide eplontersen targeting transthyretin (TTR) mRNA, fit to pooled data from two phase 1 studies in healthy volunteers (Diep 2022). First-order SC absorption with site-specific typical ka (arm vs abdomen), allometric scaling on CL by lean body mass, on Vc/Q/Vp by total body weight, and an indirect-response model with eplontersen-driven inhibition of TTR production."
  reference   <- "Diep JK, Yu RZ, Viney NJ, Schneider E, Guo S, Henry S, Monia B, Geary R, Wang Y. Population pharmacokinetic/pharmacodynamic modelling of eplontersen, an antisense oligonucleotide in development for transthyretin amyloidosis. Br J Clin Pharmacol. 2022;88(12):5389-5398. doi:10.1111/bcp.15468"
  vignette    <- "Diep_2022_eplontersen"
  units       <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight (baseline)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Power exponents on Vc (1.89), Q (2.53), and Vp (2.73); reference 72.1 kg = cohort median (Diep 2022 Table 1 and Eqs 2-4). The source paper labels this column BW.",
      source_name        = "BW"
    ),
    LBM = list(
      description        = "Lean body mass (baseline)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Power exponent 1.42 on CL; reference 51.6 kg = cohort median (Diep 2022 Table 1 and Eq 1). Lean body mass (NOT lean body weight, despite both being reported in the demographics table) is the variable used in the final model.",
      source_name        = "LBM"
    ),
    INJSITE_ARM = list(
      description        = "SC injection-site indicator: 1 = arm, 0 = abdomen",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (abdomen)",
      notes              = "Per-dose-record covariate. Diep 2022 estimated separate typical first-order absorption rate constants for arm (ka_arm = 0.217 1/h) and abdomen (ka_ab = 0.282 1/h) with a single shared eta on log(ka). The model encodes the abdomen value as the typical-value reference (consistent with abdomen being the universal SC reference site across the popPK literature) and an additive log-shift covariate effect e_injsite_arm_ka = log(ka_arm / ka_ab) on the typical log(ka) when INJSITE_ARM = 1. In multi-dose simulations a subject can switch sites between doses, so this is per-administration rather than per-subject.",
      source_name        = "Injection site (paper narrative; arm vs abdomen subgroup labels driving ka_arm vs ka_ab)"
    )
  )

  population <- list(
    n_subjects       = 55L,
    n_studies        = 2L,
    age_range        = "23 - 65 years",
    age_median       = "54 years",
    weight_range     = "50.4 - 97.0 kg",
    weight_median    = "72.1 kg",
    height_range     = "146 - 189 cm",
    bmi_range        = "18.7 - 30.7 kg/m^2",
    lbm_range        = "22.8 - 66.3 kg",
    sex_female_pct   = 34.5,
    race_ethnicity   = c(Caucasian = 30.9, `Black or African American` = 18.2, Asian = 50.9),
    disease_state    = "Healthy volunteers (no transthyretin amyloidosis); the pooled analysis included an ethnobridging cohort of Japanese descent.",
    dose_range       = "Subcutaneous eplontersen as a single 120 mg dose, single-ascending 45/60/90 mg cohorts, or 45/60/90 mg every 4 weeks for 4 doses (days 1, 29, 57, 85). Administration alternated between arm and abdomen in multi-dose cohorts.",
    regions          = "Phase 1 study NCT03728634 conducted in Canada (Western volunteers, dose escalation); phase 1 study NCT04302064 conducted in healthy volunteers of Japanese descent (single ascending dose).",
    studies          = "NCT03728634 (n = 47 enrolled; 1 single-dose 120 mg cohort and 3 multi-dose 45/60/90 mg cohorts; 10:2 randomization to active:placebo) and NCT04302064 (n = 24 enrolled; 3 single-dose 45/60/90 mg cohorts in Japanese descendants; 6:2 randomization).",
    notes            = "PK/PD analysis pooled n = 55 active-arm subjects after excluding 14 placebo subjects and 2 subjects with pre-existing antidrug antibodies (Diep 2022 Section 3.1). Final dataset: 1260 plasma eplontersen concentrations and 624 serum TTR concentrations. PK quantification used hybridization-based ECL with LLOQ 0.129 ng/mL; TTR quantification used ELISA with LLOQ 0.896 mg/dL (Section 2.1). Data below LLOQ (11.7% of PK observations) excluded per M1 method."
  )

  ini({
    # ---- Structural PK parameters (Diep 2022 Table 2 final-model column) ----
    # Reference subject: median total body weight 72.1 kg, median lean body mass
    # 51.6 kg, abdomen SC injection (INJSITE_ARM = 0). lka is the typical
    # log(ka) for abdomen; e_injsite_arm_ka is an additive log-shift covariate
    # effect that recovers the arm typical value (ka_arm = 0.217 1/h) when
    # INJSITE_ARM = 1.
    lka              <- log(0.282); label("First-order SC absorption rate constant for abdomen injection (ka_ab, 1/h; INJSITE_ARM = 0 reference)") # Diep 2022 Table 2 final-model ka_ab = 0.282 1/h
    e_injsite_arm_ka <- log(0.217 / 0.282); label("Additive log-shift covariate effect of INJSITE_ARM on lka (log(ka_arm / ka_ab), unitless)") # Diep 2022 Table 2 final-model: ka_arm / ka_ab = 0.217 / 0.282
    lcl              <- log(24.1);  label("Linear plasma clearance at reference LBM (CL, L/h)")                          # Diep 2022 Table 2 final-model CL = 24.1 L/h
    lvc              <- log(50.4);  label("Central volume of distribution at reference WT (Vc, L)")                      # Diep 2022 Table 2 final-model Vc = 50.4 L
    lq               <- log(3.64);  label("Intercompartmental clearance at reference WT (Q, L/h)")                       # Diep 2022 Table 2 final-model Q = 3.64 L/h
    lvp              <- log(2790);  label("Peripheral volume of distribution at reference WT (Vp, L)")                   # Diep 2022 Table 2 final-model Vp = 2790 L

    # ---- Covariate effects on PK (Diep 2022 Eqs 1-4) ----
    e_lbm_cl <- 1.42; label("Power exponent of LBM on CL (unitless)") # Diep 2022 Eq 1: CL = tvCL * (LBM / 51.6)^1.42
    e_wt_vc  <- 1.89; label("Power exponent of WT on Vc (unitless)")  # Diep 2022 Eq 2: Vc = tvVc * (BW / 72.1)^1.89
    e_wt_q   <- 2.53; label("Power exponent of WT on Q (unitless)")   # Diep 2022 Eq 3: Q  = tvQ  * (BW / 72.1)^2.53
    e_wt_vp  <- 2.73; label("Power exponent of WT on Vp (unitless)")  # Diep 2022 Eq 4: Vp = tvVp * (BW / 72.1)^2.73

    # ---- Structural PD parameters (Diep 2022 Table 3 final-model column) ----
    # Indirect-response model with inhibition of TTR production:
    #   d/dt(effect) = kin * (1 - Cp * Imax / (IC50 + Cp)) - kout * effect
    #   kin = rbase * kout (so effect = rbase at the no-drug steady state)
    # The PD compartment is the canonical `effect`; the observed variable
    # is the same TTR concentration assigned to a paper-named output `ttr`
    # for clarity in vignettes and PKNCA recipes.
    lrbase    <- log(31.4);    label("Baseline serum transthyretin (BL, mg/dL)")                                       # Diep 2022 Table 3 BL = 31.4 mg/dL
    lkout  <- log(0.00398); label("First-order TTR loss rate constant (kout, 1/h)")                                 # Diep 2022 Table 3 kout = 0.00398 1/h
    imax   <- 0.970;        label("Maximum fractional inhibition of TTR production by eplontersen (Imax, unitless)") # Diep 2022 Table 3 Imax = 0.970
    lic50  <- log(0.0283);  label("Plasma eplontersen concentration yielding half-maximum inhibition (IC50, ng/mL)") # Diep 2022 Table 3 IC50 = 0.0283 ng/mL

    # ---- IIV (Diep 2022 Tables 2 and 3; log-normal eta with omega^2 = log(CV^2 + 1)) ----
    # Single eta on log(ka) shared between arm and abdomen typical values
    # (paper: "interoccasion variability... included on ka to account for
    # injection site differences"; reported shrinkages 22.2% arm and 21.6%
    # abdomen are subgroup-specific, but the underlying eta is one parameter
    # paired with the reference-site lka).
    etalcl   ~ 0.04357   # Diep 2022 Table 2 IIV%CL = 21.1%   -> log(0.211^2 + 1)
    etalvc   ~ 0.24017   # Diep 2022 Table 2 IIV%Vc = 52.1%   -> log(0.521^2 + 1)
    etalq    ~ 0.12777   # Diep 2022 Table 2 IIV%Q  = 36.9%   -> log(0.369^2 + 1)
    etalvp   ~ 0.18147   # Diep 2022 Table 2 IIV%Vp = 44.6%   -> log(0.446^2 + 1)
    etalka   ~ 0.14378   # Diep 2022 Table 2 IIV%ka = 39.3%   -> log(0.393^2 + 1)
    etalrbase   ~ 0.09633   # Diep 2022 Table 3 IIV%BL  = 31.8%  -> log(0.318^2 + 1)
    etalkout ~ 0.19349   # Diep 2022 Table 3 IIV%kout = 46.2% -> log(0.462^2 + 1)
    etalic50 ~ 0.52295   # Diep 2022 Table 3 IIV%IC50 = 82.9% -> log(0.829^2 + 1)

    # ---- Residual error ----
    # PK: Diep 2022 Table 2 reports sigma^2_add = 0.0851 as the additive
    # component on log-transformed plasma concentrations. NONMEM/Phoenix
    # "additive on log scale" maps to a proportional residual in linear
    # space; propSd = sqrt(0.0851) = 0.2917 (the per-conventions translation
    # in references/naming-conventions.md NONMEM error-block table).
    # PD: Diep 2022 Table 3 reports sigma^2_prop = 0.0368 as the proportional
    # component on linear TTR concentrations; propSd_ttr = sqrt(0.0368) = 0.1918.
    propSd     <- 0.2917; label("Proportional residual error on plasma eplontersen Cc (fraction)")        # Diep 2022 Table 2 sigma^2_add = 0.0851 (log-additive form) -> sqrt(0.0851)
    propSd_ttr <- 0.1918; label("Proportional residual error on serum TTR concentration ttr (fraction)") # Diep 2022 Table 3 sigma^2_prop = 0.0368 -> sqrt(0.0368)
  })

  model({
    # ---- 1. Individual PK parameters ----
    # Typical lka switches between abdomen (INJSITE_ARM = 0; lka) and arm
    # (INJSITE_ARM = 1; lka + e_injsite_arm_ka) by an additive log shift; the
    # same etalka applies to both sites so a subject's relative ka deviation
    # is preserved across switches in injection site.
    ka <- exp(lka + e_injsite_arm_ka * INJSITE_ARM + etalka)

    cl <- exp(lcl + etalcl) * (LBM / 51.6)^e_lbm_cl
    vc <- exp(lvc + etalvc) * (WT  / 72.1)^e_wt_vc
    q  <- exp(lq  + etalq ) * (WT  / 72.1)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT  / 72.1)^e_wt_vp

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---- 2. Individual PD parameters ----
    # kin is derived so that effect = rbase is the no-drug steady state
    # (Diep 2022 Figure 1: BL = kin / kout).
    rbase   <- exp(lrbase   + etalrbase)
    kout <- exp(lkout + etalkout)
    ic50 <- exp(lic50 + etalic50)
    kin  <- rbase * kout

    # ---- 3. ODE system ----
    # depot, central, peripheral1: amounts in mg (subcutaneous dose enters
    # depot in mg). Vc and Vp in L; central / vc has units mg/L = ug/mL.
    # Multiply by 1000 to express plasma concentration in ng/mL so that the
    # paper's IC50 (0.0283 ng/mL) shares units with Cc.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    Cc <- central / vc * 1000

    # ---- 4. PD: indirect response with inhibition of TTR production ----
    # Inhibition factor 1 - Cp * Imax / (IC50 + Cp) with Cp in ng/mL
    # (Diep 2022 Figure 1 inset equation). At Cp = 0 the factor is 1
    # (uninhibited basal production); at Cp >> IC50 the factor approaches
    # 1 - Imax = 0.030 (97% inhibition).
    inh          <- 1 - Cc * imax / (ic50 + Cc)
    d/dt(effect) <- kin * inh - kout * effect
    effect(0)    <- rbase

    # ttr is the paper-named alias of the canonical `effect` PD compartment
    # (serum TTR concentration in mg/dL).
    ttr <- effect

    # ---- 5. Observations ----
    Cc  ~ prop(propSd)
    ttr ~ prop(propSd_ttr)
  })
}
