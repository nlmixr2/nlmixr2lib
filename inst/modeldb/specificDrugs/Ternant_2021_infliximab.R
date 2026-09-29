Ternant_2021_infliximab <- function() {
  description <- "Two-compartment population PK model of intravenous infliximab with double quasi-steady-state (QSS) target-mediated drug disposition (TMDD): infliximab binds TNF-alpha in both the central and the peripheral compartment, each with its own baseline target level, steady-state dissociation constant and complex elimination rate, and a shared fixed TNF-alpha elimination rate. Fitted jointly to adults with inflammatory bowel disease (Crohn's disease or ulcerative colitis; TMDD active) and adults with ankylosing spondylitis (linear reference with no target interaction). Covariates: body weight and sex on V1, sex on CL, ulcerative colitis on central baseline TNF-alpha (Ternant 2021)."
  reference <- "Ternant D, Le Tilly O, Picon L, Moussata D, Passot C, Bejan-Angoulvant T, Desvignes C, Mulleman D, Goupille P, Paintaud G. Infliximab Efficacy May Be Linked to Full TNF-alpha Blockade in Peripheral Compartment - A Double Central-Peripheral Target-Mediated Drug Disposition (TMDD) Model. Pharmaceutics. 2021;13(11):1821. doi:10.3390/pharmaceutics13111821"
  vignette <- "Ternant_2021_infliximab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  # The peripheral total-target state extends the canonical `total_target`
  # (central) with a second, peripheral pool. `targetLocationRegex` covers only
  # `target_*` / `complex_*`, so the name is declared here, following
  # Le_2015_lampalizumab_cyno.R.
  paper_specific_compartments <- c("total_target_peripheral1")

  compartmentData <- list(
    central = list(
      analyte = "infliximab (total: unbound + TNF-alpha-bound)",
      units = "mg",
      specimen = "serum",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "infliximab (total: unbound + TNF-alpha-bound)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    total_target = list(
      analyte = "TNF-alpha (total: unbound + infliximab-bound), central compartment",
      units = "nM",
      specimen = "not applicable",
      verified = TRUE
    ),
    total_target_peripheral1 = list(
      analyte = "TNF-alpha (total: unbound + infliximab-bound), peripheral compartment",
      units = "nM",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on V1 centred on the median body weight: V1 = V1_TV * (WT/65)^0.33. Ternant 2021 Methods 2.2.2 states the weight effect is 'centered on its median' but does not print the pooled median; 65 kg is the pooled median of the 158 analysed patients reconstructed from the per-cohort medians and IQRs in Table 1 (IBD 64 kg [56-72], n = 133; AS 75 kg [65-85], n = 25). See the vignette Assumptions section.",
      source_name = "BW"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (female)",
      notes = "Ternant 2021 codes SX with females as the reference category and estimates a log-additive shift for males: ln(theta_TV) = ln(theta_female) + beta_SX * SX_male (Methods 2.2.2). The model derives the male indicator as 1 - SEXF and applies exp(e_sexmale_vc * (1 - SEXF)) on V1 and exp(e_sexmale_cl * (1 - SEXF)) on CL.",
      source_name = "SX"
    ),
    DIS_CD = list(
      description = "Crohn's disease patient indicator, 1 = Crohn's disease, 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-CD subject; in Ternant 2021 the complement is ankylosing spondylitis or ulcerative colitis)",
      notes = "Together with DIS_UC this rebuilds the paper's IBD switch (IBD = DIS_CD + DIS_UC; Appendix A 'IBD = {use=regressor} ;1=IBD, 0=AS'). When IBD = 1 the central and peripheral TMDD terms are active and the target states start at their baselines; when both indicators are 0 the subject is an ankylosing-spondylitis reference patient with linear two-compartment kinetics and no target. Crohn's disease is the reference category of the ulcerative-colitis effect on central baseline TNF-alpha, so DIS_CD itself carries no coefficient.",
      source_name = "IBD / DIS"
    ),
    DIS_UC = list(
      description = "Ulcerative colitis patient indicator, 1 = ulcerative colitis, 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-UC subject; in Ternant 2021 the complement is Crohn's disease, the reference category, or ankylosing spondylitis)",
      notes = "Log-additive shift on the central baseline TNF-alpha level: R0_C = 3.3 * exp(0.57 * DIS_UC) nM (Table 2 'UC_RC0'; Results 3.2 gives R0_C = 5.8 nM in UC). Also contributes to the IBD switch (IBD = DIS_CD + DIS_UC); see DIS_CD.",
      source_name = "DIS"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 158L,
    n_studies = 2L,
    age_range = "IBD cohort median 34 years (IQR 25-41); AS cohort median 43 years (IQR 35-52)",
    weight_range = "41-110 kg (population range stated in Methods 2.3)",
    weight_median = "IBD cohort 64 kg (IQR 56-72); AS cohort 75 kg (IQR 65-85)",
    sex_female_pct = 37.3,
    race_ethnicity = "Not reported (French single-centre and bicentric cohorts)",
    disease_state = "133 adults with inflammatory bowel disease (108 Crohn's disease, 25 ulcerative colitis) treated in routine care, plus 25 adults with ankylosing spondylitis from the SPAXIM trial used as the linear-kinetics reference.",
    dose_range = "Infliximab 5 mg/kg IV infusion; IBD at weeks 0, 2, 6 and every 8 weeks thereafter, AS at weeks 0, 2, 6, 12 and 18. Starting dose median 300 mg (IBD) and 400 mg (AS).",
    regions = "France (Tours University Hospital; SPAXIM trial, NCT00607403)",
    n_observations = "1333 infliximab serum concentrations (845 IBD, 488 AS)",
    notes = "Table 1 of Ternant 2021. IBD patients with anti-drug antibodies in the first three cycles were excluded; later cycles with antibodies were discarded. Concentrations measured by a validated ELISA of unbound infliximab (LLOQ 0.103 mg/L). Estimation in MONOLIX Suite 2020 (SAEM)."
  )

  ini({
    # ---- Linear two-compartment PK (Table 2, 'Final TMDD Central + Peripheral') ----
    # Typical values are for a female patient at the 65 kg centring weight.
    lcl <- log(0.16); label("Clearance CL, female reference (L/day)") # Table 2 final model CL = 0.16 L.day-1 (RSE 9.3%)
    lvc <- log(2.6); label("Central volume V1, female reference at 65 kg (L)") # Table 2 final model V1 = 2.6 L (RSE 3.3%)
    lvp <- log(1.9); label("Peripheral volume V2 (L)") # Table 2 final model V2 = 1.9 L (RSE 8.8%)
    lq <- log(1.8); label("Intercompartmental clearance Q (L/day)") # Table 2 final model Q = 1.8 L.day-1 (RSE 2.0%)

    # ---- Central-compartment TMDD (QSS) ----
    lkss <- log(15.4); label("Central steady-state dissociation constant KSS_C (nM)") # Table 2 final model K C SS = 15.4 nM (RSE 21%)
    lrbase_target <- log(3.3); label("Central baseline TNF-alpha R0_C, Crohn's disease reference (nM)") # Table 2 final model R C 0 = 3.3 nM (RSE 28%)
    lkint <- log(0.17); label("Central infliximab-TNF-alpha complex elimination rate kint_C (1/day)") # Table 2 final model k C int = 0.17 day-1 (RSE 11%)

    # ---- Peripheral-compartment TMDD (QSS) ----
    lkss_peripheral1 <- log(0.49); label("Peripheral steady-state dissociation constant KSS_P (nM)") # Table 2 final model K P SS = 0.49 nM (RSE 11%)
    lrbase_target_peripheral1 <- log(0.46); label("Peripheral baseline TNF-alpha R0_P (nM)") # Table 2 final model R P 0 = 0.46 nM (RSE 22%)
    lkint_peripheral1 <- log(0.0079); label("Peripheral infliximab-TNF-alpha complex elimination rate kint_P (1/day)") # Table 2 final model k P int = 0.0079 day-1 (RSE 36%)

    # ---- TNF-alpha elimination, one shared value for both compartments ----
    lkdeg <- fixed(log(20)); label("TNF-alpha first-order elimination rate kout, shared by both compartments (1/day)") # Table 2 k out = 20 day-1 (fixed); Results 3.1 and Supplement Table S1 (value selected by AIC scan over 5-200 day-1)

    # ---- Covariate effects (Table 2 final model; Methods 2.2.2 forms) ----
    e_wt_vc <- 0.33; label("Power exponent of (WT/65) on V1 (unitless)") # Table 2 'BW_V 1' = 0.33 (RSE 35%)
    e_sexmale_vc <- 0.13; label("Log-additive shift on V1 for males vs female reference (unitless)") # Table 2 'SX_V 1' = 0.13 (RSE 40%)
    e_sexmale_cl <- 0.36; label("Log-additive shift on CL for males vs female reference (unitless)") # Table 2 'SX_CL' = 0.36 (RSE 26%); Results 3.2 CL_males = 0.23 L/day
    e_dis_uc_rbase_target <- 0.57; label("Log-additive shift on central baseline TNF-alpha R0_C for ulcerative colitis vs Crohn's disease (unitless)") # Table 2 'UC_R C 0' = 0.57 (RSE 47%); Results 3.2 R0_C = 5.8 nM in UC

    # ---- Inter-individual variability ----
    # Monolix reports omega as the SD of the normally distributed eta of an
    # exponential (log-normal) model; the variances below are omega^2. IIV on
    # Q, the two KSS, the two kint was set to 0 by the authors (Results 3.1).
    etalvc ~ 0.0729 # Table 2 final model omega V1 = 0.27 (SD), variance 0.27^2
    etalcl ~ 0.1225 # Table 2 final model omega CL = 0.35 (SD), variance 0.35^2
    etalvp ~ 0.1521 # Table 2 final model omega V2 = 0.39 (SD), variance 0.39^2
    etalrbase_target ~ 1.0 # Table 2 final model omega RC0 = 1.0 (SD), variance 1.0^2
    etalrbase_target_peripheral1 ~ 1.21 # Table 2 final model omega RP0 = 1.1 (SD), variance 1.1^2

    # ---- Residual error: mixed additive-proportional (Methods 2.2.2) ----
    addSd <- 1.8; label("Additive residual error (mg/L)") # Table 2 final model sigma add = 1.8 mg/L (RSE 9.8%)
    propSd <- 0.20; label("Proportional residual error (fraction)") # Table 2 final model sigma prop = 0.20 (RSE 3.0%)
  })

  model({
    # Infliximab molecular mass used by the source code to convert between mass
    # and molar units: 1 nM = 0.1442 mg/L (Appendix A:
    # 'iv(p=(1/0.144200)/V1) ;conversion concentrations mg -> nM' and
    # 'Cc = Cf*0.144200').
    mgl_per_nm <- 0.1442

    # IBD switch (Appendix A 'IBD' / 'DIS' regressors, 1 = IBD, 0 = AS).
    ibd <- DIS_CD + DIS_UC

    # ---- Individual parameters ----
    sexmale <- 1 - SEXF
    cl <- exp(lcl + etalcl + e_sexmale_cl * sexmale)
    vc <- exp(lvc + etalvc + e_sexmale_vc * sexmale) * (WT / 65)^e_wt_vc
    vp <- exp(lvp + etalvp)
    q <- exp(lq)

    kss <- exp(lkss)
    rbase_target <- exp(lrbase_target + etalrbase_target + e_dis_uc_rbase_target * DIS_UC)
    kint <- exp(lkint)
    kss_peripheral1 <- exp(lkss_peripheral1)
    rbase_target_peripheral1 <- exp(lrbase_target_peripheral1 + etalrbase_target_peripheral1)
    kint_peripheral1 <- exp(lkint_peripheral1)
    kdeg <- exp(lkdeg)

    # Zero-order TNF-alpha inputs from the baselines (Appendix A
    # 'ckin = cR0*ckout', 'pkin = pR0*pkout').
    kin <- rbase_target * kdeg
    kin_peripheral1 <- rbase_target_peripheral1 * kdeg

    # Micro-constants (Appendix A).
    k10 <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- Total infliximab concentrations in nM ----
    # The states hold total infliximab amounts in mg. The source code carries
    # both drug states as molar concentrations referenced to V1 (its peripheral
    # equation 'ddt_CP = k12*Cf - k21*Cp' makes CP = peripheral amount / V1), so
    # both amounts are divided by V1 here.
    ctot <- central / vc / mgl_per_nm
    ctot_peripheral1 <- peripheral1 / vc / mgl_per_nm

    # ---- QSS unbound concentrations (Methods 2.2.1; Appendix A) ----
    cinter <- ctot - total_target - kss
    cdelta <- max(cinter^2 + 4 * kss * ctot, 0)
    pinter <- ctot_peripheral1 - total_target_peripheral1 - kss_peripheral1
    pdelta <- max(pinter^2 + 4 * kss_peripheral1 * ctot_peripheral1, 0)
    if (ibd >= 1) {
      cfree <- 0.5 * (cinter + sqrt(cdelta))
      cfree_peripheral1 <- 0.5 * (pinter + sqrt(pdelta))
    } else {
      cfree <- ctot
      cfree_peripheral1 <- ctot_peripheral1
    }

    # ---- ODEs (Appendix A; the drug equations are scaled from nM/day to mg/day) ----
    dctot <- -k10 * cfree - k12 * cfree + k21 * cfree_peripheral1 - kint * (ctot - cfree) * ibd
    dctot_peripheral1 <- k12 * cfree - k21 * cfree_peripheral1 - kint_peripheral1 * (ctot_peripheral1 - cfree_peripheral1) * ibd
    d/dt(central) <- dctot * vc * mgl_per_nm
    d/dt(peripheral1) <- dctot_peripheral1 * vc * mgl_per_nm
    d/dt(total_target) <- (kin - kdeg * (total_target - ctot + cfree) - kint * (ctot - cfree)) * ibd
    d/dt(total_target_peripheral1) <- (kin_peripheral1 - kdeg * (total_target_peripheral1 - ctot_peripheral1 + cfree_peripheral1) - kint_peripheral1 * (ctot_peripheral1 - cfree_peripheral1)) * ibd

    # Target baselines; zero for ankylosing spondylitis (Appendix A
    # 'cRT_0 = cR0*IBD', 'pRT_0 = pR0*IBD').
    total_target(0) <- rbase_target * ibd
    total_target_peripheral1(0) <- rbase_target_peripheral1 * ibd

    # ---- Derived target readouts (Methods 2.3) ----
    target_free <- total_target - (ctot - cfree)
    target_free_peripheral1 <- total_target_peripheral1 - (ctot_peripheral1 - cfree_peripheral1)

    # ---- Observation: unbound serum infliximab (mg/L) ----
    # The ELISA measures unbound infliximab (Methods 2.1); Appendix A
    # 'Cc = Cf*0.144200'.
    Cc <- cfree * mgl_per_nm
    Cc ~ add(addSd) + prop(propSd)
  })
}
