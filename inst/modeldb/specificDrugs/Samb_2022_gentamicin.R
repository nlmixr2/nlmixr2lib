Samb_2022_gentamicin <- function() {
  description <- paste(
    "Integrated plasma-saliva population PK model for intravenous gentamicin in preterm and term",
    "neonates (Samb 2022). Plasma is the Fuchs 2014 two-compartment neonatal model with every",
    "parameter held fixed (allometric weight scaling on CL/Q and Vc/Vp, linear centred effects of",
    "gestational age on CL and Vc, postnatal age on CL and concomitant dopamine on CL). A saliva",
    "compartment is appended as a driven, non-depleting hypothetical effect compartment: first-order",
    "transfer from central (kin_saliva = 0.023 1/h) and first-order loss from saliva (kel_saliva =",
    "0.169 1/h), with the salivary concentration read on the central volume. Both saliva rate",
    "constants fall steeply with postmenstrual age (power exponents -8.8 and -5.1 referenced to",
    "244.2 days), so gentamicin appears far more readily in the saliva of premature neonates. IIV on",
    "kel_saliva only (38% CV); log-scale (exponential) residual error on saliva (49.7%).",
    sep = " "
  )
  reference <- paste(
    "Samb A, Kruizinga M, Tallahi Y, van Esdonk M, van Heel W, Driessen G, Bijleveld Y,",
    "Stuurman R, Cohen A, van Kaam A, de Haan TR, Mathot R. Saliva as a sampling matrix for",
    "therapeutic drug monitoring of gentamicin in neonates: A prospective population",
    "pharmacokinetic and simulation study. Br J Clin Pharmacol. 2022;88(4):1845-1855.",
    "doi:10.1111/bcp.15105. Plasma layer: Fuchs A, Guidi M, Giannoni E, et al. Population",
    "pharmacokinetic study of gentamicin in a large cohort of premature and term neonates.",
    "Br J Clin Pharmacol. 2014;78(5):1090-1101. doi:10.1111/bcp.12444.",
    sep = " "
  )
  vignette <- "Samb_2022_gentamicin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(
      analyte = "gentamicin",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "gentamicin",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    # Hypothetical effect compartment (Samb 2022 Methods 2.5: transport to
    # saliva is 'similar to a hypothetical effect compartment model' and its
    # mass loss from central is 'assumed to be negligible'; Figure 1 draws
    # the central-to-saliva arrow dashed). The state holds a notional amount
    # that is read on the central volume -- see the model() note.
    saliva = list(
      analyte = "gentamicin",
      units = "mg",
      specimen = "saliva",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Current body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Plasma layer only (Fuchs 2014, Samb 2022 Table S1): allometric scaling",
        "(WT/2.170 kg)^0.75 on CL and Q and (WT/2.170 kg)^1 on Vc and Vp. Table S1",
        "writes the reference as 2170 with weight in grams; the ratio is identical",
        "with WT in kg and a 2.170 kg reference. Samb 2022 cohort median 2.4 kg",
        "(range 0.7-4.3; Table 1)."
      ),
      source_name = "WT"
    ),
    GA = list(
      description = "Gestational age at birth.",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Plasma layer only: linear effects centred on 34 weeks on CL",
        "(1 + 1.87 * (GA - 34)/34) and Vc (1 - 0.922 * (GA - 34)/34), Samb 2022",
        "Table S1. Samb 2022 cohort median 34.8 weeks (range 24.3-41.7; Table 1)."
      ),
      source_name = "GA"
    ),
    PNA = list(
      description = "Postnatal age.",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Plasma layer only: linear effect on CL, 1 + 0.054 * (PNA_days - 1)/1",
        "(Samb 2022 Table S1, PNA in days, reference 1 day). The canonical PNA is",
        "in months, so model() converts with PNA_days = PNA * 30.4375. Samb 2022",
        "cohort median 1.5 days (range 0.3-6.8; Table 1)."
      ),
      source_name = "PNA"
    ),
    CONMED_DOPA = list(
      description = "Concomitant dopamine administration indicator.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant dopamine).",
      notes = paste(
        "Plasma layer only: CL multiplied by (1 - 0.120 * CONMED_DOPA) (Samb 2022",
        "Table S1, carried from Fuchs 2014). Samb 2022 does not report the",
        "number of dopamine-treated neonates; concomitant drugs were screened",
        "on the saliva parameters without effect."
      ),
      source_name = "DOPA"
    ),
    PAGE = list(
      description = "Postmenstrual age (gestational age plus postnatal age).",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Saliva layer only: power effects (PMA_days/244.2)^-8.8 on kin_saliva and",
        "(PMA_days/244.2)^-5.1 on kel_saliva (Samb 2022 Table 2 and its printed",
        "equations), where 244.2 days is the cohort median PMA (Table 1, range",
        "170.5-294.2 days). Samb 2022 reports PMA in DAYS (GA * 7 + PNA in days);",
        "the canonical PAGE is in months, so model() converts with",
        "PMA_days = PAGE * 30.4375. Supply PAGE consistent with GA and PNA."
      ),
      source_name = "PMA"
    )
  )

  covariatesDataExcluded <- list(
    BW = list(
      description = "Birth weight. Screened on the saliva parameters; not retained.",
      units = "kg",
      type = "continuous",
      notes = "Samb 2022 Methods 2.5 and Results 3.4: 'None of the other tested covariates improved the model'."
    ),
    SEXF = list(
      description = "Sex indicator (1 = female). Screened on the saliva parameters; not retained.",
      units = "(binary)",
      type = "binary",
      notes = "57.4% male in the Samb 2022 cohort (Table 1). Not retained (Results 3.4)."
    ),
    ASPHYXIA = list(
      description = "Perinatal asphyxia / controlled hypothermia indicator. Listed as candidate covariates; not tested.",
      units = "(binary)",
      type = "binary",
      notes = "Results 3.4: 'controlled hypothermia/perinatal asphyxia was not tested due to a lack of power (n = 3)'."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 54L,
    n_studies = 1L,
    n_observations = "194 analysed saliva concentrations (27 below the 0.056 mg/L LLOQ, handled with M3) and 97-99 plasma TDM concentrations",
    age_range = "GA 24.3-41.7 weeks (median 34.8); PNA 0.3-6.8 days (median 1.5); PMA 170.5-294.2 days (median 244.2)",
    weight_range = "0.7-4.3 kg current weight (median 2.4); birth weight 0.7-4.5 kg (median 2.4)",
    sex_female_pct = 42.6,
    race_ethnicity = "Not reported (two Dutch hospitals).",
    disease_state = paste(
      "Preterm and term neonates treated with intravenous gentamicin per local clinical",
      "guidelines (suspected early- or late-onset sepsis), on a NICU or a paediatric ward.",
      "21 with GA < 32 weeks, 13 with GA 32-37 weeks, 20 with GA >= 37 weeks; 3 with",
      "perinatal asphyxia / controlled hypothermia."
    ),
    dose_range = paste(
      "0.5 h IV infusion: 5 mg/kg every 48 h (GA < 32 weeks), 5 mg/kg every 36 h (GA 32-37",
      "weeks), 4 mg/kg every 24 h (GA >= 37 weeks, Emma Children's Hospital) or 5 mg/kg every",
      "36 h (GA >= 37 weeks, Juliana Children's Hospital)."
    ),
    regions = "The Netherlands (Emma Children's Hospital, Amsterdam UMC; Juliana Children's Hospital, The Hague), October 2018 - March 2020.",
    notes = paste(
      "Prospective observational study (NL7211). Saliva collected with SalivaBio Infant's",
      "Swabs, up to eight samples per neonate up to 48 h after the last dose; gentamicin",
      "C1 + C1a + C2 by LC-MS/MS. Plasma from routine TDM (1 h after the first dose and",
      "12-48 h after the first dose). Demographics from Table 1; saliva estimates from",
      "Table 2 (final model); fixed plasma model from Table S1."
    )
  )

  ini({
    # ==================================================================
    # PLASMA LAYER -- Fuchs 2014, held FIXED by Samb 2022 (Methods 2.5:
    # 'plasma PK data was described using a previously published model by
    # Fuchs et al., fixing the PK parameters'). Values from Samb 2022
    # Supporting Information Table S1, which reproduces Fuchs 2014 Table 2.
    # Reference subject: 2.170 kg, GA 34 weeks, PNA 1 day, no dopamine.
    # ==================================================================
    lcl <- fixed(log(0.089)); label("Clearance at the reference neonate (L/h)") # Table S1 theta_CL = 0.089
    lvc <- fixed(log(0.908)); label("Central volume of distribution at the reference neonate (L)") # Table S1 theta_Vc = 0.908
    lq <- fixed(log(0.157)); label("Intercompartmental clearance at the reference neonate (L/h)") # Table S1 theta_Q = 0.157
    lvp <- fixed(log(0.56)); label("Peripheral volume of distribution at the reference neonate (L)") # Table S1 theta_Vp = 0.56

    e_wt_cl_q <- fixed(0.75); label("Allometric weight exponent on CL and Q (unitless)") # Table S1 theta_CLWT = theta_QWT = 0.75
    e_wt_vc_vp <- fixed(1); label("Allometric weight exponent on Vc and Vp (unitless)") # Table S1 theta_VcWT = theta_VpWT = 1

    e_ga_cl <- fixed(1.87); label("Linear GA effect on CL centred on 34 weeks (unitless)") # Table S1 theta_CLGA = 1.87
    e_pna_cl <- fixed(0.054); label("Linear PNA effect on CL centred on 1 day (unitless)") # Table S1 theta_CLPNA = 0.054
    e_conmed_dopa_cl <- fixed(-0.12); label("Fractional change in CL with concomitant dopamine (unitless)") # Table S1 theta_CLDOPA = -0.120
    e_ga_vc <- fixed(-0.922); label("Linear GA effect on Vc centred on 34 weeks (unitless)") # Table S1 theta_VcGA = -0.922

    # IIV CL 28% and Vc 18% (CV%), correlation 87% (Table S1), converted as
    # omega^2 = log(1 + CV^2) exactly as in the shipped Fuchs_2014_gentamicin:
    #   CL: log(1 + 0.28^2) = 0.0754785; Vc: log(1 + 0.18^2) = 0.0318862;
    #   cov = 0.87 * sqrt(0.0754785 * 0.0318862) = 0.0426808
    etalcl + etalvc ~ fixed(c(
      0.0754785,
      0.0426808, 0.0318862
    )) # Table S1 IIV CL 28 pct, IIV Vc 18 pct, correlation CL-Vc 87 pct

    addSd <- fixed(0.1); label("Additive residual SD for plasma Cc (mg/L)") # Table S1 additive residual error = 0.1 mg/L
    propSd <- fixed(0.18); label("Proportional residual SD for plasma Cc (fraction)") # Table S1 proportional residual error = 18%

    # ==================================================================
    # SALIVA LAYER -- Samb 2022 Table 2, final model (OFV 738.7).
    # Paper notation k13 (central -> saliva) and k30 (loss from saliva) map
    # to the canonical kin_saliva and kel_saliva.
    # ==================================================================
    lkin_saliva <- log(0.023); label("Central-to-saliva transfer rate constant at PMA 244.2 days, k13 (1/h)") # Table 2 final model theta_k13 = 0.023 1/h (RSE 16%; bootstrap median 0.023, 95% CI 0.016-0.033)
    lkel_saliva <- log(0.169); label("Elimination rate constant from saliva at PMA 244.2 days, k30 (1/h)") # Table 2 final model theta_k30 = 0.169 1/h (RSE 15%; bootstrap median 0.171, 95% CI 0.123-0.239)

    e_page_kin_saliva <- -8.8; label("Power exponent of (PMA/244.2 days) on kin_saliva (unitless)") # Table 2 final model theta_PMA_K13 = -8.8 (RSE 16%; bootstrap median -8.7, 95% CI -11.7 to -5.7)
    e_page_kel_saliva <- -5.1; label("Power exponent of (PMA/244.2 days) on kel_saliva (unitless)") # Table 2 final model theta_PMA_K30 = -5.1 (RSE 28%; bootstrap median -4.9, 95% CI -8.1 to -2.0)

    # IIV k30 38% read as a CV: omega^2 = log(1 + 0.38^2) = 0.134886
    etalkel_saliva ~ 0.134886 # Table 2 final model IIV k30 = 38.0 pct (RSE 17%; bootstrap median 37.3, 95% CI 30.5-43.8)

    # 'Logarithmic proportional error model' on log-transformed saliva data
    # (Results 3.4): an additive error on log(Csaliva), i.e. lnorm().
    expSd_Csaliva <- 0.497; label("Log-scale residual SD for saliva Csaliva (unitless)") # Table 2 final model sigma_prop = 49.7% (RSE 7%; bootstrap median 49.0, 95% CI 40.8-56.4)
  })

  model({
    # 1. Covariate transforms. Table S1 writes the plasma equations with
    #    weight in grams (reference 2170), GA in weeks (reference 34) and
    #    PNA in days (reference 1); Table 2 writes the saliva equations with
    #    PMA in days (reference 244.2, the cohort median of Table 1).
    pna_days <- PNA * 30.4375
    pma_days <- PAGE * 30.4375

    # 2. Plasma layer (Table S1 TVCL / TVVC / Q / Vp equations).
    cl <- exp(lcl + etalcl) * (WT / 2.170)^e_wt_cl_q *
      (1 + e_ga_cl * (GA - 34) / 34) *
      (1 + e_pna_cl * (pna_days - 1) / 1) *
      (1 + e_conmed_dopa_cl * CONMED_DOPA)
    vc <- exp(lvc + etalvc) * (WT / 2.170)^e_wt_vc_vp *
      (1 + e_ga_vc * (GA - 34) / 34)
    q <- exp(lq) * (WT / 2.170)^e_wt_cl_q
    vp <- exp(lvp) * (WT / 2.170)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. Saliva layer (Table 2 printed equations k13 = theta_k13 *
    #    (PMA/244.2)^theta_PMA_k13 and k30 = theta_k30 * (PMA/244.2)^theta_PMA_k30).
    #    IIV on k30 only: a model with IIV on both k13 and k30 was rejected
    #    for eta-shrinkage of 56% and 34% (Results 3.4).
    kin_saliva <- exp(lkin_saliva) * (pma_days / 244.2)^e_page_kin_saliva
    kel_saliva <- exp(lkel_saliva + etalkel_saliva) * (pma_days / 244.2)^e_page_kel_saliva

    # 4. ODEs. Gentamicin is given as a 0.5 h IV infusion into central
    #    (supply rate or dur in the event table).
    #    SALIVA IS DRIVEN, NOT MASS-BALANCE-COUPLED: Methods 2.5 states that
    #    the 'central gentamicin mass decrease due to transport from the
    #    central compartment to the saliva compartment was assumed to be
    #    negligible', so the central equation carries no -kin_saliva*central
    #    term. There is no saliva-to-central return (Methods 2.5: oral
    #    bioavailability of gentamicin is negligible).
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(saliva) <- kin_saliva * central - kel_saliva * saliva

    # 5. Observations. No saliva volume is reported anywhere in Samb 2022
    #    (Table 2 lists only k13, k30, the two PMA exponents, sigma and IIV),
    #    and Methods 2.5 describes the saliva compartment as 'similar to a
    #    hypothetical effect compartment model', so the saliva amount is read
    #    on the central volume, as in Nguyen_2026_linezolid. This makes the
    #    late-phase saliva:plasma ratio the parameter-free
    #    kin_saliva / (kel_saliva - lambda_z), which reproduces the Figure 3
    #    individual profiles; reading saliva on an implied 1 L instead would
    #    lower the premature-neonate saliva curve roughly threefold. See the
    #    vignette Assumptions and deviations.
    Cc <- central / vc
    Csaliva <- saliva / vc

    Cc ~ add(addSd) + prop(propSd)
    Csaliva ~ lnorm(expSd_Csaliva)
  })
}
