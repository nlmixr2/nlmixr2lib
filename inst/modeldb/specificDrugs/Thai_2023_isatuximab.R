Thai_2023_isatuximab <- function() {
  description <- "Joint serum M-protein and progression-free survival (PFS) model for intravenous isatuximab alone or with dexamethasone in adults with relapsed/refractory multiple myeloma (Thai 2023). A Claret tumor-growth-inhibition model of serum M-protein, in which an isatuximab Emax effect and a linear dexamethasone K-PD effect both raise the shrinkage rate that decays with emerging resistance, is linked to a log-logistic PFS hazard through the instantaneous M-protein slope dM/dt. Isatuximab exposure comes from the Fau 2020 two-compartment population PK model, embedded here with its parameters fixed and converted to days."
  reference <- paste(
    "Thai HT, Koiwai K, Shitara Y, Kazama H, Fau JB, Semiond D,",
    "Veyrat-Follet C. Model-based simulation to support the approval of",
    "isatuximab alone or with dexamethasone for the treatment of",
    "relapsed/refractory multiple myeloma in Japanese patients.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(12):1846-1858.",
    "doi:10.1002/psp4.12947. Isatuximab PK layer: Fau JB et al. CPT",
    "Pharmacometrics Syst Pharmacol. 2020;9(11):649-658. doi:10.1002/psp4.12561.",
    sep = " "
  )
  vignette <- "Thai_2023_isatuximab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "isatuximab", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "isatuximab", units = "mg", specimen = "plasma", verified = TRUE),
    depot_kpd = list(analyte = "dexamethasone", units = "mg", specimen = "plasma", verified = TRUE),
    mprotein = list(analyte = "serum M-protein", units = "g/L", specimen = "serum", verified = TRUE),
    cumhaz_pfs = list(
      analyte = "PFS cumulative hazard",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight at baseline",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Isatuximab PK layer only (Fau 2020): power covariate on CLinf (0.621), Vc (0.472), Vp (0.719) and Q (0.477) with reference 75.6 kg. Also needed to convert the mg/kg isatuximab dose to mg.",
      source_name = "WT"
    ),
    B2M = list(
      description = "Baseline serum beta-2-microglobulin",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Two roles. (1) Isatuximab PK layer (Fau 2020): power covariate (B2M/3.90)^0.343 on CLinf. (2) PFS hazard (Thai 2023 Table 3, beta10 = -1.76e-7), encoded as exp(beta10 * log(B2M/4.53)) with 4.53 mg/L the Table 2 median; the coefficient is so small that the PFS term is inert under either a log or a linear covariate form.",
      source_name = "B2MG"
    ),
    MM_NIGG = list(
      description = "Multiple-myeloma immunoglobulin type: 1 = non-IgG MM, 0 = IgG MM",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (IgG MM)",
      notes = "Exponential effect on baseline M-protein M0 (beta1 = -0.255, Thai 2023 Table 3) and, in the Fau 2020 PK layer, on CLinf (-0.751) and KCL (-0.931).",
      source_name = "IGGTYPE"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (ALB/35)^-1.83 on M0. The paper's covariate ALBN is albumin divided by the laboratory upper limit of normal; in a ratio-to-median power form the ULN cancels, so absolute albumin centred on the Table 2 median of 35 g/L is equivalent.",
      source_name = "ALBN"
    ),
    ALPR = list(
      description = "Baseline serum alkaline phosphatase divided by the laboratory upper limit of normal",
      units = "(unitless ratio)",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (ALPR/1)^-0.233 on M0. The analysis-set median is not reported; the reference is the ULN (ALPR = 1).",
      source_name = "ALKN"
    ),
    N_PRIOR_LINES = list(
      description = "Number of prior lines of anticancer therapy",
      units = "(count)",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (N_PRIOR_LINES/4)^0.21 on M0, entered as a continuous count. The analysis-set median is not reported; the reference of 4 prior lines is a rounded value for these heavily pretreated cohorts.",
      source_name = "LINE"
    ),
    REGION_JAPAN = list(
      description = "Japanese study indicator: 1 = enrolled in the Japanese trial NCT02812706, 0 = the non-Japanese trial NCT01084252",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Japanese study)",
      notes = "Exponential effect on M0 (beta5 = -0.216). Every patient of the Japanese trial was Asian and enrolled in Japan.",
      source_name = "Japanese study"
    ),
    REFRACTORY_LEN = list(
      description = "Refractory to lenalidomide: 1 = refractory, 0 = not refractory",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (lenalidomide-refractory)",
      notes = "The paper's effect is on 'REF_Len = N' (not refractory): KL is multiplied by exp(-0.514) when REFRACTORY_LEN = 0. The typical patient of Figure S6 is lenalidomide-refractory.",
      source_name = "REF_LEN"
    ),
    MCPROT = list(
      description = "Observed baseline serum M-protein",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline value. Power effect (MCPROT/28.1)^-0.56 on the M-protein shrinkage rate KD, with 28.1 g/L the Table 2 median. Units are g/L here (the register's default is g/dL). Distinct from the model's M0 parameter, which is the estimated baseline of the M-protein state.",
      source_name = "MPROT"
    ),
    DIS_EMD = list(
      description = "Presence of plasmacytomas at baseline: 1 = yes, 0 = no",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no plasmacytoma)",
      notes = "Exponential effect on the PFS hazard (beta9 = 0.799).",
      source_name = "PCYTOMA"
    ),
    CONMED_DEXAMETHASONE = list(
      description = "Dexamethasone co-administration: 1 = isatuximab plus dexamethasone arm, 0 = isatuximab alone",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (isatuximab alone)",
      notes = "Exponential effect on the PFS hazard (beta8 = -0.363), applied in addition to the dexamethasone exposure effect on M-protein shrinkage, which comes from the doses given into depot_kpd. Time-fixed per subject (randomised arm).",
      source_name = "Dex"
    ),
    FORM_ISA_P2F2 = list(
      description = "Isatuximab drug material: 1 = P2F2 (commercial-bound), 0 = P1F1 (early-phase)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (P1F1)",
      notes = "Isatuximab PK layer only (Fau 2020): exponential effect -0.137 on Vc. P2F2 material was used in phase II stage 2 of NCT01084252 and in NCT02812706 (Appendix S1).",
      source_name = "Drug_mat"
    ),
    RACE_ASIAN = list(
      description = "Race indicator: 1 = Asian, 0 = non-Asian",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian)",
      notes = "Isatuximab PK layer only (Fau 2020): exponential effect -0.275 on Vc.",
      source_name = "Asian"
    ),
    SEXF = list(
      description = "Biological sex indicator: 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Isatuximab PK layer only (Fau 2020): exponential effect -0.126 on Vc.",
      source_name = "Sex"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 201L,
    n_observations = 1774L,
    n_studies = 2L,
    age_range = "37-84 years",
    age_median = "66 years",
    weight_range = "37.6-153 kg",
    weight_median = "70.0 kg",
    sex_female_pct = 46.3,
    race_ethnicity = c(White = 70.6, Black = 6.0, Asian = 15.9, Other = 7.5),
    disease_state = "Relapsed/refractory multiple myeloma with measurable serum M-protein (>= 2 serum M-protein values including baseline). IgG MM 69.2%; plasmacytomas 15.4%; ISS stage I/II/III 28.4/37.8/33.8%; high-risk cytogenetics 26.4%.",
    dose_range = "Isatuximab IV 3 mg/kg q2w, 10 mg/kg q2w, 10 mg/kg q2w-q4w, 10 mg/kg qw-q2w or 20 mg/kg qw-q2w (qw-q2w = weekly for 4 weeks then every 2 weeks); 38 non-Japanese patients also received dexamethasone 40 mg (20 mg if >= 75 years) on days 1, 8, 15 and 22 of each 28-day cycle.",
    regions = "Non-Japanese multinational phase II (NCT01084252, n = 170) and Japanese phase I/II (NCT02812706, n = 31).",
    notes = "Baseline demographics from Thai 2023 Table 2 (median baseline B2M 4.53 mg/L, albumin 35 g/L, serum M-protein 28.1 g/L, eGFR 76.2 mL/min/1.73 m2). Monolix 2020R1 SAEM joint fit; M-protein below 0.5 g/L (4% of data, all from the Japanese study) handled as censored."
  )

  ini({
    # --- Isatuximab PK layer: Fau 2020 Table S3 final estimates, fixed. Thai
    # 2023 used individual post-hoc PK parameters from this model as
    # regressors (Methods; Appendix S1 'Pharmacokinetic model for
    # isatuximab'). Converted from hours to days: CL, Q and Vm x 24, KCL / 24.
    lcl <- fixed(log(0.00955 * 24)); label("Steady-state linear clearance CLinf at reference covariates (L/day)") # Fau 2020 Table S3: CLinf = 0.00955 L/h
    lvc <- fixed(log(5.13)); label("Central volume of distribution Vc at reference covariates (L)") # Fau 2020 Table S3: Vc = 5.13 L
    lq <- fixed(log(0.0432 * 24)); label("Intercompartmental clearance Q at reference covariates (L/day)") # Fau 2020 Table S3: Q = 0.0432 L/h
    lvp <- fixed(log(3.62)); label("Peripheral volume of distribution Vp at reference covariates (L)") # Fau 2020 Table S3: Vp = 3.62 L
    lvmax <- fixed(log(0.136 * 24)); label("Michaelis-Menten maximum elimination rate Vm (mg/L/day)") # Fau 2020 Table S3: Vm = 0.136 mg/L/h
    lkm <- fixed(log(0.300)); label("Michaelis-Menten constant Km (mg/L)") # Fau 2020 Table S3: Km = 0.300 mg/L
    cl_time_max <- fixed(0.664); label("Log-ratio of baseline to steady-state linear clearance, CLm (unitless)") # Fau 2020 Table S3: CLm = 0.664
    lcl_t50 <- fixed(log(1055 / 24)); label("Time to half of the linear-clearance change, KCL (day)") # Fau 2020 Table S3: KCL = 1055 h
    lcl_time_hill <- fixed(log(3.91)); label("Sigmoidicity of the time-varying linear clearance, gamma (unitless)") # Fau 2020 Table S3: gamma = 3.91
    e_wt_cl <- fixed(0.621); label("Power exponent of WT/75.6 on CLinf (unitless)") # Fau 2020 Table S3: CLinf ~ Wght
    e_wt_vc <- fixed(0.472); label("Power exponent of WT/75.6 on Vc (unitless)") # Fau 2020 Table S3: Vc ~ Wght
    e_wt_vp <- fixed(0.719); label("Power exponent of WT/75.6 on Vp (unitless)") # Fau 2020 Table S3: Vp ~ Wght
    e_wt_q <- fixed(0.477); label("Power exponent of WT/75.6 on Q (unitless)") # Fau 2020 Table S3: Q ~ Wght
    e_b2m_cl <- fixed(0.343); label("Power exponent of B2M/3.90 on CLinf (unitless)") # Fau 2020 Table S3: CLinf ~ B2M
    e_nigg_cl <- fixed(-0.751); label("Exponential coefficient of non-IgG MM on CLinf (unitless)") # Fau 2020 Table S3: CLinf ~ IgType=Not_IgG
    e_nigg_cl_t50 <- fixed(-0.931); label("Exponential coefficient of non-IgG MM on KCL (unitless)") # Fau 2020 Table S3: KCL ~ IgType=Not_IgG
    e_p2f2_vc <- fixed(-0.137); label("Exponential coefficient of P2F2 drug material on Vc (unitless)") # Fau 2020 Table S3: Vc ~ Formulation=P2F2
    e_asian_vc <- fixed(-0.275); label("Exponential coefficient of Asian race on Vc (unitless)") # Fau 2020 Table S3: Vc ~ Race=Asian
    e_sexf_vc <- fixed(-0.126); label("Exponential coefficient of female sex on Vc (unitless)") # Fau 2020 Table S3: Vc ~ Sex=Female

    # Fau 2020 Table S3 reports IIV as CV; omega^2 = log(1 + CV^2) for the
    # log-normal parameters and (0.664 * CV)^2 for the normal CLm.
    etalcl ~ fixed(log(1 + 0.475^2)) # Fau 2020 Table S3: CLinf 47.5% CV
    etacl_time_max ~ fixed((0.664 * 0.972)^2) # Fau 2020 Table S3: CLm 97.2% CV (additive normal eta)
    etalcl_t50 ~ fixed(log(1 + 1.15^2)) # Fau 2020 Table S3: KCL 115% CV
    etalcl_time_hill ~ fixed(log(1 + 1.18^2)) # Fau 2020 Table S3: gamma 118% CV
    etalvc ~ fixed(log(1 + 0.257^2)) # Fau 2020 Table S3: Vc 25.7% CV
    etalq ~ fixed(log(1 + 0.858^2)) # Fau 2020 Table S3: Q 85.8% CV
    etalvp ~ fixed(log(1 + 0.456^2)) # Fau 2020 Table S3: Vp 45.6% CV
    etalvmax ~ fixed(log(1 + 0.615^2)) # Fau 2020 Table S3: Vm 61.5% CV
    etalkm ~ fixed(log(1 + 0.889^2)) # Fau 2020 Table S3: Km 88.9% CV

    # --- Dexamethasone K-PD layer (Appendix S1 mlxtran code), fixed.
    lkel <- fixed(log(2.42)); label("Dexamethasone K-PD elimination rate constant KDEd (1/day)") # Appendix S1 mlxtran: KDEd=2.42

    # --- Serum M-protein Claret TGI submodel (Thai 2023 Table 3).
    lrbase <- log(27); label("Baseline serum M-protein M0 for the reference patient (g/L)") # Table 3: M0 = 27 g/L
    e_nigg_rbase <- -0.255; label("Exponential coefficient of non-IgG MM on M0 (unitless)") # Table 3: beta1~non_IgG = -0.255
    e_alb_rbase <- -1.83; label("Power exponent of ALB/35 on M0 (unitless)") # Table 3: beta2~ALBN = -1.83
    e_alpr_rbase <- -0.233; label("Power exponent of ALPR/1 on M0 (unitless)") # Table 3: beta3~ALKN = -0.233
    e_lines_rbase <- 0.21; label("Power exponent of N_PRIOR_LINES/4 on M0 (unitless)") # Table 3: beta4~LINE = 0.21
    e_japan_rbase <- -0.216; label("Exponential coefficient of the Japanese study on M0 (unitless)") # Table 3: beta5~Japanese study = -0.216
    lp <- log(0.00531); label("Serum M-protein growth rate KL for a lenalidomide-refractory patient (1/day)") # Table 3: KL = 0.00531 day-1
    e_nonrefractory_len_p <- -0.514; label("Exponential coefficient of not being lenalidomide-refractory on KL (unitless)") # Table 3: beta6~REF_Len = N = -0.514
    lkkillmax <- log(0.0396); label("Serum M-protein shrinkage rate KD at the reference baseline M-protein (1/day)") # Table 3: KD = 0.0396
    e_mcprot_kkillmax <- -0.56; label("Power exponent of MCPROT/28.1 on KD (unitless)") # Table 3: beta7~MPROT = -0.56
    llambda <- log(0.0192); label("Rate constant of resistance appearance R (1/day)") # Table 3: R = 0.0192 day-1
    lec50 <- log(1.13); label("Isatuximab concentration giving half-maximal shrinkage stimulation EC50 (umol/L)") # Table 3: EC50 = 1.13 (umol/L; text: 170 mg/L)
    lslope_dexamethasone <- log(5.47); label("Linear coefficient k of the dexamethasone K-PD concentration on shrinkage (L/umol)") # Table 3: k = 5.47

    # Monolix omegas are standard deviations of the log-normal random effects
    # (Table 3 'Interindividual variability (%)' = 100 * omega); variance = omega^2.
    etalrbase ~ 0.398^2 # Table 3: omega_M0 = 39.8%
    etalp ~ 0.947^2 # Table 3: omega_KL = 94.7%
    etalkkillmax ~ 0.87^2 # Table 3: omega_KD = 87%
    etallambda ~ 1.13^2 # Table 3: omega_R = 113%
    etalec50 ~ 1.23^2 # Table 3: omega_EC50 = 123%

    addSd <- 0.477; label("Additive residual error on serum M-protein (g/L)") # Table 3: sigma additive = 0.477 g/L
    propSd <- 0.0721; label("Proportional residual error on serum M-protein (fraction)") # Table 3: sigma proportional = 7.21%

    # --- PFS submodel (Thai 2023 Table 3, time-to-event part).
    llambda_pfs <- log(272); label("Log-logistic scale (characteristic time) Te of the baseline PFS hazard (day)") # Table 3: Te = 272
    lgamma_pfs <- log(1.56); label("Log-logistic shape s of the baseline PFS hazard (unitless)") # Table 3: s = 1.56
    e_dexamethasone_pfs <- -0.363; label("Log-hazard coefficient of dexamethasone co-administration on PFS (unitless)") # Table 3: beta8~Dex = Y = -0.363
    e_emd_pfs <- 0.799; label("Log-hazard coefficient of plasmacytoma presence on PFS (unitless)") # Table 3: beta9~PCYTOMA = Y = 0.799
    e_b2m_pfs <- -1.76e-7; label("Log-hazard coefficient of log(B2M/4.53) on PFS (unitless)") # Table 3: beta10~B2MG = -1.76E-07
    e_mprotein_slope_pfs <- 6.65; label("Log-hazard coefficient of the serum M-protein slope dM/dt on PFS (day*L/g)") # Table 3: beta11~SlopeM = 6.65
  })

  model({
    # --- Isatuximab PK layer (Fau 2020), time in days.
    cl <- exp(lcl + etalcl) * (WT / 75.6)^e_wt_cl * (B2M / 3.90)^e_b2m_cl *
      exp(e_nigg_cl * MM_NIGG)
    cl_t50 <- exp(lcl_t50 + etalcl_t50) * exp(e_nigg_cl_t50 * MM_NIGG)
    cl_time_max_i <- cl_time_max + etacl_time_max
    cl_time_hill <- exp(lcl_time_hill + etalcl_time_hill)
    vc <- exp(lvc + etalvc) * (WT / 75.6)^e_wt_vc * exp(e_p2f2_vc * FORM_ISA_P2F2) *
      exp(e_asian_vc * RACE_ASIAN) * exp(e_sexf_vc * SEXF)
    vp <- exp(lvp + etalvp) * (WT / 75.6)^e_wt_vp
    q <- exp(lq + etalq) * (WT / 75.6)^e_wt_q
    vmax <- exp(lvmax + etalvmax)
    km <- exp(lkm + etalkm)

    # Linear CL starts at cl * exp(CLm) and decays to cl; t is time since the
    # first isatuximab dose.
    cllin <- cl * exp(cl_time_max_i * (1 - t^cl_time_hill / (cl_t50^cl_time_hill + t^cl_time_hill)))

    Cc <- central / vc
    d/dt(central) <- -(cllin / vc) * central - vc * vmax * Cc / (km + Cc) -
      (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <- (q / vc) * central - (q / vp) * peripheral1

    # --- Dexamethasone K-PD (Appendix S1): bolus doses (mg) into depot_kpd,
    # virtual volume 76.3 L.
    kel <- exp(lkel)
    d/dt(depot_kpd) <- -kel * depot_kpd
    cdex <- depot_kpd / 76.3

    # Molar concentrations in umol/L, as in the mlxtran code
    # (CiM = Cic/MW_ISA*1000 with Cic in mg/L): molecular weights 150000 g/mol
    # (isatuximab) and 392.46 g/mol (dexamethasone).
    cisa_um <- Cc / 150000 * 1000
    cdex_um <- cdex / 392.46 * 1000

    # --- Serum M-protein (Thai 2023 Methods equation; Table 3 covariates).
    # Continuous covariates enter as log(cov / reference) on the log-normal
    # parameters, categorical ones as exp(beta * indicator).
    rbase <- exp(lrbase + etalrbase) * exp(e_nigg_rbase * MM_NIGG) *
      (ALB / 35)^e_alb_rbase * (ALPR / 1)^e_alpr_rbase *
      (N_PRIOR_LINES / 4)^e_lines_rbase * exp(e_japan_rbase * REGION_JAPAN)
    p <- exp(lp + etalp) * exp(e_nonrefractory_len_p * (1 - REFRACTORY_LEN))
    kkillmax <- exp(lkkillmax + etalkkillmax) * (MCPROT / 28.1)^e_mcprot_kkillmax
    lambda <- exp(llambda + etallambda)
    ec50 <- exp(lec50 + etalec50)
    slope_dexamethasone <- exp(lslope_dexamethasone)

    mprotein(0) <- rbase
    dmdt <- p * mprotein - kkillmax * exp(-lambda * t) *
      (cisa_um / (ec50 + cisa_um) + slope_dexamethasone * cdex_um) * mprotein
    d/dt(mprotein) <- dmdt

    # --- PFS: log-logistic baseline hazard x exp(baseline covariates + link
    # coefficient * dM/dt), with hazard from time 0 (first dose).
    lambda_pfs <- exp(llambda_pfs)
    gamma_pfs <- exp(lgamma_pfs)
    tsc_pfs <- t / lambda_pfs
    hazard0_pfs <- (gamma_pfs / lambda_pfs) * tsc_pfs^(gamma_pfs - 1) / (1 + tsc_pfs^gamma_pfs)
    # The log-hazard is capped at 50 only to avoid floating-point overflow:
    # with M-protein regrowing at several g/L/day, 6.65 * dM/dt exceeds the
    # range of exp(), and at the cap the hazard (> 1e18 /day after the first
    # day) has already driven survival to zero, so the cap does not change it.
    linpred_pfs <- e_dexamethasone_pfs * CONMED_DEXAMETHASONE + e_emd_pfs * DIS_EMD +
      e_b2m_pfs * log(B2M / 4.53) + e_mprotein_slope_pfs * dmdt
    hazard_pfs <- hazard0_pfs * exp(min(linpred_pfs, 50))
    d/dt(cumhaz_pfs) <- hazard_pfs
    surv_pfs <- exp(-cumhaz_pfs)

    mprotein ~ add(addSd) + prop(propSd) + combined1()
  })
}
