Kong_2025b_piperacillin_tazobactam <- function() {
  description <- paste(
    "Joint four-compartment population PK model for intravenous",
    "piperacillin and tazobactam pooled from 13 studies spanning premature",
    "neonates to very elderly adults, including critically ill, obese,",
    "underweight and haematological-malignancy patients (Kong 2025,",
    "general-purpose model). Each drug has a two-compartment disposition",
    "model with linear elimination. Volumes scale linearly and clearance",
    "with exponent 0.75 on total body weight (70 kg reference); the",
    "inter-compartmental clearances use compartmental allometry on the",
    "individual peripheral volume. Elimination clearance of both drugs",
    "carries a shared sigmoidal postmenstrual-age maturation function, a",
    "drug-specific sigmoidal age-related decline (shared shape factor),",
    "and a shared exponential serum-creatinine effect centred on a",
    "postmenstrual-age-standardised serum creatinine. Random effects on",
    "the central volume, peripheral volume and inter-compartmental",
    "clearance are shared between the two drugs. Three observables:",
    "total piperacillin (plasma or dried blood spot), unbound piperacillin",
    "(fraction unbound 0.645) and total tazobactam (plasma or dried blood",
    "spot); study-specific effects for the Sime 2014 haematological-",
    "malignancy cohort (piperacillin CL and V2) and the Sukarnjanaset 2019",
    "cohort (unbound fraction 1)."
  )
  reference <- paste(
    "Kong D, Koomen JV, Colin PJ, et al. A Pooled Pharmacokinetic Analysis",
    "for Piperacillin/Tazobactam Across Different Patient Populations: From",
    "Premature Infants to the Elderly. Clin Pharmacokinet.",
    "2025;64(1):107-126. doi:10.1007/s40262-024-01460-6.",
    "Final NONMEM control stream and estimates: Electronic Supplementary",
    "Material 3 (ESM3).",
    sep = " "
  )
  vignette <- "Kong_2025b_piperacillin_tazobactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "piperacillin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "piperacillin", units = "mg", specimen = "plasma", verified = TRUE),
    central_taz = list(analyte = "tazobactam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_taz = list(analyte = "tazobactam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  dosing <- c("central", "central_taz")

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric size factor FSIZE = WT/70 (Eq. 16; ESM3 $PK 'FSIZE =",
        "WT/70'), exponent 1 on the volumes and 0.75 on the elimination",
        "clearances (Methods, 'Weight-Based and Compartmental Allometry').",
        "Treated as time-varying in the source analysis (Table 1",
        "footnote a). Pooled range 0.5-200 kg (Table 1)."
      ),
      source_name = "WT"
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives the maturation (Eq. 17) and decline (Eqs. 18-19) functions",
        "on the elimination clearance of both drugs and the standardised",
        "serum creatinine of Eq. 20. ESM3 works in postmenstrual YEARS",
        "throughout ('RPMA = 35+(40/52)'), so model() converts with",
        "PAGE / 12. For non-neonates (PMA >= 0.87 years) Kong 2025 took",
        "PMA as postnatal age plus 40 weeks (Methods, 'Component",
        "Datasets'; Table 1 footnote b). Time-varying in the source",
        "analysis. Pooled range 0.5-94.8 postmenstrual years (Table 1)."
      ),
      source_name = "PMA"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying exponential effect on the elimination clearance of",
        "both drugs, F_SCR = exp(-theta_SCR * (CREAT - SCR_std)) (Eq. 4),",
        "with SCR_std the postmenstrual-age-standardised value of Eq. 20.",
        "Units pinned by the dL/mg units of theta_SCR (Table 3).",
        "Missing values: the source coded subjects with no serum",
        "creatinine records as 0 and treated them as SCR_std, i.e.",
        "F_SCR = 1 (Methods; ESM3 'IF(SCR.EQ.0) FSCR = 1'); partly",
        "missing records were filled by next-value-carried-backward. To",
        "reproduce the no-record case, supply CREAT equal to SCR_std for",
        "the subject's postmenstrual age (the model() variable",
        "creat_std). Pooled range 0.17-5.82 mg/dL (Table 1)."
      ),
      source_name = "SCR"
    ),
    SAMPLE_DBS = list(
      description = "Per-observation dried-blood-spot sample indicator (1 = dried blood spot; 0 = plasma)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (plasma)",
      notes = paste(
        "Multiplies the total piperacillin prediction by f_DBS_PIP = 0.368",
        "and the total tazobactam prediction by f_DBS_TAZ = 0.448",
        "(Eqs. 24-25; ESM3 'IF(TDV.EQ.5) fDBS_PIP' / 'IF(TDV.EQ.6)",
        "fDBS_TAZ'). Dried-blood-spot samples share the residual error of",
        "the total plasma observations (Table 3, sigma_TOT&DBS). Only the",
        "Cohen-Wolkowiez 2014 infant study contributed dried blood spots",
        "(55 piperacillin and 55 tazobactam samples). Not applied to the",
        "unbound piperacillin observable Cu."
      ),
      source_name = "TDV (5 = PIP DBS, 6 = TAZ DBS)"
    ),
    STUDY_SIME2014 = list(
      description = "Sime 2014 febrile-neutropenia haematological-malignancy study indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all other pooled studies)",
      notes = paste(
        "Subject-level. 1 multiplies the piperacillin elimination",
        "clearance by theta_CL_Sime = 1.73 and the piperacillin peripheral",
        "volume by theta_V2_Sime = 0.512 (Eqs. 21-22; ESM3 'IF(SID.EQ.4)').",
        "Tazobactam is unaffected. Kong 2025 attributes the effect to the",
        "haematological-malignancy population of that study (12 patients)."
      ),
      source_name = "SID = 4"
    ),
    STUDY_SUKARNJANASET2019 = list(
      description = "Sukarnjanaset 2019 early-sepsis study indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all other pooled studies)",
      notes = paste(
        "Subject-level. 1 sets the piperacillin unbound fraction to",
        "f_UNB_Sukarnjanaset = 1 (fixed) instead of 0.645, so the unbound",
        "observable equals the total concentration (Eq. 23; ESM3",
        "'IF(SID.EQ.8)'). Kong 2025 treats this study as an outlier in",
        "protein binding (Discussion, limitations)."
      ),
      source_name = "SID = 8"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 415,
    n_studies = 13,
    age_range = "Postnatal age 1 day to 94 years; postmenstrual age 0.5-94.8 years",
    weight_range = "0.5-200 kg",
    sex_female_pct = 41,
    disease_state = paste(
      "Pooled individual data from 13 published piperacillin/tazobactam PK",
      "studies: critically ill adults with sepsis or pulmonary infection,",
      "critically ill children, premature and term infants with suspected",
      "systemic infection, infants and children aged 2 months to 6 years",
      "with infection, critically ill obese and morbidly obese adults,",
      "Australian indigenous patients with sepsis, and febrile neutropenic",
      "patients with haematological malignancy. Patients on renal",
      "replacement therapy or extracorporeal membrane oxygenation were",
      "excluded."
    ),
    dose_range = "Intravenous piperacillin/tazobactam (8:1) by intermittent, extended or continuous infusion as per the component studies",
    regions = "Multinational; pooled from 13 published studies (Table 1), countries not tabulated by Kong 2025",
    notes = paste(
      "415 patients by the Results text: 32 newborns, 127 children and",
      "adolescents, 163 adults (42 obese, 12 underweight), 74 elderly and",
      "19 very elderly. The per-study N column of Table 1 sums to 435.",
      "Female percentage derived from the per-study '% Males' column of",
      "Table 1 (about 256 males of 435). 3798 piperacillin observations",
      "(2855 total plasma, 888 unbound plasma, 55 dried blood spot) and",
      "1948 tazobactam observations (1893 total plasma, 55 dried blood",
      "spot). Estimation by FOCE-I in NONMEM 7.5."
    )
  )

  ini({
    # Structural parameters for a 70 kg, 35-year-old (35 + 40/52
    # postmenstrual years) reference with SCR = SCR_std (Table 3; ESM3
    # $THETA final estimates on the log scale).
    lvc <- log(10.4); label("Piperacillin central volume V1 at the reference (L)")                     # Table 3 theta_V1 PIP 10.4 L/70 kg; ESM3 TH1 2.34
    lcl <- log(10.6); label("Piperacillin elimination clearance CL at the reference (L/h)")             # Table 3 theta_CL PIP 10.6 L/h/70 kg; ESM3 TH2 2.36
    lvp <- log(11.6); label("Piperacillin peripheral volume V2 at the reference (L)")                  # Table 3 theta_V2 PIP 11.6 L/70 kg; ESM3 TH3 2.45
    lq <- log(15.2); label("Piperacillin inter-compartmental clearance Q2 at the reference (L/h)")     # Table 3 theta_Q2 PIP 15.2 L/h/70 kg; ESM3 TH4 2.72

    lvc_taz <- log(10.5); label("Tazobactam central volume V1 at the reference (L)")                   # Table 3 theta_V1 TAZ 10.5 L/70 kg; ESM3 TH5 2.35
    lcl_taz <- log(9.58); label("Tazobactam elimination clearance CL at the reference (L/h)")          # Table 3 theta_CL TAZ 9.58 L/h/70 kg; ESM3 TH6 2.26
    lvp_taz <- log(13.7); label("Tazobactam peripheral volume V2 at the reference (L)")                # Table 3 theta_V2 TAZ 13.7 L/70 kg; ESM3 TH7 2.62
    lq_taz <- log(16.8); label("Tazobactam inter-compartmental clearance Q2 at the reference (L/h)")   # Table 3 theta_Q2 TAZ 16.8 L/h/70 kg; ESM3 TH8 2.82

    # Allometric exponents, fixed by theory (Methods, 'Weight-Based and
    # Compartmental Allometry'; Eqs. 1, 6-9, 11-14).
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on elimination clearance (unitless)")        # Methods: exponent 0.75 for clearance terms; Eq. 7 FSIZE^0.75
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on volumes (unitless)")                          # Methods: exponent 1 for volume terms; Eq. 6 FSIZE
    e_wt_q <- fixed(0.75); label("Compartmental-allometry exponent of peripheral volume on Q2 (unitless)")         # Eq. 1 and Eq. 9 (V2/theta_V2)^0.75

    # Maturation and decline of elimination clearance (shared shape
    # factors; drug-specific decline half-age).
    lhill_mat <- log(3.35); label("Hill shape factor gamma1 for clearance maturation (log unitless)")              # Table 3 gamma1 = 3.35; ESM3 TH9 1.21
    ltmat50 <- log(54.2 / 52); label("Postmenstrual age at 50% clearance maturation MAT50 (log years)")            # Table 3 MAT50 = 54.2 weeks; /52 weeks per year as in ESM3 'RPMA = 35+(40/52)'; ESM3 TH10 0.0359 (exp = 1.037 y)
    lhill_dec <- log(1.92); label("Hill shape factor gamma2 for clearance decline (log unitless)")                 # Table 3 gamma2 = 1.92; ESM3 TH11 0.653
    ltdec50 <- log(89.1); label("Postmenstrual age at 50% decline of piperacillin clearance DEC50 (log years)")     # Table 3 DEC50 PIP = 89.1 years; ESM3 TH12 4.49
    ltdec50_taz <- log(61.6); label("Postmenstrual age at 50% decline of tazobactam clearance DEC50 (log years)")  # Table 3 DEC50 TAZ = 61.6 years; ESM3 TH13 4.12

    # Serum creatinine effect, shared between the two drugs.
    e_creat_cl <- 0.346; label("Exponential serum-creatinine coefficient on elimination clearance (dL/mg)")       # Table 3 theta_SCR = 0.346 dL/mg; ESM3 TH14 -1.06 (exp = 0.3465)

    # Piperacillin fraction unbound and dried-blood-spot ratios.
    lfu <- log(0.645); label("Piperacillin fraction unbound in plasma (log fraction)")                             # Table 3 fUNB = 0.645; ESM3 TH15 -0.439
    lfu_sukarnjanaset <- fixed(log(1)); label("Piperacillin fraction unbound in the Sukarnjanaset 2019 study (log fraction)") # Table 3 fUNB_Sukarnjanaset = 1 FIX; ESM3 TH16 '0 FIX'
    lfdbs <- log(0.368); label("Ratio of piperacillin dried-blood-spot to plasma concentration (log fraction)")    # Table 3 fDBS PIP = 0.368; ESM3 TH17 -1.00
    lfdbs_taz <- log(0.448); label("Ratio of tazobactam dried-blood-spot to plasma concentration (log fraction)")  # Table 3 fDBS TAZ = 0.448; ESM3 TH18 -0.802

    # Sime 2014 study effects on piperacillin (Eqs. 21-22).
    e_sime2014_cl <- 1.73; label("Multiplicative piperacillin clearance factor in the Sime 2014 study (unitless)")  # Table 3 theta_CL_Sime = 1.73; ESM3 TH19 0.548
    e_sime2014_vp <- 0.512; label("Multiplicative piperacillin peripheral volume factor in the Sime 2014 study (unitless)") # Table 3 theta_V2_Sime = 0.512; ESM3 TH20 -0.669

    # IIV: log-scale variances from the ESM3 final OMEGA; Table 3 reports
    # CV% = sqrt(exp(omega) - 1). etalvc, etalvp and etalq are each shared
    # by both drugs (ETA1, ETA4, ETA5 in ESM3).
    etalvc ~ 0.167     # ESM3 OMEGA ETA1 0.167, shared V1 PIP and TAZ; Table 3 IIV of V1 42.6%
    etalcl ~ 0.171     # ESM3 OMEGA ETA2 0.171; Table 3 IIV of CL PIP 43.2%
    etalcl_taz ~ 0.159 # ESM3 OMEGA ETA3 0.159; Table 3 IIV of CL TAZ 41.5%
    etalvp ~ 0.548     # ESM3 OMEGA ETA4 0.548, shared V2 PIP and TAZ; Table 3 IIV of V2 85.4%
    etalq ~ 0.358      # ESM3 OMEGA ETA5 0.358, shared Q2 PIP and TAZ; Table 3 IIV of Q2 65.6%

    # Residual error: Y = IPRED * (1 + EPS_prop) + EPS_add (ESM3 $ERROR),
    # SDs = sqrt of the ESM3 final SIGMA (Table 3 footnote c).
    propSd <- 0.302; label("Proportional residual SD, total piperacillin in plasma and DBS (fraction)")   # Table 3 sigma_TOT&DBS prop PIP 30.2%; ESM3 SIGMA EPS1 0.0910
    addSd <- 0.147; label("Additive residual SD, total piperacillin in plasma and DBS (mg/L)")            # Table 3 sigma_TOT&DBS add PIP 0.147 mg/L; ESM3 SIGMA EPS2 0.0216
    propSd_Cu <- 0.365; label("Proportional residual SD, unbound piperacillin (fraction)")                 # Table 3 sigma_UNB prop 36.5%; ESM3 SIGMA EPS3 0.133
    addSd_Cu <- 0.747; label("Additive residual SD, unbound piperacillin (mg/L)")                         # Table 3 sigma_UNB add 0.747 mg/L; ESM3 SIGMA EPS4 0.557
    propSd_taz <- 0.285; label("Proportional residual SD, total tazobactam in plasma and DBS (fraction)") # Table 3 sigma_TOT&DBS prop TAZ 28.5%; ESM3 SIGMA EPS5 0.0810 (additive TAZ 0 FIX, omitted)
  })

  model({
    # Postmenstrual age in years and the 35-year reference (ESM3 'RPMA =
    # 35+(40/52)').
    pma <- PAGE / 12
    rpma <- 35 + 40 / 52

    fsize <- WT / 70

    # Maturation (Eq. 17) and decline (Eqs. 18-19), each normalised to the
    # reference individual as in ESM3 (FMAT/RFMAT, FDEC/RFDEC). The
    # printed equations omit the normalisation; the control stream is the
    # implemented model and makes theta_CL the reference-patient value.
    hill_mat <- exp(lhill_mat)
    tmat50 <- exp(ltmat50)
    fmat <- pma^hill_mat / (pma^hill_mat + tmat50^hill_mat)
    fmat_ref <- rpma^hill_mat / (rpma^hill_mat + tmat50^hill_mat)

    hill_dec <- exp(lhill_dec)
    tdec50 <- exp(ltdec50)
    tdec50_taz <- exp(ltdec50_taz)
    fdec <- 1 - pma^hill_dec / (pma^hill_dec + tdec50^hill_dec)
    fdec_ref <- 1 - rpma^hill_dec / (rpma^hill_dec + tdec50^hill_dec)
    fdec_taz <- 1 - pma^hill_dec / (pma^hill_dec + tdec50_taz^hill_dec)
    fdec_taz_ref <- 1 - rpma^hill_dec / (rpma^hill_dec + tdec50_taz^hill_dec)

    # Standardised serum creatinine (Eq. 20; ESM3 SCRstd), PMA in years;
    # 1.42, 1.17 and 0.203 are published constants, not fitted here.
    creat_std <- exp(1.42 - (1.17 + 0.203 * log(pma / 100)) / sqrt(pma / 100))
    f_scr <- exp(-e_creat_cl * (CREAT - creat_std))

    f_sime_cl <- 1 + (e_sime2014_cl - 1) * STUDY_SIME2014
    f_sime_vp <- 1 + (e_sime2014_vp - 1) * STUDY_SIME2014

    # Piperacillin (Eqs. 6-9)
    vc <- exp(lvc + etalvc) * fsize^e_wt_vc
    cl <- exp(lcl + etalcl) * fsize^e_wt_cl * (fmat / fmat_ref) * (fdec / fdec_ref) *
      f_scr * f_sime_cl
    vp <- exp(lvp + etalvp) * fsize^e_wt_vc * f_sime_vp
    q <- exp(lq + etalq) * (vp / exp(lvp))^e_wt_q

    # Tazobactam (Eqs. 11-14)
    vc_taz <- exp(lvc_taz + etalvc) * fsize^e_wt_vc
    cl_taz <- exp(lcl_taz + etalcl_taz) * fsize^e_wt_cl * (fmat / fmat_ref) *
      (fdec_taz / fdec_taz_ref) * f_scr
    vp_taz <- exp(lvp_taz + etalvp) * fsize^e_wt_vc
    q_taz <- exp(lq_taz + etalq) * (vp_taz / exp(lvp_taz))^e_wt_q

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    kel_taz <- cl_taz / vc_taz
    k12_taz <- q_taz / vc_taz
    k21_taz <- q_taz / vp_taz

    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(central_taz) <- -(kel_taz + k12_taz) * central_taz + k21_taz * peripheral1_taz
    d/dt(peripheral1_taz) <- k12_taz * central_taz - k21_taz * peripheral1_taz

    # Observables (Eqs. 10, 15, 23-25). Dose in mg and volumes in L give
    # mg/L. Cc and Cc_taz are plasma concentrations, or dried-blood-spot
    # concentrations on records with SAMPLE_DBS = 1.
    fu <- exp(lfu) * (1 - STUDY_SUKARNJANASET2019) + exp(lfu_sukarnjanaset) * STUDY_SUKARNJANASET2019
    fdbs <- exp(lfdbs)
    fdbs_taz <- exp(lfdbs_taz)

    Cc <- central / vc * (1 + (fdbs - 1) * SAMPLE_DBS)
    Cu <- central / vc * fu
    Cc_taz <- central_taz / vc_taz * (1 + (fdbs_taz - 1) * SAMPLE_DBS)

    Cc ~ add(addSd) + prop(propSd)
    Cu ~ add(addSd_Cu) + prop(propSd_Cu)
    Cc_taz ~ prop(propSd_taz)
  })
}
