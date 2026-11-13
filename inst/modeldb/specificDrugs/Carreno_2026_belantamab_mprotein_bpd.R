Carreno_2026_belantamab_mprotein_bpd <- function() {
  description <- paste(
    "Concentration-driven tumor growth inhibition (TGI) model for serum",
    "M-protein in adults with relapsed/refractory multiple myeloma treated",
    "with belantamab mafodotin plus pomalidomide and dexamethasone (BPd;",
    "DREAMM-8): logistic growth with a modified weak Allee term, an",
    "exponential (saturating) kill term driven by the central belantamab",
    "mafodotin concentration and exponential resistance decay, with no",
    "covariates on the M-protein parameters (Carreno 2026). The",
    "belantamab mafodotin PK driver is the Papathanasiou 2025",
    "two-compartment time-varying-clearance covariate model with the",
    "values printed in the Carreno 2026 mrgsolve code (File S2).",
    sep = " "
  )
  reference <- paste(
    "Carreno F, van Noort M, Taylor A, Zeiser S, Papathanasiou T, Post TM,",
    "Lien YT(K), Ho YL, Ferron-Brady G. Simulation Framework to Investigate",
    "Efficacy and Ocular Safety of Belantamab Mafodotin Combinations in",
    "Relapsed/Refractory Multiple Myeloma. CPT Pharmacometrics Syst",
    "Pharmacol. 2026;15(6):e70276. doi:10.1002/psp4.70276.",
    "PK driver: Papathanasiou T et al. Clin Pharmacokinet.",
    "2025;64(6):925-942. doi:10.1007/s40262-025-01508-1.",
    sep = " "
  )
  vignette <- "Carreno_2026_belantamab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(analyte = "belantamab mafodotin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "belantamab mafodotin", units = "mg", specimen = "plasma", verified = TRUE),
    mprotein = list(analyte = "M-protein", units = "g/L", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    MCPROT = list(
      description = "Baseline serum monoclonal (M) protein concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed baseline in g/L (register default is g/dL; 1 g/dL = 10",
        "g/L). Seeds the initial condition of the `mprotein` state. The",
        "authors' mrgsolve code floors it at 5 g/L (`if(BMPROT < 5) BMPROTBL",
        "= 5`), the inclusion threshold of the analysis dataset; the same",
        "floor is applied here.",
        sep = " "
      ),
      source_name = "BMPROT"
    ),
    IGG = list(
      description = "Baseline serum immunoglobulin G concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed baseline. PK driver only: power exponents on CL (0.170)",
        "and Imax (0.192) with reference 15 g/L (File S2 THETA13, THETA119).",
        "Baseline IgG on the M-protein growth rate was tested in the BPd",
        "model but not retained (95% CI contained the null; Discussion).",
        sep = " "
      ),
      source_name = "IGGBL"
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "PK driver only. Power exponent 0.929 on V1 and V2 and 0.542 on CL and Q; reference 75 kg (File S2 THETA8, THETA9).",
      source_name = "WTBL"
    ),
    ALB = list(
      description = "Baseline serum albumin concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "PK driver only. Power exponents -0.698 on CL, -0.302 on V1 and +0.567 on V2; reference 40 g/L (File S2 THETA10, THETA11, THETA121).",
      source_name = "ALBBL"
    ),
    SBCMA = list(
      description = "Baseline serum soluble B-cell maturation antigen (sBCMA) concentration",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "PK driver only. Reported in ug/L (= ng/mL). Power exponents +0.113 on CL, +0.0401 on V1 and +0.160 on Imax; reference 50 (File S2 THETA12, THETA15, THETA120; SBCMAREF = 50).",
      source_name = "SBCMABL"
    ),
    BMI = list(
      description = "Baseline body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "PK driver only. Power exponent -0.459 on V1; reference 27 kg/m^2 (File S2 THETA14).",
      source_name = "IBMIBL"
    ),
    RACE_ASIAN = list(
      description = "Indicator for Asian race",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian)",
      notes = "PK driver only. Multiplicative factor 0.912 on CL (File S2 THETA116, applied when the code's FDARACEN = 1). Mutually exclusive with RACE_BLACK.",
      source_name = "FDARACEN"
    ),
    RACE_BLACK = list(
      description = "Indicator for Black / African American race",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Black)",
      notes = "PK driver only. Multiplicative factor 0.861 on CL (File S2 THETA117, applied when the code's FDARACEN = 2). Mutually exclusive with RACE_ASIAN.",
      source_name = "FDARACEN"
    ),
    COMBO_BELAMAF = list(
      description = "Belantamab mafodotin combination therapy indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (belantamab mafodotin monotherapy)",
      notes = "PK driver only. Multiplicative factor 1.44 on Imax of the time-varying clearance (File S2 THETA118, applied when the code's COMBO >= 1). Every patient this model was fitted to received BPd, so COMBO_BELAMAF = 1 for its intended use.",
      source_name = "COMBO"
    )
  )

  covariatesDataExcluded <- list(
    B2M = list(
      description = "Baseline serum beta-2-microglobulin concentration",
      units = "nmol/L",
      type = "continuous",
      notes = "Screened in the M-protein covariate analysis but not retained (Methods S1); median 231 nmol/L in the BPd analysis population (Table 1)."
    ),
    LDH = list(
      description = "Baseline serum lactate dehydrogenase",
      units = "U/L",
      type = "continuous",
      notes = "Screened in the M-protein covariate analysis but not retained (Methods S1); median 177 U/L (Table 1)."
    ),
    DIS_EMD = list(
      description = "Extramedullary disease at screening",
      units = "(binary)",
      type = "binary",
      notes = "Screened in the M-protein covariate analysis but not retained (Methods S1); present in 16% of the BPd analysis population (Table S3)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 103L,
    n_studies = 1L,
    age_range = "40-82 years",
    age_median = "67 years",
    weight_range = "46.0-117 kg",
    weight_median = "77.5 kg",
    sex_female_pct = 34,
    disease_state = paste(
      "Relapsed/refractory multiple myeloma after >= 1 prior line of",
      "therapy including lenalidomide, followed for response by serum",
      "M-protein, with baseline serum M-protein >= 5 g/L.",
      sep = " "
    ),
    dose_range = paste(
      "Belantamab mafodotin IV 2.5 mg/kg on Day 1 of Cycle 1, then 1.9",
      "mg/kg Q4W from Cycle 2, with pomalidomide and dexamethasone",
      "(Table S1).",
      sep = " "
    ),
    co_medication = "Pomalidomide + dexamethasone (BPd) in all patients.",
    biomarker_summary = paste(
      "Median (range) at baseline (Table 1): albumin 41.0 (25.0-49.0) g/L;",
      "sBCMA 40.0 (0.975-741) ug/L; IgG 18.0 (0.800-96.3) g/L;",
      "beta-2-microglobulin 231 (121-1340) nmol/L; LDH 177 (68.6-1180) U/L;",
      "M-protein 17.5 (5.30-83.7) g/L. Myeloma type IgG 69%, IgA 30%; ISS",
      "I/II/III 59/24/16% (Table S3).",
      sep = " "
    ),
    notes = paste(
      "1173 serum M-protein concentrations from 103 patients in DREAMM-8",
      "(NCT04484623) (Results; Table S2).",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Belantamab mafodotin PK driver: Papathanasiou 2025 popPK model,
    # values as printed in the Carreno 2026 mrgsolve code (File S2,
    # $PARAM block, THETA1-THETA15 and THETA116-THETA121). The M-protein
    # model was fitted to individual post hoc PK predictions from this
    # model (Methods).
    # ------------------------------------------------------------------
    lcl <- log(0.925475); label("Clearance at time zero for the reference patient (L/day)") # File S2 THETA1 = 0.925475 (TVCL)
    lvc <- log(4.21296); label("Central volume of distribution V1 for the reference patient (L)") # File S2 THETA2 = 4.21296 (TVV1)
    lq <- log(0.710541); label("Intercompartmental clearance Q for the reference patient (L/day)") # File S2 THETA3 = 0.710541 (TVQ)
    lvp <- log(6.6294); label("Peripheral volume of distribution V2 for the reference patient (L)") # File S2 THETA4 = 6.6294 (TVV2)
    cl_time_max <- -0.402758; label("Maximal log-fold change in clearance over time, Imax, monotherapy reference (unitless)") # File S2 THETA5 = -0.402758 (TVIMAX)
    lcl_t50 <- log(66.3585); label("Time of half-maximal change in clearance, TI50 (day)") # File S2 THETA6 = 66.3585 (TVTI50)
    lcl_time_hill <- log(2.86905); label("Hill exponent of time on clearance, gamma (log unitless)") # File S2 THETA7 = 2.86905 (TVGAMMA)

    e_wt_vc_vp <- 0.928527; label("Power exponent of WT on V1 and V2 (unitless)") # File S2 THETA8 = 0.928527 (WTBLonV)
    e_wt_cl_q <- 0.54204; label("Power exponent of WT on CL and Q (unitless)") # File S2 THETA9 = 0.54204 (WTBLonCL)
    e_alb_cl <- -0.698039; label("Power exponent of ALB on CL (unitless)") # File S2 THETA10 = -0.698039 (ALBBLonCL)
    e_alb_vc <- -0.301545; label("Power exponent of ALB on V1 (unitless)") # File S2 THETA11 = -0.301545 (ALBBLonV1)
    e_sbcma_cl <- 0.112764; label("Power exponent of SBCMA on CL (unitless)") # File S2 THETA12 = 0.112764 (SBCMABLonCL)
    e_igg_cl <- 0.169878; label("Power exponent of IGG on CL (unitless)") # File S2 THETA13 = 0.169878 (IGGBLonCL)
    e_bmi_vc <- -0.45933; label("Power exponent of BMI on V1 (unitless)") # File S2 THETA14 = -0.45933 (V1IBMIBL1)
    e_sbcma_vc <- 0.0400974; label("Power exponent of SBCMA on V1 (unitless)") # File S2 THETA15 = 0.0400974 (V1SBCMABL1)
    e_race_asian_cl <- 0.912487; label("Multiplicative factor on CL for Asian race (unitless)") # File S2 THETA116 = 0.912487 (CLFDARACEN1)
    e_race_black_cl <- 0.861015; label("Multiplicative factor on CL for Black race (unitless)") # File S2 THETA117 = 0.861015 (CLFDARACEN2)
    e_combo_imax <- 1.44166; label("Multiplicative factor on Imax for combination therapy (unitless)") # File S2 THETA118 = 1.44166 (IMAXCOMBO1)
    e_igg_imax <- 0.191848; label("Power exponent of IGG on Imax (unitless)") # File S2 THETA119 = 0.191848 (IMAXIGGBL1)
    e_sbcma_imax <- 0.15972; label("Power exponent of SBCMA on Imax (unitless)") # File S2 THETA120 = 0.15972 (IMAXSBCMABL1)
    e_alb_vp <- 0.566902; label("Power exponent of ALB on V2 (unitless)") # File S2 THETA121 = 0.566902 (V2ALBBL1)

    # ------------------------------------------------------------------
    # Serum M-protein TGI, BPd model (Table S4 'BPd model' column). File
    # S2 implements the BVd model; the final BVd model was the base model
    # for the BPd model (Methods), so the same parameterisation is used.
    # KGR, KDEATH and LAM are reported on the natural-log scale. E50
    # (the slope S of the kill term) is also a log-scale THETA: File S2
    # reads 'LNS = THETA69; // LNE50 in run82; S = exp(LNS + ETAS)'.
    # A50 is linear (File S2: MU_6 = log(TVAL)).
    # ------------------------------------------------------------------
    lkg <- -5.18; label("Log M-protein growth rate constant KGR (log 1/day)") # Table S4 BPd: TVKGR = -5.18 (3.94 %RSE)
    lkd <- -3.23; label("Log M-protein kill rate constant KDEATH (log 1/day)") # Table S4 BPd: TVKDEATH = -3.23 (3.24 %RSE)
    llambda <- -5.25; label("Log resistance rate constant LAM (log 1/day)") # Table S4 BPd: TVLAM = -5.25 (2.18 %RSE)
    lslope <- 0.432; label("Log slope S of the kill term 1 - exp(-S * C) (log mL/ug)") # Table S4 BPd: TVE50 = 0.432 (12.7 %RSE)
    la50 <- log(0.996); label("M-protein concentration halving the growth rate, A50 (g/L)") # Table S4 BPd: TVA50 = 0.996 g/L (7.02 %RSE)
    ltsmax <- fixed(log(150)); label("Logistic carrying capacity of serum M-protein (g/L)") # File S2 $DES: KGR * (1 - A(4)/150); not estimated

    # ------------------------------------------------------------------
    # Inter-individual variability. PK etas are not reported by Carreno
    # 2026 (the framework used individual post hoc PK estimates); they
    # are carried from the Papathanasiou 2025 model, as packaged in
    # modellib('Papathanasiou_2025_belantamab').
    # ------------------------------------------------------------------
    etalcl + etalvc ~ c(
      0.06593,
      0.0328, 0.03922
    ) # carried from Papathanasiou_2025_belantamab (Papathanasiou 2025 Table 2: CL CV 26.1%, Vc CV 20.0%, covariance 0.0328)
    etalq ~ 0.03546 # carried from Papathanasiou_2025_belantamab (Papathanasiou 2025 Table 2: Q CV 19.0%)
    etalvp ~ 0.08945 # carried from Papathanasiou_2025_belantamab (Papathanasiou 2025 Table 2: Vp CV 30.6%)
    etacl_time_max ~ 0.01366 # carried from Papathanasiou_2025_belantamab (Papathanasiou 2025 Table 2: Imax CV 29.0%, additive)
    etalcl_t50 ~ 0.39236 # carried from Papathanasiou_2025_belantamab (Papathanasiou 2025 Table 2: TI50 CV 69.3%)

    # M-protein IIV (Table S4 BPd). The 3x3 block is listed in NONMEM
    # lower-triangle order; the fourth entry is printed as a second
    # 'Covariance (KGR~KDEATH)' but sits in the KGR~LAM position.
    etalkg + etalkd + etallambda ~ c(
      1.45,
      0.533, 0.816,
      0.0457, 0.745, 0.947
    ) # Table S4 BPd: omega KGR 1.45, cov KGR-KDEATH 0.533, omega KDEATH 0.816, cov KGR-LAM 0.0457, cov KDEATH-LAM 0.745, omega LAM 0.947
    etalslope ~ fixed(0.122) # Table S4 BPd: omega E50 = 0.122, held constant
    etala50 ~ fixed(0.122) # Table S4 BPd: omega A50 = 0.122, held constant

    # Residual error. Table S4 heads the proportional row 'variance', but
    # File S2 carries it as a THETA beside the additive SD THETA, in the
    # W = sqrt(PROP^2 * IPRED^2 + ADD^2) form of its commented WBAS line,
    # so the BPd value 0.414 is likewise used as the proportional SD.
    propSd <- 0.2516; label("PK proportional residual error (fraction)") # carried from Papathanasiou_2025_belantamab (Papathanasiou 2025 Table 2: sqrt(0.0633))
    propSd_mprotein <- 0.414; label("M-protein proportional residual error (fraction)") # Table S4 BPd: residual proportional error = 0.414
    addSd_mprotein <- fixed(0.0001); label("M-protein additive residual error (g/L)") # Table S4 BPd: residual error FIXED, additive SD = 0.0001 g/L
  })

  model({
    # 1. PK driver (File S2 $PK). Imax takes its covariates before the
    #    additive eta; COMBO_BELAMAF = 1 for BPd.
    cl_time_max_i <- cl_time_max * e_combo_imax^COMBO_BELAMAF *
      (IGG / 15)^e_igg_imax * (SBCMA / 50)^e_sbcma_imax + etacl_time_max
    cl_t50_i <- exp(lcl_t50 + etalcl_t50)
    cl_time_hill <- exp(lcl_time_hill)
    cl_base <- exp(lcl + etalcl) * (WT / 75)^e_wt_cl_q * (ALB / 40)^e_alb_cl *
      (SBCMA / 50)^e_sbcma_cl * (IGG / 15)^e_igg_cl *
      e_race_asian_cl^RACE_ASIAN * e_race_black_cl^RACE_BLACK
    cl <- cl_base * exp(cl_time_max_i * t^cl_time_hill / (cl_t50_i^cl_time_hill + t^cl_time_hill))
    vc <- exp(lvc + etalvc) * (WT / 75)^e_wt_vc_vp * (ALB / 40)^e_alb_vc *
      (BMI / 27)^e_bmi_vc * (SBCMA / 50)^e_sbcma_vc
    q <- exp(lq + etalq) * (WT / 75)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 75)^e_wt_vc_vp * (ALB / 40)^e_alb_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 2. M-protein parameters (File S2 TUMOR MODEL block, without the BVd
    #    IgG covariate).
    kg <- exp(lkg + etalkg)
    kd <- exp(lkd + etalkd)
    lambda <- exp(llambda + etallambda)
    slope <- exp(lslope + etalslope)
    a50 <- exp(la50 + etala50)
    tsmax <- exp(ltsmax)

    # 3. ODEs. Dose in mg and volume in L give Cc in mg/L = ug/mL, the
    #    reciprocal unit of the slope. Growth is logistic with a weak Allee
    #    term mprotein / (mprotein + A50); the kill term is
    #    KDEATH * (1 - exp(-S * Cc)) * resistance, where File S2 integrates
    #    resistance as a state starting at 1 with d/dt = -LAM * state, i.e.
    #    exp(-LAM * t), written here in closed form.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc

    allee <- mprotein / (mprotein + a50)
    kill <- kd * (1 - exp(-slope * Cc)) * exp(-lambda * t)
    d/dt(mprotein) <- kg * (1 - mprotein / tsmax) * allee * mprotein - kill * mprotein

    # 4. Initial M-protein: baseline value floored at 5 g/L (File S2).
    mprot_bl <- MCPROT
    if (MCPROT < 5) mprot_bl <- 5
    mprotein(0) <- mprot_bl

    # 5. Observations.
    Cc ~ prop(propSd)
    mprotein ~ prop(propSd_mprotein) + add(addSd_mprotein)
  })
}
