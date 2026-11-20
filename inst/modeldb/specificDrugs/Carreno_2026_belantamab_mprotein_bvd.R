Carreno_2026_belantamab_mprotein_bvd <- function() {
  description <- paste(
    "Concentration-driven tumor growth inhibition (TGI) model for serum",
    "M-protein in adults with relapsed/refractory multiple myeloma treated",
    "with belantamab mafodotin plus bortezomib and dexamethasone (BVd;",
    "DREAMM-6 Arm B and DREAMM-7): logistic growth with a modified weak",
    "Allee term, an exponential (saturating) kill term driven by the",
    "central belantamab mafodotin concentration, exponential resistance",
    "decay, and a baseline-IgG effect on the kill rate constant (Carreno",
    "2026). The belantamab mafodotin PK driver is the Papathanasiou 2025",
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
        "Time-fixed baseline. Two roles: (1) M-protein kill rate constant,",
        "log(KDEATH) shifted by -0.0884 * (log(IGG) - 2.74), a power model",
        "centred on exp(2.74) = 15.5 g/L (Table S4; File S2 THETA71); (2)",
        "PK driver, power exponents on CL (0.170) and Imax (0.192) with",
        "reference 15 g/L (File S2 THETA13, THETA119).",
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
      notes = "PK driver only. Multiplicative factor 1.44 on Imax of the time-varying clearance (File S2 THETA118, applied when the code's COMBO >= 1). Every patient this model was fitted to received BVd, so COMBO_BELAMAF = 1 for its intended use.",
      source_name = "COMBO"
    )
  )

  covariatesDataExcluded <- list(
    B2M = list(
      description = "Baseline serum beta-2-microglobulin concentration",
      units = "nmol/L",
      type = "continuous",
      notes = "Screened in the M-protein covariate analysis but not retained (Methods S1); median 263 nmol/L in the BVd analysis population (Table 1)."
    ),
    LDH = list(
      description = "Baseline serum lactate dehydrogenase",
      units = "U/L",
      type = "continuous",
      notes = "Screened in the M-protein covariate analysis but not retained (Methods S1); median 186 U/L (Table 1)."
    ),
    DIS_EMD = list(
      description = "Extramedullary disease at screening",
      units = "(binary)",
      type = "binary",
      notes = "Screened in the M-protein covariate analysis but not retained (Methods S1); present in 9% of the BVd analysis population (Table S3)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 260L,
    n_studies = 2L,
    age_range = "32-86 years",
    age_median = "65 years",
    weight_range = "43.7-133 kg",
    weight_median = "74.7 kg",
    sex_female_pct = 45,
    disease_state = paste(
      "Relapsed/refractory multiple myeloma (>= 1 prior line of therapy)",
      "followed for response by serum M-protein, with baseline serum",
      "M-protein >= 5 g/L.",
      sep = " "
    ),
    dose_range = paste(
      "Belantamab mafodotin IV 1.9, 2.5 or 3.4 mg/kg (single, split, Q6W",
      "stretch or step-down schedules in DREAMM-6 Arm B; 2.5 mg/kg Q3W in",
      "DREAMM-7), with bortezomib and dexamethasone. Initial dose 2.5 mg/kg",
      "in 88% (Table S3).",
      sep = " "
    ),
    co_medication = "Bortezomib + dexamethasone (BVd) in all patients.",
    biomarker_summary = paste(
      "Median (range) at baseline (Table 1, N = 269): albumin 40.0 (23.0-51.0)",
      "g/L; sBCMA 40.2 (2.20-729) ug/L; IgG 20.3 (1.05-81.5) g/L;",
      "beta-2-microglobulin 263 (125-1530) nmol/L; LDH 186 (59.0-2060) U/L;",
      "M-protein 19.0 (5.00-85.0) g/L. Myeloma type IgG 75%, IgA 24%;",
      "R-ISS I/II/III 44/50/6% (Table S3).",
      sep = " "
    ),
    notes = paste(
      "4829 serum M-protein concentrations from 260 patients: 60 from",
      "DREAMM-6 Arm B (NCT03544281) and 200 from DREAMM-7 (NCT04246047)",
      "(Results; Table S2). Table 1 summarises 269 patients, including 9",
      "DREAMM-6 patients with a baseline M-protein value only.",
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
    # Serum M-protein TGI, BVd model (Table S4 'BVd model' column; the
    # parameterisation follows File S2, which implements this model).
    # KGR, KDEATH and LAM are reported on the natural-log scale. E50
    # (the slope S of the kill term) is also a log-scale THETA: File S2
    # reads 'LNS = THETA69; // LNE50 in run82; S = exp(LNS + ETAS)'.
    # A50 is linear (File S2: MU_6 = log(TVAL)).
    # ------------------------------------------------------------------
    lkg <- -4.7; label("Log M-protein growth rate constant KGR (log 1/day)") # Table S4 BVd: TVKGR = -4.7 (2.59 %RSE); File S2 THETA63 = -4.7
    lkd <- -3.06; label("Log M-protein kill rate constant KDEATH at IGG = exp(2.74) g/L (log 1/day)") # Table S4 BVd: TVKDEATH = -3.06 (2.33 %RSE); File S2 THETA64 = -3.06
    llambda <- -5.12; label("Log resistance rate constant LAM (log 1/day)") # Table S4 BVd: TVLAM = -5.12 (1.68 %RSE); File S2 THETA65 = -5.12
    lslope <- 0.398; label("Log slope S of the kill term 1 - exp(-S * C) (log mL/ug)") # Table S4 BVd: TVE50 = 0.398 (7.05 %RSE); File S2 THETA69 = 0.398 (LNS)
    la50 <- log(1.07); label("M-protein concentration halving the growth rate, A50 (g/L)") # Table S4 BVd: TVA50 = 1.07 g/L (2.38 %RSE); File S2 THETA68 = 1.07
    ltsmax <- fixed(log(150)); label("Logistic carrying capacity of serum M-protein (g/L)") # File S2 $DES: KGR * (1 - A(4)/150); not estimated
    e_igg_kd <- -0.0884; label("Coefficient of log(IGG) - 2.74 on log KDEATH (unitless)") # Table S4 BVd: IGGBL on KDEATH = -0.0884 (43 %RSE); File S2 THETA71, KDEATH_IGGBL = THETA71 * (log(IGGBL) - 2.74)

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

    # M-protein IIV (Table S4 BVd). The 3x3 block is listed in NONMEM
    # lower-triangle order; the fourth entry is printed as a second
    # 'Covariance (KGR~KDEATH)' but sits in the KGR~LAM position.
    etalkg + etalkd + etallambda ~ c(
      1.65,
      0.76, 1.03,
      0.191, 0.694, 1.16
    ) # Table S4 BVd: omega KGR 1.65, cov KGR-KDEATH 0.76, omega KDEATH 1.03, cov KGR-LAM 0.191, cov KDEATH-LAM 0.694, omega LAM 1.16
    etalslope ~ fixed(0.122) # Table S4 BVd: omega E50 = 0.122, held constant
    etala50 ~ fixed(0.122) # Table S4 BVd: omega A50 = 0.122, held constant

    # Residual error. Table S4 heads the proportional row 'variance', but
    # File S2 carries it as a THETA beside the additive SD THETA, in the
    # W = sqrt(PROP^2 * IPRED^2 + ADD^2) form of its commented WBAS line,
    # so 0.333 is used as the proportional SD.
    propSd <- 0.2516; label("PK proportional residual error (fraction)") # carried from Papathanasiou_2025_belantamab (Papathanasiou 2025 Table 2: sqrt(0.0633))
    propSd_mprotein <- 0.333; label("M-protein proportional residual error (fraction)") # Table S4 BVd: residual proportional error = 0.333; File S2 THETA67 = 0.333 (PROP ERR)
    addSd_mprotein <- fixed(0.0001); label("M-protein additive residual error (g/L)") # Table S4 BVd: residual error FIXED, additive SD = 0.0001 g/L; File S2 THETA66
  })

  model({
    # 1. PK driver (File S2 $PK). Imax takes its covariates before the
    #    additive eta; COMBO_BELAMAF = 1 for BVd.
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

    # 2. M-protein parameters (File S2 TUMOR MODEL block). The IgG effect
    #    is additive on log KDEATH, centred on log(IgG) = 2.74.
    kg <- exp(lkg + etalkg)
    kd <- exp(lkd + e_igg_kd * (log(IGG) - 2.74) + etalkd)
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
