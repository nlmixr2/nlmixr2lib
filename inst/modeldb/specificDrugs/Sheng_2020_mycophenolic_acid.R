Sheng_2020_mycophenolic_acid <- function() {
  description <- "Population PK model of unbound and total mycophenolic acid (uMPA, tMPA) and its 7-O-glucuronide metabolite (MPAG) in Chinese adult kidney transplant recipients co-treated with cyclosporine and given oral mycophenolate mofetil (MMF) (Sheng 2020). Five compartments (Figure 1): a gut depot with first-order absorption and a lag time, a two-compartment disposition of UNBOUND MPA (central + peripheral1), a one-compartment disposition of unbound MPAG (central_mpag), and a gallbladder (gallbladder_mpag) giving intermittent enterohepatic circulation. 87 percent of uMPA elimination forms MPAG (fixed); the rest is eliminated directly. Total MPA is linked to unbound MPA by a linear protein-binding constant, tMPA = uMPA * (1 + kB), with kB proportional to serum albumin; total MPAG = unbound MPAG / 0.18 (fixed unbound fraction). The gallbladder fills from central_mpag at a rate set by the estimated fraction of MPAG recycled and empties into the gut (where MPAG is assumed fully deconjugated to MPA and reabsorbed) at a fixed first-order rate during a 0.5 h window after each meal. Meal gates are read against time after the most recent dose (tad()) at 4 and 10 h, the study-1 schedule. Covariates: body weight on Q/F of uMPA, serum albumin on kB, glomerular filtration rate (CKD-EPI, mL/min) on CL/F of uMPAG. Dose in mg MMF; the model works internally in umol and reports concentrations in mg/L."
  reference <- paste(
    "Sheng C, Zhao Q, Niu W, Qiu X, Zhang M, Jiao Z.",
    "Effect of Protein Binding on Exposure of Unbound and Total Mycophenolic",
    "Acid: A Population Pharmacokinetic Analysis in Chinese Adult Kidney",
    "Transplant Recipients.",
    "Front Pharmacol. 2020;11:340.",
    "doi:10.3389/fphar.2020.00340.",
    sep = " "
  )
  vignette <- "Sheng_2020_mycophenolic_acid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(
      analyte = "mycophenolate mofetil (MMF), plus recycled MPA expressed as MMF-mass equivalents",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(analyte = "unbound mycophenolic acid (uMPA)", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(
      analyte = "unbound mycophenolic acid (uMPA)",
      units = "umol",
      specimen = "plasma",
      verified = TRUE
    ),
    central_mpag = list(
      analyte = "unbound 7-O-mycophenolic acid glucuronide (uMPAG)",
      units = "umol",
      specimen = "plasma",
      verified = TRUE
    ),
    gallbladder_mpag = list(
      analyte = "7-O-mycophenolic acid glucuronide (MPAG)",
      units = "umol",
      specimen = "bile",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight (kg).",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the apparent inter-compartmental clearance of uMPA, Q/F = 857 * (WT/70)^2.11 (Table 2 footnote c). Cohort range 40-82.5 kg (Table 1).",
      source_name = "BW"
    ),
    ALB = list(
      description = "Serum albumin (g/L).",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear (power exponent held at 1) effect on the protein-binding constant, kB = 53.4 * (ALB/40) (Table 2 footnote d; Results: 'the estimated value of the exponent for the effect of ALB on kB was quite close to 1, it was fixed at 1'). Cohort range 20-50 g/L (Table 1).",
      source_name = "ALB"
    ),
    CRCL = list(
      description = "Glomerular filtration rate estimated from serum creatinine by the CKD-EPI equation (Levey 2009), reported by the source in mL/min.",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the apparent clearance of uMPAG, CL/F = 5.71 * (GFR/80)^0.865 (Table 2 footnote e). The source labels the CKD-EPI estimate 'mL/min' (Table 1 footnote b); CKD-EPI natively returns mL/min/1.73 m^2 and the paper does not state that it was de-normalised, so supply the CKD-EPI value as the source did. Cohort range 11.2-123.8 mL/min (Table 1).",
      source_name = "GFR"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator (1 = female).",
      units = "(binary)",
      type = "binary",
      notes = "Tested on CL/F, Q/F and V/F of uMPA and uMPAG; the sex effect on Q/F entered in the forward step but was removed in backward elimination (Results)."
    ),
    AGE = list(
      description = "Age (years).",
      units = "years",
      type = "continuous",
      notes = "Screened (Methods) but not retained."
    ),
    HGB = list(
      description = "Hemoglobin (g/L).",
      units = "g/L",
      type = "continuous",
      notes = "Screened (Methods) but not retained."
    ),
    CONMED_CSA_DOSE = list(
      description = "Cyclosporine daily dose (mg/day).",
      units = "mg/day",
      type = "continuous",
      notes = "Screened (Methods) but not retained; every subject received cyclosporine."
    ),
    CONMED_ANTACID = list(
      description = "Co-administration of antacids (1 = yes).",
      units = "(binary)",
      type = "binary",
      notes = "Tested on CL/F of uMPA and on ka (Results) but not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 58L,
    n_studies = 2L,
    n_profiles = 65L,
    n_observations = "740 uMPA, 741 tMPA and 734 total MPAG concentrations (Abstract; Results).",
    age_range = "18-62 years",
    weight_range = "40-82.5 kg",
    sex_female_pct = 22.4,
    race_ethnicity = "Chinese (100%)",
    disease_state = "Adult first-time kidney transplant recipients on triple immunosuppression (MMF + cyclosporine + corticosteroids).",
    dose_range = "Oral MMF (CellCept) 750-2000 mg/day, most commonly 1000-1500 mg/day given every 12 h (Table 1).",
    regions = "China (Huashan Hospital, Fudan University, Shanghai; study 2 multicentre).",
    postoperative_time = "3-3084 days; study 1 was mostly within 3 months and study 2 mostly beyond 3 months post-transplant (Table 1; Results).",
    albumin_range = "20-50 g/L (Table 1).",
    gfr_range = "11.2-123.8 mL/min by CKD-EPI (Table 1).",
    co_medication = "Cyclosporine in all subjects (100-400 mg/day in most); corticosteroids; antacids in 11 subjects (Table 1).",
    notes = "Pooled from an early-post-transplant PK study (Jiao 2007; 20 patients, 27 full 12 h profiles) and a stable-phase bioequivalence study (Geng 2012; 38 patients, 38 profiles). Meals after the 4 and 10 h samples in study 1 and after the 3 and 9 h samples in study 2. Estimation by SAEM followed by importance sampling in NONMEM 7.4."
  )

  ini({
    # ---- Absorption (Table 2) ----
    lka    <- log(1.35);   label("First-order absorption rate constant (1/h)")        # Table 2 ka = 1.35 /h (RSE 11.1%)
    ltlag  <- log(0.447);  label("Absorption lag time (h)")                           # Table 2 Tlag = 0.447 h (RSE 16.8%)

    # ---- Unbound MPA disposition (Table 2) ----
    lcl    <- log(851);    label("Apparent clearance of unbound MPA (L/h)")           # Table 2 CLuMPA/F = 851 L/h (RSE 7.1%)
    lvc    <- log(718);    label("Apparent central volume of unbound MPA (L)")        # Table 2 VCuMPA/F = 718 L (RSE 18.5%)
    lq     <- log(857);    label("Apparent inter-compartmental clearance of unbound MPA at 70 kg (L/h)")  # Table 2 QuMPA/F = 857 L/h (RSE 11.0%)
    lvp    <- fixed(log(34300)); label("Apparent peripheral volume of unbound MPA (L)")  # Results: VPuMPA/F held at the de Winter 2009 value 34,300 L; Table 3 '34,300 FIXED'
    e_wt_q <- 2.11;        label("Power exponent of body weight on Q/F of unbound MPA (unitless)")  # Table 2 exponent BW on QuMPA/F = 2.11 (RSE 24.2%); footnote c reference 70 kg

    # ---- Metabolic branch (Methods assumption 2; Figure 1) ----
    fm     <- fixed(0.87); label("Fraction of unbound MPA elimination forming MPAG (fraction)")  # Methods: conversion ratio MPA to MPAG held at 87% (FDA CellCept label); Figure 1 k24 87%, k20 13%

    # ---- Linear protein binding of MPA (Methods Eq 1; Table 2) ----
    # tMPA = uMPA + kB * uMPA, so kB is the bound:unbound ratio
    # (dimensionless in Eq 1 although Table 2 prints the unit '/h').
    lkns       <- log(53.4); label("Linear protein-binding constant kB at albumin 40 g/L (unitless)")  # Table 2 kB = 53.4 (RSE 2.3%); footnote d reference ALB 40 g/L
    e_alb_kns  <- fixed(1);  label("Power exponent of serum albumin on kB (unitless)")                 # Results: exponent 'fixed at 1'; footnote d

    # ---- Unbound MPAG disposition (Table 2) ----
    lcl_mpag   <- log(5.71); label("Apparent clearance of unbound MPAG at GFR 80 mL/min (L/h)")  # Table 2 CLuMPAG/F = 5.71 L/h (RSE 4.4%); footnote e reference GFR 80 mL/min
    lvc_mpag   <- log(29.9); label("Apparent central volume of unbound MPAG (L)")               # Table 2 VCuMPAG/F = 29.9 L (RSE 7.7%)
    e_crcl_cl_mpag <- 0.865; label("Power exponent of GFR on CL/F of unbound MPAG (unitless)")   # Table 2 exponent GFR on CLuMPAG/F = 0.865 (RSE 11.6%)
    fu_mpag    <- fixed(0.18); label("Unbound fraction of MPAG (fraction)")                     # Methods: FUMPAG held at 18% (FDA CellCept label); Figure 1

    # ---- Intermittent enterohepatic circulation (Methods Eq 4; Table 2) ----
    # %EHC = kGG / (kGG + ke0) * 100, estimated; kGG (kbm here) is derived.
    lehcp  <- log(0.0553);  label("Fraction of MPAG recycled via the gallbladder (fraction)")  # Table 2 %EHC = 5.53 (RSE 26.2%)
    lkehc  <- fixed(log(3.708)); label("Gallbladder emptying rate constant (1/h)")           # Methods: kGB held at 3.708 /h (Guiastrennec 2016); Figure 1
    dge    <- fixed(0.5);   label("Duration of gallbladder emptying after each meal (h)")      # Methods: DGB held at 0.5 h; Figure 1
    tmeal1 <- fixed(4);     label("First meal time after each dose (h)")                      # Methods: low-fat meals after the 4 and 10 h samples (study 1)
    tmeal2 <- fixed(10);    label("Second meal time after each dose (h)")                     # Methods: low-fat meals after the 4 and 10 h samples (study 1)

    # ---- Between-subject variability ----
    # Table 2 %CV = 100 * sqrt(omega^2): the kB entry '10.0 FIXED' is the
    # variance 0.01 stated in Results, so omega^2 = (CV/100)^2 throughout.
    etalcl   ~ 0.2601        # Table 2 BSV CLuMPA/F = 51.0% -> 0.510^2
    etalq    ~ 0.2070        # Table 2 BSV QuMPA/F = 45.5% -> 0.455^2
    etalvc   ~ 0.6400        # Table 2 BSV VCuMPA/F = 80.0% -> 0.800^2
    etalka   ~ 0.2162        # Table 2 BSV ka = 46.5% -> 0.465^2
    etaltlag ~ 1.1600        # Table 2 BSV Tlag = 107.7% -> 1.077^2
    etalkns  ~ fixed(0.01)   # Table 2 BSV kB = 10.0% held by the authors; Results variance 0.01
    etalcl_mpag + etalvc_mpag ~ c(0.1011, 0.08835, 0.2343)  # Table 2 BSV CLuMPAG/F 31.8%, VCuMPAG/F 48.4%, correlation 0.574 -> cov 0.574*0.318*0.484
    etalehcp ~ 0.3795        # Table 2 BSV %EHC = 61.6% -> 0.616^2

    # ---- Residual error: exponential (Results), log-scale SD = Table 2 %CV / 100 ----
    expSd          <- 0.459; label("Log-scale residual SD of total MPA (unitless)")     # Table 2 RUV tMPA = 45.9%
    expSd_Cunbound <- 0.470; label("Log-scale residual SD of unbound MPA (unitless)")   # Table 2 RUV uMPA = 47.0%
    expSd_mpag     <- 0.220; label("Log-scale residual SD of MPAG (unitless)")          # Table 2 RUV uMPAG = 22.0%
  })

  model({
    # Molecular weights (g/mol) used by the source to convert doses and
    # concentrations to molar units (Methods, ChemIDplus).
    mw_mmf  <- 433.498
    mw_mpa  <- 320.339
    mw_mpag <- 496.462

    # ---- Individual parameters ----
    ka      <- exp(lka + etalka)
    tlag    <- exp(ltlag + etaltlag)
    cl      <- exp(lcl + etalcl)
    vc      <- exp(lvc + etalvc)
    q       <- exp(lq + etalq) * (WT / 70)^e_wt_q
    vp      <- exp(lvp)
    kns     <- exp(lkns + etalkns) * (ALB / 40)^e_alb_kns
    cl_mpag <- exp(lcl_mpag + etalcl_mpag) * (CRCL / 80)^e_crcl_cl_mpag
    vc_mpag <- exp(lvc_mpag + etalvc_mpag)
    ehcp    <- exp(lehcp + etalehcp)
    kehc    <- exp(lkehc)

    # ---- Micro-constants (Figure 1) ----
    kel  <- cl / vc              # CLuMPA/VCuMPA = k24 + k20
    k24  <- fm * kel             # uMPA -> uMPAG (87%)
    k20  <- (1 - fm) * kel       # direct uMPA elimination (13%)
    k23  <- q / vc
    k32  <- q / vp
    ke0  <- cl_mpag / vc_mpag    # uMPAG elimination
    kbm  <- ke0 * ehcp / (1 - ehcp)  # kGG from %EHC = kGG / (kGG + ke0) (Eq 4)

    # ---- Meal-triggered gallbladder emptying ----
    # The gate opens for dge hours at tmeal1 and tmeal2 after every dose.
    tdose   <- tad()
    meal    <- (tdose >= tmeal1) * (tdose < tmeal1 + dge) +
      (tdose >= tmeal2) * (tdose < tmeal2 + dge)
    kgb_act <- kehc * meal

    # ---- ODEs ----
    # depot holds mg MMF; every other state holds umol. MMF is completely
    # hydrolysed to MPA (1:1 molar), MPA -> MPAG is 1:1 molar, and biliary
    # MPAG released to the gut is completely deconjugated to MPA and
    # reabsorbed (Methods assumptions 1-3), so the gallbladder release is
    # re-expressed as MMF-mass equivalents to share the depot with the dose.
    d/dt(depot)            <- -ka * depot + kgb_act * gallbladder_mpag * mw_mmf / 1000
    d/dt(central)          <- ka * depot * 1000 / mw_mmf - (k20 + k24 + k23) * central + k32 * peripheral1
    d/dt(peripheral1)      <- k23 * central - k32 * peripheral1
    d/dt(central_mpag)     <- k24 * central - (ke0 + kbm) * central_mpag
    d/dt(gallbladder_mpag) <- kbm * central_mpag - kgb_act * gallbladder_mpag
    alag(depot) <- tlag

    # ---- Observations (mg/L) ----
    Cunbound      <- central / vc * mw_mpa / 1000                 # uMPA
    Cc            <- Cunbound * (1 + kns)                         # tMPA, Eq 1
    fu            <- 1 / (1 + kns)                                # FUMPA, Eq 2
    Cunbound_mpag <- central_mpag / vc_mpag * mw_mpag / 1000      # uMPAG
    Cc_mpag       <- Cunbound_mpag / fu_mpag                      # total MPAG

    Cc       ~ lnorm(expSd)
    Cunbound ~ lnorm(expSd_Cunbound)
    Cc_mpag  ~ lnorm(expSd_mpag)
  })
}
