Zhang_2025_nedosiran <- function() {
  description <- "Two-compartment population PK model for nedosiran (GalNAc-conjugated anti-LDHA siRNA) after subcutaneous dosing, with dual parallel transit absorption (a slow one-transit pathway carrying fraction FR1 of the dose and a fast three-transit pathway carrying the remainder), parallel linear and saturable Michaelis-Menten elimination, allometric body-weight scaling and eGFR effects on CL/F and Vc/F, and a primary-hyperoxaluria-type-1 effect on the slow absorption rate, developed from 2087 plasma concentrations in 148 healthy volunteers and patients with primary hyperoxaluria across six trials (Zhang 2025)"
  reference <- "Zhang S, Gamallo P, Rawson V. Population Pharmacokinetic and Pharmacodynamic Modelling and Simulation for Nedosiran Clinical Development and Dose Guidance in Pediatric Patients with Primary Hyperoxaluria Type 1. Clin Pharmacokinet. 2025;64(8):1213-1227. doi:10.1007/s40262-025-01540-1"
  vignette <- "Zhang_2025_nedosiran"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric size descriptor, reference 70 kg (Zhang 2025 Table 2 covariate-formula footnote). Fixed exponents 0.750 on CL/F and Q/F and 1.00 on Vc/F and Vp/F; estimated exponents -0.221 shared by ka1 and ka2 and 0.492 on Vmax. Cohort range 11.1-129.5 kg (Supplementary Table S1).",
      source_name        = "BW"
    ),
    CRCL = list(
      description        = "BSA-normalized estimated glomerular filtration rate",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 90 mL/min/1.73 m^2 (Zhang 2025 Table 2 covariate-formula footnote). Power effects on CL/F (0.969) and Vc/F (0.174). Estimating equation varies by cohort per Supplementary Table S3: CKD-EPI (adults, PHYOX2/3/5/6), MDRD (adults, PHYOX1), bedside Schwartz 2009 (children, PHYOX1/2/3), Schwartz 2012 multivariate (children, PHYOX8), Matsuo (Japanese adults, PHYOX3) and Uemura (Japanese children, PHYOX3/8). For paediatric participants a renal maturation function RMF = PMA^3.4 / (PMA^3.4 + 47.6^3.4), PMA in weeks, was applied when deriving eGFR. Cohort range 3.8-197 mL/min/1.73 m^2 (Supplementary Table S1).",
      source_name        = "eGFR"
    ),
    DIS_PH1 = list(
      description        = "Primary hyperoxaluria type 1 disease-state indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = primary hyperoxaluria type 2 patient or healthy volunteer",
      notes              = "1 for PH1 subjects, 0 for PH2 subjects and healthy volunteers (Zhang 2025 Table 2 covariate-formula footnote). Enters as the multiplicative factor 1.32^DIS_PH1 on the slow absorption rate ka1. Cohort composition: 57.4% healthy volunteers, 33.1% PH1, 9.5% PH2 (Supplementary Table S1).",
      source_name        = "PH"
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "nedosiran", units = "mg", specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "nedosiran", units = "mg", specimen = "administration site", verified = TRUE),
    depot2      = list(analyte = "nedosiran", units = "mg", specimen = "administration site", verified = TRUE),
    transit2    = list(analyte = "nedosiran", units = "mg", specimen = "administration site", verified = TRUE),
    transit3    = list(analyte = "nedosiran", units = "mg", specimen = "administration site", verified = TRUE),
    transit4    = list(analyte = "nedosiran", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "nedosiran", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "nedosiran", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 148,
    n_studies      = 6,
    n_observations = 2087,
    age_range      = "1.9-73 years",
    weight_range   = "11.1-129.5 kg",
    sex_female_pct = 35.8,
    race_ethnicity = c(White = 56.1, Asian = 19.6, `Black or African American` = 11.5, Multiple = 6.1, Other = 4.7, Unknown = 2.0),
    disease_state  = "healthy volunteers (57.4%) and patients with primary hyperoxaluria type 1 (33.1%) or type 2 (9.5%)",
    renal_function = "normal 41.9%, mild impairment 35.1%, moderate impairment 6.8%, severe impairment 4.7%, end-stage renal disease 11.5%; eGFR 3.8-197 mL/min/1.73 m^2",
    dose_range     = "0.3-12 mg/kg single SC (PHYOX1), 1.5-6 mg/kg single SC (PHYOX5, PHYOX6), 170 mg single SC (PHYOX5), 136 or 170 mg SC once-monthly and 3.5 mg/kg SC once-monthly capped at 136 or 170 mg (PHYOX2, PHYOX3, PHYOX8)",
    regions        = "multinational; PHYOX6 enrolled Japanese and Caucasian adults",
    notes          = "Pooled from PHYOX1 (NCT03392896), PHYOX2 (NCT03847909), PHYOX3 (NCT04042402), PHYOX5, PHYOX6 and PHYOX8 (NCT05001269). Age subgroups: 0 to <2 y (n=1), 2 to <6 y (n=12), 6 to <9 y (n=6), 9 to <12 y (n=5), 12 to <18 y (n=10), >=18 y (n=114). Baseline demographics are Supplementary Table S1 of Zhang 2025."
  )

  ini({
    # ---------------------------------------------------------------
    # Structural disposition (Zhang 2025 Table 2), reference subject:
    # BW 70 kg, eGFR 90 mL/min/1.73 m^2, PH type other than PH1.
    # ---------------------------------------------------------------
    lcl   <- log(6.10);   label("Apparent clearance (CL/F, L/h)")                                  # Table 2: CL/F = 6.10 L/h (RSE 11.4%, 95% CI 4.74-7.47)
    lvc   <- log(148);    label("Apparent central volume of distribution (Vc/F, L)")               # Table 2: Vc/F = 148 L (RSE 6.07%, 95% CI 130-165)
    lq    <- log(2.79);   label("Apparent inter-compartmental clearance (Q/F, L/h)")               # Table 2: Q/F = 2.79 L/h (RSE 13.0%, 95% CI 2.08-3.51)
    lvp   <- log(6560);   label("Apparent peripheral volume of distribution (Vp/F, L)")            # Table 2: Vp/F = 6560 L (RSE 22.4%, 95% CI 3690-9440)
    lvmax <- log(3.37);   label("Maximum rate of saturable elimination (Vmax, mg/h)")              # Table 2: Vmax = 3.37 mg/h (RSE 20.9%, 95% CI 1.99-4.75)
    lkm   <- log(248);    label("Michaelis-Menten constant for saturable elimination (KM, ng/mL)") # Table 2: KM = 248 ng/mL (RSE 26.7%, 95% CI 118-378)

    # ---------------------------------------------------------------
    # Dual parallel transit absorption (Zhang 2025 Fig. 1 schematic).
    # ka1 drives the slow pathway (one transit compartment), ka2 the
    # fast pathway (three transit compartments); FR1 is the fraction of
    # the dose entering the slow pathway.
    # ---------------------------------------------------------------
    lka1     <- log(0.212);   label("First-order absorption/transit rate constant, slow pathway (ka1, 1/h)") # Table 2: ka1 = 0.212 1/h (RSE 6.78%, 95% CI 0.184-0.240)
    lka2     <- log(14.9);    label("First-order absorption/transit rate constant, fast pathway (ka2, 1/h)") # Table 2: ka2 = 14.9 1/h (RSE 4.48%, 95% CI 13.6-16.2)
    logitffo <- 0.8094862;    label("Logit of the fraction of the dose absorbed via the slow pathway (FR1, logit units)") # Table 2: FR1 = 0.692 (RSE 2.20%, 95% CI 0.662-0.721); qlogis(0.692) = 0.8094862

    # ---------------------------------------------------------------
    # Covariate effects (Zhang 2025 Table 2 and its covariate-formula
    # footnote). Reference values BW = 70 kg, eGFR = 90 mL/min/1.73 m^2.
    # ---------------------------------------------------------------
    e_wt_cl_q     <- fixed(0.750); label("Allometric exponent on CL/F and Q/F (unitless)")             # Table 2: CL.BW = 0.750, shared by CL/F and Q/F per the covariate-formula footnote
    e_wt_vc_vp    <- fixed(1.00);  label("Allometric exponent on Vc/F and Vp/F (unitless)")            # Table 2: V.BW = 1.00, shared by Vc/F and Vp/F per the covariate-formula footnote
    e_wt_ka       <- -0.221;       label("Body-weight exponent shared by ka1 and ka2 (unitless)")      # Table 2: ka.BW = -0.221 (RSE 47.9%, 95% CI -0.428 to 0.0135)
    e_wt_vmax     <- 0.492;        label("Body-weight exponent on Vmax (unitless)")                    # Table 2: Vmax.BW = 0.492 (RSE 14.2%, 95% CI 0.355-0.629)
    e_crcl_cl     <- 0.969;        label("eGFR exponent on CL/F (unitless)")                           # Table 2: CL.EGFR = 0.969 (RSE 11.6%, 95% CI 0.749-1.19)
    e_crcl_vc     <- 0.174;        label("eGFR exponent on Vc/F (unitless)")                           # Table 2: Vc.EGFR = 0.174 (RSE 18.1%, 95% CI 0.112-0.236)
    e_ph1_ka1 <- 1.32;         label("Multiplicative factor on ka1 for primary hyperoxaluria type 1 (unitless)") # Table 2: ka1.PH = 1.32 (RSE 9.66%, 95% CI 1.07-1.57); enters as ka1.PH^PH

    # ---------------------------------------------------------------
    # Inter-individual variability (Zhang 2025 Table 2). The paper
    # reports log-normal IIV as CV%, so omega^2 = log(1 + CV^2); the FR1
    # IIV is logit-normal and reported directly as an SD.
    # ---------------------------------------------------------------
    etalvc       ~ 0.0749600   # Table 2: Vc/F.IIV = 27.9 CV% (shrinkage 23.7%); log(1 + 0.279^2)
    etalka1      ~ 0.2175672   # Table 2: ka1.IIV  = 49.3 CV% (shrinkage 16.8%); log(1 + 0.493^2)
    etalka2      ~ 0.2508801   # Table 2: ka2.IIV  = 53.4 CV% (shrinkage 8.94%); log(1 + 0.534^2)
    etalvmax     ~ 0.1161826   # Table 2: Vmax.IIV = 35.1 CV% (shrinkage 11.5%); log(1 + 0.351^2)
    etalogitffo  ~ 0.116964    # Table 2: FR1.IIV = 0.342 SD, additive on the logit scale (shrinkage 35.4%); 0.342^2

    # ---------------------------------------------------------------
    # Residual unexplained variability (Zhang 2025 Table 2).
    # ---------------------------------------------------------------
    expSd <- 0.3010484; label("Exponential residual error SD on the log scale (log ng/mL)") # Table 2: ExpError = 30.8 CV% (shrinkage 12.6%); sqrt(log(1 + 0.308^2))
  })

  model({
    # ==================================================================
    # 1. Individual parameters. Covariate formulas are the footnote to
    #    Zhang 2025 Table 2:
    #      CL/F  = TVCL   * (BW/70)^CL.BW  * (eGFR/90)^CL.eGFR
    #      Vc/F  = TVVc   * (BW/70)^V.BW   * (eGFR/90)^Vc.eGFR
    #      Q/F   = TVQ    * (BW/70)^CL.BW
    #      Vp/F  = TVVp   * (BW/70)^V.BW
    #      ka1   = TVka1  * (BW/70)^ka.BW  * (ka1.PH)^PH
    #      ka2   = TVka2  * (BW/70)^ka.BW
    #      Vmax  = TVVmax * (BW/70)^Vmax.BW
    #    The printed Vmax formula omits its superscript; Table 2 reports
    #    Vmax.BW = 0.492 as an estimated exponent with an RSE, so the
    #    power form is the only self-consistent reading. See the
    #    vignette's Assumptions and deviations section.
    # ==================================================================
    cl   <- exp(lcl)                 * (WT / 70)^e_wt_cl_q  * (CRCL / 90)^e_crcl_cl
    vc   <- exp(lvc + etalvc)        * (WT / 70)^e_wt_vc_vp * (CRCL / 90)^e_crcl_vc
    q    <- exp(lq)                  * (WT / 70)^e_wt_cl_q
    vp   <- exp(lvp)                 * (WT / 70)^e_wt_vc_vp
    ka1  <- exp(lka1 + etalka1)      * (WT / 70)^e_wt_ka * e_ph1_ka1^DIS_PH1
    ka2  <- exp(lka2 + etalka2)      * (WT / 70)^e_wt_ka
    vmax <- exp(lvmax + etalvmax)    * (WT / 70)^e_wt_vmax
    km   <- exp(lkm)
    ffo  <- expit(logitffo + etalogitffo)

    # ==================================================================
    # 2. Plasma concentration. Doses are in mg and vc in L, so
    #    central / vc is mg/L = ug/mL; multiply by 1000 for ng/mL, which
    #    is the scale on which KM was reported.
    # ==================================================================
    Cc <- 1000 * central / vc

    # ==================================================================
    # 3. ODE system (Zhang 2025 Fig. 1 schematic).
    #      depot    -> transit1                     -> central  (rate ka1)
    #      depot2   -> transit2 -> transit3 -> transit4 -> central (rate ka2)
    #    Elimination from central is the sum of a linear CL/F term and a
    #    saturable Michaelis-Menten term Vmax * C / (C + KM).
    # ==================================================================
    d/dt(depot)       <- -ka1 * depot
    d/dt(transit1)    <-  ka1 * depot - ka1 * transit1
    d/dt(depot2)      <- -ka2 * depot2
    d/dt(transit2)    <-  ka2 * depot2 - ka2 * transit2
    d/dt(transit3)    <-  ka2 * transit2 - ka2 * transit3
    d/dt(transit4)    <-  ka2 * transit3 - ka2 * transit4
    d/dt(central)     <-  ka1 * transit1 + ka2 * transit4 -
      (cl / vc) * central -
      (q / vc) * central + (q / vp) * peripheral1 -
      vmax * Cc / (km + Cc)
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1

    # ==================================================================
    # 4. Dose splitting. Each subcutaneous injection is entered by the
    #    USER as TWO simultaneous dose records, both carrying the full
    #    injected amount (cmt = "depot" and cmt = "depot2"); the f()
    #    multipliers below perform the FR1 / (1 - FR1) split. See the
    #    vignette for a worked event-table recipe.
    # ==================================================================
    f(depot)  <- ffo
    f(depot2) <- 1 - ffo

    # ==================================================================
    # 5. Observation and residual error.
    # ==================================================================
    Cc ~ lnorm(expSd)
  })
}
