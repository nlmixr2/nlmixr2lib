Zhang_2025_nedosiran_uoxcr <- function() {
  description <- "Population PK/PD model of the spot urine oxalate-to-creatinine ratio (Uox/Cr) for nedosiran in patients with primary hyperoxaluria type 1: the two-compartment dual-parallel-transit-absorption PK model with parallel linear and Michaelis-Menten elimination drives an effect compartment, whose concentration inhibits the zero-order production of Uox/Cr in an indirect-response turnover model through a sigmoidal Imax function, with age on the Uox/Cr baseline, developed from 668 spot Uox/Cr observations in 41 patients across three trials (Zhang 2025)"
  reference <- "Zhang S, Gamallo P, Rawson V. Population Pharmacokinetic and Pharmacodynamic Modelling and Simulation for Nedosiran Clinical Development and Dose Guidance in Pediatric Patients with Primary Hyperoxaluria Type 1. Clin Pharmacokinet. 2025;64(8):1213-1227. doi:10.1007/s40262-025-01540-1"
  vignette <- "Zhang_2025_nedosiran"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  paper_specific_compartments <- c("uoxcr")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric size descriptor for the PK layer, reference 70 kg (Zhang 2025 Table 2 covariate-formula footnote). Fixed exponents 0.750 on CL/F and Q/F and 1.00 on Vc/F and Vp/F; estimated exponents -0.221 shared by ka1 and ka2 and 0.492 on Vmax. PK/PD cohort range 11.9-115.9 kg (Supplementary Table S2).",
      source_name        = "BW"
    ),
    CRCL = list(
      description        = "BSA-normalized estimated glomerular filtration rate",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 90 mL/min/1.73 m^2 (Zhang 2025 Table 2 covariate-formula footnote). Power effects on CL/F (0.969) and Vc/F (0.174) in the PK layer. Estimating equation varies by cohort per Supplementary Table S3. PK/PD cohort range 35-197 mL/min/1.73 m^2 (Supplementary Table S2).",
      source_name        = "eGFR"
    ),
    DIS_PH1 = list(
      description        = "Primary hyperoxaluria type 1 disease-state indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = primary hyperoxaluria type 2 patient or healthy volunteer",
      notes              = "1 for PH1 subjects, 0 for PH2 subjects and healthy volunteers (Zhang 2025 Table 2 covariate-formula footnote). Enters as the multiplicative factor 1.32^DIS_PH1 on the slow absorption rate ka1 in the PK layer. Every subject contributing spot Uox/Cr data to this model had PH1, so DIS_PH1 = 1 throughout the PD analysis population.",
      source_name        = "PH"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 6 years (Zhang 2025 Table 3 covariate-formula footnote). Power effect -0.450 on the Uox/Cr baseline, so younger patients have higher baseline spot Uox/Cr. PK/PD cohort range 2-46 years (Supplementary Table S2).",
      source_name        = "AGE"
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
    peripheral1 = list(analyte = "nedosiran", units = "mg", specimen = "plasma", verified = TRUE),
    effect      = list(analyte = "nedosiran", units = "ng/mL", specimen = "not applicable", verified = TRUE),
    uoxcr       = list(analyte = "oxalate-to-creatinine ratio", units = "mmol/mol", specimen = "urine", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 41,
    n_studies      = 3,
    n_observations = 668,
    age_range      = "2-46 years",
    weight_range   = "11.9-115.9 kg",
    sex_female_pct = 46.3,
    race_ethnicity = c(White = 41.5, Multiple = 22.0, Asian = 19.5, Unknown = 14.6, `Black or African American` = 2.4),
    disease_state  = "primary hyperoxaluria type 1",
    renal_function = "normal 36.6%, mild impairment 41.5%, moderate impairment 22.0%; eGFR 35-197 mL/min/1.73 m^2",
    dose_range     = "136 or 170 mg SC once-monthly (age >= 12 years) and 3.5 mg/kg SC once-monthly capped at 136 or 170 mg (age < 12 years)",
    regions        = "multinational",
    notes          = "Patients with primary hyperoxaluria type 1 and at least one post-baseline spot Uox/Cr value from PHYOX2 (NCT03847909), PHYOX3 (NCT04042402) and PHYOX8 (NCT05001269). Age subgroups: 2 to <6 y (n=8), 6 to <9 y (n=5), 9 to <12 y (n=5), 12 to <18 y (n=6), >=18 y (n=17). The PK layer of this model was estimated on the larger 148-participant PK dataset; see modellib('Zhang_2025_nedosiran'). Baseline demographics are Supplementary Table S2 of Zhang 2025."
  )

  ini({
    # ---------------------------------------------------------------
    # PK layer (Zhang 2025 Table 2) -- identical to the standalone
    # Pop-PK model Zhang_2025_nedosiran. Reference subject: BW 70 kg,
    # eGFR 90 mL/min/1.73 m^2, PH type other than PH1.
    # ---------------------------------------------------------------
    lcl   <- log(6.10);   label("Apparent clearance (CL/F, L/h)")                                  # Table 2: CL/F = 6.10 L/h (RSE 11.4%, 95% CI 4.74-7.47)
    lvc   <- log(148);    label("Apparent central volume of distribution (Vc/F, L)")               # Table 2: Vc/F = 148 L (RSE 6.07%, 95% CI 130-165)
    lq    <- log(2.79);   label("Apparent inter-compartmental clearance (Q/F, L/h)")               # Table 2: Q/F = 2.79 L/h (RSE 13.0%, 95% CI 2.08-3.51)
    lvp   <- log(6560);   label("Apparent peripheral volume of distribution (Vp/F, L)")            # Table 2: Vp/F = 6560 L (RSE 22.4%, 95% CI 3690-9440)
    lvmax <- log(3.37);   label("Maximum rate of saturable elimination (Vmax, mg/h)")              # Table 2: Vmax = 3.37 mg/h (RSE 20.9%, 95% CI 1.99-4.75)
    lkm   <- log(248);    label("Michaelis-Menten constant for saturable elimination (KM, ng/mL)") # Table 2: KM = 248 ng/mL (RSE 26.7%, 95% CI 118-378)

    lka1     <- log(0.212);   label("First-order absorption/transit rate constant, slow pathway (ka1, 1/h)") # Table 2: ka1 = 0.212 1/h (RSE 6.78%, 95% CI 0.184-0.240)
    lka2     <- log(14.9);    label("First-order absorption/transit rate constant, fast pathway (ka2, 1/h)") # Table 2: ka2 = 14.9 1/h (RSE 4.48%, 95% CI 13.6-16.2)
    logitffo <- 0.8094862;    label("Logit of the fraction of the dose absorbed via the slow pathway (FR1, logit units)") # Table 2: FR1 = 0.692 (RSE 2.20%, 95% CI 0.662-0.721); qlogis(0.692) = 0.8094862

    e_wt_cl_q     <- fixed(0.750); label("Allometric exponent on CL/F and Q/F (unitless)")             # Table 2: CL.BW = 0.750, shared by CL/F and Q/F per the covariate-formula footnote
    e_wt_vc_vp    <- fixed(1.00);  label("Allometric exponent on Vc/F and Vp/F (unitless)")            # Table 2: V.BW = 1.00, shared by Vc/F and Vp/F per the covariate-formula footnote
    e_wt_ka       <- -0.221;       label("Body-weight exponent shared by ka1 and ka2 (unitless)")      # Table 2: ka.BW = -0.221 (RSE 47.9%, 95% CI -0.428 to 0.0135)
    e_wt_vmax     <- 0.492;        label("Body-weight exponent on Vmax (unitless)")                    # Table 2: Vmax.BW = 0.492 (RSE 14.2%, 95% CI 0.355-0.629)
    e_crcl_cl     <- 0.969;        label("eGFR exponent on CL/F (unitless)")                           # Table 2: CL.EGFR = 0.969 (RSE 11.6%, 95% CI 0.749-1.19)
    e_crcl_vc     <- 0.174;        label("eGFR exponent on Vc/F (unitless)")                           # Table 2: Vc.EGFR = 0.174 (RSE 18.1%, 95% CI 0.112-0.236)
    e_ph1_ka1 <- 1.32;         label("Multiplicative factor on ka1 for primary hyperoxaluria type 1 (unitless)") # Table 2: ka1.PH = 1.32 (RSE 9.66%, 95% CI 1.07-1.57); enters as ka1.PH^PH

    etalvc       ~ 0.0749600   # Table 2: Vc/F.IIV = 27.9 CV% (shrinkage 23.7%); log(1 + 0.279^2)
    etalka1      ~ 0.2175672   # Table 2: ka1.IIV  = 49.3 CV% (shrinkage 16.8%); log(1 + 0.493^2)
    etalka2      ~ 0.2508801   # Table 2: ka2.IIV  = 53.4 CV% (shrinkage 8.94%); log(1 + 0.534^2)
    etalvmax     ~ 0.1161826   # Table 2: Vmax.IIV = 35.1 CV% (shrinkage 11.5%); log(1 + 0.351^2)
    etalogitffo  ~ 0.116964    # Table 2: FR1.IIV = 0.342 SD, additive on the logit scale (shrinkage 35.4%); 0.342^2

    expSd <- 0.3010484; label("Exponential residual error SD on the log scale for nedosiran (log ng/mL)") # Table 2: ExpError = 30.8 CV% (shrinkage 12.6%); sqrt(log(1 + 0.308^2))

    # ---------------------------------------------------------------
    # PD layer (Zhang 2025 Table 3). Reference age 6 years. kout, the
    # Hill exponent and the effect-compartment equilibration half-life
    # were fixed to the estimates of the earlier 24-h urinary oxalate
    # PD model of Zhang et al. because the spot Uox/Cr dataset carries
    # no off-treatment observations (Zhang 2025 Sect. 3.1 and Sect. 4).
    #
    # The paper reports the PD rate constants in weeks; the model time
    # unit is hours, so each is divided by 168 h/week.
    # ---------------------------------------------------------------
    lrbase <- log(264);   label("Baseline spot urine oxalate-to-creatinine ratio at the reference age (BSL, mmol/mol)") # Table 3: BSL = 264 mmol/mol (RSE 11.9%, 95% CI 202-325)
    limax  <- log(0.687); label("Maximum inhibitory effect on Uox/Cr production (Imax, fraction)")                      # Table 3: Imax = 0.687 (RSE 3.50%, 95% CI 0.640-0.734), i.e. 68.7%
    lic50  <- log(1.68);  label("Effect-compartment concentration producing half of Imax (IC50, ng/mL)")                # Table 3: IC50 = 1.68 ng/mL (RSE 21.0%, 95% CI 0.992-2.38)
    lhill  <- fixed(log(2.56));         label("Hill exponent of the sigmoidal Imax function (unitless)")                # Table 3: Gamma = 2.56
    lkout  <- fixed(-6.208673);         label("First-order elimination rate constant of the Uox/Cr pool (kout, 1/h)")   # Table 3: kout = 0.338 1/wk; log(0.338 / 168) = -6.208673
    lke0   <- fixed(-8.576964);         label("Effect-compartment equilibration rate constant (ke0, 1/h)")              # Table 3: Lambda = 21.9 wk equilibration half-life; ke0 = log(2)/Lambda per Fig. 1; log(log(2) / (21.9 * 168)) = -8.576964

    e_age_rbase <- -0.450; label("Age exponent on the Uox/Cr baseline (unitless)") # Table 3: AGE.BSL = -0.450 (RSE 14.6%, 95% CI -0.578 to -0.321)

    etalrbase ~ 0.1638894  # Table 3: BSL.IIV  = 42.2 CV% (shrinkage 4.42%); log(1 + 0.422^2)
    etalic50  ~ 0.5457510  # Table 3: IC50.IIV = 85.2 CV% (shrinkage 34.1%); log(1 + 0.852^2)

    propSd_uoxcr <- 0.355; label("Proportional residual error on spot Uox/Cr (fraction)") # Table 3: PropError = 35.5 CV% (shrinkage 3.42%)
  })

  model({
    # ==================================================================
    # 1. PK layer -- individual parameters. Covariate formulas are the
    #    footnote to Zhang 2025 Table 2; see Zhang_2025_nedosiran for
    #    the full transcription, including the note that the printed
    #    Vmax formula omits its Vmax.BW superscript.
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
    # 2. PD layer -- individual parameters. The covariate formula is the
    #    footnote to Zhang 2025 Table 3:
    #      BSL = TVBSL * (AGE/6)^AGE.BSL
    #    kin follows from the turnover steady state at baseline in the
    #    absence of drug: kin = kout * BSL.
    # ==================================================================
    imax  <- exp(limax)
    ic50  <- exp(lic50 + etalic50)
    hill  <- exp(lhill)
    kout  <- exp(lkout)
    ke0   <- exp(lke0)
    rbase <- exp(lrbase + etalrbase) * (AGE / 6)^e_age_rbase
    kin   <- kout * rbase

    # ==================================================================
    # 3. Plasma concentration. Doses are in mg and vc in L, so
    #    central / vc is mg/L = ug/mL; multiply by 1000 for ng/mL, the
    #    scale on which KM and IC50 were reported.
    # ==================================================================
    Cc <- 1000 * central / vc

    # ==================================================================
    # 4. Drug effect. Zhang 2025 Fig. 1 prints
    #      Eff = Imax * Ceff^gamma / (Ceff^gamma + IC50)
    #    but the IC50 in the denominator has lost its gamma superscript:
    #    the printed form is dimensionally inconsistent and would not
    #    make IC50 the half-maximal concentration that Table 3 labels it
    #    as. The canonical sigmoidal Imax form is used here. See the
    #    vignette's Assumptions and deviations section.
    # ==================================================================
    eff <- imax * effect^hill / (effect^hill + ic50^hill)

    # ==================================================================
    # 5. ODE system (Zhang 2025 Fig. 1 schematic).
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
    d/dt(effect)      <-  ke0 * (Cc - effect)
    d/dt(uoxcr)       <-  kin * (1 - eff) - kout * uoxcr

    # ==================================================================
    # 6. Initial condition. Treatment starts from the untreated
    #    turnover steady state.
    # ==================================================================
    uoxcr(0) <- rbase

    # ==================================================================
    # 7. Dose splitting. Each subcutaneous injection is entered by the
    #    USER as TWO simultaneous dose records, both carrying the full
    #    injected amount (cmt = "depot" and cmt = "depot2"); the f()
    #    multipliers below perform the FR1 / (1 - FR1) split.
    # ==================================================================
    f(depot)  <- ffo
    f(depot2) <- 1 - ffo

    # ==================================================================
    # 8. Observations and residual error.
    # ==================================================================
    Cc ~ lnorm(expSd)
    uoxcr ~ prop(propSd_uoxcr)
  })
}
