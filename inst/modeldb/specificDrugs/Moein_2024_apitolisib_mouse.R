Moein_2024_apitolisib_mouse <- function() {
  description <- paste(
    "Preclinical (mouse; 786-O human renal cell adenocarcinoma xenograft).",
    "Integrated PK-PD-efficacy model for orally administered apitolisib",
    "(GDC-0980), a dual PI3K/mTOR inhibitor. Three sequentially-fitted layers",
    "in one ODE system: (1) one-compartment first-order-absorption PK",
    "(CL/F, V/F, ka carried in from the mouse PK model and held fixed);",
    "(2) an indirect-response model in which apitolisib inhibits the",
    "production rate of phosphorylated Akt (pAkt) measured in TUMOR TISSUE,",
    "expressed as %pAkt relative to a drug-free control baseline of 100%;",
    "(3) exponential tumor-volume growth at net rate p opposed by a",
    "sigmoidal shrinkage rate ks driven by the percent pAkt inhibition",
    "I = 100 - %pAkt. Imax and the pAkt Hill exponent are fixed at 1, and",
    "kout is not estimated but derived as kin/%pAkt(0) = kin/100, so the",
    "biomarker pool equilibrates essentially instantaneously and %pAkt",
    "tracks the plasma concentration directly. Only the five efficacy",
    "parameters (p, kmax, ki50, gamma2, proportional error) and the two",
    "tumor IIV terms were estimated in this integrated fit."
  )
  reference <- paste(
    "Moein A, Jin JY, Wright MR, Alicke B, Wong H. Retrospective Assessment",
    "of Translational Pharmacokinetic-Pharmacodynamic Modeling Performance:",
    "A Case Study with Apitolisib, a Dual PI3K/mTOR Inhibitor.",
    "Drugs R D. 2024;24(2):157-166. doi:10.1007/s40268-024-00459-5.",
    "PMCID PMC11315854. Structural equations from Methods Eqs. 2-5;",
    "parameter values from Tables 1 and 3; NONMEM control stream from",
    "Supplementary Information Online Resource 2 (part I).",
    "Companion human model: modellib('Moein_2024_apitolisib_human')."
  )
  vignette <- "Moein_2024_apitolisib"
  units <- list(time = "h", dosing = "mg/kg", concentration = "ug/L")

  covariateData <- list()

  compartmentData <- list(
    depot     = list(analyte = "apitolisib", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central   = list(analyte = "apitolisib", units = "mg/kg", specimen = "plasma", verified = TRUE),
    pakt      = list(analyte = "phosphorylated Akt (serine 473), percent of drug-free control", units = "%", specimen = "tumor", verified = TRUE),
    tumor_vol = list(analyte = "tumour size", units = "mm^3", specimen = "tumor", verified = TRUE)
  )

  population <- list(
    species        = "mouse (female beige nude XID, bg.nu.xid Harlan; subcutaneous 786-O human renal cell adenocarcinoma xenograft)",
    n_subjects     = 64L,
    n_studies      = 3L,
    sex_female_pct = 100,
    disease_state  = "subcutaneous 786-O human renal cell adenocarcinoma (RCC) xenograft",
    dose_range     = paste(
      "Efficacy study: vehicle plus apitolisib 0.008, 0.026, 0.085, 0.256, 1,",
      "2.5, 5, 7, 8.5, 10 and 11 mg/kg PO once daily for 17 days (n=6 per",
      "group). pAkt study: single oral dose of vehicle, 0.3, 3 or 10 mg/kg.",
      "PK study: single oral doses of 1, 5 and 10 mg/kg."
    ),
    regions        = "Preclinical (Genentech, Inc., South San Francisco, CA, USA)",
    notes          = paste(
      "The integrated PK-PD-efficacy model was fitted to n = 64 xenograft mice",
      "contributing 381 tumor-volume observations (Results Sect. 3.3). Median",
      "baseline tumor volume was 173.5 mm^3; the mean across all groups was",
      "181 mm^3 at the initiation of dosing (Methods Sect. 2.1 and Results",
      "Sect. 3.3). Tumor volumes were measured with digital calipers as",
      "TV (mm^3) = length * width^2 * 0.5 (Methods Eq. 1) pre-dose and at 72,",
      "96, 168, 240, 264, 336 and 408 h. Mice were euthanized if they lost",
      ">20% of initial body weight or if tumors exceeded 2000 mm^3.",
      "The pAkt sub-study used separate animals: tumors were collected at 0.5,",
      "2, 6, 12, 18 and 24 h post-dose (n=3 per time point; n=6 per time point",
      "for the vehicle group) and pAkt (serine 473) was quantified by Meso",
      "Scale Discovery (MSD) assay in TUMOR TISSUE -- unlike the companion",
      "human model, whose biomarker is measured in platelet-rich plasma.",
      "%pAkt was computed relative to the drug-free control (Methods Eq. 2),",
      "so the control baseline is 100% by definition. The pAkt PD analysis was",
      "a naive average pooled fit: individual %pAkt values were averaged per",
      "time point within each dose level (Online Resource 1, Sect. II).",
      "Apitolisib fraction unbound in mouse plasma was 29.4% (Results",
      "Sect. 3.2), versus 38.8% in human.",
      "The PK data (n=27 mice per dose group at 1, 5 and 10 mg/kg) were",
      "previously published and re-used here (Methods Sect. 2.1)."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Layer 1: mouse PK -- one compartment, first-order oral absorption.
    # Fitted separately to mean plasma concentration-time data pooled over
    # the 1, 5 and 10 mg/kg single-dose groups (Methods Sect. 2.3.1,
    # Results Sect. 3.1.1), then HARD-CODED as constants in the integrated
    # model's $PK block (Online Resource 2 part I: "CL = 0.387", "V = 1.13",
    # "KA = 1.1"). Encoded fixed() because they were not re-estimated in
    # this integrated fit.
    # ------------------------------------------------------------------
    lka <- fixed(log(1.1));    label("Absorption rate constant ka (1/h)")                     # Results Sect. 3.1.1: ka = 1.1 1/h (RSE 36.4%)
    lcl <- fixed(log(0.387));  label("Apparent oral clearance CL/F (L/h/kg)")                 # Results Sect. 3.1.1: CL/F = 0.387 L/h/kg (RSE 9.5%)
    lvc <- fixed(log(1.13));   label("Apparent oral volume of distribution V/F (L/kg)")       # Results Sect. 3.1.1: V/F = 1.13 L/kg (RSE 17.3%)

    # ------------------------------------------------------------------
    # Layer 2: pAkt indirect-response model (Methods Eq. 3, Table 1).
    # d(%pAkt)/dt = kin * (1 - Imax*Cp^g1/(IC50^g1 + Cp^g1)) - kout * %pAkt
    # Fitted to tumor-tissue %pAkt from the 0.3 / 3 / 10 mg/kg single-dose
    # sub-study, then held fixed in the integrated model (Online Resource 2
    # part I hard-codes IC50, KIN, KOUT, HILL and IMAX in $PK).
    # ------------------------------------------------------------------
    lic50      <- fixed(log(403));   label("Apitolisib concentration giving 50% of Imax on pAkt, IC50 (ug/L)")  # Table 1: IC50 = 403 ug/L (RSE 14.7%)
    lkin       <- fixed(log(88698)); label("pAkt zero-order production rate kin (%/h)")                         # Table 1: kin = 88,698 %/h (RSE 0.8%)
    imax       <- fixed(1);          label("Maximum fractional inhibition of pAkt production, Imax (unitless)") # Table 1: Imax = 1 (fix), to allow 100% biomarker inhibition
    lhill_pakt <- fixed(log(1));     label("Hill exponent on the apitolisib-pAkt inhibition, gamma1 (unitless)") # Table 1: gamma1 = 1 (fix); earlier iterations estimated ~1 (Online Resource 1 Sect. II)
    # kout is NOT a free parameter: Table 1 footnote and Online Resource 1
    # Sect. II state kout = kin / %pAkt(0) with %pAkt(0) = 100%, giving the
    # tabulated 886.98 1/h. It is derived in model() rather than declared here.

    # ------------------------------------------------------------------
    # Layer 3: integrated PK-PD-efficacy tumor model (Methods Eqs. 4-5,
    # Table 3). These five thetas plus the two omegas below are the ONLY
    # quantities estimated in this integrated fit.
    #   d(Tumor)/dt = p*Tumor - ks*Tumor
    #   ks = kmax * I^gamma2 / (ki50^gamma2 + I^gamma2),  I = 100 - %pAkt
    # ------------------------------------------------------------------
    lp       <- log(0.0057); label("Net tumor growth rate constant p (Kg in the paper; 1/h)")                                       # Table 3: Kg = 0.0057 1/h (RSE 3.9%)
    lkmax     <- log(0.0194); label("Maximum tumor shrinkage rate constant kmax (1/h)")                              # Table 3: Kmax = 0.0194 1/h (RSE 43.2%)
    lki50     <- log(67.2);   label("Percent pAkt inhibition giving 50% of kmax, ki50 (% pAkt inhibition)")          # Table 3: KI50 = 67.2 %pAkt inhibition (RSE 12.6%)
    lhill_ks  <- log(9.4);    label("Sigmoidicity factor of the pAkt-inhibition/shrinkage curve, gamma2 (unitless)") # Table 3: gamma2 = 9.4 (RSE 208.5%)

    # Baseline tumor volume. The authors used each mouse's OWN observed
    # baseline as the initial condition (Online Resource 2 part I:
    # "A_0(4) = TBSL ; Tumor Baseline"; Results Sect. 3.3 "Individual
    # baseline tumor data were used for model development"), so this is not
    # an estimated theta. Encoded fixed() at the reported median so the
    # model has a usable typical-value default for simulation; override per
    # subject via a tumor_vol initial-condition event.
    lrbase_tumor <- fixed(log(173.5)); label("Median baseline tumor volume (mm^3) -- data-derived, not estimated")  # Results Sect. 3.3: median baseline tumor observation = 173.5 mm^3

    # ------------------------------------------------------------------
    # Inter-individual variability (Table 3). Reported as VARIANCES on the
    # exponential (log-normal) scale -- Online Resource 1 Sect. II: "IIV of
    # parameters were modeled as exponential random-effect models ... thus
    # assumed to follow a log-normal distribution". The paper's
    # "Variability%" column is the corresponding CV%
    # (sqrt(0.0276) = 16.6%, sqrt(0.0529) = 23.0%), which confirms the
    # values are variances and not SDs.
    # ------------------------------------------------------------------
    etalp   ~ 0.0276  # Table 3: IIV Kg variance = 0.0276 (RSE 29.6%, Variability 16.6%, shrinkage 22.6%)
    etalkmax ~ 0.0529  # Table 3: IIV Kmax variance = 0.0529 (RSE 54.6%, Variability 23.0%, shrinkage 45.8%)

    # ------------------------------------------------------------------
    # Residual error.
    # ------------------------------------------------------------------
    # PK: the mouse PK model was fitted to MEAN concentration-time profiles
    # and reported "data not shown" (Results Sect. 3.1.1); no residual-error
    # magnitude is reported anywhere on disk. Encoded fixed(0) per the
    # unreported-RUV convention; see the vignette Errata.
    propSd           <- fixed(0);     label("Proportional residual error on plasma apitolisib (fraction) -- not reported")
    # pAkt: estimated in the stage-2 PK-PD model (Table 1) and not
    # re-estimated in the integrated fit, hence fixed() here.
    propSd_pakt      <- fixed(0.353); label("Proportional residual error on %pAkt (fraction)")            # Table 1: proportional error = 0.353 (RSE 5.6%)
    # Tumor volume: estimated in THIS integrated fit.
    propSd_tumor_vol <- 0.248;        label("Proportional residual error on tumor volume (fraction)")     # Table 3: proportional error = 0.248 (RSE 5.1%)
  })

  model({
    # ----- 1. Individual parameters -------------------------------------
    ka           <- exp(lka)
    cl           <- exp(lcl)
    vc           <- exp(lvc)
    ic50         <- exp(lic50)
    kin          <- exp(lkin)
    hill_pakt    <- exp(lhill_pakt)
    p            <- exp(lp + etalp)
    kmax         <- exp(lkmax + etalkmax)
    ki50         <- exp(lki50)
    hill_ks      <- exp(lhill_ks)
    rbase_tumor  <- exp(lrbase_tumor)

    # kout is derived, not estimated: Table 1 footnote gives
    # kout = kin / %pAkt(0) with %pAkt(0) = 100% by the definition of
    # %pAkt (Methods Eq. 2). 88698 / 100 = 886.98 1/h, as tabulated.
    kout <- kin / 100

    # ----- 2. PK: one compartment, first-order oral absorption ----------
    kel <- cl / vc
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central
    # Dose enters as mg/kg and vc is L/kg, so central/vc is mg/L; the
    # * 1000 converts to the ug/L scale of IC50 (Table 1) declared in
    # units$concentration. Online Resource 2 part I works in mg/L instead
    # and writes IC50 = 0.403 mg/L -- the same quantity.
    Cc <- central / vc * 1000                                       # ug/L

    # ----- 3. pAkt indirect response (Methods Eq. 3) --------------------
    # Inhibition of the pAkt PRODUCTION rate. With imax and hill_pakt both
    # fixed at 1 this collapses to Cc / (ic50 + Cc).
    inh_pakt <- imax * Cc^hill_pakt / (ic50^hill_pakt + Cc^hill_pakt)
    d/dt(pakt) <- kin * (1 - inh_pakt) - kout * pakt
    # %pAkt at the drug-free control baseline is 100% by definition
    # (Methods Eq. 2); Online Resource 2 part I sets A_0(3) = 100.
    pakt(0) <- 100

    # ----- 4. Tumor growth (Methods Eqs. 4-5) ---------------------------
    # I = 100 - %pAkt is the percent inhibition of Akt phosphorylation.
    # Online Resource 2 part I clamps it at zero ("R = 100 - A(3);
    # IF (R.LT.0) R = 0") so a transient overshoot of %pAkt above baseline
    # cannot produce a negative base for the non-integer power below.
    inhib_pakt_pct <- max(0, 100 - pakt)
    ks <- kmax * inhib_pakt_pct^hill_ks / (ki50^hill_ks + inhib_pakt_pct^hill_ks)
    d/dt(tumor_vol) <- p * tumor_vol - ks * tumor_vol
    tumor_vol(0)    <- rbase_tumor                                  # mm^3

    # ----- 5. Observations (three independent endpoints) ----------------
    # Downstream fitting with a dvid column must order the endpoints as
    # declared here (1 = Cc, 2 = pakt, 3 = tumor_vol).
    Cc        ~ prop(propSd)
    pakt      ~ prop(propSd_pakt)
    tumor_vol ~ prop(propSd_tumor_vol)
  })
}
