Moein_2024_apitolisib_human <- function() {
  description <- paste(
    "Integrated PK-PD-efficacy model for orally administered apitolisib",
    "(GDC-0980), a dual PI3K/mTOR inhibitor, in patients with advanced solid",
    "tumors or non-Hodgkin's lymphoma (two phase 1 dose-escalation studies,",
    "NCT00854152 and NCT00854126). Three sequentially-fitted layers in one",
    "ODE system: (1) two-compartment first-order-absorption population PK",
    "with an absorption lag (individual estimates carried in and held fixed);",
    "(2) an indirect-response model in which apitolisib inhibits the",
    "production rate of phosphorylated Akt (pAkt) measured in PLATELET-RICH",
    "PLASMA as a surrogate for tumor target modulation, expressed as %pAkt",
    "relative to a drug-free baseline of 100%; (3) exponential growth of the",
    "RECIST sum of longest diameters at net rate kg opposed by a sigmoidal",
    "shrinkage rate ks driven by the percent pAkt inhibition",
    "I = 100 - %pAkt. Imax and the pAkt Hill exponent are fixed at 1, and",
    "kout is not estimated but derived as kin/%pAkt(0) = kin/100, so the",
    "biomarker pool equilibrates essentially instantaneously and %pAkt",
    "tracks the plasma concentration directly. Only the five efficacy",
    "parameters (kg, kmax, ki50, gamma2, proportional error) and the four",
    "tumor IIV terms were estimated in this integrated fit.",
    "Companion preclinical model: modellib('Moein_2024_apitolisib_mouse')."
  )
  reference <- paste(
    "Moein A, Jin JY, Wright MR, Alicke B, Wong H. Retrospective Assessment",
    "of Translational Pharmacokinetic-Pharmacodynamic Modeling Performance:",
    "A Case Study with Apitolisib, a Dual PI3K/mTOR Inhibitor.",
    "Drugs R D. 2024;24(2):157-166. doi:10.1007/s40268-024-00459-5.",
    "PMCID PMC11315854. Structural equations from Methods Eqs. 2-5;",
    "parameter values from Results Sect. 3.1.2 and Tables 2 and 4; NONMEM",
    "control stream from Supplementary Information Online Resource 2",
    "(part II). Underlying phase 1 clinical data reported by Dolly SO et al.",
    "Clin Cancer Res. 2016;22(12):2874-2884."
  )
  vignette <- "Moein_2024_apitolisib"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list()

  compartmentData <- list(
    depot       = list(analyte = "apitolisib", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "apitolisib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "apitolisib", units = "mg", specimen = "plasma", verified = TRUE),
    pakt        = list(analyte = "phosphorylated Akt (serine 473) in platelet-rich plasma, percent of baseline", units = "%", specimen = "plasma", verified = TRUE),
    tumor_size  = list(analyte = "tumour size (RECIST sum of longest diameters)", units = "mm", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species       = "human",
    n_subjects    = 117L,
    n_studies     = 2L,
    disease_state = "advanced solid tumors or non-Hodgkin's lymphoma",
    dose_range    = paste(
      "Apitolisib PO: a single dose on day 1 followed by a 1-week washout,",
      "then 30, 40 or 50 mg once daily for 28 days; or 2, 4, 8, 16, 32, 40,",
      "50 or 70 mg once daily on a 3-weeks-on/1-week-off schedule; or 6, 12,",
      "25, 50, 100, 150 or 200 mg once weekly for 28 days."
    ),
    regions       = "(not reported in the modelling paper)",
    notes         = paste(
      "Two phase 1 dose-escalation studies using a 3+3 design",
      "(ClinicalTrials.gov NCT00854152 and NCT00854126). Sub-populations",
      "differ by layer: the population PK model used n = 146 patients with",
      "1835 PK observations (Results Sect. 3.1.2); the pAkt PK-PD model used",
      "the n = 37 patients in whom platelet-rich-plasma pAkt was measured",
      "(Methods Sect. 2.2); and the integrated PK-PD-efficacy model used the",
      "n = 117 evaluable patients contributing 417 tumor-size observations",
      "(Results Sect. 3.3). Patients with only a baseline tumor assessment",
      "were excluded. Median baseline tumor observation was 130.0 mm",
      "(Results Sect. 3.3).",
      "Tumor size is the sum of the longest diameters of target lesions per",
      "RECIST 1.0, or modified RECIST for malignant pleural mesothelioma,",
      "measured by CT or MRI at baseline, after the first and second cycle,",
      "and every two cycles thereafter (Methods Sect. 2.2).",
      "pAkt was collected pre-dose and 1, 3, 8 and 24 h post-dose on day 1",
      "and quantified by Meso Scale Discovery (MSD) assay; %pAkt was computed",
      "relative to baseline (Methods Eq. 2), so baseline is 100% by",
      "definition. Peak platelet pAkt suppression was 90% at doses >= 12 mg",
      "(Results Sect. 3.2). Apitolisib plasma concentrations were assayed by",
      "LC-MS/MS with an LLOQ of 0.5 ng/mL.",
      "In the integrated fit, patients WITH pAkt data used their individual",
      "IC50 estimate and patients WITHOUT pAkt data used the typical IC50",
      "(Methods Sect. 2.5). Per-tumor-type individual IC50 means ranged from",
      "4.0 to 26.1 ug/L (Supplementary Table S1); the mean for the two",
      "metastatic renal cell carcinoma patients was 16.2 ug/L, i.e. at the",
      "high end and closest to the 403 ug/L tumor-tissue estimate in the",
      "companion 786-O RCC xenograft model.",
      "Body-weight-normalized PK at a 70 kg reference weight: CL/F 0.304",
      "L/h/kg, V1/F 3.01 L/kg, V2/F 9.13 L/kg, CLd/F 0.085 L/h/kg (Results",
      "Sect. 3.1.2). Apitolisib fraction unbound in human plasma was 38.8%,",
      "versus 29.4% in mouse (Results Sect. 3.2). No covariates were retained",
      "in any layer of either model."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Layer 1: patient population PK -- two compartments, first-order oral
    # absorption with an absorption lag (Methods Sect. 2.3.2, Results
    # Sect. 3.1.2). Fitted separately; the integrated model then reads each
    # patient's INDIVIDUAL PK parameters in from the data set (Online
    # Resource 2 part II $INPUT: ICL, IV1, IV2, IQ, IKA, IALAG1). Encoded
    # fixed() at the reported typical values because they were not
    # re-estimated in this integrated fit.
    #
    # Online Resource 2 part II works on a WEEK time base and rescales
    # (CL = ICL*168, KA = IKA*168, ALAG1 = IALAG1/168). This file uses an
    # HOUR time base instead, which keeps every PK and pAkt value verbatim
    # as published and needs a conversion on only the two tumor rate
    # constants below.
    # ------------------------------------------------------------------
    lka   <- fixed(log(3.35));   label("Absorption rate constant ka (1/h)")                          # Results Sect. 3.1.2: ka = 3.35 1/h (RSE 16.8%)
    lcl   <- fixed(log(21.3));   label("Apparent oral clearance CL/F (L/h)")                         # Results Sect. 3.1.2: CL/F = 21.3 L/h (RSE 4.28%)
    lvc   <- fixed(log(210.6));  label("Apparent central volume of distribution V1/F (L)")           # Results Sect. 3.1.2: V1/F = 210.6 L (RSE 3.51%)
    lvp   <- fixed(log(639.1));  label("Apparent peripheral volume of distribution V2/F (L)")        # Results Sect. 3.1.2: V2/F = 639.1 L (RSE 15.1%)
    lq    <- fixed(log(5.93));   label("Apparent intercompartmental clearance CLd/F (L/h)")          # Results Sect. 3.1.2: CLd/F = 5.93 L/h (RSE 8.53%)
    ltlag <- fixed(log(0.465));  label("Absorption lag time (h)")                                   # Results Sect. 3.1.2: lagged time = 0.465 h (RSE 0.43%)

    # ------------------------------------------------------------------
    # Layer 2: pAkt indirect-response model (Methods Eq. 3, Table 2).
    # Fitted to platelet-rich-plasma %pAkt in the n = 37 patients with
    # biomarker data, then held fixed in the integrated model.
    # ------------------------------------------------------------------
    lic50      <- fixed(log(9.32));   label("Apitolisib concentration giving 50% of Imax on pAkt, IC50 (ug/L)")  # Table 2: IC50 = 9.32 ug/L (RSE 12.7%)
    lkin       <- fixed(log(88699));  label("pAkt zero-order production rate kin (%/h)")                         # Table 2: kin = 88,699 %/h (RSE 3.8%)
    imax       <- fixed(1);           label("Maximum fractional inhibition of pAkt production, Imax (unitless)") # Table 2: Imax = 1 (fix), to allow 100% biomarker inhibition
    lhill_pakt <- fixed(log(1));      label("Hill exponent on the apitolisib-pAkt inhibition, gamma1 (unitless)") # Table 2: gamma1 = 1 (fix); footnote ** "estimate was 0.868-0.937; therefore it was fixed to 1 in the final model"
    # kout is NOT a free parameter: Table 2 footnote * and Online Resource 1
    # Sect. II state kout = kin / %pAkt(0) with %pAkt(0) = 100%, giving the
    # tabulated 886.99 1/h. It is derived in model() rather than declared here.

    # ------------------------------------------------------------------
    # Layer 3: integrated PK-PD-efficacy tumor model (Methods Eqs. 4-5,
    # Table 4). These five thetas plus the four omegas below are the ONLY
    # quantities estimated in this integrated fit.
    #   d(Tumor)/dt = kg*Tumor - ks*Tumor
    #   ks = kmax * I^gamma2 / (ki50^gamma2 + I^gamma2),  I = 100 - %pAkt
    #
    # Table 4 reports kg and kmax per WEEK; the / 168 converts to the 1/h
    # time base of this file (168 h per week, the same factor Online
    # Resource 2 part II uses in the opposite direction).
    # ------------------------------------------------------------------
    lkg      <- log(0.0097 / 168); label("Net tumor growth rate constant kg (1/h; = 0.0097 1/week)")                  # Table 4: Kg = 0.0097 1/week (RSE 20.7%)
    lkmax    <- log(0.0142 / 168); label("Maximum tumor shrinkage rate constant kmax (1/h; = 0.0142 1/week)")         # Table 4: Kmax = 0.0142 1/week (RSE 43.9%)
    lki50    <- log(58.0);         label("Percent pAkt inhibition giving 50% of kmax, ki50 (% pAkt inhibition)")      # Table 4: KI50 = 58.0 %pAkt inhibition (RSE 34.3%)
    lhill_ks <- log(6.52);         label("Sigmoidicity factor of the pAkt-inhibition/shrinkage curve, gamma2 (unitless)") # Table 4: gamma2 = 6.52 (RSE 101.4%)

    # Baseline tumor size. The authors used each patient's OWN observed
    # baseline as the initial condition (Online Resource 2 part II $INPUT
    # "TBSL ; Individual Tumor baseline"; Results Sect. 3.3 "Individual
    # baseline tumor data were used for model development"), so this is not
    # an estimated theta. Encoded fixed() at the reported median so the
    # model has a usable typical-value default for simulation; override per
    # subject via a tumor_size initial-condition event.
    lrbase_tumor <- fixed(log(130.0)); label("Median baseline tumor size, RECIST SLD (mm) -- data-derived, not estimated")  # Results Sect. 3.3: median baseline tumor observation = 130.0 mm

    # ------------------------------------------------------------------
    # Inter-individual variability.
    # All values are VARIANCES on the exponential (log-normal) scale:
    # Methods Sect. 2.3.2 and Online Resource 1 Sect. II state "IIV in PK
    # parameters were modeled as exponential random-effect models ... thus
    # assumed to follow a log-normal distribution". Table 4's
    # "Variability%" column is the corresponding CV%
    # (sqrt(0.668) = 81.7%, sqrt(0.377) = 61.4%), which confirms variances
    # rather than SDs.
    #
    # PK and IC50 IIV are carried in from the separately-fitted upstream
    # layers and are therefore fixed(); only the four tumor IIV terms were
    # part of this integrated fit. No IIV correlations are reported.
    # ------------------------------------------------------------------
    etalcl   ~ fixed(0.17);    # Results Sect. 3.1.2: IIV CL/F variance = 0.17 (RSE 17.6%, shrinkage 13.7%)
    etalvc   ~ fixed(0.0639);  # Results Sect. 3.1.2: IIV V1/F variance = 0.0639 (RSE 23.2%, shrinkage 23.2%)
    etalvp   ~ fixed(0.881);   # Results Sect. 3.1.2: IIV V2/F variance = 0.881 (RSE 44.8%, shrinkage 40.5%)
    etalq    ~ fixed(0.208);   # Results Sect. 3.1.2: IIV CLd/F variance = 0.208 (RSE 33.5%, shrinkage 44.2%)
    etalka   ~ fixed(2.09);    # Results Sect. 3.1.2: IIV ka variance = 2.09 (RSE 16.4%, shrinkage 13.2%)
    etalic50 ~ fixed(0.296);   # Table 2: IIV IC50 variance = 0.296 (RSE 44.3%, shrinkage 15.5%)

    etalkg      ~ 0.668;          # Table 4: IIV Kg variance = 0.668 (RSE 34.9%, Variability 81.7%, shrinkage 37.9%)
    etalkmax    ~ 0.377;          # Table 4: IIV Kmax variance = 0.377 (RSE 44.6%, Variability 61.4%, shrinkage 55.8%)
    etalki50    ~ fixed(0.0225);  # Table 4: IIV KI50 variance = 0.0225 (fix) -- Online Resource 1 Sect. II: 'fixed to a small value (0.0225) to enhance the numerical estimation capability'
    etalhill_ks ~ fixed(0.0225);  # Table 4: IIV gamma2 variance = 0.0225 (fix) -- same rationale

    # ------------------------------------------------------------------
    # Residual error.
    # ------------------------------------------------------------------
    # PK: Methods Sect. 2.3.2 states "The residual error model was a
    # proportional error model" but no magnitude is reported in the paper,
    # the tables, or Online Resource 1 or 2. Encoded fixed(0) per the
    # unreported-RUV convention; see the vignette Errata.
    propSd            <- fixed(0);     label("Proportional residual error on plasma apitolisib (fraction) -- not reported")
    # pAkt: estimated in the stage-2 PK-PD model (Table 2) and not
    # re-estimated in the integrated fit, hence fixed() here.
    propSd_pakt       <- fixed(0.672); label("Proportional residual error on %pAkt (fraction)")           # Table 2: proportional error = 0.672 (RSE 9.5%)
    # Tumor size: estimated in THIS integrated fit.
    propSd_tumor_size <- 0.141;        label("Proportional residual error on tumor size (fraction)")      # Table 4: proportional error = 0.141 (RSE 2.3%)
  })

  model({
    # ----- 1. Individual parameters -------------------------------------
    ka           <- exp(lka + etalka)
    cl           <- exp(lcl + etalcl)
    vc           <- exp(lvc + etalvc)
    vp           <- exp(lvp + etalvp)
    q            <- exp(lq + etalq)
    tlag         <- exp(ltlag)
    ic50         <- exp(lic50 + etalic50)
    kin          <- exp(lkin)
    hill_pakt    <- exp(lhill_pakt)
    kg           <- exp(lkg + etalkg)
    kmax         <- exp(lkmax + etalkmax)
    ki50         <- exp(lki50 + etalki50)
    hill_ks      <- exp(lhill_ks + etalhill_ks)
    rbase_tumor  <- exp(lrbase_tumor)

    # kout is derived, not estimated: Table 2 footnote gives
    # kout = kin / %pAkt(0) with %pAkt(0) = 100% by the definition of
    # %pAkt (Methods Eq. 2). 88699 / 100 = 886.99 1/h, as tabulated.
    kout <- kin / 100

    # ----- 2. PK: two compartments, first-order oral absorption + lag ---
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    alag(depot)       <- tlag
    # Dose enters as mg and vc is L, so central/vc is mg/L; the * 1000
    # converts to the ug/L (= ng/mL) scale of IC50 (Table 2) and of the
    # 0.5 ng/mL assay LLOQ, declared in units$concentration. Online
    # Resource 2 part II instead declares AMT in ug so that central/V is
    # already ug/L -- the same quantity.
    Cc <- central / vc * 1000                                       # ug/L

    # ----- 3. pAkt indirect response (Methods Eq. 3) --------------------
    # Inhibition of the pAkt PRODUCTION rate. With imax and hill_pakt both
    # fixed at 1 this collapses to Cc / (ic50 + Cc).
    inh_pakt <- imax * Cc^hill_pakt / (ic50^hill_pakt + Cc^hill_pakt)
    d/dt(pakt) <- kin * (1 - inh_pakt) - kout * pakt
    # %pAkt at the drug-free baseline is 100% by definition (Methods
    # Eq. 2); Online Resource 2 part II sets A_0(4) = 100.
    pakt(0) <- 100

    # ----- 4. Tumor growth (Methods Eqs. 4-5) ---------------------------
    # I = 100 - %pAkt is the percent inhibition of Akt phosphorylation.
    # Online Resource 2 part II clamps it at zero ("R = 100 - A(4);
    # IF (R.LT.0) R = 0") so a transient overshoot of %pAkt above baseline
    # cannot produce a negative base for the non-integer power below.
    inhib_pakt_pct <- max(0, 100 - pakt)
    ks <- kmax * inhib_pakt_pct^hill_ks / (ki50^hill_ks + inhib_pakt_pct^hill_ks)
    d/dt(tumor_size) <- kg * tumor_size - ks * tumor_size
    tumor_size(0)    <- rbase_tumor                                 # mm (RECIST SLD)

    # ----- 5. Observations (three independent endpoints) ----------------
    # Downstream fitting with a dvid column must order the endpoints as
    # declared here (1 = Cc, 2 = pakt, 3 = tumor_size).
    Cc         ~ prop(propSd)
    pakt       ~ prop(propSd_pakt)
    tumor_size ~ prop(propSd_tumor_size)
  })
}
