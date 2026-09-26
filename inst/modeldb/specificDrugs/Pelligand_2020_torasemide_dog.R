Pelligand_2020_torasemide_dog <- function() {
  description <- "Preclinical (beagle dog). Population PK/PD model for oral torasemide in healthy male Beagle dogs: two-compartment PK with first-order absorption (F fixed at 0.98) and a urinary-excretion arm carrying 61% of total clearance. The daily amount of torasemide excreted in urine drives (i) daily diuresis through a direct power model on top of a baseline urine output and (ii) natriuresis through an indirect-response model whose zero-order input is stimulated by a sigmoid Emax function; a threshold-activated reversible 'diuretic resistance' state reduces that Emax (Pelligand 2020)"
  reference <- paste(
    "Pelligand L, Guillot E, Geneteau A, Guyonnet J, Magnier R, Elliott J,",
    "Peyrou M, Jacobs M. Population Pharmacokinetics and Pharmacodynamics",
    "Modeling of Torasemide and Furosemide After Oral Repeated Administration in",
    "Healthy Dogs. Front Vet Sci. 2020;7:151. doi:10.3389/fvets.2020.00151.",
    sep = " "
  )
  vignette <- "Pelligand_2020_torasemide_dog"
  # The PK parameters are published in absolute units for dogs of about 10 kg
  # (CL in L/h, volumes in L; Results 'Pharmacokinetic Model' rescales CL to
  # 7.7 mL/kg/h 'when scaled by an average body weight of 10 kg'), so doses are
  # absolute mg (mg/kg dose x body weight). central/vc is mg/L; the x1000 in Cc
  # reports ug/L, the unit of Figures 2-3. Urinary torasemide is reported in ug
  # (Figure 4, Table 2 additive error in ug, Table 4 EC50 in ug).
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list()

  # Paper-mechanistic PD states named after the paper's own symbols (Equations
  # 2, 6 and 7): E_Na (natriuresis rate), R_Na (resistance fraction) and Q_Na
  # (sodium amount voided in the current 24-h collection).
  paper_specific_compartments <- c("ena", "rna", "qna")

  compartmentData <- list(
    depot = list(analyte = "torasemide", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "torasemide", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "torasemide", units = "mg", specimen = "plasma", verified = TRUE),
    urine = list(analyte = "torasemide", units = "mg", specimen = "urine", verified = TRUE),
    ena = list(analyte = "sodium", units = "mEq/h", specimen = "not applicable", verified = TRUE),
    rna = list(
      analyte = "diuretic resistance (fraction)",
      units = "(fraction)",
      specimen = "not applicable",
      verified = TRUE
    ),
    qna = list(analyte = "sodium", units = "mEq", specimen = "urine", verified = TRUE)
  )

  population <- list(
    species = "beagle dog",
    n_subjects = 17L,
    n_studies = 2L,
    age_range = "1-2.1 years (study 1); 15-19 months (study 2)",
    weight_range = "9.3-11.6 kg (study 1, mean 10.1 kg); 9.9-11.2 kg (study 2, mean 10.6 kg)",
    sex_female_pct = 0,
    disease_state = "Healthy (clinical examination, haematology and biochemistry)",
    dose_range = paste(
      "Oral torasemide tablets once daily in the morning, ~30 min after food:",
      "study 1 0.1, 0.2, 0.4 (and 0.8) mg/kg/day for 14 days; study 2 0.1, 0.2,",
      "0.3 and 0.4 mg/kg/day on Day 1 and Days 5-14. The 0.8 mg/kg period was",
      "excluded from the PK/PD analysis (renal tubular findings at that dose in a",
      "safety study).",
      sep = " "
    ),
    regions = "France",
    notes = paste(
      "Two placebo-controlled crossover studies (Table 1 and Table S1): study 1,",
      "5 male Beagles, 5 periods; study 2, 12 male Beagles, 9 periods including",
      "four furosemide arms (furosemide data are not part of the model). Urine was",
      "summed over 24-h collections for the PD analysis. NONMEM 7.3, FOCE.",
      sep = " "
    )
  )

  ini({
    # ---- PK (Table 2) --------------------------------------------------------
    lfdepot <- fixed(log(0.98))
    label("Oral bioavailability (fraction)")
    # Table 2, F = 0.98 FIX; Methods 'Pharmacokinetic Modeling': taken from a
    # previous study and 'used as an input in the current model (value fixed to 0.98)'.

    lka <- log(1.66)
    label("First-order absorption rate constant (1/h)")
    # Table 2, ka = 1.66 /h (RSE 18.7%)

    lcl <- log(0.077)
    label("Total clearance (L/h)")
    # Table 2, CL = 0.077 L/h (RSE 11%)

    lfe <- log(0.611)
    label("Fraction of total clearance excreted unchanged in urine (fraction)")
    # Table 2, Fu = 0.611 (RSE 2.4%), 'Urine fraction of the total clearance'

    lvc <- log(0.145)
    label("Central volume of distribution (L)")
    # Table 2, Vd = 0.145 L (RSE 23.9%)

    lq <- log(0.262)
    label("Intercompartmental clearance (L/h)")
    # Table 2, Q = 0.262 L/h (RSE 16.8%)

    lvp <- log(0.935)
    label("Peripheral volume of distribution (L)")
    # Table 2, VdPER = 0.935 L (RSE 4.6%)

    # ---- Diuresis, direct power model (Table 3, Equation 1) -----------------
    lurprod <- log(221)
    label("Baseline urine output without treatment (mL/day)")
    # Table 3, Baseline = 221 mL/day (RSE 4.0%)

    lslope <- -4.11 * log(10)
    label("Natural log of the slope of the torasemide effect on diuresis (log mL/ug^alpha)")
    # Table 3, slope = -4.11 on the log10 scale (RSE 13%); Methods 'PK/PD modeling
    # of diuresis': '-4.11 corresponds to 10^-4.11 mL/ug or 0.0000776 mL/ug'.
    # Carried on the natural-log scale: -4.11 x ln(10) = -9.4636.

    lalpha <- log(2.05)
    label("Power exponent of the torasemide effect on diuresis (unitless)")
    # Table 3, alpha_Tora = 2.05 (RSE 9%)

    # ---- Natriuresis, indirect response with resistance (Table 4) -----------
    lkin <- log(3.67)
    label("Zero-order rate constant for production of natriuresis (mEq/h^2)")
    # Table 4, kforNa = 3.67 mEq/h2 (RSE 3.65%)

    lkout <- log(4.99)
    label("First-order rate constant for loss of natriuresis (1/h)")
    # Table 4, kabsNa = 4.99 /h (RSE 5.76%)

    lemax <- log(4.67)
    label("Maximum fractional stimulation of natriuresis production without resistance (unitless)")
    # Table 4, EmaxNa = 4.67 (RSE 10.53%)

    lec50 <- log(1080)
    label("Urinary torasemide amount producing half-maximal natriuresis stimulation (ug)")
    # Table 4, EC50Na = 1080 ug (RSE 5.97%)

    lhill <- log(2.57)
    label("Hill coefficient of the natriuresis stimulation (unitless)")
    # Table 4, alpha_ToraNa = 2.57 (RSE 12.50%)

    lkres_on <- log(0.0811)
    label("Zero-order rate constant for activation of diuretic resistance (1/h)")
    # Table 4, kresOn = 0.0811 /h (RSE 12.94%)

    lkres_off <- log(0.0894)
    label("First-order rate constant for loss of diuretic resistance (1/h)")
    # Table 4, kresOff = 0.0894 /h (RSE 11.21%)

    lres_thres <- log(0.0547)
    label("Threshold of the fractional natriuresis stimulation above which resistance is activated (unitless)")
    # Table 4, Threshold = 0.0547 (RSE 12.31%); Results: 'a fractional threshold of 0.055 or 5.5%'

    # ---- Inter-individual variability ---------------------------------------
    # Tables 2-4 report IIV as CV% of a log-normal (exponential) random effect;
    # Methods converts with CV% = 100 x sqrt(exp(omega^2) - 1), so
    # omega^2 = log(1 + CV^2). Covariances were not estimated (diagonal omega).
    etalka ~ 0.9517 # Table 2, ka IIV 126.1% CV
    etalcl ~ 0.346 # Table 2, CL IIV 64.3% CV
    etalfe ~ 0.007369 # Table 2, Fu IIV 8.6% CV
    etalurprod ~ 0.0178 # Table 3, Baseline IIV 13.4% CV
    etalslope ~ 0.1904 # Table 3, slope IIV 45.8% CV, applied multiplicatively to the linear-scale slope (see vignette)
    etalalpha ~ 0.00152 # Table 3, alpha_Tora IIV 3.9% CV
    etalkin ~ 0.03695 # Table 4, kforNa IIV 19.4% CV
    etalemax ~ 0.3289 # Table 4, EmaxNa IIV 62.4% CV

    # ---- Residual error -----------------------------------------------------
    propSd <- 0.184
    label("Proportional residual error, plasma torasemide (fraction)")
    # Table 2, 'Proportional error residual in plasma' = 18.4%

    propSd_Aurine <- 0.225
    label("Proportional residual error, urinary torasemide amount (fraction)")
    # Table 2, 'Proportional error residual in urine' = 22.5%

    addSd_Aurine <- 3.68
    label("Additive residual error, urinary torasemide amount (ug)")
    # Table 2, 'Additive error residual in urine' = 3.68 ug

    propSd_urine_vol <- 0.432
    label("Proportional residual error, daily urine volume (fraction)")
    # Table 3, 'Proportional error residual' = 43.2%

    propSd_qna <- 0.305
    label("Proportional residual error, daily urinary sodium (fraction)")
    # Table 4, 'Proportional error residual' = 30.5%
  })
  model({
    # ---- Individual parameters ----------------------------------------------
    fdepot <- exp(lfdepot)
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    fe <- exp(lfe + etalfe)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)

    urprod <- exp(lurprod + etalurprod)
    slope <- exp(lslope + etalslope)
    alpha <- exp(lalpha + etalalpha)

    kin <- exp(lkin + etalkin)
    kout <- exp(lkout)
    emax <- exp(lemax + etalemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)
    kres_on <- exp(lkres_on)
    kres_off <- exp(lkres_off)
    res_thres <- exp(lres_thres)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- PK: two compartments plus urinary excretion (Figure 1) -------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    # Renal arm = fe x CL. `urine` holds the torasemide voided since the start of
    # the current 24-h urine collection: the event table must reset it (and
    # `qna`) to zero at each collection boundary with an evid = 5, amt = 0
    # record (Results 'Pharmacodynamic Modeling': 'The bladder was considered
    # empty after each urination'; Discussion: urine volumes, sodium and
    # torasemide 'were summed over 24 h periods').
    d/dt(urine) <- fe * kel * central
    f(depot) <- fdepot

    # Qurine, the quantity of torasemide excreted in urine (ug) that drives both
    # PD sub-models (Methods 'Pharmacodynamic Modeling').
    Aurine <- 1000 * urine

    # ---- Natriuresis: indirect response with resistance (Equations 2-8) -----
    # Equation 8: resistance lowers the efficacy, Emax_apparent = Emax x (1 - R).
    emax_app <- emax * (1 - rna)
    # Equation 5 (sigmoid Emax on Qurine), with the apparent Emax.
    ena_tora <- emax_app * Aurine^hill / (ec50^hill + Aurine^hill)
    # Resistance activates only while the fractional stimulation exceeds the
    # threshold; otherwise kResON = 0 (Methods, Equation 7 text).
    kres_on_act <- 0
    if (ena_tora > res_thres) kres_on_act <- kres_on
    ena(0) <- kin / kout # Equation 4: baseline natriuresis = kforNa / kabsNa
    d/dt(ena) <- kin * (1 + ena_tora) - kout * ena # Equation 6
    d/dt(rna) <- kres_on_act - kres_off * rna # Equation 7
    d/dt(qna) <- ena # Equation 2

    # ---- Diuresis: direct power model (Equation 1) --------------------------
    # Daily urine volume (mL per 24-h collection) = Baseline + slope x Qurine^alpha,
    # meaningful at the END of each 24-h collection when Aurine is the daily
    # excreted amount.
    urine_vol <- urprod + slope * Aurine^alpha

    # ---- Observations -------------------------------------------------------
    Cc <- 1000 * central / vc # ug/L

    Cc ~ prop(propSd)
    Aurine ~ add(addSd_Aurine) + prop(propSd_Aurine)
    urine_vol ~ prop(propSd_urine_vol)
    qna ~ prop(propSd_qna)
  })
}
