Gao_2012_exenatide_glucose_insulin_rat <- function() {
  description <- paste(
    "Preclinical (rat). Integrated target-mediated PK / glucose-insulin",
    "feedback PD model for exendin-4 (exenatide) in male Sprague-Dawley rats",
    "given a 2-h continuous intravenous infusion of exendin-4 with an",
    "intravenous D-glucose challenge at 30 min (Gao 2012). The",
    "pharmacokinetic layer is the same target-mediated drug disposition",
    "(TMDD) structure as Gao_2012_exenatide_rat.R -- free drug binds GLP-1R",
    "(kon / koff) to a complex that is internalised and degraded (kint),",
    "distributes nonspecifically (k12 / k21) and is eliminated linearly",
    "(kel) -- carried at the Table 2 rat estimates and driving the",
    "pharmacodynamic layer through the free plasma concentration Cc. The",
    "subcutaneous depot of Gao_2012_exenatide_rat.R is deliberately absent:",
    "the pharmacodynamic study was intravenous-infusion only. The",
    "pharmacodynamic layer is a two-pool glucose-insulin feedback model in",
    "which glucose stimulates insulin secretion and insulin stimulates",
    "glucose disposal, each through a linear stimulation factor. Exendin-4",
    "acts by multiplying the glucose-driven insulin-secretion term by",
    "(1 + Sd), where Sd is a BIPHASIC (bell-shaped) Adair function of the",
    "free plasma concentration -- so insulin secretion rises with dose to a",
    "maximum near 7.4 nM and falls again at higher concentrations. The",
    "insulinotropic effect is glucose-dependent and is SHUT OFF whenever",
    "glucose is at or below its baseline (Gao 2012 Discussion), which is",
    "encoded by rectifying the glucose-elevation term at zero. Baseline",
    "glucose and insulin are per-dose-group covariates (FPG / INS_BL); the",
    "defaults are the saline-arm values. The paper's alternative driving",
    "function based on the drug-receptor complex (eq. 8) was rejected on AIC",
    "and is NOT extracted. Naive-pooled fit to mean profiles in ADAPT II;",
    "there is no between-subject variability and the residual-variance",
    "parameters were not reported."
  )
  reference <- paste(
    "Gao W, Jusko WJ (2012).",
    "Target-mediated Pharmacokinetic and Pharmacodynamic Model of Exendin-4",
    "in Rats, Monkeys, and Humans.",
    "Drug Metabolism and Disposition 40(5):990-997.",
    "doi:10.1124/dmd.111.042291.",
    "PMCID PMC3336795.",
    sep = " "
  )
  vignette <- "Gao_2012_exenatide"

  units <- list(
    time          = "min",
    dosing        = paste(
      "Two independent inputs. Exendin-4 is infused intravenously into",
      "`central` in nmol/min (the paper's pmol/kg/min infusion rates are",
      "converted with the animal's body weight; a 350 g rat at",
      "300 pmol/kg/min receives 0.000105 nmol/min). The D-glucose challenge",
      "is an intravenous bolus into `glucose` in mmol per kg of body weight",
      "(the paper used 5.7 mmol/kg)."
    ),
    concentration = paste(
      "Exendin-4 nmol/L (nM), converted from mass using a molecular weight of",
      "4186.6 g/mol; glucose mmol/L (mM); insulin nmol/L (nM). Figures 5 and 8",
      "of the source plot insulin in pmol/L, so multiply Ic by 1000 to compare",
      "against them. The drug and insulin scales are fixed by the units of the",
      "Table 2 and Table 4 estimates (kon in 1/(nM*min), SIns in 1/nM, SGlu in",
      "1/mM) and must not be rescaled independently."
    )
  )

  # The drug compartments are ABSOLUTE amounts (nmol) paired with the
  # absolute central volume Vc = 43.2 mL of Table 2. The `glucose`
  # compartment is instead an amount PER KILOGRAM of body weight (mmol/kg),
  # because the paper reports the glucose volume of distribution per
  # kilogram (Vg = 0.208 l/kg, Table 4) and the glucose challenge per
  # kilogram (5.7 mmol/kg, Methods). Keeping glucose on the per-kilogram
  # scale reproduces the paper's own initial condition
  # Glu(0) = Dose/Vg + Gb = 5.7/0.208 + Gb without needing a body-weight
  # covariate. `insulin` carries no volume at all: eq. 6 is written directly
  # on the insulin CONCENTRATION, so that state is nmol/L rather than an
  # amount.
  compartmentData <- list(
    central     = list(analyte = "exenatide", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "exenatide", units = "nmol", specimen = "tissue", verified = TRUE),
    complex     = list(analyte = "exenatide/GLP-1R complex", units = "nmol", specimen = "tissue", verified = TRUE),
    glucose     = list(analyte = "glucose", units = "mmol/kg body weight", specimen = "plasma", verified = TRUE),
    insulin     = list(analyte = "insulin", units = "nmol/L", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    FPG = list(
      description        = paste(
        "Baseline (predose) plasma glucose concentration, the paper's Gb.",
        "Sets both the glucose initial condition glucose(0) = Gb * Vg and,",
        "through kinG = koutG * Gb, the zero-order glucose production rate,",
        "and is the threshold above which the insulinotropic effect of",
        "exendin-4 switches on. Time-fixed per dose group."
      ),
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Gao 2012 Methods states only that 'the baseline conditions Gb and Ib",
        "were fixed as the measured predose values' -- no numeric value for",
        "Gb appears in the text or in any table. The default carried here",
        "(9.60 mmol/L) and the per-arm values used in the vignette were read",
        "off the pre-challenge plateau of Fig. 5 by pixel-calibrating the",
        "axis labels: saline 9.60, 3 pmol/kg/min 11.53, 30 pmol/kg/min 10.42,",
        "300 pmol/kg/min 9.92, 3000 pmol/kg/min 9.76 mmol/L. These are",
        "DIGITISED, not published, values -- see vignette Assumptions and",
        "deviations. The digitisation is cross-checked by the paper's own",
        "initial condition: 5.7/0.208 + 9.60 = 37.0 mmol/L, which matches the",
        "plotted post-challenge peak of the saline panel."
      ),
      source_name        = "Gb"
    ),
    INS_BL = list(
      description        = paste(
        "Baseline (predose) plasma insulin concentration, the paper's Ib.",
        "Sets the insulin initial condition insulin(0) = Ib, the zero-order",
        "insulin secretion rate kinI = koutI * Ib, and the reference level",
        "from which insulin stimulates glucose disposal. Time-fixed per dose",
        "group."
      ),
      units              = "pmol/L (divided by 1000 inside model() to the nmol/L scale on which SIns is expressed)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "As for FPG, Ib is not printed anywhere in Gao 2012. The default",
        "(304 pmol/L) and the per-arm values in the vignette were digitised",
        "from the pre-challenge plateau of Fig. 5 (right column): saline 304,",
        "3 pmol/kg/min 399, 30 pmol/kg/min 330, 300 pmol/kg/min 287,",
        "3000 pmol/kg/min 329 pmol/L. DIGITISED, not published. The canonical",
        "column is in pmol/L to match Fig. 5 and Fig. 8, which are plotted in",
        "picomoles; the model divides by 1000 because SIns is reported in",
        "1/nM (Table 4)."
      ),
      source_name        = "Ib"
    )
  )

  population <- list(
    species       = "rat (male Sprague-Dawley)",
    n_subjects    = 8L,
    n_studies     = 1L,
    weight_range  = "80-420 g",
    disease_state = "Healthy (normoglycaemic) male Sprague-Dawley rats.",
    dose_range    = paste(
      "Two-hour continuous intravenous infusion of saline or of exendin-4 at",
      "3, 30, 300 or 3000 pmol/kg/min, with an intravenous D-glucose",
      "challenge of 5.7 mmol/kg given over 2-3 min (0.5 mL/min) beginning",
      "30 min after the start of the infusion. n = 4-8 per arm."
    ),
    notes = paste(
      "Gao 2012 Materials and Methods, 'In another study, male SD rats'. No",
      "exendin-4 concentrations were measured in this study; the drug",
      "concentration profiles driving the pharmacodynamics were SIMULATED",
      "from the rat TMDD PK model of Table 2 (Gao 2012 Results,",
      "'Pharmacodynamics', and Fig. 6 top). Plasma glucose was measured by",
      "immobilised oxidase chemistry (YSI 2300 Stat Plus) and insulin by",
      "radioimmunoassay (Linco). The model was fitted to MEAN profiles by",
      "naive pooling in ADAPT II with the maximum-likelihood method, so",
      "n_subjects records the largest per-arm n rather than a pooled",
      "population size. The rats in this pharmacodynamic study (80-420 g)",
      "span a wider weight range than those in the pharmacokinetic study",
      "(350-370 g) whose parameters are reused here."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # PHARMACOKINETICS -- Gao 2012 Table 2 (rat), CV% in parentheses in
    # the source. All rates are 1/min. The rat model runs on the
    # nanomolar concentration scale (kon in 1/(nM*min), Rtot in nmol/L),
    # confirmed by koff/kon = 0.0153/0.0207 = 0.74 nM, the KD quoted in
    # Gao 2012 Results 'Rat PK'. These are the same estimates as
    # Gao_2012_exenatide_rat.R; the paper reused the PK fit unchanged to
    # generate the driving concentrations for the PD analysis.
    # ------------------------------------------------------------------
    lkel  <- log(0.0839)  ; label("Linear elimination rate constant from the central compartment (kel, 1/min)")     # Table 2: kel = 0.0839 (CV 10%)
    lk12  <- log(0.0282)  ; label("Transfer rate constant central -> peripheral1 (kpt, 1/min)")                     # Table 2: kpt = 0.0282 (CV 15%)
    lk21  <- log(0.0213)  ; label("Transfer rate constant peripheral1 -> central (ktp, 1/min)")                     # Table 2: ktp = 0.0213 (CV 5%)
    lvc   <- log(0.0432)  ; label("Central volume of distribution (Vc, L)")                                         # Table 2: Vc = 43.2 mL = 0.0432 L (CV 12%)
    lkon  <- log(0.0207)  ; label("Second-order association rate constant of exendin-4 with GLP-1R (kon, 1/(nM*min))") # Table 2: kon = 0.0207 (CV 42%)
    lkoff <- log(0.0153)  ; label("First-order dissociation rate constant of the drug-receptor complex (koff, 1/min)") # Table 2: koff = 0.0153 (CV 206%)
    lkint <- log(0.0966)  ; label("Internalisation / degradation rate constant of the drug-receptor complex (kint, 1/min)") # Table 2: kint = 0.0966 (CV 38%)
    lrtot <- log(5.21)    ; label("Total GLP-1R concentration, held constant (Rtot, nmol/L)")                       # Table 2: Rtot = 5.21 nmol/L (CV 5%)

    # ------------------------------------------------------------------
    # PHARMACODYNAMICS -- Gao 2012 Table 4, the final integrated TMDD
    # PK/PD fit to the glucose and insulin profiles of all five infusion
    # arms. CV% in parentheses in the source. The zero-order production
    # rates kinG and kinI are NOT estimated: Gao 2012 Methods gives the
    # steady-state relationships kinG = koutG * Gb and kinI = koutI * Ib,
    # so they are derived in model() from the baseline covariates.
    # ------------------------------------------------------------------
    lkout_glucose <- log(0.046)   ; label("First-order glucose elimination rate constant (koutG, 1/min)")  # Table 4: koutG = 0.046 (CV 9%)
    lkout_insulin <- log(0.483)   ; label("First-order insulin elimination rate constant (koutI, 1/min)")  # Table 4: koutI = 0.483 (CV 50%)
    lvg           <- log(0.208)   ; label("Glucose apparent volume of distribution (Vg, L/kg body weight)") # Table 4: Vg = 0.208 l/kg (CV 5%)

    # Linear inter-pool stimulation factors of the feedback loop
    # (Gao 2012 eqs. 5-6). Each enters as (1 + S * (driver - baseline)).
    # Note the two are expressed on DIFFERENT concentration scales,
    # matching the units of the pools they read: SGlu is per millimolar
    # (glucose) and SIns is per nanomolar (insulin).
    lsstim_glucose_insulin <- log(0.0684) ; label("Linear stimulation factor of glucose on insulin secretion (SGlu, 1/mM)") # Table 4: SGlu = 0.0684 (CV 20%)
    lsstim_insulin_glucose <- log(0.157)  ; label("Linear stimulation factor of insulin on glucose disposal (SIns, 1/nM)")  # Table 4: SIns = 0.157 (CV 46%)

    # Biphasic (bell-shaped) Adair drug-effect function, Gao 2012 eq. 7:
    #   Sd = Smax * C / (k1 + C + k2 * C^2)
    # The quadratic denominator term is what makes the response
    # bell-shaped; Sd is maximal at C = sqrt(k1/k2) = 7.35 nM. The two
    # constants carry reciprocal units and are NOT two instances of one
    # quantity. The paper interprets them mechanistically (Discussion) as
    # the first and second dissociation constants of a receptor that can
    # bind a second exendin-4 molecule.
    lsmax    <- log(4.67)   ; label("Maximal stimulation factor of the Adair drug-effect function (Smax, unitless)") # Table 4: Smax = 4.67 (CV 30%)
    lkadair1 <- log(0.826)  ; label("First (linear) constant of the Adair function (k1, nM)")                        # Table 4: k1 = 0.826 nM (CV 71%)
    lkadair2 <- log(0.0153) ; label("Second (quadratic, self-inhibition) constant of the Adair function (k2, 1/nM)")  # Table 4: k2 = 0.0153 1/nM (CV 69%)

    # ------------------------------------------------------------------
    # Residual error. Gao 2012 Methods states the variance model
    # Vi = (sigma1 + sigma2 * Y)^2 -- a combined additive plus
    # proportional standard deviation -- but reports NEITHER sigma1 nor
    # sigma2, for any of the three outputs. All six terms are therefore
    # encoded at zero rather than invented. See vignette Assumptions and
    # deviations.
    # ------------------------------------------------------------------
    addSd  <- fixed(0) ; label("Additive residual SD on exendin-4 concentration (nmol/L; ZERO - not reported)")     # Methods: Vi = (sigma1 + sigma2 * Y)^2; sigma1 not reported
    propSd <- fixed(0) ; label("Proportional residual SD on exendin-4 concentration (fraction; ZERO - not reported)") # Methods: Vi = (sigma1 + sigma2 * Y)^2; sigma2 not reported
    addSd_Gc  <- fixed(0) ; label("Additive residual SD on glucose concentration (mmol/L; ZERO - not reported)")       # Methods: Vi = (sigma1 + sigma2 * Y)^2; sigma1 not reported
    propSd_Gc <- fixed(0) ; label("Proportional residual SD on glucose concentration (fraction; ZERO - not reported)")   # Methods: Vi = (sigma1 + sigma2 * Y)^2; sigma2 not reported
    addSd_Ic  <- fixed(0) ; label("Additive residual SD on insulin concentration (nmol/L; ZERO - not reported)")       # Methods: Vi = (sigma1 + sigma2 * Y)^2; sigma1 not reported
    propSd_Ic <- fixed(0) ; label("Proportional residual SD on insulin concentration (fraction; ZERO - not reported)")   # Methods: Vi = (sigma1 + sigma2 * Y)^2; sigma2 not reported

    # No inter-individual variability: the analysis is a naive-pooled fit
    # to mean glucose and insulin profiles, so no OMEGA block exists.
  })

  model({
    # 1. Individual parameters (no IIV; naive-pooled mean-data fit).
    kel  <- exp(lkel)
    k12  <- exp(lk12)
    k21  <- exp(lk21)
    vc   <- exp(lvc)
    kon  <- exp(lkon)
    koff <- exp(lkoff)
    kint <- exp(lkint)
    rtot <- exp(lrtot)

    kout_glucose <- exp(lkout_glucose)
    kout_insulin <- exp(lkout_insulin)
    vg           <- exp(lvg)

    sstim_glucose_insulin <- exp(lsstim_glucose_insulin)
    sstim_insulin_glucose <- exp(lsstim_insulin_glucose)

    smax    <- exp(lsmax)
    kadair1 <- exp(lkadair1)
    kadair2 <- exp(lkadair2)

    # 2. Per-arm baselines. FPG is the paper's Gb in mmol/L and needs no
    # conversion; INS_BL is the paper's Ib supplied in the canonical
    # pmol/L and is divided by 1000 to the nmol/L scale on which SIns
    # (1/nM, Table 4) and the insulin state are expressed.
    gb <- FPG
    ib <- INS_BL / 1000

    # 3. Zero-order production rates, derived rather than estimated:
    # Gao 2012 Methods, "kinG = koutG * Gb and kinI = koutI * Ib, where Gb
    # and Ib are the basal glucose and insulin concentrations". kinG is
    # multiplied by vg because the `glucose` state is an amount per kg
    # while eq. 5 is written on the concentration.
    kin_glucose <- kout_glucose * gb * vg
    kin_insulin <- kout_insulin * ib

    # 4. Free (unoccupied) receptor concentration -- Gao 2012 eq. 3 writes
    # the binding term as kon * (Rtot - RC) * C, so the receptor pool is
    # constant and free receptor is the total pool less the occupied
    # fraction. RC is a CONCENTRATION in the source and `complex` is the
    # corresponding AMOUNT, hence the division by vc.
    rfree <- rtot - complex / vc

    # 5. Pharmacokinetics -- Gao 2012 eqs. 1-3, rewritten from the
    # published concentration form (dC/dt) to the amount form used by
    # rxode2 by multiplying through by Vc. There is no depot: the
    # pharmacodynamic study administered exendin-4 by intravenous
    # infusion only, so the subcutaneous input of eq. 4 does not apply.
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1 -
                          kon * rfree * central + koff * complex
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(complex)     <-  kon * rfree * central - (koff + kint) * complex

    # 6. Free exendin-4 concentration in plasma (nmol/L). This is the
    # driving function of the pharmacodynamic model: Gao 2012 tested the
    # drug-receptor complex RC as the driver instead (eq. 8) and rejected
    # it on AIC and on prior Goto-Kakizaki rat evidence, so eq. 8 is not
    # implemented here.
    Cc <- central / vc

    # 7. Biphasic Adair drug-effect function -- Gao 2012 eq. 7.
    sadair <- smax * Cc / (kadair1 + Cc + kadair2 * Cc * Cc)

    # 8. Glucose elevation above its own baseline, RECTIFIED AT ZERO. The
    # printed eq. 6 carries the bare difference (Glu - Gb), but Gao 2012
    # Discussion is explicit that "when glucose is not higher than basal
    # values, the effect is shut off" -- without the rectification the
    # term turns negative during the post-challenge undershoot and the
    # drug would spuriously SUPPRESS insulin secretion.
    glucose_conc <- glucose / vg
    delta_g      <- glucose_conc - gb

    # 9. Glucose-insulin feedback -- Gao 2012 eqs. 5-6. Insulin
    # stimulates glucose disposal linearly; glucose stimulates insulin
    # secretion linearly, and that glucose-driven term is amplified by
    # the drug through (1 + Sd).
    d/dt(glucose) <- kin_glucose -
                     kout_glucose * (1 + sstim_insulin_glucose * (insulin - ib)) * glucose
    d/dt(insulin) <- kin_insulin *
                     (1 + sstim_glucose_insulin * delta_g * (delta_g > 0) * (1 + sadair)) -
                     kout_insulin * insulin

    # 10. Initial conditions -- Gao 2012 eqs. 5-6. Before the challenge
    # both pools sit at their measured predose values; the D-glucose
    # bolus is supplied as a dosing record into `glucose` (mmol/kg), which
    # reproduces the printed Glu(0) = Dose/Vg + Gb.
    glucose(0) <- gb * vg
    insulin(0) <- ib

    # 11. Observations: exendin-4 (nmol/L), glucose (mmol/L) and insulin
    # (nmol/L) in plasma.
    Gc <- glucose_conc
    Ic <- insulin

    Cc ~ add(addSd) + prop(propSd)
    Gc ~ add(addSd_Gc) + prop(propSd_Gc)
    Ic ~ add(addSd_Ic) + prop(propSd_Ic)
  })
}
