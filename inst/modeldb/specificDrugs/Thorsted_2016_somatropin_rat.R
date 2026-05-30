Thorsted_2016_somatropin_rat <- function() {
  description <- paste(
    "Preclinical (hypophysectomized Sprague-Dawley rat).",
    "Mixed-effects PKPD model for recombinant human growth hormone",
    "(rhGH / somatropin) describing PK as a two-compartment model with",
    "parallel linear (CL) and Michaelis-Menten (Vmax, KM) elimination,",
    "parallel first-order subcutaneous absorption (ka1 direct path, ka2",
    "delayed through one transit compartment, with bioavailabilities",
    "F1 and F2), an indirect response model for IGF-1 induction",
    "(stimulation of kin via a three-compartment effect-delay chain",
    "feeding an Emax/EC50 stimulation), and a linear bodyweight-gain",
    "model driven by IGF-1 above baseline. Reference rat body weight",
    "is 0.1 kg (100 g) and the allometric exponents (0.75 / 1.0) are",
    "fixed."
  )
  reference <- paste(
    "Thorsted A, Thygesen P, Agerso H, Laursen T, Kreilgaard M.",
    "Translational mixed-effects PKPD modelling of recombinant human",
    "growth hormone - from hypophysectomized rat to patients.",
    "Br J Pharmacol. 2016 Jun;173(11):1742-55.",
    "doi:10.1111/bph.13473."
  )
  vignette <- "Thorsted_2016_somatropin_rat"
  units <- list(
    time = "h",
    dosing = "ug",
    concentration = "ng/mL"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling around a reference rat body weight of 0.1 kg",
        "(100 g) per Thorsted 2016 Methods. Clearance terms (CL, Q) and",
        "the Michaelis-Menten maximum elimination capacity (Vmax) scale",
        "with (WT/0.1)^0.75; distribution volumes (Vc, Vp) scale with",
        "(WT/0.1)^1.0; first-order absorption rate constants (ka1, ka2)",
        "scale with (WT/0.1)^(-0.25). KM is not scaled. Exponents are",
        "fixed at their canonical allometric values per Thorsted 2016",
        "Results / Discussion (estimating the exponents did not yield",
        "meaningful improvement over the fixed values). Baseline-only",
        "covariate for PK scaling; the bodyweight state in the model",
        "(bw) tracks rhGH-induced bodyweight gain and is reported in",
        "grams while the WT covariate is in kg.",
        sep = " "
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "rat (Sprague-Dawley, hypophysectomized; male, 90-110 g, age 6-8 weeks)",
    n_subjects     = 230L,
    n_studies      = 6L,
    age_range      = "6-8 weeks at study start (hypophysectomized at 4 weeks; acclimatized 1-2 weeks)",
    weight_range   = "90-110 g (approximately 0.1 kg)",
    sex_female_pct = 0,
    disease_state  = paste(
      "Hypophysectomized rat model of growth hormone deficiency",
      "(pituitary gland surgically removed at age 4 weeks); plasma rhGH",
      "and IGF-1 measured by LOCI immunoassay and ELISA respectively;",
      "bodyweight measured daily during the multiple-dose PD studies.",
      sep = " "
    ),
    dose_range     = paste(
      "rhGH 11-3319 ug as i.v. tail-vein bolus or s.c. injection into",
      "the scruff of the neck; six study cohorts in total (one",
      "single-dose PKPD study + five multiple-dose PD studies of up to",
      "28 daily doses; only doses through 168 h were retained in the",
      "PD model because of anti-drug-antibody formation thereafter).",
      sep = " "
    ),
    regions        = "Denmark (Novo Nordisk A/S, Maaloev)",
    notes          = paste(
      "230 male hypophysectomized Sprague-Dawley rats from Taconic M&B",
      "(Ejby, Denmark); 304 rhGH plasma concentrations, 717 IGF-1",
      "plasma concentrations, and 1248 bodyweight measurements. See",
      "Thorsted 2016 Methods (Animals, Experiments) and Table 1 for",
      "the cohort breakdown."
    )
  )

  ini({
    # PK structural parameters: Thorsted 2016 Table 2. Reference rat 0.1 kg.
    lcl     <- log(0.0285);   label("Linear clearance CL (L/h) at 0.1 kg rat")             # Table 2
    lvmax   <- log(11.5);     label("Michaelis-Menten Vmax (ug/h) at 0.1 kg rat")          # Table 2
    lkm     <- log(358);      label("Michaelis-Menten KM (ug/L = ng/mL)")                  # Table 2
    lvc     <- log(0.0069);   label("Central volume of distribution Vc (L) at 0.1 kg rat") # Table 2
    lq      <- log(0.0101);   label("Inter-compartmental clearance Q (L/h) at 0.1 kg rat") # Table 2
    lvp     <- log(0.0081);   label("Peripheral volume of distribution Vp (L) at 0.1 kg rat") # Table 2
    lka1    <- log(3.02);     label("First-order SC absorption rate ka1 (1/h) - direct path") # Table 2
    lka2    <- log(1.22);     label("First-order SC absorption rate ka2 (1/h) - delayed via transit") # Table 2
    lfdepot <- log(0.0316);   label("Bioavailability F1 through ka1 SC path (fraction)")   # Table 2
    lfdepot2 <- log(0.833);   label("Bioavailability F2 through ka2 SC path (fraction)")   # Table 2

    # PD structural parameters: Thorsted 2016 Table 2.
    # Only the single-dose PKPD-study Emax is encoded (9.88); see vignette
    # Assumptions and deviations for why the dual-Emax (PKPD vs PD studies)
    # was collapsed to the single-dose value (which is the one the paper
    # carries forward to the human projection in Table 3).
    lkout   <- log(0.0913);   label("First-order IGF-1 degradation rate kout (1/h)")       # Table 2
    lrbase     <- log(29.4);     label("Baseline IGF-1 R0 (ng/mL)")                            # Table 2
    lemax   <- log(9.88);     label("Maximum stimulation of kin relative to baseline (PKPD-study value)") # Table 2
    lec50   <- log(16.3);     label("rhGH effect-compartment concentration for 50% Emax (ug/L = ng/mL)") # Table 2
    lkcplag <- log(0.599);    label("First-order rate constant for the GH effect-delay (CPLAG) chain (1/h)") # Table 2
    lsld    <- log(0.000309); label("Slope of bodyweight gain per (IGF-1 - R0) per h (g per (ng/mL) per h)") # Table 2
    lwtbase <- log(106);      label("Estimated bodyweight at study start (g)")              # Table 2

    # Allometric exponents - fixed at canonical values per Thorsted 2016
    # Methods ("with fixed exponents of 0.75 and -0.25") and Results
    # ("scaled according to bodyweight with fixed exponents of 0.75 and
    # 1.0 for clearance processes and distribution volumes").
    e_wt_cl   <- fixed(0.75);  label("Allometric exponent on CL, Q, Vmax (unitless)")      # Methods / Results
    e_wt_vc   <- fixed(1.0);   label("Allometric exponent on Vc, Vp (unitless) - rat")     # Methods / Results
    e_wt_ka   <- fixed(-0.25); label("Allometric exponent on ka1, ka2 (unitless)")         # Methods

    # IIV - Thorsted 2016 Table 2. CV% converted to omega^2 via
    # log(CV^2 + 1). The CL-Vp block is correlated with rho = -0.568
    # (Results); the off-diagonal is rho * sqrt(var_cl) * sqrt(var_vp).
    # CL 11.6% CV -> 0.01337; Vp 18.4% CV -> 0.03330; cov -> -0.01199.
    etalcl + etalvp ~ c(0.01337,
                        -0.01199, 0.03330)                                                  # Table 2 IIV CL 11.6% CV, Vp 18.4% CV, correlation -0.568
    etalka2  ~ 0.008613                                                                     # Table 2 IIV ka2 9.3% CV: log(0.093^2 + 1)
    etalrbase   ~ 0.028504                                                                     # Table 2 IIV R0 17.0% CV: log(0.170^2 + 1)
    etalwtbase ~ 0.002697                                                                   # Table 2 IIV WT_BASE 5.2% CV: log(0.052^2 + 1)

    # Residual error - Thorsted 2016 Table 2 (PKPD-study values).
    propSd        <- 0.233;    label("Proportional residual error on rhGH plasma concentration (fraction)") # Table 2 Prop-PK
    addSd         <- 0.279;    label("Additive residual SD on rhGH plasma concentration (ng/mL = ug/L)")    # Table 2 Add-PK
    addSd_IGF1    <- 21.3;     label("Additive residual SD on IGF-1 plasma concentration (ng/mL) - PKPD study") # Table 2 Add-PD (IGF-1) PKPD
    addSd_BW      <- 2.38;     label("Additive residual SD on bodyweight (g)")                              # Table 2 Add-PD Bodyweight
  })

  model({
    # 1. Allometric scaling (WT covariate) around reference rat 0.1 kg.
    wt_norm <- WT / 0.1

    # 2. Individual PK parameters
    cl   <- exp(lcl   + etalcl) * wt_norm^e_wt_cl
    vmax <- exp(lvmax)          * wt_norm^e_wt_cl
    km   <- exp(lkm)
    vc   <- exp(lvc)            * wt_norm^e_wt_vc
    q    <- exp(lq)             * wt_norm^e_wt_cl
    vp   <- exp(lvp   + etalvp) * wt_norm^e_wt_vc
    ka1  <- exp(lka1)           * wt_norm^e_wt_ka
    ka2  <- exp(lka2  + etalka2)* wt_norm^e_wt_ka
    fdepot  <- exp(lfdepot)
    fdepot2 <- exp(lfdepot2)

    # Individual PD parameters
    kout   <- exp(lkout)
    rbase     <- exp(lrbase + etalrbase)
    emax   <- exp(lemax)
    ec50   <- exp(lec50)
    kcplag <- exp(lkcplag)
    sld    <- exp(lsld)
    wtbase <- exp(lwtbase + etalwtbase)

    # 3. Micro-constants for the 2-compartment disposition
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 4. ODE system
    # rhGH disposition: parallel linear (kel * central) and saturable
    # (Vmax * Cc / (KM + Cc)) elimination from central; central <-> peripheral1.
    # Concentration Cc has units ug/L = ng/mL.
    Cc <- central / vc

    d/dt(depot)       <- -ka1 * depot
    d/dt(depot2)      <- -ka2 * depot2
    d/dt(transit1)    <-  ka2 * depot2 - ka2 * transit1
    d/dt(central)     <-  ka1 * depot + ka2 * transit1 -
                          kel * central - k12 * central + k21 * peripheral1 -
                          (vmax * Cc) / (km + Cc)
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # GH effect-delay (CPLAG) chain: three first-order pools transferring
    # the plasma rhGH concentration to a delayed effect concentration
    # that drives the IGF-1 Emax stimulation. The final effect3
    # corresponds to CPLAG,GH in Thorsted 2016.
    d/dt(effect1) <- kcplag * (Cc - effect1)
    d/dt(effect2) <- kcplag * (effect1 - effect2)
    d/dt(effect3) <- kcplag * (effect2 - effect3)

    # IGF-1 indirect response: stimulation of kin (= kout * R0) by the
    # delayed GH effect concentration via Emax / EC50.
    stim <- emax * effect3 / (ec50 + effect3)
    d/dt(igf1) <- kout * rbase * (1 + stim) - kout * igf1
    igf1(0)    <- rbase

    # Bodyweight gain: linear in (IGF-1 - R0), only positive contributions
    # (drug-induced excursions above baseline).
    d/dt(bw) <- sld * (igf1 - rbase)
    bw(0)    <- wtbase

    # 5. Bioavailability: SC dose into depot (ka1 path) gets F1; SC dose
    # into depot2 (ka2 path) gets F2. IV doses are placed directly into
    # central with no F applied. To simulate an SC dose at time t with
    # amount X, provide TWO dose events: cmt=depot, amt=X and cmt=depot2,
    # amt=X; the f() factors split the systemic input between the two
    # parallel absorption paths.
    f(depot)  <- fdepot
    f(depot2) <- fdepot2

    # 6. Observations and residual error
    IGF1 <- igf1
    BW   <- bw

    Cc   ~ add(addSd) + prop(propSd)
    IGF1 ~ add(addSd_IGF1)
    BW   ~ add(addSd_BW)
  })
}
