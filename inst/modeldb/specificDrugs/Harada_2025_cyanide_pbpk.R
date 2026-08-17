Harada_2025_cyanide_pbpk <- function() {
  description <- "PBPK (semi-mechanistic, 6 states). Cyanide (CN) toxicokinetics during inhalation of hydrogen cyanide (HCN) gas, applied to forensic reconstruction of fire-related deaths (Harada et al. 2025, Forensic Toxicol). Inhaled HCN is absorbed across the alveolar membrane driven by a plasma:air partition coefficient, distributes from plasma into erythrocytes (saturable binding), liver, muscle and a lumped other-tissue compartment, and is detoxified in the liver by rhodanese-mediated transfer of sulfane sulfur from a dynamic endogenous sulfur-donor pool to give thiocyanate. States carry concentrations (umol/L), not amounts, reproducing the authors' mass-balance parameterisation. The model has no dosing events: exposure is driven entirely by the inhaled air concentration covariate CONC_HCN_PPM. Two observables are produced, whole-blood arterial CN (Carterial, the post-mortem LEFT cardiac blood proxy) and whole-blood mixed-venous CN (Cvenous, the RIGHT cardiac blood proxy); the arterial-venous gap is the paper's diagnostic signal for a brief high-concentration exposure. Structure and every parameter value are inherited from Stamyr et al. 2015 (Arch Toxicol 89:1287-1296) as transcribed in the Harada 2025 Supplemental Material 1 Python script, which is the authoritative on-disk implementation. Deterministic typical-value model: no IIV and no residual error are reported."
  reference <- "Harada K, Tokugawa Y, Henmi K, Miyashita Y, Sakahashi Y, Nishihori T, Sakamoto Y, Yang C, Isobe Y, Sugimoto K, Nakama K, Katada R, Matsumoto H. Analysis of cyanide exposure status in fire-related deaths using a physiologically based pharmacokinetic model. Forensic Toxicol. 2025;43(2):303-312. doi:10.1007/s11419-025-00713-8. Model structure and parameters inherited from Stamyr K, Mork AK, Johanson G. Physiologically based pharmacokinetic modeling of hydrogen cyanide levels in human breath. Arch Toxicol. 2015;89:1287-1296. doi:10.1007/s00204-014-1310-y (not open access; all values used here are taken from the Harada 2025 Supplemental Material 1 Python script, which lists them explicitly)."
  vignette <- "Harada_2025_cyanide_pbpk"
  units <- list(
    time = "min",
    dosing = "ppm (inhaled HCN air concentration; the model has no dosing events)",
    concentration = "ug/mL"
  )

  # The endogenous sulfane-sulfur donor pool that rhodanese consumes when
  # converting CN to thiocyanate. Cyanide-specific detoxification mechanism
  # with no canonical counterpart in inst/references/compartment-names.md, so
  # it is declared paper-specific rather than forced onto an unrelated name.
  paper_specific_compartments <- c("sulfur_donor")

  covariateData <- list(
    CONC_HCN_PPM = list(
      description        = "Inhaled hydrogen cyanide gas concentration in ambient air",
      units              = "ppm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "The sole exposure driver; the model has no dosing events. Converted inside model() to umol/L by dividing by 24 (the authors' conversion, Supplemental Material 1 header comment: 1 umol/L = 24 ppm = 0.027 ug/mL). Held constant from t = 0 in the published analysis, which is how Fig. 3 and the Table 2 estimates were generated; supplied as a covariate column so a user can instead pass a time-varying fire-scene profile, the fluctuation the paper's Discussion names as an unaddressed limitation. Harada 2025 estimated this quantity per case over a 0-18,000 ppm grid in 12 ppm steps.",
      source_name        = "Ca"
    )
  )

  # Issue #482. Every state holds a CONCENTRATION in umol/L rather than an
  # amount, because the source ODEs divide each mass-balance term by the
  # compartment volume; compare the `glucose` canonical, which is likewise a
  # concentration state. verified = TRUE: analyte and matrix were read off the
  # Supplemental Material 1 variable-definition comments.
  compartmentData <- list(
    plasma       = list(analyte = "cyanide", units = "umol/L", specimen = "plasma", verified = TRUE),
    erythrocytes = list(analyte = "cyanide", units = "umol/L", specimen = "blood cell", verified = TRUE),
    liver        = list(analyte = "cyanide", units = "umol/L", specimen = "tissue", verified = TRUE),
    sulfur_donor = list(analyte = "endogenous sulfane-sulfur donors", units = "umol/L", specimen = "tissue", verified = TRUE),
    muscle       = list(analyte = "cyanide", units = "umol/L", specimen = "tissue", verified = TRUE),
    other        = list(analyte = "cyanide", units = "umol/L", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 29L,
    n_studies      = 1L,
    age_range      = "31-89 years",
    age_median     = "70 years",
    sex_female_pct = 37.9,
    disease_state  = "Fire-related deaths examined at forensic autopsy. No resuscitation was attempted before death was confirmed. Carboxyhaemoglobin exceeded the 2.0% non-smoker reference in every case (range 5.4-100%); cyanide was above the 0.2 ug/mL detection limit in the left or right cardiac blood of 23 of 29 cases (79.3%), and thiocyanate was detectable in all 29 (0.92-8.5 ug/mL).",
    dose_range     = "Estimated inhaled HCN air concentrations 84-16,632 ppm with estimated exposure durations 0.05-13.65 min across the 13 cases that had usable paired left/right cardiac cyanide measurements (Table 2).",
    regions        = "Osaka, Japan (Department of Legal Medicine, Osaka University; autopsies April 2014 - March 2020)",
    notes          = "Baseline demographics, paired left/right cardiac cyanide and thiocyanate concentrations, and carboxyhaemoglobin percentages are in Harada 2025 Table 1; the per-case exposure reconstructions are in Table 2. IMPORTANT: the 29 autopsy cases are the APPLICATION dataset, not the estimation dataset. No parameter in this model was fitted to them - every structural and physiological value is inherited unchanged from Stamyr et al. 2015, whose PBPK model was built for healthy adult humans in a controlled HCN-in-breath study. Harada 2025 estimated only two per-case quantities (inhaled air concentration and exposure duration) by grid search against the paired cardiac measurements. The paper's own Discussion flags that the Stamyr physiological parameters may need adjustment for individual decedents (body weight, body-fat percentage, lung capacity) and that post-mortem redistribution of cyanide from lung into left cardiac blood is only implicitly accommodated by treating the arterial compartment as lung + arterial blood combined."
  )

  ini({
    # ----------------------------------------------------------------------
    # ALL values below are from the Harada 2025 Supplemental Material 1
    # Python script parameter block (a Springer electronic-supplementary-
    # material .docx), which names, values and unit-comments each Stamyr 2015
    # parameter explicitly. Stamyr et al. 2015 is not open access; the
    # supplement is the on-disk authority for every number here. Nothing is
    # estimated in this paper, so every parameter is fixed().
    #
    # NOT ENCODED: the supplement's parameter block also defines
    #   Clr = 0.0041  # first-order metabolic clearance of CN via remaining pathways
    # but no ODE in the published script references Clr, so it has no effect on
    # any published result and is deliberately omitted rather than carried as a
    # dangling ini() entry. See the vignette Errata.
    # ----------------------------------------------------------------------

    # Partition coefficients (unitless). lkp_<tissue> is the registered
    # tissue-to-plasma family; lppa is the plasma-to-AIR coefficient that
    # gates pulmonary uptake and is therefore not a member of that family.
    lppa       <- fixed(log(281)); label("Plasma:air partition coefficient (unitless)")            # Suppl. Material 1 (Ppa = 281)
    lkp_liver  <- fixed(log(5.1)); label("Liver:plasma partition coefficient (unitless)")           # Suppl. Material 1 (Php = 5.1)
    lkp_muscle <- fixed(log(2.8)); label("Muscle:plasma partition coefficient (unitless)")          # Suppl. Material 1 (Pmp = 2.8)
    lkp_other  <- fixed(log(5.4)); label("Other-tissue:plasma partition coefficient (unitless)")    # Suppl. Material 1 (Pop = 5.4)

    # Ventilation and blood flows (L/min). q_liver + q_muscle + q_other =
    # 1.6 + 4.9 + 4.2 = 10.7 = q_tot exactly, so the flow table closes.
    lq_alv     <- fixed(log(16.5)); label("Alveolar ventilation (L/min)")                           # Suppl. Material 1 (Qalv = 16.5)
    lq_tot     <- fixed(log(10.7)); label("Cardiac output (L/min)")                                 # Suppl. Material 1 (Qtot = 10.7)
    lq_liver   <- fixed(log(1.6));  label("Liver blood flow (L/min)")                               # Suppl. Material 1 (Qh = 1.6)
    lq_muscle  <- fixed(log(4.9));  label("Muscle blood flow (L/min)")                              # Suppl. Material 1 (Qm = 4.9)
    lq_other   <- fixed(log(4.2));  label("Other-tissue blood flow (L/min)")                        # Suppl. Material 1 (Qo = 4.2)

    # Compartment volumes (L). v_blood = 4.8 L is used only to recombine the
    # plasma and erythrocyte phases into a whole-blood observable; it is
    # deliberately NOT v_plasma + v_erythrocytes (2.9 + 1.9 = 4.8 here, which
    # happens to agree).
    lv_plasma       <- fixed(log(2.9)); label("Plasma compartment volume (L)")                      # Suppl. Material 1 (Vp = 2.9)
    lv_erythrocytes <- fixed(log(1.9)); label("Erythrocyte compartment volume (L)")                 # Suppl. Material 1 (Ve = 1.9)
    lv_liver        <- fixed(log(1.6)); label("Liver compartment volume (L)")                       # Suppl. Material 1 (Vh = 1.6)
    lv_muscle       <- fixed(log(35));  label("Muscle compartment volume (L)")                      # Suppl. Material 1 (Vm = 35)
    lv_other        <- fixed(log(18));  label("Other-tissue compartment volume (L)")                # Suppl. Material 1 (Vo = 18)
    lv_blood        <- fixed(log(4.8)); label("Whole-blood volume used to recombine plasma and erythrocyte phases (L)") # Suppl. Material 1 (Vb = 4.8)

    # Blood-phase split and plasma <-> erythrocyte exchange.
    fp_blood          <- fixed(0.6);       label("Plasma fraction of whole blood (unitless); equivalently 1 - haematocrit") # Suppl. Material 1 (Fp = 0.6)
    lcl_erythrocytes  <- fixed(log(65.9)); label("Plasma-erythrocyte exchange clearance (L/min)")   # Suppl. Material 1 (Cle = 65.9)
    emax_erythrocytes <- fixed(140);       label("Maximum cyanide binding capacity of erythrocytes (umol/L)") # Suppl. Material 1 (Emax = 140)
    kaff_erythrocytes <- fixed(1);         label("Erythrocyte binding affinity constant for cyanide from plasma (umol/L)") # Suppl. Material 1 (Kaff = 1)

    # Hepatic detoxification: rhodanese converts CN to thiocyanate using
    # sulfane sulfur drawn from a pool with zero-order supply and first-order
    # loss, so the pool depletes during a high-concentration exposure.
    lkscn <- fixed(log(0.01));   label("Second-order rate constant for thiocyanate formation (L/(umol*min))") # Suppl. Material 1 (Kscn = 0.01)
    lkfs  <- fixed(log(2.2));    label("Zero-order sulfur-donor formation rate (umol/L/min)")       # Suppl. Material 1 (Kfs = 2.2)
    lkes  <- fixed(log(0.0027)); label("First-order sulfur-donor elimination rate constant (1/min)") # Suppl. Material 1 (Kes = 0.0027)

    # Unit conversions taken verbatim from the supplement's own header
    # comment ("1 umol/L = 27 ug/L = 24 ppm"). Note the authors convert on the
    # molar mass of HCN (27.03 g/mol), not of the CN- anion (26.02 g/mol);
    # keeping their factor is what reproduces their published numbers.
    mw_hcn      <- fixed(27); label("Molar mass used by the authors for the umol/L to ug/mL conversion (g/mol)") # Suppl. Material 1 header comment (1 umol/L = 27 ug/L)
    ppm_per_umol <- fixed(24); label("Air-phase conversion factor (ppm per umol/L)")                # Suppl. Material 1 header comment (1 umol/L = 24 ppm)
  })

  model({
    # ------------------------------------------------------------------
    # Back-transform. Every parameter is fixed, so these are constants;
    # they remain overridable via rxSolve(params = c(...)) for the
    # sensitivity analyses the paper's Discussion calls for.
    # ------------------------------------------------------------------
    ppa       <- exp(lppa)
    php       <- exp(lkp_liver)
    pmp       <- exp(lkp_muscle)
    pop       <- exp(lkp_other)

    qalv      <- exp(lq_alv)
    qtot      <- exp(lq_tot)
    qh        <- exp(lq_liver)
    qm        <- exp(lq_muscle)
    qo        <- exp(lq_other)

    vp        <- exp(lv_plasma)
    ve        <- exp(lv_erythrocytes)
    vh        <- exp(lv_liver)
    vm        <- exp(lv_muscle)
    vo        <- exp(lv_other)
    vb        <- exp(lv_blood)

    fp        <- fp_blood
    cle       <- exp(lcl_erythrocytes)
    emax      <- emax_erythrocytes
    kaff      <- kaff_erythrocytes

    kscn      <- exp(lkscn)
    kfs       <- exp(lkfs)
    kes       <- exp(lkes)

    # ug/mL per umol/L (27 g/mol / 1000).
    cfUgmL    <- mw_hcn / 1000

    # ------------------------------------------------------------------
    # Exposure driver. The covariate is supplied in ppm; the ODEs work in
    # umol/L. Constant from t = 0 in the published analysis, but a
    # time-varying CONC_HCN_PPM column is honoured.
    # ------------------------------------------------------------------
    ca        <- CONC_HCN_PPM / ppm_per_umol

    # ------------------------------------------------------------------
    # Mixed venous PLASMA concentration: the flow-weighted average of the
    # three tissues' effluent plasma concentrations. Appears both in the
    # plasma mass balance (as the returning venous blood) and in the
    # venous observable.
    # ------------------------------------------------------------------
    cvPlasma  <- ((qh * liver) / php + (qm * muscle) / pmp + (qo * other) / pop) / qtot

    # ------------------------------------------------------------------
    # Erythrocyte binding. The authors write the plasma-equivalent of the
    # erythrocyte concentration as erythrocytes / (emax / (plasma + kaff)),
    # i.e. an effective erythrocyte:plasma partition coefficient
    # emax / (plasma + kaff) that FALLS as plasma cyanide rises - the
    # saturable-binding term. Kept in the authors' algebraic form.
    # ------------------------------------------------------------------
    peBind    <- emax / (plasma + kaff)

    # ------------------------------------------------------------------
    # ODE system, transcribed term-for-term from the Supplemental Material 1
    # `func` body. Each equation is already divided by its compartment
    # volume, so the states are concentrations (umol/L).
    # ------------------------------------------------------------------

    # Plasma: alveolar exchange + erythrocyte exchange + venous return.
    d/dt(plasma) <- (qalv * (ca - plasma / ppa) +
                     cle * (erythrocytes / peBind - plasma) +
                     fp * qtot * (cvPlasma - plasma)) / vp

    # Erythrocytes: exchange with plasma only.
    d/dt(erythrocytes) <- cle * (plasma - erythrocytes / peBind) / ve

    # Liver: perfusion in/out, minus rhodanese consumption, minus the
    # authors' additional -(kfs - kscn*liver*S - kes*S) term. That bracket is
    # exactly -d(sulfur_donor)/dt, so the two kscn*liver*sulfur_donor terms
    # cancel algebraically and the net liver equation reduces to
    #   (fp*qh*(plasma - liver/php) - kfs + kes*sulfur_donor) / vh.
    # This is retained verbatim because it is the published implementation
    # that generated Harada 2025 Fig. 3 and Table 2; the vignette shows the
    # simplified reading is numerically indistinguishable over the published
    # exposure range. See the vignette Errata.
    d/dt(liver) <- (fp * qh * (plasma - liver / php) -
                    kscn * liver * sulfur_donor -
                    (kfs - kscn * liver * sulfur_donor - kes * sulfur_donor)) / vh

    # Sulfane-sulfur donor pool: zero-order supply, consumption by
    # rhodanese, first-order loss. Not volume-divided in the source.
    d/dt(sulfur_donor) <- kfs - kscn * liver * sulfur_donor - kes * sulfur_donor

    # Muscle and lumped other tissue: perfusion-limited.
    d/dt(muscle) <- fp * qm * (plasma - muscle / pmp) / vm
    d/dt(other)  <- fp * qo * (plasma - other / pop) / vo

    # Initial conditions: cyanide-free body with a full sulfur-donor pool.
    # Supplemental Material 1: init = [0, 0, 0, 1, 0, 0] over
    # [Cp, Ce, Ch, S, Cm, Co], i.e. only the donor pool starts non-zero.
    sulfur_donor(0) <- 1

    # ------------------------------------------------------------------
    # Observables (ug/mL). The paper compares these against post-mortem
    # cardiac blood: arterial whole blood proxies LEFT cardiac blood and
    # mixed-venous whole blood proxies RIGHT cardiac blood.
    # ------------------------------------------------------------------
    Carterial <- (vp * plasma + ve * erythrocytes) / vb * cfUgmL
    Cvenous   <- (vp * cvPlasma + ve * erythrocytes) / vb * cfUgmL

    # Diagnostic quantities the source script also computes.
    CvenousPlasma <- cvPlasma                             # mixed venous plasma CN (umol/L)
    metScn        <- kscn * liver * sulfur_donor          # rate of CN -> SCN conversion (umol/L/min)
    avGap         <- Carterial - Cvenous                  # arterial-venous gap (ug/mL); the paper's diagnostic signal
  })
}
