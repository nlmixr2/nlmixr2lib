BlussevanOudAlblas_2019_propofol_bis <- function() {
  description <- paste(
    "Population PK-PD model for propofol in adolescents undergoing idiopathic",
    "scoliosis surgery with an intraoperative wake-up test and reinduction of",
    "anesthesia (Blusse van Oud-Alblas 2019). A two-compartment whole-blood",
    "disposition model (CL 1.37 L/min, V1 3.6 L, Q 1.15 L/min, V2 76.8 L) is",
    "linked to a TWO-compartment biophase (effect-site) distribution model and",
    "an inhibitory sigmoid Emax model for the Bispectral Index (BIS). The",
    "central biophase compartment equilibrates with whole-blood propofol at",
    "ke0 = 0.102 1/min and exchanges with a peripheral biophase compartment at",
    "ke12 = 0.121 and ke21 = 0.172 1/min; the two-compartment biophase gave a",
    "time to peak effect-site concentration of 2.63 min and was superior to a",
    "one-compartment biophase (p < 0.001). BIS baseline and Emax were both",
    "fixed to 100 because too few awake observations were available to estimate",
    "them. No covariate (bodyweight, age, gender or remifentanil infusion rate)",
    "reached significance on any PK or PD parameter. The companion model for",
    "the composite A-line ARX index endpoint of the same study is",
    "modellib('BlussevanOudAlblas_2019_propofol_caai').",
    sep = " "
  )
  reference <- paste(
    "Blusse van Oud-Alblas HJ, Brill MJE, Peeters MYM, Tibboel D, Danhof M,",
    "Knibbe CAJ. (2019). Population pharmacokinetic-pharmacodynamic model of",
    "propofol in adolescents undergoing scoliosis surgery with intraoperative",
    "wake-up test: a study using Bispectral index and composite auditory evoked",
    "potentials as pharmacodynamic endpoints. BMC Anesthesiology 19:15.",
    "doi:10.1186/s12871-019-0684-z. PMCID: PMC6343297.",
    sep = " "
  )
  vignette <- "BlussevanOudAlblas_2019_propofol"
  units <- list(time = "min", dosing = "mg", concentration = "mg/L")
  # Time is minutes because every rate constant in the paper is published in
  # min^-1 and clearances in L/min (Tables 2 and 3). Propofol was assayed in
  # WHOLE BLOOD by HPLC with fluorescence detection (Blood sampling and
  # analysis), so the model concentration is a whole-blood concentration.

  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened as a candidate covariate on the pharmacokinetic and",
        "pharmacodynamic parameters (Covariate analysis: 'Bodyweight, age,",
        "gender, and remifentanil infusion rates were tested as covariates')",
        "but NOT retained: 'None of the explored covariates significantly",
        "influence the pharmacokinetic parameters' (Results, Propofol",
        "pharmacokinetics). The Discussion notes only 'a trend towards an",
        "influence of bodyweight on clearance'. No point estimate is",
        "reported, so nothing is implementable. Cohort median 51 kg",
        "(36.6-82 kg, Table 1)."
      ),
      source_name = "bodyweight"
    ),
    AGE = list(
      description = "Chronological age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened but not retained. For the pharmacodynamics the authors",
        "report 'a negatively linear trend towards an influence of age on the",
        "EC50 (P > 0.05) of both BIS and cAAI' (Results, Propofol",
        "pharmacodynamics) and attribute the non-significance to insufficient",
        "power (Discussion). No point estimate is reported. Cohort median",
        "14.7 years (9.8-20.1 years, Table 1)."
      ),
      source_name = "age"
    ),
    SEXF = list(
      description = "Sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Screened as 'gender' (Covariate analysis) but not retained on any PK",
        "or PD parameter. Cohort was 2 male / 12 female (Table 1)."
      ),
      source_name = "gender"
    ),
    CONMED_REMIFENTANIL = list(
      description = "Concomitant remifentanil infusion rate",
      units = "ug/kg/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened as a candidate covariate ('remifentanil infusion rates were",
        "tested as covariates') but not retained. The Discussion argues that",
        "maximum propofol-remifentanil synergy was already reached at the",
        "0.2-1 ug/kg/min infusion rates used, so no contrast was available.",
        "Documentation only; not a canonical covariate column and not",
        "referenced in model()."
      ),
      source_name = "remifentanil infusion rate"
    )
  )

  compartmentData <- list(
    central = list(analyte = "propofol", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "propofol", units = "mg", specimen = "whole blood", verified = TRUE),
    effect1 = list(analyte = "propofol", units = "mg/L", specimen = "not applicable", verified = TRUE),
    effect2 = list(analyte = "propofol", units = "mg/L", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 14L,
    n_studies = 1L,
    age_range = "9.8-20.1 years",
    age_median = "14.7 years",
    weight_range = "36.6-82 kg",
    weight_median = "51 kg",
    height_range = "142-183 cm",
    height_median = "162.5 cm",
    sex_female_pct = 100 * 12 / 14,
    disease_state = paste(
      "ASA physical status I-II adolescents undergoing surgical correction of",
      "idiopathic scoliosis with an intraoperative wake-up test. Exclusion",
      "criteria were hypacusis or deafness, any neurological disease,",
      "medication affecting the central nervous system, and any",
      "contraindication to the anesthesia protocol."
    ),
    dose_range = paste(
      "Propofol 4 mg/kg IV bolus over 10 s for induction, maintenance infusion",
      "2-10 mg/kg/h, and a 3-5 mg/kg IV reinduction bolus after the",
      "intraoperative wake-up test. Concurrent remifentanil 1 ug/kg/min at",
      "induction then 0.2-1 ug/kg/min for maintenance; rocuronium 0.6 mg/kg;",
      "intrathecal morphine 5 ug/kg and IV morphine 100 ug/kg intraoperatively."
    ),
    regions = "The Netherlands (Erasmus University Medical Center, Rotterdam)",
    n_observations = paste(
      "225 whole-blood propofol samples (median 16, range 6-22 per patient)",
      "and a median of 113 (66-154) BIS observations per patient."
    ),
    notes = paste(
      "Single-centre study approved by the Institutional Ethics Committee of",
      "the Erasmus University Medical Center. Demographics from Table 1,",
      "reported as median (minimum-maximum). Median duration of propofol",
      "administration 410 min (200-460); length of the wake-up test 21.5 min",
      "(7.6-42.4); return of consciousness 52.2 min (22.6-116.33) after the",
      "end of the propofol infusion. Propofol was assayed in whole blood by",
      "HPLC with fluorescence detection, limit of quantification 0.005 mg/L,",
      "inter- and intra-assay CV below 6.7% and 3.3% over 0.5-20 mg/L.",
      "Sequential population PK then PD analysis in NONMEM VI using",
      "first-order conditional estimation with eta-epsilon interaction; the",
      "PD model was fitted to the post hoc PK parameters. Internal validation",
      "by 1000 bootstrap replicates (PK) and 250 bootstrap replicates (PD)",
      "plus normalised prediction distribution errors."
    )
  )

  ini({
    # -----------------------------------------------------------------
    # Propofol whole-blood disposition. Blusse van Oud-Alblas 2019
    # Table 2, 'PK Model Mean (CV%)' column. The parenthesised CV% in
    # Table 2 is the precision of the estimate, not a variability term.
    # NONMEM ADVAN3 TRANS4 fitted to log-transformed concentrations.
    # -----------------------------------------------------------------
    lcl <- log(1.37)
    label("Propofol clearance, CL (L/min)") # Table 2: CL = 1.37 L/min (CV 7.0%)
    lvc <- log(3.6)
    label("Propofol central volume of distribution, V1 (L)") # Table 2: V1 = 3.6 L (CV 10.2%)
    lq <- log(1.15)
    label("Propofol inter-compartmental clearance, Q (L/min)") # Table 2: Q = 1.15 L/min (CV 29.4%)
    lvp <- log(76.8)
    label("Propofol peripheral volume of distribution, V2 (L)") # Table 2: V2 = 76.8 L (CV 5.0%)

    # -----------------------------------------------------------------
    # Two-compartment biophase (effect-site) distribution model for BIS.
    # Table 3, 'BIS Mean (CV%)' column. ke0 is both the rate constant
    # from the central PK compartment into the central effect-site
    # compartment and the rate constant for drug loss from it (Fig. 1
    # caption); ke12 / ke21 connect the central and peripheral effect-site
    # compartments (Fig. 1 caption and the Pharmacodynamic model section).
    # -----------------------------------------------------------------
    lke0 <- log(0.102)
    label("Biophase equilibration rate constant, ke0 (1/min)") # Table 3: ke0 = 0.102 1/min (CV 20.4%)
    lke12 <- log(0.121)
    label("Rate constant from the central to the peripheral biophase compartment, ke12 (1/min)") # Table 3: ke12 = 0.121 1/min (CV 49.3%)
    lke21 <- log(0.172)
    label("Rate constant from the peripheral to the central biophase compartment, ke21 (1/min)") # Table 3: ke21 = 0.172 1/min (CV 35.8%)

    # -----------------------------------------------------------------
    # Inhibitory sigmoid Emax model on BIS, Eq. 4:
    #   PD_ij = PD0 - Emax_i * Ce_ij^gamma / (EC50_i^gamma + Ce_ij^gamma)
    # E0 and Emax were BOTH fixed to 100 because 'because of the
    # relatively small amount of available awake data points in the
    # present study, it was not possible to estimate baseline and Emax
    # for BIS. Therefore, in the final pharmacodynamic model both
    # parameters for BIS were fixed to 100' (Discussion).
    # -----------------------------------------------------------------
    le0 <- fixed(log(100))
    label("Baseline (drug-free) BIS, E0 (BIS units)") # Table 3: E0 = 100 Fixed
    limax <- fixed(log(100))
    label("Maximum reduction in BIS attributable to propofol, Emax (BIS units)") # Table 3: Emax (value) = 100 Fixed
    lec50 <- log(3.51)
    label("Biophase propofol concentration at half the maximum BIS effect, EC50 (mg/L)") # Table 3: EC50 = 3.51 mg/L (CV 6.61%)
    lhill <- log(1.43)
    label("Hill coefficient of the sigmoid Emax model on BIS (unitless)") # Table 3: gamma = 1.43 (CV 10.1%)

    # -----------------------------------------------------------------
    # Inter-individual variability, Eq. 1 (P_i = P_tv * exp(eta_i)),
    # log-normal with mean zero and variance omega^2.
    #
    # Table 2 reports the CL term already back-transformed, as the row
    # 'omega^2 of CL in %' = 22.1: the table footnote defines the printed
    # percentage as sqrt(exp(omega^2) - 1), so the log-scale variance is
    # omega^2 = log(1 + 0.221^2) = 0.0476857.
    #
    # Table 3 reports the PD terms as the variances themselves, which the
    # Results text confirms by back-transforming them: omega^2 = 0.059 is
    # quoted as 'a CV of 25%' (sqrt(exp(0.059) - 1) = 24.7%) and the BIS
    # Hill coefficient variance 0.117 as 'CV = 34%' (35.2%).
    # -----------------------------------------------------------------
    etalcl ~ 0.0476857 # Table 2, row 'omega^2 of CL in %' = 22.1 (CV 28.7%); log(1 + 0.221^2)
    etalec50 ~ 0.059 # Table 3, row 'omega EC50^2' BIS = 0.059 (CV 48.1%); Results quote 'CV of 25%'
    etalhill ~ 0.117 # Table 3, row 'omega gamma^2' BIS = 0.117 (CV 41.0%); Results quote 'CV = 34%'
    etalke0 ~ 0.192 # Table 3, row 'omega keo^2' BIS = 0.192 (CV 49.7%)

    # -----------------------------------------------------------------
    # Residual error.
    #
    # PK: Eq. 2 is additive on the LOG scale (Y_ij = log(cpred_ij) +
    # eps_ij), which the Methods call 'a proportional error model' and
    # which is proportional in nlmixr2's linear space. Table 2 prints the
    # term as a percentage, '19.0%', so the proportional SD is 0.19.
    #
    # PD: Eq. 5 is additive on the BIS scale (Y_ij = PDpred_ij + eps_ij)
    # and Table 3 prints the variance, so the additive SD is sqrt(59.5).
    # -----------------------------------------------------------------
    propSd <- 0.19
    label("Proportional residual error on whole-blood propofol (fraction)") # Table 2: sigma^2 = 19.0% (CV 11.9%)
    addSd_BIS <- 7.713624
    label("Additive residual error on BIS (BIS units)") # Table 3: sigma^2 BIS = 59.5 (CV 11.1%); SD = sqrt(59.5)
  })

  model({
    # 1. Individual parameters, Eq. 1: P_i = P_tv * exp(eta_i).
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)

    ke0 <- exp(lke0 + etalke0)
    ke12 <- exp(lke12)
    ke21 <- exp(lke21)
    ec50 <- exp(lec50 + etalec50)
    hill <- exp(lhill + etalhill)
    e0 <- exp(le0)
    imax <- exp(limax)

    # 2. Micro-constants of the two-compartment disposition model
    #    (NONMEM ADVAN3 TRANS4, parameterised as CL / V1 / Q / V2).
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. Propofol disposition. Propofol is given intravenously (bolus and
    #    infusion) directly into the central compartment; there is no
    #    absorption compartment.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 4. Whole-blood propofol concentration (mg/L).
    Cc <- central / vc

    # 5. Two-compartment biophase distribution model. The equations are not
    #    printed in the paper; they follow from Eq. 3 for the
    #    one-compartment case plus the Fig. 1 caption and the
    #    Pharmacodynamic model section, which state that (a) ke0 carries
    #    propofol from the central PK compartment into the central
    #    effect-site compartment AND is the rate constant for drug loss
    #    from it, and (b) ke12 / ke21 carry propofol from the central
    #    effect-site compartment to the peripheral effect-site compartment
    #    and back. Both states hold a concentration (mg/L), not an amount,
    #    exactly as in Eq. 3. This form reproduces the paper's reported
    #    time to peak effect-site concentration of 2.63 min for a bolus
    #    dose (Results, Propofol pharmacodynamics); the one-compartment
    #    biophase with the same ke0 would give 3.30 min. See the vignette.
    d/dt(effect1) <- ke0 * (Cc - effect1) - ke12 * effect1 + ke21 * effect2
    d/dt(effect2) <- ke12 * effect1 - ke21 * effect2

    # 6. Inhibitory sigmoid Emax model on BIS, Eq. 4, driven by the
    #    concentration in the CENTRAL effect-site compartment.
    BIS <- e0 - imax * effect1^hill / (ec50^hill + effect1^hill)

    Cc ~ prop(propSd)
    BIS ~ add(addSd_BIS)
  })
}
