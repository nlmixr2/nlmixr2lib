Padavia_2024_paracetamol <- function() {
  description <- "Parent-and-metabolites population PK model for intravenous paracetamol (acetaminophen) and its glucuronide, sulphate, and oxidative (cysteine + mercapturate) metabolites in extreme preterm neonates of 23-26 weeks' gestational age receiving prophylactic paracetamol for patent ductus arteriosus (Padavia 2024, TREOCAPA phase II trial). Five-compartment structure: two-compartment plasma disposition for parent paracetamol with first-order elimination, feeding three one-compartment metabolite pools whose distribution volumes are not identifiable and are therefore fixed to the parent central volume. Total parent clearance is split across four parallel routes - glucuronidation, sulphation, oxidation, and unchanged renal elimination - by an odds-style partition on three estimated ratio parameters (the paper's t1, t2, t3), each ratio being that route's clearance divided by the unchanged-renal clearance, so the four resulting fractions sum to one by construction. Body weight enters as a power-law covariate on total clearance, birth length on the peripheral volume, and gestational age both on the glucuronidation partition ratio (which shifts the whole pathway split, raising the glucuronide fraction and lowering the sulphate fraction with increasing maturity) and additionally on the sulphate metabolite's own elimination clearance."
  reference <- paste(
    "Padavia F, Treluyer JM, Cambonie G, Flamant C, Rideau A, Tauzin M,",
    "Patkai J, Gascoin G, Lumia M, Aikio O, Foissac F, Urien S, Benaboud S,",
    "Lui G, Froelicher Bournaud L, Zheng Y, Kemper R, Tortigue M,",
    "Baruteau AE, Kallio J, Hallman M, Diallo A, Levoyer L, Roze JC,",
    "Bouazza N (2024).",
    "Population pharmacokinetics of intravenous paracetamol and its",
    "metabolites in extreme preterm neonates in the context of patent",
    "ductus arteriosus treatment.",
    "Clin Pharmacokinet 63(12):1689-1700.",
    "doi:10.1007/s40262-024-01439-3.",
    sep = " "
  )
  vignette <- "Padavia_2024_paracetamol"
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Methods 2.4 states that doses and concentrations were
  # converted to micromoles using the paracetamol and metabolite molecular
  # weights, so every state is a molar amount; Methods 2.3 confirms plasma as
  # the assayed matrix for parent and all metabolites. Methods 2.4 also
  # defines the oxidative analyte as the molar sum of paracetamol-cysteine
  # and paracetamol-mercapturate, which is what the registered `cysmer`
  # suffix names.
  compartmentData <- list(
    central          = list(analyte = "paracetamol", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1      = list(analyte = "paracetamol", units = "umol", specimen = "plasma", verified = TRUE),
    central_gluc     = list(analyte = "paracetamol glucuronide", units = "umol", specimen = "plasma", verified = TRUE),
    central_sulf     = list(analyte = "paracetamol sulphate", units = "umol", specimen = "plasma", verified = TRUE),
    central_cysmer   = list(analyte = "paracetamol cysteine + mercapturate", units = "umol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight during treatment",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-law scaling on total paracetamol clearance, centred on the population median. Methods 2.2 states that 'Bodyweight and body length were also collected every day during treatment', so this is the time-varying treatment-period weight rather than the time-fixed birth weight; Results 3.2 names the retained covariate 'bodyweight'. The Table 2 footnote centres it on 800 g, which the canonical kilogram units render as 0.8 kg. Because the covariate enters as the ratio (WT / 0.8), the gram-to-kilogram change of units leaves the term numerically identical to the published (W / 800 g) form.",
      source_name        = "W"
    ),
    HT_BIRTH = list(
      description        = "Body length measured at birth",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-law scaling on the peripheral volume of distribution, centred on 33 cm per the Table 2 footnote (the Table 1 median is 32.75 cm; 33 cm is that median rounded, and 33 cm is the value the published equation uses). Time-fixed per subject. Distinct from the daily body length that Methods 2.2 also collected: Results 3.2 and the Conclusion both name the retained covariate specifically as 'birth length', while the weight covariate above is the time-varying one, so the two size descriptors in this model differ in whether they are fixed at birth.",
      source_name        = "BL"
    ),
    GA = list(
      description        = "Gestational age at birth",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-law scaling, centred on 25 weeks per the Table 2 footnote, on two parameters: the glucuronidation partition ratio (the paper's t1) and, additionally, the sulphate metabolite's own elimination clearance. Because t1 sits in the shared denominator of all four metabolisation fractions, the single GA effect on t1 moves every pathway fraction at once - this is the mechanism behind the Results 3.2 statement that 'the glucuronide pathway increased and the sulfate pathway decreased' with increasing age. Time-fixed per subject; the cohort spans only 23-26 completed weeks (Table 1), so extrapolation far outside that window is not supported by the data.",
      source_name        = "GA"
    )
  )

  # Screened in the covariate search of Methods 2.4 but NOT retained in the
  # final model, so they are documentation only and are never referenced in
  # model(). Methods 2.4 lists the tested set as "gestational age, weight,
  # height, sex, 5-min Apgar score, total and conjugated bilirubin level, and
  # blood pressure measurements (systolic and diastolic)"; Results 3.2 reports
  # that after weight, birth length and gestational age were included, "no
  # other covariate was significant". Head circumference, the 5-minute Apgar
  # score and conjugated bilirubin are described in Table 1 and were screened
  # too, but have no canonical column in inst/references/covariate-columns.md
  # and are not carried here rather than minting canonicals for covariates the
  # final model does not use.
  covariatesDataExcluded <- list(
    SEXF = list(
      description        = "Sex, 1 = female",
      units              = "(binary)",
      type               = "binary",
      reference_category = "male",
      notes              = "Screened in the Methods 2.4 stepwise covariate search but not retained. Table 1 reports 17 of 30 subjects male (56.7%). Source column was a male indicator ('Sex (M)'), so the canonical female orientation is SEXF = 1 - SEXM.",
      source_name        = "Sex (M)"
    ),
    TBILI = list(
      description        = "Total bilirubin",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in the Methods 2.4 stepwise covariate search but not retained. Table 1 median 46 umol/L (IQR 38-50.5, range 18-91).",
      source_name        = "Total bilirubin"
    ),
    SBP = list(
      description        = "Systolic arterial blood pressure",
      units              = "mmHg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in the Methods 2.4 stepwise covariate search but not retained. Table 1 median 50.5 mmHg (IQR 44-57, range 39-78).",
      source_name        = "Systolic blood pressure"
    ),
    DBP = list(
      description        = "Diastolic arterial blood pressure",
      units              = "mmHg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in the Methods 2.4 stepwise covariate search but not retained. Table 1 median 29.5 mmHg (IQR 25.25-33, range 13-47); Table 1 separately reports a minimal diastolic pressure median of 19.5 mmHg.",
      source_name        = "Diastolic blood pressure"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 30L,
    n_studies      = 1L,
    age_range      = "23-26 completed weeks' gestational age at birth; dosing began within 12 h of birth and ran 5 days, so postnatal age spans roughly 0-5 days",
    ga_range       = "23-26 weeks (2 subjects at 23, 9 at 24, 7 at 25, 12 at 26)",
    weight_range   = "0.470-0.920 kg birth weight (median 0.800 kg)",
    length_range   = "28.0-36.5 cm birth length (median 32.75 cm)",
    sex_female_pct = 43.3,
    race_ethnicity = NA_character_,
    disease_state  = "Extreme preterm neonates receiving prophylactic intravenous paracetamol for closure of the ductus arteriosus. Exclusions were birth defects or congenital anomalies, twin-to-twin transfusion syndrome, suspected pulmonary hypoplasia, and clinical instability likely to cause rapid death.",
    dose_range     = "Two dose levels of a Bayesian continual-reassessment dose-escalation trial. Level 1 (21 subjects): 20 mg/kg intravenous loading dose then 7.5 mg/kg every 6 h for 5 days. Level 2 (9 subjects): 25 mg/kg loading dose then 10 mg/kg every 6 h for 5 days. Escalation stopped at level 2 because no further efficacy gain was expected. Twenty doses total per subject.",
    regions        = "France and Finland (8 neonatal intensive care units)",
    n_observations = "121 paracetamol and 484 metabolite plasma concentrations; median 4 samples per subject (range 2-6)",
    notes          = "TREOCAPA phase II trial (NCT04459117), enrolled November 2020 to September 2021. Demographics from Table 1; 31 neonates were enrolled and 1 was excluded from the PK analysis for having no paracetamol plasma level. Estimation used Monolix 2021R2 (SAEM); confidence intervals came from 100 bootstrap replicates. Concentrations below the limit of quantification were handled as left-censored data. All mothers received antenatal steroids."
  )

  ini({
    # -----------------------------------------------------------------
    # Parent disposition. Table 2 of Padavia 2024 reports these twice:
    # once under "Parent model estimates" with RSEs and bootstrap
    # percentiles, and again under "Parent-metabolite model estimates"
    # with "-" in both uncertainty columns. Methods 2.4 explains the
    # dash: "The parent-metabolite model was developed by setting the
    # parameters of the previously constructed parent model." They are
    # therefore held fixed in the final parent-metabolite model, which
    # is the model this file encodes, and are wrapped in fixed()
    # accordingly. Because the parent-metabolite model reuses them
    # unchanged, simulating this model and looking only at Cc
    # reproduces the parent model exactly - which is why the paper's
    # two models are one file here.
    # -----------------------------------------------------------------
    lcl              <- fixed(log(0.0785))  ; label("Total paracetamol clearance at the median 0.8 kg body weight (L/h)")            # Table 2 CL = 0.0785 L/h (parent model RSE 8.49%, bootstrap 5th-95th 0.0656-0.0920); carried into the parent-metabolite model
    e_wt_cl          <- fixed(1.81)         ; label("Power exponent of body weight on total paracetamol clearance")                  # Table 2 beta_Cl_W = 1.81 (parent model RSE 0.302%, bootstrap 1.21-3.40); centred on 800 g per the Table 2 footnote
    lvc              <- fixed(log(0.124))   ; label("Paracetamol central volume of distribution (L)")                                # Table 2 V1 = 0.124 L (parent model RSE 1.11%, bootstrap 0.0181-0.575)
    lq               <- fixed(log(3.84))    ; label("Paracetamol inter-compartmental clearance (L/h)")                               # Table 2 Q = 3.84 L/h (parent model RSE 0.386%, bootstrap 0.616-32.4)
    lvp              <- fixed(log(0.672))   ; label("Paracetamol peripheral volume of distribution at the median 33 cm birth length (L)") # Table 2 V2 = 0.672 L (parent model RSE 10.8%, bootstrap 0.246-0.976)
    e_ht_birth_vp    <- fixed(6.46)         ; label("Power exponent of birth length on the peripheral volume of distribution")       # Table 2 beta_V2_BL = 6.46 (parent model RSE 24.8%, bootstrap 1.90-10.6); centred on 33 cm per the Table 2 footnote

    # -----------------------------------------------------------------
    # Pathway partition. Methods 2.4 defines the metabolisation
    # fractions from three estimated parameters t1, t2 and t3 as
    #   f_gluc      = t1 / (t1 + t2 + t3 + 1)
    #   f_sulf      = t2 / (t1 + t2 + t3 + 1)
    #   f_ox        = t3 / (t1 + t2 + t3 + 1)
    #   f_unchanged =  1 / (t1 + t2 + t3 + 1)
    # and the pathway clearances as CL times the matching fraction,
    # with the Table 2 footnote repeating all four. Dividing any
    # pathway fraction by f_unchanged cancels the denominator, so each
    # t is that route's clearance divided by the unchanged-renal
    # clearance - an odds-style ratio against the renal route, not a
    # share of total clearance. That is why these are named
    # lclrat_<metab> rather than being folded into the fm_<pathway>
    # family, whose members are shares of total clearance: an fm_ name
    # here would state the wrong denominator. The shares themselves are
    # derived in model() as fm_gluc / fm_sulf / fm_cysmer, which ARE
    # fm_<pathway> members.
    #
    # The "1" in each denominator is the unchanged-renal route's own
    # ratio against itself and is a structural constant of the
    # parameterisation, not an estimated value.
    # -----------------------------------------------------------------
    lclrat_gluc      <- log(7.53)     ; label("Glucuronidation clearance as a ratio to the unchanged-renal clearance at the median 25 weeks' gestational age (unitless)") # Table 2 t1 = 7.53 (RSE 19.4%, bootstrap 3.30-9.44)
    e_ga_clrat_gluc  <- 6.88          ; label("Power exponent of gestational age on the glucuronidation clearance ratio")            # Table 2 beta_t1_GA = 6.88 (RSE 1.33%, bootstrap 2.64-12.2); centred on 25 weeks per the Table 2 footnote
    lclrat_sulf      <- log(111)      ; label("Sulphation clearance as a ratio to the unchanged-renal clearance (unitless)")                    # Table 2 t2 = 111 (RSE 20.1%, bootstrap 36.8-128); no covariate retained
    lclrat_cysmer    <- log(4.56)     ; label("Oxidation clearance as a ratio to the unchanged-renal clearance (unitless)")                     # Table 2 t3 = 4.56 (RSE 2.13%, bootstrap 2.42-5.77); no covariate retained

    # -----------------------------------------------------------------
    # Metabolite elimination clearances. Table 2 reports these as
    # "CL E gluc", "CL E sulf" and "CL E ox"; the Results 3.2 sentence
    # "Table 2 summarises the specific elimination clearances for each
    # compound" confirms they are each metabolite's own elimination
    # clearance, not a formation clearance (the formation clearances are
    # CL times the fractions above). Gestational age carries an
    # additional effect on the sulphate route only: "An additional
    # effect of gestational age was observed on elimination clearance of
    # sulfate metabolite" (Results 3.2), matching beta_CLE_sulf_GA as
    # the only other GA coefficient in Table 2.
    # -----------------------------------------------------------------
    lcle_gluc        <- log(0.0324)   ; label("Paracetamol glucuronide elimination clearance (L/h)")                                 # Table 2 CL E gluc = 0.0324 L/h (RSE 26.3%, bootstrap 0.0168-0.0823)
    lcle_sulf        <- log(0.0173)   ; label("Paracetamol sulphate elimination clearance at the median 25 weeks' gestational age (L/h)") # Table 2 CL E sulf = 0.0173 L/h (RSE 7.66%, bootstrap 0.0132-0.0217)
    e_ga_cle_sulf    <- 6.26          ; label("Power exponent of gestational age on the sulphate elimination clearance")             # Table 2 beta_CLE_sulf_GA = 6.26 (RSE 1.78%, bootstrap 2.22-10.8); centred on 25 weeks per the Table 2 footnote
    lcle_cysmer      <- log(0.00445)  ; label("Paracetamol cysteine + mercapturate elimination clearance (L/h)")                     # Table 2 CL E ox = 0.00445 L/h (RSE 26.5%, bootstrap 0.00232-0.0162)

    # -----------------------------------------------------------------
    # Between-subject variability. Methods 2.4 gives the exponential
    # model theta_i = theta_pop * exp(eta_i), matching nlmixr2's
    # exp(l<param> + eta) form directly. The Table 1/2 shared footnote
    # states "omega inter-subject variability expressed as standard
    # deviation", so each printed omega is an SD and the variance
    # entered here is that value squared - no CV%-to-variance
    # conversion is involved. omega_CL and omega_V2 are carried from
    # the parent model (Table 2 lists them with "-" uncertainty in the
    # parent-metabolite block) and so are fixed; the five metabolite-
    # side omegas were estimated in the parent-metabolite fit.
    #
    # Table 2 reports no omega for t3, for V1 or for Q, so those three
    # parameters carry no IIV here.
    # -----------------------------------------------------------------
    etalcl            ~ fixed(0.167281)  # Table 2 omega_CL        = 0.409 SD -> 0.409^2; parent-model value (RSE 15.8%, bootstrap 0.318-0.477) carried over
    etalvp            ~ fixed(0.105625)  # Table 2 omega_V2        = 0.325 SD -> 0.325^2; parent-model value (RSE 39.1%, bootstrap 0.0985-0.847) carried over
    etalclrat_gluc    ~ 0.245025         # Table 2 omega_t1        = 0.495 SD -> 0.495^2 (RSE 31.3%, bootstrap 0.183-0.789)
    etalclrat_sulf    ~ 0.853776         # Table 2 omega_t2        = 0.924 SD -> 0.924^2 (RSE 15.7%, bootstrap 0.552-1.27)
    etalcle_gluc      ~ 0.697225         # Table 2 omega_CLE_gluc  = 0.835 SD -> 0.835^2 (RSE 24.2%, bootstrap 0.149-1.15)
    etalcle_sulf      ~ 0.077284         # Table 2 omega_CLE_sulf  = 0.278 SD -> 0.278^2 (RSE 30.2%, bootstrap 0.0915-0.472)
    etalcle_cysmer    ~ 0.857476         # Table 2 omega_CLE_ox    = 0.926 SD -> 0.926^2 (RSE 22.6%, bootstrap 0.194-1.37)

    # -----------------------------------------------------------------
    # Residual error. Methods 2.4 tested additive, proportional and
    # combined models; Results 3.2 reports "Residual variabilities were
    # best described by a proportional error model", so each of the four
    # analytes carries proportional error only. The parent value is the
    # parent model's, carried into the parent-metabolite model (Table 2
    # lists it with "-" uncertainty there) and therefore fixed.
    # -----------------------------------------------------------------
    propSd            <- fixed(0.349)  ; label("Paracetamol proportional residual SD (fraction)")                                   # Table 2 proportional error = 0.349 (parent model RSE 8.18%, bootstrap 0.290-0.381)
    propSd_gluc       <- 0.443         ; label("Paracetamol glucuronide proportional residual SD (fraction)")                       # Table 2 proportional error for glucuronide = 0.443 (RSE 12.5%, bootstrap 0.329-0.623)
    propSd_sulf       <- 0.503         ; label("Paracetamol sulphate proportional residual SD (fraction)")                          # Table 2 proportional error for sulfate = 0.503 (RSE 8.92%, bootstrap 0.420-0.593)
    propSd_cysmer     <- 0.470         ; label("Paracetamol cysteine + mercapturate proportional residual SD (fraction)")           # Table 2 proportional error for oxidative metabolites = 0.47 (RSE 10.9%, bootstrap 0.361-0.641)
  })

  model({
    # Individual parent parameters. Methods 2.4 gives the continuous
    # covariate form theta = theta_pop * (cov / Me(cov))^beta, i.e. a
    # power function of the covariate divided by its population median.
    # The medians are the Table 2 footnote's: 800 g (= 0.8 kg) for
    # weight, 33 cm for birth length, 25 weeks for gestational age.
    cl <- exp(lcl + etalcl) * (WT / 0.8)^e_wt_cl
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp + etalvp) * (HT_BIRTH / 33)^e_ht_birth_vp

    # Pathway clearance ratios against the unchanged-renal route. Only
    # the glucuronidation ratio carries a covariate (gestational age);
    # only it and the sulphation ratio carry IIV.
    clrat_gluc   <- exp(lclrat_gluc + etalclrat_gluc) * (GA / 25)^e_ga_clrat_gluc
    clrat_sulf   <- exp(lclrat_sulf + etalclrat_sulf)
    clrat_cysmer <- exp(lclrat_cysmer)

    # Metabolisation fractions, Methods 2.4. The trailing "+ 1" is the
    # unchanged-renal route's ratio against itself, so the four
    # fractions sum to exactly one and the four pathway clearances sum
    # to exactly the total clearance - the constraint the paper writes
    # as CL = CL*f_gluc + CL*f_sulf + CL*f_ox + CL*f_unchanged.
    partition_denom <- 1 + clrat_gluc + clrat_sulf + clrat_cysmer
    fm_gluc         <- clrat_gluc   / partition_denom
    fm_sulf         <- clrat_sulf   / partition_denom
    fm_cysmer       <- clrat_cysmer / partition_denom
    f_unchanged     <- 1            / partition_denom

    # Pathway clearances: "The clearance associated with each
    # elimination pathway is obtained by multiplying total clearance of
    # paracetamol (CL) by the corresponding metabolisation fraction"
    # (Methods 2.4).
    cl_gluc   <- cl * fm_gluc
    cl_sulf   <- cl * fm_sulf
    cl_cysmer <- cl * fm_cysmer
    cl_renal  <- cl * f_unchanged

    # Metabolite distribution volumes. Methods 2.4: "The metabolite
    # distribution volumes are not identifiable in parent-metabolite
    # models, so they were fixed to the parent volume." The parent
    # volume that the metabolite pools are equated to is the central
    # volume V1, the compartment the metabolites are formed from and
    # sampled in; this is a structural identifiability constraint of
    # the parameterisation rather than an estimated value, so it is
    # written here as an equality rather than as three separate fixed()
    # entries in ini() that a user could inconsistently override.
    vc_gluc   <- vc
    vc_sulf   <- vc
    vc_cysmer <- vc

    # Metabolite elimination clearances; gestational age acts on the
    # sulphate route only.
    cle_gluc   <- exp(lcle_gluc   + etalcle_gluc)
    cle_sulf   <- exp(lcle_sulf   + etalcle_sulf) * (GA / 25)^e_ga_cle_sulf
    cle_cysmer <- exp(lcle_cysmer + etalcle_cysmer)

    # Five-compartment parent-metabolite system of Results 3.2 ("a
    # five-compartment model (two compartments from the parent model and
    # three compartments for the glucuronide, sulfate, and oxidative
    # metabolites)"), drawn in supplemental Figure 1. The four parent
    # elimination arms are written out individually rather than collapsed
    # into a single cl/vc term so that the mass routed to each metabolite
    # pool is visibly the same quantity that leaves the parent; the four
    # arms sum to cl by construction (see partition_denom above).
    d/dt(central)        <- -(cl_gluc + cl_sulf + cl_cysmer + cl_renal) / vc * central -
                             q / vc * central + q / vp * peripheral1
    d/dt(peripheral1)    <-  q / vc * central - q / vp * peripheral1
    d/dt(central_gluc)   <-  cl_gluc   / vc * central - cle_gluc   / vc_gluc   * central_gluc
    d/dt(central_sulf)   <-  cl_sulf   / vc * central - cle_sulf   / vc_sulf   * central_sulf
    d/dt(central_cysmer) <-  cl_cysmer / vc * central - cle_cysmer / vc_cysmer * central_cysmer

    # Plasma concentrations. Methods 2.4 states that doses and
    # concentrations were converted to micromoles per litre using the
    # paracetamol and metabolite molecular weights, so amounts are in
    # umol, volumes in L, and each ratio is umol/L. Users dosing in mg
    # of paracetamol should convert with the paracetamol molecular
    # weight 151.16 g/mol: dose_umol = dose_mg / 0.15116. The molar
    # basis is also what makes the pathway split a mole-for-mole
    # partition of the parent - one micromole of parent leaving by the
    # glucuronidation arm arrives as one micromole of glucuronide.
    Cc        <- central        / vc
    Cc_gluc   <- central_gluc   / vc_gluc
    Cc_sulf   <- central_sulf   / vc_sulf
    Cc_cysmer <- central_cysmer / vc_cysmer

    Cc        ~ prop(propSd)
    Cc_gluc   ~ prop(propSd_gluc)
    Cc_sulf   ~ prop(propSd_sulf)
    Cc_cysmer ~ prop(propSd_cysmer)
  })
}
