Overgaard_2019_semaglutide <- function() {
  description <- "Two-compartment population PK model for subcutaneous semaglutide (GLP-1 receptor agonist) with first-order absorption and first-order elimination, pooled across nine clinical pharmacology trials in healthy volunteers and adults with type 2 diabetes (Overgaard 2019)."
  reference <- "Overgaard RV, Delff PH, Petri KCC, Anderson TW, Flint A, Ingwersen SH. Population pharmacokinetics of semaglutide for type 2 diabetes. Diabetes Therapy. 2019;10(2):649-662. doi:10.1007/s13300-019-0581-y"
  vignette <- "Overgaard_2019_semaglutide"
  units <- list(time = "hour", dosing = "nmol", concentration = "nmol/L")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL and Q (shared exponent 1.01) and on Vc and Vp (shared exponent 0.923) per Overgaard 2019 Table 4. Reference weight 85 kg per Methods (reference subject profile). Time-fixed at baseline.",
      source_name        = "WT"
    ),
    DIAB = list(
      description        = "Type 2 diabetes status (glycaemic status indicator)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normoglycaemia / healthy volunteer)",
      notes              = "Maps the glycaemic-status covariate of Overgaard 2019 (normoglycaemia vs T2D) onto the canonical DIAB binary indicator. 1 = type 2 diabetes; 0 = normoglycaemia (the reference subject profile is healthy). Multiplicative effects: 1.12 on CL and 0.544 on ka per Overgaard 2019 Table 4. The cohort comprised 277 normoglycaemic and 76 T2D subjects (Table 3); only T2D was distinguished, so the canonical DIAB (which does not separate Type 1 vs Type 2) is the appropriate column.",
      source_name        = "T2D"
    )
  )

  population <- list(
    n_subjects     = 353L,                                         # Overgaard 2019 Table 3 (two-compartment model dataset)
    n_studies      = 9L,                                           # Overgaard 2019 Table 1 (nine clinical pharmacology trials)
    age_range      = "19-70 years",                                # Overgaard 2019 Table 3
    age_mean       = "44.6 years (SD 11.8)",                       # Overgaard 2019 Table 3
    weight_range   = "51.9-121.2 kg",                              # Overgaard 2019 Table 3
    weight_mean    = "81.9 kg (SD 15.1)",                          # Overgaard 2019 Table 3
    bmi_range      = "18.7-42.8 kg/m^2",                           # Overgaard 2019 Table 3
    sex_female_pct = 36.0,                                         # 127 / 353 per Overgaard 2019 Table 3
    race_ethnicity = c(White = 92.6, Asian = 4.5, Black = 0.8,
                       OtherMissing = 2.0),                         # Overgaard 2019 Table 3 (327/16/3/7 of 353)
    diabetes_pct   = 21.5,                                         # 76 T2D / 353 total per Overgaard 2019 Table 3
    disease_state  = "Pooled clinical pharmacology cohort: healthy volunteers (277) and adults with type 2 diabetes (76); a hepatic-impairment subgroup of 25 subjects was also included.",
    dose_range     = "0.25 to 1.5 mg semaglutide once weekly subcutaneous (and 0.25 mg single intravenous in trial 7) across nine clinical pharmacology trials.",
    regions        = "Multinational; included a Japanese safety-and-PK trial (trial 1) alongside Caucasian-majority trials.",
    trials         = c("NCT02146079", "NCT02110871", "NCT02212067", "NCT02064348",
                       "NCT02147431", "NCT02079870", "NCT02231684", "NCT02022254",
                       "NCT02243098"),
    notes          = "Demographics from Overgaard 2019 Table 3 (two-compartment model column). The reference subject profile used for the typical-value parameters is a healthy, white, non-Hispanic female aged <= 65 years with body weight 85 kg, abdomen injection site, and 1.34 mg/mL drug product strength (Methods, page 654). Injection-site (abdomen / thigh) and drug-product-strength (1, 1.34, 3, 10 mg/mL) covariate effects are documented in the vignette but not encoded in this model file; see the vignette's Assumptions and deviations section for rationale."
  )

  ini({
    # Structural parameters from Overgaard 2019 Table 4 (final two-compartment model). Reference
    # subject: healthy (DIAB = 0), 85 kg body weight, abdomen injection site, 1.34 mg/mL product.
    lka     <- log(0.0253);  label("First-order SC absorption rate constant (1/hour)")             # Overgaard 2019 Table 4 (ka = 0.0253 1/h, 1.34 mg/mL)
    lcl     <- log(0.0348);  label("Clearance at reference covariates (L/hour)")                    # Overgaard 2019 Table 4 (CL = 0.0348 L/h, 85 kg, healthy)
    lvc     <- log(3.59);    label("Central volume of distribution at reference body weight (L)")   # Overgaard 2019 Table 4 (Vc = 3.59 L, 85 kg)
    lq      <- log(0.304);   label("Intercompartmental clearance at reference body weight (L/hour)") # Overgaard 2019 Table 4 (Q = 0.304 L/h, 85 kg)
    lvp     <- log(4.10);    label("Peripheral volume of distribution at reference body weight (L)") # Overgaard 2019 Table 4 (Vp = 4.10 L, 85 kg)
    lfdepot <- log(0.847);   label("Absolute SC bioavailability at abdomen, 1.34 mg/mL (fraction)")  # Overgaard 2019 Table 4 (F = 0.847)

    # Covariate-effect parameters per Overgaard 2019 Table 4, Methods page 653 (covariate model
    # specification), and page 654 (covariate functional forms).
    e_wt_cl_q  <- 1.01;      label("Power exponent of body weight on CL and Q (shared)")           # Overgaard 2019 Table 4 (BW on CL/Q = 1.01)
    e_wt_vc_vp <- 0.923;     label("Power exponent of body weight on Vc and Vp (shared)")          # Overgaard 2019 Table 4 (BW on Vc/Vp = 0.923)
    e_diab_cl  <- 1.12;      label("T2D-vs-healthy multiplier on CL (1.12^DIAB)")                  # Overgaard 2019 Table 4 (T2D on CL = 1.12)
    e_diab_ka  <- 0.544;     label("T2D-vs-healthy multiplier on ka (0.544^DIAB)")                 # Overgaard 2019 Table 4 (T2D on ka = 0.544)

    # IIV from Overgaard 2019 Table 4 (CV%): omega^2 = log(1 + CV^2). Single eta on absorption (ka).
    # Correlated block on (CL, Vc) with shared etas: eta_CL also drives Q, eta_V also drives Vp
    # (Methods page 653, "Inter-individual variability ... using a log normal distribution in a model
    # that included covariates with clear effects"). Covariance reported in Table 4 footnote.
    #   IIV CL  = 15.2% -> omega^2 = log(1 + 0.152^2) = 0.022840
    #   IIV Vc  = 15.4% -> omega^2 = log(1 + 0.154^2) = 0.023438
    #   cov(CL, Vc) = 0.0172 (Overgaard 2019 Table 4 footnote)
    #   IIV ka  = 37.9% -> omega^2 = log(1 + 0.379^2) = 0.133844
    etalcl + etalvc ~ c(0.022840,
                        0.0172, 0.023438)                                                             # Overgaard 2019 Table 4 (IIV CL 15.2%, IIV Vc 15.4%, cov 0.0172)
    etalka ~ 0.133844                                                                                  # Overgaard 2019 Table 4 (IIV ka 37.9%)

    # Residual error: additive on log-transformed observations (Overgaard 2019 Methods page 652
    # "An additive residual error model on log-transformed concentration values was assessed to
    # adequately describe the data"); maps to proportional in nlmixr2 linear space.
    propSd <- 0.103; label("Proportional residual error (SD, fraction)")                              # Overgaard 2019 Table 4 (residual error = 0.103, additive on log scale)
  })

  model({
    # Individual PK parameters per Overgaard 2019 Methods page 654 (covariate equations):
    #   CL_i = CL_ref * E_weight,CL * E_glycaemicstatus * exp(eta_CL)
    #   Q_i  = Q_ref  * E_weight,CL                     * exp(eta_CL)   [shares eta_CL]
    #   Vc_i = Vc_ref * E_weight,V                       * exp(eta_V)
    #   Vp_i = Vp_ref * E_weight,V                       * exp(eta_V)   [shares eta_V]
    #   ka_i = ka_ref * theta_ka,T2D^T2D                 * exp(eta_ka)  [for 1.34 mg/mL only]
    # All non-T2D / non-WT covariates (sex, race, ethnicity, age group, hepatic impairment,
    # injection site, drug product strength) lie within the 80-125% equivalence interval and are
    # excluded from this model; see vignette Assumptions and deviations for the omitted
    # injection-site (thigh) and drug-product-strength (1, 3, 10 mg/mL) terms.
    cl <- exp(lcl + etalcl) * (WT / 85)^e_wt_cl_q  * e_diab_cl^DIAB
    q  <- exp(lq  + etalcl) * (WT / 85)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 85)^e_wt_vc_vp
    vp <- exp(lvp + etalvc) * (WT / 85)^e_wt_vc_vp
    ka <- exp(lka + etalka) * e_diab_ka^DIAB

    # Micro-constants for the two-compartment ODE system
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment disposition with first-order SC absorption (Overgaard 2019 Methods page 653,
    # Model Description: "first-order absorption, two-compartment disposition and first-order
    # elimination, whereas intravenous data were modelled as a bolus injection into the central
    # compartment"). Bioavailability F applies to SC dosing into the depot; IV data dose directly
    # into central via the data set's cmt column.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    f(depot) <- exp(lfdepot)

    # Concentration: dose in nmol, vc in L -> central / vc in nmol/L. Semaglutide MW = 4113.6 g/mol;
    # mg-scale doses convert via 1 mg = 1000 / 4113.6 = 0.2431 umol = 243.1 nmol.
    Cc <- central / vc

    # Additive on log scale ~ proportional in linear space
    Cc ~ prop(propSd)
  })
}
