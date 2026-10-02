He_2020_phosphocreatine <- function() {
  description <- "Joint parent-metabolite population PK model for intravenous phosphocreatine (PCr) and its metabolite creatine (Cr) in children with acute myocarditis (He 2020). Four-compartment chain: two-compartment disposition for PCr (central + peripheral), of which a fixed fraction Fm = 0.75 of the PCr elimination flux forms Cr, and two-compartment disposition for Cr, with first-order elimination from both central compartments. Observed Cr is the exogenous (PCr-derived) Cr concentration plus an estimated constant endogenous baseline (66.6 umol/L). Body weight scales every clearance (exponent 0.75, fixed) and every volume (exponent 1, fixed), referenced to 20 kg; bedside-Schwartz eGFR enters Cr clearance as a power function (exponent 0.311, reference 127.78 mL/min/1.73 m^2). Amounts are in umol and time in minutes; doses in grams of phosphocreatine sodium must be converted to umol of PCr before use."
  reference <- paste(
    "He H, Zhang M, Zhao LB, Sun N, Zhang Y, Yuan Y, Wang XL.",
    "Population Pharmacokinetics of Phosphocreatine and Its Metabolite Creatine in Children With Myocarditis.",
    "Front Pharmacol. 2020;11:574141. doi:10.3389/fphar.2020.574141.",
    sep = " "
  )
  vignette <- "He_2020_phosphocreatine"

  units <- list(time = "min", dosing = "umol", concentration = "umol/L")

  covariateData <- list(
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling of every clearance (CL_PCr, Q_PCr, CL_Cr, Q_Cr; exponent 0.75) and every volume (Vc_PCr, Vp_PCr, Vc_Cr, Vp_Cr; exponent 1), both exponents fixed, referenced to 20 kg (Methods Eqs. 3-4: 'theta_CL and theta_V are the respective parameter values for a subject with a bodyweight of 20 kg'; Eq. 17 and Discussion restate the 20 kg reference). Cohort median 20.4 kg, range 7.9-86 kg (Table 1).",
      source_name = "BW"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate, BSA-normalized (creatinine-based paediatric Schwartz equation).",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (CRCL / 127.78)^0.311 on creatine clearance only (Results Eq. 17). Computed with the bedside Schwartz equation 0.413 x height (cm) / serum creatinine (mg/dL) for children aged 1-18 years (Methods Eq. 1) and the original Schwartz equation 0.45 x height (cm) / serum creatinine (mg/dL) for children under 1 year (Methods Eq. 2). Cohort median 127.7805, range 66.33-224.01 mL/min/1.73 m^2 (Table 1); children with renal insufficiency were excluded.",
      source_name = "GFR"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "A sigmoid age-maturation factor MF = AGE^Hill / (k50^Hill + AGE^Hill) on clearance (Methods Eq. 5) was tested but made the model unstable with unreasonable estimates, so it was not retained (Results, Population Pharmacokinetic Model). Age was also screened in the stepwise covariate search and not retained."
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "phosphocreatine",
      units = "umol",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "phosphocreatine",
      units = "umol",
      specimen = "plasma",
      verified = TRUE
    ),
    central_creatine = list(
      analyte = "creatine (exogenous, formed from phosphocreatine; endogenous baseline added algebraically)",
      units = "umol",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_creatine = list(
      analyte = "creatine (exogenous, formed from phosphocreatine)",
      units = "umol",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 100L,
    n_studies = 1L,
    age_range = "0.38-16.45 years (median 5.78 years)",
    weight_range = "7.9-86 kg (median 20.4 kg)",
    sex_female_pct = 44,
    race_ethnicity = "Chinese (single centre in Beijing; race not tabulated)",
    disease_state = "Children (under 18 years) with clinically diagnosed acute-stage myocarditis (onset within about half a year); renal insufficiency, prior heart disease and other causes of myocardial injury were excluded.",
    dose_range = "Single 30 +/- 2 min IV infusion of phosphocreatine sodium by age band: 0.5 g (28 days to under 1 year), 1 g (1 to under 6 years), 2 g (6 to under 18 years).",
    regions = "China (Beijing Children's Hospital, Capital Medical University)",
    renal_function = "Bedside-Schwartz eGFR median 127.78 mL/min/1.73 m^2 (range 66.33-224.01)",
    n_observations = "997 plasma concentrations (498 PCr, 499 Cr); 48.7% of PCr concentrations were below the LLOQ of 1.96 umol/L and were handled with the M3 method; all Cr concentrations were above the LLOQ of 30.53 umol/L.",
    notes = "56 males and 44 females. Samples at baseline and approximately 30, 40 or 50, 75 and 180 min after the start of infusion. Fitted in Phoenix NLME 8.2 with FOCE-ELS using a sequential approach (PCr parameters estimated first, then fixed while the Cr parameters were estimated)."
  )

  ini({
    # ---------------------------------------------------------------------
    # Phosphocreatine (parent) -- Table 2, final-model estimates. Typical
    # values for a 20 kg child (Methods Eqs. 3-4). Time unit is minutes.
    # ---------------------------------------------------------------------
    lvc <- log(8.22); label("PCr central volume Vc_PCr at 20 kg (L)")                         # Table 2: VcPCr = 8.22 L (RSE 10.1%)
    lvp <- log(3.07); label("PCr peripheral volume Vp_PCr at 20 kg (L)")                      # Table 2: VpPCr = 3.07 L (RSE 10.9%)
    lcl <- log(1.33); label("PCr clearance CL_PCr at 20 kg (L/min)")                          # Table 2: CLPCr = 1.33 L/min (RSE 3.09%)
    lq  <- log(0.136); label("PCr intercompartmental clearance Q_PCr at 20 kg (L/min)")      # Table 2: QPCr = 0.136 L/min (RSE 30%)

    # Fraction of PCr elimination that forms Cr, fixed (not identifiable
    # alongside the Cr volumes) from the approximately three-quarters
    # conversion reported in animals (Xu et al., 2014).
    fm <- fixed(0.75); label("Fraction of PCr clearance converted to Cr, Fm (fraction)")     # Methods and Results text: Fm fixed to 0.75

    # ---------------------------------------------------------------------
    # Creatine (metabolite) -- Table 2.
    # ---------------------------------------------------------------------
    lvc_creatine <- log(2.39);   label("Cr central volume Vc_Cr at 20 kg (L)")                                   # Table 2: VcCr = 2.39 L (RSE 8.17%)
    lvp_creatine <- log(2.9);    label("Cr peripheral volume Vp_Cr at 20 kg (L)")                                # Table 2: VpCr = 2.9 L (RSE 4.07%)
    lcl_creatine <- log(0.0825); label("Cr clearance CL_Cr at 20 kg and eGFR 127.78 mL/min/1.73 m^2 (L/min)")   # Table 2: CLCr = 0.0825 L/min (RSE 1.64%); Eq. 17
    lq_creatine  <- log(0.146);  label("Cr intercompartmental clearance Q_Cr at 20 kg (L/min)")                  # Table 2: QCr = 0.146 L/min (RSE 9.87%)
    lrbase_creatine <- log(66.6); label("Endogenous baseline plasma Cr concentration baseCr (umol/L)")          # Table 2: baseCr = 66.6 umol/L (RSE 3.53%)

    # ---------------------------------------------------------------------
    # Covariate effects
    # ---------------------------------------------------------------------
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent shared by all clearances (unitless)")   # Methods: power exponents fixed at 0.75 for clearances (Eq. 3)
    e_wt_vc_vp <- fixed(1);    label("Allometric exponent shared by all volumes (unitless)")      # Methods: power exponents fixed at 1.0 for volumes (Eq. 4)
    e_crcl_cl_creatine <- 0.311; label("Exponent of (eGFR/127.78) on Cr clearance (unitless)")    # Table 2: GFR on CLCr = 0.311 (RSE 30%); Eq. 17

    # ---------------------------------------------------------------------
    # Inter-individual variability -- Table 2, omega^2 (log-normal,
    # exponential IIV, Methods Eq. 6). The printed CV% equals
    # 100*sqrt(omega^2) (e.g. 0.0378 -> 19.4%). No covariances reported.
    # ---------------------------------------------------------------------
    etalcl             ~ 0.0378   # Table 2: eta CLPCr omega^2 = 0.0378 (CV 19.4%)
    etalvc_creatine    ~ 0.0882   # Table 2: eta VcCr omega^2 = 0.0882 (CV 29.7%)
    etalvp_creatine    ~ 0.0354   # Table 2: eta VpCr omega^2 = 0.0354 (CV 18.8%)
    etalcl_creatine    ~ 0.0233   # Table 2: eta CLCr omega^2 = 0.0233 (CV 15.3%)
    etalrbase_creatine ~ 0.121    # Table 2: eta baseCr omega^2 = 0.121 (CV 34.8%)

    # ---------------------------------------------------------------------
    # Residual variability -- Table 2 footnote b / Results: proportional
    # error for both analytes (Methods Eq. 7, C = IPRED x (1 + eps)).
    # Phoenix reports the residual standard deviation.
    # ---------------------------------------------------------------------
    propSd          <- 0.244;  label("PCr proportional residual SD (fraction)")   # Table 2: sigma PCr = 0.244 (RSE 6.79%)
    propSd_creatine <- 0.0519; label("Cr proportional residual SD (fraction)")    # Table 2: sigma Cr = 0.0519 (RSE 8.25%)
  })

  model({
    # ----- Reference covariate values
    ref_wt   <- 20      # kg; Methods Eqs. 3-4 and Eq. 17
    ref_crcl <- 127.78  # mL/min/1.73 m^2; Eq. 17 (cohort median, Table 1)

    allom_cl <- (WT / ref_wt)^e_wt_cl_q
    allom_v  <- (WT / ref_wt)^e_wt_vc_vp

    # ----- Individual parameters: phosphocreatine (parent)
    cl <- exp(lcl + etalcl) * allom_cl
    vc <- exp(lvc) * allom_v
    q  <- exp(lq)  * allom_cl
    vp <- exp(lvp) * allom_v

    # ----- Individual parameters: creatine (metabolite)
    # Eq. 17: CL_Cr = 0.0825 x (BW/20)^0.75 x (GFR/127.78)^0.311
    cl_creatine <- exp(lcl_creatine + etalcl_creatine) * allom_cl *
      (CRCL / ref_crcl)^e_crcl_cl_creatine
    vc_creatine <- exp(lvc_creatine + etalvc_creatine) * allom_v
    q_creatine  <- exp(lq_creatine) * allom_cl
    vp_creatine <- exp(lvp_creatine + etalvp_creatine) * allom_v
    rbase_creatine <- exp(lrbase_creatine + etalrbase_creatine)

    # ----- Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    kel_creatine <- cl_creatine / vc_creatine
    k12_creatine <- q_creatine  / vc_creatine
    k21_creatine <- q_creatine  / vp_creatine

    # ----- ODE system (Results Eqs. 11-14). The zero-order infusion rate
    # A0 enters the PCr central compartment through the dosing record. A
    # fraction Fm of the PCr elimination flux (CL_PCr / Vc_PCr x A1) forms
    # Cr mole-for-mole; the remaining (1 - Fm) leaves the system.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    d/dt(central_creatine)     <- fm * kel * central - kel_creatine * central_creatine -
      k12_creatine * central_creatine + k21_creatine * peripheral1_creatine
    d/dt(peripheral1_creatine) <- k12_creatine * central_creatine - k21_creatine * peripheral1_creatine

    # ----- Observations (Results Eqs. 15-16). Observed Cr is the
    # PCr-derived Cr concentration plus the endogenous baseline baseCr;
    # predose PCr was below the LLOQ in every subject, so PCr baseline is 0.
    Cc          <- central / vc
    Cc_creatine <- central_creatine / vc_creatine + rbase_creatine

    Cc          ~ prop(propSd)
    Cc_creatine ~ prop(propSd_creatine)
  })
}
