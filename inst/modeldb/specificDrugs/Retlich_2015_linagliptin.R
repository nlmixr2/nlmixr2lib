Retlich_2015_linagliptin <- function() {
  description <- "Two-compartment population PK model with concentration-dependent (saturable) binding of linagliptin to dipeptidyl peptidase-4 in both central and peripheral compartments, coupled with a population sigmoid Emax PK/PD model relating total linagliptin plasma concentration to plasma DPP-4 activity, in adults with type 2 diabetes mellitus (Retlich 2015 Tables 4 and 5)."
  reference <- paste(
    "Retlich S, Duval V, Graefe-Mody U, Friedrich C, Patel S, Jaehde U, Staab A.",
    "Population Pharmacokinetics and Pharmacodynamics of Linagliptin in Patients",
    "with Type 2 Diabetes Mellitus. Clin Pharmacokinet. 2015;54(7):737-750.",
    "doi:10.1007/s40262-014-0232-4.",
    "Distribution / binding parameters (V_P/F, Q/F, Kd, Amax,P/F) are fixed from the",
    "upstream popPK model in Retlich S, Duval V, Graefe-Mody U, Jaehde U, Staab A.",
    "J Clin Pharmacol. 2010;50(8):873-885."
  )
  vignette <- "Retlich_2015_linagliptin"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL (linagliptin total plasma concentration; converted from nmol/L via MW 472.54 g/mol); RFU (DPP-4 activity)")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline weight; linear-deviation effect on relative bioavailability F centred at 88 kg (population median, Retlich 2015 Table 4 footnote b).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-deviation effect on Bmax_C centred at 60 years (population median, Retlich 2015 Table 4 footnote g).",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Sex (female indicator)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Additive effect on Bmax_C (+0.457 nmol/L for females, Retlich 2015 Table 4) and on baseline DPP-4 activity BSL (+865 RFU for females, Retlich 2015 Table 5).",
      source_name        = "SEX (1 = female; 0 = male)"
    ),
    DOSE = list(
      description        = "Administered linagliptin dose at the current dosing record",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-record dose value (linagliptin doses studied: 0.5, 1, 2.5, 5, 10 mg). Used as a covariate centred at 5 mg in the empirical dose-dependent Ka and Bmax_C effects (Retlich 2015 Table 4 footnotes c and g). Pass as a constant per subject for chronic once-daily simulations.",
      source_name        = "DOSE"
    ),
    GGT = list(
      description        = "Serum gamma-glutamyltransferase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-deviation effect on linagliptin CL centred at 33 U/L (Retlich 2015 Table 4 footnote f, popPK median). The PD layer uses GGT centred at 32.3 U/L (Retlich 2015 Table 5 footnote b, popPK/PD median) with a piecewise threshold at 175 U/L.",
      source_name        = "GGT"
    ),
    DPP4_BL_RFU = list(
      description        = "Baseline plasma DPP-4 enzymatic activity in relative fluorescence units",
      units              = "RFU",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-subject baseline DPP-4 activity from the linagliptin Boehringer Ingelheim assay. Linear-deviation effect on Bmax_C centred at 12,497 RFU (Retlich 2015 Table 4 footnote g, popPK median). The popPK/PD layer estimates BSL as a per-subject parameter; the BSL_EC50 covariate effect is centred at 11,600 RFU (Retlich 2015 Table 5 footnote c) and acts on the individually predicted BSL_i.",
      source_name        = "DPP (baseline RFU)"
    ),
    ALT = list(
      description        = "Serum alanine aminotransferase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-deviation effect on baseline DPP-4 activity BSL centred at 28.8 U/L (Retlich 2015 Table 5 footnote b, popPK/PD median).",
      source_name        = "ALT"
    ),
    FPG = list(
      description        = "Fasting plasma glucose at baseline",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-deviation effect on baseline DPP-4 activity BSL centred at 8.90 mmol/L (Retlich 2015 Table 5 footnote b, popPK/PD median).",
      source_name        = "FPG"
    ),
    TRIG = list(
      description        = "Serum triglyceride concentration",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-deviation effect on baseline DPP-4 activity BSL and on EC50, both centred at 160 mg/dL (Retlich 2015 Table 5 footnotes b and c, popPK/PD median).",
      source_name        = "TG"
    ),
    TCHOL = list(
      description        = "Total serum cholesterol",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-deviation effect on baseline DPP-4 activity BSL centred at 183 mg/dL (Retlich 2015 Table 5 footnote b, popPK/PD median).",
      source_name        = "CHOL"
    ),
    CONMED_METFORMIN = list(
      description        = "Concomitant metformin co-administration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (linagliptin monotherapy)",
      notes              = "1 = Retlich 2015 Study 4 add-on-to-metformin cohort. Multiplicative effect on relative bioavailability F (+69% F for metformin co-administration, Retlich 2015 Table 4 row 'F in study 4'); attributed by the authors to a metformin-linagliptin DDI consistent with Graefe-Mody 2009.",
      source_name        = "METFORMIN"
    ),
    FORM_POWDER = list(
      description        = "Linagliptin powder-in-bottle formulation indicator (Study 1)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (tablet)",
      notes              = "1 = Retlich 2015 Study 1 powder-in-bottle formulation; 0 = tablet (formulation 1 or 2). Multiplicative effect on the absorption rate constant Ka (typical Ka for powder = 0.933 1/h vs 0.441 1/h for tablet formulation 2 = the marketed-product reference, Retlich 2015 Table 4).",
      source_name        = "FORMPOW"
    ),
    FORM_LINAG_TAB1 = list(
      description        = "Linagliptin tablet formulation 1 indicator (Study 2)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (tablet formulation 2 = marketed reference; or powder if FORM_POWDER = 1)",
      notes              = "1 = Retlich 2015 Study 2 tablet formulation 1 (development formulation, not marketed); 0 = tablet formulation 2 (the marketed linagliptin tablet) OR powder. Multiplicative effect on Ka (typical Ka for tablet formulation 1 = 0.795 1/h vs 0.441 1/h for tablet formulation 2 reference, Retlich 2015 Table 4).",
      source_name        = "FORMTAB1"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 607L,
    n_subjects_pk  = 462L,
    n_subjects_pkpd = 607L,
    n_studies      = 4L,
    studies        = "Studies 1-2: phase 1 (single + multiple oral dose, powder-in-bottle [Study 1] and tablet formulation 1 [Study 2], 12 days / 4 weeks). Studies 3-4: phase 2b (12 weeks; Study 3 linagliptin monotherapy, Study 4 add-on to metformin; both used tablet formulation 2 = marketed product).",
    age_range      = "30-78 years (pooled, Retlich 2015 Table 3)",
    age_median     = "60 years",
    weight_range   = "55-138 kg (pooled, Retlich 2015 Table 3)",
    weight_median  = "89 kg",
    bmi_median     = "30.6 kg/m^2 (range 20.4-42.2)",
    sex_female_pct = 33.9,
    race_ethnicity = c(Caucasian = 92.1, Black = 2.5, Asian = 1.8, Hispanic = 3.6),
    disease_state  = "Type 2 diabetes mellitus (T2DM). Patients with normal hepatic function and normal renal function or mild renal impairment.",
    dose_range     = "0.5-10 mg PO once daily linagliptin (Studies 3 and 4 marketed-dose range bracketed by 0.5-10 mg in phase 2b).",
    regions        = "Multiregional (Retlich 2015 Methods); ethnic origin coded as Caucasian / Black / Asian / Hispanic.",
    co_medication  = "Metformin in Study 4 (267/607 patients, 44%); linagliptin monotherapy in Studies 1, 2, 3.",
    n_observations = "PK dataset: 6907 linagliptin plasma concentrations (462 subjects). PK/PD dataset: 9674 paired linagliptin + DPP-4 activity measurements (607 subjects including placebo arms).",
    baseline_fpg_median = "9.9 mmol/L (range 5.1-20.0)",
    notes          = "Demographics summarised in Retlich 2015 Table 3."
  )

  ini({
    # ---- Structural PK parameters (Retlich 2015 Table 4) ----
    # The popPK structural parameters carried forward unchanged from the upstream
    # Retlich 2010 J Clin Pharmacol popPK model (V_P/F, Q/F, Kd, Amax,P/F) are
    # wrapped in fixed() per the Table 4 footnote d ("Parameters not estimated,
    # but fixed to estimates of the previous model").
    lka <- log(0.441);   label("Typical Ka for tablet formulation 2 (marketed linagliptin tablet); 1/h")  # Table 4 row Ka,3
    lvc      <- log(715);     label("Apparent central volume of distribution VC/F (L)")                       # Table 4 row VC/F
    lvp      <- fixed(log(1650));  label("Apparent peripheral volume of distribution VP/F (L; fixed from upstream Retlich 2010)") # Table 4 row VP/F, footnote d
    lq       <- fixed(log(412));   label("Apparent inter-compartmental clearance QP/F (L/h; fixed from upstream Retlich 2010)")    # Table 4 row QP/F, footnote d
    lcl      <- log(258);     label("Apparent clearance of the unbound linagliptin concentration CL/F (L/h)") # Table 4 row CL/F
    lbmaxc  <- log(4.97);    label("Typical Bmax_C: central-compartment DPP-4 binding-site concentration in males (nmol/L)")  # Table 4 row Bmax,C
    lamax_p  <- fixed(log(1650));  label("Apparent peripheral binding-site amount Amax_P/F (nmol; fixed from upstream Retlich 2010)") # Table 4 row Amax,P/F
    lkd      <- fixed(log(0.0652));label("DPP-4 dissociation constant Kd (nmol/L; fixed from upstream Retlich 2010)")                # Table 4 row Kd
    lfdepot  <- fixed(log(1));     label("Reference relative bioavailability (F = 1, fixed at the typical-value reference)")        # Table 4 row F, footnote a

    # ---- PK covariate effects (Retlich 2015 Table 4) ----
    e_pow_ka       <- fixed(0.7493);  label("Powder-formulation effect on log-Ka (log(0.933/0.441))")  # Table 4 rows Ka,1 vs Ka,3 (derived ratio)
    e_tab1_ka      <- fixed(0.5896);  label("Tablet-formulation-1 effect on log-Ka (log(0.795/0.441))") # Table 4 rows Ka,2 vs Ka,3 (derived ratio)
    e_dose_ka      <- -0.0651;  label("Fractional change in Ka per mg DOSE deviation from 5 mg (= -6.51 % per mg)")  # Table 4 row Dose_Ka
    e_wt_f         <- -0.00958; label("Fractional change in F per kg WT deviation from 88 kg (= -0.958 % per kg)")    # Table 4 row Weight_F
    e_metformin_f  <- 0.69;     label("Fractional increase in F for concomitant metformin co-administration (= +69 %)") # Table 4 row 'F in study 4' (169 % relative to 100 % reference)
    e_ggt_cl       <- -0.000339;label("Fractional change in CL per U/L GGT deviation from 33 U/L (= -0.0339 % per U/L)") # Table 4 row GGT_CL
    e_dpp4_bmaxc   <- 3.32e-5;  label("Fractional change in Bmax_C per RFU DPP4_BL_RFU deviation from 12,497 RFU (= 0.00332 % per RFU)") # Table 4 row DPP_Bmax,C
    e_dose_bmaxc   <- 0.0341;   label("Fractional change in Bmax_C per mg DOSE deviation from 5 mg (= 3.41 % per mg)")  # Table 4 row Dose_Bmax,C
    e_age_bmaxc    <- 0.00561;  label("Fractional change in Bmax_C per year AGE deviation from 60 years (= 0.561 % per year)") # Table 4 row Age_Bmax,C
    e_sex_bmaxc    <- 0.457;    label("Additive shift in Bmax_C for SEXF = 1 (nmol/L, female - male)")  # Table 4 row Sex_Bmax,C

    # ---- PK inter-individual variability (Retlich 2015 Table 4; CV % converted to log-normal variance via omega^2 = log(1 + CV^2)) ----
    etalfdepot + etalcl ~ c(0.20271,
                            -0.09300, 0.07291)  # omega^2_F = log(1 + 0.474^2) = 0.20271 (CV 47.4%); omega^2_CL = log(1 + 0.275^2) = 0.07291 (CV 27.5%); rho(F,CL) = -0.765, cov = -0.765 * sqrt(0.20271 * 0.07291) = -0.09300
    etalka    ~ 0.46019  # omega^2 = log(1 + 0.764^2) = 0.46019 (CV 76.4%)   Table 4 row xKa
    etalvc    ~ 0.05772  # omega^2 = log(1 + 0.244^2) = 0.05772 (CV 24.4%)   Table 4 row xVC
    etalbmaxc ~ 0.02227  # omega^2 = log(1 + 0.150^2) = 0.02227 (CV 15.0%)   Table 4 row xBmax,C

    # ---- PK residual error (Retlich 2015 Table 4 row r_prop,phase 2a; footnote h: coded as additive on log-transformed data) ----
    # For SD = 0.136 the additive-on-log error is approximately equivalent to a proportional CV of 13.6 % in linear space.
    # The paper reports a higher phase-2b residual SD (0.383); only the phase-1 (controlled-clinical-pharmacology) value is encoded here.
    # The phase-2b value is documented in the vignette's Assumptions section as a deviation.
    propSd <- 0.136;  label("Proportional residual error for linagliptin plasma concentration (fraction; phase 1 value)")  # Table 4 row r_prop,phase 2a

    # ---- PD structural parameters (Retlich 2015 Table 5) ----
    lbsl  <- log(10700);    label("Typical baseline plasma DPP-4 activity for males (BSL_male, RFU)")  # Table 5 row BSL_male
    e_sex_bsl  <- 865;           label("Additive shift in BSL for SEXF = 1 (RFU, female - male)")           # Table 5 row BSL_female: 10700 + 865 = 11565
    emax       <- 0.924;         label("Maximum fractional decrease in DPP-4 activity (Emax/100)")          # Table 5 row Emax (92.4 %)
    lec50      <- log(3.06);     label("Linagliptin concentration giving half-maximum DPP-4 inhibition (EC50, nmol/L)")  # Table 5 row EC50
    hill       <- 3.22;          label("Hill coefficient of the sigmoid Emax PD model")                      # Table 5 row HILL

    # ---- PD covariate effects (Retlich 2015 Table 5 footnotes b and c) ----
    e_bsl_ec50    <- 7.92e-5;   label("Fractional change in EC50 per RFU BSL_i deviation from 11,600 RFU (= 0.00792 % per RFU)") # Table 5 row BSL_EC50
    e_ggt_bsl     <- 0.00153;   label("Fractional change in BSL per U/L GGT deviation from 32.3 U/L when GGT <= 175 (= 0.153 % per U/L)") # Table 5 row GGT_BSL
    e_ggt_bsl_hi  <- 0.213;     label("Constant fractional shift in BSL when GGT > 175 U/L (= 21.3 %)")  # Table 5 row GGT_BSL2
    e_alt_bsl     <- 0.00175;   label("Fractional change in BSL per U/L ALT deviation from 28.8 U/L (= 0.175 % per U/L)")  # Table 5 row ALT_BSL
    e_fpg_bsl     <- 0.0146;    label("Fractional change in BSL per mmol/L FPG deviation from 8.90 mmol/L (= 1.46 % per mmol/L)")  # Table 5 row FPG_BSL
    e_trig_bsl    <- 0.000294;  label("Fractional change in BSL per mg/dL TRIG deviation from 160 mg/dL (= 0.0294 % per mg/dL)")  # Table 5 row TRIG_BSL
    e_tchol_bsl   <- 0.000261;  label("Fractional change in BSL per mg/dL TCHOL deviation from 183 mg/dL (= 0.0261 % per mg/dL)") # Table 5 row CHOL_BSL
    e_trig_ec50   <- -0.000153; label("Fractional change in EC50 per mg/dL TRIG deviation from 160 mg/dL (= -0.0153 % per mg/dL)") # Table 5 row TRIG_EC50

    # ---- PD inter-individual variability (Retlich 2015 Table 5) ----
    etalbsl   ~ 0.02808  # omega^2 = log(1 + 0.169^2) = 0.02808 (CV 16.9%)  Table 5 row xBSL
    etalec50  ~ 0.02341  # omega^2 = log(1 + 0.154^2) = 0.02341 (CV 15.4%)  Table 5 row xEC50

    # ---- PD residual error (Retlich 2015 Table 5) ----
    propSd_Dpp4Act <- 0.148;  label("Proportional residual error for plasma DPP-4 activity (fraction)")  # Table 5 row r_prop
  })

  model({
    # ---- Constants ----
    # Linagliptin molecular weight; small molecule (BI 1356, xanthine-based DPP-4 inhibitor),
    # CAS 668270-12-0, molecular formula C25H28N8O2 -> 472.54 g/mol.
    mw_linag       <- 472.54
    nmol_per_mg    <- 1e6 / mw_linag

    # ---- Population reference values for centred covariate effects (Retlich 2015 Tables 4 / 5 footnotes) ----
    ref_wt          <- 88     # kg     (Table 4 footnote b, popPK median)
    ref_dose        <- 5      # mg     (Table 4 footnotes c and g, central reference dose)
    ref_age         <- 60     # years  (Table 4 footnote g, popPK median)
    ref_ggt_pk      <- 33     # U/L    (Table 4 footnote f, popPK median)
    ref_dpp4_pk     <- 12497  # RFU    (Table 4 footnote g, popPK median DPP-4 activity)
    ref_ggt_pd      <- 32.3   # U/L    (Table 5 footnote b, popPK/PD median)
    ref_alt         <- 28.8   # U/L    (Table 5 footnote b)
    ref_fpg         <- 8.90   # mmol/L (Table 5 footnote b)
    ref_trig        <- 160    # mg/dL  (Table 5 footnotes b and c)
    ref_tchol       <- 183    # mg/dL  (Table 5 footnote b)
    ref_bsl_ec50    <- 11600  # RFU    (Table 5 footnote c)

    # ---- Individual PK parameters ----
    # Bioavailability F (multiplicative on dose-input mass; reference F = 1)
    f_indiv <- exp(lfdepot) *
               (1 + e_wt_f * (WT - ref_wt)) *
               (1 + e_metformin_f * CONMED_METFORMIN) *
               exp(etalfdepot)

    # Absorption rate constant (formulation- and dose-dependent)
    ka <- exp(lka + e_pow_ka * FORM_POWDER + e_tab1_ka * FORM_LINAG_TAB1) *
          (1 + e_dose_ka * (DOSE - ref_dose)) *
          exp(etalka)

    # Apparent volumes and clearances
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q  <- exp(lq)
    cl <- exp(lcl + etalcl) * (1 + e_ggt_cl * (GGT - ref_ggt_pk))

    # Central-compartment binding-site concentration (DPP-4 in plasma)
    bmaxc_base <- exp(lbmaxc) + e_sex_bmaxc * SEXF
    bmaxc      <- bmaxc_base *
                  (1 + e_dpp4_bmaxc * (DPP4_BL_RFU - ref_dpp4_pk)) *
                  (1 + e_dose_bmaxc * (DOSE - ref_dose)) *
                  (1 + e_age_bmaxc * (AGE - ref_age)) *
                  exp(etalbmaxc)

    # Peripheral-compartment binding-site concentration (Amax_P expressed in nmol; divide by VP to get nmol/L)
    bmaxp <- exp(lamax_p) / vp

    # DPP-4 dissociation constant (same in both compartments)
    kd <- exp(lkd)

    # ---- Quasi-equilibrium free-vs-total drug algebra ----
    # For each compartment with binding-site concentration B and dissociation constant Kd,
    # the total concentration C_total decomposes into free and bound (complex) via
    #   C_total = C_free + B * C_free / (Kd + C_free)
    # which is equivalent to the quadratic
    #   complex^2 - (C_total + B + Kd) * complex + B * C_total = 0,
    # solved with the numerically stable form
    #   complex = 0.5 * ((Ctot + B + Kd) - sqrt((Ctot + B + Kd)^2 - 4 * B * Ctot)).
    # Then C_free = C_total - complex. (Following Gibiansky 2008 / Papachristos 2020 QSS-TMDD convention.)
    ctot_c    <- central / vc
    discr_c   <- (ctot_c + bmaxc + kd)^2 - 4 * bmaxc * ctot_c
    complex_c <- 0.5 * ((ctot_c + bmaxc + kd) - sqrt(discr_c))
    cfree_c   <- ctot_c - complex_c

    ctot_p    <- peripheral1 / vp
    discr_p   <- (ctot_p + bmaxp + kd)^2 - 4 * bmaxp * ctot_p
    complex_p <- 0.5 * ((ctot_p + bmaxp + kd) - sqrt(discr_p))
    cfree_p   <- ctot_p - complex_p

    # ---- ODEs ----
    # Two-compartment with elimination on unbound drug and free-drug-driven inter-compartmental flux
    # (only unbound drug crosses membranes between central and peripheral compartments).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - cl * cfree_c - q * cfree_c + q * cfree_p
    d/dt(peripheral1) <-               q * cfree_c - q * cfree_p

    # ---- Bioavailability (mg dose to nmol amount conversion times relative F) ----
    f(depot) <- f_indiv * nmol_per_mg

    # ---- Observation: linagliptin total plasma concentration ----
    # Internal kinetics use nmol/L (matches Retlich 2015 Tables 4 and 5);
    # the observed Cc is converted to ng/mL so that the user-facing concentration
    # is dimensionally consistent with the mg dose unit (mg -> mass-based ng/mL).
    # Conversion: c_nmolL * MW_linag (g/mol) / 1000 = ng/mL
    #   1 nmol/L * MW (g/mol) = MW * 1e-9 g/L = MW ng/L = MW / 1000 ng/mL.
    Cc <- ctot_c * mw_linag / 1000
    Cc ~ prop(propSd)

    # ---- PD: sigmoid Emax of plasma DPP-4 activity vs total plasma linagliptin (Retlich 2015 Table 5) ----
    # Piecewise GGT effect on BSL: linear deviation below the threshold, constant +21.3 % above.
    ggt_lt    <- 1 + e_ggt_bsl * (GGT - ref_ggt_pd)
    ggt_hi    <- 1 + e_ggt_bsl_hi
    ggt_below <- (GGT <= 175)
    ggt_eff   <- ggt_below * ggt_lt + (1 - ggt_below) * ggt_hi

    bsl_pop <- exp(lbsl) + e_sex_bsl * SEXF
    bsl     <- bsl_pop *
               ggt_eff *
               (1 + e_alt_bsl   * (ALT   - ref_alt)) *
               (1 + e_fpg_bsl   * (FPG   - ref_fpg)) *
               (1 + e_trig_bsl  * (TRIG  - ref_trig)) *
               (1 + e_tchol_bsl * (TCHOL - ref_tchol)) *
               exp(etalbsl)

    ec50 <- exp(lec50) *
            (1 + e_bsl_ec50  * (bsl - ref_bsl_ec50)) *
            (1 + e_trig_ec50 * (TRIG - ref_trig)) *
            exp(etalec50)

    # Sigmoid Emax: DPP-4 activity = BSL * (1 - Emax * Cc_nmolL^HILL / (EC50^HILL + Cc_nmolL^HILL))
    # The PD formula uses linagliptin in nmol/L (matches Retlich 2015 Table 5 units of EC50).
    Dpp4Act <- bsl * (1 - emax * ctot_c^hill / (ec50^hill + ctot_c^hill))
    Dpp4Act ~ prop(propSd_Dpp4Act)
  })
}
