Vellinga_2025_remimazolam_tol <- function() {
  description <- paste(
    "Sequential population PK/PD logistic model of the probability of",
    "tolerance to laryngoscopy (TOL, no response to a Macintosh-blade",
    "laryngoscopy) under remimazolam with and without remifentanil in 24",
    "healthy adult volunteers (Vellinga 2025, NCT04670471).",
    "The joint remimazolam, CNS 7054 and remifentanil pharmacokinetic model of",
    "modellib('Vellinga_2025_remimazolam') is carried forward unchanged and",
    "drives the drug effect through arterial plasma concentrations (no effect",
    "compartment).",
    "The probability of tolerance is an Emax function of the remimazolam",
    "concentration with zero baseline probability, and remifentanil lowers",
    "the remimazolam C50 by the factor 1 - CR / (CR + EC50R). No CNS 7054",
    "antagonism and no acute tolerance.",
    sep = " "
  )
  reference <- paste(
    "Vellinga R, Koomen JV, Eleveld DJ, Stohr T, Pesic M, Struys MMRF, Colin PJ.",
    "Influence of Remifentanil on the Pharmacokinetics and Pharmacodynamics of",
    "Remimazolam in Healthy Volunteers. Anesthesiology. 2025;142(4):666-679.",
    "doi:10.1097/ALN.0000000000005348. PMC11892992.",
    "Remifentanil layer: Eleveld DJ, Proost JH, Vereecke H, Absalom AR,",
    "Olofsen E, Vuyk J, Struys MMRF. An Allometric Model of Remifentanil",
    "Pharmacokinetics and Pharmacodynamics. Anesthesiology. 2017;126(6):1005-1018.",
    sep = " "
  )
  vignette <- "Vellinga_2025_remimazolam"
  units <- list(
    time = "min",
    dosing = "mg (remimazolam, into central); ug (remifentanil, into central_remifentanil)",
    concentration = "ng/mL (remimazolam, CNS 7054 and remifentanil)"
  )

  # What each ODE state holds. Supplement 5 model code: the parent states and
  # the transit pool carry remimazolam mass in mg (S1 = V1/1000 makes the
  # prediction ng/mL), the CNS 7054 state carries remimazolam-mass equivalents
  # converted to CNS 7054 mass only at the output (RATIO = 425.3/439.3), and
  # the remifentanil states carry ug (S6 = V1R, so ug/L = ng/mL). All
  # concentrations are arterial plasma (Methods, Storage and Analysis of
  # Blood Samples).
  compartmentData <- list(
    central = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    transit1 = list(
      analyte = "CNS 7054 (remimazolam-mass equivalents, before reaching the metabolite pool)",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    central_cns7054 = list(
      analyte = "CNS 7054 (remimazolam-mass equivalents)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_remifentanil = list(analyte = "remifentanil", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1_remifentanil = list(analyte = "remifentanil", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral2_remifentanil = list(analyte = "remifentanil", units = "ug", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Remimazolam and CNS 7054: fixed allometric scaling to 70 kg with exponent 0.75 on CL, Q2, Q3 and CNS 7054 CL, 1 on V1, V2, V3 and the CNS 7054 volume, and -0.25 on the transit rate constant (Supplement 5 model code ALLCL / ALLV / ALLR). Remifentanil: enters the Eleveld 2017 fat-free mass, the weight maturation of clearance and the linear weight correction of V3.",
      source_name = "WEIGHT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Remifentanil layer only (Eleveld 2017): exponential ageing on V1/Q2/Q3, V2/CL and V3 centred on 35 years, the Al-Sallami fat-free-mass maturation terms, and the puberty / elderly window of the female clearance increase.",
      source_name = "AGE"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Remifentanil layer only: Al-Sallami fat-free mass in the Eleveld 2017 body-size scalar.",
      source_name = "HEIGHT"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = male",
      notes = "Remifentanil layer only. The source data code SEX as 1 = male, 2 = female (model code MAL = 2 - SEX, FEM = SEX - 1), so SEXF = SEX - 1. Selects the female Al-Sallami fat-free-mass equation and the female clearance / Q2 / V2 increase of Eleveld 2017.",
      source_name = "SEX"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 24L,
    n_studies = 1L,
    age_range = "19-70 years",
    age_mean = "43 years",
    weight_range = "51-106 kg",
    weight_mean = "74 kg",
    height_range = "159-192 cm",
    height_mean = "175 cm",
    bmi_range = "20.2-29.4 kg/m^2",
    sex_female_pct = 100 * 11 / 24,
    disease_state = "Healthy volunteers, ASA physical status I, BMI 18-30 kg/m^2, stratified by age (18-34, 35-49, 50-70 years) and sex.",
    dose_range = "Remimazolam effect-site TCI (Zhou model, ke0 0.135 1/min) stepped up and down through targets of 50 to 2000 ng/mL over three crossover periods one week apart; remifentanil plasma TCI (Eleveld model) at a constant target of 0 (period 1), 0.5 (period 2) and 0.1, 1.0, 2.0 or 4.0 ng/mL (period 3).",
    regions = "The Netherlands (single centre, University Medical Center Groningen)",
    n_observations = "374 TOL observations over the three periods, assessed only at MOAAS scores of 1 or less.",
    notes = "Three-period crossover dose-ranging trial (NCT04670471). The remimazolam / CNS 7054 model of the period-1 analysis was refitted to all three periods; remifentanil concentrations were described by the Eleveld 2017 model with only the residual error re-estimated. Below-quantification samples were excluded in the final analysis (Supplement 5). Final estimates from Supplement 5 Table 2 and the 'Model (All + interaction)' column of Supplement 5 Table 1."
  )

  ini({
    # ---------------------------------------------------------------------
    # Pharmacokinetic layer, identical to modellib('Vellinga_2025_remimazolam')
    # without its concentration residual errors. The pharmacodynamic models
    # were fit sequentially on individual post hoc PK predictions (Methods,
    # Pharmacodynamic Model Development; Supplement 6 model code reads IV1,
    # ICL, ... from the data), so these values were not re-estimated in the
    # pharmacodynamic step. They are reproduced here so the model is
    # self-contained.
    # ---------------------------------------------------------------------
    # ---------------------------------------------------------------------
    # Remimazolam, Supplement 5 Table 2 ('Parameter estimates'), all values
    # for a 70 kg individual. Where Supplement 5 Table 1 column 'Model (All +
    # interaction)' -- the same final model -- prints one more digit it is
    # used; the Table 2 value is quoted alongside.
    # ---------------------------------------------------------------------
    lvc <- log(5.31)
    label("Remimazolam central volume of distribution for 70 kg (L)") # Suppl 5 Table 1 'All + interaction' V1 = 5.31; Table 2 V1 = 5.3 L (95% CI 4.9-5.8)
    lvp <- log(12.81)
    label("Remimazolam first peripheral volume for 70 kg (L)") # Suppl 5 Table 1 'All + interaction' = 12.81; Table 2 V2 = 12.8 L (11.3-14.6)
    lvp2 <- log(40.04)
    label("Remimazolam second peripheral volume for 70 kg (L)") # Suppl 5 Table 1 'All + interaction' = 40.04; Table 2 V3 = 40.0 L (36.6-43.7)
    lcl <- log(1.15)
    label("Remimazolam elimination clearance for 70 kg (L/min)") # Suppl 5 Table 2 CL = 1.15 L/min (1.06-1.25)
    lq <- log(0.91)
    label("Remimazolam intercompartmental clearance to peripheral1 for 70 kg (L/min)") # Suppl 5 Table 2 Q2 = 0.91 L/min (0.84-0.99)
    lq2 <- log(0.29)
    label("Remimazolam intercompartmental clearance to peripheral2 for 70 kg (L/min)") # Suppl 5 Table 2 Q3 = 0.29 L/min (0.23-0.37)

    # ---------------------------------------------------------------------
    # CNS 7054 formation and disposition, Supplement 5 Table 2. The fraction
    # of remimazolam clearance forming CNS 7054 was not estimated: Table 2
    # footnote cites the Byfavo EPAR (EMA/160756/2021) for 80%.
    # ---------------------------------------------------------------------
    fm <- fixed(0.8)
    label("Fraction of remimazolam clearance forming CNS 7054 (fraction)") # Suppl 5 Table 2 FM = 0.8 (fixed), from the Byfavo EPAR
    lktr <- log(0.27)
    label("Transit rate constant of CNS 7054 formation for 70 kg (1/min)") # Suppl 5 Table 1 'All + interaction' = 0.27; Table 2 KTR = 0.3 1/min (0.2-0.4)
    lcl_cns7054 <- log(0.05)
    label("Apparent CNS 7054 clearance without remifentanil for 70 kg (L/min)") # Suppl 5 Table 2 CLM = 0.05 L/min (0.04-0.05); Table 1 'All + interaction' also 0.05
    lvc_cns7054 <- log(7.77)
    label("Apparent CNS 7054 volume of distribution for 70 kg (L)") # Suppl 5 Table 1 'All + interaction' = 7.77; Table 2 V5 = 7.8 L (7.6-7.9)

    # Remifentanil inhibition of apparent CNS 7054 clearance, main-text
    # equation 1 and Supplement 5 model code (INH = CONC^GAM / (CONC^GAM +
    # IC50^GAM); DADT(5) = -K50*A(5)*(1 - INH) + ...).
    lic50_remifentanil <- log(8.0)
    label("Remifentanil plasma concentration halving the apparent CNS 7054 clearance (ng/mL)") # Suppl 5 Table 2 IC50 = 8.0 ng/mL (5.5-13.5); main text 8.0 (5.5-13.4)
    lhill_remifentanil <- log(0.61)
    label("Hill coefficient of the remifentanil inhibition of CNS 7054 clearance (unitless)") # Suppl 5 Table 1 'All + interaction' Hill factor = 0.61; Table 2 GAM = 0.6 (0.5-0.7)

    # ---------------------------------------------------------------------
    # Remifentanil, Eleveld 2017. All values are held constant in the
    # Supplement 5 control stream ($THETA 1-14 FIX); the comments quote it.
    # THETA 1-6 are already on the log scale.
    # ---------------------------------------------------------------------
    lvc_remifentanil <- fixed(1.75911)
    label("Remifentanil central volume for the reference individual (log L)") # Suppl 5 $THETA(1) = 1.759110 FIX, V1 = 5.81 L
    lvp_remifentanil <- fixed(2.17717)
    label("Remifentanil first peripheral volume for the reference individual (log L)") # Suppl 5 $THETA(2) = 2.177170 FIX, V2 = 8.82 L
    lvp2_remifentanil <- fixed(1.61445)
    label("Remifentanil second peripheral volume for the reference individual (log L)") # Suppl 5 $THETA(3) = 1.614450 FIX, V3 = 5.03 L
    lcl_remifentanil <- fixed(0.946257)
    label("Remifentanil elimination clearance for the reference individual (log L/min)") # Suppl 5 $THETA(4) = 9.462570e-01 FIX, CL = 2.58 L/min
    lq_remifentanil <- fixed(0.540315)
    label("Remifentanil intercompartmental clearance to peripheral1 for the reference individual (log L/min)") # Suppl 5 $THETA(5) = 5.403150e-01 FIX, Q2 = 1.72 L/min
    lq2_remifentanil <- fixed(-2.08388)
    label("Remifentanil intercompartmental clearance to peripheral2 for the reference individual (log L/min)") # Suppl 5 $THETA(6) = -2.083880e+00 FIX, Q3 = 0.12 L/min
    wt50_cl_remifentanil <- fixed(2.87898)
    label("Body weight at half-maximal remifentanil clearance maturation (kg)") # Suppl 5 $THETA(7) = 2.878980e+00 FIX, 'E50 for CL maturation' (SE50)
    e_age_vc_remifentanil <- fixed(-0.00554481)
    label("Ageing coefficient on remifentanil V1, Q2 and Q3 (1/year)") # Suppl 5 $THETA(8) = -5.544810e-03 FIX, 'Aging v1/q2/q3'
    e_age_vp_remifentanil <- fixed(-0.00326985)
    label("Ageing coefficient on remifentanil V2 and CL (1/year)") # Suppl 5 $THETA(9) = -3.269850e-03 FIX, 'Aging v2/cl'
    e_age_vp2_remifentanil <- fixed(-0.0315135)
    label("Ageing coefficient on remifentanil V3 (1/year)") # Suppl 5 $THETA(10) = -3.151350e-02 FIX, 'Aging v3'
    e_sexf_cl_remifentanil <- fixed(0.470405)
    label("Fractional increase of remifentanil CL, Q2 and V2 in females aged 12-45 years (fraction)") # Suppl 5 $THETA(11) = 4.704050e-01 FIX
    e_wt_vp2_remifentanil <- fixed(-0.0260496)
    label("Linear body-weight correction on remifentanil V3 (1/kg)") # Suppl 5 $THETA(12) = -2.604960e-02 FIX, 'Weight correction v3'

    # ---------------------------------------------------------------------
    # Between-subject variability. Supplement 5 Table 2 reports CV% =
    # sqrt(exp(omega^2) - 1) * 100 (table footnote), so omega^2 =
    # log(1 + CV^2). The variances agree with the 'Between-subject
    # variability' rows of Table 1 ('All + interaction': 0.04, 0.29, 0.09,
    # 0.37). V1, V2, V3, Q2 and the CNS 7054 volume carry no IIV ($OMEGA 0.0
    # FIX).
    # ---------------------------------------------------------------------
    etalcl ~ 0.03847 # Suppl 5 Table 2 IIV CL = 19.8% CV; log(1 + 0.198^2)
    etalq2 ~ 0.28576 # Suppl 5 Table 2 IIV Q3 = 57.5% CV; log(1 + 0.575^2)
    etalcl_cns7054 ~ 0.08780 # Suppl 5 Table 2 IIV CLM = 30.3% CV; log(1 + 0.303^2)
    etalktr ~ 0.37080 # Suppl 5 Table 2 IIV KTR = 67.0% CV; log(1 + 0.670^2)

    # Remifentanil (Eleveld 2017) IIV, held constant in the Supplement 5
    # control stream ($OMEGA 1-6 FIX).
    etalvc_remifentanil ~ fixed(0.103576) # Suppl 5 $OMEGA(1) = 1.035760e-01 FIX, V1 CV 33.0%
    etalvp_remifentanil ~ fixed(0.115366) # Suppl 5 $OMEGA(2) = 1.153660e-01 FIX, V2 CV 35.0%
    etalvp2_remifentanil ~ fixed(0.809706) # Suppl 5 $OMEGA(3) = 8.097060e-01 FIX, V3 CV 111.7%
    etalcl_remifentanil ~ fixed(0.0197392) # Suppl 5 $OMEGA(4) = 1.973920e-02 FIX, CL CV 14.1%
    etalq_remifentanil ~ fixed(0.0546616) # Suppl 5 $OMEGA(5) = 5.466160e-02 FIX, Q2 CV 23.7%
    etalq2_remifentanil ~ fixed(0.285485) # Suppl 5 $OMEGA(6) = 2.854850e-01 FIX, Q3 CV 57.5%


    # ---------------------------------------------------------------------
    # Tolerance to laryngoscopy, main-text Table 1 and the Supplement 6 TOL
    # model code. Logistic model with a baseline probability of tolerance of
    # 0 (the estimated baseline was not different from the no-drug bound and
    # was held there) and a remifentanil reduction of the remimazolam C50
    # (main-text equations 3-5):
    #   C50,cor = EC50 * (1 - CR / (CR + EC50R))
    #   P(tolerant) = (CP / C50,cor) / (1 + CP / C50,cor)
    # Competitive antagonism by CNS 7054 was not supported and is absent
    # (CB = 0 in the model code).
    # ---------------------------------------------------------------------
    lec50 <- log(54128)
    label("Remimazolam plasma concentration with half-maximal probability of tolerance without remifentanil (ng/mL)") # Table 1 EC50T = 54,128 ng/mL (19,930-179,871)
    lec50_remifentanil <- log(0.02)
    label("Remifentanil plasma concentration halving the remimazolam C50 (ng/mL)") # Table 1 EC50TR = 0.02 ng/mL (0.01-0.06)

    # Exponential eta on the remimazolam C50 (Suppl 6: EXP(THETA(1) +
    # ETA(1))), printed as CV% and converted with the Supplement 5 rule
    # omega^2 = log(1 + CV^2).
    etalec50 ~ 1.848444 # Table 1 IIV on the remimazolam C50 = 231.3% CV; log(1 + 2.313^2)
  })

  model({
    # Remimazolam and CNS 7054: fixed allometric scaling to 70 kg
    # (Supplement 5 code ALLCL, ALLV, ALLR).
    allcl <- (WT / 70)^0.75
    allv <- WT / 70
    allktr <- (WT / 70)^(-0.25)

    vc <- exp(lvc) * allv
    vp <- exp(lvp) * allv
    vp2 <- exp(lvp2) * allv
    cl <- exp(lcl + etalcl) * allcl
    q <- exp(lq) * allcl
    q2 <- exp(lq2 + etalq2) * allcl
    cl_cns7054 <- exp(lcl_cns7054 + etalcl_cns7054) * allcl
    vc_cns7054 <- exp(lvc_cns7054) * allv
    ktr <- exp(lktr + etalktr) * allktr
    ic50_remifentanil <- exp(lic50_remifentanil)
    hill_remifentanil <- exp(lhill_remifentanil)

    # Remifentanil, Eleveld 2017 as coded in Supplement 5. Body size is the
    # Al-Sallami fat-free mass relative to a 70 kg, 170 cm male; the
    # maturation constants (0.88, 13.4, -12.7, 1.11, 7.1, -1.1) and the
    # fat-free-mass coefficients are part of the published model.
    ht2 <- (HT / 100)^2
    ffm_mat_male <- 0.88 + (1 - 0.88) / (1 + (AGE / 13.4)^(-12.7))
    ffm_mat_female <- 1.11 + (1 - 1.11) / (1 + (AGE / 7.1)^(-1.1))
    ffm_male <- ffm_mat_male * 42.92 * ht2 * WT / (30.93 * ht2 + WT)
    ffm_female <- ffm_mat_female * 37.99 * ht2 * WT / (35.98 * ht2 + WT)
    ffm_ref <- 42.92 * 1.7^2 * 70 / (30.93 * 1.7^2 + 70)
    bsize_remifentanil <- ((1 - SEXF) * ffm_male + SEXF * ffm_female) / ffm_ref

    kmat_remifentanil <- (WT^2 / (WT^2 + wt50_cl_remifentanil^2)) /
      (70^2 / (70^2 + wt50_cl_remifentanil^2))
    kage_v1 <- exp(e_age_vc_remifentanil * (AGE - 35))
    kage_v2 <- exp(e_age_vp_remifentanil * (AGE - 35))
    kage_v3 <- exp(e_age_vp2_remifentanil * (AGE - 35))
    puberty <- AGE^6 / (AGE^6 + 12^6)
    elderly <- AGE^6 / (AGE^6 + 45^6)
    ksex <- 1 + SEXF * puberty * (1 - elderly) * e_sexf_cl_remifentanil
    kwt_v3 <- exp(e_wt_vp2_remifentanil * (WT - 70))

    vc_remifentanil <- exp(lvc_remifentanil + etalvc_remifentanil) *
      bsize_remifentanil * kage_v1
    vp_remifentanil <- exp(lvp_remifentanil + etalvp_remifentanil) *
      bsize_remifentanil * kage_v2 * ksex
    vp2_remifentanil <- exp(lvp2_remifentanil + etalvp2_remifentanil) *
      bsize_remifentanil * kage_v3 * kwt_v3
    cl_remifentanil <- exp(lcl_remifentanil + etalcl_remifentanil) *
      bsize_remifentanil^0.75 * kage_v2 * ksex * kmat_remifentanil
    q_remifentanil <- exp(lq_remifentanil + etalq_remifentanil) *
      (vp_remifentanil / exp(lvp_remifentanil))^0.75 * kage_v1 * ksex
    q2_remifentanil <- exp(lq2_remifentanil + etalq2_remifentanil) *
      (vp2_remifentanil / exp(lvp2_remifentanil))^0.75 * kage_v1

    # Micro-constants. Remimazolam clearance is split into a non-metabolite
    # route (1 - fm) and CNS 7054 formation (fm).
    kel <- cl * (1 - fm) / vc
    kmet <- cl * fm / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2
    kel_cns7054 <- cl_cns7054 / vc_cns7054
    kel_remifentanil <- cl_remifentanil / vc_remifentanil
    k12_remifentanil <- q_remifentanil / vc_remifentanil
    k21_remifentanil <- q_remifentanil / vp_remifentanil
    k13_remifentanil <- q2_remifentanil / vc_remifentanil
    k31_remifentanil <- q2_remifentanil / vp2_remifentanil

    # Remifentanil plasma concentration (ug/L = ng/mL) and its sigmoid
    # inhibition of the apparent CNS 7054 clearance (main-text equation 1).
    Cc_remifentanil <- central_remifentanil / vc_remifentanil
    inh_cl_cns7054 <- 0
    if (Cc_remifentanil > 0) {
      inh_cl_cns7054 <- Cc_remifentanil^hill_remifentanil /
        (Cc_remifentanil^hill_remifentanil + ic50_remifentanil^hill_remifentanil)
    }

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2 - kmet * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    d/dt(transit1) <- kmet * central - ktr * transit1
    d/dt(central_cns7054) <- ktr * transit1 -
      kel_cns7054 * (1 - inh_cl_cns7054) * central_cns7054
    d/dt(central_remifentanil) <- -kel_remifentanil * central_remifentanil -
      k12_remifentanil * central_remifentanil +
      k21_remifentanil * peripheral1_remifentanil -
      k13_remifentanil * central_remifentanil +
      k31_remifentanil * peripheral2_remifentanil
    d/dt(peripheral1_remifentanil) <- k12_remifentanil * central_remifentanil -
      k21_remifentanil * peripheral1_remifentanil
    d/dt(peripheral2_remifentanil) <- k13_remifentanil * central_remifentanil -
      k31_remifentanil * peripheral2_remifentanil

    # Concentrations in ng/mL. CNS 7054 is carried in remimazolam-mass
    # equivalents and converted by the molecular-weight ratio
    # 425.3 / 439.3 printed in the Supplement 5 $ERROR block.
    Cc <- central / vc * 1000
    Cc_cns7054 <- central_cns7054 / vc_cns7054 * 1000 * 425.3 / 439.3

    # Tolerance to laryngoscopy (Supplement 6 TOL $ERROR block).
    ec50 <- exp(lec50 + etalec50)
    ec50_remifentanil <- exp(lec50_remifentanil)
    ec50_tol <- ec50 * (1 - Cc_remifentanil / (Cc_remifentanil + ec50_remifentanil))
    ca_tol <- 0
    if (Cc > 0) {
      ca_tol <- Cc / ec50_tol
    }
    prob_tol <- ca_tol / (1 + ca_tol)
    p0_tol <- 1 - prob_tol

    # Binary endpoint: simulated category 1 is a response to the stimulus
    # (DV 0) and category 2 is tolerance, no response (DV 1).
    tol ~ c(p0_tol)
  })
}
