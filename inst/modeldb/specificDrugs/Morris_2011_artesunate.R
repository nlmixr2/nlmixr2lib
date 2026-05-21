Morris_2011_artesunate <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for single-dose oral artesunate",
    "(AS) and its active metabolite dihydroartemisinin (DHA) in 26 pregnant and",
    "25 non-pregnant women with asymptomatic Plasmodium falciparum malaria in",
    "the Democratic Republic of Congo (Morris 2011). Each species has a",
    "one-compartment apparent-volume disposition, with mixed zero-order plus",
    "lagged first-order absorption of AS and complete in-vivo conversion of AS",
    "to DHA (no separate AS elimination). Pregnancy increases DHA apparent",
    "clearance by 42.3% relative to non-pregnant controls (the only retained",
    "covariate); the postpartum sub-cohort could not be characterised by a",
    "structural model and is not represented."
  )
  reference <- paste(
    "Morris CA, Onyamboko MA, Capparelli E, Koch MA, Atibu J, Lokomba V,",
    "Douoguih M, Hemingway-Foday J, Wesche D, Ryder RW, Bose C, Wright L,",
    "Tshefu AK, Meshnick S, Fleckenstein L (2011).",
    "Population pharmacokinetics of artesunate and dihydroartemisinin in",
    "pregnant and non-pregnant women with malaria.",
    "Malar J 10:114. doi:10.1186/1475-2875-10-114."
  )
  vignette <- "Morris_2011_artesunate"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  covariateData <- list(
    PREG = list(
      description        = "Pregnancy status",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = pregnant (second trimester 22-26 weeks gestation or third trimester",
        "32-36 weeks); 0 = non-pregnant control. Time-fixed per subject. Acts as",
        "a proportional effect on DHA apparent clearance via cl_dha = cl_dha_typ",
        "* (1 + e_preg_cl_dha * PREG); the structural CLM/F = 64.0 L/h is the",
        "non-pregnant reference and pregnant women have approximately 42.3%",
        "higher CLM (Morris 2011 Table 2). Postpartum data could not be",
        "characterised by any tested structural model (Results, p.123 of the",
        "trimmed source) and are not represented in this final model."
      ),
      source_name        = "PREG"
    )
  )

  population <- list(
    species           = "human",
    n_subjects        = 51L,
    n_pregnant        = 26L,
    n_nonpregnant     = 25L,
    n_studies         = 1L,
    n_observations    = "300 AS + 498 DHA quantifiable plasma concentrations used in the final fit (~41% of AS and ~2% of DHA samples were below the 1 ng/mL LLQ and excluded prior to model building; 1 AS and 1 DHA outlier additionally excluded)",
    age_range         = "18-38 years (median 23 years pregnant, 24 years non-pregnant; Table 1)",
    weight_range      = "40-84 kg (pregnant median 63 kg, non-pregnant median 52 kg; Table 1)",
    sex_female_pct    = 100,
    race_ethnicity    = "African (Democratic Republic of Congo, Kingasani Maternity Clinic, Kinshasa)",
    disease_state     = paste(
      "Asymptomatic Plasmodium falciparum parasitaemia (parasite density",
      "200-300,000 parasites/microlitre at enrollment, slide- and PCR-positive;",
      "HIV seronegative; haematocrit > 30%; no chronic hypertension or",
      "diabetes). Pregnant women in the second (22-26 weeks gestation) or third",
      "(32-36 weeks) trimester paired with a non-pregnant female control",
      "cohort. The previously enrolled pregnant subjects were also studied",
      "three months postpartum but the postpartum data are excluded from the",
      "final structural model."
    ),
    dose_range        = "Single 200 mg oral artesunate (four 50 mg tablets, Guilin Pharmaceutical Co. Ltd) at the start of an inpatient stay. The 200 mg dose was converted to molar units (~520,200 nmol using artesunate MW = 384.42 g/mol) before modelling; concentrations were similarly modelled on a molar basis in nmol/L.",
    sampling          = "Pre-dose plus 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, and 8 h after AS administration (11 nominal samples per subject).",
    regions           = "Democratic Republic of Congo (Kingasani Maternity Clinic, Kinshasa). ClinicalTrials.gov NCT00538382.",
    notes             = "Demographics and modelling cohort from Morris 2011 Table 1 (median (interquartile range))."
  )

  ini({
    # Structural parameters (Morris 2011 Table 2, "Estimate" column).
    # Concentrations were modelled on a molar basis (Methods: AS and DHA
    # concentrations were converted from ng/mL to nmol/L and the 200 mg AS
    # dose was converted to nmols before modelling). Reported typical
    # values therefore correspond to apparent clearance / volume in L/h
    # and L on the molar concentration scale.
    lka     <- log(4.28)
    label("First-order absorption rate constant K12 for the lagged first-order arm of AS absorption (1/h)")  # Morris 2011 Table 2: K12 = 4.28 1/h (%RSE 23.6)

    ldur    <- log(4.04)
    label("Duration D2 of the zero-order arm of AS absorption (h)")  # Morris 2011 Table 2: D2 = 4.04 h (%RSE 19.5)

    llagt   <- log(0.627)
    label("Lag time ALAG1 for the lagged first-order arm of AS absorption (h)")  # Morris 2011 Table 2: ALAG1 = 0.627 h (%RSE 10.9)

    lfdepot <- log(0.864)
    label("Fraction F1 of the oral AS dose absorbed by the lagged first-order arm (unitless); the remaining 1 - F1 = 0.136 is absorbed by the zero-order arm")  # Morris 2011 Table 2: F1 = 0.864 (%RSE 1.56)

    lcl     <- log(895)
    label("Apparent AS elimination clearance CL/F (L/h); equivalent to the AS-to-DHA conversion clearance under the assumption of complete in-vivo conversion")  # Morris 2011 Table 2: CL/F = 895 L/h (%RSE 5.9)

    lvc     <- log(195)
    label("Apparent AS central volume of distribution V2/F (L)")  # Morris 2011 Table 2: V2/F = 195 L (%RSE 16.4)

    lcl_dha <- log(64.0)
    label("Apparent DHA elimination clearance CLM/F in non-pregnant women (L/h); pregnant women have CLM scaled by (1 + e_preg_cl_dha)")  # Morris 2011 Table 2: CLM/F = 64.0 L/h (%RSE 6.53), non-pregnant reference

    lvc_dha <- log(91.4)
    label("Apparent DHA central volume of distribution V3/F (L)")  # Morris 2011 Table 2: V3/F = 91.4 L (%RSE 6.15)

    # Covariate effect: pregnancy on DHA apparent clearance. Encoded as a
    # proportional effect so the structural CLM/F at the reference category
    # PREG = 0 (non-pregnant) is preserved verbatim from Table 2.
    e_preg_cl_dha <- 0.423
    label("Proportional effect of pregnancy on DHA apparent clearance (unitless); cl_dha = cl_dha_typ * (1 + e_preg_cl_dha * PREG), so pregnant women have approximately 42.3% higher DHA CL/F than non-pregnant controls")  # Morris 2011 Table 2: PREG on CLM/F = 0.423 (%RSE 30.3)

    # Inter-individual variability (log-normal). Morris 2011 Table 2 reports
    # variances on the log scale; the listed "%CV" is sqrt(variance) * 100,
    # an informal approximation that holds only for small variances (true
    # log-normal %CV is sqrt(exp(var) - 1) * 100). The model file uses the
    # published variance values directly. The paper text (Results, model
    # development) reports that "the IIV values associated with CL/F and F1
    # were fixed after conclusion of covariate model building due to poor
    # precision in omega estimates", and Table 2 does not list IIV-CL/F or
    # IIV-F1: those random effects are removed from the final model
    # (interpreted here as fixed at zero, the standard NONMEM action for an
    # omega with poor precision).
    etalka     ~ 1.84    # Morris 2011 Table 2: var(eta_K12) = 1.84 (%RSE 25.3; reported as 136 %CV = 100*sqrt(1.84))
    etaldur    ~ 1.33    # Morris 2011 Table 2: var(eta_D2) = 1.33 (%RSE 22.9; reported as 115 %CV)
    etallagt   ~ 0.573   # Morris 2011 Table 2: var(eta_ALAG) = 0.573 (%RSE 20.8; reported as 75.7 %CV)
    etalvc     ~ 0.604   # Morris 2011 Table 2: var(eta_V2/F) = 0.604 (%RSE 30.1; reported as 77.7 %CV)
    etalcl_dha ~ 0.0802  # Morris 2011 Table 2: var(eta_CLM/F) = 0.0802 (%RSE 24.9; reported as 28.3 %CV)
    etalvc_dha ~ 0.0790  # Morris 2011 Table 2: var(eta_V3/F) = 0.0790 (%RSE 34.7; reported as 28.1 %CV)

    # Residual error. Morris 2011 modelled the natural-log-transformed AS
    # and DHA plasma concentrations with additive residual error on the
    # log scale (Methods: "Residual variability (RV) was modelled with an
    # additive model for log-transformed data"). NONMEM additive-on-log-
    # scale residual error maps to proportional residual error in
    # nlmixr2's linear space (per references/nonmem-translation.md and the
    # sibling Hendriksen_2013_artesunate and Birgersson_2019_artesunate
    # models); the propSd values here are the SD on the log scale
    # (sqrt of the variance reported in Table 2).
    propSd     <- sqrt(0.696)
    label("Proportional residual SD for artesunate plasma concentration (SD on the log scale)")  # Morris 2011 Table 2: var_RV(AS) = 0.696 (%RSE 11.6); SD = sqrt(0.696) ~= 0.834
    propSd_dha <- sqrt(0.174)
    label("Proportional residual SD for DHA plasma concentration (SD on the log scale)")  # Morris 2011 Table 2: var_RV(DHA) = 0.174 (%RSE 9.94); SD = sqrt(0.174) ~= 0.417
  })

  model({
    # Individual PK parameters. No allometric weight scaling: weight on
    # V2/F, CLM/F, and V3/F was tested during covariate model building
    # (Morris 2011 Results) but was not retained at the p < 0.001 backward-
    # elimination threshold. The only covariate retained in the final
    # model is pregnancy status on DHA apparent clearance.
    ka     <- exp(lka     + etalka)
    dur_zo <- exp(ldur    + etaldur)
    lagt   <- exp(llagt   + etallagt)
    fr_fo  <- exp(lfdepot)
    cl     <- exp(lcl)
    vc     <- exp(lvc     + etalvc)
    cl_dha <- exp(lcl_dha + etalcl_dha) * (1 + e_preg_cl_dha * PREG)
    vc_dha <- exp(lvc_dha + etalvc_dha)

    # Micro-constants. Under the source paper's assumption of complete,
    # irreversible in-vivo conversion of AS to DHA, all AS clearance is
    # metabolic conversion: kel = cl/vc transfers AS central -> DHA
    # central. DHA is eliminated linearly at kel_dha = cl_dha/vc_dha.
    # Doses and concentrations are tracked on a molar basis (nmol /
    # nmol/L) so no molecular-weight conversion is needed at the
    # parent-to-metabolite step.
    kel     <- cl     / vc
    kel_dha <- cl_dha / vc_dha

    # ODE system. Mixed zero-order plus lagged first-order absorption of
    # AS uses two dose records per administration: a fraction fr_fo
    # enters the depot compartment (first-order absorption with rate ka
    # and lag lagt) and the remaining (1 - fr_fo) enters the central
    # compartment directly as a zero-order input of duration dur_zo. The
    # source NONMEM compartment mapping (Morris 2011 Methods, Base model
    # development) is preserved: CMT 1 = first-order absorption depot
    # (lagged), CMT 2 = AS central (zero-order input target with
    # duration D2), CMT 3 = DHA central.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central
    d/dt(central_dha) <-  kel * central - kel_dha * central_dha

    # Bioavailability and timing for the two parallel absorption arms.
    # rxode2 applies f(), alag(), and dur() to whichever compartment a
    # given dose record targets; the user supplies two simultaneous dose
    # records per administration (one with cmt = depot for the first-
    # order arm, one with cmt = central and rate = -2 for the zero-order
    # arm) carrying the same full molar dose amount, and f() splits the
    # absorbed fraction between the two arms.
    f(depot)     <- fr_fo
    alag(depot)  <- lagt
    f(central)   <- 1 - fr_fo
    dur(central) <- dur_zo

    # Plasma concentrations (nmol/L; molar dose in nmol; volumes in L).
    Cc     <- central     / vc
    Cc_dha <- central_dha / vc_dha

    Cc     ~ prop(propSd)
    Cc_dha ~ prop(propSd_dha)
  })
}
