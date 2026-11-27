Emoto_2016_sirolimus <- function() {
  description <- paste(
    "Postmenstrual-age maturation model for apparent oral whole-blood",
    "sirolimus clearance in neonates, infants and young children (under 4",
    "years) with complicated vascular anomalies (Emoto 2016). This is a",
    "NONMEM $PRED regression, not a concentration model: the observations",
    "were maximum a posteriori Bayesian clearance estimates (one per trough",
    "sample, MW/Pharm), allometrically normalised to 70 kg with a weight",
    "exponent of 0.67, and the fitted function is a sigmoid Emax in",
    "postmenstrual age, CL_wt70 = CLmatured * PMA^hill / (tm50^hill +",
    "PMA^hill), with exponential IIV on all three parameters and a",
    "proportional residual error. The model carries no volume, absorption",
    "or dosing; the individual apparent clearance in L/h is recovered as",
    "CL_abs = CL_wt70 * (WT / 70)^0.67."
  )
  reference <- paste(
    "Emoto C, Fukuda T, Mizuno T, Schniedewind B, Christians U, Adams DM,",
    "Vinks AA.",
    "Characterizing the Developmental Trajectory of Sirolimus Clearance in",
    "Neonates and Infants.",
    "CPT Pharmacometrics Syst Pharmacol. 2016;5(8):411-417.",
    "doi:10.1002/psp4.12096.",
    sep = " "
  )
  vignette <- "Emoto_2016_sirolimus"
  units <- list(
    time = "week (unused by the model; the vignette uses postmenstrual age in weeks as the time axis)",
    dosing = "n/a (no dosing; clearance-maturation regression)",
    concentration = "L/h/70kg (allometrically size-normalised apparent oral sirolimus clearance, observation CL_wt70)"
  )

  covariateData <- list(
    PAGE = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying. Supplied in WEEKS (not the register's default",
        "months) because Emoto 2016 writes the maturation function and its",
        "TM50 = 62.9 weeks on that scale. The cohort had no premature",
        "infants, so the authors computed PMA = postnatal age (weeks) + 40",
        "(Methods, 'Population modeling analysis'). Fitted range 43-233",
        "weeks (Supplemental Table 2)."
      ),
      source_name = "PMA"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Used only to convert the size-normalised prediction CL_wt70 back",
        "to an individual clearance CL_abs = CL_wt70 * (WT / 70)^0.67. The",
        "authors normalised each observed clearance to 70 kg with this",
        "exponent BEFORE the NONMEM fit, so WT does not enter the fitted",
        "maturation function itself. Fitted range 4.0-18 kg at the start of",
        "therapy (Supplemental Table 2)."
      ),
      source_name = "Body Weight"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 24L,
    n_studies = 1L,
    age_range = paste(
      "postmenstrual age 43-233 weeks (median 120) for the maturation fit,",
      "i.e. patients younger than 4 years; the full trial cohort was 21 days",
      "to 19 years"
    ),
    weight_range = "4.0-18 kg (median 11) for the maturation fit; 4.0-101 kg in the full cohort",
    sex_female_pct = 58.3,
    race_ethnicity = "Caucasian 17, African American 4, other 2, unknown 1 (maturation-fit subset)",
    disease_state = paste(
      "Complicated vascular anomalies (including kaposiform",
      "hemangioendothelioma) in the concentration-controlled sirolimus",
      "phase II trial NCT00975819, with adequate liver and renal function and",
      "no CYP3A4 inhibitors or inducers at entry."
    ),
    dose_range = paste(
      "Oral sirolimus started at 0.8 mg/m^2 twice daily, then",
      "concentration-controlled to a trough of 10-15 ng/mL over the 4-year",
      "study."
    ),
    regions = "United States (Cincinnati Children's Hospital Medical Center)",
    notes = paste(
      "316 Bayesian clearance estimates from 24 patients younger than 4",
      "years (Supplemental Table 2). The 0.67 weight exponent was taken from",
      "a separate log-log regression of averaged first-3-month clearance on",
      "body weight in 32 patients aged 3-19 years (Supplemental Table 1 and",
      "Supplemental Figure 1A; slope 0.67, 95% CI 0.47-0.88). Seven KHE",
      "patients received tapering steroids for 1-2 months; steroid use had",
      "no significant covariate effect and their data were retained."
    )
  )

  ini({
    # Sigmoid Emax maturation of size-normalised clearance, Emoto 2016
    # Methods equation CL = CLmatured * PMA^Hill / (TM50^Hill + PMA^Hill);
    # final estimates from Table 2 (identical to Results text and the
    # Figure 1 legend). The deposited control stream (Supplement 2) carries
    # only the initial estimates (17, 60, 3).
    lcl <- log(18.1)
    label("Fully matured size-normalised apparent clearance CLmatured (L/h/70kg)") # Table 2 CLmatured = 18.1 L/h/70kg (RSE 4.7%)
    ltm50_cl <- log(62.9)
    label("Postmenstrual age at half of CLmatured, TM50 (weeks)") # Table 2 TM50 = 62.9 weeks (RSE 8.4%)
    lhill <- log(2.94)
    label("Hill coefficient of the clearance maturation function (unitless)") # Table 2 Hill = 2.94 (RSE 16.4%)

    # Allometric exponent applied by the authors before the NONMEM fit
    # (Results equation for CLsize-normalized to 70kg); not estimated in
    # the population model, so fixed here.
    e_wt_cl <- fixed(0.67)
    label("Allometric weight exponent on clearance (unitless)") # Results / Supplemental Figure 1A: log-log slope 0.67 (95% CI 0.47-0.88), used as a fixed normalisation exponent

    # IIV: Supplement 1 eq. 1, P_i = P_pop * exp(eta_i). Table 2 reports IIV
    # as CV%; converted with omega^2 = log(1 + CV^2). The vignette shows the
    # alternative omega^2 = CV^2 reading gives a near-identical VPC.
    etalcl ~ 0.012688 # Table 2 IIV CLmatured 11.3 CV%: log(1 + 0.113^2)
    etaltm50_cl ~ 0.029826 # Table 2 IIV TM50 17.4 CV%: log(1 + 0.174^2)
    etalhill ~ 0.713147 # Table 2 IIV Hill 102 CV%: log(1 + 1.02^2)

    # Proportional residual error, Supplement 1 eq. 2 Y = F * (1 + eps).
    # Table 2 prints 0.066; read as the NONMEM $SIGMA VARIANCE (SD 0.257),
    # which reproduces the Figure 2 VPC 5th/95th percentiles; the SD reading
    # (6.6%) gives a prediction interval roughly half as wide.
    propSd <- sqrt(0.066)
    label("Proportional residual SD on size-normalised clearance (fraction)") # Table 2 'Proportional residual error' = 0.066 (variance); SD = sqrt(0.066) = 0.2569
  })

  model({
    cl_mat <- exp(lcl + etalcl)
    tm50_cl <- exp(ltm50_cl + etaltm50_cl)
    hill <- exp(lhill + etalhill)

    # Emoto 2016 Methods sigmoid Emax in postmenstrual age (weeks), written
    # in the overflow-safe form 1 / (1 + (TM50/PMA)^Hill), which is
    # algebraically identical to PMA^Hill / (TM50^Hill + PMA^Hill).
    fmat <- 1 / (1 + (tm50_cl / PAGE)^hill)
    CL_wt70 <- cl_mat * fmat

    # Individual apparent oral clearance in L/h (inverse of the Results
    # normalisation CLsize-normalized = CLobserved * (70 / WT)^0.67).
    CL_abs <- CL_wt70 * (WT / 70)^e_wt_cl

    CL_wt70 ~ prop(propSd)
  })
}
