Kawata_2024_propofol_2cmt <- function() {
  description <- "New two-compartment model for propofol plasma and brain concentrations in Japanese adults and adolescents undergoing awake craniotomy (Kawata 2024; N = 29 model-building patients). Two-compartment intravenous PK estimated from the paper's own data, plus a brain compartment whose input carries a brain/plasma partition coefficient Kp: dCbrain/dt = k1e * (Kp * Cc - Cbrain). IIV on CL and Kp; proportional error on plasma and additive error on brain concentrations; no covariates retained."
  reference <- paste(
    "Kawata M, Yonezawa A, Mineharu Y, Itohara K, Mizota T, Matsui Y, Kikuchi T,",
    "Yamao Y, Hattori EY, Hamada M, Hira D, Furukawa K, Miyamoto S, Terada T,",
    "Matsubara K, Arakawa Y. Development of extended pharmacokinetic models for propofol",
    "based on measured blood and brain concentrations. Sci Rep. 2024;14:6326.",
    "doi:10.1038/s41598-024-56863-z. PMCID: PMC10943190."
  )
  vignette <- "Kawata_2024_propofol"
  units <- list(time = "min", dosing = "mg", concentration = "ug/mL")

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened; no covariate improved the new model (Results, 'PopPK modeling for elimination phase description'). Per-covariate delta OBJ values are printed only for the extended Marsh model (Supplementary Table S2)."
    ),
    HT = list(
      description = "Height",
      units = "m",
      type = "continuous",
      notes = "Paper reports height in m (Supplementary Table S1). Screened; not retained (Results)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened; not retained (Results)."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Screened as (1 + theta * sex) with sex = 1 for females; not retained (Results)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened; not retained (Results)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened; not retained (Results)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = "Paper reports albumin in g/dL (Supplementary Table S1). Screened; not retained (Results)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened; not retained (Results)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened; not retained (Results)."
    )
  )

  compartmentData <- list(
    central = list(analyte = "propofol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "propofol", units = "mg", specimen = "plasma", verified = TRUE),
    brain = list(analyte = "propofol", units = "ug/g", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 29,
    n_studies = 1,
    age_range = "13-76 years",
    age_median = "49 years",
    weight_range = "44.0-92.8 kg",
    weight_median = "62.5 kg",
    sex_female_pct = 44.8,
    race_ethnicity = c(Asian = 100),
    disease_state = "Brain tumour (glioma, 26/29) or other brain lesion (epilepsy, cavernous hemangioma, tuberculoma) undergoing awake craniotomy under propofol anesthesia",
    dose_range = "Intravenous propofol infusion (target-controlled) until awakening; plasma (radial artery) and brain (lobectomy / uncapping margin) samples immediately before and 5, 15, 25 and 75 min after cessation",
    regions = "Japan (Kyoto University Hospital)",
    notes = "Model-building group = first 29 of 57 consecutively enrolled patients (May 2019 - June 2022); the remaining 28 formed an external validation group. Demographics from Supplementary Table S1 (median (range)); 13 of 29 female. 105 plasma and 80 brain concentrations. Japanese ethnicity inferred from the single-centre Kyoto enrolment."
  )

  ini({
    lcl <- log(1.78); label("Clearance CL1 (L/min)") # Table 2: CL1 = 1.78 (5.22 %RSE)
    lvc <- log(41.1); label("Central volume V1 (L)") # Table 2: V1 = 41.1 (17.7 %RSE)
    lq <- log(0.806); label("Intercompartmental clearance CL2 (L/min)") # Table 2: CL2 = 0.806 (11.6 %RSE)
    lvp <- log(123); label("Peripheral volume V2 (L)") # Table 2: V2 = 123 (26.5 %RSE)
    lke0 <- log(0.0149); label("Plasma-brain distribution rate constant k1e = ke0' (1/min)") # Table 2: k1e = 0.0149 (35.8 %RSE)
    lkp_brain <- log(1.39); label("Brain/plasma partition coefficient Kp (unitless)") # Table 2: Kp = 1.39 (9.78 %RSE)

    # IIV: exponential. Supplementary code 2 $OMEGA carries the final variances;
    # Table 2 prints them as CV% = 100 * sqrt(omega^2)
    etalcl ~ 0.0262 # Supplementary code 2 $OMEGA 0.0262; Table 2: IIV for CL1 16.2 CV%
    etalkp_brain ~ 0.148 # Supplementary code 2 $OMEGA 0.148; Table 2: IIV for Kp 38.5 CV%

    # Residual error (Supplementary code 2 $ERROR: CMT 1 Y=F*(1+EPS(1)), else Y=F+EPS(2));
    # $SIGMA 0.045 and 0.772 are variances, Table 2 prints their square roots
    propSd <- 0.212; label("Proportional residual error on plasma Cc (fraction)") # Table 2: proportional error for plasma 21.2 CV%; sqrt($SIGMA 0.045)
    addSd_Cbrain <- 0.879; label("Additive residual error on brain Cbrain (ug/g)") # Table 2: additive error for brain 0.879 ug/g; sqrt($SIGMA 0.772)
  })

  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    ke0 <- exp(lke0)
    kp_brain <- exp(lkp_brain + etalkp_brain)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc

    # Brain state is a concentration (NONMEM A(3) observed with no scale factor);
    # brain-tissue density taken as 1 g/mL so ug/mL plasma maps to ug/g brain.
    # Methods: dCbrain/dt = k1e * Kp * Ccentral - ke0' * Cbrain, with ke0' = k1e.
    d/dt(brain) <- ke0 * kp_brain * Cc - ke0 * brain
    Cbrain <- brain

    Cc ~ prop(propSd)
    Cbrain ~ add(addSd_Cbrain)
  })
}
