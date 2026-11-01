Pereira_2024_M6495_cyno <- function() {
  description <- "Preclinical (cynomolgus monkey). Population PK/PD model for the anti-ADAMTS-5 NANOBODY M6495 (Pereira 2024): two-compartment PK with first-order SC absorption, estimated bioavailability, and parallel linear plus Michaelis-Menten (target-mediated) elimination, linked to an indirect-response model in which M6495 inhibits the production of the serum aggrecan neo-epitope ARGS. Parameters are body-weight normalised (per kg) and doses are in mg/kg."
  reference <- "Pereira JNS, Ottevaere I, Serruys B, Guehring H, Ladel C, Lindemann S. Translational pharmacokinetic and pharmacodynamic modelling of the anti-ADAMTS-5 NANOBODY(R) (M6495) using the neo-epitope ARGS as a biomarker. J Pharmacokinet Pharmacodyn. 2024;52(1):8. doi:10.1007/s10928-024-09958-z. PMCID: PMC11662058."
  vignette <- "Pereira_2024_M6495"
  units <- list(time = "h", dosing = "mg/kg", concentration = "mg/L")

  # ARGS is the serum aggrecan neo-epitope released by ADAMTS-5 cleavage; it
  # is a paper-specific biomarker turnover state with no canonical name.
  paper_specific_compartments <- c("args")

  covariateData <- list()

  covariatesDataExcluded <- list(
    STUDY = list(
      description = "Study indicator (single-dose PK/PD study versus single and multiple-dose safety study)",
      units = "(binary)",
      type = "binary",
      notes = "Tested as a covariate on baseline ARGS (Supplementary Table 4: baseline 0.427 nM, study effect 0.648, dOFV -11.286) and improved the ARGS pcVPC (Supplementary Fig. 5), but the authors kept the base model without it because the other estimates were nearly unchanged and the human projection replaces the ARGS baseline with measured human values (Results). The coding of the indicator (which study is the reference) is not stated."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "M6495", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "M6495", units = "mg/kg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "M6495", units = "mg/kg", specimen = "serum", verified = TRUE),
    # The ARGS state holds a CONCENTRATION: Eq. 5-7 write the turnover directly
    # in nM, with kin = ARGS_baseline * kout and ARGS(0) = ARGS_baseline.
    args = list(analyte = "aggrecan ARGS neo-epitope", units = "nM", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "cynomolgus monkey",
    n_subjects = 70L,
    n_studies = 2L,
    weight_range = "Not reported; Table 2 converts the per-kg estimates for a 3.9 kg monkey.",
    sex_female_pct = 47.1,
    disease_state = "Healthy cynomolgus monkeys (preclinical PK/PD and safety studies).",
    dose_range = "Single-dose PK/PD study: 0 (vehicle), 0.01, 0.1, 0.3, 1 and 6 mg/kg SC and 1 mg/kg IV. Safety study: 0 (vehicle), 6, 30 and 150 mg/kg SC once weekly for 4 doses, plus single 6 and 30 mg/kg SC doses.",
    regions = NA_character_,
    notes = "Supplementary Tables 1 and 2: 42 animals in the single-dose PK/PD study (7 groups of 3 male and 3 female) and 28 in the safety studies (4 multiple-dose groups of 3 male and 3 female plus 2 single-dose groups of 2 males), giving 33 of 70 female. 907 PK samples (621 + 286, excluding vehicle) and 1466 ARGS samples, of which 111 (7.6%) were below the ARGS LLOQ and were omitted from the fit. Mean ARGS baseline was 0.5 nM in the PK/PD study and 0.3 nM in the safety studies. PK and PD were fitted simultaneously with NONMEM 7.3, FOCE with interaction."
  )

  ini({
    # Table 1 reports the PK parameters normalised to body weight (L/kg,
    # L/kg/h, mg/kg/h); doses are therefore in mg/kg and Cc is in mg/L.
    lka <- log(0.0394); label("First-order SC absorption rate constant, Ka (1/h)")                       # Table 1: Ka 0.0394 /h (RSE 7.6%)
    lfdepot <- log(1); label("SC bioavailability, F (fraction)")                                          # Table 1: F 1 (RSE 6.3%)
    lcl <- log(0.000251); label("Linear clearance, CL (L/kg/h)")                                         # Table 1: CL 0.000251 L/kg/h (RSE 7.3%)
    lvc <- log(0.046); label("Central volume, V1 (L/kg)")                                                # Table 1: V1 0.046 L/kg (RSE 2.8%)
    lq <- log(0.00163); label("Intercompartmental clearance, Q (L/kg/h)")                                # Table 1: Q 0.00163 L/kg/h (RSE 26.9%)
    lvp <- log(0.0311); label("Peripheral volume, V2 (L/kg)")                                            # Table 1: V2 0.0311 L/kg (RSE 10.6%)
    lvmax <- log(0.000417); label("Maximum rate of non-linear elimination, Vmax (mg/kg/h)")              # Table 1: Vmax 0.000417 mg/kg/h (RSE 11.4%)
    lkm <- log(0.265); label("Concentration at half-maximal non-linear elimination, Km (mg/L)")         # Table 1: Km 0.265 mg/L (RSE 14.3%)

    lkout <- log(0.0401); label("First-order ARGS elimination rate constant, kout (1/h)")               # Table 1: Kout 0.0401 /h (RSE 10.8%)
    limax <- log(0.725); label("Maximum fractional inhibition of ARGS production, Imax (fraction)")   # Table 1: Imax 0.725 (RSE 3.4%)
    lic50 <- log(1.17); label("M6495 concentration at half-maximal inhibition, IC50 (mg/L)")            # Table 1: IC50 1.17 mg/L (RSE 27.6%)
    lrbase <- log(0.36); label("Baseline serum ARGS concentration, BASE (nM)")                         # Table 1: BASE 0.36 nM (RSE 9.8%)
    boxcox_limax <- -3.16; label("Box-Cox shape parameter for the interindividual distribution of Imax (unitless)") # Table 1: Boxcox parameter -3.16 (RSE 57%); Eq. 8-10

    # Table 1 reports IIV as CV%. The authors use CV = sqrt(omega^2): the
    # human projection (Methods, supplementary R code) converts the 46% CV of
    # the human ARGS baseline to a variance of 0.21 = 0.46^2, so the same
    # convention is applied here, omega^2 = (CV/100)^2.
    etalka ~ 0.078961 # Table 1: Ka IIV 28.1% CV -> 0.281^2
    etalfdepot ~ 0.022801 # Table 1: F IIV 15.1% CV -> 0.151^2
    etalcl ~ 0.026244 # Table 1: CL IIV 16.2% CV -> 0.162^2
    etalvmax ~ 0.037249 # Table 1: Vmax IIV 19.3% CV -> 0.193^2
    etalimax ~ 0.043264 # Table 1: Imax IIV 20.8% CV -> 0.208^2 (Box-Cox transformed in model())
    etalrbase ~ 0.293764 # Table 1: BASE IIV 54.2% CV -> 0.542^2

    propSd <- 0.152; label("Proportional residual error for M6495 (fraction)")                         # Table 1: RE1 (PK) 15.2%
    propSd_ARGS <- 0.221; label("Proportional residual error for ARGS (fraction)")                     # Table 1: RE2 (PD) 22.1%
  })

  model({
    ka <- exp(lka + etalka)
    fdepot <- exp(lfdepot + etalfdepot)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    vmax <- exp(lvmax + etalvmax)
    km <- exp(lkm)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Imax: Box-Cox transform of its eta (Eq. 8-10) keeps Imax below 1.
    phi_limax <- (exp(etalimax)^boxcox_limax - 1) / boxcox_limax
    imax <- exp(limax + phi_limax)
    ic50 <- exp(lic50)
    kout <- exp(lkout)
    rbase <- exp(lrbase + etalrbase)
    kin <- rbase * kout # Eq. 6, steady state at time 0

    Cc <- central / vc

    # Eq. 2-5. The paper writes the bioavailability into the absorption flux
    # (Ka * F * A); with linear absorption that is identical to scaling the
    # SC dose by F, so f(depot) is used. IV doses go directly into central.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - k12 * central + k21 * peripheral1 - kel * central - vmax * Cc / (km + Cc)
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(args) <- kin * (1 - imax * Cc / (ic50 + Cc)) - kout * args
    args(0) <- rbase # Eq. 7
    f(depot) <- fdepot

    ARGS <- args

    Cc ~ prop(propSd)
    ARGS ~ prop(propSd_ARGS)
  })
}
