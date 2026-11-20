Nakatani_2013_E2110_direct_rat <- function() {
  description <- paste(
    "Preclinical (rat). Direct (no hysteresis) hyperbolic Emax model of",
    "central serotonin 1A (5-HT1A) receptor occupancy (RO) measured by",
    "[11C]WAY-100635 PET in the medial prefrontal cortex (MPFC) and the dorsal",
    "raphe nucleus (DRN) versus plasma concentration of E2110, an",
    "investigational 5-HT1A receptor antagonist for overactive bladder, in",
    "male Sprague-Dawley rats (Nakatani 2013, dose-occupancy analysis: 0.3 to",
    "10 mg/kg orally, PET 4 h post-dose). Emax is 100 percent and the EC50 is",
    "region-specific (Table 1). The paper fitted this model to RO against the",
    "measured plasma concentration at the scan; to make it simulatable from a",
    "dose it is driven here by the same paper's two-compartment oral PK model",
    "of male and female rats, whose parameters the paper does not print and",
    "which were read back from the authors' drawn fitted curves in Figure 6",
    "(apparent oral values; see Nakatani_2013_E2110_rat, which carries the",
    "paper's effect-compartment time-course analysis). That PK model was",
    "fitted at 1 mg/kg only and under-predicts the 4-h plasma concentrations of",
    "the PET rats in Figure 5 (about 2-fold at 0.3-1 mg/kg, more at 3-10 mg/kg),",
    "so dose-driven occupancy is under-estimated there; where measured plasma",
    "concentrations exist, apply the Emax relationship to them directly.",
    "Doses are in mg/kg and volumes and clearances per kg body weight.",
    "Typical-value model: no interindividual variability or residual error was reported (residual",
    "terms are set to zero).",
    sep = " "
  )
  reference <- paste(
    "Nakatani Y, Suzuki M, Tokunaga M, Maeda J, Sakai M, Ishihara H,",
    "Yoshinaga T, Takenaka O, Zhang MR, Suhara T, Higuchi M. A small-animal",
    "pharmacokinetic/pharmacodynamic PET study of central serotonin 1A",
    "receptor occupancy by a potential therapeutic agent for overactive",
    "bladder. PLoS One. 2013;8(9):e75040. doi:10.1371/journal.pone.0075040.",
    "PMCID PMC3781034, PMID 24086433.",
    sep = " "
  )
  vignette <- "Nakatani_2013_E2110_rat"
  units <- list(time = "h", dosing = "mg/kg", concentration = "ng/mL")

  covariateData <- list(
    SEXF = list(
      description = "Biological sex, female indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "1 = female, 0 = male. Selects the female two-compartment PK fit of",
        "Figure 6 instead of the male fit; the female set is encoded as",
        "log-scale offsets from the male set. The direct Emax model itself",
        "was fitted in males only.",
        sep = " "
      ),
      source_name = "male / female SD rats"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "E2110", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "E2110", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "E2110", units = "mg/kg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "rat (Sprague-Dawley)",
    disease_state = "Healthy rats, isoflurane-anaesthetised during the PET scans.",
    dose_range = paste(
      "PD data: E2110 0.3, 1, 3 and 10 mg/kg orally, PET 4 h post-dose. PK",
      "data: 1 mg/kg orally, single dose.",
      sep = " "
    ),
    regions = "Japan (National Institute of Radiological Sciences and Eisai Co., Ltd.)",
    notes = paste(
      "PET dose-occupancy study: 4 male rats per dose level, each with a",
      "baseline scan, scanned 4 h after dosing with a tail-vein plasma sample",
      "taken at scan start; scans of the same rat at least 1 week apart (Figure",
      "5). PK study: 4 male and 4 female rats, 1 mg/kg orally, plasma sampled",
      "0.25-8 h; the mean profile of each sex was fitted. Plasma unbound",
      "fraction 6 percent, concentration-independent.",
      sep = " "
    )
  )

  ini({
    # PK: identical to Nakatani_2013_E2110_rat. The paper prints no PK
    # parameter values; they were read back from the drawn Figure 6 fitted
    # lines (male solid, RMSE 0.0016 log10 units; female dotted, RMSE 0.0015).
    lka <- log(2.352); label("Absorption rate constant, males (1/h)") # Figure 6 male line, read back; not printed
    lvc <- log(9.628); label("Apparent central volume V1/F per kg, males (L/kg)") # Figure 6 male line, read back; not printed
    lcl <- log(15.16); label("Apparent clearance CL/F per kg, males (L/h/kg)") # Figure 6 male line, read back; not printed
    lq <- log(5.305); label("Apparent intercompartmental clearance Q/F per kg, males (L/h/kg)") # Figure 6 male line, read back; not printed
    lvp <- log(7.742); label("Apparent peripheral volume V2/F per kg, males (L/kg)") # Figure 6 male line, read back; not printed
    e_sexf_ka <- log(1.344 / 2.352); label("Log-scale effect of female sex on ka (unitless)") # Figure 6 female line, read back: ka 1.344 1/h
    e_sexf_vc <- log(3.339 / 9.628); label("Log-scale effect of female sex on V1/F (unitless)") # Figure 6 female line, read back: V1/F 3.339 L/kg
    e_sexf_cl <- log(1.830 / 15.16); label("Log-scale effect of female sex on CL/F (unitless)") # Figure 6 female line, read back: CL/F 1.830 L/h/kg
    e_sexf_q <- log(2.161 / 5.305); label("Log-scale effect of female sex on Q/F (unitless)") # Figure 6 female line, read back: Q/F 2.161 L/h/kg
    e_sexf_vp <- log(2.537 / 7.742); label("Log-scale effect of female sex on V2/F (unitless)") # Figure 6 female line, read back: V2/F 2.537 L/kg

    # PD: direct hyperbolic Emax on plasma concentration (Equation 4).
    lemax <- fixed(log(100)); label("Maximum receptor occupancy Emax (percent)") # Results: ROmax held at 100 because 3 and 10 mg/kg gave full occupancy
    lec50_mpfc <- log(3.68); label("Plasma concentration for half-maximal RO in MPFC (ng/mL)") # Table 1 dose-response study, MPFC: 7.51 nM; Results text 3.68 ng/mL
    lec50_drn <- log(2.64); label("Plasma concentration for half-maximal RO in DRN (ng/mL)") # Table 1 dose-response study, DRN: 5.40 nM; Results text 2.64 ng/mL

    # Residual error: not reported.
    propSd <- fixed(0); label("Proportional residual SD on plasma concentration (fraction)") # not reported
    addSd_RO_MPFC <- fixed(0); label("Additive residual SD on MPFC receptor occupancy (percent)") # not reported
    addSd_RO_DRN <- fixed(0); label("Additive residual SD on DRN receptor occupancy (percent)") # not reported
  })

  model({
    ka <- exp(lka + e_sexf_ka * SEXF)
    vc <- exp(lvc + e_sexf_vc * SEXF)
    cl <- exp(lcl + e_sexf_cl * SEXF)
    q <- exp(lq + e_sexf_q * SEXF)
    vp <- exp(lvp + e_sexf_vp * SEXF)

    emax <- exp(lemax)
    ec50_mpfc <- exp(lec50_mpfc)
    ec50_drn <- exp(lec50_drn)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg/kg and volume in L/kg give mg/L; x1000 gives ng/mL.
    Cc <- 1000 * central / vc

    # Receptor occupancy (percent), Equations 3-4, direct on plasma Cp
    RO_MPFC <- emax * Cc / (ec50_mpfc + Cc)
    RO_DRN <- emax * Cc / (ec50_drn + Cc)

    Cc ~ prop(propSd)
    RO_MPFC ~ add(addSd_RO_MPFC)
    RO_DRN ~ add(addSd_RO_DRN)
  })
}
