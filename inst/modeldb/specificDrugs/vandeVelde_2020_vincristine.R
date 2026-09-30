vandeVelde_2020_vincristine <- function() {
  description <- "Two-compartment population PK model for intravenous vincristine in children with cancer (van de Velde 2020), with body-surface-area-normalised parameters and an administration-method covariate: intercompartmental clearance and peripheral volume are each exp(1.13) = 3.1-fold higher after a push injection (1-5 min, 15-min infusions pooled) than after a 1-h infusion, while clearance and central volume do not depend on administration method."
  reference <- "van de Velde ME, Panetta JC, Wilhelm AJ, van den Berg MH, van der Sluis IM, van den Bos C, Abbink FCH, van den Heuvel-Eibrink MM, Segers H, Chantrain C, van der Werff Ten Bosch J, Willems L, Evans WE, Kaspers GJL. Population Pharmacokinetics of Vincristine Related to Infusion Duration and Peripheral Neuropathy in Pediatric Oncology Patients. Cancers (Basel). 2020;12(7):1789. doi:10.3390/cancers12071789"
  vignette <- "vandeVelde_2020_vincristine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "All four structural parameters are reported per m^2 of body surface area",
        "(Table 2: Cl in L/hr/m2, V1 in L/m2, IC-Cl in L/hr/m2, V2 in L/m2), so the",
        "individual parameter is the per-m^2 value times BSA (linear scaling, no",
        "estimated exponent). The BSA formula is not stated in the source."
      ),
      source_name = "BSA"
    ),
    TINF = list(
      description = "Duration of the intravenous vincristine administration",
      units = "h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The paper's covariate is the binary administration method 'push' (Table 2",
        "footnote: PK = PKpop * exp(beta * push), push = 0 for a 1 h infusion and 1",
        "for a push injection). A push injection was defined as 1-5 min (Methods 4.2);",
        "the two 15-min infusions given off-protocol were analysed as push injections",
        "(Results 2.1); 1 h infusions lasted 60 min (38 min bag plus 22 min flush) or,",
        "at one hospital, 96 min. The model derives push = 1 when TINF < 0.5 h and",
        "push = 0 otherwise; any threshold between 0.25 h and 1 h classifies every",
        "administration in the study identically. Supply TINF as a data column in",
        "addition to the dose record's own duration (rxode2 does not expose the dose",
        "duration to model())."
      ),
      source_name = "push"
    )
  )

  covariatesDataExcluded <- list(
    CONMED_AZOLE = list(
      description = "Concomitant azole antifungal treatment in the week before or on the day of PK sampling",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant azole)",
      notes = "Tested on the population PK parameters and not significant (Results 2.2); only 7 of 70 PK occasions (6 patients) had concurrent azole treatment (Discussion).",
      source_name = "azole antifungal"
    )
  )

  compartmentData <- list(
    central = list(analyte = "vincristine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vincristine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 35L,
    n_studies = 1L,
    n_occasions = 70L,
    n_observations = 425L,
    age_range = "children and adolescents (range not reported)",
    age_mean = "10.06 years (SD 5.6)",
    sex_female_pct = 54,
    race_ethnicity = c(Caucasian = 86, `Non-Caucasian` = 14),
    disease_state = "Pediatric cancer: acute lymphoblastic leukemia 74%, Hodgkin lymphoma 17%, medulloblastoma 3%, low-grade glioma 3%, Wilms tumor 3%",
    dose_range = "Vincristine 1.5 or 2 mg/m2 IV (maximum 2 mg) as a push injection (1-5 min; n = 20) or a 1 h infusion (n = 15); dose capped at 2 mg in 20 patients (37 occasions)",
    regions = "Netherlands (4 centres), Belgium (3 centres)",
    notes = "PK substudy (35 of 90 patients) of a randomized trial of push injection versus 1 h infusion (Dutch Trial Registry NL4019); 1-5 PK occasions per patient with 1-8 samples over 24 h (Table 1, Results 2.1, Methods 4.2). Fitted in Monolix 5.1.0 (SAEM)."
  )

  ini({
    # Final 'Administration Type' model, Table 2 (reference = 1 h infusion)
    lcl <- log(30.5); label("Clearance per m^2 BSA (L/h/m^2)") # Table 2, Administration Type, Cl = 30.5 L/hr/m2 (RSE 10.0%)
    lvc <- log(20.8); label("Central volume per m^2 BSA (L/m^2)") # Table 2, Administration Type, V1 = 20.8 L/m2 (RSE 12.5%)
    lq <- log(34.2); label("Intercompartmental clearance per m^2 BSA, 1 h infusion (L/h/m^2)") # Table 2, Administration Type, IC-Cl = 34.2 L/hr/m2 (RSE 9.5%)
    lvp <- log(127.7); label("Peripheral volume per m^2 BSA, 1 h infusion (L/m^2)") # Table 2, Administration Type, V2 = 127.7 L/m2 (RSE 19.0%)

    e_push_q <- 1.13; label("Effect of push injection on log intercompartmental clearance (unitless)") # Table 2, beta on IC-CL (push) = 1.13 (RSE 4.5%); exp(1.13) = 3.1 (Results 2.2)
    e_push_vp <- 1.13; label("Effect of push injection on log peripheral volume (unitless)") # Table 2, beta on V2 (push) = 1.13 (RSE 19.1%); exp(1.13) = 3.1 (Results 2.2)

    # Monolix omega (SD of the log-normal random effect), printed under a 'CV%' header; variance = omega^2
    etalcl ~ 0.2704 # Table 2, IIV Cl = 0.52 (RSE 16.9%); 0.52^2
    etalvc ~ 0.3025 # Table 2, IIV V1 = 0.55 (RSE 23.9%); 0.55^2
    etalq ~ 0.2304 # Table 2, IIV IC-Cl = 0.48 (RSE 16.4%); 0.48^2
    etalvp ~ 0.1681 # Table 2, IIV V2 = 0.41 (RSE 17.0%); 0.41^2

    propSd <- 0.45; label("Proportional residual error (fraction)") # Table 2, Residual = 0.45 (RSE 4.3%); proportional model (Methods 4.5)
  })

  model({
    # Administration method: push injection (<= 15 min) vs 1 h infusion (Table 2 footnote)
    push <- 0
    if (TINF < 0.5) push <- 1

    cl <- exp(lcl + etalcl) * BSA
    vc <- exp(lvc + etalvc) * BSA
    q <- exp(lq + etalq + e_push_q * push) * BSA
    vp <- exp(lvp + etalvp + e_push_vp * push) * BSA

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # mg / L * 1000 = ng/mL
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
