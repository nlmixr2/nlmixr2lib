Chen_2025_hemoporfin <- function() {
  description <- paste0(
    "Three-compartment population PK model with linear elimination for ",
    "intravenous hemoporfin (a porphyrin-derivative photosensitiser used in ",
    "photodynamic therapy of port-wine stain) in a pooled cohort of 24 ",
    "Chinese pediatric patients aged 7-14 years with port-wine stain and 16 ",
    "Chinese healthy adult volunteers aged 20-43 years (Chen 2025, ",
    "NCT03125057 plus a 2012 phase I study; single 5 mg/kg 20-minute IV ",
    "infusion in both trials). All six disposition parameters are scaled by ",
    "empirical allometry on fat-free mass normalised to the adult median FFM ",
    "of 46.1 kg: clearance and both inter-compartmental clearances share one ",
    "estimated exponent (0.474) and all three volumes share an exponent ",
    "fixed at the theoretical value of 1. The estimated 0.474 exponent is ",
    "well below the theory-based 0.75 and is the paper's central finding -- ",
    "small pediatric subjects have a substantially higher clearance per ",
    "kilogram than adults, so weight-proportional down-scaling of the ",
    "approved adult 5 mg/kg dose under-exposes children. A study effect ",
    "carried on bioavailability (F1 = 1 + theta_study * STUDY_PEDIATRIC, ",
    "theta_study = -0.175) absorbs the inter-trial difference between the ",
    "pediatric and adult trials. Concentrations are in ng/mL, matching the ",
    "source control stream scaling S1 = V1/1000 with doses in mg. ",
    "Because a large proportion of the pediatric 24-hour samples were below ",
    "the 2 ng/mL LLOQ, the source was fitted with the Beal M3 method; the ",
    "combined proportional-plus-additive residual model reproduced here is ",
    "the above-LLOQ branch of that likelihood. Five companion ",
    "exposure-response models in the Chen_2025_hemoporfin_* family relate ",
    "AUC(0-30min) to the probability of a favourable efficacy outcome."
  )
  reference <- paste(
    "Chen R, Zhang B, Tao J, Yao Q, Zhou T, Ma L, Xu Z.",
    "Population Pharmacokinetics and Exposure-Response Relationship of",
    "Hemoporfin in Pediatric Patients With Port-Wine Stain.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(8):1449-1457.",
    "doi:10.1002/psp4.70050. PMCID: PMC12439283.",
    "Structure from the final NONMEM control stream in Data S2;",
    "values from Table 2.",
    sep = " "
  )
  vignette <- "Chen_2025_hemoporfin"

  # Doses are in mg and volumes in L; the source control stream sets
  # S1 = V1/1000, which converts mg/L to ng/mL, so the predicted
  # concentration -- and therefore the additive residual SD and the 2 ng/mL
  # LLOQ -- are on the ng/mL scale.
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Body-size metric for all six disposition parameters, normalised to ",
        "46.1 kg (Chen 2025 Results 3.3: body_size_adult,typical was fixed ",
        "at the median FFM of the adult healthy volunteers). Computed with ",
        "the PEDIATRIC fat-free-mass equation of Al-Sallami et al. (Clin ",
        "Pharmacokinet 2015;54:1169-1178), which is reference [20] of Chen ",
        "2025 -- not the adult Janmahasatian equation. The Al-Sallami form ",
        "applies the age-correction multiplier a + (1 - a)/(1 + (PNA/b)^-c) ",
        "to the adult fat-free mass, with a = 0.88, b = 13.4, c = 12.7 for ",
        "males and a = 1.11, b = 7.1, c = 1.1 for females; it therefore ",
        "requires WT, HT, SEXF and age. Chen 2025 selected FFM over total ",
        "body weight and over normal fat mass in the stepwise body-size ",
        "screen of Figure S1."
      ),
      source_name        = "FFM"
    ),
    STUDY_PEDIATRIC = list(
      description        = "Pediatric-trial cohort indicator (1 = the pediatric port-wine-stain trial NCT03125057, 0 = the 2012 adult healthy-volunteer phase I trial)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (adult healthy-volunteer phase I trial)",
      notes              = paste0(
        "Chen 2025 Equation 2: F1 = 1 + theta_study * STUDY. Derived in the ",
        "source control stream from the subject identifier ",
        "(IF (ID.LT.1000) STUDY=1; IF (ID.GT.1000) STUDY=0). The estimate ",
        "-0.175 means pediatric-trial exposure is 17.5% lower than the ",
        "allometric prediction; the authors state the underlying cause ",
        "(study operations, bioanalytical method, or subject physiology) is ",
        "unknown, and they deliberately did NOT correct for it in the ",
        "dose-selection simulations. Set to 0 to simulate the allometric ",
        "prediction without the inter-trial offset."
      ),
      source_name        = "STUDY"
    )
  )

  compartmentData <- list(
    central     = list(analyte = "hemoporfin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "hemoporfin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "hemoporfin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 40L,
    n_studies      = 2L,
    age_range      = "7-13 years (pediatric patients); 20-43 years (adult healthy volunteers)",
    age_median     = "9 years (pediatric); 29 years (adult)",
    weight_range   = "21-72 kg (pediatric); 50-75 kg (adult)",
    weight_median  = "32 kg (pediatric); 62 kg (adult)",
    height_range   = "120-164 cm (pediatric); 157-180 cm (adult)",
    height_median  = "139 cm (pediatric); 167 cm (adult)",
    sex_female_pct = 45,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Port-wine stain (congenital capillary malformation) of the head and neck in the pediatric cohort; healthy volunteers in the adult cohort. Pediatric patients with alanine aminotransferase, aspartate transaminase, total bilirubin, serum creatinine or blood urea nitrogen above 1.5 times the upper limit of normal were excluded.",
    dose_range     = "Single 5 mg/kg intravenous infusion over 20 minutes in both trials.",
    administration = "Intravenous infusion (20 minutes)",
    regions        = "China",
    notes          = paste0(
      "Demographics from Chen 2025 Table 1. 24 pediatric patients ",
      "(10 female, 14 male) and 16 adult healthy volunteers (8 female, ",
      "8 male); the 45% female figure is the pooled 18/40. Adult volunteers ",
      "additionally had body mass index restricted to 19.0-24.0 kg/m^2. ",
      "Pediatric sampling was sparse (0, 7 and 24 hours after the end of ",
      "infusion); adult sampling was intensive (14 post-infusion times ",
      "through 24 hours). Assay LLOQ 2 ng/mL, ULOQ 1000 ng/mL by LC-MS/MS; ",
      "a large proportion of the pediatric 24-hour samples were below the ",
      "LLOQ and were handled with the M3 method. Pediatric patients ",
      "additionally received photodynamic therapy 10-30 minutes after the ",
      "start of infusion (530 nm LED at 60 or 75 mW/cm^2); the authors argue ",
      "photobleaching has negligible effect on systemic PK. Findings apply ",
      "only to children aged 7-14 years -- no clearance maturation function ",
      "is included because maturation is essentially complete by age 7."
    )
  )

  ini({
    # ==================================================================
    # Chen 2025 Table 2, "Fix effect" block (final estimates with
    # bootstrap medians and 90% CIs). Structure -- which parameter
    # carries which exponent, the FFM normalisation, and the F1 study
    # effect -- is taken from the final NONMEM control stream in Data S2
    # ($PK block, ADVAN11 TRANS4). The $THETA / $OMEGA numbers printed
    # in that control stream are INITIAL estimates (e.g. CL 11.2, Q2
    # 0.374, theta1 0.442, theta_study -0.151) and are NOT used here;
    # every value below is the published final estimate from Table 2.
    # ==================================================================

    # ----- Typical disposition parameters at the adult reference FFM of 46.1 kg -----
    lcl  <- log(11.2)  ; label("Clearance at FFM = 46.1 kg (L/h)")                                          # Table 2, CL_adult,typical = 11.2 (RSE 4.71%; bootstrap median 11.2, 90% CI 10.4 to 12.0)
    lq   <- log(0.381) ; label("Inter-compartmental clearance Q2 between central and peripheral1 at FFM = 46.1 kg (L/h)")  # Table 2, Q2_adult,typical = 0.381 (RSE 9.46%; bootstrap median 0.385, 90% CI 0.335 to 0.465)
    lq2  <- log(0.127) ; label("Inter-compartmental clearance Q3 between central and peripheral2 at FFM = 46.1 kg (L/h)")  # Table 2, Q3_adult,typical = 0.127 (RSE 8.86%; bootstrap median 0.128, 90% CI 0.110 to 0.148)
    lvc  <- log(3.45)  ; label("Central volume of distribution V1 at FFM = 46.1 kg (L)")                    # Table 2, V1_adult,typical = 3.45 (RSE 6.42%; bootstrap median 3.46, 90% CI 3.06 to 3.80)
    lvp  <- log(0.675) ; label("First peripheral volume of distribution V2 at FFM = 46.1 kg (L)")           # Table 2, V2_adult,typical = 0.675 (RSE 6.40%; bootstrap median 0.681, 90% CI 0.600 to 0.748)
    lvp2 <- log(1.32)  ; label("Second peripheral volume of distribution V3 at FFM = 46.1 kg (L)")          # Table 2, V3_adult,typical = 1.32 (RSE 6.69%; bootstrap median 1.34, 90% CI 1.17 to 1.48)

    # ----- Allometric exponents on fat-free mass (Chen 2025 Equation 1) -----
    # Equation 1: Parameter_pediatric = Parameter_adult,typical *
    #   (body_size_pediatric / body_size_adult,typical)^theta
    # theta1 is shared by CL, Q2 and Q3; theta2 is shared by V1, V2 and V3.
    # This is the "empirical allometry" arm of the Figure S1 screen: theta1
    # was estimated rather than fixed at the theory-based 0.75, and only
    # theta2 was held at its theoretical value.
    e_ffm_cl <- 0.474        ; label("Allometric exponent on fat-free mass shared by CL, Q2 and Q3 (unitless)")  # Table 2, theta1 = 0.474 (RSE 16.3%; bootstrap median 0.460, 90% CI 0.367 to 0.609). Discussion: the CI upper bound 0.609 excludes the theory-based 0.75.
    e_ffm_vc <- fixed(1)     ; label("Allometric exponent on fat-free mass shared by V1, V2 and V3 (unitless)")  # Table 2, theta2 = 1 FIX (control stream Data S2 $THETA "(0, 1) FIX ;AS V"). Held at the theoretical value for a volume.

    # ----- Study effect on bioavailability (Chen 2025 Equation 2) -----
    e_study_f <- -0.175      ; label("Fractional change in bioavailability for the pediatric trial relative to the adult trial (unitless)")  # Table 2, theta_study = -0.175 (RSE 33.1%; bootstrap median -0.174, 90% CI -0.265 to -0.0738). Results 3.3: adding it improved fit by dOFV = -5.210 on 1 df, p < 0.05.

    # ----- Between-subject variability -----
    # IIV is exponential (Methods 2.7). Table 2 reports it as CV%, converted
    # here to the log-scale variance with omega^2 = log(CV^2 + 1).
    # The control stream carries an ETA on all six disposition parameters,
    # but $OMEGA fixes the Q3, V2 and V3 elements to 0 ("0 FIX"), which is
    # why Table 2 reports IIV for CL, Q2 and V1 only. Those three zero-
    # variance etas are omitted here rather than written as fixed(0): a
    # zero diagonal makes OMEGA singular and rxode2 then fails to simulate.
    etalcl ~ 0.0167584  # Table 2, IIV CL  = 13.0 CV% (RSE 26.2%, shrinkage 17.0%; bootstrap median 12.6, 90% CI 8.72 to 15.3); log(0.130^2 + 1) = 0.0167584
    etalq  ~ 0.0392207  # Table 2, IIV Q2  = 20.0 CV% (RSE 35.5%, shrinkage 46.9%; bootstrap median 20.0, 90% CI 16.1 to 27.6); log(0.200^2 + 1) = 0.0392207
    etalvc ~ 0.0460116  # Table 2, IIV V1  = 21.7 CV% (RSE 31.4%, shrinkage 38.3%; bootstrap median 20.8, 90% CI 14.7 to 26.8); log(0.217^2 + 1) = 0.0460116

    # ----- Residual unexplained variability -----
    # Control stream Data S2 $ERROR: SD = SQRT(PROP*PROP + ADD*ADD) with
    # PROP = THETA(9)*IPRED and ADD = THETA(10). This is exactly rxode2's
    # combined additive-plus-proportional form, on the ng/mL scale.
    propSd <- 0.267 ; label("Proportional residual error (fraction)")     # Table 2, sigma_prop = 26.7 CV% (RSE 3.92%, epsilon shrinkage 7.80%; bootstrap median 26.3, 90% CI 24.5 to 27.9). Control stream THETA(9) is the fraction multiplying IPRED.
    addSd  <- 0.240 ; label("Additive residual error (ng/mL)")            # Table 2, sigma_add = 0.240 ng/mL (RSE 9.57%, epsilon shrinkage 7.80%; bootstrap median 0.225, 90% CI 0.191 to 0.255). Control stream THETA(10).
  })

  model({
    # ----- Individual disposition parameters (Equation 1) -----
    # Control stream Data S2: CL = THETA(1)*EXP(ETA(1))*((FFM/46.1)**THETA(7))
    # with the same THETA(7) on Q2 and Q3, and THETA(8) on V1, V2 and V3.
    cl  <- exp(lcl  + etalcl) * (FFM / 46.1)^e_ffm_cl
    q   <- exp(lq   + etalq)  * (FFM / 46.1)^e_ffm_cl
    q2  <- exp(lq2)           * (FFM / 46.1)^e_ffm_cl
    vc  <- exp(lvc  + etalvc) * (FFM / 46.1)^e_ffm_vc
    vp  <- exp(lvp)           * (FFM / 46.1)^e_ffm_vc
    vp2 <- exp(lvp2)          * (FFM / 46.1)^e_ffm_vc

    # ----- Micro-constants (control stream K, K12, K21, K13, K31) -----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # ----- Three-compartment disposition (NONMEM ADVAN11 TRANS4) -----
    # Dosing is an intravenous infusion into central; the source data set
    # supplies the infusion duration per record (RATE = -2, D1 = DUR).
    d/dt(central) <- -kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # ----- Study effect on bioavailability (Equation 2) -----
    # Control stream: F1 = 1 + STUDY*THETA(11), STUDY = 1 for the pediatric
    # trial and 0 for the adult trial.
    f(central) <- 1 + e_study_f * STUDY_PEDIATRIC

    # ----- Observation and residual error -----
    # Control stream S1 = V1/1000 scales the mg amount in a volume in L to
    # ng/mL, i.e. a factor of 1000 relative to mg/L.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
