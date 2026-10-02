Assmus_2022_emodepside <- function() {
  description <- paste(
    "Three-compartment population PK model with a four-transit-compartment",
    "absorption chain (ka = ktr = 5/MTT) and linear elimination for oral",
    "emodepside in healthy male volunteers pooled from three phase I studies",
    "(single ascending dose, multiple ascending dose and relative",
    "bioavailability; Assmus 2022). CL/F, Q/F and Q2/F scale allometrically",
    "with body weight (exponent 0.75, reference 75 kg) and every volume",
    "linearly. Formulation (ASD-tablet A or B versus the oral liquid service",
    "formulation solution) and food act on both the mean transit time and the",
    "relative bioavailability, and the daily dose per kg body weight",
    "lengthens the mean transit time linearly. Inter-occasion variability on",
    "the mean transit time applies only to the two occasions of the multiple",
    "ascending dose study. Venous plasma (Cc) and dried-blood-spot (Cb)",
    "concentrations are linked by an estimated scaling factor of 0.618 and",
    "carry separate log-scale residual errors."
  )
  reference <- paste(
    "Assmus F, Hoglund RM, Monnot F, Specht S, Scandale I, Tarning J.",
    "Drug development for the treatment of onchocerciasis: Population",
    "pharmacokinetic and adverse events modeling of emodepside.",
    "PLoS Negl Trop Dis. 2022;16(3):e0010219.",
    "doi:10.1371/journal.pntd.0010219. PMCID: PMC8912909.",
    "Parameter values from Table 3 and the S1 Code NONMEM control stream",
    "(supporting information)."
  )
  vignette <- "Assmus_2022_emodepside"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling on all clearance (exponent 0.75) and volume",
        "(exponent 1) parameters, standardised to 75 kg (Methods Eq 4;",
        "S1 Code TVCL = THETA(1)*((WT/75)**0.75)). Pooled-analysis median",
        "79.1 kg, range 53.2-105 kg (Table 2)."
      ),
      source_name = "WT"
    ),
    FED = list(
      description = "Fed state at dosing; 1 = fed, 0 = fasted",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = paste(
        "S1 Code FOOD = 1 (fasted) / 2 (fed); FED = FOOD - 1. Food raises",
        "the mean transit time by 114% and lowers relative bioavailability",
        "by 24.4% (Table 3). 29 of 142 subjects were dosed fed."
      ),
      source_name = "FOOD"
    ),
    FORM_EMODEPSIDE_ASDA = list(
      description = paste(
        "Emodepside amorphous-solid-dispersion tablet A (hypromellose",
        "acetate succinate polymer); 1 = ASD-tablet A, 0 = otherwise"
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (with FORM_EMODEPSIDE_ASDB = 0: the oral liquid service formulation (LSF) solution, 1 mg/mL)",
      notes = paste(
        "S1 Code FORM2 = 3. Mutually exclusive with FORM_EMODEPSIDE_ASDB.",
        "Raises the mean transit time by 243% and lowers relative",
        "bioavailability by 31.4% versus the LSF solution (Table 3). Used",
        "only in the relative bioavailability study (35 subjects)."
      ),
      source_name = "FORM2"
    ),
    FORM_EMODEPSIDE_ASDB = list(
      description = paste(
        "Emodepside amorphous-solid-dispersion tablet B (copovidone",
        "polymer); 1 = ASD-tablet B, 0 = otherwise"
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (with FORM_EMODEPSIDE_ASDA = 0: the oral liquid service formulation (LSF) solution)",
      notes = paste(
        "S1 Code FORM2 = 4. Mutually exclusive with FORM_EMODEPSIDE_ASDA.",
        "Raises the mean transit time by 124% and lowers relative",
        "bioavailability by 20.0% versus the LSF solution (Table 3). The",
        "formulation carried forward to phase II and used in every dose",
        "finding simulation of the paper (31 subjects in the analysis)."
      ),
      source_name = "FORM2"
    ),
    DOSE_EMODEPSIDE_MGKGD = list(
      description = "Daily emodepside dose per kg body weight",
      units = "mg/kg/day",
      type = "continuous",
      reference_category = "0.08 mg/kg/day (centring value)",
      notes = paste(
        "Linear effect on the mean transit time centred at 0.08 mg/kg/day",
        "(S1 Code COV4 = 1+THETA(16)*(DOSE_KG - 0.08)): +105% per",
        "mg/kg/day (Table 3). Compute as total daily dose (mg) divided by",
        "body weight (kg); twice-daily 10 mg in a 75 kg adult is",
        "20/75 = 0.267 mg/kg/day. Observed range in the pooled data",
        "0.012-0.667 mg/kg/day (Discussion)."
      ),
      source_name = "DOSE_KG"
    ),
    OCC = list(
      description = "Occasion index for inter-occasion variability on the mean transit time",
      units = "(count)",
      type = "categorical",
      reference_category = "0 (no inter-occasion variability)",
      notes = paste(
        "Only the multiple ascending dose study carries IOV",
        "(S1 Code IF(STUDY_ID.EQ.2) IOV = ETA(10)*OCC1+ETA(11)*OCC2).",
        "OCC = 1: days 0-6 of dosing (time < 144.1 h); OCC = 2: day 7 to",
        "the last day of dosing (time >= 144.1 h). Single-dose records",
        "(single ascending dose and relative bioavailability studies) and",
        "any simulation that should omit IOV take OCC = 0, which zeroes",
        "both occasion indicators."
      ),
      source_name = "OCC"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "emodepside", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "emodepside", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "emodepside", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "emodepside", units = "mg", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "emodepside", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "emodepside", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "emodepside", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "emodepside", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 142L,
    n_studies = 3L,
    n_observations = "3,123 concentrations (2,892 venous plasma, 231 dried blood spot)",
    age_range = "18-54 years",
    age_median = "32 years",
    weight_range = "53.2-105 kg",
    weight_median = "79.1 kg",
    sex_female_pct = 0,
    race_ethnicity = c(White = 100),
    disease_state = "healthy male volunteers",
    dose_range = paste(
      "single oral doses of 1-40 mg (LSF solution) or 5-10 mg (ASD-tablet",
      "A or B), and 5 mg once daily, 10 mg once daily or 10 mg twice daily",
      "for 10 days (LSF solution)"
    ),
    regions = "United Kingdom (single phase I site, Hammersmith Medicines Research, London)",
    notes = paste(
      "Pooled from the single ascending dose (NCT02661178, n = 47 after",
      "excluding 11 subjects on the discontinued crystalline tablet X), the",
      "multiple ascending dose (NCT03383614, n = 18) and the relative",
      "bioavailability (NCT03383523, n = 77) studies (Table 1). Baseline",
      "demographics in Table 2. 113 subjects dosed fasted, 29 fed; 76",
      "received the LSF solution, 35 ASD-tablet A and 31 ASD-tablet B.",
      "Concentrations below the 1 ng/mL LLOQ (4.66%) were omitted."
    )
  )

  ini({
    # Structural parameters for a 75 kg adult dosed fasted with the LSF
    # solution at 0.08 mg/kg/day. Final estimates from Table 3; S1 Code
    # $THETA initial values are identical to them.
    # S1 Code labels the deep peripheral compartment 'Peripheral CMP1'
    # (THETA 4 Q1 = 8.45, THETA 5 V3 = 647) and the shallow one 'Peripheral
    # CMP2' (THETA 6 Q2 = 4.6, THETA 7 V4 = 44.4). Table 3 and the Discussion
    # number them the other way round (Q1/F 4.60 with Vp1/F 44.4, Q2/F 8.45
    # with Vp2/F 647). The pairings agree; the Table 3 numbering is used here.
    lmtt <- log(0.488)
    label("Mean absorption transit time, fasted LSF solution at 0.08 mg/kg/day (h)") # Table 3 'MTT(h)' 0.488 (3.5% RSE); S1 Code THETA(3)
    lcl <- log(1.29)
    label("Apparent elimination clearance CL/F at 75 kg (L/h)") # Table 3 'CL/F (L/h)' 1.29 (4.6% RSE); S1 Code THETA(1)
    lvc <- log(52.4)
    label("Apparent central volume Vc/F at 75 kg (L)") # Table 3 'Vc/F (L)' 52.4 (3.4% RSE); S1 Code THETA(2)
    lq <- log(4.60)
    label("Apparent intercompartmental clearance to the shallow peripheral compartment Q1/F at 75 kg (L/h)") # Table 3 'Q1/F (L/h)' 4.60 (5.1% RSE); S1 Code THETA(6) Q2 = 4.6
    lvp <- log(44.4)
    label("Apparent shallow peripheral volume Vp1/F at 75 kg (L)") # Table 3 'Vp1/F (L)' 44.4 (8.7% RSE); S1 Code THETA(7) V4 = 44.4
    lq2 <- log(8.45)
    label("Apparent intercompartmental clearance to the deep peripheral compartment Q2/F at 75 kg (L/h)") # Table 3 'Q2/F (L/h)' 8.45 (2.9% RSE); S1 Code THETA(4) Q1 = 8.45
    lvp2 <- log(647)
    label("Apparent deep peripheral volume Vp2/F at 75 kg (L)") # Table 3 'Vp2/F (L)' 647 (4.5% RSE); S1 Code THETA(5) V3 = 647
    lfdepot <- fixed(log(1))
    label("Relative oral bioavailability of the LSF solution, fasted (unitless)") # Table 3 'F' 1 fixed; S1 Code THETA(8) 1 FIX
    bpr <- 0.618
    label("Dried-blood-spot to venous-plasma concentration ratio (unitless)") # Table 3 'Venous plasma-DBS scaling factor (%)' 61.8 (2.9% RSE); S1 Code THETA(9) 0.618

    # Allometric exponents, fixed (Methods Eq 4).
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on CL/F, Q1/F and Q2/F (unitless)") # Methods: 'The exponent n was fixed to 0.75 for all clearance parameters'; S1 Code (WT/75)**0.75
    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on Vc/F, Vp1/F and Vp2/F (unitless)") # Methods: 'and to 1 for all volume of distribution parameters'; S1 Code (WT/75)**1.00

    # Covariate effects: fractional changes, parameter * (1 + e * covariate).
    e_form_emodepside_asda_fdepot <- -0.314
    label("Fractional change in relative bioavailability for ASD-tablet A versus the LSF solution (unitless)") # Table 3 'Formulation effect on F (%)' ASD-tablet A -31.4 (10.9% RSE); S1 Code THETA(10)
    e_form_emodepside_asdb_fdepot <- -0.200
    label("Fractional change in relative bioavailability for ASD-tablet B versus the LSF solution (unitless)") # Table 3 'Formulation effect on F (%)' ASD-tablet B -20.0 (18.8% RSE); S1 Code THETA(11)
    e_form_emodepside_asda_mtt <- 2.43
    label("Fractional change in mean transit time for ASD-tablet A versus the LSF solution (unitless)") # Table 3 'Formulation effect on MTT(%)' ASD-tablet A 243 (12.2% RSE); S1 Code THETA(12)
    e_form_emodepside_asdb_mtt <- 1.24
    label("Fractional change in mean transit time for ASD-tablet B versus the LSF solution (unitless)") # Table 3 'Formulation effect on MTT(%)' ASD-tablet B 124 (14.9% RSE); S1 Code THETA(13). Results prose says 114%; see vignette
    e_fed_fdepot <- -0.244
    label("Fractional change in relative bioavailability when dosed fed (unitless)") # Table 3 'Food intake on F (%)' -24.4 (17.3% RSE); S1 Code THETA(15)
    e_fed_mtt <- 1.14
    label("Fractional change in mean transit time when dosed fed (unitless)") # Table 3 'Food intake on MTT (%)' 114 (22.6% RSE); S1 Code THETA(14)
    e_dose_emodepside_mgkgd_mtt <- 1.05
    label("Fractional change in mean transit time per mg/kg/day of daily dose above 0.08 mg/kg/day (per mg/kg/day)") # Table 3 'Dose on MTT, %per mg/kg/day' 105 (22.3% RSE); S1 Code THETA(16), centred at 0.08

    # Inter-individual variability: S1 Code $OMEGA variances. Table 3 prints
    # the matching %CV = sqrt(exp(omega^2) - 1). No IIV on Vp1/F or on the
    # DBS scaling factor (both $OMEGA 0 FIX).
    etalfdepot ~ 0.0351 # S1 Code $OMEGA IIV_F1 0.0351; Table 3 F IIV 18.9% CV
    etalmtt ~ 0.135 # S1 Code $OMEGA IIV_MT 0.135; Table 3 MTT IIV 38.0% CV
    etalcl ~ 0.044 # S1 Code $OMEGA IIV_CL 0.044; Table 3 CL/F IIV 21.2% CV
    etalvc ~ 0.0926 # S1 Code $OMEGA IIV_V2 0.0926; Table 3 Vc/F IIV 31.1% CV
    etalq ~ 0.0843 # S1 Code $OMEGA IIV_Q2 0.0843 (on Q = 4.6); Table 3 Q1/F IIV 29.7% CV
    etalq2 ~ 0.017 # S1 Code $OMEGA IIV_Q1 0.017 (on Q = 8.45); Table 3 Q2/F IIV 13.1% CV
    etalvp2 ~ 0.0908 # S1 Code $OMEGA IIV_V3 0.0908 (on V = 647); Table 3 Vp2/F IIV 30.8% CV

    # Inter-occasion variability on MTT, multiple ascending dose study only.
    # rxode2 has no NONMEM-style occasion level, so each occasion has its own
    # eta; the second is fixed to the first ($OMEGA BLOCK(1) SAME).
    etaiov_mtt_1 ~ 0.0692 # S1 Code $OMEGA BLOCK(1) IOV MT 0.0692; Table 3 MTT IOV 26.8% CV
    etaiov_mtt_2 ~ fixed(0.0692) # S1 Code $OMEGA BLOCK(1) SAME

    # Residual error: additive on log-transformed concentrations, separate
    # for each matrix. S1 Code $SIGMA holds variances; SD = sqrt(variance).
    expSd <- 0.1425
    label("Log-scale residual SD, venous plasma (unitless)") # Table 3 'sigma, venous data' 0.0203 (5.7% RSE); S1 Code $SIGMA 0.0203; sqrt(0.0203) = 0.1425
    expSd_Cb <- 0.1913
    label("Log-scale residual SD, dried blood spot (unitless)") # Table 3 'sigma, capillary data' 0.0366 (21.4% RSE); S1 Code $SIGMA 0.0366; sqrt(0.0366) = 0.1913
  })

  model({
    # 1. Occasion indicators; OCC = 0 zeroes both (no IOV).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2

    # 2. Individual parameters (S1 Code $PK).
    cl <- exp(lcl + etalcl) * (WT / 75)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 75)^e_wt_vc
    q <- exp(lq + etalq) * (WT / 75)^e_wt_cl
    vp <- exp(lvp) * (WT / 75)^e_wt_vc
    q2 <- exp(lq2 + etalq2) * (WT / 75)^e_wt_cl
    vp2 <- exp(lvp2 + etalvp2) * (WT / 75)^e_wt_vc

    mtt <- exp(lmtt + etalmtt + iov_mtt) *
      (1 + e_form_emodepside_asda_mtt * FORM_EMODEPSIDE_ASDA) *
      (1 + e_form_emodepside_asdb_mtt * FORM_EMODEPSIDE_ASDB) *
      (1 + e_fed_mtt * FED) *
      (1 + e_dose_emodepside_mgkgd_mtt * (DOSE_EMODEPSIDE_MGKGD - 0.08))

    fdepot <- exp(lfdepot + etalfdepot) *
      (1 + e_form_emodepside_asda_fdepot * FORM_EMODEPSIDE_ASDA) *
      (1 + e_form_emodepside_asdb_fdepot * FORM_EMODEPSIDE_ASDB) *
      (1 + e_fed_fdepot * FED)

    # 3. Micro-constants. Four transit compartments (NN = 4) between the
    #    absorption compartment and central; every step uses
    #    KTR = (NN + 1) / MTT, including absorption into central (ka = ktr).
    ktr <- 5 / mtt
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # 4. ODE system (S1 Code ADVAN5: depot = COMP 1, transit1-4 = COMP 5-8).
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(central) <- ktr * transit4 - kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # 5. Bioavailability on the dosing (absorption) compartment.
    f(depot) <- fdepot

    # 6. Observations. Dose in mg and volumes in L give mg/L; x 1000 for
    #    ng/mL. Cb is the dried-blood-spot concentration (S1 Code CPC).
    Cc <- central / vc * 1000
    Cb <- Cc * bpr

    Cc ~ lnorm(expSd)
    Cb ~ lnorm(expSd_Cb)
  })
}
