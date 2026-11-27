Xu_2020_jnj54175446 <- function() {
  description <- "Population PK/PD model for the P2X7 receptor antagonist JNJ-54175446 in healthy adults after single and once-daily oral doses of a liquid formulation (Xu 2020). One-compartment disposition with first-order elimination and parallel absorption: a first-order arm (ka, lag time ALAG1) carrying 92.4 percent of the absorbed dose and a zero-order arm into the central compartment (duration TK0, lag time ALAG2) carrying the rest. Relative bioavailability falls exponentially with the dose level and is raised by a standard or a high-fat meal; body weight enters CL/F and V/F as power functions with reference 70 kg. Log-normal IIV on ka, CL/F and V/F; proportional residual error. The ex vivo LPS + BzATP stimulated IL-1beta release in whole blood is linked directly to the plasma concentration by a sigmoid inhibitory Emax model with a normally distributed baseline and additive residual error. Each oral dose is given as two event records, one to depot and one to central with rate = -2."
  reference <- paste(
    "Xu Y, Miao X, Ravenstijn P, Hijzen A, Schmidt ME, Nandy P, Zhou H.",
    "Translational Model-Informed Dose Selection for a Human Positron",
    "Emission Tomography Imaging Study of JNJ-54175446, a P2X7 Receptor",
    "Antagonist. Clin Transl Sci. 2020;13(2):309-317.",
    "doi:10.1111/cts.12711.",
    sep = " "
  )
  vignette <- "Xu_2020_jnj54175446"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")
  dosing <- c("depot", "central")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects on CL/F and V/F with reference 70 kg: CL/F = 5.10 * (WT / 70)^0.426 and V/F = 272 * (WT / 70)^1.14 (Xu 2020 Table 1 footnotes c and e; File S1 $PK). Cohort median 77.3 kg, range 46.5-113.4 kg.",
      source_name = "BWT"
    ),
    DOSE_JNJ54175446_MG = list(
      description = "JNJ-54175446 dose level of the regimen (mg per administration)",
      units = "mg",
      type = "continuous",
      reference_category = "150 mg (Fr normalisation dose)",
      notes = "Drives the empirical dose effect on relative bioavailability, Fr = exp(-0.178 * DOSE / 150) * EFOOD (Xu 2020 Eq. 1; File S1 $PK 'F1=EXP(-Beta*(DOSE/150))*(1-FK0)*FNURT'). Set it to the per-administration dose (single dose in SAD, the once-daily dose in MAD); it is constant within a regimen. The model reads the dose level from this column rather than from the dose record because the zero-order arm is an infusion record, on which podo() returns a rate.",
      source_name = "DOSE"
    ),
    FED = list(
      description = "Fed state at dosing (any meal)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = "Together with FED_HIGHFAT encodes the paper's three-level FOOD covariate: fasted = (FED 0, FED_HIGHFAT 0); standard meal = (FED 1, FED_HIGHFAT 0), Fr multiplied by 1.63; high-fat meal = (FED 1, FED_HIGHFAT 1), Fr multiplied by 1.76 (Xu 2020 Table 1 'E FOOD' rows and Results: 'subjects with high-fat and standard meals had a higher bioavailability (76% and 63% higher, respectively)'). The MAD study and the PET study dosed with a standard breakfast.",
      source_name = "FOOD"
    ),
    FED_HIGHFAT = list(
      description = "High-fat, high-calorie meal at dosing",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted or standard meal)",
      notes = "High-fat/high-calorie breakfast of SAD Part 3 (50, 300 and 600 mg; Table S1). When FED_HIGHFAT = 1 the standard-meal factor is switched off and Fr is multiplied by 1.76 alone. Xu 2020 Methods code FOOD = 1 for a high-fat meal and FOOD = 2 for a standard meal (as does File S1 $PK), while the Table 1 row labels carry the two codes the other way round; the meal-to-value pairing is the same in Table 1's descriptions and in the Results, which is what is encoded here.",
      source_name = "FOOD"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "JNJ-54175446", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "JNJ-54175446", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 107L,
    n_studies = 2L,
    n_observations = "2560 plasma concentrations; IL-1beta release samples from the same subjects",
    age_range = "18-69 years",
    age_median = "31 years",
    weight_range = "46.5-113.4 kg",
    weight_median = "77.3 kg",
    sex_female_pct = 2.8,
    race_ethnicity = c(White = 87, Other = 13),
    disease_state = "Healthy adults (young men in SAD Part 1 and the MAD study; elderly men and women in SAD Part 2).",
    dose_range = "Oral liquid formulation: single 0.5-300 mg fasted, single 50-600 mg after a high-fat meal, and 50-450 mg once daily for 10 days with a standard meal.",
    food_status = "Fasted 39 percent (42/107), high-fat meal 16 percent (17/107), standard meal 45 percent (48/107).",
    studies = "Phase 1 SAD (NCT02475148) and MAD (NCT02515955); study sites not stated in Xu 2020.",
    notes = "Pooled SAD and MAD phase 1 data (Xu 2020 Results 'Human population PK/PD model analysis results' and Table S1). PK quantified by LC-MS/MS (LLOQ 0.500 ng/mL). IL-1beta release measured ex vivo after LPS priming and BzATP stimulation, as the difference between the LPS + BzATP tube and the LPS-only tube, so single values may be negative. NONMEM 7.4, FOCE-I; the PD model was fitted sequentially on the individual empirical Bayes PK estimates."
  )

  ini({
    # Structural PK: Xu 2020 Table 1 (estimate (RSE%)); typical values for a
    # 70 kg subject. CL and V are apparent values relative to the Fr scale.
    lka <- log(0.428); label("First-order absorption rate constant (1/h)") # Table 1 'K a' = 0.428 (10.8)
    lcl <- log(5.10); label("Apparent clearance CL/F for 70 kg (L/h)") # Table 1 'CL/F' = 5.10 (4.77)
    lvc <- log(272); label("Apparent volume of distribution V/F for 70 kg (L)") # Table 1 'V/F' = 272 (3.56)
    e_wt_cl <- 0.426; label("Power exponent of (WT / 70) on CL/F (unitless)") # Table 1 'BWT on CL/F' = 0.426 (47.6); footnote c
    e_wt_vc <- 1.14; label("Power exponent of (WT / 70) on V/F (unitless)") # Table 1 'BWT on V/F' = 1.14 (12.5); footnote e

    # Parallel absorption. FK0 = 0.0760 of the absorbed dose enters by the
    # zero-order arm, so the first-order fraction is 1 - 0.0760 = 0.924
    # (Results: 'The majority of the oral dose (92.4%) was absorbed via the
    # first-order absorption route'). No IIV, so the logit encoding is exact.
    logitffo <- log(0.924 / 0.076); label("Logit of the dose fraction absorbed by the first-order arm (logit)") # Table 1 'FK0' = 0.0760 (35.7); ffo = 1 - FK0
    ld2 <- log(1.05); label("Duration of the zero-order input TK0 (h)") # Table 1 'TK0' = 1.05 (29.6); File S1 'D2=TK0'
    ltlag <- log(0.287); label("Lag time of the first-order arm ALAG1 (h)") # Table 1 'ALAG1' = 0.287 (8.58)
    ltlag2 <- log(2.03); label("Lag time of the zero-order arm ALAG2 (h)") # Table 1 'ALAG2' = 2.03 (0.0343)

    # Relative bioavailability Fr (Eq. 1; File S1 $PK). Fr itself is the
    # reference 1 (Table 1 'F r' = 1, fixed); the factors below scale it.
    e_dose_fdepot <- 0.178; label("Exponential rate of decline of Fr with DOSE / 150 mg (unitless)") # Table 1 'beta' = 0.178 (11.2); Eq. 1 Fr = exp(-beta * DOSE/150) * EFOOD
    e_fed_fdepot <- 1.63; label("Multiplicative factor on Fr for a standard meal (unitless)") # Table 1 'E FOOD ... regular food' = 1.63 (4.70)
    e_fed_highfat_fdepot <- 1.76; label("Multiplicative factor on Fr for a high-fat meal (unitless)") # Table 1 'E FOOD ... high-fat food' = 1.76 (6.10)

    # IIV: Table 1 reports CV% = sqrt(omega^2) x 100 (footnote), so omega^2 = (CV/100)^2.
    etalka ~ 1.0816 # Table 1 'IIV of K a' 104% -> 1.04^2
    etalcl ~ 0.09 # Table 1 'IIV of CL/F' 30.0% -> 0.300^2
    etalvc ~ 0.039204 # Table 1 'IIV of V/F' 19.8% -> 0.198^2

    propSd <- 0.175; label("Proportional residual error on JNJ-54175446 plasma concentration (fraction)") # Table 1 'Proportional residual error' 17.5% (3.92)

    # Ex vivo stimulated IL-1beta release: sigmoid inhibitory Emax model (Eq. 2).
    e0 <- 798; label("Baseline ex vivo stimulated IL-1beta release E0 (pg/mL)") # Table 1 'E 0' = 798 (2.50)
    limax <- log(747); label("Maximal inhibition of the IL-1beta release (pg/mL)") # Table 1 'Emax' = 747 (1.76)
    lic50 <- log(67.0); label("Plasma concentration giving half-maximal IL-1beta inhibition (ng/mL)") # Table 1 'EC50,IL1beta' = 0.0670 ug/mL (12.5) = 67.0 ng/mL; Results 'EC50,IL1beta release was 67 ng/mL'
    # Not in Table 1. Digitised by the maintainers from the Figure 5 model
    # curve (red): the curve falls from E0 to E0 - Emax between about 40 and
    # 110 ng/mL, which a simple Emax model (gamma = 1) cannot produce. A
    # logit-linear fit to seven points on the falling limb gives
    # gamma = 4.9 with the midpoint at 67.5 ng/mL, matching the tabulated
    # EC50. See the vignette Assumptions and deviations section.
    lhill <- log(4.9); label("Hill coefficient of the IL-1beta inhibition (unitless)") # Figure 5 red curve, digitised (not tabulated)

    etae0 ~ 42436 # Table 1 'IIV of E 0' = 206 pg/mL, normally distributed (Methods) -> SD 206, variance 206^2
    etalic50 ~ 0.201601 # Table 1 'IIV of EC50,IL1beta' 44.9% -> 0.449^2

    addSd_il1b <- 199; label("Additive residual error on IL-1beta release (pg/mL)") # Table 1 'Additive residual error (pg/mL)' = 199 (0.966)
  })

  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    kel <- cl / vc

    ffo <- expit(logitffo)
    d2 <- exp(ld2)
    tlag <- exp(ltlag)
    tlag2 <- exp(ltlag2)

    # Relative bioavailability (Eq. 1). FED_HIGHFAT = 1 replaces the
    # standard-meal factor by the high-fat factor. As printed and as coded
    # in File S1, Fr = exp(-beta * DOSE / 150), which is exp(-0.178) = 0.837
    # (not 1) at the 150 mg reference dose.
    f_food <- e_fed_fdepot^(FED * (1 - FED_HIGHFAT)) * e_fed_highfat_fdepot^FED_HIGHFAT
    fr <- exp(-e_dose_fdepot * DOSE_JNJ54175446_MG / 150) * f_food

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # NONMEM ADVAN2 with two dose records per administration: CMT 1 (depot,
    # F1, ALAG1) and CMT 2 (central, RATE = -2, D2 = TK0, F2, ALAG2).
    f(depot) <- fr * ffo
    alag(depot) <- tlag
    f(central) <- fr * (1 - ffo)
    alag(central) <- tlag2
    dur(central) <- d2

    # File S1 'S2 = V/1000': mg / L x 1000 = ng/mL.
    Cc <- 1000 * central / vc

    # Ex vivo IL-1beta release (Eq. 2), driven directly by Cc.
    base_il1b <- e0 + etae0
    imax <- exp(limax)
    ic50 <- exp(lic50 + etalic50)
    hill <- exp(lhill)
    il1b <- base_il1b - imax * Cc^hill / (ic50^hill + Cc^hill)

    Cc ~ prop(propSd)
    il1b ~ add(addSd_il1b)
  })
}
