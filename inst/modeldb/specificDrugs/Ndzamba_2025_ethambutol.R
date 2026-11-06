Ndzamba_2025_ethambutol <- function() {
  description <- paste(
    "Two-compartment population PK model for oral ethambutol in plasma in",
    "South African adults with drug-susceptible tuberculosis and HIV (TB-HAART",
    "study; Ndzamba 2025). Absorption runs through four explicit transit",
    "compartments (NN = 3 fixed, ktr = 4/MTT, MTT 0.717 h) into an absorption",
    "compartment emptied at ka 0.877 1/h; elimination is first order.",
    "Clearance (34.8 L/h) and central volume (156 L) are allometrically scaled",
    "on fat-free mass (computed in the model with the Janmahasatian formula",
    "from weight, height and sex) against a 42 kg reference, and",
    "intercompartmental clearance (34.4 L/h) and peripheral volume (564 L) on",
    "total body weight against 57 kg, with fixed 0.75 / 1 exponents. The",
    "Antib-4 fixed-dose combination has 26.0% higher bioavailability and",
    "37.1% longer MTT than e-275 Rifafour. Between-subject variability on",
    "clearance and seven-occasion between-occasion variability on ka, F and",
    "MTT, inflated 2.20-fold for doses taken unobserved at home; combined",
    "proportional (24.8%) plus fixed additive (0.02 mg/L) residual error.",
    "Doses are mg of ethambutol base (0.737 x ethambutol hydrochloride)."
  )
  reference <- paste(
    "Ndzamba B, Denti P, McIlleron H, Smith P, Mthiyane T, Rustomjee R,",
    "Onyebujoh P, Resendiz-Galvan JE, Wasserman S (2025).",
    "Pharmacokinetics of ethambutol and weight banded dosing in South African",
    "adults newly diagnosed with tuberculosis and HIV.",
    "Antimicrob Agents Chemother 69(2):e01200-24. doi:10.1128/aac.01200-24",
    sep = " "
  )
  vignette <- "Ndzamba_2025_ethambutol"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  # Oral doses enter the first transit compartment, not `depot`: `depot` here
  # is the supplementary control stream's ABS compartment, which sits
  # downstream of the four transit compartments and empties into `central` at
  # KA ($MODEL COMP=(TRANSIT1,DEFDOSE) ... COMP=(ABS); K12 = K23 = K34 = K45 =
  # KTR, K56 = KA).
  dosing <- "transit1"

  compartmentData <- list(
    transit1 = list(
      analyte = "ethambutol",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit2 = list(
      analyte = "ethambutol",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit3 = list(
      analyte = "ethambutol",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit4 = list(
      analyte = "ethambutol",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    depot = list(
      analyte = "ethambutol",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "ethambutol",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "ethambutol",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight; scales the peripheral parameters and is an input to the fat-free mass formula.",
      units = "kg",
      type = "continuous",
      source_name = "WEIGHT",
      notes = paste(
        "Q and Vp are scaled by (WT/57)^0.75 and (WT/57) (Table 2 footnote d;",
        "control stream TVWT = 57). Recorded at baseline only (Discussion).",
        "Cohort median (IQR) 55.2 (50.0-62.4) kg, Table 1."
      )
    ),
    HT = list(
      description = "Body height, an input to the fat-free mass formula.",
      units = "cm",
      type = "continuous",
      source_name = "HEIGHT",
      notes = paste(
        "The control stream's HEIGHT column is already in metres (HTM =",
        "HEIGHT); the model converts the canonical cm column with HT / 100.",
        "Cohort median (IQR) 1.59 (1.54-1.68) m, Table 1."
      )
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male; selects the sex-specific constants of the fat-free mass formula.",
      units = "(binary)",
      type = "binary",
      source_name = "SEX_F1",
      reference_category = "0 (male)",
      notes = paste(
        "Control stream 'IF (SEX_F1.EQ.1) THEN ; female', the same coding as",
        "the canonical. Sex has no effect on any PK parameter other than",
        "through FFM: the Results report that scaling CL by FFM explained a",
        "previously significant lower CL in women. 54.1% female, Table 1."
      )
    ),
    FORM_EMB_ANTIB4 = list(
      description = "1 = Antib-4 (Rusan Pharma) four-drug fixed-dose combination, 0 = e-275 Rifafour (Sanofi-aventis), the reference product.",
      units = "(binary)",
      type = "binary",
      source_name = "FORM",
      reference_category = "0 (e-275 Rifafour)",
      notes = paste(
        "Control stream FORM = 1 Rifafour, FORM = 2 Antib-4, so",
        "FORM_EMB_ANTIB4 = as.integer(FORM == 2). Both products contain",
        "275 mg ethambutol hydrochloride with isoniazid, pyrazinamide and",
        "rifampicin. Antib-4 multiplies F by 1.260 and MTT by 1.371. Patients",
        "were not switched between products; Antib-4 was given to 7 of 61",
        "(11.5%) participants when Rifafour was out of stock."
      )
    ),
    OCC = list(
      description = "Integer dosing-occasion index (1-7) for the between-occasion random effects.",
      units = "(count)",
      type = "categorical",
      source_name = "OCC",
      notes = paste(
        "The control stream multiplexes ETA(17)-ETA(23) (ka),",
        "ETA(24)-ETA(30) (F) and ETA(31)-ETA(37) (MTT) with IF (OCC==k)",
        "blocks for k = 1..7. Table 2 footnote b defines the dose taken",
        "at home before a clinic visit, with its pre-dose sample, as an",
        "occasion separate from the observed clinic dose. With PK visits on",
        "days 0, 7, 14 and 28 and pre-dose samples on days 7, 14 (plus the day",
        "13 admission sample) and 28, this gives 1 = day 0 observed dose,",
        "2/3 = day 7 unobserved/observed, 4/5 = day 14 unobserved/observed,",
        "6/7 = day 28 unobserved/observed. Records with OCC outside 1-7 carry",
        "no between-occasion variability."
      )
    ),
    SELFADMIN = list(
      description = "1 = the dose was taken at home without observation by study staff, 0 = directly observed dose in the clinic.",
      units = "(binary)",
      type = "binary",
      source_name = "OBSERVED",
      reference_category = "0 (directly observed dose)",
      notes = paste(
        "The control stream's OBSERVED column has the opposite polarity",
        "(1 = observed), so SELFADMIN = 1 - OBSERVED. Carries the",
        "'IF (OBSERVED==0)' branch that multiplies the between-occasion etas",
        "of ka, F and MTT by ExtraBOV = THETA(13) (Table 2 'Scaling factor",
        "for BOV of unobserved doses' 2.20). Set SELFADMIN = 0 to simulate",
        "observed dosing."
      )
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age in years",
      units = "years",
      type = "continuous",
      notes = "Explored as a covariate (Methods) and not retained. Cohort median (IQR) 32 (27-37) years, Table 1."
    ),
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault equation",
      units = "mL/min",
      type = "continuous",
      notes = paste(
        "Tested and not retained (Results); the Discussion attributes this",
        "to stable renal function across the cohort. Table 1 prints the unit",
        "as mL/h, a typographical error for mL/min (median 97)."
      )
    ),
    CD4_ABS = list(
      description = "CD4+ T-cell count",
      units = "cells/uL",
      type = "continuous",
      notes = "Tested and not retained (Results). Cohort median (IQR) 254 (159-323) cells/uL, Table 1."
    ),
    CONMED_EFV = list(
      description = "Efavirenz-based antiretroviral therapy started on day 13",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "ART (efavirenz, lamivudine, zidovudine) from day 13 was tested and",
        "not retained (Results); 41 of 61 participants (67.2%) started ART",
        "on day 13, Table 1."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 61L,
    n_studies = 1L,
    age_range = "27-37 years (IQR)",
    age_median = "32 years",
    weight_range = "50.0-62.4 kg (IQR)",
    weight_median = "55.2 kg",
    ffm_median = "42 kg (IQR 36.3-45.5)",
    height_median = "1.59 m (IQR 1.54-1.68)",
    sex_female_pct = 54.1,
    race_ethnicity = "Not tabulated; adults enrolled in Durban, South Africa.",
    disease_state = paste(
      "Newly diagnosed, sputum-positive, drug-susceptible pulmonary",
      "tuberculosis in ART-naive adults living with HIV (CD4 median 254",
      "cells/uL); 41 (67.2%) started efavirenz-based ART on day 13."
    ),
    dose_range = paste(
      "Ethambutol hydrochloride 550-1375 mg once daily in the four-drug",
      "fixed-dose combination (275 mg ethambutol hydrochloride per tablet)",
      "by WHO weight band: 550 mg <37.9 kg (n = 1), 825 mg 38-54.9 kg",
      "(n = 31), 1100 mg 55-70 kg (n = 19), 1375 mg >70 kg (n = 10)."
    ),
    regions = "South Africa (Durban)",
    notes = paste(
      "TB-HAART study, 2007-2008. PK visits on days 0, 7, 14 and 28 after an",
      "overnight fast with samples at 1, 1.5, 2, 4, 8 and 12 h post-dose,",
      "pre-dose samples on days 7, 14 and 28, and one more on day 13. 1542",
      "samples, 10 (0.65%) below the 0.1 mg/L LLOQ. 54 participants (88.5%)",
      "took e-275 Rifafour and 7 (11.5%) Antib-4. Baseline characteristics",
      "are Table 1."
    )
  )

  ini({
    # Disposition. Table 2 typical values refer to FFM 42 kg for CL and Vc
    # (footnote c) and to total body weight 57 kg for Q and Vp (footnote d).
    # Values are apparent and relative to the ethambutol BASE dose: the
    # 275 mg ethambutol hydrochloride tablet was taken as 202.68 mg base
    # (salt factor 0.737, Methods). The supplementary control stream lists
    # initial values carried from an earlier run, so Table 2 governs.
    lcl <- log(34.8)
    label("Apparent clearance CL/F at FFM 42 kg (L/h)") # Table 2 Clearance 34.8 L/h (32.9-36.7); THETA 1
    lvc <- log(156)
    label("Apparent central volume Vc/F at FFM 42 kg (L)") # Table 2 Central volume 156 L (139-172); THETA 2
    lq <- log(34.4)
    label("Apparent intercompartmental clearance Q/F at WT 57 kg (L/h)") # Table 2 Intercompartmental clearance 34.4 L/h (31.9-37.1); THETA 9
    lvp <- log(564)
    label("Apparent peripheral volume Vp/F at WT 57 kg (L)") # Table 2 Peripheral volume 564 L (489-663); THETA 8

    e_ffm_cl <- fixed(0.75)
    label("Allometric exponent on fat-free mass for CL/F (unitless)") # Methods: exponents fixed to 0.75 for clearances; control stream ALLMCL_FFM = (FFM/TVFFM)**0.75
    e_ffm_vc <- fixed(1)
    label("Allometric exponent on fat-free mass for Vc/F (unitless)") # Methods: exponents fixed to 1 for volumes; control stream ALLMV_FFM = (FFM/TVFFM)
    e_wt_q <- fixed(0.75)
    label("Allometric exponent on total body weight for Q/F (unitless)") # Methods: exponents fixed to 0.75 for clearances; control stream ALLMCL_WT = (WT/TVWT)**0.75
    e_wt_vp <- fixed(1)
    label("Allometric exponent on total body weight for Vp/F (unitless)") # Methods: exponents fixed to 1 for volumes; control stream ALLMV_WT = (WT/TVWT)

    # Absorption: four explicit transit compartments at ktr = (NN + 1)/MTT
    # feeding a first-order absorption compartment.
    lka <- log(0.877)
    label("Absorption rate constant ka from the absorption compartment (1/h)") # Table 2 Absorption rate constant 0.877 1/h (0.755-1.02); THETA 3
    lmtt <- log(0.717)
    label("Mean transit time MTT with e-275 Rifafour (h)") # Table 2 Mean transit time 0.717 h (0.641-0.788); THETA 7
    lnn <- fixed(log(3))
    label("Number of transit compartments NN; the chain has NN + 1 compartments (unitless)") # Table 2 Transit compartments (NN) 3 (FIXED); THETA 12
    lfdepot <- fixed(log(1))
    label("Oral bioavailability F with e-275 Rifafour (unitless)") # Table 2 Bioavailability (F) 1 (FIXED); THETA 4

    # Formulation effects, control stream TVBIO = THETA(4)*(1 + THETA(14)) and
    # TVMTT = THETA(7)*(1 + THETA(15)) when FORM = 2 (Antib-4).
    e_form_antib4_fdepot <- 0.260
    label("Fractional change in F with Antib-4 vs e-275 Rifafour (unitless)") # Table 2 Antib-4 formulation effect on bioavailability +26.0% (12.0-41.3); THETA 14
    e_form_antib4_mtt <- 0.371
    label("Fractional change in MTT with Antib-4 vs e-275 Rifafour (unitless)") # Table 2 Antib-4 formulation effect on mean transit time +37.1% (10.7-76.5); THETA 15

    e_selfadmin_iov <- 2.20
    label("Multiplier on the between-occasion etas of ka, F and MTT for an unobserved dose (-fold)") # Table 2 Scaling factor for BOV of unobserved doses 2.20 (1.92-2.49); THETA 13

    # Random effects. Table 2 footnote e prints 'CV% = omega2 x 100', but the
    # control stream $OMEGA initials reproduce the tabulated percentages as
    # sqrt(omega^2) (BSV CL sqrt(0.0391723) = 19.8% vs 19.9%), so the
    # percentages are the log-scale standard deviation and the variances
    # below are (percentage / 100)^2. BSV on Vc, ka, F, Vp, Q and MTT, BOV on
    # CL and between-visit variability on CL are '0 FIX' in the control
    # stream, so those etas are absent.
    etalcl ~ 0.039601 # Table 2 Clearance BSV 19.9% (16.7-23.5); 0.199^2

    # BOV on ka and F: Table 2 attributes 27.8% to ka and 52.8% to F, but
    # these follow the $OMEGA comment labels ('17 BOVBIO', '24 BOVKA') and
    # not the $PK code, which assigns ETA(17)-ETA(23) to BOVKA and
    # ETA(24)-ETA(30) to BOVBIO. The code is what NONMEM estimated, so the
    # 52.8% variance is carried on ka and 27.8% on F. The paper's own
    # Figure 2 supports this reading: its simulated steady-state AUC0-24
    # has a log-scale SD of about 0.31 within each weight band, close to
    # sqrt(0.199^2 + 0.278^2) = 0.34 and far below sqrt(0.199^2 + 0.528^2)
    # = 0.56 (see the vignette).
    etaiov_ka_1 ~ 0.278784 # Table 2 BOV 52.8% (45.4-61.3) printed on the F row; control stream ETA(17) = BOVKA; 0.528^2
    etaiov_ka_2 ~ fixed(0.278784) # equal to occasion 1 per control stream $OMEGA BLOCK(1) SAME
    etaiov_ka_3 ~ fixed(0.278784) # equal to occasion 1 per control stream $OMEGA BLOCK(1) SAME
    etaiov_ka_4 ~ fixed(0.278784) # equal to occasion 1 per control stream $OMEGA BLOCK(1) SAME
    etaiov_ka_5 ~ fixed(0.278784) # equal to occasion 1 per control stream $OMEGA BLOCK(1) SAME
    etaiov_ka_6 ~ fixed(0.278784) # equal to occasion 1 per control stream $OMEGA BLOCK(1) SAME
    etaiov_ka_7 ~ fixed(0.278784) # equal to occasion 1 per control stream $OMEGA BLOCK(1) SAME
    etaiov_fdepot_1 ~ 0.077284 # Table 2 BOV 27.8% (25.2-30.3) printed on the ka row; control stream ETA(24) = BOVBIO; 0.278^2
    etaiov_fdepot_2 ~ fixed(0.077284) # equal to occasion 1 per control stream $OMEGA BLOCK(1) SAME
    etaiov_fdepot_3 ~ fixed(0.077284) # equal to occasion 1 per control stream $OMEGA BLOCK(1) SAME
    etaiov_fdepot_4 ~ fixed(0.077284) # equal to occasion 1 per control stream $OMEGA BLOCK(1) SAME
    etaiov_fdepot_5 ~ fixed(0.077284) # equal to occasion 1 per control stream $OMEGA BLOCK(1) SAME
    etaiov_fdepot_6 ~ fixed(0.077284) # equal to occasion 1 per control stream $OMEGA BLOCK(1) SAME
    etaiov_fdepot_7 ~ fixed(0.077284) # equal to occasion 1 per control stream $OMEGA BLOCK(1) SAME
    etaiov_mtt_1 ~ 0.439569 # Table 2 Mean transit time BOV 66.3% (60.3-74.0); control stream ETA(31) = BOVMTT; 0.663^2
    etaiov_mtt_2 ~ fixed(0.439569) # equal to occasion 1 per control stream $OMEGA BLOCK(1) SAME
    etaiov_mtt_3 ~ fixed(0.439569) # equal to occasion 1 per control stream $OMEGA BLOCK(1) SAME
    etaiov_mtt_4 ~ fixed(0.439569) # equal to occasion 1 per control stream $OMEGA BLOCK(1) SAME
    etaiov_mtt_5 ~ fixed(0.439569) # equal to occasion 1 per control stream $OMEGA BLOCK(1) SAME
    etaiov_mtt_6 ~ fixed(0.439569) # equal to occasion 1 per control stream $OMEGA BLOCK(1) SAME
    etaiov_mtt_7 ~ fixed(0.439569) # equal to occasion 1 per control stream $OMEGA BLOCK(1) SAME

    # Residual error: W = SQRT(ADD**2 + PROP**2) with ADD = THETA(6) +
    # 0.2 * LLOQ, THETA(6) FIX 4.4E-7 and LLOQ = 0.1 mg/L.
    propSd <- 0.248
    label("Proportional residual error (fraction)") # Table 2 Proportional error 24.8% (23.8-25.9); THETA 5
    addSd <- fixed(0.02)
    label("Additive residual error (mg/L)") # Table 2 Additive error 0.02 mg/L (FIXED), footnote f: 20% of the 0.1 mg/L LLOQ
  })

  model({
    # 1. Occasion indicators for the seven-occasion BOV.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)
    oc7 <- (OCC == 7)
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2 + oc3 * etaiov_ka_3 +
      oc4 * etaiov_ka_4 + oc5 * etaiov_ka_5 + oc6 * etaiov_ka_6 +
      oc7 * etaiov_ka_7
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 +
      oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4 +
      oc5 * etaiov_fdepot_5 + oc6 * etaiov_fdepot_6 +
      oc7 * etaiov_fdepot_7
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2 + oc3 * etaiov_mtt_3 +
      oc4 * etaiov_mtt_4 + oc5 * etaiov_mtt_5 + oc6 * etaiov_mtt_6 +
      oc7 * etaiov_mtt_7

    # Unobserved doses inflate the between-occasion etas of all absorption
    # parameters, per $PK 'IF (OBSERVED==0) THEN BOVKA = BOVKA*ExtraBOV ...'.
    iov_scale <- 1 + (e_selfadmin_iov - 1) * SELFADMIN

    # 2. Fat-free mass (Janmahasatian 2005), WT in kg and height in metres;
    # the constants select on sex as in the control stream.
    htm <- HT / 100
    ffm_whsmax <- 42.92 + (37.99 - 42.92) * SEXF
    ffm_whs50 <- 30.93 + (35.98 - 30.93) * SEXF
    ffm <- (ffm_whsmax * htm^2 * WT) / (ffm_whs50 * htm^2 + WT)

    # 3. Formulation factors (1 for e-275 Rifafour).
    f_form <- 1 + e_form_antib4_fdepot * FORM_EMB_ANTIB4
    mtt_form <- 1 + e_form_antib4_mtt * FORM_EMB_ANTIB4

    # 4. Individual parameters. FFM scales the central parameters and total
    # body weight the peripheral ones.
    cl <- exp(lcl + etalcl) * (ffm / 42)^e_ffm_cl
    vc <- exp(lvc) * (ffm / 42)^e_ffm_vc
    q <- exp(lq) * (WT / 57)^e_wt_q
    vp <- exp(lvp) * (WT / 57)^e_wt_vp
    ka <- exp(lka + iov_ka * iov_scale)
    mtt <- exp(lmtt + iov_mtt * iov_scale) * mtt_form
    nn <- exp(lnn)
    fdepot <- exp(lfdepot + iov_fdepot * iov_scale) * f_form

    ktr <- (nn + 1) / mtt
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 5. ODE system: control stream ADVAN5 compartments TRANSIT1-4, ABS,
    # CENTRAL and PERI1.
    d/dt(transit1) <- -ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(depot) <- ktr * transit4 - ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(transit1) <- fdepot

    # 6. Observation and residual error.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
