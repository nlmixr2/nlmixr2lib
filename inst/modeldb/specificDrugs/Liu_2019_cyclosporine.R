Liu_2019_cyclosporine <- function() {
  description <- paste(
    "Two-compartment population PK model for whole-blood cyclosporine A after a",
    "single 200 mg oral dose of the Sandimmun Neoral microemulsion in healthy",
    "Chinese adults, given alone or with 500 mg oral baicalin (Liu 2019).",
    "Savic transit-compartment absorption (non-integer number of transit",
    "compartments N and mean transit time MTT) into an absorption compartment",
    "emptied at a first-order rate ka, linear elimination, and a proportional",
    "residual error. Between-subject variability on CL/F, ka, Q/F and N, and",
    "between-occasion (study period) variability on MTT and N. No covariates",
    "were retained, and baicalin co-administration had no clinically relevant",
    "effect on any parameter, so the final model carries no baicalin term."
  )
  reference <- paste(
    "Liu R, Li X, Wei J, Liu S, Chang Y, Zhang J, Zhang J, Zhang X, Fuhr U,",
    "Taubert M, Tian X (2019). A Single Dose of Baicalin Has No Clinically Significant",
    "Effect on the Pharmacokinetics of Cyclosporine A in Healthy Chinese",
    "Volunteers. Frontiers in Pharmacology 10:518.",
    "doi:10.3389/fphar.2019.00518."
  )
  vignette <- "Liu_2019_cyclosporine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "cyclosporine A", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "cyclosporine A", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "cyclosporine A", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    OCC = list(
      description = "Study period (occasion) of the two-period fixed-sequence crossover",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Values 1 and 2. Period 1 (day 1) is cyclosporine alone and period 2 (day 15, after a 2-week washout) is cyclosporine with 500 mg baicalin (Methods 'Study Design'). The paper estimates inter-occasional variability on MTT and N (Results 'Model Development'; Table 3 'IOV MTT' and 'IOV N'), one occasion per period. Decomposed inside model() into the binary indicators oc1 and oc2 that multiplex the per-occasion etas; the occasion-2 variances are fixed equal to the occasion-1 variances (NONMEM $OMEGA BLOCK(1) SAME). A single-dose simulation should still supply OCC = 1.",
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    CONMED_BAICALIN = list(
      description = "Baicalin co-administration indicator (1 = cyclosporine given with a single 500 mg oral baicalin dose, 0 = cyclosporine alone)",
      units = "(binary)",
      type = "binary",
      notes = "The paper's 'test' variable (Eq. 1: PAR = TVPAR * theta_test^test). The baicalin factor was tested on N, MTT, ka, CL/F, Vc/F, Vp/F and Q/F one parameter at a time, each with a 1000-sample bootstrap (Methods 'Effect of Baicalin on the Pharmacokinetics of Cyclosporine A'). Every 95% CI included 1, so no effect was retained in the final model (Table 3). Bootstrap medians (95% CI) digitised by the maintainers from Figure 4: CL 0.99 (0.95-1.04), ka 0.98 (0.29-1.49), N 1.02 (0.89-1.17), MTT 1.02 (0.98-1.05), Vc 0.97 (0.90-1.03), Vp 1.07 (0.98-1.15), Q 1.07 (0.95-1.26)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested in the stepwise covariate model and not retained (Results 'Model Development')."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Tested in the stepwise covariate model and not retained (Results 'Model Development')."
    ),
    SEXF = list(
      description = "Biological sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Tested in the stepwise covariate model and not retained (Results 'Model Development')."
    ),
    HCT = list(
      description = "Hematocrit",
      units = "fraction",
      type = "continuous",
      notes = "Tested in the stepwise covariate model and not retained (Results 'Model Development'); the paper does not report the units used."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Tested in the stepwise covariate model and not retained (Results 'Model Development')."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Tested in the stepwise covariate model and not retained (Results 'Model Development')."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Tested in the stepwise covariate model and not retained (Results 'Model Development'); the paper does not report the units used."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Tested in the stepwise covariate model and not retained (Results 'Model Development'); the paper does not report the units used."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 16L,
    n_studies = 1L,
    age_range = "19-34 years",
    bmi_range = "19.4-25.6 kg/m^2",
    sex_female_pct = 50,
    race_ethnicity = "Chinese",
    disease_state = "Healthy volunteers",
    dose_range = "Single oral 200 mg cyclosporine A (eight 25 mg Sandimmun Neoral soft capsules) after an overnight fast, alone (period 1) and with a single oral 500 mg baicalin dose (period 2)",
    regions = "China (First Affiliated Hospital of Zhengzhou University, Henan)",
    notes = "Single-centre, open-label, two-period, fixed-sequence study with a 2-week washout (Methods 'Study Design'). 8 men and 8 women (Methods 'Study Population'); body weight is not reported. Whole-blood cyclosporine A was sampled pre-dose and at 0.5, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 10, 12, 16, 24, 36 and 48 h after each dose and assayed by UHPLC-MS/MS (LLOQ 10 ng/mL). All 16 subjects completed both periods."
  )

  ini({
    # Table 3 reports the bootstrap median of each parameter (1000 bootstrap
    # samples of the final NONMEM 7.4.1 model); the paper prints no separate
    # point estimate, so the bootstrap medians are encoded. All clearances and
    # volumes are apparent (divided by the unknown oral bioavailability F).
    lcl <- log(62.8)
    label("Apparent clearance CL/F (L/h)") # Table 3 'CL / F (l / h)' = 62.8 (54.4-71.2)
    lvc <- log(254)
    label("Apparent central volume of distribution Vc/F (L)") # Table 3 'Vc / F (l)' = 254 (226-281)
    lvp <- log(388)
    label("Apparent peripheral volume of distribution Vp/F (L)") # Table 3 'Vp / F (l)' = 388 (344-456)
    lq <- log(23.6)
    label("Apparent intercompartmental clearance Q/F (L/h)") # Table 3 'Q (l / h)' = 23.6 (19.3-29.5)
    lka <- log(12.4)
    label("First-order absorption rate constant from the absorption compartment (1/h)") # Table 3 'Ka(h-1)' = 12.4 (6.57-33.9)
    lmtt <- log(0.812)
    label("Absorption mean transit time MTT (h)") # Table 3 'MTT (h)' = 0.812 (0.797-0.831)
    lntr <- log(20.2)
    label("Number of absorption transit compartments N (unitless, non-integer)") # Table 3 'N' = 20.2 (16.6-25.6)

    # IIV and IOV are reported as CV% = sqrt(exp(omega^2) - 1) * 100 (Table 3
    # footnote), so omega^2 = log(1 + CV^2).
    etalcl ~ 0.0157513 # Table 3 'IIV CL (CV%)' = 12.6 -> log(1 + 0.126^2)
    etalka ~ 1.22451 # Table 3 'IIV Ka (CV%)' = 155 -> log(1 + 1.55^2)
    etalq ~ 0.0639593 # Table 3 'IIV Q (CV%)' = 25.7 -> log(1 + 0.257^2)
    etalntr ~ 0.0399935 # Table 3 'IIV N (CV%)' = 20.2 -> log(1 + 0.202^2)

    # Inter-occasional variability, one occasion per study period. The
    # occasion-2 variance is fixed equal to the occasion-1 variance (NONMEM
    # $OMEGA BLOCK(1) SAME).
    etaiov_mtt_1 ~ 0.00192536 # Table 3 'IOV MTT (CV%)' = 4.39 -> log(1 + 0.0439^2)
    etaiov_mtt_2 ~ fixed(0.00192536) # Table 3 'IOV MTT (CV%)' = 4.39, occasion 2 (SAME)
    etaiov_ntr_1 ~ 0.0305059 # Table 3 'IOV N (CV%)' = 17.6 -> log(1 + 0.176^2)
    etaiov_ntr_2 ~ fixed(0.0305059) # Table 3 'IOV N (CV%)' = 17.6, occasion 2 (SAME)

    # Proportional residual error (Results 'Model Development'). Table 3
    # footnote: 'CV% for PRV computed as sqrt(exp(sigma^2) - 1)', so
    # sigma^2 = log(1 + 0.194^2) and the proportional SD is its square root.
    propSd <- 0.192211
    label("Proportional residual error (fraction)") # Table 3 'PRV (CV%)' = 19.4 (16.0-22.9) -> sqrt(log(1 + 0.194^2))
  })
  model({
    # Occasion indicators multiplex the per-period IOV etas (OCC = 1 is the
    # cyclosporine-alone period, OCC = 2 the baicalin co-administration period).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2
    iov_ntr <- oc1 * etaiov_ntr_1 + oc2 * etaiov_ntr_2

    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq + etalq)
    ka <- exp(lka + etalka)
    mtt <- exp(lmtt + iov_mtt)
    ntr <- exp(lntr + etalntr + iov_ntr)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Savic 2007 transit-compartment absorption (Methods 'Absorption Model
    # Selection'; Table 2 model M8). transit() is rxode2's analytical form of a
    # chain of ntr transit compartments, ktr = (ntr + 1) / mtt, feeding the
    # absorption compartment, which empties into central at ka. The oral
    # bioavailability is not identifiable (all clearances and volumes are
    # apparent), so the whole dose enters the chain.
    d/dt(depot) <- transit(ntr, mtt, 1) - ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # transit() reads the dose amount from podo(depot) and delivers it through
    # the chain, so the dose record must not also deposit a bolus in depot.
    f(depot) <- 0

    # Dose in mg, volume in L -> mg/L; x 1000 gives ng/mL (whole blood).
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
