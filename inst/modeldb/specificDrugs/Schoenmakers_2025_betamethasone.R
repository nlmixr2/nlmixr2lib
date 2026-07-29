Schoenmakers_2025_betamethasone <- function() {
  description <- "Two-compartment population PK model with first-order absorption (no lag time) for intramuscular betamethasone in pregnant women admitted with imminent preterm birth, including early-onset pre-eclampsia (eoPE; diagnosed before 34 weeks gestation). Apparent clearance is multiplied by 0.617 (a 38% reduction, or ~60% of the non-eoPE clearance) when eoPE is present; this is the only retained covariate after backward elimination at P < 0.01. Body weight, BMI, lean body weight, age, gestational age, number of foetuses, white blood cell counts and CRP were screened but did not retain after backward elimination."
  reference <- paste(
    "Schoenmakers S, Li L, Kluivers ACM, Broekhuizen M, Harhangi MS, Ronde E,",
    "DeKoninck PLJ, Reiss I, Danser AHJ, Allegaert K, van den Berg SAA, van Zelst BD,",
    "van Schaik RHN, Simons SHP, Koch BCP, Sassen SDT. (2025).",
    "Pharmacokinetics of betamethasone in pre-eclampsia: An in vivo and ex vivo study.",
    "Br J Clin Pharmacol 91(11):2327-2339. doi:10.1002/bcp.70035.",
    sep = " "
  )
  vignette <- "Schoenmakers_2025_betamethasone"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    DIS_EOPE = list(
      description        = "Binary indicator of early-onset pre-eclampsia (eoPE; diagnosed before 34 weeks gestation); 1 = eoPE present, 0 = not eoPE.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no eoPE; pregnant woman with imminent preterm birth and no pre-eclampsia at the time of betamethasone administration)",
      notes              = "Time-fixed per subject. The Schoenmakers 2025 final model multiplies CL/F by ThetaPE = 0.617 (Eq 2 with the categorical-indicator form) when eoPE is present, giving a 38% reduction in apparent clearance (non-eoPE CL/F = 15.6 L/h; eoPE CL/F = 15.6 x 0.617 ~= 9.6 L/h, matching the abstract's reported 9.35 L/h vs 15.78 L/h median individual estimates). The source dataset cohort was 23 non-eoPE women + 5 eoPE women drawn from a single Dutch obstetric centre between January and October 2021.",
      source_name        = "DIS_EOPE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 28L,
    n_studies      = 1L,
    n_observations = 194L,
    age_range      = "21-41 years",
    age_median     = "31 years",
    weight_range   = "50-121 kg",
    weight_median  = "68 kg",
    sex_female_pct = 100,
    race_ethnicity = NULL,
    ga_range       = "23+5 to 33+6 weeks (gestational age at first betamethasone dose)",
    ga_median      = "27+4 weeks",
    disease_state  = paste(
      "Pregnant women admitted with imminent preterm birth and treated with",
      "intramuscular betamethasone for foetal lung maturation. Subgroups:",
      "23 women without pre-eclampsia and 5 women with early-onset",
      "pre-eclampsia (eoPE; pre-eclampsia diagnosed < 34 weeks gestation)."
    ),
    dose_range     = paste(
      "Celestone Chronodose 11.4 mg intramuscular betamethasone (6 mg betamethasone",
      "sodium phosphate + 5.4 mg betamethasone acetate, an instant-release plus",
      "depot-prodrug combination) given once daily for 2 consecutive days (2 doses",
      "total). The model lumps the two prodrug forms into a single first-order",
      "absorption Ka because active betamethasone, not the prodrug esters, was",
      "measured."
    ),
    regions        = "Netherlands (Erasmus University Medical Center, Rotterdam; single-centre study, MEC-2019-0650)",
    notes          = paste(
      "Demographics from Schoenmakers 2025 Table 1 (n = 28 women contributing",
      "194 maternal serum samples; samples taken pre-dose and at 0-0.5, 1-3,",
      "5-8, 10-12, and 20-24 h after each of the two doses). 8 samples were",
      "excluded for abnormal concentrations (potential sampling errors / time",
      "recording errors). Two of the 30 initially enrolled women withdrew",
      "(Figure S1). Twin pregnancies n = 2 (both in the non-eoPE subgroup);",
      "singleton n = 26. Active betamethasone serum concentrations quantified",
      "by LC-MS/MS. Race / ethnicity not reported. The model file encodes the",
      "final-model parameters from Table 2 ('Final model' column with SIR-",
      "confirmed point estimates); the in-file source-trace comments cite the",
      "specific Table 2 row for each value. The two betamethasone metabolites",
      "(11-keto-betamethasone and 6-beta-hydroxybetamethasone) reported in",
      "Table S3 are NOT modelled here; the supplement is not on disk and the",
      "main paper only describes the parent betamethasone fit. See vignette",
      "Assumptions and deviations for the IIV / residual-error reporting-",
      "convention interpretation."
    )
  )

  ini({
    # Structural parameters - typical values reference a pregnant woman
    # WITHOUT early-onset pre-eclampsia (DIS_EOPE = 0). The Schoenmakers 2025
    # final model fitted a 2-compartment open model with first-order
    # absorption (no lag time) and first-order elimination, with a
    # proportional residual error and IIV on CL and Vc (Methods 2.2 and
    # Results 3.1).
    lka  <- log(1.67);  label("First-order absorption rate constant after IM betamethasone (1/h)")           # Schoenmakers 2025 Table 2 final, Ka = 1.67 (RSE 27%)
    lcl  <- log(15.6);  label("Apparent clearance for a non-eoPE pregnant woman (CL/F, L/h)")                  # Schoenmakers 2025 Table 2 final, CL/F = 15.6 (RSE 5%)
    lvc  <- log(46.1);  label("Apparent central volume of distribution (Vc/F, L)")                            # Schoenmakers 2025 Table 2 final, Vc/F = 46.1 (RSE 34%)
    lvp  <- log(109);   label("Apparent peripheral volume of distribution (Vp/F, L)")                         # Schoenmakers 2025 Table 2 final, Vp/F = 109 (RSE 7%)
    lq   <- log(99.8);  label("Apparent inter-compartmental clearance (Q/F, L/h)")                            # Schoenmakers 2025 Table 2 final, Q/F = 99.8 (RSE 23%); Table 2 unit label is 'L' but the row is for an inter-compartmental clearance, so units are L/h (confirmed by CL/F row units L/h above and by the parameter's role in the rate equation).

    # Covariate effect: early-onset pre-eclampsia (eoPE) on CL.
    # The source paper's Eq. 2 (categorical covariate form):
    #   theta_j = theta_jTPV * (theta_COVi)^FLAG_COVi
    # is implemented in nlmixr2 as a log-additive shift on lcl:
    #   cl = exp(lcl + etalcl + e_eope_cl * DIS_EOPE)
    # so DIS_EOPE = 0 yields the typical-value CL, and DIS_EOPE = 1 multiplies CL
    # by exp(log(0.617)) = 0.617 (a 38% reduction).
    e_eope_cl <- log(0.617);  label("Log-multiplicative effect of early-onset pre-eclampsia on CL/F")         # Schoenmakers 2025 Table 2 final, ThetaPE = 0.617 (RSE 12%; SIR 95% CI 0.28-0.76); paper Methods Eq 2.

    # IIV on the log scale. The source paper reports IIV-CL = 21.3% and
    # IIV-Vc = 89.5% in Table 2 as a percentage without an explicit
    # conversion formula. We interpret these as approximate CV% on the
    # natural (linear) scale and back-transform via the exact log-normal
    # relation
    #   omega^2 = log(1 + CV^2)
    # so that the etas are normal(0, omega^2) and the parameter is
    # exp(eta) times the typical value. See vignette Assumptions and
    # deviations for the alternative sqrt(omega^2)*100 reading.
    # CV-CL = 21.3% -> omega^2 = log(1 + 0.213^2) = 0.0444
    # CV-Vc = 89.5% -> omega^2 = log(1 + 0.895^2) = 0.5887
    etalcl ~ 0.0444  # Schoenmakers 2025 Table 2 final, IIV-CL = 21.30% (RSE 27%, shrinkage 24%); omega^2 = log(1 + 0.213^2).
    etalvc ~ 0.5887  # Schoenmakers 2025 Table 2 final, IIV-Vc = 89.50% (RSE 30%, shrinkage 20%); omega^2 = log(1 + 0.895^2).

    # Proportional residual error. The source paper reports
    # 'Residual variability = 0.224' in Table 2 without a unit annotation
    # (whereas IIV values carry an explicit '%' sign). We interpret 0.224
    # as the SD of the proportional residual (i.e., propSd = 22.4% CV),
    # consistent with conventional NONMEM popPK reporting of
    # log-normal proportional error. See vignette Assumptions and
    # deviations for the alternative variance reading.
    propSd <- 0.224;  label("Proportional residual SD (fraction)")                                            # Schoenmakers 2025 Table 2 final, residual variability = 0.224 (RSE 11%); interpreted as SD of the proportional residual.
  })

  model({
    # Individual PK parameters. The Schoenmakers 2025 final model retains
    # only the eoPE indicator on CL after backward elimination at
    # P < 0.01 (Results 3.1: 'keeping the weight-related covariates in
    # the model did not result in a significantly better model and
    # therefore was not added to the model'). Vc carries IIV but no
    # covariates; Vp, Q and Ka carry neither covariates nor IIV in the
    # published final model.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl + e_eope_cl * DIS_EOPE)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q  <- exp(lq)

    # Two-compartment open model with first-order absorption and
    # first-order elimination, no lag time. Bioavailability is not
    # identifiable from IM-only data, so all volumes and clearances
    # are apparent (X/F). The depot compartment absorbs the lumped
    # betamethasone-sodium-phosphate + betamethasone-acetate prodrug
    # mixture (Methods 2.2: 'we were unable to model the prodrug esters
    # independently').
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Plasma concentration of active betamethasone. Units: dose in mg,
    # volumes in L -> concentration in mg/L (equivalent to ug/mL). The
    # source paper plots concentrations in ug/L (ng/mL), which is the
    # same value as mg/L * 1000.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
