Yang_2026_lurasidone <- function() {
  description <- "One-compartment population PK model with first-order absorption for lurasidone in Chinese psychiatric inpatients, with a linear age effect and a concomitant-valproate effect on apparent clearance"
  reference <- paste(
    "Yang Y, Lu H, Xiao T, Ni X, Wang Z, Chen Y, Dai L, Song E, Su F and Wen Y (2026).",
    "Impact of valproate co-medication and age on lurasidone exposure:",
    "a population pharmacokinetic study and real-world evaluation in Chinese psychiatric inpatients.",
    "Front Pharmacol 17:1810528. doi:10.3389/fphar.2026.1810528",
    sep = " "
  )
  vignette <- "Yang_2026_lurasidone"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    AGE = list(
      description        = "Age at the pharmacokinetic observation",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters CL/F as a LINEAR centred effect, not the power model that Yang 2026 Section 2.5",
        "describes generically for continuous covariates: the final model equation printed in",
        "Section 3.3 is CL/F = 339 * [1 - 0.0125 * (AGE - 22)] * (1 + 0.477 * VPA) * exp(eta).",
        "The centring value of 22 years is the cohort median age (Table 1: median 22.00,",
        "range 13.00-70.00), which is a free transcription check on the constant.",
        "CAUTION: a linear (rather than power or exponential) covariate form is not bounded below.",
        "The bracketed multiplier reaches zero at AGE = 102 years and turns negative above it,",
        "which would give a non-physical negative clearance. The model is only interpretable over",
        "the cohort age range of 13-70 years; the paper itself flags that the elderly arm rests on",
        "n = 4 subjects aged >= 65 and is 'hypothesis-generating' (Section 4.5). Simulations",
        "outside 13-70 years should be avoided rather than clamped, since clamping would depart",
        "from the published model."
      ),
      source_name        = "AGE"
    ),
    CONMED_VPA = list(
      description        = "Concomitant valproate (valproic acid / sodium valproate) co-medication",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant valproate)",
      notes              = paste(
        "Yang 2026 Section 3.3: 'VPA denotes valproate comedication status (binary categorical",
        "variable): 0 = no comedication, 1 = comedication.' Enters CL/F multiplicatively as",
        "(1 + 0.477 * CONMED_VPA), i.e. a 47.7% INCREASE in apparent clearance and therefore",
        "lower total serum concentrations. 59/156 patients (28%) were on valproate (Table 1).",
        "The authors justify the binary (rather than concentration-driven) encoding with a",
        "Spearman correlation of the final-model CL/F ETAs against paired valproate serum",
        "concentrations in a 39-observation subset (r = 0.12, p = 0.519; Supplementary Figure S1).",
        "Direction note: the increase is counter-intuitive for a drug classically described as a",
        "CYP inhibitor. Yang 2026 Section 4.3 proposes PXR-mediated CYP3A4/MDR1 induction and,",
        "more importantly for interpretation, protein-binding displacement -- lurasidone is >99%",
        "and valproate ~90% protein bound, so displacement raises the free fraction and hence the",
        "APPARENT clearance of the TOTAL drug that the assay measures. Users should therefore not",
        "read this coefficient as a change in unbound exposure."
      ),
      source_name        = "VPA"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "lurasidone", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "lurasidone", units = "mg",
      specimen = "serum", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 156,
    n_observations = 212,
    n_studies      = 1,
    age_range      = "13-70 years",
    age_median     = "22 years",
    weight_range   = "36-138 kg",
    weight_median  = "65.5 kg",
    sex_female_pct = 69.9,
    race_ethnicity = c(Asian = 100),
    disease_state  = "hospitalised psychiatric patients (ICD-10 diagnoses; schizophrenia and bipolar disorder)",
    dose_range     = "20-120 mg once daily oral (median 60 mg/day)",
    regions        = "China (Guangdong)",
    age_groups     = "60 adolescents (13-17 years), 92 adults (18-64 years), 4 elderly (>= 65 years)",
    co_medication  = "valproate 59/156 (28%); lithium carbonate 93/156 (44%)",
    notes          = paste(
      "Retrospective analysis of routine therapeutic drug monitoring (TDM) data from Han Chinese",
      "psychiatric inpatients at The Affiliated Brain Hospital of Guangzhou Medical University.",
      "Baseline demographics are Yang 2026 Table 1; sex is reported there as 47 male / 109 female",
      "(69.9% female; the Results text rounds this to 69%). Sampling was sparse and predominantly",
      "steady-state troughs, which is why Ka could not be estimated and why the structural model",
      "was restricted to one compartment. Observed concentrations had a median of 6.30 ng/mL",
      "(range 1.00-44.53), largely BELOW the AGNP therapeutic reference range of 15-40 ng/mL;",
      "the authors attribute this chiefly to real-world meals falling short of the >= 350 kcal",
      "needed for full lurasidone bioavailability (Section 4.4) rather than to intrinsic",
      "hyper-clearance. Patients with severe hepatic impairment (Child-Pugh >= 10) or severe",
      "renal impairment (CrCL <= 30 mL/min) were excluded."
    )
  )

  ini({
    # Structural parameters -- Yang 2026 Table 3 (Final model, Estimated column) and the
    # final model equation printed in Section 3.3.
    lka <- fixed(log(0.679)); label("Absorption rate constant (Ka, 1/h)")                                              # Table 3 "Ka (h-1) 0.679,FIX"; fixed from Hu 2017 healthy-Chinese PK because the sparse trough design could not identify it (Sections 2.5, 3.2)
    lcl <- log(339); label("Apparent clearance at the reference age of 22 years without valproate (CL/F, L/h)")        # Table 3 CL/F = 339 (RSD 7%; bootstrap median 337, 95% CI 283-394)
    lvc <- log(13600); label("Apparent volume of distribution (V/F, L)")                                               # Table 3 V/F = 13600 (RSD 20%; bootstrap median 13531, 95% CI 8802-19807)

    # Covariate effects on CL/F -- both enter multiplicatively; see model() for the exact form.
    e_age_cl <- 0.0125; label("Fractional decrease in CL/F per year of age above the 22-year reference (1/year)")      # Table 3 theta CL-AGE = 0.0125 (RSD 16%; bootstrap median 0.0122, 95% CI 0.0062-0.0172). Table 3 reports the positive magnitude; the minus sign is carried by the Section 3.3 equation.
    e_conmed_vpa_cl <- 0.477; label("Fractional increase in CL/F with concomitant valproate (fraction)")               # Table 3 theta CL-VPA = 0.477 (RSD 29%; bootstrap median 0.478, 95% CI 0.109-0.880) = the +47.7% clearance increase quoted in the Abstract and Section 3.3

    # IIV -- exponential on CL/F and V/F (Section 2.5). Table 3 "Random effects" block reports
    # NONMEM OMEGA variances on the log scale; CV% = sqrt(exp(omega^2) - 1).
    etalcl ~ 0.505                                                                                                    # Table 3 random effect on CL/F = 0.505 (RSD 20%, eta-shrinkage 19%; bootstrap median 0.480, 95% CI 0.307-0.673) -> 81% CV
    etalvc ~ 1.24                                                                                                     # Table 3 random effect on V/F = 1.24 (RSD 19%, eta-shrinkage 41%; bootstrap median 1.259, 95% CI 0.539-1.974) -> 157% CV

    # Residual error -- exponential/log-normal (Section 2.5: "RUV was evaluated using exponential
    # error models"). Table 3 reports the SIGMA variance 0.0779; the SD is its square root.
    expSd <- sqrt(0.0779); label("Log-normal residual error (log-scale SD)")                                           # Table 3 "Residual error 0.0779" (RSD 24%, epsilon-shrinkage 39%; bootstrap median 0.0757, 95% CI 0.0411-0.1220); sqrt(0.0779) = 0.2791
  })

  model({
    # Individual parameters. Yang 2026 Section 3.3 prints the final model as
    #   CL/F = 339 * [1 - 0.0125 * (AGE - 22)] * (1 + 0.477 * VPA) * exp(eta_i)
    #   V/F  = 13600
    #   Ka   = 0.679 /h
    # The age term is LINEAR and centred on the cohort median age of 22 years (Table 1) -- it is
    # NOT the median-normalised power model that Section 2.5 describes generically for continuous
    # covariates. Where the generic Methods text and the printed final equation disagree, the
    # printed equation governs.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      (1 - e_age_cl * (AGE - 22)) *
      (1 + e_conmed_vpa_cl * CONMED_VPA)
    # The printed V/F line carries no eta, but Section 2.5 states IIV on both CL/F and V/F was
    # modelled with exponential error models, and Table 3 reports the V/F random effect (1.24)
    # with a 41% eta-shrinkage -- a shrinkage figure only exists for an estimated eta.
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # `central` is in mg and `vc` in L, so central/vc is mg/L; x 1000 converts to ng/mL, the
    # unit in which the TDM assay and every concentration in the paper are reported.
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd)
  })
}
