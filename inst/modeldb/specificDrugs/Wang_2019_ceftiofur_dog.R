Wang_2019_ceftiofur_dog <- function() {
  description <- paste(
    "Preclinical (beagle dog).",
    "Two-compartment mammillary population PK model for ceftiofur",
    "equivalents (ceftiofur plus its desfuroylceftiofur metabolites that",
    "retain an intact beta-lactam ring, assayed together after",
    "derivatisation to desfuroylceftiofur acetamide) after a single",
    "2.2 mg/kg dose of ceftiofur sodium given intravenously or",
    "subcutaneously to healthy beagles. First-order elimination from the",
    "central compartment and first-order absorption from the subcutaneous",
    "depot; absolute bioavailability is estimated because the same animals",
    "supplied both routes. Female sex lowers the subcutaneous absorption",
    "rate roughly two-fold. Inter-individual variability on clearance and",
    "central volume is almost perfectly correlated (r = 0.999); the",
    "authors' search drove the absorption-rate and inter-compartmental",
    "clearance variances to zero and they were fixed there. Fitted by SAEM",
    "in Monolix 2018R2 with M4 handling of below-limit-of-quantification",
    "records.",
    "UNITS: the model is coded in ABSOLUTE units (L, L/h) with the dose in",
    "mg, not per kilogram. Wang 2019 Table 1 labels the disposition",
    "parameters 'L/kg' and 'L/h/kg', but that suffix is falsified by the",
    "paper's own simulations and figures -- see the vignette's",
    "'Assumptions and deviations' section. A 2.2 mg/kg dose in a 10 kg",
    "beagle is amt = 22 mg.",
    sep = " "
  )
  reference <- paste(
    "Wang J, Schneider BK, Xue J, Sun P, Qiu J, Mochel JP, Cao X.",
    "Pharmacokinetic Modeling of Ceftiofur Sodium Using Non-linear",
    "Mixed-Effects in Healthy Beagle Dogs.",
    "Front Vet Sci. 2019;6:363.",
    "doi:10.3389/fvets.2019.00363.",
    sep = " "
  )
  vignette <- "Wang_2019_ceftiofur_dog"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ug/mL"
  )

  covariateData <- list(
    SEXF = list(
      description = "Biological sex of the dog; 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Retained on the subcutaneous absorption rate only (Wang 2019",
        "Results, 'Pharmacokinetic Model'). The source writes the",
        "relationship as log(ka_i) = log(ka_pop) + beta * sex_{i=f} + eta_i",
        "with sex_{i=f} = 1 for a female, so ka_pop is the MALE typical",
        "value and the canonical SEXF orientation needs no value",
        "transformation and no sign flip. Six males and six females were",
        "enrolled and the route assignment was blocked on sex (3 of each",
        "per route). Body weight and age were also screened by the",
        "automated Monolix search and were not retained; they are recorded",
        "in covariatesDataExcluded."
      ),
      source_name = "sex"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Screened by the automated Pearson correlation test in Monolix",
        "2018R2 after median-normalisation and log-transformation (Wang",
        "2019 Methods, 'Inclusion of Covariate Relationships') and not",
        "retained at P < 0.05. The cohort spanned only 9-12 kg, so the",
        "study had little power to resolve a size effect."
      )
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = paste(
        "Screened alongside body weight after median-normalisation and",
        "log-transformation and not retained. The cohort spanned only",
        "1.5-2.5 years."
      )
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "ceftiofur",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "ceftiofur equivalents",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "ceftiofur equivalents",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "beagle dog",
    n_subjects = 12L,
    n_studies = 1L,
    age_range = "1.5-2.5 years",
    weight_range = "9-12 kg",
    sex_female_pct = 50,
    disease_state = "Healthy (screened by physical examination, haematology, clinical chemistry and coagulation time)",
    dose_range = paste(
      "Single 2.2 mg/kg dose of ceftiofur sodium, reconstituted from",
      "sterile powder in 20 mL bacteriostatic water per 1 g vial.",
      "Six dogs received it intravenously (cephalic vein) and six",
      "subcutaneously (behind the shoulders), assigned by a block design",
      "on sex so that 3 males and 3 females went to each route."
    ),
    regions = "China (China Agricultural University, Beijing)",
    notes = paste(
      "Plasma sampled at 0, 0.08 (intravenous arm only), 0.25, 0.5, 0.75,",
      "1, 1.5, 2, 3, 4, 6, 8, 12, 24, 36, 48 and 72 h post dose; 198",
      "concentrations from both routes were pooled and fitted",
      "simultaneously. Ceftiofur and every desfuroylceftiofur metabolite",
      "retaining an intact beta-lactam ring were cleaved with",
      "dithioerythritol, derivatised with iodoacetamide to",
      "desfuroylceftiofur acetamide and quantified by UPLC-MS/MS; the",
      "model therefore describes TOTAL ceftiofur equivalents, not parent",
      "ceftiofur, and free drug is roughly 10% of that total. LLOQ",
      "100 ng/mL, calibration range 100-5000 ng/mL; 8 of 198 records",
      "(4.0%) were below it and were handled by the M4 likelihood",
      "method. See Wang 2019 Methods, 'Drug Supply and Animals' through",
      "'Handling of Below Limit of Quantification (BLQ) Data'."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Wang 2019 Table 1, 'Point estimate' column.
    #
    # UNITS. Table 1's unit column reads 'L/h/kg' and 'L/kg'. Those values
    # are used here as ABSOLUTE L/h and L for a typical study beagle, with
    # the dose supplied in mg. The '/kg' suffix cannot be right: read as
    # per-kg, a 2.2 mg/kg subcutaneous dose peaks at 0.86 ug/mL and stays
    # above the MIC50 of 0.5 ug/mL for only 4.9 h, whereas the Abstract
    # reports ~30 h, Table 2 reports target attainment out to an MIC of
    # 4 ug/mL, and Figure 2's concentration axis runs to 20,000 ng/mL.
    # Read as absolute, the same numbers reproduce every one of those.
    # The numeric values below are exactly as printed; only the unit label
    # is treated as a typographic error. See the vignette's 'Assumptions
    # and deviations' section for the full arbitration.
    # ------------------------------------------------------------------
    lka <- log(1.43)
    label("Absorption rate constant after subcutaneous dosing, MALE reference (log 1/h)") # Table 1, Absorption (S.C): Ka = 1.43 1/h (RSE 11.9%)
    lcl <- log(0.25)
    label("Clearance (log L/h)") # Table 1, Clearance: CL = 0.25 (RSE 8.29%)
    lvc <- log(1.69)
    label("Central volume of distribution (log L)") # Table 1, Central compartment volume of distribution: V1 = 1.69 (RSE 6.9%)
    lvp <- log(1.28)
    label("Peripheral volume of distribution (log L)") # Table 1, Peripheral compartment volume of distribution: V2 = 1.28 (RSE 12.9%)
    lq <- log(0.16)
    label("Inter-compartmental clearance (log L/h)") # Table 1, Inter-compartmental clearance: Q = 0.16 (RSE 13.6%)
    lfdepot <- log(0.937)
    label("Absolute bioavailability after subcutaneous dosing (log fraction)") # Table 1, Bioavailability (S.C): F = 93.7% (RSE 11.4%)

    # ------------------------------------------------------------------
    # Covariate effect -- Wang 2019 Results, 'Pharmacokinetic Model':
    #   log(ka_i) = log(ka_pop) + beta * sex_{i=f} + eta_i
    # with sex_{i=f} = 1 for a female and 0 otherwise, so the effect is
    # additive on the log scale and ka_pop above is the male value. The
    # implied male:female ratio is exp(0.643) = 1.90, matching the paper's
    # 'two times greater in male vs. female dogs'.
    # ------------------------------------------------------------------
    e_sexf_ka <- -0.643
    label("Effect of female sex on log absorption rate (additive on log scale)") # Table 1, Coefficient (Ka and sex): beta_sex = -0.643 (RSE 20.1%)

    # ------------------------------------------------------------------
    # IIV -- Wang 2019 Table 1, 'IIV (%)' column. The table footnote
    # states the column is expressed as CV%, and Methods writes every
    # individual parameter as phi_i = mu * exp(eta_i), so the log-scale
    # variance is omega^2 = log(1 + CV^2):
    #   CL  24.0% -> log(1 + 0.240^2) = 0.05600219
    #   V1  32.4% -> log(1 + 0.324^2) = 0.09982362
    #   V2  25.7% -> log(1 + 0.257^2) = 0.06395929
    #   F   52.0% -> log(1 + 0.520^2) = 0.23933181
    #
    # Ka and Q carry no eta: the Table 1 footnote reads 'Model parameter
    # estimated to converge to a null value and fixed to 0', so those two
    # variances were driven to zero during estimation and held there. They
    # are omitted rather than written as fixed(0), which would make the
    # omega matrix singular and break rxode2 simulation.
    #
    # The CL-V1 block reproduces Table 1's corr(cl_v1) = 99.9%, which
    # Supplemental Figure 1C confirms as a correlation of the random
    # effects themselves. The off-diagonal is
    # 0.999 * sqrt(0.05600219 * 0.09982362) = 0.07469381. The block is
    # positive definite but close to singular (eigenvalues 0.1558 and
    # 7.18e-05).
    # ------------------------------------------------------------------
    etalcl + etalvc ~ c(
      0.05600219,
      0.07469381, 0.09982362
    ) # Table 1: IIV CL 24%, IIV V1 32.4%, corr(cl_v1) 99.9% (RSE 6.24%)
    etalvp ~ 0.06395929 # Table 1: IIV V2 25.7%
    etalfdepot ~ 0.23933181 # Table 1: IIV F 52%

    # ------------------------------------------------------------------
    # Residual error -- Wang 2019 Results: 'A log-normal error model best
    # captured the residual variability in the model'. The MAGNITUDE is
    # not printed anywhere in the paper: Table 1 has no residual-error
    # row, and Supplemental Figure 1 shows only IWRES scatter, eta
    # boxplots and the random-effect correlation matrix.
    #
    # NON-PAPER PROVENANCE -- digitised from Figure 2. That panel plots
    # observations against individual predictions with 'dotted black
    # lines: 90% prediction interval'. For a log-normal residual the
    # bands sit at f * exp(+/- 1.6449 * expSd), so the band-to-identity
    # ratio reads the parameter off directly. Measured on the 300 dpi
    # render of the left (I.V) panel over 239 independent pixel columns:
    # a constant 37.5 px offset on a 257.5 px/decade axis = 0.1456
    # decades = a ratio of 1.3984, giving expSd = ln(1.3984) / 1.6449 =
    # 0.204 (about 20.6% CV). The offset being CONSTANT across the whole
    # concentration range is itself the confirmation that the residual is
    # pure log-scale with no additive component. Treat as +/- 0.01.
    # ------------------------------------------------------------------
    expSd <- 0.204
    label("Log-normal residual error, SD on the log scale") # Figure 2, 90% prediction-interval bands (digitised; not printed in Table 1)
  })

  model({
    # 1. Individual parameters. Log-normal IIV per Wang 2019 Methods,
    #    phi_i = mu * exp(eta_i). Ka and Q have no eta (fixed to 0).
    ka <- exp(lka + e_sexf_ka * SEXF)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q <- exp(lq)

    # 2. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. Two-compartment mammillary system with a first-order
    #    subcutaneous depot -- Wang 2019 Figure 1.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 4. Absolute bioavailability applies to the subcutaneous depot only;
    #    intravenous doses go straight into central and carry F = 1.
    f(depot) <- exp(lfdepot + etalfdepot)

    # 5. Observation. Amount in mg over volume in L gives mg/L, which is
    #    ug/mL -- the units the paper's MIC thresholds use (Figure 2 plots
    #    the same quantity as ng/mL, i.e. 1000-fold larger numbers).
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
