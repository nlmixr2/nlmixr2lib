`Buil-Bruna_2015_lanreotide` <- function() {
  description <- "One-compartment population PK model with parallel first- and zero-order subcutaneous absorption for lanreotide Autogel/Depot in patients with gastroenteropancreatic neuroendocrine tumors (Buil-Bruna 2015). A linear effect of body weight on apparent clearance and a small categorical effect of sex on the first-order absorbed fraction are retained; absolute bioavailability F is not identifiable and is structurally anchored at 1, so apparent CL/F and Vd/F are reported. Concentrations are predicted in ng/mL; residual error is additive on the log-transformed observations (LTBS), mapped to proportional in linear space."
  reference <- "Buil-Bruna N, Garrido MJ, Dehez M, Manon A, Nguyen TXQ, Gomez-Panzani EL, Troconiz IF. Population Pharmacokinetic Analysis of Lanreotide Autogel/Depot in the Treatment of Neuroendocrine Tumors: Pooled Analysis of Four Clinical Trials. Clin Pharmacokinet. 2016;55(4):461-473. doi:10.1007/s40262-015-0329-4"
  vignette <- "Buil-Bruna_2015_lanreotide"
  units <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear effect on apparent total serum clearance per Buil-Bruna 2015 Eq. 2: CL/F = theta_CL * [1 + theta_BW * (BW - 74)]. The reference value 74 kg is the population median (Buil-Bruna 2015 Results Sect. 3.3.3). Source column 'BW' in the paper.",
      source_name        = "BW"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female; reference is female because Buil-Bruna 2015 Eq. 1 sets theta_SEX = 0 for females and -0.024 for males, so the published theta_F1 = 0.994 is the female typical value)",
      notes              = "Buil-Bruna 2015 Eq. 1: F1 = theta_F1 * (1 + theta_SEX) with theta_SEX = 0 (females), -0.024 (males). Mapped to canonical SEXF (1 = female). Implemented as F1 = exp(lfdepot + etalfdepot) * (1 + e_sexf_fdepot * (1 - SEXF)), so the published coefficient -0.024 multiplies the male indicator (1 - SEXF) and the typical value lfdepot exactly equals the female value 0.994 from Table 2. Source uses categorical SEX with female as reference (theta_SEX = 0 for females); the canonical SEXF reference category is reversed (male = 0), so the coefficient is carried on the male indicator (1 - SEXF) to preserve source fidelity.",
      source_name        = "SEX"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 290L,
    n_observations  = 1541L,
    n_studies       = 4L,
    age_range       = "pooled mean 60.7 years (CV 18.2%); per-study means 58.5-63.3 years; median 62 years (Buil-Bruna 2015 Eq. 2 reference)",
    weight_range    = "pooled mean 75.1 kg (CV 22.2%); per-study means 69.3-78.0 kg; median 74 kg (Buil-Bruna 2015 Eq. 2 reference); BW also discussed up to 89 kg in subgroup analyses",
    sex_female_pct  = NA_real_,
    race_ethnicity  = "Predominantly White; RACE was not formally tested as a covariate because the majority of patients were White. A 90% prediction interval was retrospectively superimposed on data from 10 Asian and 10 Black/African American patients (Buil-Bruna 2015 Fig. 5b) to confirm no apparent ethnicity effect.",
    disease_state   = "Functioning and non-functioning gastroenteropancreatic neuroendocrine tumors (GEP-NETs). Primary tumor location distribution in the pooled dataset (Table 1): foregut 10, midgut 87, hindgut 13, other 13, unknown 121.",
    dose_range      = "60, 90, or 120 mg lanreotide Autogel/Depot every 4 weeks (deep SC injection); Study 4 dose-titration regimens 60/90/120 mg",
    regions         = "Pooled phase III + phase II studies: CLARINET (multi-regional), ELECT (multi-regional), Study 3 (Spain, 17 centres), Study 4 dose-titration (multi-regional)",
    renal_function  = "Cockcroft-Gault CLCR (mL/min): normal >90, n = 130; mild 60-89, n = 100; moderate 30-59, n = 58; severe <30, n = 2 (Buil-Bruna 2015 Results Sect. 3.2)",
    notes           = "Pooled population PK dataset of 1541 serum lanreotide concentrations from 290 GEP-NET patients across four clinical trials (Study 1 CLARINET n = 101, Study 2 ELECT n = 104, Study 3 n = 30, Study 4 dose-titration n = 71). Serum lanreotide was quantified by a validated radioimmunoassay (LLOQ 0.078 ng/mL; intra-/inter-assay CV 2.3-13.6%; accuracy >89%). Eleven samples from six patients showed anti-lanreotide antibodies at lanreotide concentrations of 4.7-8.19 ng/mL; ADA status had no PK effect (Fig. 5a)."
  )

  ini({
    # Structural PK parameters - Buil-Bruna 2015 Table 2 final-model estimates.
    # CL/F and Vd/F are apparent parameters; absolute bioavailability F was not
    # identifiable and is structurally anchored at 1 (Fig. 2 caption: "F
    # absolute bioavailability (not known and arbitrarily set to 1)"). The
    # paper writes CL/F = 513 L/day, Vd/F = 18.3 L, ka = 1.59e-2 /day, D0 =
    # 2.96 day, theta_F1 = 0.994 (female reference).
    lcl     <- log(513);     label("Apparent total serum clearance CL/F at WT = 74 kg (L/day)") # Buil-Bruna 2015 Table 2: theta_CL = 513 L/day
    lvc     <- log(18.3);    label("Apparent volume of distribution Vd/F (L)")                  # Buil-Bruna 2015 Table 2: Vd/F = 18.3 L
    lka     <- log(0.0159);  label("First-order absorption rate ka (1/day)")                    # Buil-Bruna 2015 Table 2: ka = 1.59e-2 /day
    lfdepot <- log(0.994);   label("First-order absorbed fraction F1 (female reference; unitless)") # Buil-Bruna 2015 Table 2: theta_F1 = 0.994
    lD0     <- log(2.96);    label("Zero-order absorption duration D0 (day)")                   # Buil-Bruna 2015 Table 2: D0 = 2.96 day

    # Covariate effects - Buil-Bruna 2015 Table 2.
    # Body weight on CL/F: linear deviation from the median (Eq. 2)
    #   CL/F = theta_CL * [1 + theta_BW * (BW - 74)]
    # where 74 kg is the population median. theta_BW = 9.77e-3 /kg (95% CI of
    # the magnitude printed in Table 2 with sign ambiguity; magnitude and sign
    # consistent with the +14.7% CL/F at 89 kg vs 74 kg reported in Sect. 3.5).
    e_wt_cl       <-  9.77e-3; label("Linear coefficient of (WT - 74 kg) on CL/F (per kg)")     # Buil-Bruna 2015 Table 2: theta_BW = 9.77e-3 /kg

    # Sex on F1: paper Eq. 1
    #   F1 = theta_F1 * (1 + theta_SEX)
    # with theta_SEX = 0 (females) or -0.024 (males). Carried on the male
    # indicator (1 - SEXF) so the typical value lfdepot stays exactly at the
    # published female value 0.994.
    e_sexf_fdepot <- -0.024;   label("Linear coefficient of male indicator (1 - SEXF) on F1 (unitless)") # Buil-Bruna 2015 Table 2: theta_SEX (males) = -0.024

    # Inter-patient variability (IPV). Table 2 reports IPV as CV%; convert to
    # log-normal variance via omega^2 = log(1 + CV^2). Paper Methods Sect. 2.4:
    # "IPV was modeled exponentially" (theta_i = theta * exp(eta_i)).
    #   IPV CL/F = 27%  -> log(1 + 0.27^2)   = 0.07031
    #   IPV Vd/F = 150% -> log(1 + 1.50^2)   = 1.17865
    #   IPV ka   = 61%  -> log(1 + 0.61^2)   = 0.31641
    #   IPV F1   = 1.05%-> log(1 + 0.0105^2) = 0.0001103
    # The non-diagonal elements of OMEGA were tested and did not improve fit
    # significantly (Sect. 3.3.1), so IPVs are independent here.
    etalcl     ~ 0.07031     # Buil-Bruna 2015 Table 2: IPV CL/F = 27%
    etalvc     ~ 1.17865     # Buil-Bruna 2015 Table 2: IPV Vd/F = 150%
    etalka     ~ 0.31641     # Buil-Bruna 2015 Table 2: IPV ka = 61%
    etalfdepot ~ 0.0001103   # Buil-Bruna 2015 Table 2: IPV F1 = 1.05% (variability on log-F1; the paper's exponential parameterization can numerically produce F1 > 1 for a small fraction of subjects given F1 typical = 0.994; see vignette Assumptions and deviations)

    # Residual error - Buil-Bruna 2015 Table 2: 0.275 [log(ng/mL)]. The paper
    # Methods Sect. 2.4 state: "Data were transformed logarithmically [...]
    # Residual variability was modeled considering an additive error in the
    # logarithmic domain of the transformed data." This log-transform-both-
    # sides (LTBS) parameterization, Y = log(F) + EPS(1) in NONMEM, maps to
    # a proportional error in linear-space nlmixr2 per the standard
    # NONMEM-to-nlmixr2 translation: Cc ~ prop(propSd) with propSd = sigma_log.
    propSd <- 0.275; label("LTBS-equivalent proportional residual SD (fraction)") # Buil-Bruna 2015 Table 2: sigma_log = 0.275 log(ng/mL)
  })

  model({
    # Individual PK parameters. CL/F has a linear deviation from the median
    # body weight (74 kg) per Buil-Bruna 2015 Eq. 2. The bioavailability for
    # the first-order absorption pathway, F1, is the female typical value
    # adjusted downward by 2.4% in males (Eq. 1) and individualized via
    # etalfdepot.
    cl  <- exp(lcl + etalcl) * (1 + e_wt_cl * (WT - 74))
    vc  <- exp(lvc + etalvc)
    ka  <- exp(lka + etalka)
    F1  <- exp(lfdepot + etalfdepot) * (1 + e_sexf_fdepot * (1 - SEXF))
    D0  <- exp(lD0)
    kel <- cl / vc

    # Parallel first- and zero-order subcutaneous absorption (Buil-Bruna 2015
    # Fig. 2 + Methods Sect. 2.4.2). The total prescribed dose enters the
    # 'depot' bookkeeping compartment; internally it splits into two arms:
    #   * the F1 fraction follows first-order absorption with rate ka, and
    #     contributes ka * depot to central;
    #   * the F2 = 1 - F1 fraction is delivered into 'central' at a constant
    #     rate (1 - F1) * Dose / D0 over the first D0 days after each dose.
    # The default-and-conditional-upgrade pattern for kzero (rather than
    # compute-then-zero) keeps kzero finite before any depot dose has been
    # given, because podo(depot) returns NA in that pre-dose window and an
    # NA seed would propagate through the ODE solver.
    kzero <- 0.0
    if (tad(depot) <= D0) kzero <- (1 - F1) * podo(depot) / D0

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot + kzero - kel * central

    # Route only the F1 fraction into the depot compartment; the remaining
    # (1 - F1) fraction reaches central via kzero, avoiding the need for a
    # second user-side dose record per administration.
    f(depot) <- F1

    # Concentration: dose in mg and vc in L give central / vc in mg/L =
    # ug/mL. Multiply by 1000 to express Cc in ng/mL, matching Buil-Bruna
    # 2015 Table 2 (sigma in log(ng/mL)) and the y-axis units of Fig. 1.
    Cc <- (central / vc) * 1000
    Cc ~ prop(propSd)
  })
}
