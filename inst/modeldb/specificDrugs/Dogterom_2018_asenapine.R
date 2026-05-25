Dogterom_2018_asenapine <- function() {
  description <- "Two-compartment population PK model with first-order sublingual absorption for asenapine in pediatric patients (aged 10-17 years) with schizophrenia, bipolar I disorder, or other psychiatric disorders (Dogterom 2018 Drug Design, Development and Therapy). Central / peripheral volumes and absorption-rate constant were fixed from a Phase I-only fit; no intrinsic covariate (age, BMI, race, sex) was retained in the final model. Residual error switches per observation between intensive Phase I PK sampling (27.8% CV) and sparse Phase III efficacy sampling (56.0% CV), with an additional between-subject scaling of the residual SD (19.2% CV)."
  reference <- paste(
    "Dogterom P, Riesenberg R, de Greef R, Dennie J, Johnson M, Pilla Reddy V,",
    "Miltenburg AMM, Findling RL, Jakate A, Carrothers TJ, Troyer MD.",
    "Asenapine pharmacokinetics and tolerability in a pediatric population.",
    "Drug Des Devel Ther. 2018;12:2677-2693.",
    "doi:10.2147/DDDT.S171475.",
    sep = " "
  )
  vignette <- "Dogterom_2018_asenapine"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    SAMPLE_INTENSIVE = list(
      description        = "Per-observation indicator of sampling intensity: 1 = intensive (rich post-dose Phase I PK profile), 0 = sparse (Phase III efficacy-study population sample).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (sparse Phase III efficacy sampling)",
      notes              = "Record-level indicator that switches the proportional residual-error magnitude per observation: intensive Phase I sampling 27.8% CV, sparse Phase III sampling 56.0% CV (Dogterom 2018 Table 3 final model). The pooled analysis carried 2,451 observations from 561 pediatric patients across two Phase I PK studies (A7501022 and P06522, intensive sampling between 15 min and 72 h postdose) and two Phase III efficacy studies (P05896 schizophrenia, P06107 bipolar I disorder, sparse population samples). Set SAMPLE_INTENSIVE = 1 on observations from a Phase I PK profile, 0 on observations from a Phase III efficacy study.",
      source_name        = "SAMPLE_INTENSIVE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 561L,
    n_observations = 2451L,
    n_studies      = 4L,
    age_range      = "10-17 years",
    weight_range   = "37 kg minimum (Study 1 inclusion criterion); Study 2 had no body-weight cutoff",
    sex_female_pct = 43,
    race_ethnicity = c(White = NA_real_, Black = NA_real_, Asian = 0, Other = 0),
    disease_state  = "Schizophrenia, bipolar I disorder, or other psychiatric disorders (autism, conduct disorder, oppositional defiant disorder) requiring chronic antipsychotic medication",
    dose_range     = "1-10 mg sublingual twice daily; Phase I cohorts at 1, 3, 5, 10 mg BID (Study 1) and 2.5, 5, 10 mg BID (Study 2); Phase III cohorts at 2.5, 5, 10 mg BID",
    regions        = "United States",
    notes          = "Pooled population PK analysis of two Phase I multiple ascending-dose studies (A7501022 / Study 1 in 40 patients; P06522 / Study 2 / NCT01206517 in 30 patients) and two Phase III fixed-dose efficacy studies (P06107 / NCT01244815 bipolar I disorder for 3 weeks; P05896 / NCT01190254 schizophrenia for 8 weeks). Inclusion limited to White and Black/African-American patients across all four studies (Discussion: 'demographics among patients who met inclusion criteria only included white and black/African-American racial groups'). Sex balance approximated from Study 1 (23M/17F) and Study 2 (17M/13F); the pooled Phase III demographic breakdown is not reported in the main text. Prespecified covariates age, BMI, body weight, sex, race, and dose were tested by stepwise covariate selection (forward p<0.01, backward p<0.001); none were retained in the final model. Asenapine is administered as rapidly-dissolving sublingual tablets; the model treats absorption as first-order from a depot compartment."
  )

  ini({
    # Structural parameters -- Dogterom 2018 Table 3 final-model population estimates.
    # Per Table 3 footnote b, V2/F, V3/F, and Ka were fixed in the final pooled-data
    # model to the values obtained when fitting the two Phase I PK studies alone;
    # CL/F and Q/F were estimated in the pooled fit. CL/F and Q/F are reported with
    # RSEs (3.09% and 14.3% respectively) and bootstrap CIs; the fixed parameters
    # have no RSE or bootstrap entry. (Note: Table 3 labels the Ka units as 'hours'
    # but Ka is a first-order rate constant; the implied units are 1/h.
    # Verification: ka = 2.98/h with kel = 296/2740 = 0.108/h gives Tmax = ln(ka/kel)/(ka-kel) =
    # 1.15 h, matching the paper's 'Tmax ~ 1 hour' and Table 1 medians 0.5-1.8 h.)
    lka     <- fixed(log(2.98));   label("Sublingual absorption rate constant Ka (1/h)")              # Table 3: Ka = 2.98 (fixed from Phase I-only fit)
    lcl     <- log(296);           label("Apparent oral clearance CL/F (L/h)")                        # Table 3: CL/F = 296 L/h (RSE 3.09%)
    lvc     <- fixed(log(2740));   label("Apparent central volume of distribution V2/F (L)")          # Table 3: V2/F = 2,740 L (fixed from Phase I-only fit)
    lq      <- log(120);           label("Apparent intercompartmental clearance Q/F (L/h)")           # Table 3: Q/F = 120 L/h (RSE 14.3%)
    lvp     <- fixed(log(2490));   label("Apparent peripheral volume of distribution V3/F (L)")       # Table 3: V3/F = 2,490 L (fixed from Phase I-only fit)
    # Bioavailability is not separately identifiable when CL/F and V/F are
    # reported as apparent values; F is anchored at 1 so the IIV(F) below
    # captures between-subject variability in the depot bioavailability.
    lfdepot <- fixed(log(1));      label("Sublingual bioavailability F (fixed anchor; CL and V are apparent)")  # structural anchor
    # Residual-variability scaling anchor (log scale; anchored at 1). The
    # paper's 'IIV (rV)' in Table 3 is encoded as etalrv below, scaling the
    # per-observation residual SD by exp(lrv + etalrv); the lrv anchor exists
    # so the eta pairs with a typical-value fixed effect per the nlmixr2lib
    # naming convention.
    lrv     <- fixed(log(1));      label("Residual-variability scaling anchor rV (fixed anchor)")  # structural anchor

    # Inter-individual variability (exponential model). Paper reports %CV;
    # the log-scale variance is omega^2 = log(1 + CV^2):
    #   IIV CL/F = 66.2%  -> omega^2 = log(1 + 0.662^2) = 0.36355
    #   IIV V2/F = 113%   -> omega^2 = log(1 + 1.13^2)  = 0.82309
    #   corr(CL/F, V2/F) = 0.921 -> cov = 0.921 * sqrt(0.36355 * 0.82309) = 0.50388
    #   IIV Ka   = 68.9%  -> omega^2 = log(1 + 0.689^2) = 0.38824 (fixed; no RSE
    #     reported in Table 3, consistent with Ka having been fixed structurally)
    #   IIV F    = 54.0%  -> omega^2 = log(1 + 0.540^2) = 0.25599
    #   IIV rV   = 19.2%  -> omega^2 = log(1 + 0.192^2) = 0.03620
    # Table 3 footnote also defines the correlation as cov / sqrt(var_cl * var_v2).
    etalcl + etalvc ~ c(0.36355, 0.50388, 0.82309)                                                    # Table 3: IIV CL/F = 66.2% (RSE 19.5%); IIV V2/F = 113% (RSE 21.3%); corr 0.921 (RSE 20.5%)
    etalka          ~ fixed(0.38824)                                                                  # Table 3: IIV Ka = 68.9% (no RSE; shrinkage 73.5%; treated as fixed alongside the structural Ka fix)
    etalfdepot      ~ 0.25599                                                                         # Table 3: IIV F = 54.0% (RSE 22.7%)
    etalrv          ~ 0.03620                                                                         # Table 3: IIV rV = 19.2% (RSE 29.8%) -- between-subject scaling of residual SD

    # Residual error. The paper used a 'log-additive error model'; following
    # the Macpherson 2015 convention this is encoded as nlmixr2's prop() in
    # linear space. Two magnitudes are estimated, switched per observation by
    # SAMPLE_INTENSIVE in model().
    propSdPhaseI   <- 0.278; label("Proportional residual error, intensive Phase I PK sampling (fraction)")        # Table 3: PK studies = 27.8% (RSE 5.29%)
    propSdPhaseIII <- 0.560; label("Proportional residual error, sparse Phase III efficacy sampling (fraction)")   # Table 3: Efficacy studies = 56.0% (RSE 5.05%)
  })

  model({
    # Individual PK parameters (exponential IIV).
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # Sublingual bioavailability on the depot compartment; structural F = 1
    # with exponential between-subject variability.
    f(depot) <- exp(lfdepot + etalfdepot)

    # Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system: sublingual depot -> central -> peripheral1, first-order
    # absorption and elimination, with reversible distribution to peripheral1.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # Plasma asenapine concentration. Dose is in mg, vc is in L, so central/vc
    # has units mg/L = ug/mL; multiply by 1000 to obtain ng/mL (the bioanalytical
    # reporting unit in both Phase I studies: LLOQ 0.025 ng/mL, ULOQ 20 ng/mL).
    Cc <- (central / vc) * 1000

    # Residual-error magnitude switched per observation by sampling intensity,
    # plus a between-subject log-normal scaling factor capturing IIV on residual
    # variability (Dogterom 2018 'IIV (rV)' row of Table 3):
    #   SAMPLE_INTENSIVE = 1 -> Phase I intensive sampling, 27.8% CV
    #   SAMPLE_INTENSIVE = 0 -> Phase III sparse sampling,   56.0% CV
    propSd <- (propSdPhaseI * SAMPLE_INTENSIVE + propSdPhaseIII * (1 - SAMPLE_INTENSIVE)) * exp(lrv + etalrv)
    Cc ~ prop(propSd)
  })
}
