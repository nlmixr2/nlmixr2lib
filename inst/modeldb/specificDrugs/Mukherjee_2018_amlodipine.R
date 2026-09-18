Mukherjee_2018_amlodipine <- function() {
  description <- paste(
    "Two-compartment oral pharmacokinetic reduction of the Simcyp",
    "minimal-PBPK-with-single-adjusting-compartment (SAC) model for",
    "amlodipine, linked to the paper's own systolic-blood-pressure (SBP)",
    "pharmacodynamic model (Mukherjee 2018). The PBPK model was built in",
    "Simcyp V15R1 and its whole-body mass-balance equations are not",
    "published, but the amlodipine compound layer is reported in full and",
    "is sufficient to rebuild the disposition as an ordinary compartmental",
    "model: first-order absorption into a depot with a lag time,",
    "distribution between a systemic compartment and the SAC (encoded as",
    "the canonical peripheral1 with the reported Q_SAC as q), and",
    "first-order elimination. Volumes are the reported per-kilogram values",
    "(systemic = Vss - V_SAC, peripheral1 = V_SAC) so they scale linearly",
    "with body weight. The reduction reproduces the paper's own predicted",
    "single-dose oral Cmax, Tmax, AUC-infinity and terminal half-life and",
    "the intravenous AUC-infinity to within 2.6 percent at the",
    "back-solved Simcyp population-representative weight of 82.2 kg.",
    "The PD layer is transcribed from the Lua script the authors deposited",
    "as supplementary Figure S2, which is the code that generated the",
    "published figures: SBP is a cosine baseline plus a linear",
    "concentration effect that is delayed by an effect-compartment",
    "build-up factor 1 - exp(-keo * time-after-first-dose). Note that the",
    "main-text Equations 1 and 3 print this factor as exp(-keo * t), which",
    "would abolish the drug effect within days and contradicts the",
    "sustained day-43 effect the paper shows in Figure 4; the deposited",
    "code is taken as authoritative.",
    "Ritonavir co-administration enters as the empirical binary covariate",
    "CONMED_RTV acting on relative bioavailability and on clearance,",
    "calibrated to the paper's own predicted Cmax and AUC ratios for",
    "ritonavir 100 mg once daily (Table 4, Menon row). The mechanistic",
    "time-dependent CYP3A4 inhibition and induction that produced those",
    "ratios lives in the separate Simcyp ritonavir compound model of",
    "Shebley 2017 and is NOT reproducible here, so the time course of",
    "onset and washout of the interaction (the paper's Figure 3) is not",
    "captured; see the vignette for the full list of deviations.",
    "This is a typical-value simulation model: the source reports no",
    "inter-individual variance components and no pharmacokinetic residual",
    "error, so there are no etas and propSd is fixed at zero.",
    sep = " "
  )
  reference <- paste(
    "Mukherjee D, Zha J, Menon RM, Shebley M. (2018).",
    "Guiding dose adjustment of amlodipine after co-administration with",
    "ritonavir containing regimens using a physiologically-based",
    "pharmacokinetic/pharmacodynamic model.",
    "J Pharmacokinet Pharmacodyn 45(3):443-456.",
    "doi:10.1007/s10928-018-9574-0.",
    sep = " "
  )
  vignette <- "Mukherjee_2018_amlodipine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed. Mukherjee 2018 Table 1 reports both distribution volumes",
        "per kilogram (Vd 21.4 L/kg, V_SAC 11 L/kg), so the systemic and",
        "peripheral volumes scale linearly with body weight; the exponent is",
        "therefore structurally 1 rather than an estimated allometric power.",
        "The reference weight 82.2 kg is not printed anywhere in the paper -",
        "it is the Simcyp 'population representative' healthy-volunteer weight,",
        "back-solved from the paper's own predicted terminal half-life of",
        "39.9 h (Table 4) and audited against the four other printed",
        "predictions; see the vignette Errata. Clearance and inter-compartmental",
        "clearance are reported in absolute L/h and are NOT weight-scaled.",
        sep = " "
      ),
      source_name = "not reported"
    ),
    CONMED_RTV = list(
      description = "Concomitant ritonavir co-administration indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no ritonavir)",
      notes = paste(
        "1 = ritonavir 100 mg once daily at steady state, the perpetrator",
        "regimen used in every simulation scenario of Mukherjee 2018. The two",
        "coefficients are NOT mechanistic: they are back-solved so that a",
        "single 5 mg oral amlodipine dose reproduces the paper's own predicted",
        "Cmax ratio of 1.42 and AUC-infinity ratio of 2.28 against the",
        "ritonavir-free reference (Table 4, Menon et al. row). The underlying",
        "mechanism - reversible plus mechanism-based CYP3A4 inhibition together",
        "with CYP3A4 induction - is carried by the separate Simcyp ritonavir",
        "compound model of Shebley 2017, which is not on disk and is a platform",
        "database artefact, so no time-varying interaction term can be encoded.",
        "The indicator is therefore a steady-state switch: it cannot reproduce",
        "the onset and 5-day washout of the interaction shown in Figure 3.",
        "The paper's other reported interaction arm (indinavir 800 mg twice",
        "daily plus ritonavir 100 mg twice daily; predicted Cmax ratio 1.74 and",
        "AUC24 ratio 1.89, Glesby et al. row of Table 4) is a different",
        "perpetrator combination and is deliberately NOT folded into this",
        "column.",
        sep = " "
      ),
      source_name = "not reported"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "amlodipine",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "amlodipine",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "amlodipine",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 12,
    n_studies = 11,
    age_range = "23-64 years across the contributing studies (Table S2)",
    disease_state = paste(
      "Healthy adult volunteers for the pharmacokinetic layer; adults with",
      "essential hypertension for the systolic-blood-pressure layer.",
      sep = " "
    ),
    dose_range = "2.5-20 mg oral once daily; 10 mg intravenous over 10 min",
    notes = paste(
      "The PK model was optimised against the 12-subject intravenous and oral",
      "single-dose study of Faulkner 1986 and the 18-subject",
      "indinavir/ritonavir interaction study of Glesby 2005, and verified",
      "against four further studies (Table 2). All simulations in the paper",
      "use the Simcyp 'population representative' virtual healthy volunteer,",
      "i.e. a single typical subject rather than a sampled population. The PD",
      "parameters were fitted to the mean systolic blood pressure of the 12",
      "hypertensive patients of Donnelly 1993 (aged 25-64 years, amlodipine",
      "5 mg once daily for 6 weeks) - only mean observations were published,",
      "so no inter-individual variability could be estimated even though the",
      "paper notes it was high for m and keo.",
      "Reported fractional contributions to amlodipine clearance, retained",
      "here as provenance because the pathway split is kinetically inert in a",
      "compartmental reduction: CYP3A4 intrinsic clearance 170 L/h, CYP3A5",
      "intrinsic clearance 43.5 L/h, biliary clearance 12 L/h, non-specific",
      "systemic clearance 16 L/h and renal clearance 1.8 L/h (6 percent of",
      "total), all Table 1.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Nothing below is estimated in this reduction. Every value is either
    # a Mukherjee 2018 Table 1 / Table 3 / Table 4 entry, a verbatim
    # constant from the deposited Lua script (supplementary Figure S2),
    # or an arithmetic consequence of one of those.
    #
    # The single quantity the paper does not print is the body weight of
    # the Simcyp 'population representative'. It is needed because
    # Table 1 gives both volumes per kilogram. It was obtained by
    # requiring the reduction to return the paper's own predicted
    # terminal half-life of 39.9 h (Table 4), which gives 82.2 kg, and
    # the four remaining printed predictions were then held out as an
    # audit: oral Cmax -1.6%, oral Tmax -0.2%, oral AUC-infinity -2.3%,
    # intravenous AUC-infinity -2.6%.
    # ------------------------------------------------------------------

    lka <- fixed(log(0.75))
    label("First-order absorption rate constant ka (1/h)")
    # Table 1, 'Optimized parameters': ka = 0.75 1/h, optimised against the
    # oral plasma profile of Faulkner 1986 and close to the 0.8 1/h initial
    # estimate of Flynn et al.

    ltlag <- fixed(log(3.2))
    label("Absorption lag time (h)")
    # Table 1, 'Optimized parameters': absorption lag time = 3.2 h.

    lfdepot <- fixed(log(0.666))
    label("Oral bioavailability (unitless)")
    # Table 4, Faulkner et al. (oral) row: model-predicted F = 66.6%
    # (observed 64%). Table 1 records that the fraction absorbed fa = 1, so
    # this is entirely first-pass extraction.

    lcl <- fixed(log(33.9))
    label("Systemic plasma clearance CL (L/h)")
    # Table 1, 'Initial estimates': 33.9 L/h intravenous clearance from
    # Faulkner 1986. Consistent with the final model's own prediction:
    # Table 4 gives a predicted intravenous AUC-infinity of 303 ng*h/mL
    # after 10 mg, i.e. 10 mg / 303 ng*h/mL = 33.0 L/h.

    lvc <- fixed(log(854.88))
    label("Systemic compartment volume vc (L) at the 82.2 kg reference weight")
    # Table 1: Vd = 21.4 L/kg and V_SAC = 11 L/kg. In the Simcyp minimal-PBPK
    # layout the systemic compartment holds the balance of the steady-state
    # volume, so vc = (21.4 - 11) L/kg * 82.2 kg = 854.88 L. The liver and
    # portal-vein volumes of the platform layout are not reported and are not
    # subtracted; together they are under 0.3% of this volume.

    lvp <- fixed(log(904.2))
    label("Single adjusting compartment volume vp (L) at the 82.2 kg reference weight")
    # Table 1, 'Optimized parameters': V_SAC = 11 L/kg; 11 * 82.2 = 904.2 L.

    lq <- fixed(log(90))
    label("Inter-compartmental clearance to the single adjusting compartment q (L/h)")
    # Table 1, 'Optimized parameters': Q_SAC = 90 L/h. Reported in absolute
    # L/h, so unlike the volumes it is not weight-scaled.

    e_wt_vc_vp <- fixed(1)
    label("Exponent of body weight on vc and vp (unitless)")
    # Structural, not estimated: Table 1 reports both volumes in L/kg, which
    # is exactly linear scaling.

    e_conmed_rtv_fdepot <- fixed(1.3783)
    label("Ritonavir multiplier on oral bioavailability (unitless)")
    # Back-solved together with e_conmed_rtv_cl so that a 5 mg single oral
    # dose reproduces both ritonavir ratios predicted in Table 4, Menon
    # et al. row: Cmax ratio 1.42 and AUC-infinity ratio 2.28. Empirical,
    # not mechanistic - see covariateData[['CONMED_RTV']].

    e_conmed_rtv_cl <- fixed(0.6045)
    label("Ritonavir multiplier on systemic clearance (unitless)")
    # Partner of e_conmed_rtv_fdepot; 1.3783 / 0.6045 = 2.28, the predicted
    # AUC-infinity ratio of Table 4, Menon et al. row.

    # ---- Systolic blood pressure model -------------------------------
    # Main-text Equation 3 with the parameter values of Table 3, encoded
    # as the authors implemented it in the Lua script of supplementary
    # Figure S2 (Po 148.84, A 8.245, om 0.463, m -1285.93, keo 0.049).

    le0 <- fixed(log(148.84))
    label("Start-of-day baseline systolic blood pressure P0 (mmHg)")
    # Figure S2 Lua script, 'Po'; Table 3 prints the same value rounded to
    # 148.8 mmHg.

    lamp <- fixed(log(8.245))
    label("Amplitude of the baseline systolic blood pressure oscillation (mmHg)")
    # Figure S2 Lua script, 'A'; Table 3 prints 8.25 mmHg.

    lfcirc <- fixed(log(1.76))
    label("Frequency of the baseline systolic blood pressure oscillation (cycles/day)")
    # Table 3: circadian frequency f = 1.76 1/day. Cross-checks against the
    # Lua script's angular frequency om = 0.463 1/h, since
    # 2 * pi * 1.76 / 24 = 0.4608 1/h.

    lke0 <- fixed(log(0.049))
    label("Effect-delay rate constant keo (1/h)")
    # Table 3 and Figure S2 Lua script: keo = 0.049 1/h.

    slope_drug <- fixed(-3.145)
    label("Linear effect of amlodipine on systolic blood pressure (mmHg per ng/mL)")
    # Table 3: m = -3.145 mmHg*mL/ng. Cross-checks exactly against the Lua
    # script's molar form m = -1285.93 mmHg/uM, since
    # 3.145 * 408.88 g/mol = 1285.93.

    propSd <- fixed(0)
    label("Proportional residual error on plasma concentration (fraction)")
    # Not reported: the source is a deterministic platform simulation of a
    # single population representative, with no residual-error model.

    addSd_SBP <- fixed(5.04)
    label("Additive residual error on systolic blood pressure (mmHg)")
    # Table 3, 'Drug effect model' block, residual-error column: 5.04.
    # Footnote c defines it as sqrt(sum((y - y')^2) / N), i.e. a root mean
    # squared error in mmHg. The corresponding baseline-only figure is 2.89.
  })

  model({
    # Reference body weight for the per-kilogram volumes of Table 1 (kg).
    # Not printed in the paper; back-solved from the predicted terminal
    # half-life - see the note at the top of ini().
    wtref <- 82.2

    # Molecular weight of amlodipine, used only to state the unit chain
    # (g/mol, Table 1); the concentration-effect slope is already carried
    # in mmHg per ng/mL so no conversion is applied here.

    ka <- exp(lka)
    tlag <- exp(ltlag)
    q <- exp(lq)
    vc <- exp(lvc) * (WT / wtref)^e_wt_vc_vp
    vp <- exp(lvp) * (WT / wtref)^e_wt_vc_vp
    cl <- exp(lcl) * e_conmed_rtv_cl^CONMED_RTV
    fdepot <- exp(lfdepot) * e_conmed_rtv_fdepot^CONMED_RTV

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central +
      k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    alag(depot) <- tlag
    f(depot) <- fdepot

    # central is in mg and vc in L, so the factor 1000 converts mg/L to
    # ng/mL, the units the paper reports plasma amlodipine in.
    Cc <- 1000 * central / vc

    # Systolic blood pressure, Equation 3 as implemented in Figure S2.
    # tday is clock time within the day, the Lua script's t - 24*floor(t/24);
    # the oscillation restarts at each midnight boundary. The drug term uses
    # time after the first dose, so the effect builds up over the first days
    # of therapy and is essentially complete by day 43.
    e0 <- exp(le0)
    amp <- exp(lamp)
    fcirc <- exp(lfcirc)
    ke0 <- exp(lke0)
    tday <- t - 24 * floor(t / 24)
    SBP <- e0 + amp * cos(2 * pi * fcirc / 24 * tday) +
      slope_drug * Cc * (1 - exp(-ke0 * max(0, tafd())))

    Cc ~ prop(propSd)
    SBP ~ add(addSd_SBP)
  })
}
