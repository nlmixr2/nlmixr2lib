Morse_2012_ghb_rbc_invitro <- function() {
  description <- paste(
    "In vitro (rat Sprague-Dawley erythrocytes).",
    "Saturable-plus-linear membrane transport model for uptake of the drug",
    "of abuse gamma-hydroxybutyrate (GHB) into freshly isolated rat red",
    "blood cells, characterised at two buffer pH values. Total unidirectional",
    "uptake is the sum of two parallel arms (source eq. 1): a saturable",
    "Michaelis-Menten arm carried by monocarboxylate transporter 1 (MCT1),",
    "vmax_rbc * C / (km_rbc + C), and a linear arm carried by the anion",
    "exchanger band 3 (AE1) together with passive diffusion, kinf_rbc * C.",
    "The two arms were resolved experimentally by selective inhibition: 1 mM",
    "p-chloromercuribenzene sulfonate (pCMBS) abolishes the saturable arm and",
    "10 uM 4,4'-diisothiocyanostilbene-2,2'-disulfonate (DIDS) abolishes the",
    "linear arm. Affinity for MCT1 is strongly pH dependent - the Michaelis",
    "constant is 17.0 mM at pH 7.4 but only 2.2 mM at pH 6.5 - so the whole",
    "parameter set is switched on PH_MEDIUM rather than interpolated; only",
    "those two pH levels were studied. A third arm of the paper titrated the",
    "competing MCT substrate L-lactate against uptake of 10 mM GHB at pH 7.4",
    "(source eq. 2) and is carried here as a sigmoidal Imax factor on total",
    "uptake, driven by the LACT covariate.",
    "The model is a UNIDIRECTIONAL INITIAL-RATE law: uptake was measured over",
    "5 s (pH 6.5) and 15 s (pH 7.4) windows chosen to lie in the linear range,",
    "so it carries no efflux, no trans-stimulation and no approach to",
    "equilibrium, and the authors state explicitly that equilibrium-exchange",
    "km and vmax would be higher than these unidirectional values. States are",
    "normalised per mg of red-cell protein, the units the assay reports, so",
    "rbc_ghb holds nmol/mg protein rather than the concentration units usual",
    "for the rbc_<analyte> family.",
    "Companion in vivo result, not part of this model: blood/plasma",
    "partitioning of GHB in rats was LINEAR across 400-1500 mg/kg with a B/P",
    "ratio of 0.75, and was unchanged by co-administered L-lactate (0.76),",
    "even though the same L-lactate dose significantly increased GHB renal",
    "and total clearance. The in vitro parameters here explain that",
    "dissociation: at blood pH the L-lactate IC50 is 19.1 mM, far above the",
    "2-5 mM plasma L-lactate actually reached in vivo.",
    sep = " "
  )
  reference <- paste(
    "Morse BL, Felmlee MA, Morris ME.",
    "gamma-Hydroxybutyrate blood/plasma partitioning: effect of physiologic",
    "pH on transport by monocarboxylate transporters.",
    "Drug Metab Dispos. 2012;40(1):64-69.",
    "doi:10.1124/dmd.111.041285. PMID: 21976621. PMCID: PMC3250051.",
    "Structural equations: eq. 1 (saturable plus linear uptake) and eq. 2",
    "(sigmoidal Imax inhibition by L-lactate), Materials and Methods,",
    "'Data and Statistical Analysis', p. 66.",
    "Parameter values for eq. 1 at both pH values and the L-lactate IC50:",
    "Results, 'In Vitro GHB Uptake', p. 67.",
    "The Imax and Hill coefficient of eq. 2 are reported nowhere in the",
    "paper and were recovered by digitising the fitted curve of Figure 5;",
    "see the vignette 'Assumptions and deviations'.",
    "No supplementary material accompanies this article and no erratum was",
    "located.",
    sep = " "
  )
  vignette <- "Morse_2012_ghb_rbc_invitro"

  units <- list(time = "min", dosing = "mmol/L (bath concentration)", concentration = "mmol/L")

  # The bath state is dosed directly, following the in vitro convention of
  # HernandezLozano_2025_apramycin_invitro.R: there is no depot and no
  # central compartment, so the "dose" is the extracellular GHB
  # concentration placed into the buffer at time 0.
  dosing <- c("ghb")

  # `ghb` holds the static extracellular (bath) GHB concentration and is not
  # a canonical compartment role; `rbc_ghb` is a regex-validated member of
  # the rbc_<analyte> family (see inst/references/compartment-names.md).
  paper_specific_compartments <- c("ghb")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Checked against Morse 2012 Materials and Methods,
  # 'In Vitro Uptake' (p. 66) and eq. 1.
  compartmentData <- list(
    ghb = list(
      analyte = "gamma-hydroxybutyrate", units = "mmol/L", specimen = "administration site",
      verified = TRUE
    ),
    rbc_ghb = list(
      analyte = "gamma-hydroxybutyrate", units = "nmol/mg red-cell protein",
      specimen = "blood cell", verified = TRUE
    )
  )

  covariateData <- list(
    PH_MEDIUM = list(
      description = paste(
        "pH of the HEPES-buffered uptake medium in which the erythrocyte",
        "suspension was incubated. Selects the whole eq. 1 parameter set",
        "(vmax_rbc, km_rbc, kinf_rbc), because the paper fitted pH 7.4 and",
        "pH 6.5 independently rather than estimating a pH effect."
      ),
      units = "pH units",
      type = "continuous",
      reference_category = NULL,
      source_name = "pH",
      notes = paste(
        "Uptake buffer was 137 mM NaCl, 5.4 mM KCl, 2.8 mM CaCl2, 1.2 mM",
        "MgCl2 and 10 mM HEPES titrated to the stated pH (Materials and",
        "Methods, 'In Vitro Uptake').",
        "Only two levels were studied, 7.4 (blood) and 6.5 (chosen to",
        "approximate the more acidic renal tubule, the other MCT-mediated",
        "site in GHB disposition). The model dispatches on the midpoint",
        "threshold 6.95 rather than interpolating, because no intermediate",
        "pH was characterised; a value below 6.95 is treated as the pH 6.5",
        "condition and a value at or above it as the pH 7.4 condition.",
        "The direction is load-bearing for the paper's argument: lower pH",
        "gives HIGHER MCT1 affinity (km_rbc 2.2 mM at pH 6.5 vs 17.0 mM at",
        "pH 7.4), which is why GHB transport saturates and is inhibitable by",
        "L-lactate in the kidney but not at the red-cell membrane."
      )
    ),
    LACT = list(
      description = paste(
        "L-lactate concentration in the uptake buffer, acting as a competing",
        "MCT substrate. Drives the sigmoidal Imax inhibition of total GHB",
        "uptake in source eq. 2. Set to 0 for the uninhibited condition."
      ),
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      source_name = "L-lactate",
      notes = paste(
        "In vitro range studied: 0 to 150 mM (Figure 5), titrated against",
        "uptake of 10 mM GHB at pH 7.4 only.",
        "The inhibition parameters were therefore characterised at a single",
        "GHB concentration and a single pH; the model applies the factor",
        "multiplicatively to total uptake at any GHB concentration, which is",
        "the form of eq. 2 but is an extrapolation beyond the fitted",
        "condition. It is NOT gated on PH_MEDIUM, because no L-lactate",
        "titration was run at pH 6.5; a user simulating the pH 6.5 condition",
        "with non-zero LACT is extrapolating and should say so.",
        "This column is the same analyte as the clinical serum-lactate",
        "covariate and deliberately reuses that canonical: the paper's own",
        "argument turns on comparing the in vitro IC50 of 19.1 mM against",
        "the 2-5 mM PLASMA L-lactate its rats actually reached in vivo",
        "(Figure 3), so one column must span both matrices. Note the",
        "in vitro concentrations are far above any physiologic value."
      )
    )
  )

  population <- list(
    species = "in vitro (rat Sprague-Dawley erythrocytes)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    disease_state = "not applicable (in vitro transport assay)",
    dose_range = "0-150 mM GHB in the uptake buffer; 0-150 mM L-lactate in the inhibition arm",
    notes = paste(
      "Erythrocytes were freshly isolated from male Sprague-Dawley rats",
      "(Harlan), washed three times to remove extracellular lactate, and",
      "resuspended to a 40 percent suspension in uptake buffer. Uptake was",
      "stopped by centrifugation through silicone oil. Unless stated",
      "otherwise, experiments were run at ROOM TEMPERATURE (not 37 C) and",
      "in triplicate, repeated three times; the L-lactate titration is from",
      "two experiments in triplicate. The authors note that km and IC50 for",
      "MCT substrates rise with temperature, so the in vivo values at 37 C",
      "are expected to exceed these.",
      "Number of animals contributing erythrocytes is not reported, hence",
      "n_subjects = NA.",
      "The accompanying in vivo arm - which this model does not describe -",
      "used male Sprague-Dawley rats of 270-330 g given 400-1500 mg/kg GHB",
      "intravenously, with or without L-lactate, n = 74 rats in total across",
      "serial and sacrificial sampling (Materials and Methods, 'In Vivo",
      "Experimental Protocol'; Table 1)."
    )
  )

  ini({
    # ---- Source eq. 1, pH 7.4 (Results, 'In Vitro GHB Uptake', p. 67) ----
    # "Total transport and transport in the presence of pCMBS at pH 7.4 were
    #  simultaneously fitted with nonlinear regression, resulting in values of
    #  17.0 mM, 20.9 nmol/mg per minute, and 0.24 ul/mg per minute for Km,
    #  Vmax, and P at pH 7.4, respectively."
    # All six eq. 1 values are point estimates from WinNonlin nonlinear
    # regression; the paper reports no standard errors or confidence
    # intervals for any of them, so none is carried here.
    lvmax_rbc_ph74 <- log(20.9); label("Log maximum rate of MCT1-mediated GHB influx into erythrocytes at pH 7.4 (nmol/mg protein/min)")  # Results p.67: Vmax = 20.9 nmol/mg per minute at pH 7.4
    lkm_rbc_ph74 <- log(17.0); label("Log Michaelis constant for MCT1-mediated GHB influx at pH 7.4 (mmol/L)")                            # Results p.67: Km = 17.0 mM at pH 7.4
    lkinf_rbc_ph74 <- log(0.24); label("Log linear (band 3 plus passive) GHB influx coefficient at pH 7.4 (uL/mg protein/min)")           # Results p.67: P = 0.24 ul/mg per minute at pH 7.4

    # ---- Source eq. 1, pH 6.5 (Results, 'In Vitro GHB Uptake', p. 67) ----
    # "Fitting the total and pCMBS-inhibited transport at pH 6.5 resulted in
    #  a much lower Km value of 2.2 mM, with Vmax and P values of 5.3
    #  nmol/mg per minute and 0.11 ul/mg per minute, respectively."
    lvmax_rbc_ph65 <- log(5.3); label("Log maximum rate of MCT1-mediated GHB influx into erythrocytes at pH 6.5 (nmol/mg protein/min)")   # Results p.67: Vmax = 5.3 nmol/mg per minute at pH 6.5
    lkm_rbc_ph65 <- log(2.2); label("Log Michaelis constant for MCT1-mediated GHB influx at pH 6.5 (mmol/L)")                             # Results p.67: Km = 2.2 mM at pH 6.5
    lkinf_rbc_ph65 <- log(0.11); label("Log linear (band 3 plus passive) GHB influx coefficient at pH 6.5 (uL/mg protein/min)")           # Results p.67: P = 0.11 ul/mg per minute at pH 6.5

    # ---- Source eq. 2: L-lactate inhibition at pH 7.4 -------------------
    # IC50 is the only one of the three eq. 2 parameters the paper prints.
    lic50_lact <- log(19.1); label("Log L-lactate concentration giving half-maximal inhibition of GHB uptake at pH 7.4 (mmol/L)")         # Results p.67: "Fitting GHB transport in the presence of L-lactate resulted in an IC50 value for L-lactate of 19.1 mM (Fig. 5)."

    # FIGURE-DERIVED, NOT PRINTED. Imax and the Hill coefficient (the paper's
    # gamma) appear in eq. 2 but are reported in neither the text, Table 1,
    # Table 2, nor the Figure 5 caption. Both were recovered by digitising
    # the fitted solid line of Figure 5 (page 68) at 300 dpi, calibrating on
    # the axis frame and label positions, and least-squares fitting eq. 2 to
    # 650 traced curve points with IC50 held at the published 19.1 mM. The
    # fit is essentially exact (RMSE 0.004 nmol/mg/min, below the 0.014-unit
    # pixel resolution) and the digitisation is corroborated by re-fitting
    # with IC50 FREE, which returns 19.7 mM against the published 19.1 mM
    # -- a 3 percent recovery of a number that was not used to fit.
    # Fixing the Hill coefficient to 1 instead degrades the fit 20-fold
    # (RMSE 0.082), so the sigmoidicity is a real feature of the published
    # curve and not an artefact of the tracing. See the vignette
    # 'Assumptions and deviations'.
    limax_lact <- log(0.613); label("Log maximum fractional inhibition of total GHB uptake by L-lactate at pH 7.4 (unitless) -- figure-derived")  # Digitised from Figure 5 fitted curve; asymptotes 7.34 (control) and 2.84 (plateau)
    lhill_lact <- log(1.53); label("Log Hill coefficient of the L-lactate inhibition curve (unitless) -- figure-derived")                         # Digitised from Figure 5 fitted curve; the paper's eq. 2 symbol is gamma

    # ---- Residual error ---------------------------------------------------
    # The in vitro data were fitted by nonlinear regression in WinNonlin
    # (Materials and Methods, 'Data and Statistical Analysis'); no residual
    # error model and no sigma are reported anywhere in the paper. Held at 0
    # rather than invented. Figures 4 and 5 plot mean +/- S.E.M. of
    # triplicate determinations. See the vignette Errata.
    addSd <- fixed(0); label("Additive residual SD on accumulated red-cell GHB (nmol/mg protein); not reported by the source")  # Materials and Methods p.66: nonlinear regression in WinNonlin; no residual error model reported
  })

  model({
    # ---- 1. Select the pH-specific eq. 1 parameter set -------------------
    # The paper fitted pH 7.4 and pH 6.5 as independent datasets and shares
    # no parameter between them, so this is a switch and not a covariate
    # slope. Only the two studied levels are characterised; 6.95 is the
    # midpoint dispatch threshold.
    if (PH_MEDIUM < 6.95) {
      vmax_rbc <- exp(lvmax_rbc_ph65)
      km_rbc <- exp(lkm_rbc_ph65)
      kinf_rbc <- exp(lkinf_rbc_ph65)
    } else {
      vmax_rbc <- exp(lvmax_rbc_ph74)
      km_rbc <- exp(lkm_rbc_ph74)
      kinf_rbc <- exp(lkinf_rbc_ph74)
    }

    imax_lact <- exp(limax_lact)
    ic50_lact <- exp(lic50_lact)
    hill_lact <- exp(lhill_lact)

    # ---- 2. Source eq. 2: sigmoidal inhibition by L-lactate --------------
    # Uptake = Control * (1 - Imax * C^gamma / (IC50^gamma + C^gamma)).
    # Written as a dimensionless multiplier on total uptake, so that the
    # role played by "Control" in the published equation is filled by the
    # eq. 1 prediction. LACT = 0 gives inhib = 1 and recovers eq. 1 exactly.
    inhib <- 1 - imax_lact * LACT^hill_lact / (ic50_lact^hill_lact + LACT^hill_lact)

    # ---- 3. Source eq. 1: the two parallel transport arms ----------------
    # `ghb` holds the extracellular bath concentration in mmol/L directly,
    # so no volume term is needed. Bath concentration is treated as constant
    # over the 5-15 s initial-rate window, which is what the assay measures.
    uptake_mct <- vmax_rbc * ghb / (km_rbc + ghb)
    uptake_lin <- kinf_rbc * ghb
    uptake <- (uptake_mct + uptake_lin) * inhib

    # Fraction of total uptake that is MCT-mediated, i.e. the fraction
    # abolished by 1 mM pCMBS. This is the quantity tabulated in Table 2 and
    # is exposed as an output so it can be checked directly against it.
    fmct <- uptake_mct / (uptake_mct + uptake_lin)

    # Uptake in the presence of pCMBS (saturable arm abolished) and of DIDS
    # (linear arm abolished), the two inhibitor conditions plotted in
    # Figure 4.
    uptake_pcmbs <- uptake_lin * inhib
    uptake_dids <- uptake_mct * inhib

    # ---- 4. ODE system ----------------------------------------------------
    d/dt(ghb) <- 0
    d/dt(rbc_ghb) <- uptake

    # ---- 5. Observation and error ----------------------------------------
    # The assay datum is the radiolabel accumulated inside the cells over a
    # fixed short window (5 s at pH 6.5, 15 s at pH 7.4), which the paper
    # divides by that window to report a rate. The observation is therefore
    # the accumulated state, and `uptake` above is its time derivative.
    # Because the bath is constant, uptake is constant in time and the two
    # are related exactly by rbc_ghb(t) = uptake * t, so solving to t = 1 min
    # reads the published uptake rate off the state directly.
    rbc_ghb ~ add(addSd)
  })
}
