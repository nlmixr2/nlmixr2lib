Kroemer_2024_ceftazidime_avibactam_fosfomycin_tkc <- function() {
  description <- "In vitro (Escherichia coli, clinical MDR isolate; CTX-M-15 + TEM-4 ESBL and OXA-244 carbapenemase). Semi-mechanistic PK/PD model of the static time-kill experiments for ceftazidime, avibactam and fosfomycin: a susceptible and a joint resistant bacterial subpopulation under a shared capacity limit, with sigmoidal Emax or power-model mono-drug kill, Bliss independence as the additivity criterion, and general pharmacodynamic interaction (GPDI) terms describing the >99% reduction of the ceftazidime EC50 by avibactam and of the fosfomycin EC50 by ceftazidime. Drug concentrations are static covariates (no dosing events). Sibling model: Kroemer_2024_ceftazidime_avibactam_fosfomycin_hfim."
  reference <- "Kroemer N, Amann LF, Farooq A, Pfaffendorf C, Martens M, Decousser JW, Gregoire N, Nordmann P, Wicha SG. Pharmacokinetic/pharmacodynamic analysis of ceftazidime/avibactam and fosfomycin combinations in an in vitro hollow fiber infection model against multidrug-resistant Escherichia coli. Microbiol Spectr. 2024 Jan;12(1):e03318-23. doi:10.1128/spectrum.03318-23. Structural equations: supplemental material Text S3 Eqs 1-8. Parameter estimates: supplemental material Table S3 (with 95% confidence intervals from sampling importance resampling)."
  vignette <- "Kroemer_2024_ceftazidime_avibactam_fosfomycin"
  units <- list(time = "h", dosing = "mg/L", concentration = "mg/L")

  # Bacterial subpopulations of the static time-kill model (Text S3 Eqs 1-2).
  # Neither maps onto a canonical PK compartment name.
  paper_specific_compartments <- c("bact_susceptible", "bact_resistant")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    bact_susceptible = list(analyte = "Escherichia coli (susceptible bacteria)", units = NA_character_, specimen = "bile", verified = FALSE),
    bact_resistant   = list(analyte = "Escherichia coli (resistant bacteria)", units = NA_character_, specimen = "bile", verified = FALSE)
  )

  covariateData <- list(
    CONC_CAZ_MGL = list(
      description = "Static (time-invariant) ceftazidime concentration in the time-kill tube. Drives the sigmoidal Emax kill of both bacterial subpopulations and, as perpetrator, the GPDI shift of the fosfomycin EC50 on the resistant subpopulation.",
      units = "mg/L",
      type = "continuous",
      reference_category = "n/a -- 0 mg/L is the drug-free growth control; the static time-kill experiments spanned 0.002-128 mg/L (Results, Static time-kill experiments; Fig. 1)",
      notes = "Held constant for the 30-h experiment; the static model has no PK component. Ceftazidime MIC of the isolate was 16 mg/L alone and 0.125 mg/L in the presence of 4 mg/L avibactam (Results, Bacterial isolate and susceptibility).",
      source_name = "CAZ"
    ),
    CONC_AVI_MGL = list(
      description = "Static (time-invariant) avibactam concentration in the time-kill tube. Drives its own (weak) kill of both bacterial subpopulations and, as perpetrator, the GPDI shift of the ceftazidime EC50 on both subpopulations.",
      units = "mg/L",
      type = "continuous",
      reference_category = "n/a -- 0 mg/L is the avibactam-free arm; concentrations up to ~64 mg/L were applied (Text S5: the avibactam effect at the highest applied concentration was 0.294 1/h, which the power model slope_avi_r = 0.0787 and hill_avi_r = 0.317 reach at ~64 mg/L)",
      notes = "Held constant for the 30-h experiment. The bioanalysis found no degradation of ceftazidime at the avibactam concentrations studied, so a permanent (rather than concentration-dependent) inhibition of the beta-lactamases was assumed and no PK drug interaction is modelled (Discussion; Text S5).",
      source_name = "AVI"
    ),
    CONC_FOF_MGL = list(
      description = "Static (time-invariant) fosfomycin concentration in the time-kill tube. Drives a power-model kill of the susceptible subpopulation and a sigmoidal Emax kill of the resistant subpopulation.",
      units = "mg/L",
      type = "continuous",
      reference_category = "n/a -- 0 mg/L is the drug-free growth control; the static time-kill experiments spanned 2-16 mg/L (Results, Static time-kill experiments; Fig. 1)",
      notes = "Held constant for the 30-h experiment. Fosfomycin MIC of the isolate was 16/0.5 mg/L (microdilution/agar dilution). All fosfomycin-containing media were supplemented with 25 mg/L glucose-6-phosphate per EUCAST recommendations.",
      source_name = "FOF"
    )
  )

  population <- list(
    species = "in vitro (Escherichia coli, single clinical multidrug-resistant isolate)",
    n_subjects = 1L,
    n_studies = 1L,
    disease_state = "Clinical E. coli isolated from a rectal swab in a multidrug-resistance screening at Henri-Mondor Hospital, Paris. Whole-genome sequencing identified CTX-M-15 and TEM-4 (extended-spectrum beta-lactamase) and OXA-244 (OXA-48-like carbapenemase) genes. MICs: ceftazidime 16 mg/L, ceftazidime/avibactam 0.125 mg/L, fosfomycin 16/0.5 mg/L (microdilution/agar dilution); susceptible to ceftazidime/avibactam and fosfomycin, resistant to ceftazidime alone (EUCAST). The same isolate is E. coli YAL_AMA of Kroemer et al. 2023.",
    model_system = "Static time-kill experiments over 30 h at 37 degrees C in ambient air, 10 mL cation-adjusted Mueller-Hinton broth 2, inoculated at 10^6 CFU/mL. Drugs were added after 2 h of preincubation and bacterial counts quantified at 0, 2, 4, 8, 24 and 30 h after drug addition by serial dilution, plating and manual colony counting. Experiments were performed in duplicate.",
    initial_inoculum = "10^6 CFU/mL targeted; 10^6.86 CFU/mL estimated for the susceptible subpopulation, 10^3.14 CFU/mL for the resistant subpopulation",
    dose_range = "Static concentrations: ceftazidime 0.002-128 mg/L, avibactam up to ~64 mg/L, fosfomycin 2-16 mg/L, alone and in combination",
    notes = paste(
      "Estimated in NONMEM 7.5.0 with second-order conditional estimation with interaction (LAPLACIAN-I) and the CVODES solver (ADVAN14).",
      "Built sequentially: mono-drug effects first, then the interaction parameters with the mono-drug parameters fixed, then all parameters estimated simultaneously.",
      "Parameter uncertainty from sampling importance resampling in PsN 5.0, using the covariance-step relative standard errors as the proposal distribution.",
      "The lower limit of quantification of the plating method was about 10^1-10^2 CFU/mL; non-quantifiable counts were empirically set to 1 CFU/mL and included in the model-based data evaluation.",
      "Inter-experimental variability is carried as interindividual variability (exponential) on the inoculum and the growth rate of the resistant subpopulation.",
      "Sibling model fitted to the dynamic hollow fiber data: Kroemer_2024_ceftazidime_avibactam_fosfomycin_hfim."
    )
  )

  ini({
    # --- Structural growth model (Table S3, "Structural model parameters") ---
    log10inoc_s <- 6.86; label("Inoculum of the susceptible subpopulation (log10 CFU/mL)")   # Table S3: 6.86 [6.74-6.96]
    log10inoc_r <- 3.14; label("Inoculum of the resistant subpopulation (log10 CFU/mL)")     # Table S3: 3.14 [2.81-3.43]
    log10bmax   <- 8.92; label("Maximum bacterial capacity Bmax (log10 CFU/mL)")             # Table S3: 8.92 [8.76-9.06]
    kgs         <- 1.81; label("Growth rate of the susceptible subpopulation (1/h)")         # Table S3: 1.81 [1.61-2.15]
    kgr         <- 0.45; label("Growth rate of the resistant subpopulation (1/h)")           # Table S3: 0.45 [0.41-0.52]

    # --- Mono-drug PD, susceptible subpopulation (Table S3, "Mono drug PD parameters") ---
    # Ceftazidime and avibactam use the sigmoidal Emax form (Eq 3); fosfomycin uses the power form (Eq 4).
    emax_caz_s  <- 3.40;  label("Emax of ceftazidime on the susceptible subpopulation (1/h)")   # Table S3: 3.40 [3.08-3.87]
    ec50_caz_s  <- 5.31;  label("EC50 of ceftazidime on the susceptible subpopulation (mg/L)")  # Table S3: 5.31 [4.23-6.54]
    hill_caz_s  <- 2.32;  label("Hill factor of ceftazidime on the susceptible subpopulation")  # Table S3: 2.32 [1.77-2.88]
    emax_avi_s  <- 3.3;   label("Emax of avibactam on the susceptible subpopulation (1/h)")     # Table S3: 3.3 [2.76-3.98]
    ec50_avi_s  <- 22.3;  label("EC50 of avibactam on the susceptible subpopulation (mg/L)")    # Table S3: 22.3 [16.50-29.10]
    hill_avi_s  <- 1.13;  label("Hill factor of avibactam on the susceptible subpopulation")    # Table S3: 1.13 [0.92-1.39]
    slope_fof_s <- 2.71;  label("Power-model slope of fosfomycin on the susceptible subpopulation (L/(mg*h))") # Table S3: 2.71 [2.51-3.09]
    hill_fof_s  <- 0.333; label("Power-model exponent of fosfomycin on the susceptible subpopulation")         # Table S3: 0.333 [0.292-0.377]

    # --- Mono-drug PD, resistant subpopulation (Table S3) ---
    emax_caz_r  <- 0.659;  label("Emax of ceftazidime on the resistant subpopulation (1/h)")   # Table S3: 0.659 [0.592-0.746]
    ec50_caz_r  <- 74.40;  label("EC50 of ceftazidime on the resistant subpopulation (mg/L)")  # Table S3: 74.40 [64.25-92.04]
    hill_caz_r  <- 8.45;   label("Hill factor of ceftazidime on the resistant subpopulation")  # Table S3: 8.45 [5.98-10.95]
    emax_fof_r  <- 0.635;  label("Emax of fosfomycin on the resistant subpopulation (1/h)")    # Table S3: 0.635 [0.582-0.718]
    ec50_fof_r  <- 4.70;   label("EC50 of fosfomycin on the resistant subpopulation (mg/L)")   # Table S3: 4.70 [3.67-5.72]
    hill_fof_r  <- 4.08;   label("Hill factor of fosfomycin on the resistant subpopulation")   # Table S3: 4.08 [2.76-5.79]
    slope_avi_r <- 0.0787; label("Power-model slope of avibactam on the resistant subpopulation (L/(mg*h))") # Table S3: 0.0787 [0.0554-0.1101]
    hill_avi_r  <- 0.317;  label("Power-model exponent of avibactam on the resistant subpopulation")         # Table S3: 0.317 [0.194-0.420]

    # --- GPDI: avibactam shifts the ceftazidime EC50 (Table S3, "Interaction model: avibactam affecting ceftazidime") ---
    # Table S3 footnote 1: the INT parameters were estimated on a log scale, back-transformed as TV = exp(theta) - 1.
    # Table S3 footnote 2: the interaction EC50s were estimated on a log scale, back-transformed as TV = exp(theta).
    lint_avi_caz_s     <- -6.70;  label("Log-scale INT: maximum change of the ceftazidime EC50 on (S) mediated by avibactam; INT = exp(theta) - 1")  # Table S3: -6.70 [-8.13 to -5.65]
    lec50int_avi_caz_s <- -16.20; label("Log EC50 of avibactam in the interaction on the ceftazidime EC50 on (S) (log mg/L)")                        # Table S3: -16.20 [-18.93 to -13.99]
    hillint_avi_caz_s  <- 0.266;  label("Hill factor of avibactam in the interaction on the ceftazidime EC50 on (S)")                                # Table S3: 0.266 [0.226-0.312]
    lint_avi_caz_r     <- -13.50; label("Log-scale INT: maximum change of the ceftazidime EC50 on (R) mediated by avibactam; INT = exp(theta) - 1")  # Table S3: -13.50 [-19.43 to -9.69]
    lec50int_avi_caz_r <- -5.23;  label("Log EC50 of avibactam in the interaction on the ceftazidime EC50 on (R) (log mg/L)")                        # Table S3: -5.23 [-5.53 to -5.04]
    hillint_avi_caz_r  <- fixed(1); label("Hill factor of avibactam in the interaction on the ceftazidime EC50 on (R)")                      # Table S3: 1, footnote 3 "parameter was fixed to a constant"

    # --- GPDI: ceftazidime shifts the fosfomycin EC50 on (R) (Table S3, "Interaction model: ceftazidime affecting fosfomycin") ---
    # Text S5: a monodirectional interaction (ceftazidime affecting fosfomycin) was best (dAIC -271.991 vs Bliss independence).
    lint_caz_fof_r     <- -7.85;  label("Log-scale INT: maximum change of the fosfomycin EC50 on (R) mediated by ceftazidime; INT = exp(theta) - 1") # Table S3: -7.85 [-10.08 to -5.58]
    lec50int_caz_fof_r <- -12.40; label("Log EC50 of ceftazidime in the interaction on the fosfomycin EC50 on (R) (log mg/L)")                       # Table S3: -12.40 [-15.01 to -10.50]
    hillint_caz_fof_r  <- 0.239;  label("Hill factor of ceftazidime in the interaction on the fosfomycin EC50 on (R)")                               # Table S3: 0.239 [0.179-0.311]

    # --- Variability model (Table S3, "Variability model") ---
    # Table S3 footnote 4: %CV = sqrt(exp(omega^2) - 1) * 100%, i.e. omega^2 = log(CV^2 + 1).
    etalog10inoc_r ~ 0.108787   # Table S3: inter-experimental variability on the inoculum of (R) = 33.9 %CV [27.8-38.9]; log(0.339^2 + 1) = 0.108787
    etakgr         ~ 0.049810   # Table S3: inter-experimental variability on the growth rate of (R) = 22.6 %CV [19.9-25.8]; log(0.226^2 + 1) = 0.049810

    addSd <- 1.30; label("Additive residual SD on the log10 bacterial count (log10 CFU/mL)") # Table S3: additive residual variability sigma = 1.30 [1.20-1.38]
  })

  model({
    # ---- 0. Guard against log10() of a non-positive solver excursion ----
    eps <- 1e-30

    # ---- 1. Back-transform the log-scale growth parameters (Table S3 units) ----
    # Inter-experimental variability enters as an exponential coefficient of
    # variation on the tabulated log10 inoculum and on the growth rate of (R)
    # (Text S5: "implemented on both parameters as exponential coefficient of
    # variation"). See the vignette Assumptions section: the paper does not
    # state whether the exponential eta multiplies the log10-scale inoculum or
    # the linear CFU count; the log10-scale reading is used because only it
    # reproduces the ~24-h spread in the timing of resistance emergence that the
    # variability model was introduced to capture.
    inoc_s <- 10^log10inoc_s                          # CFU/mL
    inoc_r <- 10^(log10inoc_r * exp(etalog10inoc_r))  # CFU/mL
    bmax   <- 10^log10bmax                            # CFU/mL
    kgr_i  <- kgr * exp(etakgr)                       # 1/h

    # ---- 2. GPDI terms (Text S3 Eq 8) ----
    #   Theta_GPDI = Theta * (1 + INT * C^H_INT / (EC50_INT^H_INT + C^H_INT))
    # Both INT and EC50_INT were estimated on the log scale (Table S3 footnotes 1-2).
    int_avi_caz_s     <- exp(lint_avi_caz_s) - 1
    ec50int_avi_caz_s <- exp(lec50int_avi_caz_s)
    int_avi_caz_r     <- exp(lint_avi_caz_r) - 1
    ec50int_avi_caz_r <- exp(lec50int_avi_caz_r)
    int_caz_fof_r     <- exp(lint_caz_fof_r) - 1
    ec50int_caz_fof_r <- exp(lec50int_caz_fof_r)

    # Avibactam potentiates ceftazidime by shifting its EC50 on both subpopulations.
    ec50_caz_s_eff <- ec50_caz_s *
      (1 + int_avi_caz_s * CONC_AVI_MGL^hillint_avi_caz_s /
         (ec50int_avi_caz_s^hillint_avi_caz_s + CONC_AVI_MGL^hillint_avi_caz_s))
    ec50_caz_r_eff <- ec50_caz_r *
      (1 + int_avi_caz_r * CONC_AVI_MGL^hillint_avi_caz_r /
         (ec50int_avi_caz_r^hillint_avi_caz_r + CONC_AVI_MGL^hillint_avi_caz_r))
    # Ceftazidime potentiates fosfomycin by shifting its EC50 on the resistant subpopulation.
    ec50_fof_r_eff <- ec50_fof_r *
      (1 + int_caz_fof_r * CONC_CAZ_MGL^hillint_caz_fof_r /
         (ec50int_caz_fof_r^hillint_caz_fof_r + CONC_CAZ_MGL^hillint_caz_fof_r))

    # ---- 3. Mono-drug effects (Text S3 Eq 3 sigmoidal Emax, Eq 4 power) ----
    e_caz_s <- emax_caz_s * CONC_CAZ_MGL^hill_caz_s /
      (ec50_caz_s_eff^hill_caz_s + CONC_CAZ_MGL^hill_caz_s)
    e_avi_s <- emax_avi_s * CONC_AVI_MGL^hill_avi_s /
      (ec50_avi_s^hill_avi_s + CONC_AVI_MGL^hill_avi_s)
    e_fof_s <- slope_fof_s * CONC_FOF_MGL^hill_fof_s

    e_caz_r <- emax_caz_r * CONC_CAZ_MGL^hill_caz_r /
      (ec50_caz_r_eff^hill_caz_r + CONC_CAZ_MGL^hill_caz_r)
    e_fof_r <- emax_fof_r * CONC_FOF_MGL^hill_fof_r /
      (ec50_fof_r_eff^hill_fof_r + CONC_FOF_MGL^hill_fof_r)
    e_avi_r <- slope_avi_r * CONC_AVI_MGL^hill_avi_r

    # ---- 4. Combined effects (Text S3 Eqs 5-7) ----
    # Susceptible subpopulation: fosfomycin is a power model, so the probabilistic
    # correction terms of Bliss independence are negligible and the criterion
    # collapses to simple addition of the three effects (Eq 7; Text S5).
    e_s <- e_caz_s + e_avi_s + e_fof_s

    # Resistant subpopulation: full three-drug Bliss independence (Eq 6),
    # normalised by the largest mono-drug Emax (Eq 5). Text S5 restricts that
    # normalisation to ceftazidime and fosfomycin: avibactam is described by a
    # power model on (R) and its effect at the highest applied concentration
    # (0.294 1/h) is small relative to Emax_CAZ (0.659) and Emax_FOF (0.635).
    emax_r_norm <- max(emax_caz_r, emax_fof_r)
    ea <- e_caz_r / emax_r_norm
    eb <- e_avi_r / emax_r_norm
    ec <- e_fof_r / emax_r_norm
    e_r <- (ea + eb + ec - ea * eb - ea * ec - eb * ec + ea * eb * ec) * emax_r_norm

    # ---- 5. Bacterial ODEs (Text S3 Eqs 1-2) ----
    # NOTE ON EQ 2: the supplement prints the kill term of the resistant
    # compartment as "- S * E_R" (and repeats it in the HFIM Eq 10). Taken
    # literally, the resistant subpopulation would be killed in proportion to
    # the SUSCEPTIBLE count, which at the estimated inocula (10^6.86 vs 10^3.14)
    # drives R negative within a fraction of an hour and leaves R uncontrollable
    # once S is eradicated. Fig. 3 of the main text shows the drug effects acting
    # ON the resistant bacteria ("dashed arrow, drug effects on the resistant
    # bacteria"), and the standard semi-mechanistic form (Nielsen 2011, ref 6 of
    # the supplement) is "- R * E_R". The typographical "S" is therefore read as
    # "R" here. Documented in the vignette Assumptions section.
    d/dt(bact_susceptible) <-
      bact_susceptible * kgs *
      (1 - (bact_susceptible + bact_resistant) / bmax) -
      bact_susceptible * e_s
    d/dt(bact_resistant) <-
      bact_resistant * kgr_i *
      (1 - (bact_susceptible + bact_resistant) / bmax) -
      bact_resistant * e_r

    bact_susceptible(0) <- inoc_s
    bact_resistant(0)   <- inoc_r

    # ---- 6. Output: total viable count on the log10 scale ----
    Cc <- log10(bact_susceptible + bact_resistant + eps)
    Cc ~ add(addSd)
  })
}
