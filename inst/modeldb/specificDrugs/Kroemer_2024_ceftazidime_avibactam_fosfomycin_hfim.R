Kroemer_2024_ceftazidime_avibactam_fosfomycin_hfim <- function() {
  description <- "In vitro (Escherichia coli, clinical MDR isolate; CTX-M-15 + TEM-4 ESBL and OXA-244 carbapenemase). Semi-mechanistic PK/PD model of the dynamic hollow fiber infection model for ceftazidime/avibactam and fosfomycin: the static time-kill model (susceptible plus joint resistant subpopulation, sigmoidal Emax or power-model kill, Bliss independence, GPDI synergy) extended with two phenotypically less susceptible subpopulations growing at 3x MIC, whose emergence is suppressed by ceftazidime and fosfomycin through subpopulation synergy. Drug concentrations are dosable one-compartment states declining with the nominal hollow fiber half-life. Sibling model: Kroemer_2024_ceftazidime_avibactam_fosfomycin_tkc."
  reference <- "Kroemer N, Amann LF, Farooq A, Pfaffendorf C, Martens M, Decousser JW, Gregoire N, Nordmann P, Wicha SG. Pharmacokinetic/pharmacodynamic analysis of ceftazidime/avibactam and fosfomycin combinations in an in vitro hollow fiber infection model against multidrug-resistant Escherichia coli. Microbiol Spectr. 2024 Jan;12(1):e03318-23. doi:10.1128/spectrum.03318-23. Structural equations: supplemental material Text S3 Eqs 3-13. Parameter estimates: supplemental material Table S4 (with 95% confidence intervals from sampling importance resampling); parameters marked 'FIX to TKC parameter' are inherited from Table S3. Nominal hollow fiber pharmacokinetics: main-text Table 1."
  vignette <- "Kroemer_2024_ceftazidime_avibactam_fosfomycin"
  units <- list(time = "h", dosing = "mg/L", concentration = "mg/L")

  # Four bacterial subpopulations (Text S3 Eqs 9-12) plus the three antibiotic
  # concentration states of the hollow fiber central compartment. None maps onto
  # a canonical PK compartment name.
  paper_specific_compartments <- c(
    "bact_susceptible", "bact_resistant", "bact_resistant_cza", "bact_resistant_fof",
    "conc_caz", "conc_avi", "conc_fof"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    bact_susceptible   = list(analyte = "Escherichia coli (susceptible bacteria)", units = NA_character_, specimen = "not applicable", verified = FALSE),
    bact_resistant     = list(analyte = "Escherichia coli (resistant bacteria)", units = NA_character_, specimen = "not applicable", verified = FALSE),
    bact_resistant_cza = list(analyte = "Escherichia coli (resistant to ceftazidime bacteria)", units = NA_character_, specimen = "not applicable", verified = FALSE),
    bact_resistant_fof = list(analyte = "Escherichia coli (resistant to fosfomycin bacteria)", units = NA_character_, specimen = "not applicable", verified = FALSE),
    conc_caz           = list(analyte = "ceftazidime", units = NA_character_, specimen = "administration site", verified = FALSE),
    conc_avi           = list(analyte = "avibactam", units = NA_character_, specimen = "administration site", verified = FALSE),
    conc_fof           = list(analyte = "fosfomycin", units = NA_character_, specimen = "administration site", verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species = "in vitro (Escherichia coli, single clinical multidrug-resistant isolate)",
    n_subjects = 1L,
    n_studies = 1L,
    disease_state = "Clinical E. coli isolated from a rectal swab in a multidrug-resistance screening at Henri-Mondor Hospital, Paris. Whole-genome sequencing identified CTX-M-15 and TEM-4 (extended-spectrum beta-lactamase) and OXA-244 (OXA-48-like carbapenemase) genes. MICs: ceftazidime 16 mg/L, ceftazidime/avibactam 0.125 mg/L, fosfomycin 16/0.5 mg/L (microdilution/agar dilution). Retention of the CTX-M-15 and OXA-244 enzymes in the hollow fiber cartridge was confirmed by antigen testing after 2 h of preincubation.",
    model_system = "Dynamic hollow fiber infection model over 72 h at 37 degrees C in ambient air (FiberCell C2011 polysulfone cartridge, 200 mL central compartment, 40 mL extra-capillary space inoculated at 10^6 CFU/mL). Fresh Mueller-Hinton broth 2 was pumped in at 0.834-1.28 mL/min with matched removal to control the drug pharmacokinetics; the first dose was given by syringe driver after 2 h of preincubation. Total and 3x MIC phenotypically resistant bacterial counts were quantified by sampling, serial dilution and plating. Each hollow fiber experiment was performed once.",
    initial_inoculum = "10^6 CFU/mL targeted; susceptible subpopulation fixed at 10^6.86 CFU/mL from the static time-kill model, resistant subpopulation estimated at 10^2.89 CFU/mL",
    dose_range = "Ceftazidime/avibactam 0.06/0.015 to 2/0.5 g q8h and fosfomycin 0.125 to 6 g q8h, alone and in combination, mimicked as bolus injections reproducing the simulated Cmax and Cmin of Table 1",
    notes = paste(
      "Estimated in NONMEM 7.5.0 (LAPLACIAN-I, CVODES/ADVAN14); parameter uncertainty from sampling importance resampling in PsN 5.0.",
      "All drug-effect parameters and the growth constants of the susceptible and resistant subpopulations were fixed to the static time-kill estimates (Table S3); the hollow fiber data could not fully inform them and the estimates tended towards the static values.",
      "Pharmacokinetic profiles were not estimated: they were simulated from published clinical population PK models of ceftazidime/avibactam and fosfomycin (refs 21, 22 of the main text) and reproduced in the hollow fiber system as bolus injections matching Cmax and Cmin, with a joint elimination half-life of about 2 h across all three drugs (Table 1: 1.81-3.03 h). Bioanalysis by UHPLC-MS/MS confirmed the nominal profiles, so the nominal pharmacokinetics were used for modelling.",
      "Assumed unbound fractions were 92% for avibactam, 85% for ceftazidime and effectively 100% for fosfomycin (Methods, Pharmacokinetics); the model works with the total concentrations mimicked in the hollow fiber system.",
      "The two less susceptible subpopulations are the colonies growing on agar supplemented with 3x MIC of the respective drug; their inocula and growth are the only route to regrowth once the susceptible and resistant subpopulations are killed.",
      "Sibling model fitted to the static time-kill data: Kroemer_2024_ceftazidime_avibactam_fosfomycin_tkc."
    )
  )

  ini({
    # --- Structural growth model (Table S4, "Structural model parameters") ---
    log10inoc_s <- fixed(6.86); label("Inoculum of the susceptible subpopulation (log10 CFU/mL; the static time-kill estimate)") # Table S4: 6.86 FIX to TKC parameter
    log10inoc_r <- 2.89;        label("Inoculum of the resistant subpopulation (log10 CFU/mL)")                                          # Table S4: 2.89 [2.68-2.99]
    log10bmax   <- 9.73;        label("Maximum bacterial capacity Bmax (log10 CFU/mL)")                                                  # Table S4: 9.73 [9.50-9.91]
    kgs         <- fixed(1.81); label("Growth rate of the susceptible subpopulation (1/h; the static time-kill estimate)")      # Table S4: 1.81 FIX to TKC parameter
    kgr         <- fixed(0.45); label("Growth rate of the resistant subpopulation (1/h; the static time-kill estimate)")        # Table S4: 0.45 FIX to TKC parameter

    # --- Mono-drug PD, susceptible subpopulation (Table S4; all FIXED to Table S3) ---
    emax_caz_s  <- fixed(3.40);  label("Emax of ceftazidime on the susceptible subpopulation (1/h)")   # Table S4: 3.40 FIX to TKC parameter
    ec50_caz_s  <- fixed(5.31);  label("EC50 of ceftazidime on the susceptible subpopulation (mg/L)")  # Table S4: 5.31 FIX to TKC parameter
    hill_caz_s  <- fixed(2.32);  label("Hill factor of ceftazidime on the susceptible subpopulation") # Table S4: 2.32 FIX to TKC parameter
    emax_avi_s  <- fixed(3.3);   label("Emax of avibactam on the susceptible subpopulation (1/h)")     # Table S4: 3.3 FIX to TKC parameter
    ec50_avi_s  <- fixed(22.3);  label("EC50 of avibactam on the susceptible subpopulation (mg/L)")    # Table S4: 22.3 FIX to TKC parameter
    hill_avi_s  <- fixed(1.13);  label("Hill factor of avibactam on the susceptible subpopulation")   # Table S4: 1.13 FIX to TKC parameter
    slope_fof_s <- fixed(2.71);  label("Power-model slope of fosfomycin on the susceptible subpopulation (L/(mg*h))") # Table S4: 2.71 FIX to TKC parameter
    hill_fof_s  <- fixed(0.333); label("Power-model exponent of fosfomycin on the susceptible subpopulation")        # Table S4: 0.333 FIX to TKC parameter

    # --- Mono-drug PD, resistant subpopulation (Table S4; all FIXED to Table S3) ---
    emax_caz_r  <- fixed(0.659);  label("Emax of ceftazidime on the resistant subpopulation (1/h)")   # Table S4: 0.659 FIX to TKC parameter
    ec50_caz_r  <- fixed(74.40);  label("EC50 of ceftazidime on the resistant subpopulation (mg/L)")  # Table S4: 74.40 FIX to TKC parameter
    hill_caz_r  <- fixed(8.45);   label("Hill factor of ceftazidime on the resistant subpopulation") # Table S4: 8.45 FIX to TKC parameter
    emax_fof_r  <- fixed(0.635);  label("Emax of fosfomycin on the resistant subpopulation (1/h)")    # Table S4: 0.635 FIX to TKC parameter
    ec50_fof_r  <- fixed(4.70);   label("EC50 of fosfomycin on the resistant subpopulation (mg/L)")   # Table S4: 4.70 FIX to TKC parameter
    hill_fof_r  <- fixed(4.08);   label("Hill factor of fosfomycin on the resistant subpopulation")  # Table S4: 4.08 FIX to TKC parameter
    slope_avi_r <- fixed(0.0787); label("Power-model slope of avibactam on the resistant subpopulation (L/(mg*h))") # Table S4: 0.0787 FIX to TKC parameter
    hill_avi_r  <- fixed(0.317);  label("Power-model exponent of avibactam on the resistant subpopulation")        # Table S4: 0.317 FIX to TKC parameter

    # --- GPDI: avibactam shifts the ceftazidime EC50 (Table S4; all FIXED to Table S3) ---
    # Table S4 footnote 1: INT parameters estimated on a log scale, TV = exp(theta) - 1.
    # Table S4 footnote 2: interaction EC50s estimated on a log scale, TV = exp(theta).
    lint_avi_caz_s     <- fixed(-6.70);  label("Log-scale INT: maximum change of the ceftazidime EC50 on (S) mediated by avibactam; INT = exp(theta) - 1") # Table S4: -6.70 FIX to TKC parameter
    lec50int_avi_caz_s <- fixed(-16.20); label("Log EC50 of avibactam in the interaction on the ceftazidime EC50 on (S) (log mg/L)")                        # Table S4: -16.20 FIX to TKC parameter
    hillint_avi_caz_s  <- fixed(0.266);  label("Hill factor of avibactam in the interaction on the ceftazidime EC50 on (S)")                               # Table S4: 0.266 FIX to TKC parameter
    lint_avi_caz_r     <- fixed(-13.50); label("Log-scale INT: maximum change of the ceftazidime EC50 on (R) mediated by avibactam; INT = exp(theta) - 1") # Table S4: -13.50 FIX to TKC parameter
    lec50int_avi_caz_r <- fixed(-5.23);  label("Log EC50 of avibactam in the interaction on the ceftazidime EC50 on (R) (log mg/L)")                        # Table S4: -5.23 FIX to TKC parameter
    hillint_avi_caz_r  <- fixed(1);      label("Hill factor of avibactam in the interaction on the ceftazidime EC50 on (R)")                               # Table S4: 1 FIX to TKC parameter

    # --- GPDI: ceftazidime shifts the fosfomycin EC50 on (R) (Table S4; all FIXED to Table S3) ---
    lint_caz_fof_r     <- fixed(-7.85);  label("Log-scale INT: maximum change of the fosfomycin EC50 on (R) mediated by ceftazidime; INT = exp(theta) - 1") # Table S4: -7.85 FIX to TKC parameter
    lec50int_caz_fof_r <- fixed(-12.40); label("Log EC50 of ceftazidime in the interaction on the fosfomycin EC50 on (R) (log mg/L)")                        # Table S4: -12.40 FIX to TKC parameter
    hillint_caz_fof_r  <- fixed(0.239);  label("Hill factor of ceftazidime in the interaction on the fosfomycin EC50 on (R)")                               # Table S4: 0.239 FIX to TKC parameter

    # --- Less susceptible (3x MIC) subpopulation model (Table S4, "Less susceptible subpopulation model") ---
    log10inoc_rcza <- fixed(-18); label("Inoculum of the ceftazidime/avibactam less susceptible subpopulation (log10 CFU/mL; the final estimate)") # Table S4: -18, footnote 3 "parameter was fixed to final estimate"
    log10inoc_rfof <- -2.15;      label("Inoculum of the fosfomycin less susceptible subpopulation (log10 CFU/mL)")                                         # Table S4: -2.15 [-2.98 to -1.52]
    kgr2           <- 2.37;       label("Growth rate of both less susceptible subpopulations (1/h)")                                                        # Table S4: 2.37 [2.09-2.68]; merged to one parameter for both subpopulations (Text S5)

    ec50_cza_rcza  <- 0.576;      label("EC50 of ceftazidime suppressing the CZA less susceptible subpopulation (mg/L)")  # Table S4: 0.576 [0.441-0.765]
    hill_cza_rcza  <- fixed(1);   label("Hill factor of ceftazidime suppressing the CZA less susceptible subpopulation")   # Table S4: 1, footnote 4 "parameter was fixed to a constant"
    ec50_fof_rcza  <- 1.38;       label("EC50 of fosfomycin suppressing the CZA less susceptible subpopulation (mg/L)")   # Table S4: 1.38 [1.00-2.49]
    hill_fof_rcza  <- fixed(20);  label("Hill factor of fosfomycin suppressing the CZA less susceptible subpopulation")    # Table S4: 20, footnote 4 "parameter was fixed to a constant" (empirically fixed for very steep concentration-effect relations)
    ec50_fof_rfof  <- 6.84;       label("EC50 of fosfomycin suppressing the FOF less susceptible subpopulation (mg/L)")   # Table S4: 6.84 [6.48-7.17]
    hill_fof_rfof  <- fixed(20);  label("Hill factor of fosfomycin suppressing the FOF less susceptible subpopulation")    # Table S4: 20, footnote 4 "parameter was fixed to a constant"
    ec50_cza_rfof  <- 0.049;      label("EC50 of ceftazidime suppressing the FOF less susceptible subpopulation (mg/L)")  # Table S4: 0.049 [0.040-0.057]
    hill_cza_rfof  <- 2.49;       label("Hill factor of ceftazidime suppressing the FOF less susceptible subpopulation")  # Table S4: 2.49 [1.76-4.20]

    # --- Nominal hollow fiber pharmacokinetics (main-text Table 1; not estimated) ---
    thalf <- fixed(1.81); label("Joint elimination half-life of ceftazidime, avibactam and fosfomycin in the hollow fiber system (h)") # Table 1: modal simulated half-life across the 14 hollow fiber experiments (range 1.81-3.03 h); Discussion: "a then-joint elimination half-life of approximately 2 h"

    # --- Variability model (Table S4, "Variability model") ---
    # Table S4 footnote 5: %CV = sqrt(exp(omega^2) - 1) * 100%, i.e. omega^2 = log(CV^2 + 1).
    # Text S5: the static-model IIV on the growth rate of (R) was not needed for the hollow fiber data.
    etalog10inoc_r    ~ 0.208826  # Table S4: inter-experimental variability on the inoculum of (R) = 48.2 %CV [40.4-56.9]; log(0.482^2 + 1) = 0.208826
    etalog10inoc_rcza ~ 0.204927  # Table S4: inter-experimental variability on the inoculum of the CZA less susceptible subpopulation = 47.7 %CV [28.6-65.9]; log(0.477^2 + 1) = 0.204927

    addSd          <- 3.28;  label("Additive residual SD on the log10 total bacterial count (log10 CFU/mL)")               # Table S4: additive residual variability on the total bacterial count sigma = 3.28 [2.98-3.67]
    addSd_CFUrcza  <- 0.906; label("Additive residual SD on the log10 CZA less susceptible count (log10 CFU/mL)")          # Table S4: additive residual variability on less susceptible sigma = 0.906 [0.837-1.00]; ONE shared parameter was estimated for both less susceptible subpopulations (Text S5), reproduced here as two identical entries because rxode2 requires one residual parameter per endpoint
    addSd_CFUrfof  <- 0.906; label("Additive residual SD on the log10 FOF less susceptible count (log10 CFU/mL)")          # Table S4: additive residual variability on less susceptible sigma = 0.906 [0.837-1.00]; same shared estimate as addSd_CFUrcza
  })

  model({
    # ---- 0. Guard against log10() of a non-positive solver excursion ----
    eps <- 1e-30

    # ---- 1. Back-transform the log-scale growth parameters ----
    # Inter-experimental variability enters as an exponential coefficient of
    # variation on the tabulated log10 inocula (Text S5). See the vignette
    # Assumptions section for why the eta is applied to the log10-scale value
    # rather than to the linear CFU count.
    inoc_s    <- 10^log10inoc_s                                  # CFU/mL
    inoc_r    <- 10^(log10inoc_r * exp(etalog10inoc_r))          # CFU/mL
    inoc_rcza <- 10^(log10inoc_rcza * exp(etalog10inoc_rcza))    # CFU/mL
    inoc_rfof <- 10^log10inoc_rfof                               # CFU/mL
    bmax      <- 10^log10bmax                                    # CFU/mL

    # ---- 2. Nominal hollow fiber pharmacokinetics (main-text Table 1) ----
    # A single central compartment per drug: doses are given as bolus increments
    # of concentration (mg/L) and decline with the joint half-life.
    kel <- log(2) / thalf                                        # 1/h

    # ---- 3. GPDI terms (Text S3 Eq 8) ----
    #   Theta_GPDI = Theta * (1 + INT * C^H_INT / (EC50_INT^H_INT + C^H_INT))
    int_avi_caz_s     <- exp(lint_avi_caz_s) - 1
    ec50int_avi_caz_s <- exp(lec50int_avi_caz_s)
    int_avi_caz_r     <- exp(lint_avi_caz_r) - 1
    ec50int_avi_caz_r <- exp(lec50int_avi_caz_r)
    int_caz_fof_r     <- exp(lint_caz_fof_r) - 1
    ec50int_caz_fof_r <- exp(lec50int_caz_fof_r)

    ec50_caz_s_eff <- ec50_caz_s *
      (1 + int_avi_caz_s * conc_avi^hillint_avi_caz_s /
         (ec50int_avi_caz_s^hillint_avi_caz_s + conc_avi^hillint_avi_caz_s))
    ec50_caz_r_eff <- ec50_caz_r *
      (1 + int_avi_caz_r * conc_avi^hillint_avi_caz_r /
         (ec50int_avi_caz_r^hillint_avi_caz_r + conc_avi^hillint_avi_caz_r))
    ec50_fof_r_eff <- ec50_fof_r *
      (1 + int_caz_fof_r * conc_caz^hillint_caz_fof_r /
         (ec50int_caz_fof_r^hillint_caz_fof_r + conc_caz^hillint_caz_fof_r))

    # ---- 4. Mono-drug kill of the susceptible and resistant subpopulations ----
    #        (Text S3 Eq 3 sigmoidal Emax, Eq 4 power)
    e_caz_s <- emax_caz_s * conc_caz^hill_caz_s /
      (ec50_caz_s_eff^hill_caz_s + conc_caz^hill_caz_s)
    e_avi_s <- emax_avi_s * conc_avi^hill_avi_s /
      (ec50_avi_s^hill_avi_s + conc_avi^hill_avi_s)
    e_fof_s <- slope_fof_s * conc_fof^hill_fof_s

    e_caz_r <- emax_caz_r * conc_caz^hill_caz_r /
      (ec50_caz_r_eff^hill_caz_r + conc_caz^hill_caz_r)
    e_fof_r <- emax_fof_r * conc_fof^hill_fof_r /
      (ec50_fof_r_eff^hill_fof_r + conc_fof^hill_fof_r)
    e_avi_r <- slope_avi_r * conc_avi^hill_avi_r

    # ---- 5. Combined kill (Text S3 Eqs 5-7) ----
    # Susceptible: Bliss independence collapses to effect addition (Eq 7).
    e_s <- e_caz_s + e_avi_s + e_fof_s
    # Resistant: three-drug Bliss independence (Eq 6) normalised by the largest
    # mono-drug Emax (Eq 5), restricted by Text S5 to ceftazidime and fosfomycin.
    emax_r_norm <- max(emax_caz_r, emax_fof_r)
    ea <- e_caz_r / emax_r_norm
    eb <- e_avi_r / emax_r_norm
    ec <- e_fof_r / emax_r_norm
    e_r <- (ea + eb + ec - ea * eb - ea * ec - eb * ec + ea * eb * ec) * emax_r_norm

    # ---- 6. Growth inhibition of the less susceptible subpopulations ----
    #        (Text S3 Eq 13: the sigmoidal Emax model simplified to Emax = 1,
    #        i.e. full inhibition of growth at saturating concentration)
    # Text S5: the inhibitory effects of ceftazidime and avibactam were merged for
    # the subpopulation synergy and estimated as a function of the CEFTAZIDIME
    # concentration alone, so conc_avi does not enter these four terms.
    i_cza_rcza <- conc_caz^hill_cza_rcza / (ec50_cza_rcza^hill_cza_rcza + conc_caz^hill_cza_rcza)
    i_fof_rcza <- conc_fof^hill_fof_rcza / (ec50_fof_rcza^hill_fof_rcza + conc_fof^hill_fof_rcza)
    i_fof_rfof <- conc_fof^hill_fof_rfof / (ec50_fof_rfof^hill_fof_rfof + conc_fof^hill_fof_rfof)
    i_cza_rfof <- conc_caz^hill_cza_rfof / (ec50_cza_rfof^hill_cza_rfof + conc_caz^hill_cza_rfof)

    # ---- 7. Bacterial ODEs (Text S3 Eqs 9-12) ----
    # NOTE ON EQ 10: as in the static model (Eq 2), the supplement prints the kill
    # term of the resistant compartment as "- S * E_R". Read as "- R * E_R" here:
    # taken literally the resistant subpopulation would be killed in proportion to
    # the susceptible count, which drives R negative within a fraction of an hour
    # and leaves R uncontrollable once S is eradicated. Fig. 3 shows the drug
    # effects acting on the resistant bacteria and the standard semi-mechanistic
    # form (Nielsen 2011, ref 6 of the supplement) is "- R * E_R". Documented in
    # the vignette Assumptions section.
    cfu_all <- bact_susceptible + bact_resistant + bact_resistant_cza + bact_resistant_fof
    cap <- 1 - cfu_all / bmax

    d/dt(bact_susceptible) <- bact_susceptible * kgs * cap - bact_susceptible * e_s
    d/dt(bact_resistant)   <- bact_resistant   * kgr * cap - bact_resistant   * e_r
    # Eqs 11-12: the less susceptible subpopulations are not killed, only
    # growth-inhibited, multiplicatively by both drugs (subpopulation synergy).
    d/dt(bact_resistant_cza) <- bact_resistant_cza * kgr2 * cap * (1 - i_cza_rcza) * (1 - i_fof_rcza)
    d/dt(bact_resistant_fof) <- bact_resistant_fof * kgr2 * cap * (1 - i_fof_rfof) * (1 - i_cza_rfof)

    # ---- 8. Antibiotic concentration states (mg/L) ----
    d/dt(conc_caz) <- -kel * conc_caz
    d/dt(conc_avi) <- -kel * conc_avi
    d/dt(conc_fof) <- -kel * conc_fof

    # ---- 9. Initial conditions ----
    bact_susceptible(0)   <- inoc_s
    bact_resistant(0)     <- inoc_r
    bact_resistant_cza(0) <- inoc_rcza
    bact_resistant_fof(0) <- inoc_rfof

    # ---- 10. Outputs (log10 CFU/mL) ----
    # Total count = S + R + R_CZA + R_FOF (Fig. 3 caption).
    Cc <- log10(cfu_all + eps)
    # Text S5: "To account for the lower limit of quantification of the less
    # susceptible subpopulations a baseline was set to 0 log10(CFU/mL)". The
    # +1 CFU/mL inside the log10 realises that baseline smoothly: the output sits
    # at 0 while the subpopulation is below 1 CFU/mL and rises once it emerges,
    # so the timing of emergence is driven by the inoculum and the suppressing
    # drug concentrations, exactly as the paper describes.
    CFUrcza <- log10(bact_resistant_cza + 1)
    CFUrfof <- log10(bact_resistant_fof + 1)

    Cc      ~ add(addSd)
    CFUrcza ~ add(addSd_CFUrcza)
    CFUrfof ~ add(addSd_CFUrfof)
  })
}
