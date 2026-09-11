Olivo_2025_vancomycin_invitro <- function() {
  description <- "In vitro (Staphylococcus aureus ATCC 43300, methicillin-resistant). Semi-mechanistic time-kill pharmacodynamic model of vancomycin against MRSA with adaptive resistance, in the lineage of Vera-Yunca 2020. Bacteria occupy an active, drug-susceptible state (gro) and a dormant, drug-insusceptible state (pers); transfer into the dormant state is density dependent at rate kad = (kgrow - kdeath) * (gro + pers) / nmax, and dormant bacteria return at kda. Vancomycin adds a sigmoidal Emax kill term to the death rate of the active state only. Adaptive resistance is carried by two states, aroff and aron, converted at rate kon * conc, and aron inflates the potency term linearly as ec50 = (1 + slopeAr * aron) * ec500, which is what produces the regrowth observed at 1 x MIC after 12 h. ec500 is proportional to the strain MIC (source Equation 4), so the fitted 1.05 mg/L at MIC 2 mg/L reproduces the paper's own 2.1 and 4.2 mg/L at MIC 4 and 8 mg/L. IMPORTANT: the bacterial states carry log10 CFU/mL, not CFU/mL -- see the extended note in ini() and the vignette Errata; nmax is therefore used directly as 8.94 rather than exponentiated. The paper's whole-body PBPK component was built in PK-Sim and is NOT reproduced here (only compound properties are tabulated; organ volumes, blood flows and partition coefficients are internal to the platform), so vancomycin exposure enters as the covariate CONC_VAN_MGL, which may be held static to replicate the time-kill experiment or driven with a tissue concentration-time profile to replicate the paper's coupled PBPK/PD simulations."
  reference <- "Ben Olivo L, Silva de Lemos JL, Rodrigues VJ, Kretschmer DB, Cruz WdA, Staudt KJ, Annaert P, Verlindo de Araujo B. PBPK/PD Model of Vancomycin in Sepsis: Linking Interstitial Exposure in Perfusion-Limited Tissues to MRSA Infection. Pharmaceutics. 2025 Aug 26;17(9):1111. doi:10.3390/pharmaceutics17091111. PMCID: PMC12473409. PD model structure: Supplementary File Section S1, Equations S1-S7, plus the drug-effect Equation (3) in Materials and methods 2.2 and the EC50-MIC relationship Equation (4) in Results 3.3. Parameter estimates with RSE and sampling-importance-resampling 95% CIs: Table 4. Observed time-kill data and model fit: Figure 3. Visual predictive check: Supplementary Figure S4. The adaptive-resistance structure is adopted from Vera-Yunca D, Girard P, Parra-Guillen ZP, Munafo A, Ottinger S, Terranova N. Machine learning and quantitative systems pharmacology to predict the effect of vancomycin. (Source reference 24.) The EC50-MIC interrelationship is taken from Schmidt S, Barbour A, Sahre M, Rand KH, Derendorf H. PK/PD: new insights for antibacterial and antiviral applications. Curr Opin Pharmacol. 2008 (source reference 25)."
  vignette <- "Olivo_2025_vancomycin"
  units <- list(
    time = "h",
    dosing = "mg/L (static or time-varying covariate, not an administered event)",
    concentration = "log10 CFU/mL for the model observation Cc"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. All four states are latent PD bookkeeping states of an
  # in vitro time-kill system, so none has a biological specimen.
  compartmentData <- list(
    gro = list(analyte = "Staphylococcus aureus ATCC 43300, active (growing, vancomycin-susceptible) subpopulation", units = "log10 CFU/mL", specimen = "not applicable", verified = TRUE),
    pers = list(analyte = "Staphylococcus aureus ATCC 43300, dormant (non-replicating, vancomycin-insusceptible) subpopulation", units = "log10 CFU/mL", specimen = "not applicable", verified = TRUE),
    aroff = list(analyte = "adaptive-resistance OFF subpopulation fraction", units = "fraction", specimen = "not applicable", verified = TRUE),
    aron = list(analyte = "adaptive-resistance ON subpopulation fraction", units = "fraction", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    CONC_VAN_MGL = list(
      description = "Vancomycin concentration driving the bactericidal effect and the adaptive-resistance switch (mg/L).",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "In the source time-kill experiment this is a static bath concentration: Materials and methods 2.2 states that 'the concentrations used during time-kill experiments represent 0.25, 0.5, 1, 2, 4, 6, and 8 times the MIC for the antimicrobial, alongside a growth control', i.e. 0 (control), 0.5, 1, 2, 4, 8, 12 and 16 mg/L at the fitted strain's MIC of 2 mg/L.",
        "Unlike the other members of the CONC_<drug>_MGL family this covariate is NOT restricted to a static value. Section 2.3 couples the same PD model to PBPK-predicted unbound interstitial concentration-time profiles in kidney, liver, lung and subcutis (Figures 4 and 5), so a time-varying series is the intended second use. Both uses are exercised in the validation vignette.",
        "This is the free (unbound) concentration. In the protein-free Mueller-Hinton broth of the time-kill assay the total bath concentration is the free concentration; in the coupled PBPK/PD simulations the driver is explicitly the unbound interstitial concentration (Figure 1 legend: 'C is the unbound drug concentration in interstitial compartment')."
      ),
      source_name = "Conc (Supplementary Equations S4-S5); C (Equation 3)"
    )
  )

  population <- list(
    species = "in vitro (Staphylococcus aureus ATCC 43300, methicillin-resistant)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    organism = "Methicillin-resistant Staphylococcus aureus ATCC 43300; broth-microdilution MIC of vancomycin 2 mg/L, classified susceptible (Results 3.3)",
    system = "Static time-kill curves in sterile flasks containing 20 mL Mueller-Hinton broth inoculated with 100 uL of bacterial suspension, pre-incubated 3 h 30 min to reach exponential (log) phase; viable counts at 0, 1, 2, 4, 6, 8, 10, 12 and 24 h, in triplicate per concentration",
    medium = "Mueller-Hinton broth",
    temperature = "35 C",
    duration = "24 h",
    starting_inoculum = "Approximately 10^7.2 CFU/mL; the paper does not tabulate the inoculum, so this value was read from the earliest observations across the eight arms of Figure 3 (see the bact0 note in ini())",
    mic_values = c(`Staphylococcus aureus ATCC 43300` = "2 mg/L"),
    concentration_range = "0 (growth control) and 0.25, 0.5, 1, 2, 4, 6, 8 x MIC",
    disease_state = "not applicable (in vitro)",
    notes = paste(
      "Model development was done in NONMEM 7.4 with PsN 4.9.0; robustness was assessed by sampling importance resampling (n = 1000), whose medians and 95% CIs are the last column of Table 4.",
      "The paper additionally reports simulations for hypothetical MRSA strains with MICs of 4 and 8 mg/L, obtained by scaling EC50 through Equation (4); those are reproduced here by changing the `mic` parameter rather than by refitting.",
      "The whole-body PBPK half of the paper (PK-Sim, Open Systems Pharmacology Suite 11.0, healthy volunteers and a 100-subject virtual septic population) is deliberately not part of this model file; see the description field and the vignette Errata."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # SCALE OF THE BACTERIAL STATES -- load-bearing, read before editing.
    #
    # Equations S1-S3 are written as though the states A and D were counts in
    # CFU/mL, which is the usual convention for this model family (compare
    # HernandezLozano_2025_apramycin_invitro.R, which exponentiates its
    # tabulated log10 Bmax). Under that reading the parameters of Table 4 are
    # not merely a poor fit, they are structurally impossible: the maximum
    # attainable kill rate Emax = 0.061 /h is SMALLER than the net growth rate
    # kgrow - kdeath = 1.96 - 1.85 = 0.11 /h, so d(gro)/dt > 0 at every
    # concentration and the model can never kill anything -- yet Figure 3
    # shows bacterial decline in five of its eight arms.
    #
    # The states are instead in log10 CFU/mL, which is what the unit column of
    # Table 4 says for the capacity parameter ("Maximum bacteria, Log CFU/mL,
    # 8.94"): nmax is used directly, not exponentiated, and the system is
    # logistic in log space. Three independent checks pass under this reading
    # and fail under the CFU reading:
    #   1. Growth control curvature. Observed control rises steeply early
    #      (about 0.175 log10/h between 2 and 6 h) and flattens late (about
    #      0.025 log10/h between 12 and 24 h). A CFU-scale model with net
    #      growth 0.11 /h cannot exceed 0.11/ln(10) = 0.048 log10/h at any
    #      time. The log-scale form gives 0.11 * A * (1 - A/8.94), i.e. 0.133
    #      log10/h at A = 7.5 and 0.026 log10/h at A = 8.7.
    #   2. Kill plateau. Saturating concentrations drive the system to
    #      nmax * (1 - Emax/(kgrow - kdeath)) = 8.94 * 0.4455 = 3.98 log10
    #      CFU/mL; the 4, 6 and 8 x MIC arms of Figure 3 converge on about 4.6
    #      and flatten.
    #   3. Regrowth at 1 x MIC. Reproduced at the observed time (dip to about
    #      6.2 by 8-12 h, recovery to about 7.5 by 24 h) with no retuning.
    # Every value below is the printed value; nothing was fitted to Figure 3.
    # ---------------------------------------------------------------------

    # ---- Bacterial system (Table 4) --------------------------------------
    lkgrow <- log(1.96); label("Log growth rate constant of the active bacterial state (1/h)")   # Table 4: k growth = 1.96 1/h, RSE 1.8%, SIR 1.95 (1.89-2.01)
    lkdeath <- log(1.85); label("Log natural death rate constant, both bacterial states (1/h)")  # Table 4: k death = 1.85 1/h, RSE 1.4%, SIR 1.85 (1.80-1.89). Methods 2.2: "a natural death rate in both states"
    nmax <- 8.94; label("Maximum bacterial density (log10 CFU/mL)")                              # Table 4: Maximum bacteria = 8.94 Log CFU/mL, RSE 1.8%, SIR 8.93 (8.66-9.21)
    lkda <- log(0.022); label("Log transfer rate constant from the dormant to the active state (1/h)") # Table 4: k DA = 0.022 1/h, RSE 3.7%, SIR 0.022 (0.020-0.023)

    # ---- Vancomycin effect (Table 4, Equation 3) -------------------------
    lemax <- log(0.061); label("Log maximum vancomycin kill rate constant (1/h)")                # Table 4: Maximum effect = 0.061 1/h, RSE 0.5%, SIR 0.062 (0.060-0.065). Discussion: "an Emax of 0.061 h-1 ... resulted in a half-kill time of 11 h"
    lec50Ref <- log(1.05); label("Log EC50 in the absence of adaptive resistance, at the reference MIC (mg/L)") # Table 4: EC50 in absence of AR = 1.05 mg/L, RSE 0.7%, SIR 1.05 (1.03-1.06)
    hill <- 5.74; label("Hill coefficient of the sigmoidal vancomycin effect (unitless)")        # Table 4: Hill factor = 5.74, RSE 3.4%, SIR 5.72 (5.41-6.05)

    # ---- Adaptive resistance (Table 4, Equations S4-S7) ------------------
    lkon <- log(0.021); label("Log adaptive-resistance activation rate constant (L/(mg*h))")     # Table 4: k ON = 0.021 1/h, RSE 3.7%, SIR 0.022 (0.021-0.023). Equations S4-S5 multiply kON by a concentration, so the dimension is L/(mg*h); Table 4's "h-1" omits the concentration term
    slopeAr <- 3.24; label("Linear slope of adaptive resistance on EC50 (unitless)")             # Table 4: Slope = 3.24, RSE 1.4%, SIR 3.23 (3.17-3.31). Equation S6: AReff = 1 + Slope * ARon

    # ---- Strain MIC (Results 3.3, Equation 4) ----------------------------
    # Equation (4) is MIC = (d/(Emax - d))^(1/gamma) * EC50, in which d, Emax
    # and gamma are strain-independent assay constants; EC50 is therefore
    # exactly proportional to MIC. Holding the ratio at the fitted strain and
    # rescaling reproduces the paper's own numbers: 1.05 * 4/2 = 2.1 and
    # 1.05 * 8/2 = 4.2 mg/L, which are the values quoted in Results 3.3 for
    # MICs of 4 and 8 mg/L. `mic` is exposed so a user can switch strain with
    # rxSolve(params = c(mic = 4)) without refitting.
    micRef <- fixed(2); label("MIC of vancomycin for the fitted strain, ATCC 43300 (mg/L)")      # Results 3.3: "The MIC of VAN for the MRSA strain was determined to be 2 ug/mL"
    mic <- fixed(2); label("MIC of vancomycin for the strain being simulated (mg/L)")            # Results 3.3: simulations reported for MICs of 2, 4 and 8 ug/mL

    # ---- Initial condition -----------------------------------------------
    # FIGURE-DERIVED, NOT TABULATED. Table 4 has no inoculum row and the
    # Methods give only the preparation ("100 uL of bacterial suspension" into
    # 20 mL of broth, 3 h 30 min to log phase). 7.2 log10 CFU/mL is the mean
    # of the earliest observed counts across the eight arms of Figure 3, whose
    # panels start between about 6.8 and 7.5. See the vignette Errata.
    bact0 <- 7.2; label("Initial density of the active bacterial state (log10 CFU/mL) -- figure-derived")  # Figure 3 (digitised); not reported in Table 4 or the Supplementary File

    # ---- Residual error ---------------------------------------------------
    # Table 4's last row is "Proportional error | % | 0.16 | RSE 13.6%". Read
    # literally as a 0.16% proportional error on a log10 count of about 7 it
    # would imply a residual SD of 0.011 log10 CFU/mL, and read as 16% it
    # would imply 1.1 log10 CFU/mL; the Figure S4 visual predictive check
    # shows a 10th-to-90th percentile band about 1.0 log10 CFU/mL wide, i.e.
    # a residual SD of 1.0/(2*1.2816) = 0.39. The value is therefore the
    # NONMEM $SIGMA VARIANCE, 0.16, on the log10 CFU/mL observation, giving
    # SD = sqrt(0.16) = 0.4 -- which matches the VPC band and the visible
    # spread of the Figure 3 triplicates. Encoded as additive because the
    # observation is already a log10 count. See the vignette Errata.
    addSd <- 0.4; label("Additive residual SD on the log10 bacterial count (log10 CFU/mL)")      # Table 4: "Proportional error" = 0.16 read as the $SIGMA variance; SD = sqrt(0.16) = 0.4
  })

  model({
    # ---- 1. Back-transform ------------------------------------------------
    kgrow <- exp(lkgrow)
    kdeath <- exp(lkdeath)
    kda <- exp(lkda)
    kon <- exp(lkon)
    emax <- exp(lemax)

    # ---- 2. Potency at this strain's MIC (Equation 4) ---------------------
    ec500 <- exp(lec50Ref) * mic / micRef

    # ---- 3. Density-dependent transfer into the dormant state (Eq S3) -----
    kad <- (kgrow - kdeath) * (gro + pers) / nmax

    # ---- 4. Adaptive resistance on potency (Eq S6, S7) --------------------
    areff <- 1 + slopeAr * aron
    ec50 <- areff * ec500

    # ---- 5. Vancomycin effect (Equation 3) --------------------------------
    effect <- emax * CONC_VAN_MGL^hill / (ec50^hill + CONC_VAN_MGL^hill)

    # ---- 6. Bacterial system (Eq S1, S2) ----------------------------------
    # Equation S2 as printed reads "kAD * A - kd * A - kDA * D", i.e. the
    # natural death term is applied to A a second time rather than to D. That
    # is a typographical error: Materials and methods 2.2 states the model has
    # "a natural death rate in both states", and with the printed form the
    # total population obeys d(gro+pers)/dt = (kgrow - 2*kdeath) * gro < 0, so
    # even the untreated growth control would decay monotonically instead of
    # rising to nmax as it does in Figure 3. The corrected term used here also
    # places the stationary state at gro + pers = nmax * (kdeath + kda)/kdeath,
    # i.e. 1.012 * nmax, which is the carrying capacity the parameter is named
    # for. See the vignette Errata.
    d/dt(gro) <- kgrow * gro - kad * gro - (kdeath + effect) * gro + kda * pers
    d/dt(pers) <- kad * gro - kdeath * pers - kda * pers

    # ---- 7. Adaptive-resistance states (Eq S4, S5) ------------------------
    d/dt(aroff) <- -kon * CONC_VAN_MGL * aroff
    d/dt(aron) <- kon * CONC_VAN_MGL * aroff

    # ---- 8. Initial conditions --------------------------------------------
    # All bacteria start active; the dormant state fills through kad. The
    # adaptive-resistance pool starts wholly OFF so that Equation S6 gives
    # AReff = 1 and the potency starts at ec500.
    gro(0) <- bact0
    aroff(0) <- 1

    # ---- 9. Observation ---------------------------------------------------
    # The observed quantity is the total viable count in log10 CFU/mL, which
    # is the sum of the two states on the model's own scale -- the same sum
    # Equation S3 uses for the density-dependent transfer.
    Cc <- gro + pers
    Cc ~ add(addSd)
  })
}
