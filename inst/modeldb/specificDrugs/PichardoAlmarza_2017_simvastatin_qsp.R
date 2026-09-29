PichardoAlmarza_2017_simvastatin_qsp <- function() {
  description <- "QSP (pharmacokinetic/pharmacodynamic module). Simvastatin parent/metabolite pharmacokinetics coupled to an inhibitory turnover model of circulating LDL, forming the drug module of the hybrid multiscale quantitative-systems-pharmacology model of atherosclerosis progression and statin adherence of Pichardo-Almarza and Diaz-Zuccarini (2017). Simvastatin (lactone) is absorbed first order from a gut compartment into a single simvastatin compartment whose total apparent clearance CL2 splits 30 percent to formation of the active metabolite simvastatin acid and 70 percent to other elimination; simvastatin acid is eliminated first order from its own compartment. Circulating LDL is a zero-order-in, first-order-out turnover pool whose production Kin is inhibited by simvastatin acid through a fractional Emax function, with Kout defined as Kin divided by the LDL baseline so the drug-free pool sits exactly at baseline. Parameters are the typical-patient values of Table 1, originally estimated by Kim 2011 in healthy male volunteers. The multiscale arterial-wall module of the source (hemodynamics, endothelial LDL transport, oxidised LDL, monocytes, macrophages, foam cells and plaque growth) is NOT encoded here because its wall geometry, endothelial surface area, LDL wall diffusivity and per-species volume constants are not reported in any available source; see the validation vignette for the full gap list. The vignette also reproduces the paper's two-state Markov chain for medication adherence, which is a dosing-schedule construct rather than part of the differential-equation system."
  reference <- paste(
    "Pichardo-Almarza C, Diaz-Zuccarini V.",
    "Understanding the Effect of Statins and Patient Adherence in Atherosclerosis",
    "via a Quantitative Systems Pharmacology Model Using a Novel, Hybrid, and",
    "Multi-Scale Approach.",
    "Front Pharmacol. 2017;8:635. doi:10.3389/fphar.2017.00635.",
    "The 30/70 split of CL2 between metabolite formation and other elimination is",
    "read from the predecessor paper by the same authors, which prints the",
    "percentages on the same model diagram:",
    "Pichardo-Almarza C, Metcalf L, Finkelstein A, Diaz-Zuccarini V.",
    "Using a Systems Pharmacology Approach to Study the Effect of Statins on the",
    "Early Stage of Atherosclerosis in Humans.",
    "CPT Pharmacometrics Syst Pharmacol. 2015;4(1):41-50. doi:10.1002/psp4.7.",
    "Both papers attribute the PKPD parameter estimates to",
    "Kim J, Ahn BJ, Chae HS, Han S, Doh K, Choi J, et al.",
    "A population pharmacokinetic-pharmacodynamic model for simvastatin that",
    "predicts low-density lipoprotein-cholesterol reduction in patients with",
    "primary hyperlipidaemia.",
    "Basic Clin Pharmacol Toxicol. 2011;109(3):156-163.",
    "doi:10.1111/j.1742-7843.2011.00700.x.",
    sep = " "
  )
  vignette <- "PichardoAlmarza_2017_simvastatin_qsp"

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Pichardo-Almarza 2017 Figure 4 and
  # Pichardo-Almarza 2015 Figure 4 (which labels the same three PK
  # compartments 'Gut (Cmpt 1)', 'Simvastatin (cmpt 2)' and 'Simvastatin acid
  # (cmpt 3)') and Table 1.
  compartmentData <- list(
    depot = list(analyte = "simvastatin (lactone)", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "simvastatin (lactone)", units = "mg", specimen = "plasma", verified = TRUE),
    central_acid = list(analyte = "simvastatin acid", units = "mg", specimen = "plasma", verified = TRUE),
    ldl = list(analyte = "low-density lipoprotein particles", units = "nmol/L", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    # The source reports a single set of typical-patient PKPD values
    # (Pichardo-Almarza 2017 Table 1, 'parameter values for a Typical
    # Patient') with no covariate relationships on any PK or PD parameter.
    # Between-individual variability in the source's virtual population is
    # introduced on arterial-wall quantities (lumen radius, blood viscosity)
    # and on the circulating LDL level, not on the PKPD parameters; those
    # wall quantities belong to the module that is not encoded here.
  )

  population <- list(
    species = "human",
    n_subjects = 27L,
    disease_state = "Healthy adult male volunteers. Pichardo-Almarza 2015 Methods ('PKPD model of statins') states that the PKPD model 'was developed using data collected from 27 healthy male volunteers with a daily dose of 40 mg of simvastatin given for 14 d'; Pichardo-Almarza 2017 describes the same model as fitted 'in male volunteers'. Note that the originating publication (Kim 2011) is titled for prediction in patients with primary hyperlipidaemia, so the estimation cohort (healthy volunteers) and the intended prediction population (hyperlipidaemia) differ; only the two Pichardo-Almarza papers were available to confirm the estimation cohort.",
    sex_female_pct = 0,
    dose_range = "40 mg simvastatin orally once daily for 14 days.",
    sampling = "Simvastatin lactone (parent) and simvastatin acid (metabolite) plasma concentrations at 0, 0.5, 1, 1.5, 2, 3, 3.5, 4, 5, 6, 8, 10, 12 and 24 h after the dose on days 1, 7 and 14; LDL measured daily after an overnight fast (Pichardo-Almarza 2015 Methods).",
    notes = "The QSP paper's own simulated cohort (Pichardo-Almarza 2015 Table 3: 1,000 virtual subjects with log-normal lumen radius mean 0.03 dm, blood viscosity mean 0.004 Pa s, and circulating LDL mean 1,402 nmol/L, SD 141, 5th-95th percentile 1,190-1,644) varies only inputs to the arterial-wall module, which is not encoded here; no PKPD parameter is varied between individuals and none carries a covariate effect. The PKPD module is therefore deterministic: the source reports a single typical-patient parameter set with no inter-individual variability and no residual-error model."
  )

  ini({
    # Notes on parameter values:
    # * Every value is the typical-patient value of Pichardo-Almarza 2017
    #   Table 1 (identical to Pichardo-Almarza 2015 Table 1), which both
    #   papers attribute to Kim 2011. The source reports point values only --
    #   no standard errors, no confidence intervals, no inter-individual
    #   variances and no residual-error estimates -- so every parameter is
    #   encoded with fixed() and the residual error is fixed at zero.
    # * FM is NOT printed in Pichardo-Almarza 2017: its Figure 4 labels the
    #   simvastatin-to-simvastatin-acid arm 'CL23' and gives no value. The
    #   predecessor figure (Pichardo-Almarza 2015 Figure 4) draws the same two
    #   arms out of compartment 2 as 'CL2 (30%)' into simvastatin acid and
    #   'CL2 (70%)' to elimination, which identifies CL2 as the TOTAL
    #   clearance of compartment 2 and fixes FM = 0.30.
    # * EC50 is reported in ng/mL, so the metabolite concentration driving the
    #   PD must be in ng/mL; model() converts mg/L to ng/mL explicitly.

    # --- Pharmacokinetics: simvastatin (parent) ---
    # Ka is tabulated and the model diagram draws a gut compartment feeding
    # the simvastatin compartment through it, so the structure is kept. Note
    # that the source's own published PK figure (2015 Figure 2b) was NOT
    # generated through that path: its simvastatin curve rises vertically to
    # Dose/V2 at the dose time, i.e. the dose was placed directly in the
    # simvastatin compartment. The vignette quantifies both readings.
    lka <- fixed(log(2.76))
    label("Absorption rate constant from the gut compartment, Ka (1/h)") # 2017 Table 1: Ka = 2.76 1/hr
    lcl <- fixed(log(1740))
    label("Apparent total clearance of simvastatin, CL2 (L/h)") # 2017 Table 1: CL2 = 1,740 L/hr ('Clearance compartment 2')
    lvc <- fixed(log(8980))
    label("Apparent central volume of simvastatin, V2 (L)") # 2017 Table 1: V2 = 8,980 L ('Volume compartment 2')

    # --- Pharmacokinetics: simvastatin acid (active metabolite) ---
    logitfm <- fixed(qlogis(0.30))
    label("Logit of the fraction of simvastatin clearance forming simvastatin acid, FM (unitless)") # 2015 Figure 4: 'CL2 (30%)' into simvastatin acid vs 'CL2 (70%)' to other elimination
    lcle_acid <- fixed(log(383))
    label("Apparent elimination clearance of simvastatin acid, CL3 (L/h)") # 2017 Table 1: CL3 = 383 L/hr ('Clearance compartment 3')
    lvc_acid <- fixed(log(1190))
    label("Apparent central volume of simvastatin acid, V3 (L)") # 2017 Table 1: V3 = 1,190 L ('Volume compartment 3')

    # --- Pharmacodynamics: LDL turnover inhibited by simvastatin acid ---
    lkin <- fixed(log(29.52))
    label("Zero-order production rate of circulating LDL, Kin (nmol/L/h)") # 2017 Table 1: Kin = 29.52 nmol/L/hr
    lrbase <- fixed(log(1400))
    label("Baseline circulating LDL, LDLbaseline (nmol/L)") # 2017 Table 1: LDLbaseline = 1,400 nmol/L
    limax <- fixed(log(0.489))
    label("Maximum fractional inhibition of LDL production, Emax (unitless)") # 2017 Table 1: Emax = 0.489
    lic50 <- fixed(log(0.0868))
    label("Simvastatin acid concentration at half-maximal inhibition, EC50 (ng/mL)") # 2017 Table 1: EC50 = 0.0868 ng/mL

    # --- Residual error ---
    # The source is a simulation model: it reports no residual-error estimate
    # for any of the three outputs, so each is fixed at zero rather than
    # invented.
    propSd <- fixed(0)
    label("Proportional residual SD for simvastatin (fraction; ZERO - not reported in source)")
    propSd_acid <- fixed(0)
    label("Proportional residual SD for simvastatin acid (fraction; ZERO - not reported in source)")
    propSd_ldl <- fixed(0)
    label("Proportional residual SD for circulating LDL (fraction; ZERO - not reported in source)")
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    fm <- expit(logitfm)
    cle_acid <- exp(lcle_acid)
    vc_acid <- exp(lvc_acid)

    kin <- exp(lkin)
    rbase <- exp(lrbase)
    imax <- exp(limax)
    ic50 <- exp(lic50)

    # Kout is not an independent parameter. Pichardo-Almarza 2015 Table 1
    # lists it as 'Kin/LDLbaseline', i.e. it is defined so that the drug-free
    # turnover pool sits exactly at the reported baseline.
    kout <- kin / rbase

    # CL2 is the TOTAL clearance out of the simvastatin compartment; the
    # fraction fm of that flux appears as simvastatin acid and the remaining
    # (1 - fm) leaves the system (2015 Figure 4).
    kel <- cl / vc
    kel_acid <- cle_acid / vc_acid

    # Amounts are in mg and volumes in L, giving mg/L; the factor 1000
    # converts to the ng/mL scale on which EC50 is reported.
    Cc <- 1000 * central / vc
    Cc_acid <- 1000 * central_acid / vc_acid

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    d/dt(central_acid) <- fm * kel * central - kel_acid * central_acid

    ldl(0) <- rbase
    d/dt(ldl) <- kin * (1 - imax * Cc_acid / (ic50 + Cc_acid)) - kout * ldl

    Cc ~ prop(propSd)
    Cc_acid ~ prop(propSd_acid)
    ldl ~ prop(propSd_ldl)
  })
}
