Li_2025_modakafuspAlfa <- function() {
  description <- "Two-compartment population PK model for modakafusp alfa (TAK-573), an anti-CD38 IgG4 antibody fused to attenuated interferon alfa-2b, in adults with relapsed or refractory multiple myeloma (Li 2025). Unbound drug leaves the central compartment by parallel linear and Michaelis-Menten pathways, the latter a Michaelis-Menten approximation of target-mediated disposition. Independently of that, unbound drug binds reversibly to an anti-drug-antibody (ADA) pool and the resulting drug-ADA complex is eliminated. Both the total ADA pool and the drug-ADA dissociation rate constant are driven by the observed log3-transformed ADA titer, so ADA-mediated elimination grows as titers rise. Time-varying body weight acts on central volume only."
  reference <- paste(
    "Li C, Santulli A, Van Wart S, Yang L, Suryanarayan K, Cook SF, Parot X,",
    "Mager DE, Gupta N. Population Pharmacokinetic and",
    "Pharmacokinetic-Pharmacodynamic Modeling of Serum M-Protein Response for",
    "Modakafusp Alfa in a Phase 1/2 Study of Patients With Relapsed or",
    "Refractory Multiple Myeloma. Clin Transl Sci. 2025;18(7):e70296.",
    "doi:10.1111/cts.70296.",
    "Parameter values from Table 2 and its footnotes b, c and d; ODE structure",
    "from Figure 1 and from the NONMEM control stream reproduced in the",
    "Supporting Information. Companion model:",
    "modellib('Li_2025_modakafuspAlfa_mprotein').",
    sep = " "
  )
  vignette <- "Li_2025_modakafuspAlfa"
  units <- list(time = "day", dosing = "nmol", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight, time-varying. Enters the typical value of the central volume as a power function centred on 80.8 kg.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Li 2025 Methods used the CURRENT body weight at the start of each treatment cycle rather than the baseline weight, because weight-based doses were re-calculated each cycle; over the study individual weight changed by -14.1 kg (-14.3%) to +16.0 kg (+21.7%). The reference weight of 80.8 kg comes from a prior interim analysis, not from the final analysis population (whose mean weight was 79.9 kg). Body weight was tested on CL, Vc, Q, Vp and Vmax with both estimated and theoretically fixed exponents (0.75 for clearances, 1 for volumes) and was statistically significant only on Vc.",
      source_name        = "WTKG"
    ),
    ADA_TITER = list(
      description        = "Anti-drug-antibody titer carried on the log3 scale and rounded to the nearest integer, time-varying. Drives both the total ADA target pool and the drug-ADA dissociation rate constant.",
      units              = "log3(reciprocal dilution), rounded to the nearest integer (a value of 4 is a reciprocal titer of 81)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "ZERO-ENCODING: ADA-negative is encoded as 0 on this log3 scale, and the model gates every ADA term on ADA_TITER > 3 (the NONMEM stream's L3OFFTITER = 3.0), so the ADA pool is exactly zero for ADA-negative records and for reciprocal titers at or below 27. Li 2025 Methods: 'ADA titer values were transformed to log 3 scale and rounded to the nearest integer for use as a time-varying covariate.' Observed positive values ran from 4 to 13 (Figure S5); the maximum-titer quartiles of ADA-positive patients were 75 to 675, 675 to 54700, 54700 to 164000 and above 164000 on the reciprocal scale (Figure S4 legend). Titers were measured pre-dose within each cycle and carried forward between measurements. In the PK-evaluable population 12.5% were ADA-positive at baseline and 52.1% became ADA-positive on treatment.",
      source_name        = "L3TITER"
    )
  )

  compartmentData <- list(
    central     = list(analyte = "modakafusp alfa (unbound; the species measured by the ELISA)", units = "nmol", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "modakafusp alfa (unbound)", units = "nmol", specimen = "serum", verified = TRUE),
    complex     = list(analyte = "modakafusp alfa bound to anti-drug antibody", units = "nM", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 96L,
    n_studies      = 1L,
    age_range      = "34-84 years",
    age_mean       = "63.1 years (SD 10.5)",
    weight_range   = "60-169 kg",
    weight_mean    = "79.9 kg (SD 23.0)",
    sex_female_pct = 42.7,
    race_ethnicity = c(Caucasian = 78.1, Black = 14.6, Asian = 4.2, Other_or_missing = 3.1),
    disease_state  = "Relapsed or refractory multiple myeloma with disease progression after at least three prior lines of therapy; ECOG performance status 0 (14.6%), 1 (82.3%) or 2 (3.1%).",
    dose_range     = "0.001 to 6.0 mg/kg intravenously, initially over a 4-h ramped infusion and after a protocol amendment over 1 h (or more than 2 h at the highest dose level), on weekly, every-2-week, every-3-week and every-4-week schedules.",
    regions        = "Multicenter; Phase 1/2 iinnovate-1 trial (NCT03215030), Parts 1 (dose escalation) and 2 (dose expansion), data cutoff 30 May 2022.",
    renal_function = "Creatinine clearance mean 59.8 mL/min/1.73 m^2 (SD 19.4), range 24.1-101.",
    n_observations = "2297 quantifiable serum modakafusp alfa concentrations plus 806 records below the 6.25 ng/mL LLOQ, which were handled as censored data by the Beal M3 likelihood method.",
    notes          = "Demographics from Li 2025 Table 1, PK analysis population column. 100 patients were enrolled (56 in Part 1, 44 in Part 2); the 96 with at least one quantifiable concentration form the PK analysis population. Estimation used NONMEM 7.4.4 with Monte Carlo importance sampling and Laplacian conditional estimation with eta-epsilon interaction."
  )

  ini({
    # ---- Structural PK disposition (Li 2025 Table 2).
    # The model time base here is DAYS. Li 2025 reports CL and Q in L/day and
    # everything else on a per-hour basis; the source NONMEM stream works
    # entirely in hours (its $PK carries the comment "dose in nmol, conc in
    # nmol/L, time in hr"). Per-hour rate constants are therefore multiplied
    # by 24 below and the L/day clearances are used exactly as printed.
    # Amounts are nmol and concentrations nM, as in the source stream.
    lcl     <- log(0.0553)   ; label("Linear clearance of unbound drug from the central compartment (L/day)")        # Li 2025 Table 2: CL = 0.0553 L/day (%RSE 8.97)
    lvc     <- log(5.08)     ; label("Central volume of distribution for a 80.8 kg patient (L)")                     # Li 2025 Table 2: Vc coefficient for a typical 80.8 kg patient = 5.08 L (%RSE 6.44)
    e_wt_vc <- 0.509         ; label("Power exponent on (WT/80.8) for the central volume (unitless)")                # Li 2025 Table 2: Vc-weight exponent = 0.509 (%RSE 30.4); Table 2 footnote b: TVVc (L) = 5.08 * (Weight/80.8)^0.509
    lq      <- log(0.137)    ; label("Intercompartmental clearance (L/day)")                                         # Li 2025 Table 2: Q = 0.137 L/day (%RSE 16.4)
    lvp     <- log(4.01)     ; label("Peripheral volume of distribution (L)")                                        # Li 2025 Table 2: Vp = 4.01 L (%RSE 16.7)

    # ---- Saturable (Michaelis-Menten approximation) elimination pathway.
    lvmax   <- log(4.77 * 24); label("Maximum rate of Michaelis-Menten elimination (nmol/day)")                      # Li 2025 Table 2: Vmax = 4.77 nmol/h (%RSE 13.6); x24 for the day time base
    lkm     <- log(1.26)     ; label("Unbound drug concentration producing half of Vmax (nM)")                       # Li 2025 Table 2: KM = 1.26 nM (%RSE 40.0)

    # ---- Anti-drug-antibody binding (a TMDD binding model applied to ADA,
    # independent of the target-mediated disposition captured by the
    # Michaelis-Menten arm). Li 2025 Table 2 and its footnotes c and d.
    lrtot_ada            <- log(1.35)  ; label("Total ADA target pool at a log3 ADA titer of 4, i.e. a reciprocal titer of 81 (nM)")  # Li 2025 Table 2: Rtot,ADA coefficient for ADA titer of 81 = 1.35 nM (%RSE 118)
    e_ada_titer_rtot_ada <- 0.412      ; label("Power exponent on (titer/81) for the total ADA pool (unitless)")                      # Li 2025 Table 2: Rtot,ADA-Titer power = 0.412 (%RSE 4.29). Table 2 footnote c prints the exponent as 4.12, which is a decimal-point error: the typical-value line of Figure S5 runs from 1.35 nM at log3 titer 4 to about 79 nM at log3 titer 13, which is 1.35*(3^9)^0.412 = 79.3 and rules 4.12 out by many orders of magnitude. The source NONMEM stream also hardcodes 0.412.
    lkon_ada             <- fixed(log(3.6 * 24))    ; label("ADA association rate constant (1/(nM*day))")                            # Li 2025 Table 2: kon,ADA = 3.6 1/(nM*h), literature value from Chen 2016 (reference 14); x24 for the day time base
    lkoff_ada_max        <- fixed(log(3600 * 24))   ; label("ADA dissociation rate constant at the reference log3 titer of 4 (1/day)")  # Li 2025 Table 2: koff,ADA,max = 3600 1/h, literature value from Chen 2016 (reference 14); x24 for the day time base
    lkoff_ada_min        <- fixed(log(0.036 * 24))  ; label("ADA dissociation rate constant as the ADA titer approaches infinity (1/day)")  # Li 2025 Table 2: koff,ADA,min = 0.036 1/h, literature value from Chen 2016 (reference 14); x24 for the day time base
    kdecay_ada           <- 0.0275     ; label("Rate of decline of koff,ADA per unit log3 ADA titer (1/log3 titer unit)")            # Li 2025 Table 2: kdec = 0.0275 (%RSE 5.23). Dimensionless in time: it acts on the log3 titer covariate, not on time.
    lkel_ada             <- log(0.210 * 24)         ; label("Elimination rate constant of the modakafusp alfa-ADA complex (1/day)")  # Li 2025 Table 2: kel,ADA = 0.210 1/h (%RSE 10.8); x24 for the day time base

    # ---- Inter-individual variability. Li 2025 Table 2 reports the omega^2
    # values below together with the coefficient of variation implied by the
    # log-normal relationship CV = sqrt(exp(omega^2) - 1); e.g. CL
    # sqrt(exp(0.307) - 1) = 60.0%, matching the printed "60.0% CV". The
    # source $OMEGA was a BLOCK(5) over CL, Vc, Q, Vp and Vmax plus a
    # separate BLOCK(1) for Rtot,ADA, but only the diagonal variances are
    # published, so the five correlated etas are carried here as independent
    # etas. See the vignette Errata.
    etalcl       ~ 0.307  ; label("IIV variance on log CL")                # Li 2025 Table 2: omega^2 for CL = 0.307 (60.0% CV, %RSE 85.1)
    etalvc       ~ 0.391  ; label("IIV variance on log Vc")                # Li 2025 Table 2: omega^2 for Vc = 0.391 (69.2% CV, %RSE 28.7)
    etalq        ~ 2.84   ; label("IIV variance on log Q")                 # Li 2025 Table 2: omega^2 for Q = 2.84 (402% CV, %RSE 33.0)
    etalvp       ~ 1.04   ; label("IIV variance on log Vp")                # Li 2025 Table 2: omega^2 for Vp = 1.04 (135% CV, %RSE 36.9)
    etalvmax     ~ 1.14   ; label("IIV variance on log Vmax")              # Li 2025 Table 2: omega^2 for Vmax = 1.14 (146% CV, %RSE 36.0)
    etalrtot_ada ~ 1.63   ; label("IIV variance on log Rtot,ADA")          # Li 2025 Table 2: omega^2 for Rtot,ADA = 1.63 (203% CV, %RSE 50.4)

    # ---- Residual error. The source $ERROR is
    # Y = IPRED + IPRED*EPS(1) + EPS(2) with the additive $SIGMA slot fixed
    # to zero, so only the proportional component is active.
    propSd <- 0.3949684   ; label("Proportional residual error on serum modakafusp alfa (fraction)")  # Li 2025 Table 2: residual variability sigma^2 = 0.156 (39.5% CV); SD = sqrt(0.156) = 0.394968
  })

  model({
    # ---- Molecular weight conversion. The source stream carries drug as
    # nmol / nM throughout and converts to the assay scale in $ERROR with
    # "IPRED = (A(1)/S1)*186.0", i.e. a molecular weight of 186 kDa, for
    # which 1 nM = 186 ng/mL.
    ngmlPerNm <- 186

    # ---- Individual disposition parameters.
    cl      <- exp(lcl + etalcl)
    vc      <- exp(lvc + etalvc) * (WT / 80.8)^e_wt_vc
    q       <- exp(lq + etalq)
    vp      <- exp(lvp + etalvp)
    vmax    <- exp(lvmax + etalvmax)
    km      <- exp(lkm)
    kel     <- cl / vc
    k12     <- q / vc
    k21     <- q / vp

    # ---- ADA binding parameters, all driven by the time-varying log3 titer.
    # The source gates every ADA term on L3TITER > L3OFFTITER = 3.0, so a
    # log3 titer of 3 or below (a reciprocal titer of 27 or below, and in
    # particular the ADA-negative encoding of 0) contributes no ADA pool.
    adaOn      <- (ADA_TITER > 3)
    kon_ada    <- exp(lkon_ada)
    kel_ada    <- exp(lkel_ada)
    koff_max   <- exp(lkoff_ada_max)
    koff_min   <- exp(lkoff_ada_min)

    # Total ADA pool, Table 2 footnote c:
    #   TVRtot,ADA (nM) = 1.35 * (Titer/81)^0.412
    # written here on the log3 titer that the data actually carry, so that
    # (3^ADA_TITER)/81 is the reciprocal-titer ratio the footnote uses.
    rtot_ada   <- adaOn * exp(lrtot_ada + etalrtot_ada) * (3^ADA_TITER / 81)^e_ada_titer_rtot_ada

    # Titer-dependent dissociation rate constant, Table 2 footnote d:
    #   koff,ADA (1/h) = 3.6 * exp(ln(0.01) - ln(0.01/1000) * exp(-kdec*(log3(Titer)-4)))
    # With kon,ADA = 3.6 1/(nM*h) the two constants inside the logarithms are
    # the minimum and maximum ADA dissociation equilibrium constants, KD,min =
    # 0.01 nM and KD,max = 1000 nM, so the expression is exactly the
    # interpolation between koff,ADA,max = 3600 1/h at a log3 titer of 4 and
    # koff,ADA,min = 0.036 1/h as the titer grows without bound. Rewritten
    # below in terms of the two published koff values, which is algebraically
    # identical because koff,max/koff,min = 3600/0.036 = KD,max/KD,min.
    # Checked against Figure S6, which plots 3.6e3 1/h at log3 titer 4 and
    # decays to the 3.6e-2 1/h asymptote.
    koff_ada   <- koff_min * exp(log(koff_max / koff_min) * exp(-kdecay_ada * (ADA_TITER - 4)))

    # ---- Unbound drug concentration and the saturable elimination rate.
    cp         <- central / vc
    mmElim     <- vmax * cp / (km + cp)

    # Free (unoccupied) ADA is the total pool less the complex; both are
    # concentrations in nM. Li 2025 Discussion: "Rtot,ADA was assumed to be
    # constant at a given time-varying titer level (i.e., no turnover of
    # ADA)", so the pool is replenished as complex is cleared.
    adaFree    <- adaOn * (rtot_ada - complex)

    # ---- ODE system. central and peripheral1 hold AMOUNTS in nmol; complex
    # holds a CONCENTRATION in nM, exactly as A(3) does in the source $DES
    # (its DADT(3) is written in concentration units and the return flux into
    # the central compartment is multiplied back up by Vc).
    d/dt(central)     <- k21 * peripheral1 - (kel + k12) * central - mmElim -
      kon_ada * central * adaFree + koff_ada * complex * vc
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(complex)     <- kon_ada * cp * adaFree - koff_ada * complex - kel_ada * complex

    # ---- Observation: unbound modakafusp alfa, the species the ELISA
    # measures, reported on the assay scale.
    Cc <- cp * ngmlPerNm
    Cc ~ prop(propSd)
  })
}
