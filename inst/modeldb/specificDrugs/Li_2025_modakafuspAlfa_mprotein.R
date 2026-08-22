Li_2025_modakafuspAlfa_mprotein <- function() {
  description <- "Sequential population PK-PD model for the serum M-protein response to modakafusp alfa (TAK-573) in adults with relapsed or refractory multiple myeloma (Li 2025). The PK layer is the Michaelis-Menten approximation model with anti-drug-antibody (ADA) binding described in modellib('Li_2025_modakafuspAlfa'); unbound drug concentration drives a Claret tumor-growth-inhibition model in which serum M-protein grows exponentially and is killed by a saturable Emax function whose potency decays exponentially as resistance develops. Baseline serum M-protein depends on baseline serum albumin."
  reference <- paste(
    "Li C, Santulli A, Van Wart S, Yang L, Suryanarayan K, Cook SF, Parot X,",
    "Mager DE, Gupta N. Population Pharmacokinetic and",
    "Pharmacokinetic-Pharmacodynamic Modeling of Serum M-Protein Response for",
    "Modakafusp Alfa in a Phase 1/2 Study of Patients With Relapsed or",
    "Refractory Multiple Myeloma. Clin Transl Sci. 2025;18(7):e70296.",
    "doi:10.1111/cts.70296.",
    "PD parameter values from Table 3 and its footnotes a and b; PK parameter",
    "values from Table 2; ODE structure from Equation 1, Figure 1 and the",
    "PK-PD NONMEM control stream reproduced in the Supporting Information.",
    "Tumor-dynamics form follows Claret L, Girard P, Hoff PM, et al.",
    "J Clin Oncol. 2009;27(25):4103-4108. Companion model:",
    "modellib('Li_2025_modakafuspAlfa').",
    sep = " "
  )
  vignette <- "Li_2025_modakafuspAlfa"
  units <- list(
    time          = "day",
    dosing        = "nmol",
    concentration = "g/L (serum M-protein, the fitted endpoint; the derived unbound modakafusp alfa concentration Cc is in ng/mL)"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight, time-varying. Enters the typical value of the central volume as a power function centred on 80.8 kg.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Carried unchanged from the population PK layer; see modellib('Li_2025_modakafuspAlfa'). Li 2025 Methods used the CURRENT body weight at the start of each treatment cycle rather than the baseline weight, because weight-based doses were re-calculated each cycle.",
      source_name        = "WTKG"
    ),
    ADA_TITER = list(
      description        = "Anti-drug-antibody titer carried on the log3 scale and rounded to the nearest integer, time-varying. Drives both the total ADA target pool and the drug-ADA dissociation rate constant in the PK layer.",
      units              = "log3(reciprocal dilution), rounded to the nearest integer (a value of 4 is a reciprocal titer of 81)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "ZERO-ENCODING: ADA-negative is encoded as 0 on this log3 scale, and the model gates every ADA term on ADA_TITER > 3 (the NONMEM stream's L3OFFTITER = 3.0), so the ADA pool is exactly zero for ADA-negative records and for reciprocal titers at or below 27. Carried unchanged from the population PK layer; see modellib('Li_2025_modakafuspAlfa').",
      source_name        = "L3TITER"
    ),
    ALB = list(
      description        = "Baseline serum albumin. Enters the typical value of the baseline serum M-protein concentration as a power function.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "TIME-FIXED baseline value, not a time-varying laboratory value: the source data column is ALB0 and the effect is on the M-protein baseline, not on a disposition parameter. UNIT CONVERSION: Li 2025 reports albumin in g/dL (mean 3.75 g/dL, SD 0.57, range 1.9-4.6 in the M-protein analysis population; Table 1) and Table 3 footnote a centres the effect on 3.6 g/dL. The nlmixr2lib canonical unit for ALB is g/L, so the reference used here is 36 g/L and any ALB column must be supplied in g/L (multiply a g/dL value by 10). The power ratio is unchanged by the conversion because it is a ratio of two albumin values. Albumin on baseline M-protein was the only covariate retained in the PK-PD model (p < 0.00001).",
      source_name        = "ALB0"
    )
  )

  compartmentData <- list(
    central     = list(analyte = "modakafusp alfa (unbound; the species measured by the ELISA)", units = "nmol", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "modakafusp alfa (unbound)", units = "nmol", specimen = "serum", verified = TRUE),
    complex     = list(analyte = "modakafusp alfa bound to anti-drug antibody", units = "nM", specimen = "serum", verified = TRUE),
    mprotein    = list(analyte = "serum M-protein (monoclonal immunoglobulin, a marker of multiple myeloma tumor burden)", units = "g/L", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 60L,
    n_studies      = 1L,
    age_range      = "34-84 years",
    age_mean       = "63.1 years (SD 10.4)",
    weight_range   = "40-168 kg",
    weight_mean    = "80.3 kg (SD 22.5)",
    sex_female_pct = 42.4,
    race_ethnicity = c(Caucasian = 77.8, Black = 15.2, Asian = 4.0, Other_or_missing = 3.0),
    disease_state  = "Relapsed or refractory multiple myeloma with a baseline serum M-protein concentration of at least 5 g/L, after at least three prior lines of therapy.",
    dose_range     = "0.001 to 6.0 mg/kg intravenously on weekly, every-2-week, every-3-week and every-4-week schedules.",
    regions        = "Multicenter; Phase 1/2 iinnovate-1 trial (NCT03215030), Parts 1 (dose escalation) and 2 (dose expansion), data cutoff 30 May 2022.",
    baseline_albumin = "3.75 g/dL (SD 0.57), range 1.9-4.6, i.e. 37.5 g/L (SD 5.7), range 19-46.",
    n_observations = "492 serum M-protein observations, 86 (17.4%) of which were below the 0.1 g/L SPEP limit of quantification and were handled as censored data.",
    notes          = "The serum M-protein evaluable population is the 60 patients with a baseline serum M-protein of at least 5 g/L. Demographics are quoted from Li 2025 Table 1, 'Serum M-protein analysis population' column; note that the sex, ECOG and race counts in that column sum to 99 rather than 60, so the column appears to describe a wider PK-PD-evaluable set than the 60 patients the Results text reports for the final M-protein fit (see the vignette Errata). Estimation used NONMEM 7.4.4 with the first-order conditional estimation method with interaction and the Laplacian option; individual post hoc PK parameters from the population PK run were fixed as the drug-exposure driver."
  )

  ini({
    # ---- PK layer. Identical to modellib('Li_2025_modakafuspAlfa'); see that
    # file for the full source trace. The source fitted the PD sequentially
    # with each patient's post hoc PK parameters held fixed, so carrying the
    # published typical values and their IIV here is what makes this file
    # simulatable on its own. Time base is DAYS; per-hour rate constants from
    # Li 2025 Table 2 are multiplied by 24. Table 2 heads the CL and Q rows
    # "(L/day)", but both values are per HOUR; see the annotated derivation in
    # modellib('Li_2025_modakafuspAlfa') and the vignette Errata.
    lcl     <- log(0.0553 * 24); label("Linear clearance of unbound drug from the central compartment (L/day)")      # Li 2025 Table 2: CL = 0.0553 (%RSE 8.97), tabulated as "L/day" but in fact L/h; x24 for the day time base
    lvc     <- log(5.08)     ; label("Central volume of distribution for a 80.8 kg patient (L)")                     # Li 2025 Table 2: Vc coefficient for a typical 80.8 kg patient = 5.08 L (%RSE 6.44)
    e_wt_vc <- 0.509         ; label("Power exponent on (WT/80.8) for the central volume (unitless)")                # Li 2025 Table 2: Vc-weight exponent = 0.509 (%RSE 30.4); Table 2 footnote b
    lq      <- log(0.137 * 24); label("Intercompartmental clearance (L/day)")                                        # Li 2025 Table 2: Q = 0.137 (%RSE 16.4), tabulated as "L/day" but in fact L/h; x24 for the day time base
    lvp     <- log(4.01)     ; label("Peripheral volume of distribution (L)")                                        # Li 2025 Table 2: Vp = 4.01 L (%RSE 16.7)
    lvmax   <- log(4.77 * 24); label("Maximum rate of Michaelis-Menten elimination (nmol/day)")                      # Li 2025 Table 2: Vmax = 4.77 nmol/h (%RSE 13.6); x24 for the day time base
    lkm     <- log(1.26)     ; label("Unbound drug concentration producing half of Vmax (nM)")                       # Li 2025 Table 2: KM = 1.26 nM (%RSE 40.0). The PK-PD control stream hardcodes LTVKM = 0.229, i.e. exp(0.229) = 1.257 nM, confirming the Table 2 value.

    lrtot_ada            <- log(1.35)  ; label("Total ADA target pool at a log3 ADA titer of 4, i.e. a reciprocal titer of 81 (nM)")  # Li 2025 Table 2: Rtot,ADA coefficient for ADA titer of 81 = 1.35 nM (%RSE 118)
    e_ada_titer_rtot_ada <- 0.412      ; label("Power exponent on (titer/81) for the total ADA pool (unitless)")                      # Li 2025 Table 2: Rtot,ADA-Titer power = 0.412 (%RSE 4.29), also hardcoded as 0.412 in the PK-PD control stream. Table 2 footnote c prints 4.12, a decimal-point error falsified by Figure S5.
    lkon_ada             <- fixed(log(3.6 * 24))    ; label("ADA association rate constant (1/(nM*day))")                            # Li 2025 Table 2: kon,ADA = 3.6 1/(nM*h), literature value from Chen 2016 (reference 14); x24
    lkoff_ada_max        <- fixed(log(3600 * 24))   ; label("ADA dissociation rate constant at the reference log3 titer of 4 (1/day)")  # Li 2025 Table 2: koff,ADA,max = 3600 1/h, literature value from Chen 2016 (reference 14); x24
    lkoff_ada_min        <- fixed(log(0.036 * 24))  ; label("ADA dissociation rate constant as the ADA titer approaches infinity (1/day)")  # Li 2025 Table 2: koff,ADA,min = 0.036 1/h, literature value from Chen 2016 (reference 14); x24
    kdecay_ada           <- 0.0275     ; label("Rate of decline of koff,ADA per unit log3 ADA titer (1/log3 titer unit)")            # Li 2025 Table 2: kdec = 0.0275 (%RSE 5.23). The PK-PD control stream instead hardcodes KD_TITERDEC = 0.210, which is kel,ADA's value and is contradicted by Figure S6; see the vignette Errata.
    lkel_ada             <- log(0.210 * 24)         ; label("Elimination rate constant of the modakafusp alfa-ADA complex (1/day)")  # Li 2025 Table 2: kel,ADA = 0.210 1/h (%RSE 10.8); x24

    # ---- Serum M-protein (Claret tumor-growth-inhibition) layer,
    # Li 2025 Equation 1 and Table 3. Rate constants are reported per week
    # and are divided by 7 for the day time base, mirroring the source
    # stream's own "/24/7" conversion from weeks to its hourly time base.
    lrbase      <- log(17.4)        ; label("Baseline serum M-protein at a serum albumin of 36 g/L (g/L)")                    # Li 2025 Table 3: MP0 coefficient = 17.4 g/L (%RSE 4.40); Table 3 footnote a
    e_alb_rbase <- -2.74            ; label("Power exponent on (ALB/36) for the baseline serum M-protein (unitless)")         # Li 2025 Table 3: ALB0 effect coefficient on MP0 = -2.74 (%RSE 28.7); Table 3 footnote a: TVMP0 = 17.4 * (ALB/3.6)^-2.74 with albumin in g/dL, equivalently (ALB/36) with albumin in g/L
    lp          <- log(0.0303 / 7)  ; label("First-order serum M-protein growth rate constant (1/day)")                       # Li 2025 Table 3: kg = 0.0303 1/week (%RSE 18.5); /7 for the day time base
    llambda     <- log(0.0568 / 7)  ; label("First-order rate constant for the appearance of resistance to treatment (1/day)")# Li 2025 Table 3: kr = 0.0568 1/week (%RSE 61.2); /7 for the day time base
    lkkillmax   <- log(0.816 / 7)   ; label("Maximal serum M-protein kill rate constant (1/day)")                             # Li 2025 Table 3: kkill,max = 0.816 1/week (%RSE 42.6); /7 for the day time base
    lec50       <- log(447)         ; label("Unbound modakafusp alfa concentration producing half of kkill,max (ng/mL)")      # Li 2025 Table 3: EC50 = 447 ng/mL (%RSE 48.0), on the unbound-drug driver

    # ---- PK-layer inter-individual variability (Li 2025 Table 2), carried
    # as independent etas because only the diagonal of the source BLOCK(5)
    # OMEGA is published. See the vignette Errata.
    etalcl       ~ 0.307  ; label("IIV variance on log CL")                # Li 2025 Table 2: omega^2 for CL = 0.307 (60.0% CV)
    etalvc       ~ 0.391  ; label("IIV variance on log Vc")                # Li 2025 Table 2: omega^2 for Vc = 0.391 (69.2% CV)
    etalq        ~ 2.84   ; label("IIV variance on log Q")                 # Li 2025 Table 2: omega^2 for Q = 2.84 (402% CV)
    etalvp       ~ 1.04   ; label("IIV variance on log Vp")                # Li 2025 Table 2: omega^2 for Vp = 1.04 (135% CV)
    etalvmax     ~ 1.14   ; label("IIV variance on log Vmax")              # Li 2025 Table 2: omega^2 for Vmax = 1.14 (146% CV)
    etalrtot_ada ~ 1.63   ; label("IIV variance on log Rtot,ADA")          # Li 2025 Table 2: omega^2 for Rtot,ADA = 1.63 (203% CV)

    # ---- PD-layer inter-individual variability (Li 2025 Table 3). EC50 and
    # kkill,max were estimated as a 2x2 block; Table 3 footnote b gives the
    # covariance, which reproduces the printed correlation coefficient:
    # -0.374 / sqrt(1.57 * 1.11) = -0.283.
    etalrbase    ~ 0.343  ; label("IIV variance on log MP0")               # Li 2025 Table 3: omega^2 for MP0 = 0.343 (63.9% CV, %RSE 28.7)
    etalp        ~ 1.01   ; label("IIV variance on log kg")                # Li 2025 Table 3: omega^2 for kg = 1.01 (132% CV, %RSE 29.9)
    etallambda   ~ 1.25   ; label("IIV variance on log kr")                # Li 2025 Table 3: omega^2 for kr = 1.25 (157% CV, %RSE > 100)
    etalec50 + etalkkillmax ~ c(1.57,
                                -0.374, 1.11)                              # Li 2025 Table 3: omega^2 for EC50 = 1.57 (195% CV) and for kkill,max = 1.11 (143% CV); footnote b covariance = -0.374 (174% RSE, correlation -0.282)

    # ---- Residual error on serum M-protein. The source $ERROR is
    # Y = IPRED + IPRED*EPS(1) + EPS(2), a combined proportional-plus-additive
    # model, and Table 3 reports both components as variances.
    propSd <- 0.1268858  ; label("Proportional residual error on serum M-protein (fraction)")  # Li 2025 Table 3: proportional sigma^2 = 0.0161 (%RSE 5.63); SD = sqrt(0.0161) = 0.126886
    addSd  <- 2.5651511  ; label("Additive residual error on serum M-protein (g/L)")           # Li 2025 Table 3: additive sigma^2 = 6.58 (%RSE 13.3); SD = sqrt(6.58) = 2.565151
  })

  model({
    # ---- Molecular weight conversion, 186 kDa: 1 nM = 186 ng/mL. The source
    # $DES converts the driver with "CP = C*186".
    ngmlPerNm <- 186

    # ---- PK layer (see modellib('Li_2025_modakafuspAlfa') for the annotated
    # version of these lines).
    cl      <- exp(lcl + etalcl)
    vc      <- exp(lvc + etalvc) * (WT / 80.8)^e_wt_vc
    q       <- exp(lq + etalq)
    vp      <- exp(lvp + etalvp)
    vmax    <- exp(lvmax + etalvmax)
    km      <- exp(lkm)
    kel     <- cl / vc
    k12     <- q / vc
    k21     <- q / vp

    adaOn      <- (ADA_TITER > 3)
    kon_ada    <- exp(lkon_ada)
    kel_ada    <- exp(lkel_ada)
    koff_max   <- exp(lkoff_ada_max)
    koff_min   <- exp(lkoff_ada_min)
    rtot_ada   <- adaOn * exp(lrtot_ada + etalrtot_ada) * (3^ADA_TITER / 81)^e_ada_titer_rtot_ada
    koff_ada   <- koff_min * exp(log(koff_max / koff_min) * exp(-kdecay_ada * (ADA_TITER - 4)))

    cp         <- central / vc
    mmElim     <- vmax * cp / (km + cp)
    adaFree    <- adaOn * (rtot_ada - complex)

    d/dt(central)     <- k21 * peripheral1 - (kel + k12) * central - mmElim -
      kon_ada * central * adaFree + koff_ada * complex * vc
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(complex)     <- kon_ada * cp * adaFree - koff_ada * complex - kel_ada * complex

    # Unbound modakafusp alfa on the assay scale; this is the PD driver.
    # Li 2025 Results preferred the unbound driver over total drug and over a
    # neutralising-antibody 'switch' model because all three fitted similarly.
    Cc <- cp * ngmlPerNm

    # ---- Serum M-protein layer, Li 2025 Equation 1:
    #   dMP/dt = kg*MP - kkill,max * exp(-kr*t) * (Cp/(EC50 + Cp)) * MP
    # The exponential term is the Claret time-dependent loss of drug effect as
    # resistance develops, so t is time since the start of treatment.
    mp0      <- exp(lrbase + etalrbase) * (ALB / 36)^e_alb_rbase
    kg       <- exp(lp + etalp)
    kr       <- exp(llambda + etallambda)
    kkillmax <- exp(lkkillmax + etalkkillmax)
    ec50     <- exp(lec50 + etalec50)

    mprotein(0) <- mp0
    d/dt(mprotein) <- kg * mprotein -
      kkillmax * exp(-kr * t) * (Cc / (ec50 + Cc)) * mprotein

    mprotein ~ add(addSd) + prop(propSd)
  })
}
