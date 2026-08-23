DeJongh_2025_azd7648_mouse <- function() {
  description <- paste(
    "Preclinical (mouse, female SCID and athymic nude).",
    "Population PK model for the DNA-PK inhibitor AZD7648 after oral",
    "dosing, with first-order absorption, two-compartment disposition",
    "constrained to Vp = Vc, and parallel linear (CL) plus saturable",
    "Michaelis-Menten (Vmax / Km) elimination from the central",
    "compartment. Mouse strain enters as a categorical covariate on two",
    "parameters: absorption rate (estimated in SCID, fixed to a high",
    "value in nude mice for lack of absorption-phase samples) and",
    "relative bioavailability. A single joint inter-individual random",
    "effect acts on both elimination routes together and is correlated",
    "with the random effect on Vc. Amounts are molar (umol/kg) and",
    "concentrations are uM. Parameter values from DeJongh 2025 Table 1."
  )
  reference <- paste(
    "DeJongh J, Cadogan E, Davies M, Ramos-Montoya A, Smith A,",
    "van Steeg T, Richards R. (2025).",
    "Defining preclinical efficacy with the DNAPK inhibitor AZD7648",
    "in combination with olaparib: a minimal systems",
    "pharmacokinetic-pharmacodynamic model.",
    "J Pharmacokinet Pharmacodyn 52:17.",
    "doi:10.1007/s10928-025-09962-x.",
    sep = " "
  )
  vignette <- "DeJongh_2025_azd7648_olaparib_xenograft"

  units <- list(
    time          = "h",
    dosing        = "umol/kg",
    concentration = "uM"
  )

  compartmentData <- list(
    depot       = list(analyte = "AZD7648", units = "umol/kg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "AZD7648", units = "umol/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "AZD7648", units = "umol/kg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    STRAIN_NUDE = list(
      source_name        = "STR",
      description        = paste(
        "Mouse strain indicator: 1 = Hsd:Athymic Nude-Foxn1nu, 0 =",
        "C.B-17/IcrHan(R)Hsd-Prkdcscid (SCID). Gates both the absorption",
        "rate constant and relative bioavailability."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "NONMEM control stream (DeJongh 2025 supplementary file 5,",
        "10928_2025_9962_MOESM5_ESM.txt): `NUDE=0; IF (STR.EQ.1) NUDE=1`",
        "with `F1 = (1 + NUDE*THETA(8))` and strain-specific KA."
      )
    )
  )

  population <- list(
    species        = "mouse (female SCID C.B-17/IcrHan(R)Hsd-Prkdcscid and female Hsd:Athymic Nude-Foxn1nu)",
    n_studies      = 8L,
    age_range      = "at least 6 weeks of age",
    weight_range   = "(not reported in the source publication)",
    sex_female_pct = 100,
    disease_state  = "non-tumour-bearing and FaDu ATM-knockout xenograft-bearing satellite PK groups",
    dose_range     = paste(
      "Oral AZD7648 18.75-150 mg/kg as single and multiple doses,",
      "alone and in combination with oral olaparib. Doses enter the",
      "model in umol/kg (100 mg/kg = 262.8812 umol/kg, i.e. a molar",
      "mass of 380.4 g/mol, read from the AMT and DOSE_AZD columns of",
      "the supplied NONMEM dataset, supplementary file 7)."
    ),
    regions        = "preclinical (AstraZeneca, Cambridge UK)",
    notes          = paste(
      "Plasma concentration data after single or multiple oral dosing of",
      "AZD7648 in SCID and nude mice, pooled from eight studies into one",
      "population PK dataset; the underlying PK study designs are",
      "reported in Fok et al. (DeJongh 2025 reference 3). Blood sampled",
      "by tail-vein venipuncture (20 uL per time point, max 100 uL over",
      "24 h), mixed 1:5 with PBS, and analysed by UPLC-MS/MS against an",
      "11-point calibration curve (1-10,000 nM); LLOQ 0.009 nM.",
      "Estimation was FOCEI in NONMEM 7.4.3 with ADVAN6. Condition",
      "number 243. See DeJongh 2025 Materials and methods and Table 1."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural PK -- DeJongh 2025 Table 1, cross-checked against the
    # theta list carried into the PK-PD control stream (supplementary
    # file 4, 10928_2025_9962_MOESM4_ESM.txt lines 78-92), which reports
    # the same estimates to five significant figures.
    # ------------------------------------------------------------------
    lka_scid <- log(2.7726);   label("Absorption rate constant, SCID mice (ka, 1/h)")            # Table 1 Ka* = 2.77 (RSE 7.83%); control stream THETA(1) = 2.7726
    lka_nude <- fixed(log(9.9)); label("Absorption rate constant, nude mice (ka, 1/h)")          # Table 1: 9.90, no informative absorption-phase samples in nude mice; validated by likelihood profiling
    lcl      <- log(0.25119);  label("Linear clearance (CL, L/kg/h)")                            # Table 1 CL = 0.251 (RSE 7.69%); control stream THETA(3) = 0.25119
    lvmax    <- log(4.2741);   label("Maximum saturable elimination rate (Vmax, umol/kg/h)")     # Table 1 Vmax = 4.27 (RSE 13.1%); control stream THETA(4) = 4.2741
    lkm      <- log(3.7268);   label("Concentration at half-maximal elimination (Km, uM)")       # Table 1 Km = 3.73 (RSE 8.32%); control stream THETA(5) = 3.7268
    lvc      <- log(3.4475);   label("Central volume of distribution (Vc, L/kg)")                # Table 1 Vc = 3.45 (RSE 17.2%); control stream THETA(6) = 3.4475
    lq       <- log(0.93210);  label("Intercompartmental clearance (Q, L/kg/h)")                 # Table 1 Q = 0.932 (RSE 23.5%); control stream THETA(7) = 0.93210

    # Vp is not a separate parameter: DeJongh 2025 Table 1 records
    # "Assumption: Vp = Vp" for the peripheral volume and the Results
    # text states AZD7648 "was described by a reduced two-compartment
    # model (Vc = Vp)". The control stream sets `Vp = Vc` after the Vc
    # random effect is applied, so the IIV on Vc propagates to Vp.

    # ------------------------------------------------------------------
    # Covariate effect -- strain on relative bioavailability.
    # Control stream: F1 = (1 + NUDE*THETA(8)); THETA(8) = -0.51998, so
    # nude mice have F1 = 0.480 relative to SCID mice (F1 = 1).
    # ------------------------------------------------------------------
    e_strain_nude_fdepot <- -0.51998; label("Fractional change in relative bioavailability in nude vs SCID mice (unitless)")  # Table 1 "F1 (relative to SCID mice)" = -0.520 (RSE -35.3%)

    # ------------------------------------------------------------------
    # Inter-individual variability -- DeJongh 2025 Table 1, OMEGA BLOCK(2).
    # etalcl is the control stream's `IIVcl`: ONE joint random effect that
    # multiplies the summed elimination rate (linear CL term plus
    # Michaelis-Menten term), not a random effect on CL alone. The paper
    # states "Random effects for inter-individual variability were
    # identified on Vc and on the sum of both elimination routes. A
    # non-diagonal element was included in the Omega matrix to account for
    # their correlation." Values are variances / covariance on the
    # log scale, as estimated by NONMEM.
    # ------------------------------------------------------------------
    etalcl + etalvc ~ c(0.0631,
                       -0.0291, 0.193)                                                          # Table 1: Om_1 (CL+Vmax) = 0.0631, Om_1x2 = -0.0291, Om_2 (V) = 0.193

    # ------------------------------------------------------------------
    # Residual error. Table 1 reports the proportional residual error as a
    # VARIANCE (0.304); the control stream confirms the scale with
    # `W = IPRED*SQRT(THETA(9))` and `$SIGMA 1 FIX`, so the residual SD is
    # sqrt(0.304) = 0.5514 (55.1% CV).
    # ------------------------------------------------------------------
    propSd <- 0.55136; label("Proportional residual error (fraction)")                           # sqrt(Table 1 proportional residual error variance 0.304)
  })

  model({
    # --- Strain-dependent parameters ------------------------------------
    # Written as a linear blend of the two log-scale strata rather than an
    # if() so the expression stays differentiable for re-fitting.
    ka     <- exp(lka_scid * (1 - STRAIN_NUDE) + lka_nude * STRAIN_NUDE)
    fdepot <- 1 + e_strain_nude_fdepot * STRAIN_NUDE

    # --- Individual structural parameters --------------------------------
    cl   <- exp(lcl)
    vmax <- exp(lvmax)
    km   <- exp(lkm)
    vc   <- exp(lvc + etalvc)
    vp   <- vc
    q    <- exp(lq)

    # Joint random effect on both elimination routes (control stream
    # `IIVcl = EXP(ETA(1))`, applied to the whole elimination expression).
    iivEl <- exp(etalcl)

    Cc <- central / vc

    # --- Disposition (DeJongh 2025 Equations 1-3) ------------------------
    # The Michaelis-Menten term is transcribed from the NONMEM $DES block
    # (supplementary files 4 and 5): `Vmax*A(2)/(A(2)/Vc + Km)`, i.e. the
    # numerator is the central AMOUNT. Printed Equation 1 of the paper
    # omits that factor and reads `Vmax/(A1/Vc + Km)`, which is
    # dimensionally a clearance rather than a rate; see the vignette
    # Errata. The control stream is authoritative.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <- ka * depot -
                         (q / vc) * central + (q / vp) * peripheral1 -
                         ((cl / vc) * central + vmax * central / (Cc + km)) * iivEl
    d/dt(peripheral1) <- (q / vc) * central - (q / vp) * peripheral1

    f(depot) <- fdepot

    Cc ~ prop(propSd)
  })
}
