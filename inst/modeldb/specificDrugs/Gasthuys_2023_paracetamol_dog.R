Gasthuys_2023_paracetamol_dog <- function() {
  description <- "Preclinical (beagle dog). Two-compartment population PK model for paracetamol (acetaminophen) in healthy adult beagle dogs given single oral and intravenous doses under six prandial conditions (Gasthuys 2023). Absorption is sequential zero-order then first-order: an oral dose enters the depot compartment as a zero-order input of duration D1, scaled by the absolute oral bioavailability F, and leaves the depot into the central compartment by a first-order rate Ka; an intravenous dose enters the central compartment directly, which is what makes F identifiable. Body weight acts on CL, Vd, Vp and Q as a power function with exponents fixed at the values in Gasthuys 2023 Table 3. Fed state and caloric content were screened as covariates but no food effect was retained, so the model returns a single typical profile for all six dosing conditions; the treatment-related differences the authors could still discern by eye were absorbed into inter-occasion variability on F, Ka and D1."
  reference <- "Gasthuys E, Sandra L, Statelova M, Vertzoni M, Vermeulen A. The Use of Population Pharmacokinetics to Extrapolate Food Effects from Human Adults and Beagle Dogs to the Pediatric Population Illustrated with Paracetamol as a Test Case. Pharmaceuticals. 2024;17(1):53. doi:10.3390/ph17010053 (published online 2023-12-28). Beagle dog study design and data originally reported in Statelova et al. 2023 (cited as ref [14])."
  vignette <- "Gasthuys_2023_paracetamol"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Subject-level, time-fixed. Applied as a power function with exponents fixed at",
        "0.75 for CL and 1 for Vd, Vp and Q (Gasthuys 2023 Table 3, 'Covariate model' rows).",
        "The reference weight is 9.70 kg, the median body weight of the six beagle dogs",
        "(Gasthuys 2023 Methods 4.1.1), NOT the 70 kg written into the paper's general",
        "allometric equation (Eq. 1). Eq. 1 is a generic statement of the covariate form",
        "that is shared with the companion human-adult model; taking it literally for the dogs",
        "is falsified by the paper's own data. See vignette 'Assumptions and deviations'."
      ),
      source_name = "WT"
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "paracetamol", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "paracetamol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "paracetamol", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "beagle dog",
    n_subjects     = 6L,
    n_studies      = 1L,
    n_observations = "575 plasma paracetamol concentrations from six beagle dogs across eight occasions (Gasthuys 2023 Results 2.1). Concentrations below the 7.5 ng/mL limit of quantification (1.19% of dog records) were excluded before fitting.",
    age_range      = "1.7-4.0 years (median 1.9 years)",
    age_median     = "1.9 years",
    weight_range   = "7.90-13.3 kg",
    weight_median  = "9.70 kg",
    disease_state  = "clinically healthy adult beagle dogs (Marshall BioResources, Lyon, France)",
    dose_range     = "168 mg paracetamol per occasion: 7 mL of Panadol suspension by oral gavage followed by 10 mL tap water, and a single 16.8 mL / 168 mg intravenous injection on day 14",
    regions        = "Belgium (Janssen R&D animal facility; ethical committee approval 512)",
    notes          = paste(
      "Sequential (not crossover) design over 50 days, one prandial condition per occasion",
      "(Gasthuys 2023 Figure 5): day 0 fasted, day 14 intravenous, day 22 reference meal 200",
      "(100 g, 200 kcal), day 28 infant formula 100 (150 mL, 100 kcal), day 35 reference meal 100",
      "(50 g, 100 kcal) and infant formula 200 (300 mL, 200 kcal), day 50 fasted after 0.1 M",
      "HCl/KCl (pH 1.6) pretreatment. Plasma assayed by validated HPLC-UV, LLOQ 7.5 ng/mL.",
      "Adult dogs were used deliberately; the authors note in the Discussion that physiological",
      "maturation is still ongoing below 2 years of age and that juvenile dogs might better mimic",
      "the paediatric population this model was built to extrapolate to."
    )
  )

  ini({
    # ======================================================================
    # Structural model -- Gasthuys 2023 Table 3, 'Estimate Beagle Dogs'
    # column. All structural values are estimated point estimates reported
    # with an RSE%; none carries a FIX flag.
    #
    # Disposition parameters are the typical values at the 9.70 kg
    # reference weight (see covariateData$WT notes and the vignette).
    # ======================================================================
    logitfdepot <- qlogis(0.80); label("Absolute oral bioavailability F (unitless)")                                # Gasthuys 2023 Table 3: F = 0.80, RSE 2.56% -- held on the logit scale so the simulated F cannot leak above 1
    lka <- log(2.86);            label("First-order absorption rate constant Ka from depot to central (1/h)")       # Gasthuys 2023 Table 3: ka = 2.86 1/h, RSE 14.2%
    lvc <- log(9.53);            label("Central volume of distribution Vd (L)")                                     # Gasthuys 2023 Table 3: Vd = 9.53 L, RSE 1.55%
    lcl <- log(9.29);            label("Total body clearance CL (L/h)")                                             # Gasthuys 2023 Table 3: CL = 9.29 L/h, RSE 1.92%
    lvp <- log(34.9);            label("Peripheral volume of distribution Vp (L)")                                  # Gasthuys 2023 Table 3: Vp = 34.9 L, RSE 4.56%
    lq  <- log(2.82);            label("Inter-compartmental flow Q (L/h)")                                          # Gasthuys 2023 Table 3: Q = 2.82 L/h, RSE 0.835%
    ld1 <- log(0.64);            label("Duration of the zero-order input into the depot compartment dT1 (h)")       # Gasthuys 2023 Table 3: dT1 = 0.64 h, RSE 16.1%

    # ======================================================================
    # Covariate model -- Gasthuys 2023 Table 3, 'Covariate model' rows.
    # All four exponents are reported as [FIX].
    #
    # Note that Table 3 gives beta_WTonQ = 1 while Methods 4.2 states the
    # exponent was fixed to 0.75 "for the clearance and intercompartmental
    # flow (CL, Q)". Table 3 is the record of the final fit and is used
    # here; the disagreement moves Q by at most 8% across the observed
    # 7.90-13.3 kg weight range. Documented in vignette Errata.
    # ======================================================================
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL (unitless)")  # Gasthuys 2023 Table 3: beta_WTonCL = 0.75 [FIX]
    e_wt_vc <- fixed(1);    label("Allometric exponent on Vd (unitless)")  # Gasthuys 2023 Table 3: beta_WTonVd = 1 [FIX]
    e_wt_vp <- fixed(1);    label("Allometric exponent on Vp (unitless)")  # Gasthuys 2023 Table 3: beta_WTonVp = 1 [FIX]
    e_wt_q  <- fixed(1);    label("Allometric exponent on Q (unitless)")   # Gasthuys 2023 Table 3: beta_WTonQ  = 1 [FIX]

    # ======================================================================
    # Random effects -- Gasthuys 2023 Table 3.
    #
    # The dog column reports NO between-subject variability at all (both
    # 'BSV Vd' and 'BSV CL' rows are '-'); the only random effects are
    # inter-occasion variability on F, ka and dT1. nlmixr2lib has no
    # idiomatic encoding for inter-occasion variability separate from
    # between-subject variability, so per the convention already used in
    # Bienczak_2016_efavirenz.R, Bienczak_2016_nevirapine.R and
    # Svensson_2018_bedaquiline.R, IOV reported on a parameter that carries
    # no BSV term is folded in as a BSV-equivalent. Here that is all three.
    #
    # Monolix (version 2019R2, per Gasthuys 2023 Methods 4.5) reports omega
    # as the standard deviation of the random effect, so the tabulated
    # values are squared to give the log-scale (or logit-scale) variances
    # nlmixr2 expects.
    # ======================================================================
    etalogitfdepot ~ 0.11^2  # Gasthuys 2023 Table 3: IOV F = 0.11 (RSE 14.3%), folded as a BSV-equivalent; variance = 0.11^2 = 0.0121 on the logit scale
    etalka         ~ 0.63^2  # Gasthuys 2023 Table 3: IOV ka = 0.63 (RSE 18.0%), folded as a BSV-equivalent; variance = 0.63^2 = 0.3969
    etald1         ~ 0.86^2  # Gasthuys 2023 Table 3: IOV dT1 = 0.86 (RSE 15.2%), folded as a BSV-equivalent; variance = 0.86^2 = 0.7396

    # ======================================================================
    # Residual error -- Gasthuys 2023 Table 3, 'Error model parameters'
    # rows, footnote *1 "combined proportional and additive error model"
    # with a = additive term and b = proportional term.
    # ======================================================================
    addSd  <- 0.055; label("Additive residual error on plasma paracetamol (ug/mL)")  # Gasthuys 2023 Table 3: a = 0.055, RSE 9.30%
    propSd <- 0.13;  label("Proportional residual error on plasma paracetamol (fraction)")  # Gasthuys 2023 Table 3: b = 0.13, RSE 7.13%
  })

  model({
    # Allometric reference body weight: the median of the six beagle dogs
    # (Gasthuys 2023 Methods 4.1.1, "a median body weight of 9.70 kg").
    wtref <- 9.70  # kg

    # Individual parameters. F is carried on the logit scale; keep the
    # fixed effect and eta on their own line so the term stays
    # mu-referenced.
    logitfdepot_ind <- logitfdepot + etalogitfdepot
    fdepot <- expit(logitfdepot_ind)
    ka <- exp(lka + etalka)
    d1 <- exp(ld1 + etald1)

    # Disposition parameters carry the fixed-exponent allometry only; the
    # dog fit estimated no between-subject variability on them.
    cl <- exp(lcl) * (WT / wtref)^e_wt_cl
    vc <- exp(lvc) * (WT / wtref)^e_wt_vc
    vp <- exp(lvp) * (WT / wtref)^e_wt_vp
    q  <- exp(lq)  * (WT / wtref)^e_wt_q

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment disposition with a sequential zero-order (into
    # depot) then first-order (depot -> central) absorption input. An
    # intravenous dose is given directly to `central`.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Oral dose records target `depot`: the bioavailable fraction F is
    # delivered uniformly over [0, dT1].
    f(depot)   <- fdepot
    dur(depot) <- d1

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
