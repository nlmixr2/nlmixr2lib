SoriaChacartegui_2026_tramadol <- function() {
  description <- "Two-compartment population PK model for oral tramadol in European healthy volunteers (Soria-Chacartegui 2026), with a Savic transit-compartment chain delivering the dose into a depot compartment, first-order absorption from the depot into the central compartment, first-order elimination, a linear body-weight effect on the central volume, and a CYP2D6 intermediate-metabolizer effect that lowers clearance by 19.8% relative to the pooled normal-plus-ultrarapid-metabolizer reference group."
  reference <- "Soria-Chacartegui P, Wurthwein G, Zubiaur P, Almenara S, Ochoa D, Abad-Santos F, Hempel G. Role of pharmacogenetics on tramadol pharmacokinetics: a population pharmacokinetic model. Eur J Drug Metab Pharmacokinet. 2026. doi:10.1007/s13318-026-00986-3"
  vignette <- "SoriaChacartegui_2026_tramadol"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot = list(
      analyte  = "tramadol",
      units    = "mg",
      specimen = "administration site",
      verified = TRUE,
      notes    = paste(
        "A_dep of Supplementary Figure S1: the absorption compartment that the",
        "transit chain drains into and that empties into central at rate Ka.",
        "The transit chain itself (A_0 ... A_n) is NOT carried as explicit ODE",
        "states -- the number of transit compartments is estimated as a real",
        "number (NCMT = 17.6, Table 2), so the chain is evaluated through its",
        "Savic 2007 closed-form gamma density, which is what NONMEM fits for a",
        "non-integer chain length. The dose record therefore targets depot, and",
        "f(depot) <- 0 suppresses the bolus so the whole dose arrives through",
        "the gamma input instead."
      )
    ),
    central = list(
      analyte  = "tramadol",
      units    = "mg",
      specimen = "plasma",
      verified = TRUE,
      notes    = "A_c of Supplementary Figure S1; Cc = 1000 * central / vc is the plasma tramadol concentration in ng/mL."
    ),
    peripheral1 = list(
      analyte  = "tramadol",
      units    = "mg",
      specimen = "plasma",
      verified = TRUE,
      notes    = "A_p of Supplementary Figure S1; exchanges with central through Q (K23 = Q/Vc, K32 = Q/Vp)."
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      source_name        = "Weight",
      notes              = paste(
        "Enters Vc/F -- and ONLY Vc/F -- through the PsN stepwise-covariate-",
        "modelling linear function for a continuous covariate,",
        "Vc = Vc_typical * (1 + e_wt_vc * (WT - 70)), centred on the cohort",
        "median weight of 70.0 kg (Table 1, 'Total' row). Soria-Chacartegui",
        "2026 Results 3.2 reports that the linear weight effect on CL entered",
        "the forward step of the SCM but was dropped in backward elimination,",
        "and that allometric scaling of CL and Vc was rejected because the RSEs",
        "became implausible (120% for Vc, 113% for Ka), so NO weight effect on",
        "clearance is encoded here. The cohort weight range is narrow (IQR",
        "66.2-74.2 kg; trial eligibility required BMI 18.5-30 kg/m2), so the",
        "linear extrapolation of this effect far outside that range is not",
        "supported by the data -- see the vignette Assumptions and deviations."
      )
    ),
    CYP2D6_IM = list(
      description        = "CYP2D6 intermediate-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = CYP2D6 normal or ultrarapid metabolizer (the pooled reference group)",
      source_name        = "CYP2D6 phenotype (IM versus NM + UM)",
      notes              = paste(
        "Enters CL/F through the PsN stepwise-covariate-modelling linear",
        "function for a bivariate categorical covariate (Methods 2.4),",
        "CL = CL_typical * (1 + e_cyp2d6_im_cl * CYP2D6_IM). Phenotypes were",
        "inferred from the CPIC and DPWG guidelines after CYP2D6 genotyping",
        "plus a copy-number-variation assay (Methods 2.3). The cohort held 8",
        "IMs, 14 NMs and 2 UMs, and NO poor metabolizers; the authors pooled",
        "NM with UM to gain power (Methods 2.4), so the reference category here",
        "is the pooled NM + UM group and the register's paired CYP2D6_PM",
        "indicator is deliberately NOT used -- there were no PMs to estimate a",
        "PM effect from. Setting CYP2D6_PM aside means this model must not be",
        "used to predict poor-metabolizer exposure, in whom clearance would be",
        "expected to fall further (Discussion, and CPIC 2021)."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 24L,
    n_studies      = 1L,
    age_range      = "Median 24 years (IQR 22-27)",
    age_median     = "24 years",
    weight_range   = "Median 70.0 kg (IQR 66.2-74.2); men 72.5 kg (69.2-77.8), women 57.6 kg (55.4-68.6)",
    weight_median  = "70.0 kg",
    sex_female_pct = 29.2,
    race_ethnicity = "European (self-reported biogeographic origin); all 24 volunteers",
    disease_state  = "Healthy volunteers (phase I bioequivalence trial, EUDRA-CT 2013-000196-32)",
    dose_range     = "Single oral 37.5 mg tramadol hydrochloride oral drops (Adolonta), co-administered with 400 mg ibuprofen arginine oral solution (Espidifen), under fasting conditions",
    administration = "Oral solution (drops)",
    regions        = "Spain (single centre: Hospital Universitario de La Princesa, Madrid)",
    notes          = paste(
      "Median height 1.72 m (IQR 1.67-1.77); median BMI 22.6 kg/m2 (IQR",
      "21.8-25.5); trial eligibility required BMI 18.5-30 kg/m2 (Table 1 and",
      "Discussion). 18 post-dose plasma samples per volunteer from 0.25 h to",
      "24 h plus a pre-dose baseline; LC-MS/MS with a 0.5-250 ng/mL calibration",
      "range and an LLOQ of 0.5 ng/mL, and every post-dose sample was above the",
      "LLOQ (Supplementary Material). Estimation by FOCEi in NONMEM 7.3.",
      "",
      "PHARMACOGENETIC COVARIATE SCREEN. Three phenotype covariates were tested",
      "on clearance (Methods 2.4): CYP2D6 (IM n = 8 versus pooled NM n = 14 +",
      "UM n = 2), CYP2B6 (PM n = 2 versus pooled NM n = 11 + IM n = 11) and",
      "CYP3A4 (IM n = 3 versus NM n = 20; one volunteer's phenotype was missing",
      "because of the rs4986910 SNV, whose functional effect is unclear). Only",
      "the CYP2D6 effect was retained by the stepwise covariate modelling",
      "(Results 3.2, Supplementary Figure S2). The CYP2B6 and CYP3A4 phenotype",
      "indicators are recorded here rather than in covariatesDataExcluded",
      "because no canonical covariate column exists for a CPIC-style CYP2B6",
      "poor-metabolizer indicator or for a CYP3A4 intermediate-metabolizer",
      "indicator, and minting canonical names for columns that no model file",
      "references would be speculative registration. The authors themselves",
      "attribute both negative results to the small phenotype counts and call",
      "for larger studies (Discussion, Limitations).",
      "",
      "The remaining genotyped pharmacogenes (ABCB1, ABCC2, CYP3A5, SLC22A1,",
      "UGT1A8, UGT2B7; Supplementary Table S1) were excluded from covariate",
      "testing altogether because no allele definitions or phenotype inference",
      "rules exist for them (Methods 2.4).",
      "",
      "The clearance and the volumes are apparent (/F) quantities: only oral",
      "data were collected, so bioavailability is not identifiable and no",
      "bioavailability parameter is estimated."
    )
  )

  ini({
    # -----------------------------------------------------------------------
    # Structural disposition parameters. Table 2 ('Parameter estimates of the
    # final population pharmacokinetic model'). All are apparent (/F) values;
    # Vc/F is the typical value at the 70.0 kg centring weight.
    # -----------------------------------------------------------------------
    lcl <- log(51.1); label("Apparent clearance CL/F (L/h)")                                # Table 2: CL = 51.1 L/h (RSE 6.5%; bootstrap 51.3, 95% CI 44.8-58.3)
    lvc <- log(126); label("Apparent central volume Vc/F (L) at WT = 70 kg")                # Table 2: Vc = 126 L (RSE 14.3%; bootstrap 127, 95% CI 84.0-173)
    lq <- log(175); label("Apparent intercompartmental clearance Q/F (L/h)")                # Table 2: Q = 175 L/h (RSE 14.3%; bootstrap 171, 95% CI 133-232)
    lvp <- log(171); label("Apparent peripheral volume Vp/F (L)")                           # Table 2: Vp = 171 L (RSE 8.2%; bootstrap 169, 95% CI 141-204)

    # -----------------------------------------------------------------------
    # Absorption: Savic transit chain into the depot, then first-order ka from
    # the depot into central (Results 3.2, Figure 2, Supplementary Figure S1).
    # -----------------------------------------------------------------------
    lka <- log(3.09); label("First-order absorption rate constant, depot to central (1/h)") # Table 2: Ka = 3.09 1/h (RSE 20%; bootstrap 3.15, 95% CI 2.07-5.01)
    lmtt <- log(0.24); label("Mean transit time of the absorption chain (h)")               # Table 2: MTT = 0.24 h (RSE 5.4%; bootstrap 0.24, 95% CI 0.21-0.27)
    lntr <- log(17.6); label("Number of transit compartments (continuous, dimensionless)")  # Table 2: NCMT = 17.6 (RSE 18.7%; bootstrap 17.5, 95% CI 11.3-90.1)

    # -----------------------------------------------------------------------
    # Covariate effects. Both are PsN SCM 'linear' functions, so each is a
    # FRACTIONAL shift applied multiplicatively to the typical value.
    # -----------------------------------------------------------------------
    e_cyp2d6_im_cl <- -0.198; label("Fractional change in CL/F for a CYP2D6 intermediate metabolizer (unitless)") # Table 2: 'Being CYP2D6 IM as covariate for CL' = -0.198 (RSE 30.7%; bootstrap -0.205, 95% CI -0.329 to -0.044); Results 3.2 states the 19.8% CL reduction
    e_wt_vc <- 0.0311; label("Fractional change in Vc/F per kg of body weight above 70 kg (1/kg)")                # Table 2: 'Weight as covariate for Vc' = 0.0311 (RSE 31.3%; bootstrap 0.0303, 95% CI 0.0081-0.0486)

    # -----------------------------------------------------------------------
    # Interindividual variability. Table 2 reports each IIV as a percentage
    # computed from the estimated OMEGA by footnote (3),
    #   IIV% = 100 * sqrt(exp(OMEGA^2) - 1),
    # so the variance restored here is OMEGA^2 = log(1 + (IIV%/100)^2).
    #
    # OFF-DIAGONAL: Table 2 footnote (2) and Results 3.2 state that the final
    # model carries an OMEGA BLOCK between CL and Vc (dOFV = -33 versus the
    # diagonal model; Supplementary Table S2 model 2312), but NEITHER the
    # covariance NOR the correlation is reported anywhere in the paper or the
    # supplement. The block structure is therefore preserved so that a re-fit
    # estimates the covariance, but the off-diagonal STARTING value is 0 --
    # inventing a correlation would be fabricating a parameter. Simulation
    # from this file consequently treats etalcl and etalvc as uncorrelated;
    # see the vignette Assumptions and deviations.
    # -----------------------------------------------------------------------
    # Diagonal 1 -- Table 2: IIV CL = 30.5% -> log(1 + 0.305^2) = 0.088949
    # Off-diagonal -- UNREPORTED, started at 0 (see the note above)
    # Diagonal 2 -- Table 2: IIV Vc = 68.4% -> log(1 + 0.684^2) = 0.383803
    # (Comments must stay ABOVE this block: a comment inside an omega c(...)
    #  is a known rxode2 parse failure.)
    etalcl + etalvc ~ c(
      0.088949,
      0, 0.383803
    )
    etalka ~ 0.773067                                                                       # Table 2: IIV Ka = 108% -> log(1 + 1.08^2) = 0.773067
    etalmtt ~ 0.041956                                                                      # Table 2: IIV MTT = 20.7% -> log(1 + 0.207^2) = 0.041956

    # -----------------------------------------------------------------------
    # Residual error (Results 3.2 and Table 2). The additive component was
    # non-identifiable (RSE ~48%, minimisation and boundary problems), so the
    # authors FIXED it to a negligible 0.01 ng/mL, which left the OFV and the
    # fit unchanged but kept the model stable.
    # -----------------------------------------------------------------------
    propSd <- 0.0983; label("Proportional residual error (fraction)")                        # Table 2: 'Proportional residual error' = 9.83% (RSE 5.5%; bootstrap 9.60, 95% CI 8.60-10.8)
    addSd <- fixed(0.01); label("Additive residual error (ng/mL)")                           # Table 2 footnote (1): 'Additive residual error fixed to 0.01'; Results 3.2 gives the units as ng/mL
  })

  model({
    # -----------------------------------------------------------------------
    # 1. Individual parameters.
    #
    #    Both covariate effects use the PsN 4 stepwise-covariate-modelling
    #    'linear' function (Methods 2.4), i.e. a fractional shift multiplying
    #    the typical value:
    #      * categorical (bivariate): TV * (1 + theta * I{non-reference})
    #        -- CYP2D6_IM = 1 gives CL * (1 - 0.198), the 19.8% reduction
    #        reported in Results 3.2 and reproduced by the simulated CL of
    #        41.3 versus 51.0 L/h in Table 4.
    #      * continuous: TV * (1 + theta * (COV - median(COV)))
    #        -- centred on the cohort median weight of 70.0 kg (Table 1). The
    #        centring value is confirmed by the paper's own arithmetic: the
    #        Discussion sums the typical Vc and Vp to the reported 297 L
    #        steady-state volume, which requires Vc = 126 L exactly, and the
    #        simulations in Section 2.4 are run at weight = 70 kg.
    # -----------------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (1 + e_cyp2d6_im_cl * CYP2D6_IM)
    vc <- exp(lvc + etalvc) * (1 + e_wt_vc * (WT - 70))
    q <- exp(lq)
    vp <- exp(lvp)

    ka <- exp(lka + etalka)
    mtt <- exp(lmtt + etalmtt)
    ntr <- exp(lntr)

    # Transit rate constant. Supplementary Figure S1 writes the chain as
    #   dA_0/dt   = Dose input - Ktr*A_0
    #   dA_i/dt   = Ktr*A_(i-1) - Ktr*A_i,   i = 1..n
    #   dA_dep/dt = Ktr*A_n - Ka*A_dep
    # so a dose molecule crosses n + 1 first-order steps of rate Ktr between
    # the dose record and the depot. The time to the depot is therefore
    # Erlang/gamma with shape n + 1 and mean (n + 1)/Ktr, and that mean IS the
    # reported MTT -- hence Ktr = (ntr + 1)/mtt, the Savic 2007 convention that
    # rxode2's transit() also uses. At the typical values this is
    # (17.6 + 1)/0.24 = 77.5 1/h.
    ktr <- (ntr + 1) / mtt

    # -----------------------------------------------------------------------
    # 2. Micro-constants, exactly as defined in Supplementary Figure S1.
    # -----------------------------------------------------------------------
    kel <- cl / vc
    k23 <- q / vc
    k32 <- q / vp

    # -----------------------------------------------------------------------
    # 3. Absorption input.
    #
    #    NCMT is estimated as a real number (17.6), so the transit chain cannot
    #    be written as an integer number of ODE states; the Savic 2007 gamma
    #    density that the chain converges to is evaluated instead. The density
    #    is written out explicitly rather than through rxode2's transit() macro
    #    because transit() rescales its internal dose lookup by bioavailability
    #    and therefore delivers ZERO dose when combined with the conventional
    #    f(depot) <- 0 bolus suppression in a model in nlmixr2 UI form (see the
    #    same note in Jiang_2024_empagliflozin.R and Courlet_2023_cabamiquine.R).
    #    podo(depot) is not bioavailability-adjusted, so it returns the full
    #    administered amount. Soria-Chacartegui 2026 estimates no bioavailability
    #    term -- CL and the volumes are apparent /F quantities -- so nothing
    #    scales podo() here.
    #
    #    tad() and podo() return NA until the first dose record, so the input is
    #    gated on a dose having been given.
    # -----------------------------------------------------------------------
    tdose <- tad(depot)
    if (tdose > 0) {
      transit_in <- exp(log(podo(depot)) + log(ktr) + ntr * log(ktr * tdose) -
                          ktr * tdose - lgamma(ntr + 1))
    } else {
      transit_in <- 0
    }

    # -----------------------------------------------------------------------
    # 4. ODE system (Supplementary Figure S1).
    # -----------------------------------------------------------------------
    d/dt(depot) <- transit_in - ka * depot
    d/dt(central) <- ka * depot + k32 * peripheral1 - (k23 + kel) * central
    d/dt(peripheral1) <- k23 * central - k32 * peripheral1

    # The bolus into depot is suppressed; the gamma input above delivers the
    # whole dose. MTT = 0.24 h is three orders of magnitude shorter than any
    # clinically used tramadol dosing interval, so the chain has emptied long
    # before the next dose and the most-recent-dose semantics of tad()/podo()
    # lose no mass on multiple dosing (asserted in the vignette).
    f(depot) <- 0

    # -----------------------------------------------------------------------
    # 5. Observation and error. Dose in mg over volume in L gives mg/L; the
    #    paper reports tramadol plasma concentrations in ng/mL, which is 1000x
    #    larger (1 mg/L = 1000 ng/mL). The additive term is on that ng/mL scale.
    # -----------------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
