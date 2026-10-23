Yip_2021_carbamazepine <- function() {
  description <- "Simultaneous parent + three-metabolite population PK model for carbamazepine (CBZ) and its metabolites carbamazepine-10,11-epoxide (CBZE), 2-hydroxycarbamazepine (2OH-CBZ) and 3-hydroxycarbamazepine (3OH-CBZ) in healthy adult volunteers and adults with epilepsy (Yip 2021). One-compartment parent with first-order absorption and elimination feeding three metabolite compartments that share one apparent volume of distribution (an identifiability assumption); total apparent parent clearance partitioned via a fixed metabolite fraction into CBZE, 2OH-CBZ and 3OH-CBZ formation. Covariates: binary enzyme-autoinduction status (>= 14 days of dosing) on total clearance and on 2OH-CBZ elimination, phenytoin coadministration and total daily carbamazepine dose on total clearance, and the EPHX1 rs2234922 homozygous-variant (416G/G) genotype on CBZE elimination. Concentrations are molar (mmol/L); all clearances and volumes are apparent (CL/F, V/F) because no intravenous data were available."
  reference <- "Yip VLM, Pertinez H, Meng X, Maggs JL, Carr DF, Park BK, Marson AG, Pirmohamed M. Evaluation of clinical and genetic factors in the population pharmacokinetics of carbamazepine. Br J Clin Pharmacol. 2021;87(6):2572-2588. doi:10.1111/bcp.14667"
  vignette <- "Yip_2021_carbamazepine"
  units <- list(time = "h", dosing = "mmol", concentration = "mmol/L")

  # Issue #482: what each ODE state holds, amount units, biological matrix.
  # All analytes assayed in plasma by LC-MS/MS (Yip 2021 Methods 2.4).
  # Amounts are molar (mmol); a CBZ dose in mg is converted with the CBZ
  # molecular weight 236.27 g/mol before entering the depot.
  compartmentData <- list(
    depot = list(
      analyte = "carbamazepine (CBZ)",
      units = "mmol",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "carbamazepine (CBZ)",
      units = "mmol",
      specimen = "plasma",
      verified = TRUE
    ),
    central_cbze = list(
      analyte = "carbamazepine-10,11-epoxide (CBZE)",
      units = "mmol",
      specimen = "plasma",
      verified = TRUE
    ),
    central_2oh = list(
      analyte = "2-hydroxycarbamazepine (2OH-CBZ)",
      units = "mmol",
      specimen = "plasma",
      verified = TRUE
    ),
    central_3oh = list(
      analyte = "3-hydroxycarbamazepine (3OH-CBZ)",
      units = "mmol",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    AUTOIND_STATUS = list(
      description = "Binary carbamazepine enzyme-autoinduction status for the datapoint: 1 = induced period (at least ~14 days of continuous CBZ dosing, so self-induction of metabolism is established), 0 = non-induced / first-exposure period.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-induced / first-exposure period)",
      notes = "Per-record flag reflecting the duration of prior CBZ dosing (Yip 2021 Methods 2.8 / Eqs 10-11, 17-18, and Results text). Data from PICME I (single 400 mg dose in healthy volunteers) and the first PICME II autoinduction visit were coded 0; data from the second PICME II autoinduction visit and the entire PICME II maintenance group (chronic dosing) were coded 1. The authors state a continuous autoinduction enzyme-turnover model (as used for rifampicin) failed to minimise, so a binary step function was retained. Time-varying within a subject who transitions from naive to chronic dosing. Enters as a power-of-binary multiplier on CLTOT (e_autoind_cl = 1.73) and on 2OH-CBZ elimination clearance (e_autoind_cl_2oh = 1.63).",
      source_name = "IND"
    ),
    CONMED_PHT = list(
      description = "Concomitant phenytoin at the PK observation: 1 = taking phenytoin, 0 = not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant phenytoin)",
      notes = "4 of 69 maintenance-group patients received concomitant phenytoin (Yip 2021 Table 3). Enters as a power-of-binary multiplier on total apparent CBZ clearance (e_pht_cl = 1.93, i.e. +93%), reflecting phenytoin's induction of CBZ-metabolising CYP isoforms (Yip 2021 Table 5, Eq 17).",
      source_name = "PHT"
    ),
    DOSE_CBZ_MGD = list(
      description = "Subject's own total daily carbamazepine dose in mg/day, used to model a dose-dependent increase in apparent total clearance (greater autoinduction at higher doses).",
      units = "mg/day",
      type = "continuous",
      reference_category = "normalised to 400 mg/day, the most common daily dose (Yip 2021 Methods 2.8, Eq 13); the term vanishes at 400 mg/day",
      notes = "Here carbamazepine is the index drug, so DOSE_CBZ_MGD is the modelled drug's OWN daily dose rather than a concomitant-drug dose. Enters CLTOT ADDITIVELY on the natural-log scale as ln((DOSE_CBZ_MGD/400)^e_dose_cl) with e_dose_cl = 1.36 (Yip 2021 Eqs 13 and 17), not as a multiplicative power term. Maintenance-group daily doses were 819 +/- 367 mg (range 200-1600 mg); PICME I was a single 400 mg dose. Must be strictly positive.",
      source_name = "DOSE"
    ),
    SNP_EPHX1_RS2234922_HOM = list(
      description = "EPHX1 rs2234922 (c.416A>G, Y113H) homozygous-variant indicator: 1 = 416G/G homozygous variant, 0 = otherwise (the union of 416A/A wild-type and 416A/G heterozygous). Recessive-model encoding.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (416A/A wild-type or 416A/G heterozygous)",
      notes = "EPHX1 encodes microsomal epoxide hydrolase, which hydrolyses CBZE to the dihydrodiol. The 416G/G genotype (3 of 80 subjects; Yip 2021 Table 4) reduces epoxide-hydrolase activity. Enters as a power-of-binary multiplier on CBZE elimination clearance (e_ephx1_hom_cl_cbze = 0.502, i.e. CBZE clearance reduced ~50% in 416G/G carriers; Yip 2021 Table 5, Discussion). The published Eq 18 misprints this covariate on CL3E, but Table 5, Table 6, Methods 2.8 and the Discussion all place it on CBZE elimination (CL2E) -- see vignette Errata.",
      source_name = "EPHXAGHOM"
    )
  )

  # Demographics and candidate SNPs / comedications that Yip 2021 screened
  # but did not retain in the final model. Listed for provenance only; none
  # is referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age (years).",
      units = "years",
      type = "continuous",
      notes = "Screened as a continuous covariate (centred at the population mean 45.8 years) on CLTOT, V1 and the metabolite clearances (Yip 2021 Methods 2.8, Eq 12); not retained. No demographic covariate was significant in the final model (Results 3.1)."
    ),
    WT = list(
      description = "Total body weight (kg).",
      units = "kg",
      type = "continuous",
      notes = "Screened as a continuous covariate (centred at the population mean 80.8 kg) on CLTOT, V1 and the metabolite clearances (Yip 2021 Methods 2.8, Eq 12); not retained."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male).",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not retained. Cohort 46/80 male overall (Yip 2021 Table 3)."
    ),
    SMOKE = list(
      description = "Current-smoker status (1 = smoker, 0 = non-smoker).",
      units = "(binary)",
      type = "binary",
      notes = "Screened as a categorical effect on CLTOT (Yip 2021 Methods 2.8, Eq 14); not retained. 20/69 maintenance-group patients were smokers (Table 3)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 80,
    n_studies = 2,
    n_assay_sets = 248,
    age_range = "21.7-75.8 years",
    age_mean = "45.8 years (population mean used for covariate centring)",
    weight_range = "39.0-140.5 kg",
    weight_mean = "80.8 kg (population mean used for covariate centring)",
    sex_female_pct = 42.5,
    disease_state = "Healthy adult male volunteers (PICME I, n = 8) and adults with epilepsy (PICME II: autoinduction group n = 3, maintenance group n = 69).",
    dose_range = "PICME I: single 400 mg oral immediate-release dose. PICME II autoinduction: 200 mg once daily titrated to 200 mg AM / 400 mg PM over 42 days. PICME II maintenance: 819 +/- 367 mg/day (range 200-1600 mg/day).",
    regions = "United Kingdom (Liverpool / Walton Centre).",
    notes = "Baseline demographics in Yip 2021 Table 3; SNP distribution in Table 4. Model fit in NONMEM 7.3 (FOCE-I) using closed-form analytical solutions with DOWHILE dose superposition; maintenance-group subjects were modelled as 8 weeks of BID dosing to approximate steady state. Twenty SNPs in 11 genes (ABCB1, ABCC2, UGT2B7, EPHX1, CYP2C8, CYP2C19, CYP2B6, CYP3A4, CYP3A5, MPO, POR; Table 2) and 15 comedications were screened; only EPHX1 rs2234922 416G/G (on CBZE clearance) and phenytoin and total daily dose (on CLTOT) were retained. Interoccasion variability was tested but failed to produce acceptable fits."
  )

  ini({
    # Structural PK -- parent carbamazepine (CBZ), 1-compartment with
    # first-order absorption. Point estimates from Yip 2021 Table 5.
    # lcl is the typical TOTAL apparent clearance CLTOT, which is
    # partitioned in model() into the three metabolite-formation routes.
    lcl <- log(1.96); label("Typical total apparent carbamazepine clearance CLTOT/F (L/h)") # Yip 2021 Table 5: CLTOT = 1.96 L/h
    lvc <- log(161); label("Apparent central volume of distribution V1/F for CBZ (L)") # Yip 2021 Table 5: V1 = 161 L
    lkastar <- log(0.409); label("Absorption-rate offset KAstar above the elimination rate (1/h)") # Yip 2021 Table 5: KASTAR = 0.409 /h; KA = CLTOT/V1 + KAstar (Eq 9)

    # Fraction of parent clearance forming the two hydroxy metabolites
    # (2OH-CBZ + 3OH-CBZ). CL2 (to CBZE) = CLTOT*(1-FM1); CL3 = CL4
    # (to 2OH, 3OH) = CLTOT*FM1/2 (Eqs 6-8). FIXED at 0.64 for
    # identifiability (36% of CBZ -> CBZE by prior mass balance); a small
    # IIV was still estimated on it (Table 5).
    lfm1 <- fixed(log(0.64)); label("Log fraction of parent clearance to the 2OH + 3OH hydroxy metabolites FM1 (unitless)") # Yip 2021 Table 5 + Methods 2.8: FM1 fixed to 0.64

    # Structural PK -- metabolite elimination clearances (apparent), all
    # sharing one apparent metabolite volume lvc_met (= V2; shared across
    # the three metabolites for identifiability, Methods 2.8).
    lcl_cbze <- log(9.71); label("Apparent CBZE elimination clearance CL2E/F (L/h)") # Yip 2021 Table 5: CL2E = 9.71 L/h
    lcl_2oh <- log(726); label("Apparent 2OH-CBZ elimination clearance CL3E/F (L/h)") # Yip 2021 Table 5: CL3E = 726 L/h
    lcl_3oh <- log(543); label("Apparent 3OH-CBZ elimination clearance CL4E/F (L/h)") # Yip 2021 Table 5: CL4E = 543 L/h
    lvc_met <- log(39.2); label("Shared apparent volume of distribution V2/F for all three metabolites (L)") # Yip 2021 Table 5: V2 = 39.2 L

    # Covariate effects (power-of-binary multipliers unless noted).
    e_autoind_cl <- 1.73; label("Autoinduction multiplier on CLTOT (unitless)") # Yip 2021 Table 5: CLTOT induction effect = 1.73; Eq 17
    e_autoind_cl_2oh <- 1.63; label("Autoinduction multiplier on 2OH-CBZ elimination clearance (unitless)") # Yip 2021 Table 5: CL3E induction effect = 1.63; Eq 18
    e_pht_cl <- 1.93; label("Phenytoin-coadministration multiplier on CLTOT (unitless)") # Yip 2021 Table 5: Cov.Eff.(PHT on CLTOT) = 1.93; Eq 17
    e_dose_cl <- 1.36; label("Daily-dose exponent on CLTOT, additive ln term (unitless)") # Yip 2021 Table 5: Cov.Eff.(DOSE on CLTOT) = 1.36; Eqs 13, 17
    e_ephx1_hom_cl_cbze <- 0.502; label("EPHX1 416G/G multiplier on CBZE elimination clearance (unitless)") # Yip 2021 Table 5: Cov.Eff.(EPHXAGHOM on CL2E) = 0.502; see Errata re Eq 18

    # Interindividual variability. Yip 2021 Table 5 reports IIV as % (CV);
    # converted to log-normal variance via omega^2 = log(1 + CV^2).
    # CLTOT 25.3%: log(1 + 0.253^2) = 0.06205
    # FM1    6.3%: log(1 + 0.063^2) = 0.003961
    # CL2E  27.9%: log(1 + 0.279^2) = 0.07497
    # CL3E  42.3%: log(1 + 0.423^2) = 0.16459
    # CL4E  34.5%: log(1 + 0.345^2) = 0.11244
    # KASTAR 75.3%: log(1 + 0.753^2) = 0.44919
    # V1    13.6%: log(1 + 0.136^2) = 0.018327
    # V2   149.7%: log(1 + 1.497^2) = 1.17595
    etalcl ~ 0.062044 # Yip 2021 Table 5: 25.3% IIV on CLTOT
    etalfm1 ~ 0.003961 # Yip 2021 Table 5: 6.3% IIV on FM1 (ETA-shrinkage 65%)
    etalcl_cbze ~ 0.074960 # Yip 2021 Table 5: 27.9% IIV on CL2E (CBZE elimination)
    etalcl_2oh ~ 0.164606 # Yip 2021 Table 5: 42.3% IIV on CL3E (2OH-CBZ elimination)
    etalcl_3oh ~ 0.112458 # Yip 2021 Table 5: 34.5% IIV on CL4E (3OH-CBZ elimination)
    etalkastar ~ 0.449169 # Yip 2021 Table 5: 75.3% IIV on KASTAR
    etalvc ~ 0.018327 # Yip 2021 Table 5: 13.6% IIV on V1
    etalvc_met ~ 1.175885 # Yip 2021 Table 5: 149.7% IIV on V2 (ETA-shrinkage 64%)

    # Residual error -- exponential (log-normal) error, one SD per analyte
    # (Yip 2021 Methods 2.8, Table 5 'resid err.Sd').
    expSd <- 0.221; label("CBZ exponential residual SD (log-scale)") # Yip 2021 Table 5: CBZ resid err.Sd = 0.221
    expSd_cbze <- 0.310; label("CBZE exponential residual SD (log-scale)") # Yip 2021 Table 5: CBZE resid err.Sd = 0.310
    expSd_2oh <- 0.297; label("2OH-CBZ exponential residual SD (log-scale)") # Yip 2021 Table 5: 2OH-CBZ resid err.Sd = 0.297
    expSd_3oh <- 0.362; label("3OH-CBZ exponential residual SD (log-scale)") # Yip 2021 Table 5: 3OH-CBZ resid err.Sd = 0.362
  })

  model({
    # 1. Individual total apparent parent clearance (Yip 2021 Eq 17). The
    #    autoinduction and phenytoin effects enter multiplicatively; the
    #    daily-dose effect enters ADDITIVELY as an ln term (Eq 13), and the
    #    whole typical value is scaled by exp(etalcl).
    cltot <- (exp(lcl) * e_autoind_cl^AUTOIND_STATUS * e_pht_cl^CONMED_PHT +
      log((DOSE_CBZ_MGD / 400)^e_dose_cl)) * exp(etalcl)

    # 2. Remaining individual parameters.
    vc <- exp(lvc + etalvc)
    vmet <- exp(lvc_met + etalvc_met)
    fm1 <- exp(lfm1 + etalfm1)
    kastar <- exp(lkastar + etalkastar)

    # Metabolite elimination clearances with their covariate effects.
    cl_cbze <- exp(lcl_cbze + etalcl_cbze) * e_ephx1_hom_cl_cbze^SNP_EPHX1_RS2234922_HOM
    cl_2oh <- exp(lcl_2oh + etalcl_2oh) * e_autoind_cl_2oh^AUTOIND_STATUS
    cl_3oh <- exp(lcl_3oh + etalcl_3oh)

    # 3. Micro-constants. KA constrained above the parent elimination rate
    #    kel via the KAstar offset (Yip 2021 Eq 9), which prevents flip-flop.
    kel <- cltot / vc
    ka <- kel + kastar
    kel_cbze <- cl_cbze / vmet
    kel_2oh <- cl_2oh / vmet
    kel_3oh <- cl_3oh / vmet

    # 4. ODE system (Yip 2021 Eqs 1-5). Parent elimination kel*central is
    #    fully partitioned into the three metabolite-formation fluxes:
    #    (1-fm1) -> CBZE, fm1/2 -> 2OH-CBZ, fm1/2 -> 3OH-CBZ. Mole-for-mole
    #    conversion (concentrations are molar).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    d/dt(central_cbze) <- (1 - fm1) * kel * central - kel_cbze * central_cbze
    d/dt(central_2oh) <- (fm1 / 2) * kel * central - kel_2oh * central_2oh
    d/dt(central_3oh) <- (fm1 / 2) * kel * central - kel_3oh * central_3oh

    # 5. Observations (mmol/L) and exponential (log-normal) residual error.
    Cc <- central / vc
    Cc_cbze <- central_cbze / vmet
    Cc_2oh <- central_2oh / vmet
    Cc_3oh <- central_3oh / vmet

    Cc ~ lnorm(expSd)
    Cc_cbze ~ lnorm(expSd_cbze)
    Cc_2oh ~ lnorm(expSd_2oh)
    Cc_3oh ~ lnorm(expSd_3oh)
  })
}
