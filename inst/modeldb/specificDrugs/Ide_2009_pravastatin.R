Ide_2009_pravastatin <- function() {
  description <- "Population PK model for orally administered pravastatin with enterohepatic circulation (Ide 2009) in healthy Japanese male volunteers. Absorption is described by an Erlang chain of 8 transit compartments (N_depot = 8); disposition is one-compartment central with a gallbladder recirculation compartment whose release is gated by the gallbladder-emptying time tg (continuous filling from central via k12 for t < tg, gated release to central via k21 for t >= tg) producing the characteristic second-peak phenomenon. SLCO1B1 *15 haplotype carrier status (paired heterozygote / homozygote indicators) increases relative oral bioavailability Frel multiplicatively (1.50x and 1.95x respectively). Gastric conversion of pravastatin to its inactive 3'alpha-isopravastatin (RMS-416) is highly variable; the source paper corrected for this by using an apparent dose (actual dose x Fa, where Fa = AUCpra / (AUCpra + AUCrms)) as the model input, so the packaged model fixes the depot bioavailability anchor at the population-mean Fa = 0.571 derived from Table II mean AUC values."
  reference <- paste(
    "Ide T, Sasaki T, Maeda K, Higuchi S, Sugiyama Y, Ieiri I.",
    "Quantitative Population Pharmacokinetic Analysis of Pravastatin",
    "Using an Enterohepatic Circulation Model Combined With",
    "Pharmacogenomic Information on SLCO1B1 and ABCC2 Polymorphisms.",
    "J Clin Pharmacol. 2009;49(11):1309-1317.",
    "doi:10.1177/0091270009341960.",
    sep = " "
  )
  vignette <- "Ide_2009_pravastatin"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    SLCO1B1_HAP15_HET = list(
      description        = "SLCO1B1 *15 haplotype heterozygote indicator: 1 if subject carries exactly one *15 allele (diplotype *1a/*15 or *1b/*15), 0 otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (*15-noncarrier: *1a/*1a, *1a/*1b, or *1b/*1b)",
      notes              = "Time-fixed germline haplotype. The *15 haplotype is defined by the cis combination of the 388A>G (rs2306283, N130D) and 521T>C (rs4149056, V174A) variants. Paper covariate equation Eq. on p. 1311: Pi = TVP * theta1^HT * theta2^HM where HT = SLCO1B1_HAP15_HET and HM = SLCO1B1_HAP15_HOM. Paired with SLCO1B1_HAP15_HOM. Distribution in Ide 2009 (Table I, n = 57): 28 noncarriers (49%), 23 heterozygotes (40%), 6 homozygotes (11%).",
      source_name        = "HT"
    ),
    SLCO1B1_HAP15_HOM = list(
      description        = "SLCO1B1 *15 haplotype homozygote indicator: 1 if subject carries two *15 alleles (diplotype *15/*15), 0 otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (*15-noncarrier: *1a/*1a, *1a/*1b, or *1b/*1b)",
      notes              = "Time-fixed germline haplotype. See SLCO1B1_HAP15_HET notes for the *15 haplotype definition and the paper's covariate equation. Paired with SLCO1B1_HAP15_HET to encode a three-level haplotype categorical (noncarrier / heterozygote / homozygote) with *15-noncarrier as the implicit reference (both indicators = 0). Distribution in Ide 2009 (Table I): 6 of 57 (10.5%) homozygotes.",
      source_name        = "HM"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 57L,
    n_observations = 636L,
    n_studies      = 1L,
    age_range      = "20-40 years",
    weight_range   = "49.7-97.8 kg",
    sex_female_pct = 0,
    race_ethnicity = "Japanese (100%)",
    disease_state  = "Healthy volunteers; pooled from three previous studies by the same group (Nishizato 2003, Maeda 2006, Suwannakul 2008) prescreened by SLCO1B1 haplotype to enrich for *15 carriers.",
    dose_range     = "Single 10 mg oral dose (Daiichi-Sankyo Co Ltd, Tokyo) with 150 mL water after overnight fast; food given 4 hours post-dose.",
    regions        = "Japan (Kyushu University, Kyushu Pharmacology Research Clinic, Tottori University, University of Tokyo).",
    genotype_distribution = "SLCO1B1 (Table I): 7 *1a/*1a, 6 *1a/*1b, 15 *1b/*1b, 7 *1a/*15, 16 *1b/*15, 6 *15/*15. ABCC2 -24C>T: 40 C/C, 17 C/T, 0 T/T. ABCC2 1249G>A: 43 G/G, 12 G/A, 2 A/A. ABCC2 1446C>G: 57 C/C, 0 C/G, 0 G/G. ABCC2 3972C>T: 37 C/C, 19 C/T, 1 T/T.",
    sampling_window = "Rich sampling at 0, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6, 8, 12, 24 h post-dose in 33 subjects and 0, 0.25, 0.5, 0.75, 1, 2, 4, 6, 8, 12, 24 h post-dose in 24 subjects (Methods).",
    notes          = "Haplotype frequencies were not natural Japanese-population frequencies; subjects were recruited from a pool of approximately 100 prescreened Japanese volunteers to enrich for *15 carriers. BQL rates: 52 of 636 (8.2%) pravastatin observations and 92 of 636 (14.5%) RMS-416 observations were below LLOQ and omitted. The paper analysed pravastatin parent kinetics only; RMS-416 (3'alpha-isopravastatin, the inactive gastric biotransformation product) was used to construct the apparent-bioavailability correction Fa but was not separately modelled. Body weight was tested as a covariate and rejected (no significant effect on disposition)."
  )

  # Implementation notes (see vignette 'Assumptions and deviations' for the
  # full justification):
  # * Apparent-dose framework. Ide 2009 used "apparent dose" = actual dose * Fa
  #   as the NONMEM model input, where Fa = AUCpra / (AUCpra + AUCrms) is a
  #   per-subject correction for the highly variable gastric conversion of
  #   pravastatin to RMS-416 (the inactive 3'alpha-isopravastatin isomer).
  #   The packaged model fixes lfdepot at log(0.571), the population-mean Fa
  #   computed from Table II mean AUC values (56.2 / (56.2 + 42.3) = 0.571).
  #   This reproduces the paper's typical-value predictions when a user
  #   simulates the actual 10 mg oral dose. Users with per-subject AUC ratios
  #   can override f(depot) externally to replicate individual-level fits.
  # * Gallbladder gating. The paper's K21 (gallbladder -> central rate) is
  #   conceptually gated by tg (gallbladder-emptying time): the gallbladder
  #   accumulates drug via K12 continuously from t = 0 and begins releasing
  #   only at t >= tg. The packaged model encodes this gate via the
  #   multiplicative (t >= tg) factor on the K21 term. This is faithful for
  #   single-dose simulation (the paper's design). For multi-dose simulation
  #   with inter-dose intervals >> tg, the gallbladder release is continuous
  #   after the first dose; the model is therefore best suited to single-dose
  #   simulation. Note that nlmixr2 / rxode2 expose absolute simulation time
  #   as the lowercase variable `t`.
  # * Rate-constant parameterisation. The central-to-gallbladder
  #   bidirectional flow is parameterised as the rate constants k12 and k21
  #   (paper's notation), not as Q / Vp -- the gallbladder is mechanism-
  #   defined (a delayed-release container) rather than a volume-defined
  #   peripheral compartment, so the canonical lq / lvp form does not apply.
  #   `lk12` and `lk21` deviate from the canonical lq / lvp pattern; the
  #   deviation is paper-faithful and consistent with prior EHC popPK
  #   extractions (see Benkali_2010_tacrolimus which uses lkcp / lkpc for
  #   similar reasons).
  # * `tg` parameter. Encoded as `ltg` (log-transformed gallbladder-emptying
  #   time). Not in the canonical PK-parameter register; paper-named.
  ini({
    # Erlang absorption (Ide 2009 Table II, Final Model).
    # Mean transit time MTT and number of transit compartments N_depot = 8
    # (text p. 1312: "the EHC model with Erlang's distribution for the
    # absorption phase (Ndepot = 8) presented the best description of the
    # data"). The transit rate constant ktr is recovered inside model() as
    # ktr = 8 / mtt.
    lmtt <- log(0.660); label("Mean transit time MTT (h)")                                          # Table II Final Model MTT = 0.660 h (RSE 4.9%)

    # Apparent clearance and central volume of distribution (Ide 2009 Table II).
    # The reported CL/F and V1/F are the values that map the "apparent dose"
    # (actual dose * Fa) onto observed concentrations -- see lfdepot below
    # for the population-mean Fa anchor.
    lcl <- log(139); label("Apparent clearance CL/F (L/h)")                                          # Table II Final Model CL/F = 139 L/h (RSE 6.4%)
    lvc <- log(253); label("Apparent central volume V1/F (L)")                                       # Table II Final Model V1/F = 253 L (RSE 8.1%)

    # Enterohepatic circulation rate constants.
    # k12: continuous biliary excretion from central to gallbladder.
    # k21: post-emptying release from gallbladder to central (gated by t >= tg
    #      in the model() body, matching the paper's gallbladder-emptying
    #      mechanism).
    lk12 <- log(0.183); label("Central -> gallbladder rate constant k12 (1/h)")                      # Table II Final Model K12 = 0.183 1/h (RSE 10.3%)
    lk21 <- log(0.436); label("Gallbladder -> central rate constant k21 (1/h, gated by t >= tg)")    # Table II Final Model K21 = 0.436 1/h (RSE 6.8%)
    ltg  <- log(3.61);  label("Gallbladder-emptying time tg (h)")                                    # Table II Final Model tg = 3.61 h (RSE 1.6%)

    # Bioavailability anchor: population-mean apparent bioavailability Fa
    # from gastric conversion. Fa = AUCpra / (AUCpra + AUCrms); using Table II
    # cohort means (56.2 + 42.3 ng*h/mL) -> Fa = 56.2 / 98.5 = 0.5706.
    # FIXED because the paper did not estimate Fa as a model parameter; it
    # was computed externally per subject from observed AUCs. See vignette
    # Assumptions and deviations.
    lfdepot <- fixed(log(0.571)); label("Population-mean Fa bioavailability anchor (unitless)")      # Derived from Table II mean AUCpra = 56.2 and AUCrms = 42.3 ng*h/mL (p. 1311)

    # SLCO1B1 *15 haplotype effects on relative bioavailability Frel.
    # Paper covariate equation (p. 1311):
    #   Frel = 1 * theta1^HT * theta2^HM   (Pi = TVP * theta1^HT * theta2^HM)
    # Inside model() this becomes:
    #   Frel = e_slco1b1_hap15_het_fdepot^SLCO1B1_HAP15_HET *
    #          e_slco1b1_hap15_hom_fdepot^SLCO1B1_HAP15_HOM
    e_slco1b1_hap15_het_fdepot <- 1.50; label("SLCO1B1 *15 heterozygote fold-effect on Frel (unitless)") # Table II Final Model Frel SLCO1B1*15 heterozygote = 1.50 (RSE 8.8%)
    e_slco1b1_hap15_hom_fdepot <- 1.95; label("SLCO1B1 *15 homozygote fold-effect on Frel (unitless)")   # Table II Final Model Frel SLCO1B1*15 homozygote = 1.95 (RSE 15.6%)

    # Inter-individual variability (exponential IIV models per Methods p. 1311:
    # "Exponential and combined (ie, proportional and additive) variance
    # models were used to describe the interindividual and residual
    # variability, respectively."). Paper-reported %CV translates to log-scale
    # variance via omega^2 = log(1 + CV^2):
    #   MTT  CV 35.2% -> log(1 + 0.352^2) = 0.1168
    #   CL/F CV 27.9% -> log(1 + 0.279^2) = 0.0749
    #   V1/F CV 34.4% -> log(1 + 0.344^2) = 0.1119
    #   K12  CV 51.2% -> log(1 + 0.512^2) = 0.2328
    # Final Model retained IIV on MTT, CL/F, V1/F, K12; IIV on K21 and tg was
    # tested and not retained (Results p. 1312: "Introduction of the
    # interindividual variability in K21 and tg did not improve the model").
    etalmtt ~ 0.1168                                                                                  # Table II Final Model IIV MTT  = 35.2% CV
    etalcl  ~ 0.0749                                                                                  # Table II Final Model IIV CL/F = 27.9% CV
    etalvc  ~ 0.1119                                                                                  # Table II Final Model IIV V1/F = 34.4% CV
    etalk12 ~ 0.2328                                                                                  # Table II Final Model IIV K12  = 51.2% CV

    # Residual error (combined proportional + additive, Ide 2009 Table II).
    # Concentration units are ng/mL; the additive error magnitude is reported
    # on that scale.
    propSd <- 0.188;  label("Proportional residual error (fraction)")                                # Table II Final Model proportional error = 18.8% (RSE 14.6%)
    addSd  <- 0.176;  label("Additive residual error (ng/mL)")                                       # Table II Final Model additive error = 0.176 ng/mL (RSE 30.7%)
  })
  model({
    # Individual PK parameters (exponential IIV)
    mtt <- exp(lmtt + etalmtt)
    cl  <- exp(lcl  + etalcl)
    vc  <- exp(lvc  + etalvc)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21)
    tg  <- exp(ltg)

    # Erlang transit absorption: chain of N_depot = 8 compartments linked
    # by the common rate constant ktr = N_depot / MTT (Savic 2007 formulation,
    # adopted by Ide 2009 Eq. on p. 1310: Ktr = Ndepot / MTT). Implemented as
    # depot + transit1 + ... + transit7 = 8 compartments before central.
    ktr <- 8 / mtt

    # Micro-constant for elimination from central
    kel <- cl / vc

    # SLCO1B1 *15 haplotype multiplicative effect on relative bioavailability.
    # Paper covariate equation (p. 1311): Frel = 1 * 1.50^HT * 1.95^HM.
    frel <- e_slco1b1_hap15_het_fdepot^SLCO1B1_HAP15_HET *
            e_slco1b1_hap15_hom_fdepot^SLCO1B1_HAP15_HOM

    # Gallbladder release is gated by the gallbladder-emptying time tg
    # (paper's "tg" parameter; Final Model tg = 3.61 h). The gallbladder
    # accumulates from t = 0 via k12 and begins releasing via k21 only at
    # t >= tg. `t` is rxode2's absolute time variable.
    flow_c_to_g <- k12 * central
    flow_g_to_c <- k21 * gallbladder * (t >= tg)

    # ODE system
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2 - ktr * transit3
    d/dt(transit4)    <-  ktr * transit3 - ktr * transit4
    d/dt(transit5)    <-  ktr * transit4 - ktr * transit5
    d/dt(transit6)    <-  ktr * transit5 - ktr * transit6
    d/dt(transit7)    <-  ktr * transit6 - ktr * transit7
    d/dt(central)     <-  ktr * transit7 - kel * central - flow_c_to_g + flow_g_to_c
    d/dt(gallbladder) <-  flow_c_to_g - flow_g_to_c

    # Bioavailability anchor (population-mean Fa) multiplied by the
    # SLCO1B1 *15 Frel multiplier. The amount entering depot for an
    # actual oral dose D is D * exp(lfdepot) * frel.
    f(depot) <- exp(lfdepot) * frel

    # Pravastatin plasma concentration. central is in mg (dose units); vc is in L.
    # central/vc gives mg/L; 1 mg/L = 1000 ng/mL, so multiply by 1000 to match
    # the source paper's ng/mL bioanalytical units.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd) + add(addSd)
  })
}
