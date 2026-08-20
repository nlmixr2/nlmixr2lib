PillayFuentesLorente_2024_albendazole <- function() {
  description <- "Joint population PK model for the two main albendazole metabolites after a single oral 400 mg albendazole dose (with or without co-administered ivermectin) in adolescents infected with Trichuris trichiura in Tanzania and Cote d'Ivoire: Savic transit-compartment absorption feeding a two-compartment albendazole sulfoxide disposition that converts quantitatively to a one-compartment albendazole sulfone, with a study-population (country) effect on both apparent clearances"
  reference <- paste(
    "Pillay-Fuentes Lorente V, Nwogu-Attah JN, Steffens B, Bram D, Sprecher V,",
    "Hofmann D, Buettcher M, Pillai G, Mouksassi S, Coulibaly J, Pfister M, Keiser J.",
    "Understanding Drug Exposure and Trichuris trichiura Cure Rates:",
    "A Pharmacometric Approach for Albendazole-Ivermectin Co-medication",
    "in Tanzania and Cote d'Ivoire.",
    "Drugs R D. 2024;24(2):205-215.",
    "doi:10.1007/s40268-024-00476-4.",
    sep = " "
  )
  vignette <- "PillayFuentesLorente_2024_albendazole"

  # Albendazole sulfone is a paper-specific second analyte rather than a
  # generalisable payload / metabolite class, so its compartment is declared
  # here instead of promoting an `abzson` suffix into
  # R/conventions.R::registeredMetabolites. Same mechanism (and same
  # rationale) as the `central_cms` declaration in Mohamed_2012_colistin.R.
  paper_specific_compartments <- c("central_abzson")

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The parent drug albendazole is NOT a state in this model: up to 35% of the
  # albendazole concentrations were below the limit of quantification and
  # including albendazole did not give a good fit (Results 3.2), so the
  # published final model carries only the two metabolites. The authors assumed
  # 100% conversion of albendazole to albendazole sulfoxide and of albendazole
  # sulfoxide to albendazole sulfone (Methods 2.2), so the administered
  # albendazole dose enters the absorption chain unchanged and every unit of
  # sulfoxide cleared becomes a unit of sulfone. No molar-mass correction is
  # applied because the published apparent volumes and clearances were
  # estimated on that basis.
  compartmentData <- list(
    depot = list(
      analyte = "albendazole", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "albendazole sulfoxide", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "albendazole sulfoxide", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    central_abzson = list(
      analyte = "albendazole sulfone", units = "mg",
      specimen = "plasma", verified = TRUE
    )
  )

  covariateData <- list(
    REGION_TANZANIA = list(
      description        = "Study population / enrollment-country indicator (1 = Tanzania, 0 = Cote d'Ivoire)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Cote d'Ivoire)",
      notes              = paste(
        "The only covariate retained in the published final model (Results 3.3).",
        "Applied as an exponential categorical effect on the apparent clearance",
        "of BOTH metabolites: exp(0.56) = 1.75 (75% higher albendazole",
        "sulfoxide CL/F in Tanzania) and exp(0.38) = 1.46 (46% higher",
        "albendazole sulfone CL/F in Tanzania), matching the percentages quoted",
        "in the Abstract, Key Points, Results 3.3 and Discussion.",
        "Time-constant per subject."
      ),
      source_name        = "Country"
    )
  )

  # Covariates screened during model development but NOT retained in the final
  # model (Methods 2.3, Results 3.3). Documented here so the paper's covariate
  # screen is preserved without declaring covariates that model() never uses.
  covariatesDataExcluded <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Tested both as a free power-function covariate on CL and V and as",
        "theory-based allometric scaling with exponents fixed to 0.75 (CL) and",
        "1 (V). Neither improved the fit, so allometric scaling was removed",
        "from the final model (Results 3.3). Means 51.8 kg (SD 8.18) in",
        "Tanzania and 44.4 kg (SD 9.52) in Cote d'Ivoire (Table 1)."
      ),
      source_name        = "WT"
    ),
    BMI = list(
      description        = "Body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened on CL and V of both metabolites; did not improve the model",
        "fit (Results 3.3). Means 19.9 kg/m^2 (SD 2.19) in Tanzania and",
        "18.9 kg/m^2 (SD 2.76) in Cote d'Ivoire (Table 1)."
      ),
      source_name        = "BMI"
    ),
    SEXF = list(
      description        = "Sex (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Screened as a categorical covariate on CL and V; did not improve the",
        "model fit (Results 3.3). 44.4% female in Tanzania and 47.1% female in",
        "Cote d'Ivoire (Table 1)."
      ),
      source_name        = "Sex"
    ),
    CONMED_IVERMECTIN = list(
      description        = "Co-administered ivermectin 200 ug/kg single oral dose (1 = albendazole + ivermectin, 0 = albendazole alone)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (albendazole alone)",
      notes              = paste(
        "The paper's primary co-medication hypothesis. Tested on the apparent",
        "clearance of both metabolites and REJECTED: inclusion raised the",
        "corrected Bayesian information criterion by 9.07, so it was not",
        "retained (Results 3.3). 59.3% of the Tanzanian and 47.1% of the",
        "Ivorian participants received the combination (Table 1). Ivermectin",
        "concentrations were not available at the time of model development,",
        "so ivermectin PK is not part of this model."
      ),
      source_name        = "Medication"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 44L,
    n_studies      = 2L,
    age_range      = "12-19 years",
    age_median     = "mean 15.8 years (SD 1.37) in Tanzania; mean 14.6 years (SD 1.86) in Cote d'Ivoire (Table 1)",
    weight_range   = "mean 51.8 kg (SD 8.18) in Tanzania; mean 44.4 kg (SD 9.52) in Cote d'Ivoire (Table 1)",
    sex_female_pct = 45.5,
    disease_state  = "adolescents with confirmed Trichuris trichiura (human whipworm) infection",
    dose_range     = "single oral albendazole 400 mg, alone or co-administered with a single oral ivermectin 200 ug/kg dose, taken after a high-fat meal",
    regions        = "Tanzania (Pemba Island) and Cote d'Ivoire",
    co_medication  = "ivermectin 200 ug/kg single oral dose in 24 of 44 participants; no effect on either apparent clearance was retained",
    notes          = paste(
      "Pooled pharmacokinetic data from two phase III trials, 227 concentration",
      "records from 44 adolescents (Results 3.1). Sampling schedules differed by",
      "country and are strongly informative for which phase of the profile each",
      "cohort contributes: Tanzania n = 27, 108 samples at 6, 21, 27 and 48 h",
      "post-dose (elimination phase); Cote d'Ivoire n = 17, 119 samples at 1, 2,",
      "4, 6, 8, 22 and 27 h post-dose (absorption phase) (Methods 2.1).",
      "Albendazole, albendazole sulfoxide and albendazole sulfone were assayed",
      "in duplicate by validated LC-MS/MS (Schulz 2019, doi:10.1128/AAC.02489-18)."
    )
  )

  ini({
    # ==================================================================
    # Absorption. Savic (2007) transit-compartment chain feeding a depot
    # that absorbs into the albendazole sulfoxide central compartment at
    # rate ka. Monolix parameterises the chain by the mean transit time
    # Mtt and the transit rate constant Ktr; the implied (non-integer)
    # number of transit compartments is nn = Mtt * Ktr - 1, derived in
    # model() rather than stored here, because the paper reports Mtt and
    # Ktr and not nn.
    #
    # lktr and lmtt are wrapped in fixed() even though the paper ESTIMATED
    # both (RSE 6.01% and 8.72%). At the published values the transit chain
    # is degenerate and is encoded at its limiting form (see the model()
    # comment and the vignette), so Ktr and Mtt do not enter any ODE here:
    # they are carried for provenance only. Leaving them estimable would
    # ship a structurally singular refit (two parameters with identically
    # zero gradient). This is a deliberate, documented deviation - the
    # published values and their estimated status are preserved here and in
    # the vignette's Assumptions and deviations section.
    # ==================================================================
    lka  <- fixed(log(0.2))    ; label("Absorption rate constant ka (1/h)")     # Table 2 'Absorption rate constant, [Ka (h-1)] = 0.2 (fix)'; Results 3.2 'fixed to 0.2 based on the estimate from a previous model (unpublished data)'
    lktr <- fixed(log(0.12))   ; label("Transit rate constant Ktr (1/h)")       # Table 2 'Transit rate constant, [Ktr (h-1)] = 0.12' (RSE 6.01%); estimated in the paper, fixed here (inert in this encoding)
    lmtt <- fixed(log(0.095))  ; label("Mean transit time Mtt (h)")             # Table 2 'Mean transit time, [Mtt (h)] = 0.095' (RSE 8.72%); estimated in the paper, fixed here (inert in this encoding)

    # ==================================================================
    # Albendazole sulfoxide disposition (two compartments). All volumes
    # and clearances are apparent (F-relative) because albendazole is
    # given orally and the fraction absorbed and metabolised is not
    # separately identifiable.
    # ==================================================================
    lcl  <- log(2.09)          ; label("Apparent clearance of albendazole sulfoxide Clox/F (L/h)")                        # Table 2 'Clearance of albendazole sulfoxide, [Clox (L h-1)] = 2.09' (RSE 5.62%)
    lvc  <- log(81.81)         ; label("Apparent central volume of albendazole sulfoxide Vox/F (L)")                      # Table 2 'Volume of albendazole sulfoxide in the CC, [Vox (L)] = 81.81' (RSE 10.8%)
    lq   <- log(20.85)         ; label("Apparent inter-compartmental clearance of albendazole sulfoxide Qox/F (L/h)")     # Table 2 'Intercompartmental clearance for sulfoxide, [Qox (L h-1)] = 20.85' (RSE 11.30%)
    lvp  <- log(1931.53)       ; label("Apparent peripheral volume of albendazole sulfoxide Pvox/F (L)")                  # Table 2 'Volume of sulfoxide in the PC, [Pvox (L)] = 1931.53' (RSE 5.61%)

    # ==================================================================
    # Albendazole sulfone disposition (one compartment). Formation is
    # assumed quantitative from albendazole sulfoxide (Methods 2.2), so
    # the sulfone input rate is the full sulfoxide elimination rate.
    # ==================================================================
    lcl_abzson <- log(22.87)   ; label("Apparent clearance of albendazole sulfone Clon/F (L/h)")                          # Table 2 'Clearance of albendazole sulfone, [Clon (L h-1)] = 22.87' (RSE 5.38%)
    lvc_abzson <- log(11.51)   ; label("Apparent central volume of albendazole sulfone Von/F (L)")                        # Table 2 'Volume of albendazole sulfone in the CC, [Von (L)] = 11.51' (RSE 26.9%)

    # ==================================================================
    # Study-population (country) effect on the two apparent clearances.
    # Methods 2.3 Eq. 2 gives the categorical covariate model as
    # P_i = Ppop * exp(beta_cat_1) for the non-reference category, so the
    # tabulated coefficients are on the log scale with Cote d'Ivoire as
    # the reference. exp(0.56) = 1.751 and exp(0.38) = 1.462 reproduce
    # the "75% higher" and "46% higher" statements in the Abstract, Key
    # Points, Results 3.3 and Discussion.
    # ==================================================================
    e_region_tanzania_cl        <- 0.56 ; label("Log-scale effect of Tanzanian study population on albendazole sulfoxide apparent clearance (unitless)")  # Table 2 'Clox-Country (Tanzania) = 0.56' (RSE 18.2%)
    e_region_tanzania_cl_abzson <- 0.38 ; label("Log-scale effect of Tanzanian study population on albendazole sulfone apparent clearance (unitless)")    # Table 2 'Clon-Country (Tanzania) = 0.38' (RSE 17.5%)

    # ==================================================================
    # Inter-individual variability. Monolix reports omega as the standard
    # deviation of the normally distributed random effect on the log
    # scale for log-normally distributed parameters (Methods 2.2,
    # "Population modeling was performed utilizing a lognormal
    # distribution on PK parameters"), so the nlmixr2 variance is
    # omega^2. The paper reports no IIV on Ktr, Mtt or the albendazole
    # sulfone clearance, so none is encoded.
    # ==================================================================
    etalvc        ~ 0.2704   # omega(Vox)  = 0.52 (RSE 16.0%), Table 2 'BSV on volume of ALB-OX in the CP';                 var = 0.52^2  = 0.2704
    etalcl        ~ 0.0121   # omega(Clox) = 0.11 (RSE 33.3%), Table 2 'Variability on clearance of ALB-OX in the CP';       var = 0.11^2  = 0.0121
    etalq         ~ 0.4761   # omega(Qox)  = 0.69 (RSE 12.3%), Table 2 'Variability on intercompartmental clearance of ALB-OX'; var = 0.69^2 = 0.4761
    etalvp        ~ 0.0256   # omega(Pvox) = 0.16 (RSE 30.3%), Table 2 'Variability on volume of ALB-OX in the PC';          var = 0.16^2  = 0.0256
    etalvc_abzson ~ 0.81     # omega(Von)  = 0.90 (RSE 22.4%), Table 2 'Variability on volume of ALB-ON in the CP';          var = 0.90^2  = 0.81

    # ==================================================================
    # Residual error. "The proportional error model best described the
    # residual errors" (Results 3.2). Table 2 labels BOTH residual rows
    # "Sulfoxide, proportional"; the second row is a transcription slip
    # for the sulfone, since the model has exactly two outputs and every
    # other block of Table 2 lists the sulfoxide row before the sulfone
    # row. See the vignette Errata section.
    # ==================================================================
    propSd        <- 0.43 ; label("Albendazole sulfoxide proportional residual error (fraction)")  # Table 2 'Residual variability / Sulfoxide, proportional = 0.43' (RSE 6.18%)
    propSd_Cc_abzson <- 0.22 ; label("Albendazole sulfone proportional residual error (fraction)")    # Table 2 second 'Residual variability' row = 0.22 (RSE 6.75%); row label repeats 'Sulfoxide' but denotes the sulfone
  })

  model({
    # 1. Absorption parameters. Monolix parameterises the Savic (2007)
    #    transit chain by Mtt and Ktr, with the implied (non-integer)
    #    number of transit compartments nn = Mtt * Ktr - 1 and a transit
    #    input that is a Gamma(shape = nn + 1, rate = Ktr) density scaled
    #    by the dose. At the published estimates the shape parameter is
    #    Mtt * Ktr = 0.095 * 0.12 = 0.0114, i.e. nn = -0.9886, just inside
    #    the nn > -1 support boundary. That Gamma(0.0114, 0.12) input is
    #    overwhelmingly front-loaded: 88.5% of the dose has arrived by
    #    0.0001 h (0.36 s), 95.7% by 0.1 h, 98.1% by 1 h - the earliest
    #    observation in either cohort - and 99.7% by 8 h. The chain is
    #    therefore effectively an instantaneous input on the timescale the
    #    data resolve.
    #
    #    The ODE below uses the exact nn -> -1 limiting form (the dose
    #    lands directly in `depot`) for two verified reasons. First,
    #    rxode2's transit() builtin cannot integrate the t^nn singularity
    #    at the published nn: sweeping nn at the published Mtt, the solve
    #    succeeds down to nn = -0.60 and fails at nn <= -0.65, well above
    #    nn = -0.9886. Second, the substitution is quantitatively
    #    negligible: against a mass-conserving discretisation of the exact
    #    Gamma(0.0114, 0.12) input at the published Mtt AND Ktr, the
    #    limiting form differs by at most 2.8% (sulfoxide) and 3.2%
    #    (sulfone) at t >= 1 h, by 0.9% in Cmax, and by 0.07% in
    #    AUC(0-48 h) for both analytes. Absorption into plasma is
    #    rate-limited by the fixed ka = 0.2 1/h either way. The vignette
    #    reruns that comparison as a live check.
    #
    #    `nn` is derived here so the reported Mtt and Ktr stay traceable,
    #    but it does not enter any ODE - which is why lktr and lmtt are
    #    fixed() in ini() above.
    ka  <- exp(lka)
    ktr <- exp(lktr)
    mtt <- exp(lmtt)
    nn  <- ktr * mtt - 1

    # 2. Individual disposition parameters. Study population enters as an
    #    exponential categorical effect (Methods 2.3 Eq. 2) with
    #    Cote d'Ivoire (REGION_TANZANIA = 0) as the reference.
    cl        <- exp(lcl + etalcl) * exp(e_region_tanzania_cl * REGION_TANZANIA)
    vc        <- exp(lvc + etalvc)
    q         <- exp(lq  + etalq)
    vp        <- exp(lvp + etalvp)
    cl_abzson <- exp(lcl_abzson) * exp(e_region_tanzania_cl_abzson * REGION_TANZANIA)
    vc_abzson <- exp(lvc_abzson + etalvc_abzson)

    # 3. ODE system (Fig. 2 structural scheme). The albendazole dose lands
    #    in `depot` (the degenerate transit chain of step 1 delivers it
    #    there instantaneously) and is absorbed at rate ka into the
    #    albendazole sulfoxide central compartment (100% conversion
    #    assumed, Methods 2.2). Sulfoxide has a two-compartment
    #    disposition; everything cleared from the sulfoxide central
    #    compartment becomes albendazole sulfone (again 100% conversion),
    #    which has a one-compartment disposition with linear elimination.
    d/dt(depot)          <- -ka * depot
    d/dt(central)        <- ka * depot -
                              cl * central / vc -
                              q  * central / vc + q * peripheral1 / vp
    d/dt(peripheral1)    <- q  * central / vc - q * peripheral1 / vp
    d/dt(central_abzson) <- cl * central / vc -
                              cl_abzson * central_abzson / vc_abzson

    # 4. Observations. Amount (mg) / apparent volume (L) gives mg/L;
    #    multiply by 1000 to compare with the ng/mL figures and exposure
    #    thresholds quoted in the paper.
    Cc        <- central        / vc
    Cc_abzson <- central_abzson / vc_abzson

    Cc        ~ prop(propSd)
    Cc_abzson ~ prop(propSd_Cc_abzson)
  })
}
