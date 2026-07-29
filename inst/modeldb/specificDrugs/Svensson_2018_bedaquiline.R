Svensson_2018_bedaquiline <- function() {
  description <- "Three-compartment population PK model for the antimycobacterial bedaquiline (BDQ) and a two-compartment N-desmethyl metabolite M2 in healthy adult volunteers following single 400 mg oral doses, with four-transit-compartment first-order absorption (rate of absorption from the last transit compartment fixed equal to the inter-transit transfer rate, i.e. KA = KTR) and a multiplicative formulation effect adding 23% to the typical mean absorption time when the four 100 mg tablets are suspended in water before swallowing relative to swallowing the tablets whole."
  reference <- paste(
    "Svensson E. M., du Bois J., Kitshoff R., de Jager V. R., Wiesner L.,",
    "Norman J., Nachman S., Smith B., Diacon A. H., Hesseling A. C.,",
    "Garcia-Prats A. J. (2018).",
    "Relative bioavailability of bedaquiline tablets suspended in water:",
    "Implications for dosing in children.",
    "British Journal of Clinical Pharmacology 84(10):2384-2392.",
    "doi:10.1111/bcp.13696.",
    "Structural PK starting point: Svensson 2016 (DDMODEL00000219;",
    "see modellib('Svensson_2016_bedaquiline')).",
    sep = " "
  )
  vignette <- "Svensson_2018_bedaquiline"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    FORM_SUSPENSION = list(
      description        = "Suspended-tablet formulation indicator (1 = four 100 mg bedaquiline tablets suspended in approximately 30 mL water at bedside immediately before swallowing; 0 = four 100 mg tablets swallowed whole with water)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (tablets swallowed whole; the typical-value MAT reference in Svensson 2018 Table 2)",
      notes              = "Per-dose-occasion (not per-subject) formulation indicator. Multiplicative effect on the typical mean absorption time MAT relative to the whole-tablet reference: MAT_susp = MAT_whole * (1 + e_susp_mat * FORM_SUSPENSION) with the +23% point estimate from Svensson 2018 Table 2 ('Effect of suspending on MAT (%) = 23'). The same paper found no statistically significant difference in relative bioavailability (95% nonparametric CI 94-108%, predefined bioequivalence criteria 80-125%), so F = 1 is held identically across formulations; only MAT shifts. Sibling formulation indicator under the FORM_* family alongside FORM_TABLET (Kyhl 2016 / Tikiso 2021 tablet vs liquid solution), FORM_CAPSULE, FORM_POWDER, and the various drug-product-version indicators (FORM_SAR_DP2, FORM_ISA_P2F2, FORM_LINAG_TAB1, FORM_VISMO_PHASEI).",
      source_name        = "FORM"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 24L,
    n_studies      = 1L,
    age_range      = "19-37 years",
    age_median     = "23.5 years",
    weight_range   = "45.6-88.5 kg",
    weight_median  = "63.4 kg",
    sex_female_pct = 62.5,
    race_ethnicity = c(Black = 87.5, MixedRace = 12.5),
    disease_state  = "Healthy adult volunteers (no clinical evidence of QT prolongation, dysrhythmia, significant cardiac conditions, liver/kidney disease, HIV infection, hepatitis B/C, hypothyroidism, suspected or documented active TB, or recent household TB exposure; on no QT-prolonging medications and no CYP3A4 inducers or inhibitors).",
    dose_range     = "Single 400 mg oral dose of bedaquiline administered twice 14 days apart (open-label two-period crossover); on each occasion as either four 100 mg tablets swallowed whole with 240 mL water, or four 100 mg tablets suspended in 30 mL water with two further rinses (20 mL + 10 mL). Both doses given within 30 min after a standardized 670 kcal breakfast with at least 33% fat content.",
    regions        = "Cape Town, South Africa (single site, November-December 2016).",
    notes          = "Baseline demographics from Svensson 2018 Table 1. Twelve participants were randomized to whole-tablets-first then suspension at the second occasion; twelve to suspension first then whole tablets. All 24 participants completed the study. There were 552 concentration observations each for bedaquiline and M2; below-LLOQ observations (LLOQ 0.01 ug/mL for BDQ and M2) were excluded from the analysis (5 postdose BDQ and 81 postdose M2 observations were below LLOQ; 81 of the 81 M2 BLOQ samples occurred within 4 h of dose administration). The pediatric clinical motivation of this study (using suspended tablets when a pediatric formulation is not yet available) is not reflected in the cohort itself, which is healthy adult volunteers."
  )

  ini({
    # Structural absorption: four-transit-compartment first-order chain with
    # the rate of absorption from the last transit compartment fixed equal
    # to the inter-transit transfer rate KTR (Svensson 2018 Results paragraph
    # 1 of the Pharmacokinetic analysis subsection: 'The absorption model was
    # simplified without a statistically significant loss of fit by making
    # the rate of absorption from the last transit compartment the same as
    # the rate of transfer between the transit compartments.'). Under this
    # simplification, the typical mean absorption time MAT relates to KTR by
    # MAT = (NN + 1) / KTR, where NN is the number of transit compartments.
    lmat   <- log(2.63)       ; label("Mean absorption time MAT (h, on log scale)")  # Svensson 2018 Table 2 MAT = 2.63 h (RSE 5.0%)
    nn_fix <- fixed(4)        ; label("Number of transit compartments (Savic-style integer, unitless)") # Svensson 2018 Table 2 NN = 4.00 (RSE 10.9%); held at the integer value 4 to allow an explicit four-transit ODE chain

    # Structural disposition: 3-compartment bedaquiline (parent) +
    # 2-compartment M2 metabolite. Apparent volumes and clearances (i.e.,
    # CL/F, V/F for the parent and CL/(F*fm), V/(F*fm) for the metabolite,
    # where fm is the fraction of bedaquiline metabolised to M2 -- not
    # separately identifiable from the F-relative parameters).
    lcl     <- log(5.67)      ; label("Apparent bedaquiline clearance CL/F (L/h)")                        # Svensson 2018 Table 2 'CL_BDQ/F = 5.67 L/h' (RSE 10.1%)
    lvc     <- log(130)       ; label("Apparent bedaquiline central volume Vc/F (L)")                     # Svensson 2018 Table 2 'V_BDQ/F = 130 L' (RSE 6.1%)
    lq      <- log(6.33)      ; label("Apparent bedaquiline first-peripheral inter-compartmental clearance Q1/F (L/h)") # Svensson 2018 Table 2 'Q_BDQ,1/F = 6.33 L/h' (RSE 9.6%)
    lvp     <- log(3020)      ; label("Apparent bedaquiline first-peripheral volume Vp1/F (L)")            # Svensson 2018 Table 2 'VP_BDQ,1/F = 3020 L' (RSE 28.0%)
    lq2     <- log(4.83)      ; label("Apparent bedaquiline second-peripheral inter-compartmental clearance Q2/F (L/h)") # Svensson 2018 Table 2 'Q_BDQ,2/F = 4.83 L/h' (RSE 15.5%)
    lvp2    <- log(64.5)      ; label("Apparent bedaquiline second-peripheral volume Vp2/F (L)")           # Svensson 2018 Table 2 'VP_BDQ,2/F = 64.5 L' (RSE 13.1%)
    lcl_m2  <- log(17.2)      ; label("Apparent M2 metabolite clearance CLM2/(F*fm) (L/h)")                # Svensson 2018 Table 2 'CL_M2/F/fm = 17.2 L/h' (RSE 11.8%)
    lvc_m2  <- log(1380)      ; label("Apparent M2 metabolite central volume Vc_m2/(F*fm) (L)")            # Svensson 2018 Table 2 'V_M2/F/fm = 1380 L' (RSE 9.2%)
    lq_m2   <- log(126)       ; label("Apparent M2 metabolite inter-compartmental clearance Q_m2/(F*fm) (L/h)") # Svensson 2018 Table 2 'Q_M2/F/fm = 126 L/h' (RSE 12.9%)
    lvp_m2  <- log(3450)      ; label("Apparent M2 metabolite peripheral volume Vp_m2/(F*fm) (L)")         # Svensson 2018 Table 2 'VP_M2/F/fm = 3450 L' (RSE 11.7%)

    # Bioavailability is fixed to F = 1 because CL and V are reported as
    # apparent F-relative values (CL/F, V/F). The relative-bioavailability
    # study found no statistically significant difference between the two
    # formulations (95% nonparametric CI of relative F 94-108%, predefined
    # BE criteria 80-125%), so F is identical across whole and suspended
    # tablets; only MAT shifts. See vignette's Assumptions and deviations.
    lfdepot <- fixed(log(1))  ; label("Bioavailability F (fixed at 1 because CL and V are apparent F-relative values; relative F whole vs suspension was bioequivalent within 80-125% per Svensson 2018 Results)")

    # Formulation covariate effect. Encoded as additive-fractional on MAT:
    # MAT_typ = MAT_whole * (1 + e_susp_mat * FORM_SUSPENSION). At
    # FORM_SUSPENSION = 1 the typical MAT is 23% longer than at
    # FORM_SUSPENSION = 0.
    e_susp_mat <- 0.23        ; label("Effect of suspending tablets on typical MAT (fraction; +23% relative to whole tablets)") # Svensson 2018 Table 2 'Effect of suspending on MAT (%) = 23' (RSE 43.0%; 95% CI 2.1-48%, P = 0.03)

    # Inter-individual variability (between-subject, BSV). Final estimates
    # converted from the paper's CV% notation via omega^2 = log(1 + CV^2)
    # for log-normal IIV. The full paper reports an additional set of
    # between-occasion (BOV) terms on F (9.1% CV) and MAT (66.3% CV); the
    # BOV-F term is dropped here as nlmixr2lib has no idiomatic encoding
    # for between-occasion variability separate from BSV, and the BOV-MAT
    # term is folded in as the BSV-MAT-equivalent below (paper reports no
    # separate BSV on MAT). See vignette's Assumptions and deviations.
    etalcl + etalcl_m2 ~ c(0.02891, 0.002935, 0.04122)  # Svensson 2018 Table 2 IIV CL_BDQ = 17.1% CV, IIV CL_M2 = 20.5% CV, correlation 8.5%; omega^2_cl = log(1+0.171^2) = 0.02891, omega^2_clm2 = log(1+0.205^2) = 0.04122, cov = 0.085 * sqrt(0.02891 * 0.04122) = 0.002935
    etalvc      ~ 0.07682   # Svensson 2018 Table 2 IIV V_BDQ = 28.3% CV; omega^2 = log(1 + 0.283^2) = 0.07682
    etalq       ~ 0.02957   # Svensson 2018 Table 2 IIV Q_BDQ,1 = 17.3% CV; omega^2 = log(1 + 0.173^2) = 0.02957
    etalvc_m2   ~ 0.06677   # Svensson 2018 Table 2 IIV V_M2 = 26.3% CV; omega^2 = log(1 + 0.263^2) = 0.06677
    etalvp_m2   ~ 0.04863   # Svensson 2018 Table 2 IIV VP_M2 = 22.3% CV; omega^2 = log(1 + 0.223^2) = 0.04863
    etalfdepot  ~ 0.04983   # Svensson 2018 Table 2 IIV F = 22.6% CV; omega^2 = log(1 + 0.226^2) = 0.04983 (IOV F = 9.1% CV is dropped; see vignette Assumptions and deviations)
    etalmat     ~ 0.36764   # Svensson 2018 Table 2 IOV MAT = 66.3% CV folded in as BSV-MAT equivalent (no separate BSV on MAT reported); omega^2 = log(1 + 0.663^2) = 0.36764

    # Residual error.  The paper reports proportional residual errors on
    # both bedaquiline and M2 concentrations and a correlation between the
    # BDQ and M2 residuals (53.1% per Table 2; RSE 11.7%); the correlation
    # is dropped here because nlmixr2lib has no idiomatic encoding for
    # cross-output residual correlation. The paper also reports a 1.67-fold
    # weighting of the residual error during the first 6 h after dose
    # ('Weighting residual error samples 0-6 h = 1.67' in Table 2,
    # introduced to absorb the unmodelled dual-peak behaviour described in
    # the Discussion); this time-varying multiplier on the residual SD is
    # dropped here for similar lack of idiomatic encoding. See vignette
    # Assumptions and deviations for both deviations.
    propSd     <- 0.231       ; label("Bedaquiline residual error (proportional SD as a fraction)") # Svensson 2018 Table 2 'Proportional error BDQ = 23.1% CV'
    propSd_m2  <- 0.114       ; label("M2 metabolite residual error (proportional SD as a fraction)") # Svensson 2018 Table 2 'Proportional error M2 = 11.4% CV'
  })

  model({
    # 1. Derived terms.  Formulation effect on the typical MAT.  IIV on MAT
    #    (etalmat) is applied around the formulation-adjusted typical.
    mat_form <- 1 + e_susp_mat * FORM_SUSPENSION

    # 2. Individual PK parameters.
    mat   <- exp(lmat + etalmat) * mat_form
    nn    <- nn_fix
    ktr   <- (nn + 1) / mat   # transit-chain rate (1/h); under KA = KTR with NN+1 first-order steps from depot to central, MAT = (NN+1)/KTR

    cl    <- exp(lcl   + etalcl)
    vc    <- exp(lvc   + etalvc)
    q     <- exp(lq    + etalq)
    vp    <- exp(lvp)
    q2    <- exp(lq2)
    vp2   <- exp(lvp2)
    cl_m2 <- exp(lcl_m2 + etalcl_m2)
    vc_m2 <- exp(lvc_m2 + etalvc_m2)
    q_m2  <- exp(lq_m2)
    vp_m2 <- exp(lvp_m2 + etalvp_m2)

    # 3. ODE system.
    #    Four explicit transit compartments (NN = 4) with shared rate KTR
    #    on depot-to-transit1, transit1-to-transit2, transit2-to-transit3,
    #    transit3-to-transit4, and transit4-to-central (the KA = KTR
    #    simplification of Svensson 2018 Pharmacokinetic-analysis paragraph).
    #    Three-compartment bedaquiline disposition (central + peripheral1
    #    + peripheral2) with linear elimination feeding the metabolite
    #    central_m2 compartment; M2 has its own two-compartment disposition
    #    (central_m2 + peripheral1_m2) and is eliminated from central_m2.
    d/dt(depot)            <- -ktr * depot
    d/dt(transit1)         <-  ktr * depot     - ktr * transit1
    d/dt(transit2)         <-  ktr * transit1  - ktr * transit2
    d/dt(transit3)         <-  ktr * transit2  - ktr * transit3
    d/dt(transit4)         <-  ktr * transit3  - ktr * transit4
    d/dt(central)          <-  ktr * transit4 -
                               cl  * central / vc -
                               q   * central / vc + q  * peripheral1 / vp -
                               q2  * central / vc + q2 * peripheral2 / vp2
    d/dt(peripheral1)      <-  q   * central / vc - q  * peripheral1 / vp
    d/dt(peripheral2)      <-  q2  * central / vc - q2 * peripheral2 / vp2
    d/dt(central_m2)       <-  cl  * central / vc -
                               cl_m2 * central_m2 / vc_m2 -
                               q_m2  * central_m2 / vc_m2 + q_m2 * peripheral1_m2 / vp_m2
    d/dt(peripheral1_m2)   <-  q_m2  * central_m2 / vc_m2 - q_m2 * peripheral1_m2 / vp_m2

    # 4. Bioavailability.  F is fixed at 1 with etalfdepot carrying BSV on F
    #    (paper reports BSV F = 22.6% CV in Table 2; BOV F = 9.1% CV is
    #    dropped per the ini() comment above).
    f(depot) <- exp(lfdepot + etalfdepot)

    # 5. Observed plasma concentrations.  Dose mg / volume L = mg/L = ug/mL.
    Cc    <- central    / vc
    Cc_m2 <- central_m2 / vc_m2

    Cc    ~ prop(propSd)
    Cc_m2 ~ prop(propSd_m2)
  })
}
