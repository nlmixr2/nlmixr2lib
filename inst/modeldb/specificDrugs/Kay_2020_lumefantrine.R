# Population PK model of lumefantrine in HIV-infected and HIV-uninfected
# Ugandan children with uncomplicated malaria receiving artemether-lumefantrine
# alone or with concomitant antiretroviral therapy (efavirenz, lopinavir/ritonavir,
# or nevirapine). Source: Kay K, Goodwin J, Mwebaza N, et al. ASTMH 2020
# poster 2167 (Metrum Research Group / Yale / IDRC).

Kay_2020_lumefantrine <- function() {
  description <- paste(
    "Two-compartment population PK model for oral lumefantrine in 277 HIV-",
    "infected and HIV-uninfected Ugandan children (3 months to ~10 years)",
    "with uncomplicated malaria receiving artemether-lumefantrine alone or",
    "with concomitant ART (efavirenz, lopinavir/ritonavir, or nevirapine)",
    "(Kay 2020, ASTMH poster 2167). First-order absorption (depot ->",
    "central) feeds a 2-compartment lumefantrine disposition (central +",
    "peripheral1). Body-weight allometric scaling enters on all clearance",
    "and volume terms with a fixed volume exponent of 1 and a piecewise",
    "age-dependent clearance exponent (0.75 for age >60 mo, 0.9 for",
    ">24-60 mo, 1.0 for >3-24 mo, 1.2 for <=3 mo). Age also enters as a",
    "covariate on relative bioavailability F (younger children have",
    "reduced F). Three ART drug-drug interactions are encoded as linear-",
    "deviation effects on apparent oral clearance CL/F and on first-order",
    "absorption KA: efavirenz, lopinavir/ritonavir, and nevirapine.",
    "Diagonal IIV is retained on CL/F, V2/F, Q/F, V3/F, and KA. The",
    "NONMEM additive-on-log-scale residual is encoded as a proportional",
    "residual in linear concentration space (consistent with the related",
    "Hoglund 2015 Ugandan-adult lumefantrine model).",
    sep = " "
  )
  reference <- paste(
    "Kay K, Goodwin J, Mwebaza N, Ruiz A, Ehrlich H, Ou J, Freeman T,",
    "Wade M, Huang L, Wang K, Li F, Aweeka FT, Riggs M, Kajubi R, Parikh S",
    "(2020). Modeling and simulation of lumefantrine pharmacokinetics in",
    "HIV-infected and HIV-uninfected children with malaria and the role of",
    "lumefantrine exposure as a potential driver of drug resistance.",
    "American Society of Tropical Medicine and Hygiene (ASTMH) Annual",
    "Meeting, 16-18 November 2020, poster 2167.",
    "Source PDF:",
    "https://metrumrg.com/wp-content/uploads/2020/11/KAY_2167_ASTMH.pdf.",
    sep = " "
  )
  vignette <- "Kay_2020_lumefantrine"
  units    <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight (baseline; treated as time-fixed within an episode)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Body weight enters on all clearance and volume terms via",
        "allometric power scaling (V exponent 1; CL exponent piecewise on",
        "AGE in months). Reference weight WT_REF = 14 kg approximates the",
        "study-population pediatric median used by the source poster; the",
        "reference value is not reported explicitly in Kay 2020 and is",
        "documented as an assumption in the vignette Errata.",
        sep = " "
      ),
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age (used both for the age-dependent allometric exponent on CL/F and as a continuous covariate on relative bioavailability F)",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "AGE is supplied in years (canonical units); the model converts",
        "to months internally (AGE * 12) to evaluate the piecewise",
        "allometric exponent breakpoints stated by Kay 2020 (>60, >24-60,",
        ">3-24, <=3 months). AGE also enters on relative bioavailability F",
        "as a centered linear-deviation effect with reference age AGE_REF",
        "= 5 years. The reference age is not reported explicitly in",
        "Kay 2020 and is documented as an assumption in the vignette",
        "Errata.",
        sep = " "
      ),
      source_name        = "AGE"
    ),
    CONMED_EFV = list(
      description        = "Concomitant efavirenz-based ART indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = HIV-infected child receiving efavirenz (EFV)-based ART at",
        "the time of antimalarial dosing; 0 = no concomitant EFV. The",
        "Kay 2020 cohort had 13 HIV+ children on EFV-based ART (of 38",
        "HIV+ children genotyped for the pfcrt K76 sub-population).",
        "Time-fixed per episode. EFV is a CYP3A4 inducer and enters the",
        "model as a linear-deviation multiplier on apparent oral",
        "lumefantrine clearance (CL/F *= 1 + e_efv_cl * CONMED_EFV with",
        "e_efv_cl = +0.982) and on the first-order absorption rate",
        "constant (KA *= 1 + e_efv_ka * CONMED_EFV with e_efv_ka =",
        "+0.484). The CL/F coefficient corresponds to an approximate 2x",
        "increase in clearance, consistent with the 2.1- to 3.4-fold",
        "EFV-driven reduction in LF exposure reported by Parikh 2016",
        "(reference [1] of the poster).",
        sep = " "
      ),
      source_name        = "EFV"
    ),
    CONMED_LPV = list(
      description        = "Concomitant lopinavir/ritonavir-based ART indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = HIV-infected child receiving lopinavir/ritonavir (LPV/r)-",
        "based ART at the time of antimalarial dosing; 0 = no concomitant",
        "LPV/r. The Kay 2020 cohort had 11 HIV+ children on LPV/r-based",
        "ART (of 38 HIV+ children genotyped for the pfcrt K76",
        "sub-population). Time-fixed per episode. LPV/r is a CYP3A4",
        "inhibitor (ritonavir-driven) and enters the model as a linear-",
        "deviation multiplier on apparent oral lumefantrine clearance",
        "(CL/F *= 1 + e_lpv_cl * CONMED_LPV with e_lpv_cl = -0.514) and",
        "on the first-order absorption rate constant (KA *= 1 + e_lpv_ka",
        "* CONMED_LPV with e_lpv_ka = -0.212). The CL/F coefficient",
        "corresponds to an approximate 2x increase in LF exposure,",
        "consistent with the 2.1-fold LPV/r-driven increase in LF",
        "exposure reported by Parikh 2016 (reference [1] of the poster).",
        sep = " "
      ),
      source_name        = "LPVr"
    ),
    CONMED_NVP = list(
      description        = "Concomitant nevirapine-based ART indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = HIV-infected child receiving nevirapine (NVP)-based ART at",
        "the time of antimalarial dosing; 0 = no concomitant NVP. The",
        "Kay 2020 cohort had 14 HIV+ children on NVP-based ART (of 38",
        "HIV+ children genotyped for the pfcrt K76 sub-population).",
        "Time-fixed per episode. NVP is a CYP3A4 inducer and enters the",
        "model as a linear-deviation multiplier on apparent oral",
        "lumefantrine clearance (CL/F *= 1 + e_nvp_cl * CONMED_NVP with",
        "e_nvp_cl = +0.0191) and on the first-order absorption rate",
        "constant (KA *= 1 + e_nvp_ka * CONMED_NVP with e_nvp_ka =",
        "-0.0589). Both NVP coefficients are small and statistically",
        "non-significant (95% CIs span zero) but the indicators are",
        "retained in the final model as reported in Table 1.",
        sep = " "
      ),
      source_name        = "NVP"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 277L,
    n_studies       = 1L,
    n_episodes      = 364L,
    n_hiv_neg       = 161L,
    n_hiv_pos       = 116L,
    age_range       = "3 months to ~10 years (paediatric; piecewise allometric exponent breakpoints at 3, 24, and 60 months)",
    disease_state   = paste(
      "Uncomplicated malaria; HIV-infected children additionally on daily",
      "ART (EFV, LPV/r, or NVP) and on trimethoprim-sulfamethoxazole (TS)",
      "prophylaxis against opportunistic infections."
    ),
    dose_range      = paste(
      "Artemether-lumefantrine (AL) standard pediatric weight-band",
      "dosing per WHO guidelines, BID for 3 days (6 doses)."
    ),
    regions         = "Uganda (high-malaria-transmission area of eastern Uganda)",
    notes           = paste(
      "Demographics summarized from Kay 2020 ASTMH poster 2167 Methods.",
      "Sub-population genotyped for pfcrt K76 status: 102 HIV- children",
      "(119 episodes) + 38 HIV+ children (57 episodes: 13 EFV, 11 LPV/r,",
      "14 NVP). Companion adult lumefantrine extraction from the same",
      "Aweeka / Kampala lineage: modellib('Hoglund_2015_lumefantrine');",
      "the present model extends Hoglund 2015 to a pediatric cohort with",
      "body-weight + age-on-F covariates and piecewise age-dependent",
      "allometric clearance scaling."
    )
  )

  ini({
    # Structural population mean parameters from Kay 2020 ASTMH poster
    # 2167 Table 1 ("Lumefantrine parameter summary", "Estimate" column).
    # Values are apparent (relative to F = 1) and given on the linear
    # scale; log() is applied here for the nlmixr2 internal log-scale.
    # The reference body weight WT_REF = 14 kg and reference age
    # AGE_REF = 5 yr are not reported in Kay 2020 and are documented in
    # the vignette Errata.

    lcl <- log(1.20)
    label("Apparent lumefantrine clearance CL/F at WT_REF = 14 kg (L/h)")
    # Kay 2020 Table 1: CL/F = theta_1 = 1.20 L/h (95% CI 0.952, 1.45)

    lvc <- log(24.1)
    label("Apparent lumefantrine central volume of distribution V2/F at WT_REF = 14 kg (L)")
    # Kay 2020 Table 1: V2/F = theta_2 = 24.1 L (95% CI 19.0, 29.2)

    lq <- log(0.380)
    label("Apparent lumefantrine inter-compartmental clearance Q/F at WT_REF = 14 kg (L/h)")
    # Kay 2020 Table 1: Q/F = theta_3 = 0.380 L/h (95% CI 0.258, 0.501)

    lvp <- log(767)
    label("Apparent lumefantrine peripheral volume of distribution V3/F at WT_REF = 14 kg (L)")
    # Kay 2020 Table 1: V3/F = theta_4 = 767 L (95% CI 174, 1.36e+03)

    lka <- log(0.0215)
    label("First-order lumefantrine absorption rate constant KA (1/h)")
    # Kay 2020 Table 1: KA = theta_5 = 0.0215 1/h (95% CI 0.0197, 0.0234)

    lfdepot <- fixed(log(1))
    label("Relative bioavailability F (unitless; fixed at 1 at AGE_REF = 5 yr in the no-ART reference)")
    # Anchor: F = 1 by convention in the relative-bioavailability
    # parameterization. The AGE and ART effects are deviations from this
    # anchor (see e_age_f, e_*_cl, e_*_ka entries below). Anchor is fixed
    # in the standard apparent-CL/F lumefantrine convention (cf.
    # Hoglund_2015_lumefantrine.R lfdepot <- fixed(log(1))).

    # Covariate effects on bioavailability F. Linear-deviation form
    # extended from binary to continuous: F = F_anchor * (1 + e_age_f *
    # (AGE/AGE_REF - 1)) with AGE_REF = 5 yr. This continuous-deviation
    # form is the natural extension of the binary linear-deviation form
    # used by the Hoglund 2015 Ugandan-adult lumefantrine model from the
    # same Aweeka / Kampala lineage. The exact functional form is not
    # printed in the Kay 2020 poster (Table 1 lists only the coefficient
    # value); the assumption is documented in the vignette Errata.
    e_age_f <- 0.204
    label("Linear-deviation effect of (AGE/AGE_REF - 1) on relative bioavailability F (fractional)")
    # Kay 2020 Table 1: AGE_F = theta_6 = 0.204 (95% CI -0.0586, 0.467).
    # 95% CI spans zero -> effect retained in model but not significant.

    # ART covariate effects on CL/F: linear-deviation form CL/F =
    # TVCL * (1 + e * IND), consistent with Hoglund 2015 (same group,
    # same drug, Ugandan adults). Table 1 reports each coefficient on
    # the same linear scale as the CL/F apparent value.
    e_efv_cl <- 0.982
    label("Linear-deviation effect of efavirenz on lumefantrine CL/F (fractional)")
    # Kay 2020 Table 1: EFV_CL/F = theta_7 = 0.982 (95% CI 0.163, 1.80)

    e_lpv_cl <- -0.514
    label("Linear-deviation effect of lopinavir/ritonavir on lumefantrine CL/F (fractional)")
    # Kay 2020 Table 1: LPV/r_CL/F = theta_8 = -0.514 (95% CI -0.696, -0.332)

    e_nvp_cl <- 0.0191
    label("Linear-deviation effect of nevirapine on lumefantrine CL/F (fractional)")
    # Kay 2020 Table 1: NVP_CL/F = theta_9 = 0.0191 (95% CI -0.324, 0.362).
    # 95% CI spans zero -> retained but not significant.

    # ART covariate effects on KA: same linear-deviation form as CL/F.
    e_efv_ka <- 0.484
    label("Linear-deviation effect of efavirenz on lumefantrine KA (fractional)")
    # Kay 2020 Table 1: EFV_KA = theta_10 = 0.484 (95% CI 0.282, 0.685)

    e_lpv_ka <- -0.212
    label("Linear-deviation effect of lopinavir/ritonavir on lumefantrine KA (fractional)")
    # Kay 2020 Table 1: LPV/r_KA = theta_11 = -0.212 (95% CI -0.305, -0.120)

    e_nvp_ka <- -0.0589
    label("Linear-deviation effect of nevirapine on lumefantrine KA (fractional)")
    # Kay 2020 Table 1: NVP_KA = theta_12 = -0.0589 (95% CI -0.207, 0.0891).
    # 95% CI spans zero -> retained but not significant.

    # Diagonal IIV. Kay 2020 Table 1 footnote: "CV% of omegas =
    # sqrt(exp(estimate) - 1) * 100", which corresponds to the log-
    # normal CV% recovered from a NONMEM variance on the internal log
    # scale. The "Variance" column values are therefore the log-scale
    # omega^2 directly (e.g. omega^2 = 0.735 -> CV% =
    # sqrt(exp(0.735) - 1) * 100 = 104.2%, matching the "CV%=104"
    # printed in the Estimate column). No off-diagonal covariances were
    # reported, so a diagonal Omega is encoded.
    etalcl ~ 0.735
    # Kay 2020 Table 1: IIV-CL/F = Omega_{1,1} = 0.735 (CV%=104; 95% CI 0.526, 0.945; shrinkage 5.86%)

    etalvc ~ 0.813
    # Kay 2020 Table 1: IIV-V2/F = Omega_{2,2} = 0.813 (CV%=112; 95% CI 0.504, 1.12; shrinkage 32.6%)

    etalq ~ 0.0250
    # Kay 2020 Table 1: IIV-Q/F = Omega_{3,3} = 0.0250 (CV%=15.9; 95% CI 0.0250, 0.0250; shrinkage 67.5%).
    # The 95% CI collapsed to the point estimate suggests the variance
    # was at its lower bound; high shrinkage (67.5%) indicates the data
    # do not strongly inform Q/F IIV.

    etalvp ~ 0.956
    # Kay 2020 Table 1: IIV-V3/F = Omega_{4,4} = 0.956 (CV%=127; 95% CI 0.590, 1.32; shrinkage 36.4%)

    etalka ~ 0.0280
    # Kay 2020 Table 1: IIV-KA = Omega_{5,5} = 0.0280 (CV%=16.9; 95% CI 0.00705, 0.0490; shrinkage 58.2%)

    # Residual variance. Kay 2020 Table 1 footnote: "CV% of sigma =
    # sqrt(estimate) * 100", so the reported Sigma_{1,1} = 0.200 is the
    # variance of a proportional residual in linear concentration space
    # (CV% = sqrt(0.200) * 100 = 44.7%). The propSd entry in nlmixr2
    # parameterizes the SD; SD = sqrt(0.200) = 0.4472.
    propSd <- sqrt(0.200)
    label("Proportional residual SD for lumefantrine plasma concentration (SD on linear scale)")
    # Kay 2020 Table 1: Proportional = Sigma_{1,1} = 0.200 (CV%=44.7; 95% CI 0.173, 0.226)
  })

  model({
    # Reference centering values. WT_REF and AGE_REF are not reported
    # explicitly in Kay 2020 Table 1 and are documented as assumptions
    # in the vignette Errata. WT_REF = 14 kg is the pediatric median
    # weight plausible for a 3-mo-to-10-yr Ugandan cohort (cross-check:
    # back-extrapolation to a 70-kg adult with the >60-month allometric
    # exponent 0.75 yields CL = 1.20 * (70/14)^0.75 = 4.0 L/h, matching
    # published adult lumefantrine CL/F ~4-5 L/h). AGE_REF = 5 yr is the
    # upper boundary of one allometric age bin and is plausible as the
    # study-population median age.
    WT_REF  <- 14
    AGE_REF <- 5

    # Age-dependent allometric exponent on CL/F. Kay 2020 Final Results
    # text: "Fixed effects on volume used an exponent of 1 while the
    # fixed effects on clearance used an exponent of 0.75, 0.9, 1.0 or
    # 1.2 for children aged >60 months, >24 to 60 months, >3 to 24
    # months and 3 months, respectively." Encoded as a sum of indicator
    # products on AGE expressed in months (AGE * 12 to convert from the
    # canonical years).
    AGE_MO <- AGE * 12
    allo_cl_eff <-
      0.75 * (AGE_MO >  60) +
      0.9  * (AGE_MO >  24) * (AGE_MO <= 60) +
      1.0  * (AGE_MO >   3) * (AGE_MO <= 24) +
      1.2  * (AGE_MO <=  3)

    # Individual PK parameters. Body-weight allometric scaling on all
    # clearance and volume terms (V exponent 1.0; CL exponent piecewise
    # via allo_cl_eff). ART DDIs on CL/F and KA as linear-deviation
    # multipliers (one effect per ART; the trial assigned mutually
    # exclusive regimens so at most one of CONMED_EFV / CONMED_LPV /
    # CONMED_NVP is 1 at any time, and the multiplicative form is the
    # canonical encoding for independently estimated categorical
    # effects, mirroring Hoglund_2015_lumefantrine.R).
    cl <- exp(lcl + etalcl) *
            (WT / WT_REF)^allo_cl_eff *
            (1 + e_efv_cl * CONMED_EFV) *
            (1 + e_lpv_cl * CONMED_LPV) *
            (1 + e_nvp_cl * CONMED_NVP)

    vc <- exp(lvc + etalvc) * (WT / WT_REF)

    q  <- exp(lq + etalq)  *
            (WT / WT_REF)^allo_cl_eff

    vp <- exp(lvp + etalvp) * (WT / WT_REF)

    ka <- exp(lka + etalka) *
            (1 + e_efv_ka * CONMED_EFV) *
            (1 + e_lpv_ka * CONMED_LPV) *
            (1 + e_nvp_ka * CONMED_NVP)

    # Typical relative bioavailability with the AGE-on-F covariate
    # effect. Continuous-covariate linear-deviation form anchored at
    # AGE = AGE_REF (centered), bounded to remain positive. Functional
    # form not printed in the source poster; documented as an
    # assumption in the vignette Errata.
    fdepot_typ <- exp(lfdepot) * (1 + e_age_f * (AGE / AGE_REF - 1))

    # Disposition micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system: first-order absorption -> 2-compartment lumefantrine
    # disposition.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot   - kel * central -
                          k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Relative bioavailability applied to the depot.
    f(depot) <- fdepot_typ

    # Plasma concentration in ng/mL. Dose in mg, vc in L gives
    # central / vc in mg/L = ug/mL; multiplying by 1000 yields ng/mL,
    # the unit Kay 2020 reports throughout (Figure 1 y-axis: LF
    # concentration ng/mL).
    Cc <- 1000 * central / vc

    Cc ~ prop(propSd)
  })
}
