Cleary_2023_risdiplam <- function() {
  description <- "Population PK model for risdiplam (Evrysdi) in healthy adults and patients with spinal muscular atrophy aged 2 months to 61 years (Cleary 2023 Clin Pharmacokinet Table 2, final PPK model). Three transit absorption compartments feed a linear two-compartment disposition model. Apparent clearance and intercompartmental clearance scale allometrically with time-varying body weight (estimated exponent 0.276, reference 33 kg); apparent volumes use a separate estimated exponent (0.860). Sigmoidal maturation functions of time-varying postnatal age act on CL/F (Age50 0.877 y) and Vc/F (Age50 0.322 y), and healthy adults carry a higher CL/F than patients with SMA. The proportional residual error switches between venous and capillary blood samples."
  reference <- paste(
    "Cleary Y, Kletzl H, Grimsey P, Heinig K, Ogungbenro K, Silber Baumann HE,",
    "Frey N, Aarons L, Galetin A, Gertz M.",
    "Estimation of FMO3 ontogeny by mechanistic population pharmacokinetic",
    "modelling of risdiplam and its impact on drug-drug interactions in children.",
    "Clin Pharmacokinet. 2023;62(6):891-904.",
    "doi:10.1007/s40262-023-01241-7.",
    sep = " "
  )
  vignette <- "Cleary_2023_risdiplam"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight (time-varying over the observation period; median PK observation period 358 days in patients with SMA).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric power covariate with a median reference weight of 33 kg. Two separate exponents were estimated: 0.276 shared by CL/F and Q/F, and 0.860 shared by Vc/F and Vp/F (Cleary 2023 Table 2 footnote: '[WT/33]^0.276 for CL/F and [WT/33]^0.86 for Vc/F'). Time-varying body weight was used so that growth-driven changes in disposition over the multi-year observation window are captured. The CL/F exponent of 0.276 is far below the theory-based 0.75; the authors interpret this as evidence of higher metabolic activity per gram of liver in children (Discussion, ESM Fig. S25-S26).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Postnatal age (time-varying over the observation period).",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Drives two saturable maturation functions, Age/(Age + Age50), applied to CL/F (Age50 = 0.877 y) and to Vc/F (Age50 = 0.322 y) per the Cleary 2023 Table 2 footnote. Postnatal (not postmenstrual) age was used because most infants in the risdiplam clinical programme were born full term (Sect. 2.3). Time-varying within subject: the analysis used all available time-varying age, body weight and height records. The maturation function is a bare Emax form with no Hill coefficient and no adult normalisation, so it asymptotes to 1 rather than equalling 1 at a reference age.",
      source_name        = "Age"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant cohort indicator (1 = healthy adult, 0 = patient with spinal muscular atrophy).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (patient with spinal muscular atrophy)",
      notes              = "Cleary 2023 Sect. 3.1 introduced 'a factor for CL/F of healthy adults (n = 61)'. All 61 healthy participants in the pooled dataset are adults (26 from NCT02633709 and 35 from NCT03988907); the remaining 464 subjects are patients with SMA types 1, 2 or 3. Table 2 reports the effect only as 'Factor 0.524' and never states how it enters the model, so the encoding is an inferred assumption: CL/F x (1 + 0.524 x DIS_HEALTHY), i.e. healthy adults clear risdiplam 1.524-fold faster than patients with SMA at the same weight and age. The literal multiplicative reading (CL/F x 0.524) is excluded because it would place healthy adults 48% BELOW patients, contradicting ESM Fig. S4 (post-hoc CL/F geometric mean 5.60 L/h healthy vs 3.52 L/h adult SMA), the Discussion ('approximately 30% lower CL/F' in adult SMA) and the Mech-PPK Table 3 adult intrinsic clearances (healthy:SMA ratio 1.374 for CYP3A and 1.375 for FMO3). The linear form was preferred over exp(0.524) = 1.689 because 1.524 sits closer to those in-paper contrasts (1.374-1.59); see the vignette's Assumptions and deviations section. Operator-ratified in sidecar request-002 q2 (option A).",
      source_name        = "HEALTHY"
    ),
    SAMPLE_CAPILLARY = list(
      description        = "Per-observation blood sampling-site indicator (1 = capillary blood sample, 0 = venous blood sample).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (venous sample)",
      notes              = "Record-level indicator that switches the proportional residual-error magnitude per observation: venous 23.4% CV, capillary 34.2% CV (Cleary 2023 Table 2, sigma1 and sigma2). Cleary 2023 Sect. 3.1: 'a separate residual error term for capillary blood samples (3% of the dataset) were introduced.' Both matrices were assayed by the same validated LC-MS/MS method (LLOQ 0.25 ng/mL, ESM 'Bioanalytical method'), so the contrast is the sample matrix / collection site rather than the assay. May vary within a subject. Set SAMPLE_CAPILLARY = 1 on observations drawn as capillary (finger-prick / heel-prick) samples, 0 otherwise.",
      source_name        = "SAMPLE_CAPILLARY"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 525L,
    n_observations = 10205L,
    n_studies      = 5L,
    age_range      = "2 months-61 years",
    age_median     = "7.5 years (model development set), 13.0 years (evaluation/validation set)",
    weight_range   = "4.1-109 kg",
    weight_median  = "20.9 kg (model development set), 33.7 kg (evaluation/validation set); 33 kg allometric reference",
    sex_female_pct = 47,
    disease_state  = "Spinal muscular atrophy types 1, 2 and 3 (n = 464) pooled with healthy adults (n = 61)",
    dose_range     = "0.00106-18 mg oral solution, single or once-daily. Approved regimens: 0.2 mg/kg (2 months to <2 years), 0.25 mg/kg (>=2 years and <20 kg), 5 mg (>=2 years and >=20 kg)",
    regions        = "Multinational (NCT02633709, NCT03032172, NCT02908685, NCT02913482, NCT03988907)",
    notes          = "Pooled analysis of five clinical studies (Cleary 2023 ESM Table S1). Demographics from ESM Table S2: 130 subjects (26 healthy adults, 104 patients with SMA; 2492 observations) formed the model-development set and a further 395 subjects (35 healthy adults, 360 patients with SMA; 7713 observations) formed the evaluation/validation set; the final PPK model reported in Table 2 was fitted to all 525 subjects and 10,205 observations. Sex 278 male : 247 female across both sets. Median PK observation period 358 days in patients with SMA (439 days and up to 3 years in the 382 paediatric patients). Race/ethnicity was not reported. Estimation was FOCE-I in NONMEM 7.4; OFV = 64499; every parameter had RSE < 26% with bootstrap 95% CIs (200 replicates stratified by study, 87.5% converged). Eta-shrinkage was 5.43% (CL/F), 22.9% (ktr) and 10.1% (Vc/F). This file encodes the final PPK model only; the companion mechanistic PPK (Mech-PPK) model of Table 3, which carries the paper's novel in vivo FMO3 ontogeny function, is not yet in nlmixr2lib because it fixes the hepatic CYP3A ontogeny to the closed-access Upreti & Wahlstrom (2016) function, whose coefficients appear in no on-disk source."
  )

  ini({
    # Structural parameters -- Cleary 2023 Table 2 ("Final population PK model
    # parameters of risdiplam in subjects >= 2 months of age"), Fixed effects.
    # Reference subject: WT = 33 kg, fully mature age, patient with SMA
    # (DIS_HEALTHY = 0). All values are apparent (oral, /F) parameters.
    lcl  <- log(2.64);  label("Apparent clearance CL/F (L/h)")                        # Table 2: CL/F = 2.64 L/h (RSE 2.13%, bootstrap 95% CI 2.52-2.76)
    lktr <- log(5.18);  label("Transit absorption rate constant ktr (1/h)")           # Table 2: ktr = 5.18 1/h (RSE 2.74%, bootstrap 95% CI 4.84-5.52)
    lvc  <- log(98.0);  label("Apparent central volume Vc/F (L)")                     # Table 2: Vc/F = 98.0 L (RSE 1.80%, bootstrap 95% CI 93.8-103)
    lq   <- log(0.682); label("Apparent intercompartmental clearance Q/F (L/h)")      # Table 2: Q/F = 0.682 L/h (RSE 10.5%, bootstrap 95% CI 0.589-1.50)
    lvp  <- log(92.9);  label("Apparent peripheral volume Vp/F (L)")                  # Table 2: Vp/F = 92.9 L (RSE 25.8%, bootstrap 95% CI 49.6-133)

    # Allometric exponents -- Cleary 2023 Table 2, Covariate effects. Both were
    # ESTIMATED (not fixed to the theory-based 0.75 / 1.0), so they are not
    # wrapped in fixed(). Reference weight 33 kg (Table 2 footnote).
    e_wt_cl_q  <- 0.276; label("Allometric exponent on CL/F and Q/F (unitless)")      # Table 2: "Effect of WT on CL/F and Q/F" = 0.276 (RSE 11.8%, bootstrap 95% CI 0.167-0.341)
    e_wt_vc_vp <- 0.860; label("Allometric exponent on Vc/F and Vp/F (unitless)")     # Table 2: "Effect of WT on Vc/F and Vp/F" = 0.860 (RSE 3.34%, bootstrap 95% CI 0.792-0.915)

    # Maturation half-ages for the saturable Age/(Age + Age50) functions
    # (Table 2 footnote). Applied to CL/F and to Vc/F only -- Q/F and Vp/F
    # carry allometry but no maturation term.
    age50_cl <- 0.877; label("Age at 50% of mature CL/F (years)")                     # Table 2: "Age50-CL/F" = 0.877 y (RSE 17.1%, bootstrap 95% CI 0.640-1.26)
    age50_vc <- 0.322; label("Age at 50% of mature Vc/F (years)")                     # Table 2: "Age50-Vc/F" = 0.322 y (RSE 21.1%, bootstrap 95% CI 0.226-0.653)

    # Healthy-adult effect on CL/F. Table 2 reports only "Healthy subjects on
    # CL/F | Factor | 0.524" and never states the functional form; the linear
    # form CL/F x (1 + 0.524 x DIS_HEALTHY) is an INFERRED encoding, not a
    # sourced one (operator-ratified, sidecar request-002 q2 option A). See
    # covariateData$DIS_HEALTHY$notes and the vignette's Assumptions and
    # deviations section for the full derivation and rejected alternatives.
    e_dis_healthy_cl <- 0.524; label("Fractional increase in CL/F for healthy adults (unitless)") # Table 2: "Healthy subjects on CL/F" = 0.524 (RSE 13.1%, bootstrap 95% CI 0.392-0.751)

    # IIV -- Cleary 2023 Table 2, Random effects. The paper reports variances
    # with the corresponding CV% in parentheses, and CV% = sqrt(variance)
    # exactly reproduces every printed value: sqrt(0.0678) = 26.0%,
    # sqrt(0.272) = 52.2%, sqrt(0.0651) = 25.5%. So these are variances on the
    # log scale. No correlations between etas were reported.
    etalcl  ~ 0.0678; label("IIV on CL/F (log-scale variance)")                       # Table 2: CL/F (CV) = 0.0678 (26.0%), RSE 8.21%, bootstrap 95% CI 0.0574-0.0810
    etalktr ~ 0.272;  label("IIV on ktr (log-scale variance)")                        # Table 2: ktr (CV) = 0.272 (52.2%), RSE 11.4%, bootstrap 95% CI 0.211-0.329
    etalvc  ~ 0.0651; label("IIV on Vc/F (log-scale variance)")                       # Table 2: Vc/F (CV) = 0.0651 (25.5%), RSE 8.85%, bootstrap 95% CI 0.0535-0.0790

    # Residual error -- Cleary 2023 Table 2, Error model. Two proportional
    # magnitudes were estimated and switched per observation by sampling site.
    # As with the etas, the tabulated numbers are variances and the printed CV%
    # equals their square root: sqrt(0.0546) = 23.4%, sqrt(0.117) = 34.2%.
    propSdVenous    <- sqrt(0.0546); label("Proportional residual error, venous samples (fraction)")    # Table 2: sigma1 proportional-venous (CV) = 0.0546 (23.4%), RSE 3.17%, bootstrap 95% CI 0.0512-0.0575
    propSdCapillary <- sqrt(0.117);  label("Proportional residual error, capillary samples (fraction)") # Table 2: sigma2 proportional-capillary (CV) = 0.117 (34.2%), RSE 15.7%, bootstrap 95% CI 0.0842-0.167
  })

  model({
    # 1. Derived covariate terms
    #
    # Allometric scaling on the median body weight of 33 kg, and the saturable
    # postnatal-age maturation functions (Cleary 2023 Table 2 footnote):
    #   [WT/33]^0.276 for CL/F (and Q/F), [WT/33]^0.86 for Vc/F (and Vp/F)
    #   Age (y)/[Age (y) + 0.877 (y)] for CL/F
    #   Age (y)/[Age (y) + 0.322 (y)] for Vc/F
    # Both maturation functions are bare Emax forms that asymptote to 1 with no
    # adult normalisation, so the Table 2 typical values are the fully mature
    # (adult) values.
    allo_cl_q     <- (WT / 33)^e_wt_cl_q
    allo_vc_vp    <- (WT / 33)^e_wt_vc_vp
    maturation_cl <- AGE / (AGE + age50_cl)
    maturation_vc <- AGE / (AGE + age50_vc)

    # Healthy-adult clearance factor. INFERRED functional form -- see ini()
    # comment and covariateData$DIS_HEALTHY$notes.
    healthy_cl <- 1 + e_dis_healthy_cl * DIS_HEALTHY

    # 2. Individual parameters
    ktr <- exp(lktr + etalktr)
    cl  <- exp(lcl + etalcl) * allo_cl_q  * maturation_cl * healthy_cl
    vc  <- exp(lvc + etalvc) * allo_vc_vp * maturation_vc
    q   <- exp(lq)           * allo_cl_q
    vp  <- exp(lvp)          * allo_vc_vp

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system -- three transit compartments for absorption connected to a
    #    linear two-compartment disposition model (Cleary 2023 Sect. 2.1 and
    #    ESM Sect. 2). The dose enters `depot` and moves through transit1 ->
    #    transit2 -> transit3 -> central, every transfer sharing the single
    #    estimated rate constant ktr, so the mean transit time is 4/ktr.
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2 - ktr * transit3
    d/dt(central)     <-  ktr * transit3 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central  - k21 * peripheral1

    # 5. Observation and error.
    #    Doses are in mg and Vc/F in L, so central/vc is mg/L; the factor of
    #    1000 converts to the ng/mL reported throughout the paper (LLOQ
    #    0.25 ng/mL, calibration range to 250 ng/mL; ESM "Bioanalytical
    #    method"). The proportional residual SD switches per observation
    #    between venous and capillary samples.
    Cc <- 1000 * central / vc
    propSd <- propSdCapillary * SAMPLE_CAPILLARY + propSdVenous * (1 - SAMPLE_CAPILLARY)
    Cc ~ prop(propSd)
  })
}
