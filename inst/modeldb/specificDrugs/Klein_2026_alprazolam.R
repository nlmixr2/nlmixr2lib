Klein_2026_alprazolam <- function() {
  description <- paste(
    "Two-compartment population pharmacokinetic model for inhaled alprazolam delivered by the",
    "Staccato(R) hand-held thermal-aerosol device, pooled across three trials in adults and",
    "adolescents with epilepsy (Klein 2026, N = 99). Absorption is modelled as two parallel",
    "processes whose fractions sum to 1: an immediate-release fraction that enters the central",
    "compartment as a bolus at the dose time (the paper's F2 = 0.369 'fast absorption fraction'),",
    "and the complementary fraction that enters a depot compartment and is absorbed first-order",
    "with ka = 5.37 /h (absorption half-life 7.74 min). Both fractions are supplied from the same",
    "administered dose, so an event table must carry TWO dose records at each administration time,",
    "one to depot and one to central; the f() terms then split the dose. The dose fraction is",
    "carried on the logit scale (the library's logitfdepot idiom) so individual values stay inside",
    "(0, 1) -- the paper's own 95% CI for F2 is exactly symmetric on the logit scale, which proves",
    "the authors estimated it there. All disposition parameters are apparent (/F): the inhaled",
    "bioavailability is not separately identifiable. No relationship between body weight and",
    "clearance or intercompartmental clearance was detectable, so allometric scaling is applied",
    "only to the central and peripheral volumes, through a single freely estimated exponent shared",
    "by both. Clearance carries a multiplicative increase with concomitant strong hepatic",
    "enzyme-inducing antiseizure medications (carbamazepine, phenobarbital, phenytoin, or",
    "primidone). Interindividual variability is log-normal on ka, CL, Vc, Q and Vp (with CL and Vc",
    "correlated) and normal on the logit dose fraction; residual error is proportional."
  )
  reference <- paste(
    "Klein P, Aungaroon G, Biton V, Liow KK, Phillips S, Wychowski T, et al.",
    "Pharmacokinetics and tolerability of single-dose Staccato(R) alprazolam in adolescents",
    "with epilepsy, and population pharmacokinetic analysis to support dose selection in",
    "adolescents. Epilepsia. 2026;67(1):109-119. doi:10.1111/epi.18643.",
    "Population PK parameters are from Supporting Information Table S3."
  )
  vignette <- "Klein_2026_alprazolam"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Scales the central and peripheral volumes through a single freely estimated exponent",
        "(e_wt_vc_vp = 0.470; Table S3 'Allometric scaling Vc and Vp'). Weight was tested on CL",
        "and Q and NO relationship could be detected, so CL and Q are NOT weight-scaled (Methods",
        "section 2.2 and Table S3 title). The paper does not state the normalising reference",
        "weight; 70 kg is assumed here as the rounded standard -- see the vignette Assumptions",
        "section, where the assumption is corroborated against the observed adolescent terminal",
        "half-lives in Table 1."
      ),
      source_name        = "WT"
    ),
    CONMED_EIAED = list(
      description        = "Concomitant strong hepatic enzyme-inducing antiseizure medication",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no inducer antiseizure medication)",
      notes              = paste(
        "Klein 2026 counts carbamazepine, phenobarbital, phenytoin and primidone as strong",
        "hepatic enzyme-inducing antiseizure medications (Methods section 2.2 and Introduction).",
        "15/99 subjects (15.2%) were on one at baseline, all in the adult ENGAGE-E-001 trial;",
        "no adolescent received one (Table S1). Applied on the log scale as",
        "cl = exp(lcl + etalcl + e_conmed_eiaed_cl * CONMED_EIAED), giving a 95.5% increase in",
        "clearance. Patients on strong CYP3A4 inducers were excluded from the adolescent trial",
        "UP0100 by protocol (Methods section 2.1.2)."
      ),
      source_name        = "Inducer ASMs"
    )
  )

  # Covariates screened in the stepwise covariate analysis but NOT retained in the
  # final model (Methods section 2.2: "Covariates at baseline considered in the
  # analysis included height, sex, age, concomitant ASMs, benzodiazepine use,
  # creatinine clearance, and alanine aminotransferase (Table S1)"). No point
  # estimates are reported for any of them, so they are documented here rather
  # than encoded. Chronic benzodiazepine use (27/99 subjects, Table S1) was also
  # screened; it has no canonical column and is recorded in CONMED_AED's notes.
  covariatesDataExcluded <- list(
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened in the covariate analysis; not retained. Overall mean 1.66 m (SD 0.103), Table S1."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened in the covariate analysis; not retained. 70/99 (70.7%) female, Table S1."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened in the covariate analysis; not retained. The Results report the fold change in",
        "alprazolam clearance across 95% of the observed age range as 1.40 (95% CI 0.944-2.09),",
        "i.e. not significantly different from 1. Overall mean age 32.2 years (SD 13.5), Table S1."
      )
    ),
    CONMED_AED = list(
      description = "Concomitant antiseizure medication (any, other than the strong inducers)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened in the covariate analysis; not retained. Only the strong enzyme-inducing subset",
        "(CONMED_EIAED) survived. Chronic benzodiazepine use was screened as a separate covariate",
        "and likewise not retained (27/99, 27.3%, Table S1); it has no canonical covariate column",
        "and is recorded here rather than as its own entry."
      )
    ),
    CRCL = list(
      description = "Body-surface-area-normalized creatinine clearance",
      units       = "mL/min/1.73m2",
      type        = "continuous",
      notes       = "Screened in the covariate analysis; not retained. Overall mean 117 (SD 21.8), Table S1."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened in the covariate analysis; not retained. Overall mean 17.8 IU/L (SD 11.9), Table S1."
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "alprazolam", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "alprazolam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "alprazolam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 99,
    n_studies      = 3,
    age_range      = "12-17 years (adolescents, UP0100); adults in AMDC-002-202 and ENGAGE-E-001",
    age_median     = "32.2 years (mean, SD 13.5) overall; 15.1 years (mean, SD 1.90) in the adolescent trial",
    weight_range   = "33.5-81.2 kg in the adolescent trial UP0100; overall mean 75.5 kg (SD 22.5)",
    weight_median  = "75.5 kg (mean) overall; 54.4 kg (mean) in the adolescent trial",
    sex_female_pct = 70.7,
    race_ethnicity = c(White = 78.6, `Pacific Islander` = 7.1, `Other/mixed` = 14.3),
    disease_state  = "focal, generalized, or focal and generalized epilepsy; photosensitive epilepsy in AMDC-002-202",
    dose_range     = "single inhaled dose of Staccato alprazolam 0.5, 1 or 2 mg",
    regions        = "United States (UP0100); regions not reported for AMDC-002-202 or ENGAGE-E-001",
    co_medication  = "15/99 (15.2%) on a strong enzyme-inducing antiseizure medication; 27/99 (27.3%) chronic benzodiazepine users",
    notes          = paste(
      "Pooled analysis dataset of three trials (Table S1): UP0100 (Phase 1, adolescents,",
      "N = 14, NCT04857307), AMDC-002-202 (Phase 2a, adults with photosensitive epilepsy,",
      "N = 5, NCT02351115) and ENGAGE-E-001 (Phase 2b, adults with epilepsy, N = 80,",
      "NCT03478982). Race percentages are from the adolescent trial only (Table S2); race is",
      "not reported for the pooled dataset. The main text reports N = 84 enrolled in",
      "ENGAGE-E-001 while Table S1 reports N = 80 in the PK modelling dataset."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Structural parameters -- Table S3 (Supporting Information).
    # Every 95% CI in Table S3 is Wald-symmetric on the scale used
    # below, which is what identifies each estimation scale:
    #   ka   exp((log 3.88 + log 7.45)/2) = 5.38  (linear mean 5.67)
    #   CL   exp((log 7.23 + log 9.38)/2) = 8.24  (linear mean 8.31)
    #   Vc   exp((log 23.1 + log 59.2)/2) = 37.0  (linear mean 41.2)
    #   Q    exp((log 131  + log 271 )/2) = 188.5 (linear mean 201)
    #   Vp   exp((log 47.2 + log 77.0)/2) = 60.3  (linear mean 62.1)
    # so all five are estimated on the log scale.
    # ---------------------------------------------------------------
    lka <- log(5.37);  label("Absorption rate constant of the slow (first-order) arm (1/h)")            # Table S3: ka 5.37 (95% CI 3.88-7.45); Results: absorption half-life 7.74 (5.58-10.7) min = log(2)/5.37 h
    lcl <- log(8.23);  label("Apparent clearance CL/F (L/h)")                                           # Table S3: CL/F 8.23 (95% CI 7.23-9.38)
    lvc <- log(36.9);  label("Apparent central volume Vc/F at 70 kg (L)")                               # Table S3: Vc/F 36.9 (95% CI 23.1-59.2)
    lq  <- log(189);   label("Apparent intercompartmental clearance Q/F (L/h)")                         # Table S3: Q/F 189 (95% CI 131-271); not weight-scaled (Table S3 title)
    lvp <- log(60.3);  label("Apparent peripheral volume Vp/F at 70 kg (L)")                            # Table S3: Vp/F 60.3 (95% CI 47.2-77.0)

    # Dose split between the two parallel absorption arms. The paper reports
    # F2 = the FAST (immediate-release into central) fraction; fdepot below is
    # its complement, the fraction routed through the first-order depot arm.
    # F2's 95% CI is exactly symmetric on the LOGIT scale and on no other:
    #   logit: expit((qlogis(0.156) + qlogis(0.649))/2) = 0.3689  <- matches 0.369
    #   linear:      (0.156 + 0.649)/2                  = 0.4025
    #   log:    exp((log 0.156 + log 0.649)/2)          = 0.3182
    # so the fraction was estimated on the logit scale, and the IIV below is a
    # normal deviate on that same scale. logit(1 - p) = -logit(p), so taking the
    # complement leaves the variance unchanged.
    logitfdepot <- qlogis(1 - 0.369); label("Logit of the fraction of the dose entering the first-order depot arm (logit units)")  # Table S3: F2 (fast absorption fraction) 0.369 (95% CI 0.156-0.649), so fdepot = 1 - 0.369 = 0.631

    # ---------------------------------------------------------------
    # Covariate effects
    # ---------------------------------------------------------------
    e_wt_vc_vp <- 0.470; label("Allometric exponent on (WT/70) shared by Vc and Vp (unitless)")         # Table S3: 'Allometric scaling Vc and Vp' 0.470 (95% CI 0.341-0.598); CI is symmetric on the natural scale ((0.341 + 0.598)/2 = 0.4695), so this is a plain THETA

    # Reported in Table S3 as a "% change" of 95.5 (95% CI 50.5-154.0). That CI
    # is symmetric on the log scale and not on the linear one:
    #   (log 1.505 + log 2.540)/2 = 0.67048 = log(1.955)   -> 95.5% increase
    #   (50.5 + 154.0)/2          = 102.25                 -> would be 102.3%
    # so the effect multiplies CL by exp(coefficient), as encoded in model().
    e_conmed_eiaed_cl <- log(1.955); label("Log fold-change in CL/F with a concomitant enzyme-inducing antiseizure medication (unitless)")  # Table S3: effect of inducer ASMs on CL, +95.5% (95% CI 50.5-154.0%)

    # ---------------------------------------------------------------
    # Interindividual variability -- Table S3 "IIV" column.
    # The IIV column is the omega STANDARD DEVIATION on each parameter's own
    # estimation scale, expressed as a percentage, NOT the exact log-normal
    # CV% (sqrt(exp(omega^2) - 1)). Figure 4 settles this: the supplement
    # states AUCinf was computed as dose / CL, so the simulated AUCinf
    # distribution is exactly log-normal with median 2000/8.23 = 243 ng*h/mL,
    # and its published 5th-95th band can only match one reading --
    #   omega = 0.586 (SD reading):     93 - 637 ng*h/mL
    #   omega = 0.543 (exact CV% read): 100 - 593 ng*h/mL
    #   Figure 4, no-inducer panel:    ~93 - ~645 ng*h/mL
    # and again on the inducer panel (median 2000/(8.23*1.955) = 124):
    #   SD reading 47-326, CV% reading 51-303, figure ~45 - ~330.
    # The vignette scores both readings against the digitised figure.
    # ---------------------------------------------------------------
    etalka ~ 0.730^2                                                          # Table S3: IIV on ka 73.0% (shrinkage 26.0%)
    # Table S3: IIV on CL/F 58.6% (shrinkage 5.3%) and on Vc/F 76.3%
    # (shrinkage 8.7%); footnote "Correlation between CL and Vc: 0.342".
    # (The source trace sits ABOVE the block: a trailing comment on a
    # multi-line omega c(...) trips an rxode2 comment-to-label parse bug.)
    etalcl + etalvc ~ c(0.586^2,
                        0.342 * 0.586 * 0.763, 0.763^2)
    etalq  ~ 0.396^2                                                          # Table S3: IIV on Q/F 39.6% (shrinkage 44.7%)
    etalvp ~ 0.209^2                                                          # Table S3: IIV on Vp/F 20.9% (shrinkage 37.9%)
    etalogitfdepot ~ 1.078^2                                                  # Table S3: IIV on F2 107.8% (shrinkage 51.6%); a normal deviate on the logit scale, per the CI-symmetry evidence above

    # ---------------------------------------------------------------
    # Residual error
    # ---------------------------------------------------------------
    propSd <- 0.100; label("Proportional residual error (fraction)")          # Table S3: proportional RUV 10.0% (95% CI 8.30-11.7), shrinkage 30.3%
  })

  model({
    # 1. Individual parameters. CL and Q are deliberately NOT weight-scaled:
    #    the paper found no detectable weight relationship for either.
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl + e_conmed_eiaed_cl * CONMED_EIAED)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp

    # Fraction of the dose taking the slow, first-order route through the depot.
    fdepot <- expit(logitfdepot + etalogitfdepot)

    # 2. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. ODE system
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 4. Dose split across the two parallel absorption arms. Both arms are fed
    #    from the same administered dose, so the event table must carry two dose
    #    records per administration -- one to depot, one to central -- and these
    #    f() terms apportion them. The central-compartment record is the
    #    immediate-release ("fast") fraction the paper reports as F2.
    f(depot)   <- fdepot
    f(central) <- 1 - fdepot

    # 5. Observation. Dose is in mg and vc in L, so central/vc is mg/L; multiply
    #    by 1000 to report ng/mL as the paper does.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
