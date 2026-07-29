Jonsson_2011_ethambutol <- function() {
  description <- "Two-compartment population PK model for oral ethambutol in adult South African pulmonary tuberculosis patients (Jonsson 2011), with one transit compartment preceding first-order absorption, allometric scaling on clearance (3/4) and volume (1) terms relative to a 50 kg reference, an HIV-status effect on bioavailability (15.4% reduction), and 4-occasion inter-occasion variability on apparent oral clearance. Parameter values are taken from the publication's Table 2 (NONMEM final estimates column); see inst/modeldb/ddmore/Jonsson_2011_ethambutol_ddmore.R for the DDMoRE-bundle replicate of the same fit."
  reference <- paste(
    "Jonsson S, Davidse A, Wilkins J, Van der Walt JS, Simonsson USH,",
    "Karlsson MO, Smith P, McIlleron H. (2011). Population pharmacokinetics",
    "of ethambutol in South African tuberculosis patients.",
    "Antimicrob Agents Chemother 55(9):4230-4237.",
    "doi:10.1128/AAC.00274-11.",
    "See modellib('Jonsson_2011_ethambutol_ddmore') for the DDMoRE-bundle",
    "(DDMODEL00000220) replicate of the same fit.",
    sep = " "
  )
  vignette <- "Jonsson_2011_ethambutol"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")
  replicate_of <- "inst/modeldb/ddmore/Jonsson_2011_ethambutol_ddmore.R"

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling with theory-based fixed exponents on a 50 kg reference: 3/4 on CL/F and Q/F, 1 on V1/F and V2/F. Jonsson 2011 Results paragraph 'The initial inclusion of allometrically scaled body weight with fixed exponents on all volume and clearance terms resulted in a drop in OFV of 24 units (P < 0.001) and body weight was retained in the model.' Cohort mean WT 47 kg, range 29-86 kg (Table 1).",
      source_name        = "WT"
    ),
    HIV_POS = list(
      description        = "HIV-1 antibody-positive comorbidity indicator (1 = HIV-positive, 0 = HIV-negative).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (HIV-negative)",
      notes              = "13% HIV-positive in the Jonsson 2011 combined cohort (24 of 189 patients per Table 1). Multiplicative shift on bioavailability `f(transit1) <- 1 + e_hiv_pos_f * HIV_POS`; the Table 2 estimate of -0.154 corresponds to a 15.4% reduction in ethambutol bioavailability for HIV-positive subjects (Table 2 footnote a: 'Typical value of F = 1 - 0.154 x HIV, with HIV being 0 and 1 for negative and positive, respectively').",
      source_name        = "HIV"
    ),
    OCC = list(
      description        = "Integer-valued occasion / period indicator for inter-occasion-variability multiplexing.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Values 1, 2, 3, 4 identify the dosing / sampling occasion within subject. The DP Marais SANTA Centre cohort (60 subjects) was sampled on four occasions over a 2-week window (paper Methods 'Patients' paragraph 2: 'Patients were sampled over a 2-week period on four occasions at least 2 weeks after the start of therapy'); the Brewelskloof Hospital cohort (129 subjects) was sampled on a single occasion. Decomposed inside `model()` into binary indicators `oc1` .. `oc4` that multiplex the 4 IOV etas on log-CL. For single-occasion records pass OCC = 1 so the first IOV eta applies.",
      source_name        = "OCC"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 189L,
    n_studies      = 2L,
    age_range      = "16-72 years (median 36)",
    age_median     = "36 years",
    weight_range   = "29-86 kg (median 47)",
    weight_median  = "47 kg",
    sex_female_pct = 46,
    race_ethnicity = c(Black = 16, Coloured = 83, White = 1),
    hiv_positive_pct = 13,
    disease_state  = "Adults with pulmonary tuberculosis pooled across two South African centers (DP Marais SANTA Centre and Brewelskloof Hospital). 13% HIV-positive (24 of 189; HIV is a within-cohort comorbidity rather than the primary indication). 18% had BMIs < 16 kg/m^2 (severe malnutrition); 47% had BMIs >= 18.5 kg/m^2 (not malnourished).",
    dose_range     = "Oral ethambutol 800-1500 mg daily, multiple-dose at steady state, combined with a standard antitubercular backbone (isoniazid, pyrazinamide, rifampin; streptomycin in some). DP Marais SANTA Centre: 5 days/week (drug-free Saturday/Sunday). Brewelskloof Hospital: daily.",
    regions        = "South Africa (two centers near Cape Town and Worcester).",
    notes          = "1,869 plasma ethambutol concentrations across 189 patients analyzed by NONMEM VI (FOCEI with interaction). Estimated baseline creatinine clearance 79 mL/min, range 23-150 (Cockcroft-Gault, truncated at 150); renal function was tested as a covariate on CL/F but not retained in the final model. Demographics from Table 1; final parameter estimates from Table 2."
  )

  ini({
    # Structural PK parameters -- Jonsson 2011 Table 2 'Final parameter estimates by
    # NONMEM, together with bootstrap estimates presented as medians' (NONMEM Estimate
    # column). Reference weight is 50 kg per Table 2 footnote a: 'typical value of
    # CL/F = 39.7 * (body weight/50)^3/4; typical value of central volume of distribution
    # = 82.4 * (body weight/50); ...' (footnote uses bootstrap median 39.7 for the
    # closed-form CL/F formula but the NONMEM final estimate is 39.9).
    lcl    <- log(39.9)  ; label("Apparent oral clearance CL/F at WT = 50 kg (L/h)")               # Table 2 row 'CL/F (liter/h)' NONMEM Estimate = 39.9 (% RSE 3.1)
    lvc    <- log(82.4)  ; label("Apparent central volume of distribution V1/F at WT = 50 kg (L)") # Table 2 row 'V1/F (liters)' NONMEM Estimate = 82.4 (% RSE 42)
    lka    <- log(0.474) ; label("Absorption rate constant from depot to central, ka (1/h)")       # Table 2 row 'ka (h^-1)' NONMEM Estimate = 0.474 (% RSE 24)
    lvp    <- log(623)   ; label("Apparent peripheral volume of distribution V2/F at WT = 50 kg (L)") # Table 2 row 'V2/F (liters)' NONMEM Estimate = 623 (% RSE 22)
    lq     <- log(34.3)  ; label("Apparent inter-compartmental clearance Q/F at WT = 50 kg (L/h)") # Table 2 row 'Q (liters/h)' NONMEM Estimate = 34.3 (% RSE 10)
    lmtt   <- log(0.789) ; label("Mean transit time MTT through the absorption transit compartment (h)") # Table 2 row 'MTT (h)' NONMEM Estimate = 0.789 (% RSE 17)

    # Allometric exponents -- fixed by theory at the canonical values (Anderson and
    # Holford; references 3, 34, 35 of the paper). Methods 'Pharmacokinetic analysis'
    # paragraph 'Initially, allometric scaling by body weight was introduced on all
    # clearance and volume terms with powers of 3/4 and 1, respectively (3, 34, 35).'
    # Discussion 'The effect of body size was modeled by means of allometrically scaled
    # body weight on clearance and volume terms.'
    e_wt_cl_q  <- fixed(0.75) ; label("Allometric exponent on CL/F and Q/F with body weight (unitless)")  # Methods 'Pharmacokinetic analysis' paragraph 'with powers of 3/4 and 1, respectively'
    e_wt_vc_vp <- fixed(1.0)  ; label("Allometric exponent on V1/F and V2/F with body weight (unitless)") # Methods 'Pharmacokinetic analysis' paragraph 'with powers of 3/4 and 1, respectively'

    # HIV covariate effect on bioavailability. Multiplicative fractional change:
    # F = 1 + e_hiv_pos_f * HIV_POS (HIV_POS = 1 gives 15.4% lower bioavailability).
    e_hiv_pos_f <- -0.154 ; label("HIV-positive multiplicative effect on bioavailability (fractional change)") # Table 2 row 'Effect of presence of HIV on F (fractional change)' NONMEM Estimate = -0.154 (% RSE 40); Table 2 footnote a explicit closed form

    # Inter-individual variability (IIV). The paper reports IIV as % CV (Table 2
    # rows 'IIV on ka', 'IIV on MTT', 'IIV on CL/F'). The NONMEM internal variance
    # follows the paper's reported convention CV% ~ sqrt(omega^2) * 100 (i.e. the
    # variance equals (CV/100)^2 to within rounding of the back-transformed CV).
    # The same paper also reports four-occasion IOV on CL/F = 36% CV via the same
    # convention (Methods 'Pharmacokinetic analysis' paragraph: 'Interoccasion
    # variability was evaluated for oral clearance').
    #
    # IIV is included only on ka, MTT, and CL/F (Results 'Population
    # pharmacokinetics' paragraph: 'IIV terms were included on the absorption rate
    # coefficient, CL/F, and the mean transit time, and an IOV term was included on
    # CL/F.'). No IIV on V1, V2, or Q.
    etalka  ~ 0.1521  # Table 2 'IIV on ka (% CV)' = 39 (eta shrinkage 40%); omega^2 = (0.39)^2 = 0.1521
    etalmtt ~ 0.8649  # Table 2 'IIV on MTT (% CV)' = 93 (eta shrinkage 22%); omega^2 = (0.93)^2 = 0.8649
    etalcl  ~ 0.0400  # Table 2 'IIV on CL/F (% CV)' = 20 (eta shrinkage 49%); omega^2 = (0.20)^2 = 0.0400

    # Inter-occasion variability (IOV) on log-CL across four occasions. The paper's
    # IOV form (Methods 'Pharmacokinetic analysis' paragraph): (CL/F)_ij = TV(CL/F) *
    # exp(eta_i + kappa_ij) with eta_i ~ N(0, omega^2) (BSV) and kappa_ij ~ N(0, pi^2)
    # (IOV; same pi^2 across all occasions). nlmixr2 has no NONMEM `$OMEGA BLOCK(1)
    # SAME` shortcut, so each occasion gets its own eta with the variance fixed to
    # the shared value after the first occasion (matching the convention in
    # Jonsson_2011_ethambutol_ddmore.R, Aregbe_2012_alvespimycin.R, and
    # Xie_2019_agomelatine.R).
    etaiov_cl_1 ~ 0.1296         # Table 2 'IOV on CL/F (% CV)' = 36; pi^2 = (0.36)^2 = 0.1296 (estimated)
    etaiov_cl_2 ~ fix(0.1296)    # SAME-equivalent: fixed equal to occasion-1 IOV variance
    etaiov_cl_3 ~ fix(0.1296)    # SAME-equivalent: fixed equal to occasion-1 IOV variance
    etaiov_cl_4 ~ fix(0.1296)    # SAME-equivalent: fixed equal to occasion-1 IOV variance

    # Combined residual error on the linear (mg/L) scale. The paper used log-
    # transformed observations with combined additive + proportional error terms
    # (Methods 'Pharmacokinetic analysis' paragraph: 'Single-compartment and
    # multicompartment models with first-order absorption and elimination were
    # fitted to log-transformed data ... Residual variability was described with
    # both additive and proportional error terms.'). On the back-transformed
    # linear scale this is the combined-error variance var(eps) = (propSd * Cc)^2
    # + addSd^2, i.e. nlmixr2's default `combined2` (Pythagorean SD) form.
    propSd <- 0.318  ; label("Proportional residual error (fraction)") # Table 2 row 'Proportional residual error (% CV)' Estimate = 31.8 (% RSE 4.4); reported on the linear scale
    addSd  <- 0.107  ; label("Additive residual error (mg/L)")         # Table 2 row 'Additive residual error (mg/liter)' Estimate = 0.107 (% RSE 20); reported on the linear scale
  })

  model({
    # Decompose the integer-valued occasion column into binary indicators for IOV
    # multiplexing on log-CL (matches the paper's IOV formulation across 4
    # occasions). For single-occasion subjects (the Brewelskloof Hospital arm)
    # pass OCC = 1 so the first IOV eta applies.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)

    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 + oc4 * etaiov_cl_4

    # Individual PK parameters with allometric scaling on a 50 kg reference weight
    # (theory-based fixed exponents 0.75 on CL/F and Q/F, 1.0 on V1/F and V2/F per
    # Table 2 footnote a's closed-form covariate equations).
    cl  <- exp(lcl  + etalcl  + iov_cl) * (WT / 50)^e_wt_cl_q
    vc  <- exp(lvc)                     * (WT / 50)^e_wt_vc_vp
    ka  <- exp(lka  + etalka)
    vp  <- exp(lvp)                     * (WT / 50)^e_wt_vc_vp
    q   <- exp(lq)                      * (WT / 50)^e_wt_cl_q
    mtt <- exp(lmtt + etalmtt)
    ktr <- 1 / mtt

    # Two-compartment oral PK with one transit compartment preceding first-order
    # absorption. Dose lands in `transit1` (matching the paper's structural
    # description 'one transit compartment prior to first-order absorption'); the
    # transit compartment empties into the absorption compartment `depot` at rate
    # `ktr = 1 / MTT`, which absorbs into `central` at rate `ka`. The central /
    # peripheral block uses elimination CL/Vc and inter-compartmental Q/Vc -> Q/Vp.
    d/dt(transit1)    <- -ktr * transit1
    d/dt(depot)       <-  ktr * transit1 - ka * depot
    d/dt(central)     <-  ka  * depot   - cl / vc * central - q / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <-  q   / vc      * central - q / vp * peripheral1

    # Bioavailability applied to the dosing compartment (transit1). The paper's
    # Table 2 footnote a closed form: F = 1 - 0.154 * HIV_POS; HIV-negative
    # reference keeps F = 1, HIV-positive reduces F by 15.4%.
    f(transit1) <- 1 + e_hiv_pos_f * HIV_POS

    # Concentration in plasma. Dose units mg, Vc units L -> Cc units mg/L (= ug/mL).
    Cc <- central / vc

    # Combined proportional + additive residual error on the linear scale (default
    # Pythagorean / combined2 form). The paper reports the additive component in
    # mg/L (Table 2 row 'Additive residual error') and the proportional component
    # as % CV (Table 2 row 'Proportional residual error').
    Cc ~ add(addSd) + prop(propSd)
  })
}
