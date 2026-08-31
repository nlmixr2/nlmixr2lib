Kamal_2014_oseltamivir <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for oral oseltamivir",
    "(prodrug, OP) and its active metabolite oseltamivir carboxylate",
    "(OC) in 133 infants aged 2 weeks through 11 months with confirmed",
    "influenza, pooled from the CASG114 (NCT00391768, USA) and WP22849",
    "(NCT01286142, Europe) prospective PK/PD studies. Oseltamivir is",
    "described by a two-compartment model with first-order absorption",
    "and first-order conversion to OC; OC is described by a",
    "one-compartment model with first-order elimination (three",
    "disposition compartments in total). Complete (100%) conversion of",
    "oseltamivir to OC was assumed because the OP-to-OC conversion",
    "fraction and the OC volume are not simultaneously identifiable",
    "without OC administration data, so the parent elimination clearance",
    "CL/F is the OP-to-OC formation clearance and all OC terms are",
    "apparent (conditioned on oral bioavailability F and on fm = 1).",
    "Covariates: body weight on all three clearance terms (CL/F, Q/F,",
    "CLM/F) and all three volume terms (V2/F, V3/F, VM/F) by allometry",
    "with exponents fixed mechanistically at 0.75 and 1 and a reference",
    "weight of 8 kg; and postnatal age as a linear (not power) effect on",
    "both OC parameters CLM/F and VM/F, referenced to 24 weeks, using a",
    "single coefficient constrained equal for the two parameters",
    "(VM,AGE = CLM,AGE), corresponding to an increase of about 72% per",
    "year of age. Inter-individual variability is exponential on CL/F,",
    "V2/F, and CLM/F as a fully correlated 3x3 block; ka, Q/F, V3/F, and",
    "VM/F carry no random effect in the final model. Residual error is",
    "proportional only for oseltamivir (46.7% CV) and combined",
    "proportional plus additive for OC (11.8% CV + 39.5 ng/mL). The",
    "model supported the US Food and Drug Administration approval of a",
    "3 mg/kg twice-daily oseltamivir dose for infants aged 2 weeks",
    "through 11 months.")
  reference <- paste(
    "Kamal MA, Acosta EP, Kimberlin DW, Gibiansky L, Jester P,",
    "Niranjan V, Rath B, Clinch B, Sanchez PJ, Ampofo K, Whitley R,",
    "Rayner CR. The posology of oseltamivir in infants with influenza",
    "infection using a population pharmacokinetic approach.",
    "Clin Pharmacol Ther. 2014;96(3):380-389.",
    "doi:10.1038/clpt.2014.120. PMID 24865390.",
    sep = " ")
  vignette <- "Kamal_2014_oseltamivir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Both analytes were measured in plasma by a
  # validated HPLC/MS-MS method (Kamal 2014 Methods, first paragraph and
  # "Limits of quantitation for oseltamivir (parent drug) and OC were 1
  # and 50 ng/ml"), so analyte and specimen are confirmed against the
  # source.
  compartmentData <- list(
    depot            = list(analyte = "oseltamivir", units = "mg", specimen = "administration site", verified = TRUE),
    central          = list(analyte = "oseltamivir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1      = list(analyte = "oseltamivir", units = "mg", specimen = "plasma", verified = TRUE),
    central_oselcarb = list(analyte = "oseltamivir carboxylate", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight (baseline).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling of every clearance term (CL/F, Q/F, CLM/F) as (WT/8)^0.75 and every volume term (V2/F, V3/F, VM/F) as (WT/8), per the Kamal 2014 Table 3 footnote (page 383) and the Figure 2 caption. Both exponents were FIXED mechanistically rather than estimated: Kamal 2014 Methods, 'Allometric scaling and covariate model development' -- 'Clearance and volume parameters were scaled allometrically with fixed power coefficients of 0.75 and 1, respectively, using body weight'; and 'Body weight and age were strongly correlated, and their effects on model parameters could not be estimated simultaneously. Body weight dependencies were therefore fixed mechanistically, while the remaining parameter dependencies on body size and age were explained by age effects.' The 8 kg reference is stated twice: Table 3 footnote ('(WT/8)0.75 and (WT/8)') and the Figure 2 caption ('CLM*(WT/8)0.75, where WT is infant body weight in kilograms and the median body weight is 8 kg'). Note the 8 kg reference is NOT the cohort mean, which Table 2 gives as 6.5 kg (SD 2.1, range 2.9-12.4); see the vignette Assumptions and deviations. Kamal 2014 Discussion page 384 notes 'Scaling to an adult weight (70 kg) instead of the median weight would not have altered the model because the power of the weight dependence was fixed', so the reference weight is a re-centering constant only.",
      source_name        = "WT"
    ),
    PNA = list(
      description        = "Postnatal age (chronological age since birth) at baseline.",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear (additive-in-the-multiplier, NOT power) effect on the two OC parameters: CLM/F ~ 1 + CLM,AGE*(AGE/24 - 1) and VM/F ~ 1 + VM,AGE*(AGE/24 - 1) with VM,AGE = CLM,AGE, per the Kamal 2014 Table 3 footnote (page 383). A SINGLE estimated coefficient (theta11 = 0.33) therefore drives both OC clearance and OC volume; the equality is a structural constraint of the published model, not two coincidentally equal estimates. The source column is the infant's chronological age in WEEKS with a reference of 24 weeks; the canonical PNA carries months, so the reference is reparameterised inside model() as ref_pna = 24 * 7 / 30.4375 = 5.5195 months (the same days-to-months reparameterisation the register documents for Zhao_2018_omeprazole). Mapped to canonical PNA rather than canonical AGE because Kamal 2014 explicitly distinguishes 'age' from 'postconceptual age' in its covariate list (Methods, 'Allometric scaling and covariate model development': 'Covariates were weight, age or postconceptual age, study, and sex'), which is exactly the register's PNA (postnatal) versus PAGE (postmenstrual) distinction; and because the cohort spans only 1.9-49.9 weeks, where a years scale loses resolution. Kamal 2014 Results page 381 reports the effect as OC clearance and volume 'increasing linearly with age (~72% per year)', which reproduces from 0.33*(52.1775/24) = 0.717. Baseline-only here: the analysis used a single baseline age per subject over a 5- to 10-day treatment course. Kamal 2014 Discussion page 384 cautions that 'This dependence should not be extrapolated beyond the range of the analysis data (<1 year).'",
      source_name        = "AGE"
    )
  )

  # Covariates that Kamal 2014 screened in its full covariate model but
  # did NOT retain in the final model. Documented here (rather than in
  # covariateData) so the provenance of the covariate screen is preserved
  # without declaring covariates that model() never references.
  covariatesDataExcluded <- list(
    PAGE = list(
      description        = "Postmenstrual (the paper's 'postconceptual') age.",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as an alternative to postnatal age and not retained: Kamal 2014 Results page 381 -- 'Use of postconceptual age instead of age did not improve the model (it should be noted that infants who were premature but of postconceptual age >=36 weeks were included in both studies, but no infants were diagnosed or started on treatment until after attainment of ages beyond full term).' Table 2 reports it in WEEKS: mean 61.6 (SD 15.3), range 38.4-90.0 across all 133 subjects. No point estimate is available because the effect was rejected.",
      source_name        = "PCA"
    ),
    GA = list(
      description        = "Gestational age at birth.",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reported as a baseline descriptor in Kamal 2014 Table 2 (all subjects: mean 38.1 weeks, SD 4.5, range 24.0-43.0) and used to derive postconceptual age, but not itself listed among the screened covariates in Methods, 'Allometric scaling and covariate model development'. Retained here for provenance only. Preterm neonates of postmenstrual age <36 weeks were excluded from both trials (Kamal 2014 Discussion page 385).",
      source_name        = "GA"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Screened in the full covariate model (Kamal 2014 Methods, 'Allometric scaling and covariate model development': 'Covariates were weight, age or postconceptual age, study, and sex') and not retained. Results page 381 states the final model contained exactly 'two significant covariate effects, with OC apparent clearance and volume increasing linearly with age'. Cohort split per Table 2: 59 female (44.4%), 74 male (55.6%). No point estimate is available because the effect was rejected.",
      source_name        = "SEX"
    ),
    STUDY_CASG114 = list(
      description        = "Indicator for the CASG114 study arm (1 = CASG114, 0 = WP22849).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (WP22849)",
      notes              = "Screened as a relative-bioavailability (F1) effect to test for a formulation difference between the two trials -- CASG114 used a proprietary 12 mg/mL powder-for-oral-suspension formulation and WP22849 used a 10 mg/mL suspension compounded from 75 mg capsules (Kamal 2014 Methods, first paragraph). Not retained: Results page 381 -- 'Addition of the effect of study CASG114 to F1 (relative bioavailability) showed no evidence that the differing formulations used in the two studies (CASG114 and WP22849) had any effect on the pharmacokinetics of the parent drug or its metabolite.' Cohort split per Table 2: CASG114 68 subjects (51.1%), WP22849 65 subjects (48.9%). No point estimate is available because the effect was rejected; relative bioavailability is therefore 1 for both studies and no F term appears in model().",
      source_name        = "STUDY"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 133L,
    n_studies      = 2L,
    age_range      = "1.9-49.9 weeks postnatal age (0-12 months); label indication is 2 weeks through 11 months",
    age_mean       = "23.5 weeks (SD 14.9)",
    weight_range   = "2.9-12.4 kg",
    weight_mean    = "6.5 kg (SD 2.1)",
    sex_female_pct = 44.4,
    race_ethnicity = c(White = 78.9, Black = 10.5, Other = 9.0, Unknown = 1.5),
    disease_state  = "Infants aged 0-12 months with influenza confirmed by rapid diagnostic testing, RT-PCR, or viral culture, with influenza-like symptoms of <=96 h duration at enrolment.",
    dose_range     = "Oral suspension twice daily, weight-based and age-stratified: 2 mg/kg (0-1 month), 2.5 mg/kg (1-3 months), 3 mg/kg (3-12 months) in WP22849; 3 mg/kg (3-8 months) and 3 or 3.5 mg/kg (>=9 months) in CASG114. Ten doses (5 days) total, extendable by a further 5 days if clinically indicated.",
    regions        = "United States (CASG114, 2006-2010) and Europe, 11 academic centres (WP22849, 2011 and 2012 influenza seasons).",
    ga_range       = "24.0-43.0 weeks gestational age at birth (mean 38.1, SD 4.5); preterm neonates of postmenstrual age <36 weeks were excluded from both trials",
    pca_range      = "38.4-90.0 weeks postconceptual age (mean 61.6, SD 15.3)",
    n_observations = "604 oseltamivir and 648 oseltamivir carboxylate plasma samples (Kamal 2014 Results, first paragraph).",
    notes          = "Cohort statistics from Kamal 2014 Tables 1 and 2 (page 382). Study identifiers: CASG114 = NCT00391768 (NIAID Collaborative Antiviral Study Group, United States); WP22849 = NCT01286142 (Europe). Age-band enrolment per Table 1: <=30 days 13 subjects (9.8%), 31-90 days 33 (24.8%), 91-180 days 23 (17.3%), 181-270 days 35 (26.3%), >=271 days 29 (21.8%). NONMEM 7.2.0 was used for the population analysis (Methods, 'Population pharmacokinetic analysis'). Assay limits of quantitation were 1 ng/mL for oseltamivir and 50 ng/mL for OC; assay coefficients of variation were about 3% and 6% respectively. CAUTION on the Table 2 ethnicity rows as printed: Hispanic 34 (25.6%), Non-Hispanic/Latino 93 (69.9%), Other 12 (9.0%), Unknown 6 (4.5%) sum to 145 subjects and 109%, which exceeds the analysis N of 133; the values are recorded verbatim from the source and the discrepancy is unresolved in the publication. The race rows are internally consistent (105 + 14 + 12 + 2 = 133) and are the ones carried into race_ethnicity above. Ethnicity and race were not screened as covariates."
  )

  ini({
    # Structural fixed effects -- Kamal 2014 Table 3, page 383, "Estimate
    # (95% CI)" column (values cross-checked against the bootstrap-median
    # column, which agrees to within rounding for every parameter). All
    # clearance and volume terms are APPARENT: they are conditioned on
    # oral bioavailability F, and the OC terms are additionally
    # conditioned on the fraction of oseltamivir metabolised to OC, which
    # the authors fixed at 1 (Methods, "Population pharmacokinetic
    # analysis": "In the absence of data on OC administration, the
    # oseltamivir to OC conversion fraction and OC volumes could not be
    # identified. Therefore, complete (100%) oseltamivir to OC metabolism
    # was assumed and OC central volume was estimated.").
    # Reference subject: WT = 8 kg, postnatal age = 24 weeks.

    lka <- log(0.905); label("Oseltamivir first-order absorption rate constant ka (1/h)")                            # Table 3 theta1: ka = 0.905 1/h (95% CI 0.691-1.12), 12.1% RSE
    lcl <- log(80.4);  label("Oseltamivir apparent clearance CL/F at WT 8 kg (L/h)")                                 # Table 3 theta2: CL/F = 80.4 L/h (72.6-88.2), 4.96% RSE
    lvc <- log(165);   label("Oseltamivir apparent central volume V2/F at WT 8 kg (L)")                              # Table 3 theta3: V2/F = 165 L (124-209), 13.0% RSE
    lq  <- log(19.6);  label("Oseltamivir apparent intercompartmental clearance Q/F at WT 8 kg (L/h)")               # Table 3 theta4: Q/F = 19.6 L/h (17-22.1), 6.59% RSE
    lvp <- log(348);   label("Oseltamivir apparent peripheral volume V3/F at WT 8 kg (L)")                           # Table 3 theta5: V3/F = 348 L (152-544), 28.7% RSE

    lcl_oselcarb <- log(4.75); label("OC apparent clearance CLM/F at WT 8 kg, postnatal age 24 weeks (L/h)")         # Table 3 theta6: CLM/F = 4.75 L/h (4.34-5.16), 4.40% RSE
    lvc_oselcarb <- log(40.2); label("OC apparent central volume VM/F at WT 8 kg, postnatal age 24 weeks (L)")       # Table 3 theta7: VM/F = 40.2 L (36.6-43.8), 4.55% RSE

    # Allometric exponents. Kamal 2014 Table 3 footnote (page 383):
    # "Clearance (CL/F, Q/F, and CLM/F) and volume (V2/F, V3/F, and VM/F)
    # parameters were allometrically scaled as (WT/8)0.75 and (WT/8),
    # respectively, where WT = body weight." Both are fixed(), not
    # estimated -- Methods, "Allometric scaling and covariate model
    # development": "Clearance and volume parameters were scaled
    # allometrically with fixed power coefficients of 0.75 and 1" and
    # "Body weight dependencies were therefore fixed mechanistically".
    # Neither appears as a theta in Table 3, which is the confirming
    # signal: only theta1-theta11 were estimated.
    e_wt_cl <- fixed(0.75); label("Allometric exponent of (WT/8) on every clearance term CL/F, Q/F, CLM/F (unitless)")
    e_wt_vc <- fixed(1);    label("Allometric exponent of (WT/8) on every volume term V2/F, V3/F, VM/F (unitless)")

    # Postnatal-age effect on the two OC parameters. Kamal 2014 Table 3
    # footnote: "Age dependence was described by the linear function:
    # CLM/F~ 1 + CLM,AGE (AGE/24-1); VM/F ~1+VM,AGE (AGE/24-1),
    # VM,AGE = CLM,AGE, where AGE is subject's age in weeks."
    # ONE estimated coefficient drives BOTH OC parameters because the
    # authors constrained VM,AGE = CLM,AGE, so a single parameter is used
    # in both places rather than two numerically identical ones. The
    # effect is LINEAR in age, not a power function.
    e_age_cl_oselcarb <- 0.33; label("Linear coefficient of (postnatal age / 24 weeks - 1) on both OC CLM/F and OC VM/F (unitless)")  # Table 3 theta11: CLM,AGE = 0.33 (0.265-0.396), 10.1% RSE

    # Inter-individual variability. Kamal 2014 Table 3 reports the
    # interindividual covariance matrix Omega directly, as VARIANCES on
    # the eta (log) scale for the diagonal and as raw covariances for the
    # off-diagonal entries. The table footnote defines the printed
    # correlations as R = Omega(2,1)/sqrt(Omega(1,1)*Omega(2,2)), and all
    # three reproduce exactly from the tabulated entries
    # (0.157/sqrt(0.154*0.662) = 0.4917 vs printed 0.491;
    #  0.0638/sqrt(0.154*0.127) = 0.4562 vs printed 0.455;
    #  0.169/sqrt(0.662*0.127) = 0.5828 vs printed 0.584),
    # which confirms the entries are variances/covariances rather than
    # SDs or %CVs. No omega^2 = log(CV^2 + 1) conversion is applicable or
    # needed here. The block is positive definite as published (leading
    # principal minors 0.154, 0.0773, 0.00611), so no numerical nudge is
    # applied. IIV is present ONLY on CL/F, V2/F and CLM/F: Table 3 lists
    # no Omega for ka, Q/F, V3/F, or VM/F in the final model, even though
    # the base model started with "Random effects added to all
    # parameters" (Methods, "Population pharmacokinetic analysis").
    #
    # Per-entry source trace (all from Kamal 2014 Table 3, page 383;
    # trailing comments cannot be used inside the c() below because
    # rxode2's comment-to-label pass fails to parse them):
    #   Omega(1,1) omega^2 CL  = 0.154  (95% CI 0.101-0.208),   17.6% RSE, printed CV = 39.5%, shrinkage 7.2%
    #   Omega(2,1) cov(CL, V2) = 0.157  (95% CI 0.0576-0.256),  32.3% RSE, printed R = 0.491
    #   Omega(2,2) omega^2 V2  = 0.662  (95% CI 0.294-1.03),    28.4% RSE, printed CV = 81.4%, shrinkage 13.3%
    #   Omega(3,1) cov(CL, CLM)= 0.0638 (95% CI 0.0284-0.0991), 28.3% RSE, printed R = 0.455
    #   Omega(3,2) cov(V2, CLM)= 0.169  (95% CI 0.0875-0.251),  24.7% RSE, printed R = 0.584
    #   Omega(3,3) omega^2 CLM = 0.127  (95% CI 0.1-0.154),     10.7% RSE, printed CV = 35.7%, shrinkage 2.7%
    etalcl + etalvc + etalcl_oselcarb ~ c(
      0.154,
      0.157,  0.662,
      0.0638, 0.169,  0.127
    )

    # Residual error. The Kamal 2014 Table 3 footnote defines all three
    # sigma rows as STANDARD DEVIATIONS -- "sigmaOP,PROP, standard
    # deviation (oseltamivir), proportional; sigmaOC,PROP, standard
    # deviation (oseltamivir carboxylate), proportional; sigmaOC,ADD,
    # standard deviation (oseltamivir carboxylate), additive" -- and the
    # Variability column restates theta8 = 0.467 as "CV = 46.7%",
    # theta9 = 0.118 as "CV = 11.8%", and theta10 = 39.5 as "SD = 39.5".
    # The tabulated numbers are therefore already SDs on the linear
    # concentration scale; no square root is taken.
    # Note that the final model has NO additive residual term for the
    # parent -- Table 3 enumerates theta1 through theta11 and contains no
    # sigmaOP,ADD row -- even though Methods states that "combined
    # proportional and additive error models were used for oseltamivir
    # and OC" for the BASE model. See the vignette Assumptions and
    # deviations.
    propSd          <- 0.467;  label("Oseltamivir proportional residual SD (fraction)")  # Table 3 theta8: sigmaOP,PROP = 0.467 (0.420-0.515), 5.15% RSE, CV = 46.7%, shrinkage 10.0%
    propSd_oselcarb <- 0.118;  label("OC proportional residual SD (fraction)")           # Table 3 theta9: sigmaOC,PROP = 0.118 (0.094-0.141), 10.2% RSE, CV = 11.8%, shrinkage 12.1%
    addSd_oselcarb  <- 0.0395; label("OC additive residual SD (mg/L)")                   # Table 3 theta10: sigmaOC,ADD = 39.5 ng/mL (32.9-46.1), 8.53% RSE = 0.0395 mg/L
  })

  model({
    # Reference covariate values. The 8 kg weight reference is stated in
    # the Kamal 2014 Table 3 footnote and again in the Figure 2 caption.
    # The age reference is 24 WEEKS (Table 3 footnote, "AGE/24 ... where
    # AGE is subject's age in weeks"); the canonical PNA column carries
    # MONTHS, so the reference is written as the source's 24 weeks
    # converted in place (24 * 7 / 30.4375 = 5.5195 months) rather than as
    # a pre-computed constant, so the published number stays visible.
    ref_wt  <- 8
    ref_pna <- 24 * 7 / 30.4375

    # Allometric size factors, shared by every clearance term and every
    # volume term respectively (Table 3 footnote).
    allom_cl <- (WT / ref_wt)^e_wt_cl
    allom_v  <- (WT / ref_wt)^e_wt_vc

    # Linear postnatal-age factor, applied to BOTH OC parameters with the
    # single published coefficient (the VM,AGE = CLM,AGE constraint).
    age_oselcarb <- 1 + e_age_cl_oselcarb * (PNA / ref_pna - 1)

    # Individual parameters -- oseltamivir (parent, OP). IIV is on CL/F
    # and V2/F only; ka, Q/F and V3/F carry no eta in the final model.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * allom_cl
    vc <- exp(lvc + etalvc) * allom_v
    q  <- exp(lq)           * allom_cl
    vp <- exp(lvp)          * allom_v

    # Individual parameters -- oseltamivir carboxylate (OC, metabolite).
    # IIV is on CLM/F only; VM/F carries no eta in the final model.
    cl_oselcarb <- exp(lcl_oselcarb + etalcl_oselcarb) * allom_cl * age_oselcarb
    vc_oselcarb <- exp(lvc_oselcarb)                   * allom_v  * age_oselcarb

    # Micro-constants. kel is simultaneously the parent elimination rate
    # constant and the OC FORMATION rate constant: under the paper's
    # complete-conversion assumption the entire parent central-compartment
    # elimination flux enters the OC compartment.
    kel          <- cl / vc
    k12          <- q  / vc
    k21          <- q  / vp
    kel_oselcarb <- cl_oselcarb / vc_oselcarb

    # ODE system -- two-compartment oseltamivir with first-order
    # absorption plus one-compartment OC (Kamal 2014 Results, "The
    # population pharmacokinetic model", and Methods, "Population
    # pharmacokinetic analysis": "oseltamivir was described by a linear
    # two-compartment model with first-order absorption, and OC was
    # described by a one-compartment model ... First-order absorption and
    # direct first-order conversion of parent drug to metabolite were
    # assumed").
    d/dt(depot)            <- -ka * depot
    d/dt(central)          <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1)      <-  k12 * central - k21 * peripheral1
    d/dt(central_oselcarb) <-  kel * central - kel_oselcarb * central_oselcarb

    # Observations. Dose in mg divided by apparent volume in L gives
    # concentration in mg/L (the source reports ng/mL; 1 mg/L = 1000
    # ng/mL, and the OC additive residual SD has been converted
    # accordingly in ini()).
    Cc          <- central          / vc
    Cc_oselcarb <- central_oselcarb / vc_oselcarb

    Cc          ~ prop(propSd)
    Cc_oselcarb ~ add(addSd_oselcarb) + prop(propSd_oselcarb)
  })
}
