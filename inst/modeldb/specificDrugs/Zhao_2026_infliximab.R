Zhao_2026_infliximab <- function() {
  description <- "Two-compartment population PK model of intravenous infliximab in young paediatric patients with inflammatory bowel disease aged 10 years or younger, with allometric body-weight scaling (exponents fixed at 0.75 on clearance terms and 1 on volume terms) and power effects of serum albumin and C-reactive protein on clearance (Zhao 2026). Developed on 640 serum concentrations from 104 children across 14 European and Canadian centres; 14.4 percent of measurements were below the limit of quantification and were handled with the M3 method. Peripheral volume V2 and inter-compartmental clearance Q could not be estimated from the sparse trough-dominated data and are FIXED to the values of the published Chung and Clemente-Bautista paediatric models. Inter-individual variability sits on clearance only. Residual variability is proportional with six separate magnitudes selected per observation by which commercial infliximab ELISA measured the sample, the paper's central methodological contribution. Typical clearance in this cohort (0.779 L/day/65 kg) is roughly twice that reported for older children and adults."
  reference <- "Zhao Q, Jongsma MME, Vuijk SA, de Winter BCM, Martinez-Vinson C, Kolho KL, Norsa L, Hussey S, Wine E, Cohen S, Shouval DS, Assa A, Lev-Tzion R, de Meij T, Wolters VM, Huynh HQ, Preijers T, de Ridder L. Population Pharmacokinetics Analysis of Infliximab in up to 10-Year-Old Patients with Paediatric Inflammatory Bowel Disease: Label-Recommended Dose Fails to Achieve Therapeutic Target Concentration. Clin Pharmacokinet. 2026;65(1):79-94. doi:10.1007/s40262-025-01565-6"
  vignette <- "Zhao_2026_infliximab"
  paper_specific_residual_sds <- c(
    "propSdSanquin", "propSdImmundiagnostik", "propSdCaltag",
    "propSdMatrixBiotek", "propSdBenHorin", "propSdPromonitor"
  )
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Zhao 2026 Results 3.3 ("A
  # two-compartmental popPK model best described the IFX disposition in this
  # population") and Methods 2.1 / Supplementary Table 1, which describe serum
  # infliximab measured by commercial ELISA kits. Infliximab is given as a
  # short intravenous infusion, so there is no absorption compartment.
  compartmentData <- list(
    central     = list(analyte = "infliximab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "infliximab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on every disposition parameter, normalised to 65 kg, with the exponents FIXED a priori rather than estimated: 0.75 on CL and Q, 1 on V1 and V2 (Zhao 2026 Eq. 3 and Eq. 4, both of which carry (body weight/65 kg)^0.75 explicitly, plus Results 3.3 'The allometric exponent was fixed at 0.75 for all CL terms and 1 for all volume of distribution terms a priori' and 'Estimating the exponents for allometric scaling using body weight did not improve the model fit'). The 65 kg reference is an adult-sized normalising constant that the authors used so their estimates could be compared with the published adult and older-paediatric models in Tables 3 and 4; it is far OUTSIDE this cohort's weight range of 9.5-40.9 kg (median 25 kg, Table 1), so the reported CL and V1 values are extrapolations and should not be read as predictions for a 65 kg patient. Zhao 2026 Discussion makes this point directly: 'it could not be ignored that the covariate distribution impacts scaling of CL and V1, ultimately impacting the estimation of the typical values (Eq. 3)'. Time-varying in principle (children grow over a multi-year treatment course), but Table 1 tabulates weight at the first infusion only.",
      source_name        = "body weight"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on clearance only, normalised to the cohort median: (ALB/40.5)^-0.758 (Zhao 2026 Eq. 3 with theta_cov from Table 2 'ALB on CL' = -0.758, RSE 32%, bootstrap median -0.754 with 95% CI -1.270 to -0.246). The 40.5 g/L centring value is the Table 1 total-cohort median and is named as the reference in Results 3.3: 'CL decreased by approximately 17% for every 10 g/L increase in albumin from the reference value of 40.5 g/L'. That printed 17% figure is reproduced by a 35 to 45 g/L step, (45/35)^-0.758 = 0.827, rather than by a step starting at the 40.5 g/L reference itself, (50.5/40.5)^-0.758 = 0.846 or a 15.4% decrease -- the exponent is what is load-bearing and both readings confirm its sign and order of magnitude. Zhao 2026 reports albumin in SI g/L already (Table 1 median 40.5, range 19.6-50, mean 39.6 +/- 4.8), so no unit conversion is applied and the canonical SI column is used directly. Time-varying: albumin is a routine laboratory value drawn alongside each therapeutic-drug-monitoring sample.",
      source_name        = "ALB"
    ),
    CRP = list(
      description        = "C-reactive protein",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on clearance only, normalised to the cohort median: (CRP/3)^0.0545 (Zhao 2026 Eq. 3 with theta_cov from Table 2 'CRP on CL' = 0.0545, RSE 36%, bootstrap median 0.0541 with 95% CI 0.016-0.093). The 3 mg/L centring value is the Table 1 total-cohort MEDIAN, named as the reference in Results 3.3, and Methods 2.3.2 confirms the convention ('For the continuous covariates, values were normalized to the population median') -- do NOT substitute the Table 1 mean of 10.4 mg/L, which is far higher because CRP is a strongly right-skewed acute-phase reactant. The centring value and the exponent are pinned exactly by the printed sensitivity statement 'CL increased by approximately 1.6% for every 1 mg/L increase in CRP from the reference value of 3 mg/L': (4/3)^0.0545 = 1.0158, i.e. 1.58%. Zhao 2026 reports CRP in SI mg/L (Table 1 median 3, range 0.02-302), matching the canonical column, so no conversion is applied. CAUTION: the power form is undefined at CRP = 0 and collapses clearance toward zero as CRP approaches zero; the observed minimum in this cohort is 0.02 mg/L, at which the multiplier is 0.762, so users simulating a CRP of exactly 0 must floor the value (the assay's lower reporting limit is the natural floor). Time-varying: CRP is drawn alongside each therapeutic-drug-monitoring sample.",
      source_name        = "CRP"
    ),
    ASSAY_IFX = list(
      description        = "Integer (1-6) per-observation code naming which commercial infliximab ELISA measured the sample, used to select the proportional residual-error magnitude. 1 = Sanquin (Amsterdam), 2 = Immundiagnostik AG (Germany), 3 = Caltag Laboratories, 4 = Matrix Biotek, 5 = Shomron Ben-Horin in-house assay (Sheba Medical Center), 6 = Promonitor.",
      units              = "(integer 1-6)",
      type               = "categorical",
      reference_category = "None -- the residual-error model is a pure six-way stratification with no reference level; every observation belongs to exactly one assay and each assay selects its own estimated proportional residual SD. Code 2 (Immundiagnostik) is the sensible default for a user who does not know which kit produced a measurement: it is the assay the paper itself pivots on in the Table 5 pairwise comparisons, it is one of the two commercial kits the Discussion identifies as mutually harmonised, and its estimate sits mid-range.",
      notes              = "Per-observation (record-level) code. Stratifying the residual error by assay is the paper's central methodological contribution: Zhao 2026 Methods 2.3.1 states that an additive, proportional and combined error model was evaluated 'for each of the six different [IFX] assays', and the Discussion argues that the published models' failure to do so is why they predicted this cohort so poorly ('in the external validation of the published models, the data from different IFX assays were lumped into one residual error model ... Not surprisingly, the published models did not predict sufficiently'). The six code levels are the six assays tabulated in Zhao 2026 Table 1 and given one Table 2 residual row each. Supplementary Table 1 maps 14 contributing centres onto assay vendors and lists two labels that do NOT appear among the six modelled strata: 'apDIA' (Universitair Ziekenhuis Gent) contributed no patients to the PK analysis, and 'Sanquin (LC-MS/MS)' (Utrecht) is subsumed into the Sanquin stratum; 'Immune diagnostic' (Dublin) and 'Promonitor/Sanquin' (Tampere) are the Immundiagnostik and Promonitor strata under variant spellings. An integer stratum column rather than a family of paired binary indicators follows the ASSAY_CMIA register entry's own guidance for a dataset mixing more than two assays and the STUDY_TACRO_FRANCKE precedent, and the model file derives the (ASSAY_IFX == n) indicators inline. Auto-approved member of the ASSAY_<method> canonical family. Note that the assay enters ONLY the residual error, never a structural parameter or a scaling factor on Cc, so predicted exposure is by construction independent of the code -- which is the analytic form of the paper's own conclusion that 'the differences obtained across the IFX assays did not affect overall drug exposure' (Table 5, AUC week 6-14 P > 0.05 for all 15 pairwise comparisons).",
      source_name        = "IFX assay"
    )
  )

  # Covariates that Zhao 2026 screened in the univariate and stepwise
  # multivariate covariate search (Methods 2.3.2) but did NOT retain in the
  # final model: only albumin and CRP survived forward inclusion (P < 0.05) and
  # backward elimination (P < 0.01), per Results 3.3 "Albumin and CRP were
  # selected as covariate relationships associated with CL". The paper reports
  # no point estimate, no objective-function change and no confidence interval
  # for any rejected covariate, so none of these can be encoded. Recorded here
  # so the provenance of the covariate screen is not lost; documentation only,
  # never referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened on CL; not retained. Cohort median 8.2 years (range 1.2-10.0), mean 7.6 +/- 2.0 per Table 1. Age is strongly collinear with body weight in this cohort, and Results 3.5 attributes the apparent CL difference between the under-6 and 6-to-10 age strata to body size rather than to age itself: 'this observed difference was primarily attributable to differences in body weight between the two groups'."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened on CL; not retained. Cohort median 124.5 cm (range 74-150), mean 123.5 +/- 14.0 per Table 1."
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = "Screened on CL; not retained. Cohort median 1.0 m^2 (range 0.4-1.9), mean 1.0 +/- 0.2 per Table 1. Body weight was retained instead, as the a priori allometric scaling term."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL as a categorical covariate (Methods 2.3.2 'gender'); not retained. Table 1 reports 54/104 (52%) for the total cohort. NOTE: Table 1 prints two consecutive rows both labelled 'Female, n (%)', 54 (52%) and 50 (48%), which sum to 104; the second row is evidently male and the duplicated label is a typesetting error."
    ),
    DIS_CD = list(
      description = "Crohn's disease indicator (versus ulcerative colitis)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL as a categorical covariate (Methods 2.3.2 'diagnosis type (Crohn's disease or ulcerative colitis)'); not retained. Table 1 reports 59/104 (57%) Crohn's disease and 45/104 (43%) ulcerative colitis. The complementary DIS_UC canonical describes the same dichotomy from the other side."
    ),
    CONMED_IMMUNOMOD = list(
      description = "Concomitant immunomodulator indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL as a categorical covariate (Methods 2.3.2 'immunomodulatory use (with or without)'); not retained. Table 1 reports 46 patients on an immunomodulator and 58 not; the accompanying percentages, 59.0% and 41.0%, are transposed relative to the counts (46/104 = 44.2% and 58/104 = 55.8%)."
    ),
    ESR = list(
      description = "Erythrocyte sedimentation rate",
      units       = "mm/h",
      type        = "continuous",
      notes       = "Screened on CL; not retained. Cohort median 18.8 mm/h (range 0.1-124), mean 19.6 +/- 13.6 per Table 1. Retained on CL by four of the six published comparator models in Table 3 (Bauman, Xiong, Clemente-Bautista, Colman), which makes its rejection here a genuine finding rather than an omission. No canonical column is registered for this concept because it is screened and rejected with no reported point estimate and is never referenced in model(); the name here is the paper's own abbreviation."
    ),
    PGA = list(
      description = "Physician's Global Assessment of disease activity, a four-level ordinal (remission / mild / moderate / severe)",
      units       = "(4-level ordinal)",
      type        = "categorical",
      notes       = "Screened on CL as a categorical covariate (Methods 2.3.2); not retained. Table 1's PGA block is internally inconsistent: the total-cohort column reads remission 2 (2%), mild 12 (11.5%), moderate 54 (52.0%), severe 36 (34.6%), summing to 104, but the under-6 column reads remission 10 (63%) with moderate 0, which cannot be reconciled with a total-cohort remission count of 2. No canonical column is registered for this concept: the existing BLPHYVAS canonical is a 0-100 visual-analogue scale rather than this four-level ordinal, and PGA was rejected here with no reported point estimate and is never referenced in model(), so the name used is the paper's own abbreviation."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 104L,
    n_studies      = 1L,
    n_centres      = 14L,
    age_range      = "1.2-10.0 years",
    age_median     = "8.2 years",
    weight_range   = "9.5-40.9 kg",
    weight_median  = "25 kg",
    sex_female_pct = 52,
    race_ethnicity = NULL,
    disease_state  = "paediatric inflammatory bowel disease (57% Crohn's disease, 43% ulcerative colitis); baseline C-reactive protein mean 10.4 +/- 22.3 mg/L and erythrocyte sedimentation rate mean 19.6 +/- 13.6 mm/h, both higher than published values for children over 10 years, indicating a heavier inflammatory burden",
    dose_range     = "intravenous infliximab, median 6 mg/kg (range 3.5-15) per infusion; label induction at weeks 0, 2 and 6 followed by maintenance every 8 weeks, with intensified regimens in clinical use",
    regions        = "Europe (Netherlands, France, Finland, Italy, Ireland, Czech Republic, Israel, Belgium) and Canada",
    n_observations = 640L,
    notes          = "Baseline demographics are Zhao 2026 Table 1 (total-cohort column); covariate values are recorded at the first infusion. Retrospective therapeutic-drug-monitoring data, patients treated 2004-2016 with data collected 2015-2019, from 14 centres (Supplementary Table 1). Every child was under 10 years old at the start of infliximab therapy, and infliximab is off-label below 6 years -- 16 of the 104 patients were under 6. Observation counts differ between sections of the paper: the Abstract and Results 3.1 report 2150 measured concentrations (the full therapeutic-drug-monitoring collection, of which only four were peak samples), while Results 3.3 states that the model-development dataset comprised 640 measurements from the same 104 patients, of which 95 (14.4%) were below the limit of quantification and were handled with the M3 method. n_observations records the 640 actually fitted. The stated 14.4% does not equal 95/640 = 14.8%, a minor arithmetic slip in the source. Race and ethnicity are not reported."
  )

  ini({
    # Structural parameters, all reported for a 65 kg reference subject.
    # Zhao 2026 Table 2, "Final model estimates" column.
    lcl <- log(0.779); label("Clearance (L/day) for a 65 kg subject")                                  # Table 2 final model: CL = 0.779 L/day/65 kg (RSE 7%); bootstrap median 0.775, 95% CI 0.615-0.943
    lvc <- log(17.2);  label("Central volume of distribution (L) for a 65 kg subject")                 # Table 2 final model: V1 = 17.2 L/65 kg (RSE 16%); bootstrap median 17.1, 95% CI 10.831-23.629

    # V2 and Q could not be estimated from this trough-dominated dataset (only
    # four peak samples were available) and were FIXED to values taken from the
    # published paediatric models of Chung et al. and Clemente-Bautista et al.
    # Zhao 2026 Results 3.3: "Owing to the low number of peak concentrations,
    # the volume of distribution of peripheral compartment (V2) and
    # inter-compartment clearance (Q) could not be estimated and, therefore,
    # fixed to values from the published models".
    #
    # REFERENCE-WEIGHT AMBIGUITY: Zhao 2026 Table 3 prints the donor values as
    # V2 = 1.21 L/37.4 kg and Q = 0.0697 L/day/37.4 kg in the Chung column
    # (Clemente-Bautista's are 1.19 L/46.4 kg and 0.0696 L/day/46.4 kg), while
    # the refined-model column and Table 2 print the SAME NUMBERS against a
    # 65 kg reference. Zhao therefore adopted the numeric values without
    # renormalising them from the donor's reference weight. These are encoded as
    # PRINTED for this model -- 65 kg, matching Table 2 -- because the paper's
    # own table for its own model is the authority. The alternative reading
    # (Chung's 37.4 kg reference) would make V2 and Q 1.7x larger at any weight
    # and lengthen the terminal half-life by 13.2%; it fits Zhao 2026 Table 6
    # very slightly worse (median absolute difference 3.45 vs 3.11 percentage
    # points over the 18 cells), which is not a decisive arbitration because Q
    # is an order of magnitude below CL and the peripheral compartment barely
    # influences troughs. See the vignette's Assumptions section.
    lvp <- fixed(log(1.21));   label("Peripheral volume of distribution (L) for a 65 kg subject")      # Table 2 final model: V2 = 1.21 L/65 kg (FIX), also Table 3 refined-model column
    lq  <- fixed(log(0.0697)); label("Inter-compartmental clearance (L/day) for a 65 kg subject")      # Table 2 final model: Q = 0.0697 L/day/65 kg (FIX), also Table 3 refined-model column

    # Allometric exponents, fixed a priori rather than estimated.
    e_wt_cl_q   <- fixed(0.75); label("Allometric body-weight exponent on CL and Q (unitless)")        # Eq. 3 and Eq. 4 carry (body weight/65 kg)^0.75 for clearance; Results 3.3 'The allometric exponent was fixed at 0.75 for all CL terms'
    e_wt_vc_vp  <- fixed(1);    label("Allometric body-weight exponent on V1 and V2 (unitless)")       # Eq. 3 and Eq. 4 text: 'if the parameter is volume of distribution, the exponential for body weight/65 kg is 1'; Results 3.3 'and 1 for all volume of distribution terms a priori'

    # Covariate effects on clearance, both power exponents on a
    # median-normalised covariate (Eq. 3).
    e_alb_cl <- -0.758; label("Power exponent for serum albumin on CL (unitless)")                     # Table 2 final model: 'ALB on CL' = -0.758 (RSE 32%); bootstrap median -0.754, 95% CI -1.270 to -0.246
    e_crp_cl <- 0.0545; label("Power exponent for C-reactive protein on CL (unitless)")                # Table 2 final model: 'CRP on CL' = 0.0545 (RSE 36%); bootstrap median 0.0541, 95% CI 0.016-0.093

    # Inter-individual variability on clearance only. Table 2 reports it as
    # 41.1% CV, converted to a log-normal variance by omega^2 = log(1 + CV^2)
    # = log(1 + 0.411^2) = 0.156.
    etalcl ~ 0.156  # Table 2 final model row 'CL (%CV)' = 41.1% (RSE 12%, shrinkage 15%); bootstrap median 39.8%, 95% CI 26.88%-52.35%. Base model was 46.3%, so the two covariates explained 5.2 percentage points

    # Proportional residual error, one magnitude per infliximab ELISA. Table 2
    # prints these six values as bare NONMEM $SIGMA VARIANCES; the SDs that
    # nlmixr2's prop() expects are their square roots. See the vignette's
    # "Assumptions and deviations" section for the %RSE argument that
    # establishes the variance scale.
    propSdSanquin         <- 0.975; label("Proportional residual error, Sanquin assay (fraction)")             # Table 2 final model 'Prop. RUV for Sanquin' = 0.951 (RSE 27%, eps-shrinkage 33%); bootstrap 0.941, 95% CI 0.387-1.516. sqrt(0.951) = 0.975
    propSdImmundiagnostik <- 0.884; label("Proportional residual error, Immundiagnostik assay (fraction)")     # Table 2 final model 'Prop. RUV for Immundiagnostik' = 0.781 (RSE 9%, eps-shrinkage 14%); bootstrap 0.777, 95% CI 0.603-0.959. sqrt(0.781) = 0.884
    propSdCaltag          <- 0.785; label("Proportional residual error, Caltag assay (fraction)")              # Table 2 final model 'Prop. RUV for Caltag' = 0.617 (RSE 9%, eps-shrinkage 10%); bootstrap 0.617, 95% CI 0.509-0.724. sqrt(0.617) = 0.785
    propSdMatrixBiotek    <- 0.769; label("Proportional residual error, Matrix Biotek assay (fraction)")       # Table 2 final model 'Prop. RUV for Matrix biotek' = 0.591 (RSE 20%, eps-shrinkage 6%); bootstrap 0.584, 95% CI 0.320-0.862. sqrt(0.591) = 0.769
    propSdBenHorin        <- 0.943; label("Proportional residual error, Ben-Horin in-house assay (fraction)")  # Table 2 final model 'Prop. RUV for Shomron Ben-Horin' = 0.889 (RSE 25%, eps-shrinkage 25%); bootstrap 0.848, 95% CI 0.088-1.690. sqrt(0.889) = 0.943
    propSdPromonitor      <- 0.956; label("Proportional residual error, Promonitor assay (fraction)")          # Table 2 final model 'Prop. RUV for Promonitor' = 0.914 (RSE 31%, eps-shrinkage 11%); bootstrap 0.861, 95% CI -0.021-1.848. sqrt(0.914) = 0.956
  })

  model({
    # 1. Individual parameters. Eq. 3 multiplies the typical value by each
    #    median-normalised covariate raised to its power exponent, by
    #    exp(eta), and by the allometric weight term. The albumin and CRP
    #    centring constants are the Table 1 cohort medians named in
    #    Results 3.3 (40.5 g/L and 3 mg/L); the 65 kg weight reference is the
    #    normalising constant in Eq. 3 and Eq. 4 themselves.
    cl <- exp(lcl + etalcl) * (ALB / 40.5)^e_alb_cl * (CRP / 3)^e_crp_cl * (WT / 65)^e_wt_cl_q
    vc <- exp(lvc) * (WT / 65)^e_wt_vc_vp
    vp <- exp(lvp) * (WT / 65)^e_wt_vc_vp
    q  <- exp(lq)  * (WT / 65)^e_wt_cl_q

    # 2. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. ODE system. Infliximab is administered as a short intravenous
    #    infusion, so the dose enters the central compartment directly and
    #    there is no absorption or bioavailability term.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 4. Observation. Doses are in mg and volumes in L, so Cc is in mg/L,
    #    which is numerically identical to the ug/mL used for the paper's
    #    5 mg/L target trough.
    Cc <- central / vc

    # 5. Residual error. The proportional magnitude is selected per
    #    observation by which ELISA measured the sample (Methods 2.3.1:
    #    the error model was evaluated "for each of the six different [IFX]
    #    assays"), so it is assembled as a model variable and passed to
    #    prop(). Exactly one indicator is 1 for any valid ASSAY_IFX code.
    propSdAssay <-
      propSdSanquin         * (ASSAY_IFX == 1) +
      propSdImmundiagnostik * (ASSAY_IFX == 2) +
      propSdCaltag          * (ASSAY_IFX == 3) +
      propSdMatrixBiotek    * (ASSAY_IFX == 4) +
      propSdBenHorin        * (ASSAY_IFX == 5) +
      propSdPromonitor      * (ASSAY_IFX == 6)

    Cc ~ prop(propSdAssay)
  })
}
