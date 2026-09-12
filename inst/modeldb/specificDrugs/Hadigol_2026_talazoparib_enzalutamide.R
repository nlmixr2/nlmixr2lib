Hadigol_2026_talazoparib_enzalutamide <- function() {
  description <- paste(
    "Coupled three-analyte population PK model for oral talazoparib",
    "(substrate of interest, bare canonical names), oral enzalutamide",
    "(perpetrator, sibling-drug suffix _enz) and its active N-desmethyl",
    "metabolite (suffix _ndmenz) in 811 men with metastatic",
    "castration-resistant prostate cancer from the Phase 3 TALAPRO-2 trial",
    "(NCT03395197, Part 1 plus Part 2 Cohort 1, unselected for homologous",
    "recombination repair deficiency). Each analyte is two-compartment;",
    "talazoparib and enzalutamide have first-order absorption and the",
    "metabolite is formed directly into its own central compartment. The",
    "enzalutamide / N-desmethyl parent-metabolite pair is unidentifiable in",
    "Fmet versus Vcn, so Fmet was FIXED to 0.634 from a published",
    "physiologically based PK model rather than estimated (Hadigol 2026",
    "Methods, Parent-Metabolite Modeling); enzalutamide leaves its central",
    "compartment at CLe/Fe, of which a fraction Fmet is routed to the",
    "metabolite and the remainder is true elimination. Talazoparib is a P-gp",
    "and BCRP substrate and enzalutamide plus its N-desmethyl metabolite are",
    "P-gp inhibitors, so talazoparib apparent clearance is INHIBITED by the",
    "summed plasma concentration of the two perpetrator species through the",
    "linear relationship CLt/Ft = CLt0/Ft * (1 - theta_slope * (Ce + Cn))",
    "(equation 1); the two species were treated as equipotent on the basis",
    "of in vitro data. At the typical steady-state enzalutamide 160 mg QD",
    "exposure this reduces talazoparib CL/F by about 20%, which is why the",
    "combination dose is talazoparib 0.5 mg QD rather than the 1 mg QD used",
    "as monotherapy. Baseline body weight and age are covariates on CLe/Fe",
    "and Vce/Fe, baseline body weight on CLn and Vcn, and baseline",
    "creatinine clearance on talazoparib base clearance CLt0/Ft via a power",
    "function (equation 6) - the only clinically relevant covariate the",
    "analysis retained, reducing CLt0/Ft by 8% for mild, 27% for moderate",
    "and 47% for severe renal impairment. The three sub-models were fit",
    "SEQUENTIALLY (enzalutamide, then metabolite with enzalutamide",
    "parameters fixed to empirical Bayes estimates, then talazoparib with",
    "both fixed), because a simultaneous fit of parent and metabolite was",
    "not stable; there is therefore no cross-analyte random-effect",
    "covariance to carry. Residual error is additive on log-transformed",
    "concentrations for all three analytes (Hadigol 2026)."
  )
  reference <- paste(
    "Hadigol M, Williams JH, Shi H, Yang DZ, Hoffman J, Wang DD.",
    "Population Pharmacokinetics Analysis of Talazoparib and Enzalutamide",
    "Combination Therapy for Patients With Metastatic Castration-Resistant",
    "Prostate Cancer. J Clin Pharmacol. 2026;66(2):e70125.",
    "doi:10.1002/jcph.70125."
  )
  vignette <- "Hadigol_2026_talazoparib_enzalutamide"

  # Amounts are carried in mg and volumes in L, so every clearance and
  # volume below is the value printed in Hadigol 2026 Tables 1-3 with no
  # rescaling. The three observed concentrations are converted to ng/mL
  # (1 mg/L = 1000 ng/mL) inside model() so that (a) they match the
  # bioanalytical assay units and the units of Table 5, and (b) the
  # drug-drug-interaction slope theta_slope can be used exactly as printed
  # in mL/ng without a hidden conversion factor.
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed baseline body weight; reference 80.1 kg, stated by",
        "Hadigol 2026 to be the population median (Table 5 Note: 'Baseline",
        "body weight and age were considered at the median values of",
        "population; that is, baseline body weight = 80.1 kg and age = 71",
        "years'). Enters in TWO different functional forms in one model, so",
        "the reference cannot be read off a single equation: a POWER form on",
        "enzalutamide CLe/Fe and Vce/Fe (equations 2 and 3, exponents 0.549",
        "and 3.495, both ESTIMATED rather than fixed to allometric values),",
        "and a LINEAR-deviation form on N-desmethyl enzalutamide CLn and Vcn",
        "(equations 4 and 5, coefficients 0.006 and 0.027 per kg from 80.1",
        "kg). Baseline body weight was also tested on talazoparib Vct/Ft and",
        "CLt0/Ft and was NOT retained (removed in the first backward",
        "elimination step, then re-tested and rejected). One patient had a",
        "missing baseline body weight, imputed with the population median",
        "(Methods, Data for Analysis)."
      ),
      source_name        = "BWT"
    ),
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference 71 years, the population median (Hadigol 2026 Table 5",
        "Note). Enters in TWO different functional forms: a LINEAR-deviation",
        "form on enzalutamide CLe/Fe, '1 - 0.003 * (Age - 71)' (equation 2),",
        "and a POWER form on Vce/Fe, '(Age / 71)^1.223' (equation 3). Age was",
        "also tested on talazoparib CLt0/Ft and was not a significant",
        "covariate. Age is highly correlated with baseline creatinine",
        "clearance in this cohort, which is why the two were not entered into",
        "the same stepwise search for enzalutamide (Results, Population PK",
        "Model for Enzalutamide)."
      ),
      source_name        = "Age"
    ),
    CRCL = list(
      description        = paste(
        "Baseline creatinine clearance, RAW (not BSA-normalized) mL/min"
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Raw un-normalised creatinine clearance in mL/min, NOT the",
        "mL/min/1.73 m^2 BSA-normalised default of the CRCL canonical - the",
        "same per-model variant recorded for Delattre 2010 amikacin and Chen",
        "2023 nemonoxacin. Hadigol 2026 calls the column BCCL (baseline",
        "creatinine clearance) and bands it by the standard renal-disease",
        "cutoffs in mL/min (Table 4: normal > 90, mild 60-89, moderate",
        "30-59, severe 15-29), which are raw-mL/min cutoffs. The paper does",
        "NOT state which estimating equation produced BCCL (Cockcroft-Gault",
        "is the usual choice for these cutoffs but is not named). Enters",
        "talazoparib base clearance as the POWER form '(BCCL / 86.85)^0.455'",
        "(equation 6); a linear form was tried first and replaced by the",
        "power form because the linear model gave a high condition number",
        "(Results, Population PK Model for Talazoparib). The reference value",
        "86.85 mL/min is the value Hadigol 2026 names for the typical",
        "patient; the paper does NOT state that it is the cohort median, and",
        "Table S2 (demographics) is not on disk, so it must not be reported",
        "as one. Baseline creatinine clearance was deliberately NOT entered",
        "into the enzalutamide or metabolite covariate search because it is",
        "highly correlated with both baseline body weight and age. Two",
        "patients had a missing baseline creatinine clearance, imputed with",
        "the population median (Methods, Data for Analysis)."
      ),
      source_name        = "BCCL"
    )
  )

  covariatesDataExcluded <- list(
    CONMED_CYP3A4_MOD = list(
      description = "Concomitant moderate CYP3A4 inhibitor",
      units       = NA_character_,
      type        = "binary",
      notes       = paste(
        "Screened and REJECTED twice. On enzalutamide CLe/Fe it entered the",
        "full model in the forward stepwise search but was eliminated in the",
        "single backward elimination step. On N-desmethyl enzalutamide CLn it",
        "survived backward elimination but was then removed because its",
        "relative standard error was 72%, which made the estimate unreliable",
        "and the model unstable (Results, Population PK Model for",
        "N-Desmethyl Enzalutamide). No point estimate is reported for either",
        "effect in any on-disk source, so neither can be encoded."
      )
    ),
    CONMED_PGP_INHIB_MOD = list(
      description = "Concomitant moderate P-glycoprotein inhibitor",
      units       = NA_character_,
      type        = "binary",
      notes       = paste(
        "Explored as a categorical covariate on talazoparib relative",
        "bioavailability Ft and on kat, and not retained. The effect of",
        "STRONG P-gp inhibitors was not evaluated at all, because too few",
        "patients took one (Results, Population PK Model for Talazoparib in",
        "Combination with Enzalutamide). The talazoparib MONOTHERAPY model",
        "did carry a potent-P-gp-inhibitor effect on relative",
        "bioavailability (Methods, Prior Knowledge and Modeling",
        "Experience), so its absence here is a property of this dataset, not",
        "of the drug."
      )
    ),
    RACE_ASIAN = list(
      description = "Asian race (versus non-Asian)",
      units       = NA_character_,
      type        = "binary",
      notes       = paste(
        "Tested on talazoparib CLt/Ft because race (Asian versus non-Asian)",
        "was a significant covariate on CLt/Ft in the talazoparib",
        "monotherapy model; not retained here. No point estimate is",
        "reported."
      )
    ),
    REGION_CHINA = list(
      description = "Chinese region (versus non-Chinese)",
      units       = NA_character_,
      type        = "binary",
      notes       = paste(
        "Tested on talazoparib CLt0/Ft and Vct/Ft and found not to be a",
        "significant covariate on either (Results and Conclusions). No point",
        "estimate is reported."
      )
    )
  )

  compartmentData <- list(
    depot                = list(analyte = "talazoparib",                  units = "mg", specimen = "administration site", verified = TRUE),
    central              = list(analyte = "talazoparib",                  units = "mg", specimen = "plasma",              verified = TRUE),
    peripheral1          = list(analyte = "talazoparib",                  units = "mg", specimen = "plasma",              verified = TRUE),
    depot_enz            = list(analyte = "enzalutamide",                 units = "mg", specimen = "administration site", verified = TRUE),
    central_enz          = list(analyte = "enzalutamide",                 units = "mg", specimen = "plasma",              verified = TRUE),
    peripheral1_enz      = list(analyte = "enzalutamide",                 units = "mg", specimen = "plasma",              verified = TRUE),
    central_ndmenz       = list(analyte = "N-desmethyl enzalutamide",     units = "mg", specimen = "plasma",              verified = TRUE),
    peripheral1_ndmenz   = list(analyte = "N-desmethyl enzalutamide",     units = "mg", specimen = "plasma",              verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 811L,
    n_studies      = 1L,
    age_median     = "71 years",
    weight_median  = "80.1 kg",
    sex_female_pct = 0,
    disease_state  = paste(
      "Metastatic castration-resistant prostate cancer (mCRPC). All",
      "patients are men. 792 patients came from TALAPRO-2 Part 2 Cohort 1,",
      "which enrolled patients UNSELECTED for homologous recombination",
      "repair (HRR) gene deficiencies; the HRR-deficient Cohort 2 was NOT",
      "pooled into this analysis. The remaining 19 patients came from the",
      "open-label, non-randomized Part 1, whose primary objective was to",
      "choose the Part 2 starting dose of talazoparib. 414 patients were",
      "randomized to talazoparib plus enzalutamide and 397 to placebo plus",
      "enzalutamide, so all 811 contribute enzalutamide and N-desmethyl",
      "enzalutamide data while only the 412 of 414 with at least one",
      "non-BLQ talazoparib observation contribute talazoparib data."
    ),
    dose_range     = paste(
      "Enzalutamide 160 mg orally once daily, with protocol-permitted",
      "reductions to 120 mg and 80 mg once daily. Talazoparib 0.5 mg orally",
      "once daily, reduced to a 0.35 mg once-daily starting dose for",
      "patients with moderate renal impairment. Part 1 initially used",
      "talazoparib 1 mg once daily in combination, which produced roughly",
      "two-fold higher than expected talazoparib trough concentrations at",
      "Week 5 and a higher than expected frequency of Grade 3/4",
      "haematologic adverse events, and was reduced to 0.5 mg once daily.",
      "Only the hard-capsule talazoparib formulation was used."
    ),
    renal_function = paste(
      "Baseline creatinine clearance spans normal through severe",
      "impairment. Hadigol 2026 Table 4 evaluates the covariate at the",
      "cutoff midpoints 90 (normal), 75 (mild), 45 (moderate) and 22.5",
      "(severe) mL/min; the Table 4 Note records that there were NO",
      "patients with end-stage renal disease in the dataset."
    ),
    notes          = paste(
      "Data are pooled from Part 1 and Part 2 Cohort 1 of the Phase 3",
      "TALAPRO-2 trial (NCT03395197). Dense plasma sampling was designed",
      "for Part 1 and only SPARSE sampling for Part 2, which is the main",
      "reason absorption parameters are imprecisely informed relative to",
      "disposition. The dataset comprised 6030 plasma observations for",
      "enzalutamide and for N-desmethyl enzalutamide and 2961 for",
      "talazoparib. Concentrations below the limit of quantification were",
      "set to missing rather than handled by an M3/M4 likelihood method,",
      "because the BLQ fraction was small. Assay ranges were 20 to 20,000",
      "ng/mL for enzalutamide and for N-desmethyl enzalutamide (LLOQ 20",
      "ng/mL) and 25 to 25,000 pg/mL for talazoparib (LLOQ 25 pg/mL).",
      "NONMEM 7.4.3 with FOCE-I; stepwise covariate model building in PsN",
      "4.9.0 at forward-addition p = 0.05 and backward-deletion p = 0.001.",
      "Parameter uncertainty was assessed by PsN sampling importance",
      "resampling (2000 samples / 1000 resamples), reported here as the SIR",
      "median and 95% CI in the ini() source-trace comments. Tables S1 and",
      "S2 (study-population summary and baseline demographics) are NOT on",
      "disk, so age, weight and creatinine-clearance RANGES and the",
      "race/ethnicity and regional breakdown cannot be recorded; only the",
      "medians named in the main text are given above."
    )
  )

  ini({
    # =====================================================================
    # TALAZOPARIB -- the substrate of interest, so it carries the BARE
    # canonical parameter names.
    # Source: Hadigol 2026 Table 3, 'Talazoparib Final Model
    # Pharmacokinetics Parameters Summary'. All clearances and volumes are
    # APPARENT (divided by the relative bioavailability Ft, which was fixed
    # to 1). The final model has NO absorption lag time, unlike the
    # talazoparib monotherapy model described in Methods, Prior Knowledge
    # and Modeling Experience.
    # =====================================================================
    lcl <- log(5.078)
    label("Talazoparib apparent BASE clearance CLt0/Ft, i.e. clearance in the absence of enzalutamide and N-desmethyl enzalutamide exposure (L/h)") # Table 3 theta CLt0/Ft = 5.078 L/h, SE 0.197, RSE 3.887%, SIR median 5.072 [4.821, 5.344]
    lvc <- log(14.389)
    label("Talazoparib apparent central volume of distribution Vct/Ft (L)") # Table 3 theta Vct/Ft = 14.389 L, SE 0.066, RSE 0.462%, SIR median 14.382 [14.067, 14.735]
    lvp <- log(382.135)
    label("Talazoparib apparent peripheral volume of distribution Vpt/Ft (L)") # Table 3 theta Vpt/Ft = 382.135 L, SE 18.094, RSE 4.735%, SIR median 382.208 [364.472, 399.944]
    lq <- log(15.568)
    label("Talazoparib apparent inter-compartmental clearance Qt/Ft (L/h)") # Table 3 theta Qt/Ft = 15.568 L/h, SE 0.008, RSE 0.053%, SIR median 15.568 [15.536, 15.598]
    lka <- log(0.157)
    label("Talazoparib first-order absorption rate constant kat (1/h)") # Table 3 theta kat = 0.157 1/h, SE 1.11e-4, RSE 0.071%, SIR median 0.157 [0.156, 0.157]
    lfdepot <- fixed(log(1))
    label("Talazoparib relative bioavailability Ft (unitless, held at 1)") # Table 3 row 'theta Ft (Fixed)' = 1.000, no SE reported; this is what makes CLt/Ft and Vct/Ft apparent rather than absolute

    # Talazoparib covariate effect: baseline creatinine clearance on the
    # apparent BASE clearance, as a power function. Hadigol 2026
    # equation (6): CLt0/Ft = 5.078 * (BCCL / 86.85)^0.455.
    e_crcl_cl <- 0.455
    label("Power exponent for baseline creatinine clearance on talazoparib apparent base clearance CLt0/Ft (unitless)") # Table 3 row 'BCCL effect on theta CLt0/Ft' = 0.455, SE 0.051, RSE 11.111%, SIR median 0.452 [0.38, 0.523]; equation (6)

    # Drug-drug interaction: talazoparib is a P-gp / BCRP substrate and
    # enzalutamide plus its N-desmethyl metabolite inhibit P-gp, so the
    # apparent clearance of talazoparib falls linearly with the SUMMED
    # plasma concentration of the two perpetrator species. Hadigol 2026
    # equation (1): CLt/Ft = CLt0/Ft * (1 - theta_slope * (Ce + Cn)).
    # The two species are summed with EQUAL weight because they were
    # 'considered similar in terms of their potency, which is supported by
    # in vitro studies' (Methods, Drug-Drug Interactions). Printed in
    # mL/ng, so it multiplies a concentration in ng/mL -- which is why
    # model() converts all three concentrations to ng/mL.
    e_cenz_cl <- 6.58e-6
    label("Linear slope of the inhibitory effect of summed enzalutamide plus N-desmethyl enzalutamide plasma concentration on talazoparib apparent clearance (mL/ng)") # Table 3 theta slope = 6.58e-6 mL/ng, SE 1.22e-6, RSE 18.514%, SIR median 6.56e-6 [4.98e-6, 8.21e-6]; equation (1). The base talazoparib model gave 7.43e-6 mL/ng, slightly reduced to this value in the final model

    # =====================================================================
    # ENZALUTAMIDE -- the perpetrator. Sibling-drug suffix _enz.
    # Source: Hadigol 2026 Table 1, 'Enzalutamide Final Model
    # Pharmacokinetics Parameters Summary'. Fe was fixed to 1, so every
    # enzalutamide clearance and volume is apparent.
    # =====================================================================
    lcl_enz <- log(0.425)
    label("Enzalutamide apparent clearance CLe/Fe (L/h)") # Table 1 theta CLe/Fe = 0.425 L/h, SE 0.003, RSE 0.74%, SIR median 0.425 [0.419, 0.431]
    lvc_enz <- log(25.293)
    label("Enzalutamide apparent central volume of distribution Vce/Fe (L)") # Table 1 theta Vce/Fe = 25.293 L, SE 2.862, RSE 11.31%, SIR median 25.155 [21.773, 28.655]
    lvp_enz <- log(45.928)
    label("Enzalutamide apparent peripheral volume of distribution Vpe/Fe (L)") # Table 1 theta Vpe/Fe = 45.928 L, SE 2.836, RSE 6.17%, SIR median 46.001 [43.1, 48.806]
    lq_enz <- log(20.644)
    label("Enzalutamide apparent inter-compartmental clearance Qe/Fe (L/h)") # Table 1 theta Qe/Fe = 20.644 L/h, SE 1.505, RSE 7.29%, SIR median 20.655 [18.274, 23.322]
    lka_enz <- log(3.431)
    label("Enzalutamide first-order absorption rate constant kae (1/h)") # Table 1 theta kae = 3.431 1/h, SE 0.475, RSE 13.85%, SIR median 3.423 [2.79, 4.552]
    lfdepot_enz <- fixed(log(1))
    label("Enzalutamide relative bioavailability Fe (unitless, held at 1)") # Table 1 row 'theta Fe (Fixed)' = 1.000, no SE reported

    # Enzalutamide covariate effects. Hadigol 2026 equations (2) and (3):
    #   CLe/Fe = 0.425 * (BWT / 80.1)^0.549 * [1 - 0.003 * (Age - 71)]
    #   Vce/Fe = 25.293 * (BWT / 80.1)^3.495 * (Age / 71)^1.223
    # Note the two covariates enter CLe/Fe and Vce/Fe in DIFFERENT
    # functional forms: body weight is a power term on both parameters,
    # whereas age is a linear-deviation term on clearance but a power term
    # on volume. The body-weight exponent on Vce/Fe is 3.495, far above any
    # allometric value; it is estimated, and the reported SE and SIR
    # interval both confirm the magnitude, so it is transcribed as printed.
    e_wt_cl_enz <- 0.549
    label("Power exponent for baseline body weight on enzalutamide apparent clearance CLe/Fe (unitless)") # Table 1 row 'Body weight effect on theta CLe/Fe' = 0.549, SE 0.035, RSE 6.42%, SIR median 0.55 [0.485, 0.619]; equation (2)
    e_age_cl_enz <- -0.003
    label("Linear coefficient for age on enzalutamide apparent clearance CLe/Fe, per year from the 71-year reference (1/year)") # Table 1 row 'Age effect on theta CLe/Fe' = -0.003, SE 0.001, RSE 27.53%, SIR median -0.003 [-0.005, -0.002]; equation (2)
    e_wt_vc_enz <- 3.495
    label("Power exponent for baseline body weight on enzalutamide apparent central volume Vce/Fe (unitless)") # Table 1 row 'Body weight effect on theta Vce/Fe' = 3.495, SE 0.279, RSE 7.99%, SIR median 3.505 [3.105, 3.904]; equation (3)
    e_age_vc_enz <- 1.223
    label("Power exponent for age on enzalutamide apparent central volume Vce/Fe (unitless)") # Table 1 row 'Age effect on theta Vce/Fe' = 1.223, SE 0.271, RSE 22.17%, SIR median 1.238 [0.647, 1.842]; equation (3)

    # =====================================================================
    # N-DESMETHYL ENZALUTAMIDE -- the active metabolite of enzalutamide.
    # Metabolite suffix _ndmenz, following the ndm<drug> family.
    # Source: Hadigol 2026 Table 2, 'N-Desmethyl Enzalutamide Final Model
    # Pharmacokinetics Parameters Summary'. These clearances and volumes
    # are NOT divided by a bioavailability term, because the metabolite is
    # formed systemically and the upstream Fe is fixed to 1.
    # =====================================================================
    lcl_ndmenz <- log(0.286)
    label("N-desmethyl enzalutamide clearance CLn (L/h)") # Table 2 theta CLn = 0.286 L/h, SE 0.003, RSE 0.89%, SIR median 0.286 [0.281, 0.291]
    lvc_ndmenz <- log(44.944)
    label("N-desmethyl enzalutamide central volume of distribution Vcn (L)") # Table 2 theta Vcn = 44.944 L, SE 4.189, RSE 9.32%, SIR median 45.032 [41.657, 48.433]
    lvp_ndmenz <- log(52.373)
    label("N-desmethyl enzalutamide peripheral volume of distribution Vpn (L)") # Table 2 theta Vpn = 52.373 L, SE 3.788, RSE 7.23%, SIR median 52.341 [49.486, 54.751]
    lq_ndmenz <- log(9.443)
    label("N-desmethyl enzalutamide inter-compartmental clearance Qn (L/h)") # Table 2 theta Qn = 9.443 L/h, SE 1.349, RSE 14.28%, SIR median 9.433 [8.157, 10.909]

    # Fraction of enzalutamide routed to the N-desmethyl metabolite. FIXED,
    # not estimated: with plasma observations for parent and metabolite
    # only, Vcn and Fmet are unidentifiable because an increase in Vcn can
    # be compensated by a decrease in Fmet as long as the total apparent
    # clearance of enzalutamide stays equal to CLe/Fe. Hadigol 2026 fixed
    # Fmet to 0.634 from a published enzalutamide physiologically based PK
    # model, which fit better than the alternative literature approach of
    # assuming equal central volumes for parent and metabolite (Methods,
    # Parent-Metabolite Modeling).
    fm <- fixed(0.634)
    label("Fraction of enzalutamide clearance routed to the N-desmethyl metabolite Fmet (unitless)") # Table 2 row 'theta Fmet (Fixed)' = 0.634, no SE reported; fixed from a published enzalutamide PBPK model, not estimated

    # N-desmethyl enzalutamide covariate effects. Hadigol 2026 equations
    # (4) and (5), both LINEAR deviations from the 80.1 kg reference (not
    # power terms, unlike the parent):
    #   CLn = 0.286  * [1 + 0.006 * (BWT - 80.1)]
    #   Vcn = 44.944 * [1 + 0.027 * (BWT - 80.1)]
    e_wt_cl_ndmenz <- 0.006
    label("Linear coefficient for baseline body weight on N-desmethyl enzalutamide clearance CLn, per kg from the 80.1 kg reference (1/kg)") # Table 2 row 'Effect of body weight theta CLn' = 0.006, SE 0.001, RSE 9.42%, SIR median 0.006 [0.005, 0.007]; equation (4)
    e_wt_vc_ndmenz <- 0.027
    label("Linear coefficient for baseline body weight on N-desmethyl enzalutamide central volume Vcn, per kg from the 80.1 kg reference (1/kg)") # Table 2 row 'Effect of body weight on theta Vcn' = 0.027, SE 0.001, RSE 2.52%, SIR median 0.027 [0.026, 0.028]; equation (5)

    # =====================================================================
    # INTER-INDIVIDUAL VARIABILITY. Exponential (multiplicative) random
    # effects, so each omega below is a VARIANCE on the log scale. All
    # three tables report the same footnote, 'CV, approximate percent
    # coefficient of variation calculated as sqrt(omega^2) x 100%', which
    # confirms the tabulated numbers are variances and not CVs and lets
    # each value be used directly. Enzalutamide and the metabolite each
    # have a full block Omega; talazoparib has no reported off-diagonal,
    # so its two etas stay diagonal. The three sub-models were fit
    # SEQUENTIALLY with upstream parameters fixed to empirical Bayes
    # estimates, so there is no cross-analyte covariance to carry.
    # =====================================================================
    etalcl ~ 0.073       # Table 3 row 'omega^2 theta CLt0/Ft' = 0.073, SE 0.007, RSE 10.207%, CV 26.984%, shrinkage 10.20%, SIR median 0.074 [0.062, 0.087]
    etalvc ~ 2.582       # Table 3 row 'omega^2 Vct/Ft' = 2.582, SE 0.319, RSE 12.364%, CV 160.674%, shrinkage 46.81%, SIR median 2.55 [2.023, 3.325]

    etalcl_enz + etalvc_enz ~ c(0.037,
                                -0.015, 0.350) # Table 1 full block Omega: 'omega^2 CLe/Fe' = 0.037 (CV 19.32%, shrinkage 4.82%); 'omega Vce/Fe omega CLe/Fe' = -0.015 (RSE 65.35%); 'omega^2 Vce/Fe' = 0.350 (CV 59.18%, shrinkage 39.53%). Implied correlation -0.13

    etalcl_ndmenz + etalvc_ndmenz ~ c(0.057,
                                      0.006, 0.342) # Table 2 full block Omega: 'omega^2 CLn' = 0.057 (CV 23.83%, shrinkage 5.36%); 'omega Vcn omega CLn' = 0.006 (RSE 91.42%); 'omega^2 Vcn' = 0.342 (CV 58.44%, shrinkage 14.75%). Implied correlation 0.04

    # =====================================================================
    # RESIDUAL ERROR. All three analytes were fit to LOG-TRANSFORMED
    # observations with an additive residual, which nlmixr2 expresses as
    # `lnorm(expSd)` where expSd is the log-scale residual SD. Each table
    # reports the SD directly under the row name 'Thetarized sigma', the
    # NONMEM idiom in which $SIGMA is fixed to 1 and the residual SD is
    # carried as an estimated THETA -- which is why these are SDs and not
    # variances, and why they have their own RSE and shrinkage.
    # =====================================================================
    expSd <- 0.354
    label("Talazoparib log-scale additive residual SD (unitless)") # Table 3 row 'Thetarized sigma' = 0.354, SE 0.013, RSE 3.676%, shrinkage 8.50%, SIR median 0.354 [0.344, 0.365]; Results: 'The additive residual variability for log-transformed observation data was 0.354'
    expSd_enz <- 0.145
    label("Enzalutamide log-scale additive residual SD (unitless)") # Table 1 row 'Thetarized sigma' = 0.145, SE 0.010, RSE 6.86%, shrinkage 8.99%, SIR median 0.146 [0.143, 0.149]
    expSd_ndmenz <- 0.125
    label("N-desmethyl enzalutamide log-scale additive residual SD (unitless)") # Table 2 row 'Thetarized sigma' = 0.125, SE 0.007, RSE 5.20%, shrinkage 11.65%, SIR median 0.126 [0.123, 0.128]
  })

  model({
    # -----------------------------------------------------------------
    # 1. Enzalutamide individual parameters. Hadigol 2026 equations (2)
    #    and (3). Body weight is a power term on both CLe/Fe and Vce/Fe;
    #    age is a LINEAR deviation on CLe/Fe but a POWER term on Vce/Fe.
    # -----------------------------------------------------------------
    cl_enz <- exp(lcl_enz + etalcl_enz) * (WT / 80.1)^e_wt_cl_enz * (1 + e_age_cl_enz * (AGE - 71))
    vc_enz <- exp(lvc_enz + etalvc_enz) * (WT / 80.1)^e_wt_vc_enz * (AGE / 71)^e_age_vc_enz
    vp_enz <- exp(lvp_enz)
    q_enz  <- exp(lq_enz)
    ka_enz <- exp(lka_enz)

    # -----------------------------------------------------------------
    # 2. N-desmethyl enzalutamide individual parameters. Hadigol 2026
    #    equations (4) and (5) -- linear body-weight deviations, not
    #    power terms.
    # -----------------------------------------------------------------
    cl_ndmenz <- exp(lcl_ndmenz + etalcl_ndmenz) * (1 + e_wt_cl_ndmenz * (WT - 80.1))
    vc_ndmenz <- exp(lvc_ndmenz + etalvc_ndmenz) * (1 + e_wt_vc_ndmenz * (WT - 80.1))
    vp_ndmenz <- exp(lvp_ndmenz)
    q_ndmenz  <- exp(lq_ndmenz)

    # -----------------------------------------------------------------
    # 3. Perpetrator plasma concentrations, in ng/mL. Amounts are in mg
    #    and volumes in L, so central/vc is mg/L and 1 mg/L = 1000
    #    ng/mL. Converting here means theta_slope below is used exactly
    #    as printed in mL/ng, and means Cc_enz / Cc_ndmenz come out in
    #    the assay's own units (LLOQ 20 ng/mL).
    # -----------------------------------------------------------------
    Cc_enz    <- 1000 * central_enz    / vc_enz
    Cc_ndmenz <- 1000 * central_ndmenz / vc_ndmenz

    # -----------------------------------------------------------------
    # 4. Talazoparib individual parameters. Two stages:
    #    (a) equation (6), the renal-function power effect giving the
    #        apparent BASE clearance cl0 = CLt0/Ft; then
    #    (b) equation (1), the linear P-gp-inhibition term giving the
    #        actual, TIME-VARYING apparent clearance
    #          CLt/Ft = CLt0/Ft * (1 - theta_slope * (Ce + Cn)).
    #    The quantity that actually drives elimination must be the
    #    variable named `cl`, so the interaction is folded into `cl`
    #    itself rather than into a separately-named total.
    #    When the enzalutamide depot is never dosed, central_enz and
    #    central_ndmenz stay 0, the bracket evaluates to exactly 1, and
    #    cl falls back to the uninhibited base clearance.
    # -----------------------------------------------------------------
    cl0 <- exp(lcl + etalcl) * (CRCL / 86.85)^e_crcl_cl
    cl  <- cl0 * (1 - e_cenz_cl * (Cc_enz + Cc_ndmenz))
    vc  <- exp(lvc + etalvc)
    vp  <- exp(lvp)
    q   <- exp(lq)
    ka  <- exp(lka)

    # -----------------------------------------------------------------
    # 5. ODE system. Talazoparib first so that the bare canonical
    #    `depot` / `central` / `peripheral1` occupy the leading
    #    compartment slots, matching the substrate-of-interest role that
    #    gives them the unsuffixed names.
    #
    #    Enzalutamide leaves its central compartment at the full
    #    apparent clearance cl_enz; a fraction `fm` of that flux is
    #    routed into the metabolite's central compartment and the
    #    remaining (1 - fm) is true elimination. Writing it this way
    #    keeps the total apparent clearance of enzalutamide equal to
    #    CLe/Fe regardless of fm, which is exactly the constraint
    #    Hadigol 2026 identifies as the reason fm and Vcn cannot both be
    #    estimated. The metabolite has no depot: it is formed directly
    #    into plasma.
    # -----------------------------------------------------------------
    d/dt(depot)              <- -ka * depot
    d/dt(central)            <-  ka * depot - (cl / vc) * central - (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1)        <-  (q / vc) * central - (q / vp) * peripheral1

    d/dt(depot_enz)          <- -ka_enz * depot_enz
    d/dt(central_enz)        <-  ka_enz * depot_enz - (cl_enz / vc_enz) * central_enz - (q_enz / vc_enz) * central_enz + (q_enz / vp_enz) * peripheral1_enz
    d/dt(peripheral1_enz)    <-  (q_enz / vc_enz) * central_enz - (q_enz / vp_enz) * peripheral1_enz

    d/dt(central_ndmenz)     <-  fm * (cl_enz / vc_enz) * central_enz - (cl_ndmenz / vc_ndmenz) * central_ndmenz - (q_ndmenz / vc_ndmenz) * central_ndmenz + (q_ndmenz / vp_ndmenz) * peripheral1_ndmenz
    d/dt(peripheral1_ndmenz) <-  (q_ndmenz / vc_ndmenz) * central_ndmenz - (q_ndmenz / vp_ndmenz) * peripheral1_ndmenz

    # -----------------------------------------------------------------
    # 6. Relative bioavailability. Both Ft and Fe were fixed to 1, so
    #    these evaluate to exactly 1 and are present only to make the
    #    apparent-parameter convention explicit and machine-readable.
    # -----------------------------------------------------------------
    f(depot)     <- exp(lfdepot)
    f(depot_enz) <- exp(lfdepot_enz)

    # -----------------------------------------------------------------
    # 7. Observations. All three in ng/mL; all three fit to
    #    log-transformed data with an additive residual.
    # -----------------------------------------------------------------
    Cc <- 1000 * central / vc

    Cc        ~ lnorm(expSd)
    Cc_enz    ~ lnorm(expSd_enz)
    Cc_ndmenz ~ lnorm(expSd_ndmenz)
  })
}
