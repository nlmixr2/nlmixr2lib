Li_2017_nhl_pfs_mbma <- function() {
  description <- paste(
    "MBMA. Model-based meta-analysis of PROGRESSION-FREE SURVIVAL (PFS) in",
    "non-Hodgkin lymphoma (NHL), fitted to study-arm-level summary PFS",
    "Kaplan-Meier curves digitised from 112 published clinical trials",
    "(155 cohorts, 11,824 patients, 3,098 observations, reported 1993-2015).",
    "PFS time follows a two-parameter Weibull distribution",
    "S(t) = exp(-(t / lambda)^k) with scale lambda = 110 months and shape",
    "(Weibull slope) k = 0.79 in the reference cohort, so the median PFS",
    "time is lambda * log(2)^(1/k) = 69.1 months. The reference cohort is",
    "rituximab as the only treatment, 100% mantle cell lymphoma (MCL) and",
    "100% treatment-naive. Three covariate blocks act log-linearly on",
    "lambda and none was retained on k: cohort regimen (rituximab,",
    "CHOP/CHOP-like, bendamustine, other non-chemotherapy drugs; bortezomib",
    "and other chemotherapy were eliminated), the cohort's NHL-histology",
    "mix (follicular lymphoma, DLBCL, other NHL against an MCL reference),",
    "and the cohort's fraction of treatment-experienced patients. Because",
    "the median PFS time is directly proportional to lambda, each",
    "coefficient exponentiates to a median-PFS-time ratio: bendamustine",
    "4.0-fold, rituximab 3.1-fold, CHOP/CHOP-like 2.3-fold, other drugs",
    "2.3-fold, follicular lymphoma 1.6-fold, DLBCL 0.24-fold and an",
    "all-treatment-experienced cohort 0.16-fold. Variability is",
    "MBMA-specific and between-STUDY only: inter-study variance 1.4 on",
    "log-lambda and 0.17 on log-k. There is no ODE, no dosing and no",
    "PK layer; time is in months and the model outputs the",
    "progression-free probability sur, the cumulative hazard cumhaz and",
    "the median PFS time tmed. Suitable simulation scope is study-arm",
    "PFS curves, median PFS times and between-arm hazard ratios, NOT",
    "individual-patient event times. Parameter values are Supplemental",
    "Content Table S1 (Final Model column); the Weibull form is",
    "Supplemental Content Model Description equation 1, whose published",
    "rendering carries a typesetting error corrected here (see the",
    "vignette Errata).",
    sep = " "
  )
  reference <- paste(
    "Li M, Dave N, Salem AH, Freise KJ.",
    "Model-based meta-analysis of progression-free survival in",
    "non-Hodgkin lymphoma patients.",
    "Medicine (Baltimore). 2017;96(35):e7988.",
    "doi:10.1097/MD.0000000000007988.",
    "The model equation, the covariate-model narrative and the final",
    "parameter estimates are in the Supplemental Content",
    "(Model Description; Table S1), http://links.lww.com/MD/B853;",
    "the main text carries no parameter table.",
    sep = " "
  )
  vignette <- "Li_2017_nhl_pfs"

  units <- list(
    time = "month",
    dosing = "n/a (no dosing events; the cohort regimen enters only through the binary treatment-indicator covariates)",
    concentration = paste(
      "unitless (the outputs are the progression-free probability sur in",
      "0-1, the cumulative hazard cumhaz and the median PFS time tmed in",
      "months; there is no drug concentration in this model)"
    )
  )

  covariateData <- list(
    CONMED_RITUX = list(
      description = "Study-arm regimen indicator: 1 = the cohort's regimen includes rituximab (alone or with any chemotherapy or non-chemotherapy backbone), 0 = the cohort received no rituximab.",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (rituximab-containing regimen) for THIS model -- see notes; the canonical register defines the column itself with 0 = no rituximab",
      notes = paste(
        "MBMA study-arm-level indicator (a property of the trial cohort,",
        "not of an individual patient). Rituximab was used in 98 of the",
        "155 cohorts (63.2%; Li 2017 Table 1).",
        "POLARITY IS LOAD-BEARING AND INVERTED RELATIVE TO THE PAPER'S",
        "PRINTED COEFFICIENT. Table S1 tabulates the effect as",
        "'No rituximab on lambda' = -1.13, i.e. the coefficient applies",
        "to the ABSENCE of rituximab, because the Supplemental Content",
        "Model Description defines the reference population as",
        "'rituximab as the only treatment'. The canonical column here",
        "keeps the register's 1 = on-rituximab convention, so model()",
        "multiplies the coefficient by (1 - CONMED_RITUX). A downstream",
        "user who flips this gets a 3.1-fold error in median PFS in the",
        "wrong direction, and no mass-balance or AUC check exists in a",
        "survival model to catch it -- only the sign of the",
        "median-PFS-time ratio against Li 2017 Figure 3A does.",
        "The canonical entry's own Notes anticipate exactly this case:",
        "document the value transformation rather than registering a",
        "second reverse-coded canonical."
      ),
      source_name = "Rituximab (Li 2017 Table 1 'Type of treatment'); 'No rituximab on lambda' (Table S1 coefficient row)"
    ),
    CONMED_CHOP = list(
      description = "Study-arm regimen indicator: 1 = the cohort's regimen includes CHOP (cyclophosphamide, doxorubicin/hydroxydaunorubicin, vincristine, prednisone) or a CHOP-like regimen, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no CHOP or CHOP-like backbone)",
      notes = paste(
        "MBMA study-arm-level indicator. Li 2017 Methods defines",
        "'CHOP-like' as a regimen sharing 3 of the 4 CHOP drugs, so this",
        "single column pools CHOP with its close variants exactly as the",
        "source model did; a finer decomposition is not recoverable from",
        "the paper. Present in 41 of 155 cohorts (26.4%; Table 1).",
        "Not mutually exclusive with CONMED_RITUX: R-CHOP cohorts carry",
        "both columns at 1 and receive both coefficients additively on",
        "log-lambda."
      ),
      source_name = "CHOP/CHOP-like (Li 2017 Table 1 'Type of treatment' and Table S1 coefficient row)"
    ),
    CONMED_BENDAMUSTINE = list(
      description = "Study-arm regimen indicator: 1 = the cohort's regimen includes bendamustine, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no bendamustine)",
      notes = paste(
        "MBMA study-arm-level indicator. The smallest treatment stratum",
        "in the database: 11 of 155 cohorts (7%; Li 2017 Table 1), which",
        "is why its coefficient carries the largest treatment-effect",
        "shift between the final model (1.39) and the",
        "complete-covariate sensitivity analysis (0.526, RSE 54%).",
        "Bendamustine carries the largest single median-PFS-time ratio",
        "in the model, exp(1.39) = 4.0-fold (Li 2017 Figure 3A)."
      ),
      source_name = "Bendamustine (Li 2017 Table 1 'Type of treatment' and Table S1 coefficient row)"
    ),
    CONMED_NONCHEMO_OTHER = list(
      description = "Study-arm regimen indicator: 1 = the cohort's regimen includes at least one non-chemotherapy anticancer drug other than rituximab and other than bortezomib, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no other non-chemotherapy drug)",
      notes = paste(
        "MBMA study-arm-level residual-bucket indicator, and the second",
        "most common treatment level in the database: 60 of 155 cohorts",
        "(38.7%; Li 2017 Table 1). Li 2017 Methods enumerates the pooled",
        "agents as temsirolimus, GM-CSF, enzastaurin, celecoxib,",
        "tositumomab, alisertib, lenalidomide, pidilizumab, idelalisib,",
        "flavopiridol, cladribine, ibrutinib, everolimus, SAM486A,",
        "dacetuzumab, bevacizumab, galiximab, obinutuzumab, inotuzumab,",
        "thalidomide, epratuzumab, dexamethasone and ofatumumab. Because",
        "one coefficient is shared across that whole list, the column",
        "cannot distinguish a BTK inhibitor from a growth factor, and it",
        "must not be read as a class effect for any single agent.",
        "The Abstract's phrase is 'other nonchemotherapy drugs, aside",
        "from bortezomib', which is what fixes bortezomib OUT of this",
        "bucket and into its own (eliminated) level; see",
        "covariatesDataExcluded$CONMED_BORTEZOMIB."
      ),
      source_name = "Other drugs (Li 2017 Table 1 'Type of treatment' and Table S1 'Other drugs on lambda')"
    ),
    TUMTP_FL_PCT = list(
      description = "Study-arm-level percentage (0-100) of the cohort whose NHL histology is follicular lymphoma (FL).",
      units = "%",
      type = "continuous",
      reference_category = "0% (no FL patients); the model's reference histology is 100% mantle cell lymphoma, i.e. all three of TUMTP_FL_PCT, TUMTP_DLBCL_PCT and TUMTP_OTHER_NHL_PCT at 0",
      notes = paste(
        "MBMA study-arm-level covariate; the aggregate counterpart of the",
        "individual-level binary TUMTP_FL. Li 2017 screened 'percentage of",
        "patients with each NHL subtype', so a cohort with a mixed",
        "histology takes intermediate values -- 25 of 155 cohorts (16.1%)",
        "are classified 'Mixed' in Table 1 and are exactly why the",
        "covariate is a continuous prevalence rather than a per-cohort",
        "histology label. Cohorts predominantly FL: 44 of 155 (28.4%).",
        "SCALE IS LOAD-BEARING. The canonical _PCT family is in percent",
        "(0-100) while Li 2017's fitted coefficients are per unit",
        "FRACTION -- Table 1 prints the sibling aggregate covariates as",
        "fractions ('performance status >= 2, median 0.08, range 0-0.6'),",
        "and Figure 3B places an all-FL cohort at a median-PFS-time ratio",
        "of exp(0.46) = 1.58, which is only reproducible if the",
        "coefficient multiplies 1 rather than 100. model() therefore",
        "divides by 100, following the RACE_ASIAN_PCT / PS_ECOG_0_PCT",
        "precedent in Franzese_2026_pdl1_nsclc_mbma. Supplying a fraction",
        "in this column instead of a percent shrinks every histology",
        "effect to a 100th of its size and is invisible in any summary",
        "statistic of the column itself."
      ),
      source_name = "Percentage of patients with follicular lymphoma (Li 2017 Methods; Table 1 'NHL subtype'; Table S1 'Follicular lymphoma on lambda')"
    ),
    TUMTP_DLBCL_PCT = list(
      description = "Study-arm-level percentage (0-100) of the cohort whose NHL histology is diffuse large B-cell lymphoma (DLBCL).",
      units = "%",
      type = "continuous",
      reference_category = "0% (no DLBCL patients); the reference histology is 100% MCL",
      notes = paste(
        "MBMA study-arm-level covariate; the aggregate counterpart of the",
        "individual-level binary TUMTP_DLBCL. The largest histology",
        "stratum, 49 of 155 cohorts (31.6%; Li 2017 Table 1), and the",
        "largest negative effect in the model: an all-DLBCL cohort has",
        "exp(-1.41) = 0.24 of the reference median PFS time (Figure 3B),",
        "consistent with DLBCL being the aggressive high-grade subtype.",
        "Scale is percent (0-100) and model() divides by 100 -- see the",
        "TUMTP_FL_PCT notes for why the fraction-versus-percent reading",
        "is settled by Figure 3B rather than assumed."
      ),
      source_name = "Percentage of patients with DLBCL (Li 2017 Methods; Table 1 'NHL subtype'; Table S1 'DLBCL on lambda')"
    ),
    TUMTP_OTHER_NHL_PCT = list(
      description = "Study-arm-level percentage (0-100) of the cohort whose NHL histology is neither follicular lymphoma, nor DLBCL, nor mantle cell lymphoma (e.g. peripheral T-cell lymphoma, marginal zone lymphoma).",
      units = "%",
      type = "continuous",
      reference_category = "0% (no other-NHL patients); the reference histology is 100% MCL",
      notes = paste(
        "MBMA study-arm-level covariate; the aggregate counterpart of the",
        "individual-level binary TUMTP_OTHER_NHL, but with a narrower",
        "membership: here MCL is EXCLUDED from the bucket because MCL is",
        "the model's reference histology and carries its own (zero-by-",
        "construction) level, whereas the individual-level canonical",
        "pools MCL into the residual. Li 2017 Introduction names",
        "peripheral T-cell lymphoma and marginal zone lymphoma as the",
        "members. Only 4 of 155 cohorts (2.5%) are predominantly",
        "other-NHL (Table 1), which is why this is the one effect the",
        "paper reports as not distinguishable from the MCL reference:",
        "-0.17 with a 113% RSE and a 95% CI of -0.54 to 0.21 that",
        "brackets zero. It survived backward elimination only as part of",
        "the tumor-subtype block. Treat exp(-0.17) = 0.84 as",
        "indistinguishable from 1, matching the paper's own statement",
        "that other subtypes were 'similar to the MCL patient",
        "population'. Scale is percent (0-100); model() divides by 100."
      ),
      source_name = "Percentage of patients with other NHL subtypes (Li 2017 Methods; Table 1 'NHL subtype' 'Other'; Table S1 'Other lymphomas on lambda')"
    ),
    LINE_1L_PCT = list(
      description = "Study-arm-level percentage (0-100) of the cohort that is treatment-naive, i.e. receiving the trial regimen as first-line therapy with no prior anticancer treatment for NHL.",
      units = "%",
      type = "continuous",
      reference_category = "100% (an all-treatment-naive cohort) -- this is the model's reference level, NOT 0%",
      notes = paste(
        "MBMA study-arm-level covariate; the aggregate counterpart of the",
        "individual-level binary LINE_1L, whose 1 = first-line /",
        "treatment-naive polarity it keeps.",
        "POLARITY IS INVERTED RELATIVE TO THE PAPER'S PRINTED",
        "COEFFICIENT, for the same reason as CONMED_RITUX. Table S1",
        "tabulates '% experienced patients on lambda' = -1.85, so the",
        "coefficient applies to the fraction of TREATMENT-EXPERIENCED",
        "patients, while the Supplemental Content Model Description sets",
        "the reference population at '100% treatment naive patients'.",
        "model() therefore forms the experienced fraction as",
        "(100 - LINE_1L_PCT) / 100. Li 2017 Figure 3C is the arbiter and",
        "confirms both ends: 'All naive patients' sits exactly on the",
        "reference line at 1.0 and 'All experienced patients' at",
        "exp(-1.85) = 0.16.",
        "This is the single largest covariate effect in the model -- a",
        "6.4-fold swing in median PFS time from one end of the column to",
        "the other, larger than any treatment or histology effect.",
        "Li 2017 Results report a SUBANALYSIS (not in Table S1, and not",
        "encoded here) finding treatment-naive and one-prior-line cohorts",
        "statistically indistinguishable while two-or-more-prior-line",
        "cohorts had under a tenth the median PFS time; so this linear",
        "column approximates what is really a threshold at two prior",
        "lines, and the approximation is the source model's, not this",
        "extraction's. Scale is percent (0-100)."
      ),
      source_name = "Percentage of patients with different numbers of prior treatments (naive, 1, 2+) (Li 2017 Methods); '% experienced patients on lambda' (Table S1 coefficient row)"
    )
  )

  covariatesDataExcluded <- list(
    TUMTP_MCL_PCT = list(
      description = "Study-arm-level percentage (0-100) of the cohort whose NHL histology is mantle cell lymphoma (MCL).",
      units = "%",
      type = "continuous",
      notes = "The model's REFERENCE histology, represented implicitly by TUMTP_FL_PCT, TUMTP_DLBCL_PCT and TUMTP_OTHER_NHL_PCT all being 0, so it is deliberately not referenced in model(). Declared here to record the reference level explicitly, following the Lu_2017_polatuzumab_neuropathy treatment of TUMTP_FL. 33 of 155 cohorts (21.3%) are predominantly MCL (Li 2017 Table 1)."
    ),
    CONMED_BORTEZOMIB = list(
      description = "Study-arm regimen indicator: 1 = the cohort's regimen includes bortezomib, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      notes = "Screened as its own treatment level (22 of 155 cohorts, 14.2%; Li 2017 Table 1) but ELIMINATED in backward elimination at P < 0.001 with a 95% CI including 0, on both lambda and k. No point estimate is published for the final model (Table S1 prints '---'), so no coefficient can be transcribed. Li 2017 Figure 3A places bortezomib exactly on the median-PFS-time-ratio = 1.0 reference line, and the Abstract states bortezomib did not prolong median PFS. A cohort receiving bortezomib is represented in this model by its other regimen components only."
    ),
    CONMED_CHEMO_OTHER = list(
      description = "Study-arm regimen indicator: 1 = the cohort's regimen includes a chemotherapy regimen that is neither CHOP/CHOP-like nor bendamustine, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      notes = "Screened (36 of 155 cohorts, 23.2%; Li 2017 Table 1) but ELIMINATED in backward elimination; Table S1 prints '---' for the final model and Figure 3A places it on the 1.0 reference line. Li 2017 Methods enumerates the pooled regimens (vorinostat, fludarabine, MCP, DHAP, ICE, BEAM, pentostatin, GemOx, NIMP, THP, gemcitabine, CMD, FM, FC, pirarubicin, PEPC, IFE, AD+C, cytarabine, mitoxantrone, FMD, COP-X). Distinct from the registered CONMED_CHEMO, which means any chemotherapy backbone rather than this residual bucket."
    ),
    AGE = list(
      description = "The cohort's MEDIAN age, not an individual patient's age.",
      units = "year",
      type = "continuous",
      notes = "Screened on both lambda and k and not retained (P > 0.01). Cohort medians: 63 years overall, range 47-83 across the 155 cohorts; not reported in 8 cohorts (5.1%) and imputed as the database median there (Li 2017 Table 1 and Methods). Li 2017 Discussion: 'Impact of age and sex was not found to be significant.' An aggregate median rather than the individual-level AGE the canonical register defines."
    ),
    SEXF_PCT = list(
      description = "Study-arm-level percentage (0-100) of the cohort that is female.",
      units = "%",
      type = "continuous",
      notes = "Screened and not retained on lambda; entered forward selection on k and was then ELIMINATED in backward elimination. The source column is the percentage of MALE patients (median 59%, range 29-91; not reported in 14 cohorts (9%) and imputed as the median), so the canonical female-referenced column is SEXF_PCT = 100 - %male. No point estimate is published for either polarity."
    ),
    PRIOR_RITUX_PCT = list(
      description = "Study-arm-level percentage (0-100) of the cohort previously administered rituximab before enrolling in the trial.",
      units = "%",
      type = "continuous",
      notes = "Screened and not retained (P > 0.01). DISTINCT FROM CONMED_RITUX, which records rituximab in the cohort's own trial regimen and IS retained -- prior exposure and current exposure are separate covariates with opposite fates in this model. Missing in 12 studies and imputed as the median. Li 2017 Discussion notes the null result agrees with Johnston et al., who found no PFS difference between rituximab-naive and rituximab-retreated patients."
    ),
    PS_ECOG_2_PCT = list(
      description = "Study-arm-level percentage (0-100) of the cohort with ECOG performance status 2 at baseline.",
      units = "%",
      type = "continuous",
      notes = "Li 2017's screened covariate is the percentage with performance status >= 2, which spans ECOG 2 and ECOG 3; per the register's PS_ECOG_2_PCT notes the two are kept as separate columns and a model that needs the composite sums them inside model(). Entered forward selection on both lambda and k and was then ELIMINATED in backward elimination, so no coefficient exists and no summation is needed here. Cohort values (as fractions in the source): median 0.08, range 0-0.6; not reported in 40 cohorts (25.8%) (Li 2017 Table 1)."
    ),
    PS_ECOG_3_PCT = list(
      description = "Study-arm-level percentage (0-100) of the cohort with ECOG performance status 3 at baseline.",
      units = "%",
      type = "continuous",
      notes = "The second half of Li 2017's screened performance-status >= 2 composite; see PS_ECOG_2_PCT. Screened and eliminated, so no coefficient exists. Li 2017 does not report the ECOG 2 and ECOG 3 components separately."
    ),
    STAGE_GE3_PCT = list(
      description = "Study-arm-level percentage (0-100) of the cohort with Ann Arbor disease stage III or IV (advanced-stage disease) at baseline.",
      units = "%",
      type = "continuous",
      notes = "Screened on both lambda and k and ELIMINATED in backward elimination, so no coefficient exists. Cohort values (as fractions in the source): median 0.86, range 0-1; not reported in 33 cohorts (21.3%) (Li 2017 Table 1). Documented here only; no register entry is minted for a covariate that never reaches model()."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 11824L,
    n_studies = 112L,
    n_cohorts = 155L,
    n_randomized_studies = 22L,
    n_observations = "3,098 study-arm-level PFS observations digitised from published Kaplan-Meier curves",
    age_range = "cohort median ages 47-83 years; median across the 155 cohorts 63 years (not reported in 8 cohorts)",
    sex_female_pct = 41,
    disease_state = "non-Hodgkin lymphoma: follicular lymphoma, diffuse large B-cell lymphoma, mantle cell lymphoma, other subtypes (peripheral T-cell and marginal zone lymphoma) and mixed-histology cohorts; both treatment-naive and relapsed/refractory",
    dose_range = "n/a -- no dose information was modelled. Li 2017 Discussion states that treatment duration, drug potency and dose 'were not considered in the model development, as there were not sufficient details available for many of the trials', so the treatment covariates are presence/absence indicators only",
    regions = "not reported; trials identified from PubMed, FDA reviews, clinicaltrials.gov, and ASCO / ASH proceedings, published 1993-2015",
    notes = paste0(
      "Trial selection followed the Cochrane Handbook and PRISMA: 513 NHL ",
      "trials were screened, 179 met the inclusion criteria, 17 were ",
      "removed as duplicates, and 50 of the remaining 162 were excluded ",
      "for not reporting a PFS Kaplan-Meier curve, leaving 112 studies / ",
      "155 cohorts. Trials with fewer than 25 subjects were excluded. ",
      "sex_female_pct is derived as 100 - 59, the complement of the ",
      "median 59% male reported in Table 1, and is a median across ",
      "cohorts rather than a pooled patient-level percentage. ",
      "Treatment coverage across the 155 cohorts (cohorts can appear in ",
      "several levels): rituximab 98 (63.2%), other drugs 60 (38.7%), ",
      "CHOP/CHOP-like 41 (26.4%), other chemotherapy 36 (23.2%), ",
      "bortezomib 22 (14.2%), bendamustine 11 (7%). ",
      "Fitted with NONMEM 7.2.0. Internal validation on 1,000 simulated ",
      "trials: the observed median PFS time fell inside the predicted ",
      "90% confidence interval for 92% of cohorts, and the predicted ",
      "hazard ratio was within 25% of the observed value for 75% of ",
      "studies and within 0.5- to 2-fold for 100% (median-PFS database ",
      "80 studies / 107 cohorts / 17 randomised studies; hazard-ratio ",
      "database 11 randomised studies / 24 cohorts). ",
      "A complete-covariate SENSITIVITY ANALYSIS removing the 21 studies ",
      "with imputed covariates (12 missing prior-rituximab, 6 missing ",
      "median age, 7 missing sex) is reported in the second half of ",
      "Table S1. It is deliberately NOT extracted as a separate model ",
      "per the library's replicate-author-structure policy, which ",
      "excludes robustness checks the authors did not report as final; ",
      "its estimates are tabulated in the vignette for reference. It ",
      "changed k and lambda by under 15% but did retain treatment and ",
      "tumor-subtype effects on k that the final model eliminated."
    )
  )

  ini({
    # ==================================================================
    # All values are Supplemental Content Table S1, "Final Model and
    # Sensitivity Analysis Parameter Estimates", FINAL MODEL column
    # (the Sensitivity Analysis column is not used -- see population
    # notes). The main text of Li 2017 contains no parameter table;
    # Results direct the reader to the Supplemental Content for both
    # the model description and Table S1.
    #
    # STRUCTURAL FORM. Supplemental Content Model Description, equation
    # (1), states that PFS time follows a Weibull distribution, that
    # lambda is the scale parameter and that k is the shape parameter
    # 'also known as the Weibull slope'. The equation as PUBLISHED
    # renders as ln(-ln(PFS(t))) = k*ln(t) - k*ln(t), i.e. identically
    # zero: the lambda in the second term was lost in typesetting (the
    # supplement ships the equation as a 241x19 px raster, so no
    # higher-resolution original exists; the .docx carries no OMML or
    # MathType object). The intended equation is the standard Weibull
    # probability-plot linearisation
    #   ln(-ln(PFS(t))) = k*ln(t) - k*ln(lambda),
    # equivalently PFS(t) = exp(-(t/lambda)^k). Three independent
    # checks fix this reading and exclude the mirror-image
    # k*ln(lambda) - k*ln(t):
    #   (a) the mirror image gives a survivor function INCREASING in t
    #       (PFS -> 1 as t -> infinity), which is inadmissible and
    #       contradicts the monotone-decreasing visual predictive check
    #       in Li 2017 Figure 2;
    #   (b) 'Weibull slope' names the coefficient of ln(t) in exactly
    #       this linearisation, so k must multiply +ln(t);
    #   (c) under this reading the median PFS time is
    #       lambda*log(2)^(1/k), directly proportional to lambda, which
    #       is what makes every Table S1 coefficient exponentiate to
    #       the median-PFS-time RATIO plotted in Figure 3 -- all nine
    #       reproduce to the printed precision (see the vignette source
    #       trace).
    # See the vignette Errata; nothing here is tuned or back-solved.
    # ==================================================================

    lk <- log(0.79)
    label("Weibull shape parameter k, the Weibull slope (log scale; unitless)") # Table S1 Final Model, row 'k' = 0.79 (RSE 4.3%, 95% CI 0.72, 0.85). Log-transformed here to keep k > 0; the source estimated it on the natural scale.

    llambda <- log(110)
    label("Weibull scale parameter lambda for the reference cohort (log scale; log-months)") # Table S1 Final Model, row 'lambda' = 110 (RSE 26.5%, 95% CI 52.8, 167.2). Months, per the Figure 2 visual-predictive-check x-axis 'Time(month)'.

    # ------------------------------------------------------------------
    # Covariate effects on log-lambda. All are log-linear (the
    # Supplemental Content Model Description: 'the exponential model was
    # centered by the median covariate value'), but the centering that
    # matters for reproduction is the stated REFERENCE POPULATION --
    # 'Patient population with rituximab as the only treatment, MCL as
    # the tumor subtype, and 100% treatment naive patients were the
    # reference population' -- not the database medians. Li 2017
    # Figure 3 is the arbiter: MCL and 'All naive patients' sit exactly
    # on the median-PFS-ratio = 1.0 line, which only holds under
    # reference-population centering.
    #
    # NO covariate was retained on k in the final model: Table S1 prints
    # '---' for every 'on k' row. The Supplemental Content records that
    # treatment on k and % male patients on k entered forward selection
    # and were removed in backward elimination.
    # ------------------------------------------------------------------

    e_norituximab_llambda <- -1.13
    label("Effect on log-lambda of the cohort NOT receiving rituximab (log scale; unitless)") # Table S1 Final Model, row 'No rituximab on lambda' = -1.13 (RSE 9.4%, 95% CI -1.34, -0.92). Applies to (1 - CONMED_RITUX): exp(1.13) = 3.10-fold longer median PFS with rituximab, the 'Rituximab' point in Figure 3A.

    e_chop_llambda <- 0.819
    label("Effect on log-lambda of a CHOP or CHOP-like regimen (log scale; unitless)") # Table S1 Final Model, row 'CHOP/CHOP-like on lambda' = 0.819 (RSE 26.7%, 95% CI 0.39, 1.25); exp(0.819) = 2.27-fold, the 'CHOP/CHOP-like' point in Figure 3A.

    e_nonchemoother_llambda <- 0.813
    label("Effect on log-lambda of any other non-chemotherapy drug in the regimen (log scale; unitless)") # Table S1 Final Model, row 'Other drugs on lambda' = 0.813 (RSE 8%, 95% CI 0.69, 0.94); exp(0.813) = 2.25-fold, the 'Other drugs' point in Figure 3A.

    e_bendamustine_llambda <- 1.39
    label("Effect on log-lambda of bendamustine in the regimen (log scale; unitless)") # Table S1 Final Model, row 'Bendamustine on lambda' = 1.39 (RSE 15.7%, 95% CI 0.96, 1.82); exp(1.39) = 4.01-fold, the largest treatment effect and the 'Bendamustine' point in Figure 3A.

    e_fl_llambda <- 0.46
    label("Effect on log-lambda per unit fraction of follicular-lymphoma patients in the cohort (log scale; unitless)") # Table S1 Final Model, row 'Follicular lymphoma on lambda' = 0.46 (RSE 25%, 95% CI 0.23, 0.68); exp(0.46) = 1.58, the 60% longer median PFS versus MCL quoted in the Abstract and plotted in Figure 3B.

    e_dlbcl_llambda <- -1.41
    label("Effect on log-lambda per unit fraction of DLBCL patients in the cohort (log scale; unitless)") # Table S1 Final Model, row 'DLBCL on lambda' = -1.41 (RSE 8.5%, 95% CI -1.64, -1.18); exp(-1.41) = 0.244, the '25% of MCL' median PFS quoted in the Abstract and plotted in Figure 3B.

    e_othernhl_llambda <- -0.17
    label("Effect on log-lambda per unit fraction of other-NHL patients in the cohort (log scale; unitless)") # Table S1 Final Model, row 'Other lymphomas on lambda' = -0.17 (RSE 113%, 95% CI -0.54, 0.21 -- brackets zero); exp(-0.17) = 0.84, reported as 'similar to the MCL patient population' (Figure 3B).

    e_experienced_llambda <- -1.85
    label("Effect on log-lambda per unit fraction of treatment-experienced patients in the cohort (log scale; unitless)") # Table S1 Final Model, row '% experienced patients on lambda' = -1.85 (RSE 12.8%, 95% CI -2.32, -1.39). Applies to (100 - LINE_1L_PCT)/100; exp(-1.85) = 0.157, the 'All experienced patients' point in Figure 3C.

    # ==================================================================
    # BETWEEN-STUDY variability. This is an MBMA: the random effects are
    # inter-STUDY, not inter-subject, so they are named eta_study_* per
    # the library's MBMA convention and must not be read as popPK BSV.
    # The Table S1 footnote reads 'omega^2_k is the variance of the
    # inter-study variability on k; omega^2_lambda is the variance of
    # the inter-study variability on lambda'.
    #
    # The footnote says VARIANCE, so the tabulated numbers are used
    # as-is with no CV back-transformation. The source does not state
    # whether the etas are additive or exponential; both are encoded
    # here as EXPONENTIAL (log-scale), which the magnitudes force:
    #   - lambda: an additive eta with SD sqrt(1.4) = 1.18 MONTHS on a
    #     110-month scale is a 1% between-study spread, which cannot
    #     coexist with the cohort-level 10th-to-90th-percentile spread
    #     in the Figure 2 visual predictive check (from ~0.02 to ~0.85
    #     progression-free at 30 months);
    #   - k: an additive eta with SD sqrt(0.17) = 0.41 on k = 0.79 puts
    #     ~3% of studies at k < 0, which is not an admissible Weibull
    #     shape.
    # Encoding both on the log scale also matches the exponential
    # covariate model the same parameters carry. See the vignette
    # Assumptions and deviations.
    # ==================================================================

    eta_study_lk ~ 0.17 # Table S1 Final Model, row 'omega^2_k' = 0.17 (RSE 13.4%, shrinkage 5%); variance on log-k, SD = 0.412
    eta_study_llambda ~ 1.4 # Table S1 Final Model, row 'omega^2_lambda' = 1.4 (RSE 14.5%, shrinkage 2.3%); variance on log-lambda, SD = 1.183

    # ==================================================================
    # RESIDUAL ERROR. Supplemental Content Model Description: 'residual
    # error was entered as the standard error of ln(-ln(PFS(t))), scaled
    # by an estimated constant', citing Arends 2008 on meta-analysis of
    # summary survival curves. The per-observation residual SD is
    # therefore expSd * SE_i, where SE_i is the standard error of the
    # cohort's own transformed Kaplan-Meier estimate at that time point
    # (a property of the digitised dataset -- it depends on the number
    # still at risk -- and not a model parameter; Li 2017 does not print
    # the formula it used).
    #
    # The value below is the UNWEIGHTED scalar. Per-observation
    # weighting by SE_i is left to downstream simulation code, which is
    # the library's ratified convention when the observation weight
    # touches ONLY the residual (see the N_ARM register entry and the
    # Mercier_2014_tramadol_tapentadol_mbma /
    # Chen_2025_methotrexate_*_mbma / Asiimwe_2025_trastuzumab*_mbma /
    # Hanan_2026_peginterferon_alfa_*_mbma precedent). Li 2017's
    # between-study etas are NOT SE-weighted, so no data column is
    # needed inside model().
    #
    # expSd (log-normal, '~ lnorm()') rather than addSd is the exact
    # encoding: a residual that is additive on ln(-ln(PFS)) = ln(cumhaz)
    # is by definition log-normal on cumhaz itself.
    # ==================================================================

    expSd <- 1.0954451
    label("Residual-error scalar, as an SD on the ln(-ln(PFS)) = ln(cumhaz) scale, for an observation whose own standard error is 1 (unitless); the per-observation SD is expSd * SE_i and the SE_i weighting is applied downstream, not here") # Table S1 Final Model, row 'sigma^2' (residual error scalar variance) = 1.2 (RSE 14.0%, shrinkage 3.4%); NONMEM variance -> SD = sqrt(1.2) = 1.0954451
  })

  model({
    # ------------------------------------------------------------------
    # Study-arm inputs supplied per row (all cohort-level aggregates,
    # not patient characteristics):
    #   CONMED_RITUX          -- 1 = regimen includes rituximab
    #   CONMED_CHOP           -- 1 = regimen includes CHOP / CHOP-like
    #   CONMED_BENDAMUSTINE   -- 1 = regimen includes bendamustine
    #   CONMED_NONCHEMO_OTHER -- 1 = regimen includes another
    #                            non-chemotherapy drug
    #   TUMTP_FL_PCT          -- % of the cohort with FL histology
    #   TUMTP_DLBCL_PCT       -- % of the cohort with DLBCL histology
    #   TUMTP_OTHER_NHL_PCT   -- % of the cohort with other-NHL histology
    #   LINE_1L_PCT           -- % of the cohort that is treatment-naive
    #
    # There is no ODE, no dosing and no PK layer: PFS is algebraic in
    # time. Every histology / line-of-therapy column is a PERCENT and
    # the source coefficients are per unit FRACTION, so each is divided
    # by 100 (the Franzese_2026_pdl1_nsclc_mbma convention).
    # ------------------------------------------------------------------

    # 1. Weibull shape for this cohort. No covariate was retained on k;
    #    only the between-study random effect acts here.
    k <- exp(lk + eta_study_lk)

    # 2. Weibull scale for this cohort, in months. The two INVERTED
    #    polarities are the easiest thing to get wrong in this model and
    #    are spelled out separately: Table S1 prints the coefficients
    #    for the ABSENCE of rituximab and for the fraction of
    #    treatment-EXPERIENCED patients, while the canonical columns
    #    carry rituximab PRESENCE and the treatment-NAIVE fraction.
    noRituximab <- 1 - CONMED_RITUX
    fracExperienced <- (100 - LINE_1L_PCT) / 100

    covLambda <-
      e_norituximab_llambda * noRituximab +
      e_chop_llambda * CONMED_CHOP +
      e_bendamustine_llambda * CONMED_BENDAMUSTINE +
      e_nonchemoother_llambda * CONMED_NONCHEMO_OTHER +
      e_fl_llambda * (TUMTP_FL_PCT / 100) +
      e_dlbcl_llambda * (TUMTP_DLBCL_PCT / 100) +
      e_othernhl_llambda * (TUMTP_OTHER_NHL_PCT / 100) +
      e_experienced_llambda * fracExperienced

    lambda <- exp(llambda + covLambda + eta_study_llambda)

    # 3. Cumulative hazard of progression or death, (t/lambda)^k. This
    #    is the observed quantity's model prediction: cumhaz = -ln(PFS),
    #    so ln(cumhaz) is exactly the ln(-ln(PFS(t))) scale the source
    #    fitted on.
    #
    #    NO epsilon offset is added to t. It is tempting to write
    #    ((t + 1e-6)/lambda)^k so that the log below stays finite at
    #    t = 0, but that is wrong here: a small offset raised to a SMALL
    #    POWER is not small. With k as low as 0.24 (well inside the
    #    between-study distribution, k = 0.79 * exp(eta) with variance
    #    0.17) and lambda = 10 months, (1e-6/10)^0.24 = 0.021, i.e. a 2%
    #    drop in the survivor function at time zero. The exact form
    #    gives 0^k = 0 and sur(0) = 1 identically for any k > 0.
    #
    #    The consequence is that cumhaz = 0 at t = 0, so the observation
    #    scale ln(cumhaz) is -Inf there. That is faithful rather than a
    #    defect: ln(-ln(PFS)) is undefined wherever PFS = 1, so the
    #    source's own transformed dataset cannot contain a t = 0 record
    #    either. Place observation records at t > 0.
    cumhaz <- (t / lambda)^k

    # 4. Progression-free probability, the quantity Li 2017 plots as
    #    '% PFS' in the Figure 2 visual predictive check.
    sur <- exp(-cumhaz)

    # 5. Median PFS time in months, lambda * log(2)^(1/k). Directly
    #    proportional to lambda, which is why exponentiating any
    #    Table S1 coefficient gives the median-PFS-time ratio plotted
    #    in Figure 3. For the reference cohort:
    #    110 * log(2)^(1/0.79) = 69.1 months.
    tmed <- lambda * log(2)^(1 / k)

    # 6. Observation. The residual is additive on
    #    ln(-ln(PFS)) = ln(cumhaz), which is a log-normal residual on
    #    cumhaz itself, so '~ lnorm(expSd)' is the source error model
    #    rather than an approximation of it. expSd is the UNWEIGHTED
    #    scalar; multiply it by the observation's own standard error
    #    SE_i downstream to reproduce the fitted weighting.
    cumhaz ~ lnorm(expSd)
  })
}
