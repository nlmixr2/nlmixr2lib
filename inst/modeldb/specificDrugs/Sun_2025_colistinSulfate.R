Sun_2025_colistinSulfate <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous colistin sulfate in",
    "critically ill adults with carbapenem-resistant organism infections",
    "(Sun 2025; n = 178 Chinese ICU patients, 364 sparse therapeutic-drug-",
    "monitoring plasma samples spanning 0.16-4.91 mg/L; a further 26 patients",
    "and 57 samples were held out for external validation). Linear elimination",
    "from the central compartment with intravenous-infusion input.",
    "Cockcroft-Gault creatinine clearance enters clearance as a power function",
    "centred on 71.40 mL/min (exponent 0.456) and body weight enters the",
    "central volume as a power function centred on 67.89 kg (exponent 1.2);",
    "age, sex and albumin were screened but not retained. Exponential",
    "inter-individual variability was estimated on all four structural",
    "parameters. The residual-error MAGNITUDE is not reported by the paper and",
    "is fixed at zero here (see the ini() comments and the vignette Errata).",
    "Colistin sulfate is administered as the active drug and must not be",
    "confused with colistimethate sodium (CMS), the inactive prodrug modelled",
    "in Plachouras 2009, Mohamed 2012, Jacobs 2016 and Karaiskos 2015."
  )
  reference <- paste(
    "Sun Q, Li X, Wang G, Wang X, Xing B, Xun Z, Lu N, Li Z (2025).",
    "Population pharmacokinetics of colistin sulfate in critically ill",
    "patients based on NONMEM.",
    "Sci Rep 15:18295.",
    "doi:10.1038/s41598-025-03503-9.",
    sep = " "
  )
  vignette <- "Sun_2025_colistinSulfate"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    # Methods "Blood sample collection and determination of plasma colistin
    # sulfate concentration": samples were centrifuged to separate PLASMA,
    # stored at -80 C and assayed by a validated LC-MS/MS method. The assayed
    # analyte is colistin itself, reported as "the total concentration of
    # colistin sulfate ... calculated as the sum of colistin A and colistin
    # B" -- not a colistimethate metabolite.
    central     = list(analyte = "colistin sulfate", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "colistin sulfate", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = paste(
        "Creatinine clearance estimated with the Cockcroft-Gault equation,",
        "reported as RAW mL/min and NOT normalised to 1.73 m^2 body surface",
        "area."
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject (a single baseline value). Power effect on CL:",
        "CL = 2.66 * (CrCL / 71.40)^0.456 * exp(eta) per Sun 2025 Eq. 1, with",
        "the exponent 0.456 also tabulated as the 'dCLdCrCL' row of Table 2",
        "(RSE 8.6%, 95% CI 0.379-0.533; bootstrap median 0.466, 95% CI",
        "0.39-0.59). The centring constant 71.40 mL/min appears ONLY inside",
        "Eq. 1 -- it is not tabulated, and it is NOT the cohort median (Table",
        "1 gives CrCL 51.9 [6.7, 271.8] mL/min); it is most consistent with",
        "the cohort arithmetic mean, which the paper does not report. Unlike",
        "most papers the Cockcroft-Gault form IS written out, as Eq. 5:",
        "CrCL (mL/min) = (140 - AGE) * WT (kg) / (0.818 * Scr (umol/L)), with",
        "the result multiplied by 0.85 for females (Clinical data collection).",
        "The 0.818 denominator is the standard 72 * Scr(mg/dL) form rewritten",
        "for Scr in umol/L (72 / 88.4 = 0.814). Fitted over an observed range",
        "of 6.7-271.8 mL/min; the paper's Monte Carlo simulations span 10-120",
        "mL/min. Patients receiving continuous renal replacement therapy",
        "during colistin sulfate treatment were EXCLUDED, so the model",
        "carries no information about renal replacement therapy."
      ),
      source_name        = "CrCL"
    ),
    WT = list(
      description        = "Total body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Power effect on the central volume:",
        "V1 = 49.70 * (WT / 67.89)^1.2 * exp(eta) per Sun 2025 Eq. 2, with the",
        "exponent 1.2 also tabulated as the 'dVdWT' row of Table 2 (legend",
        "'dV1 dWT, exponential parameter coefficient of WT to V1'; RSE 22.2%,",
        "95% CI 0.679-1.721; bootstrap median 1.28, 95% CI 0.60-1.75). As with",
        "CRCL the centring constant 67.89 kg appears ONLY inside Eq. 2 and is",
        "not the cohort median (Table 1 gives WT 70 [40, 100] kg); it is most",
        "consistent with the unreported cohort arithmetic mean. NOTE that 1.2",
        "is NOT the allometric 1.0 for a volume -- it is an estimated",
        "exponent, and its confidence interval (0.679-1.721) contains 1.0, so",
        "the data do not exclude simple proportionality. The exponent is not",
        "shared with clearance: CL carries no weight term at all in this",
        "model, which is why this is `e_wt_vc` and not a `wt_<param>`",
        "allometric pair. Weight also enters CRCL through the Cockcroft-Gault",
        "equation (Eq. 5), so body size affects CL indirectly via the renal",
        "term; the paper does not discuss this partial collinearity.",
        "Fitted over an observed range of 40-100 kg."
      ),
      source_name        = "WT"
    )
  )

  # Screened during stepwise covariate modelling but not retained in the final
  # model. Results "PPK analysis": "The final results indicated that age, sex,
  # albumin, and other clinical variables had no statistically significant
  # correlation with the pharmacokinetic (PK) parameters." The paper reports
  # neither a forward-inclusion dOFV nor any coefficient for these, so nothing
  # is encoded. The remaining screened laboratory variables (urea, ALT, AST,
  # prothrombin time) are named only in the Table 1 demographics and are not
  # individually reported as covariate candidates, so they are not listed here.
  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Age at baseline.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened in the stepwise covariate model and rejected (Results, 'PPK",
        "analysis'). Table 1 gives a median of 78.5 [31, 98] years -- an",
        "unusually elderly ICU cohort. No point estimate is reported, so",
        "nothing is encoded. Age DOES enter the model indirectly, through the",
        "Cockcroft-Gault equation that generates CRCL (Eq. 5)."
      ),
      source_name        = "Age"
    ),
    SEXF = list(
      description        = "Female sex indicator.",
      units              = "unitless",
      type               = "categorical",
      reference_category = "male",
      notes              = paste(
        "Screened in the stepwise covariate model and rejected (Results, 'PPK",
        "analysis'). Table 1 gives 117 male / 61 female. No point estimate is",
        "reported, so nothing is encoded. Sex DOES enter the model indirectly,",
        "through the 0.85 female multiplier of the Cockcroft-Gault equation",
        "that generates CRCL (Eq. 5)."
      ),
      source_name        = "Gender"
    ),
    ALB = list(
      description        = "Serum albumin.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened in the stepwise covariate model and rejected (Results, 'PPK",
        "analysis'). Table 1 gives 30.97 +/- 4.43 g/L, i.e. the cohort was",
        "uniformly hypoalbuminaemic with little spread, which plausibly limits",
        "the power to detect an albumin effect. No point estimate is reported,",
        "so nothing is encoded."
      ),
      source_name        = "ALB"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 178L,
    n_studies      = 1L,
    n_observations = 364L,
    age_range      = "31-98 years (inclusion criterion age >= 18 years)",
    age_median     = "78.5 years",
    weight_range   = "40-100 kg",
    weight_median  = "70 kg",
    sex_female_pct = 100 * 61 / 178,
    race_ethnicity = c(Asian = 100),
    disease_state  = paste(
      "Critically ill intensive-care patients treated with intravenous",
      "colistin sulfate for carbapenem-resistant organism (CRO) infections.",
      "Infection sites were pulmonary (93.26%), urinary tract (10.11%),",
      "bloodstream (3.37%), abdomen (2.81%) and endocardium (0.56%).",
      "Isolates were Klebsiella pneumoniae (51.69%), Acinetobacter baumannii",
      "(39.32%), Pseudomonas aeruginosa (8.43%) and Escherichia coli",
      "(3.93%); 13.48% were treated empirically with no isolate. Almost all",
      "patients received a concomitant antibacterial, most often meropenem",
      "(61.80%), tigecycline (24.72%), fosfomycin (10.67%) or",
      "cefoperazone-sulbactam (9.55%). The 14-day mortality rate was 27.53%",
      "and 11.8% developed acute kidney injury by RIFLE criteria (urinary",
      "output criterion excluded). Patients who died within 24 h of the first",
      "dose, who received nebulised or intrathecal colistin sulfate, or who",
      "received continuous renal replacement therapy during treatment were",
      "excluded."
    ),
    dose_range     = paste(
      "Intravenous infusion of colistin sulfate (marketed specification 0.5",
      "MU; Asia Pioneer Pharmaceutical, Shanghai). Table 1 gives a median",
      "daily dose of 1.5 [1.0, 4.0] MU. The product label recommends 1.0-1.5",
      "MU/day divided into 2-3 maintenance doses; the study hospital",
      "additionally recommended a 1.0-2.0 MU loading dose. Treatment lasted",
      "11.21 +/- 6.39 days. IMPORTANT: the paper expresses every dose in",
      "million international units (MU) and NEVER states an MU-to-mg",
      "conversion, while its parameters (CL in L/h, V in L) and observations",
      "(mg/L) are on a mass basis. This model therefore takes dose in mg. See",
      "`notes` below and the vignette Errata for the back-solved conversion."
    ),
    sampling       = paste(
      "Sparse, opportunistic therapeutic drug monitoring retrieved",
      "retrospectively from computerised records; the precise dose and",
      "sampling clock times were recorded. Of the 364 concentrations, 52",
      "(14.29%) were first-dose peaks, 154 (42.31%) steady-state peaks and",
      "158 (43.41%) troughs. The paper does not state the nominal offsets at",
      "which 'peak' and 'trough' were drawn. The Discussion flags the small",
      "first-dose-peak fraction as a limitation: 'there was insufficient",
      "information regarding the initial rapid distribution phase'."
    ),
    renal_function = paste(
      "Cockcroft-Gault creatinine clearance 51.9 [6.7, 271.8] mL/min, raw and",
      "not BSA-normalised (Table 1). Baseline serum creatinine 79 [19, 753]",
      "umol/L. Patients on continuous renal replacement therapy during",
      "treatment were excluded."
    ),
    regions        = "People's Republic of China (single centre; Beijing Electric Power Hospital, Beijing).",
    notes          = paste(
      "Baseline demographics from Sun 2025 Table 1 and Results 'Baseline",
      "characteristics of patients'. Retrospective single-centre cohort",
      "collected May 2022 to May 2024 (ethics approval n090327). The final",
      "model was fit in NONMEM 7.5 via Pirana 3.0 and evaluated with a",
      "1000-sample nonparametric bootstrap (Table 2), a prediction- and",
      "variability-corrected VPC (Supplementary Figure S1) and an external",
      "validation in 26 further patients / 57 concentrations (Supplementary",
      "Tables S1-S3).",
      "",
      "DOSE UNITS. The paper's dosing is entirely in MU and no MU-to-mg",
      "potency is given anywhere in the article or supplement, so the",
      "conversion had to be back-solved from the paper's own output. Fitting",
      "the 18 digitised hourly-average peak and trough values of Figure 3",
      "panel A (1.0 MU q8h, 1 h infusion) with this model gives 1 MU ~ 46 mg",
      "(RMS log-error 0.07 over a 5-fold concentration range; the fit is",
      "flat between roughly 42 and 48 mg/MU, which is finer than the",
      "digitisation supports). Two independent routes agree: repeating the",
      "paper's external validation with typical-value predictions against",
      "the 57 observed concentrations of Supplementary Table S2 minimises at",
      "~45 mg/MU, and the sibling model modellib('Ma_2026_colistinSulfate')",
      "recovers 44 mg per 10^6 IU from an unrelated cohort by an unrelated",
      "route. That range is consistent with the pharmacopoeial potency of the",
      "colistin SULFATE salt (~22,000 IU/mg) rather than with colistin BASE",
      "activity (30,000 IU/mg, i.e. 33.3 mg/MU); the distinction matters",
      "because the two differ by ~35%. NOTE that the paper's Figure 2",
      "(probability of target attainment) cannot be reconciled with this",
      "conversion under a steady-state identity -- but Figure 2 refutes that",
      "reading itself, since it draws four regimens sharing one maintenance",
      "dose as four distinct curves, which is only possible if its AUC0-24 is",
      "a FIRST-24-hour quantity. The conversion is NOT encoded as a model",
      "parameter -- it is a property of the marketed product, not of the PK",
      "-- so a user supplying mg-denominated doses is unaffected by any",
      "residual uncertainty in it. The vignette derives it and shows every",
      "supporting fit.",
      "",
      "The unbound fraction of 0.5 and the fAUC(0-24)/MIC >= 10 target used",
      "for the paper's probability-of-target-attainment simulations are",
      "literature assumptions, not fitted parameters, and are therefore not",
      "encoded here.",
      "",
      "Two internal text-versus-table inconsistencies in the demographics do",
      "not affect the model: Results names 'ceftazidime-avibactam (9.55%)'",
      "where Table 1 has 'Cefoperazone/sulbactam 17 (9.55%)', and Results",
      "names 'Enterobacter cloacae (3.93%)' where Table 1 has 'E. coli 7",
      "(3.93%)'. Table 1 is followed above in both cases."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters. Typical values refer to the covariate
    # reference subject: CrCL 71.40 mL/min and WT 67.89 kg (the centring
    # constants of Sun 2025 Eqs. 1-2).
    # ------------------------------------------------------------------

    # ERRATUM. Eq. 1 as typeset reads "CL(L/h) = 266 (CrCL/71.40)^0.456
    # exp(etaCL)" -- the decimal point of the leading coefficient is
    # missing. Table 2 gives tvCL = 2.66 L/h (RSE 7%, 95% CI 2.293-3.027;
    # bootstrap median 2.67, 95% CI 2.30-3.04), so 2.66 is the estimate and
    # the "266" in Eq. 1 is a typographical error -- the table carries a
    # standard error and a bootstrap interval and the equation does not, and
    # 2.66 is the value the rest of the paper is consistent with. Arithmetic
    # settles it independently: since Css,avg = daily dose / (CL * 24), a CL
    # of 266 L/h would put the typical steady-state concentration a
    # hundredfold below the 0.93 +/- 0.35 mg/L that Table 1 reports and two
    # orders of magnitude below the 0.16-4.91 mg/L range of the 364 assayed
    # samples -- i.e. below the assay's own lower limit of quantification.
    # (Note that the usual "the half-life would be absurd" argument does NOT
    # work here: because Q = 1.63 L/h and V2 = 109 L make the peripheral
    # compartment rate-limiting, CL = 266 L/h still yields a terminal
    # half-life of ~47 h. Only the exposure argument discriminates.)
    lcl <- log(2.66);  label("Clearance (CL, L/h) at CrCL 71.40 mL/min")     # Table 2, tvCL 2.66 (RSE 7%; bootstrap median 2.67, 95% CI 2.30-3.04); Eq. 1 leading coefficient, mis-typeset as "266"
    lvc <- log(49.70); label("Central volume of distribution (V1, L) at WT 67.89 kg")  # Table 2, tvV1 49.7 (RSE 5.8%, 95% CI 44.036-55.364; bootstrap median 50.04, 95% CI 44.16-55.26); also the Eq. 2 leading coefficient
    lq  <- log(1.63);  label("Inter-compartmental clearance (L/h)")       # Table 2, tvQ 1.63 (RSE 19.8%, 95% CI 0.999-2.261; bootstrap median 1.67, 95% CI 0.44-3.67); also Eq. 3
    lvp <- log(109);   label("Peripheral volume of distribution (L)")    # Table 2, tvV2 109 (RSE 24.8%, 95% CI 56.08-161.92; bootstrap median 103.65, 95% CI 58.18-161.31); also Eq. 4

    # CAVEAT on the distribution parameters. Q = 1.63 L/h against V2 = 109 L
    # makes the peripheral compartment equilibrate very slowly, so these
    # values imply a terminal half-life of about 80 h (Vss = 159 L) and
    # continued accumulation for several days -- much longer than the 9-18 h
    # usually reported for colistin. The paper's own Discussion names the
    # cause: only 14.29% of the 364 samples were first-dose peaks, so "there
    # was insufficient information regarding the initial rapid distribution
    # phase of the drug in the body", and Q and V2 are correspondingly the
    # least precisely estimated parameters in Table 2 (RSE 19.8% and 24.8%;
    # the tvQ bootstrap interval 0.44-3.67 spans more than eightfold).
    # Neither parameter is altered here -- the published values are packaged
    # as published -- but users simulating beyond a few days should treat the
    # terminal phase as poorly identified. The vignette quantifies this.

    # ------------------------------------------------------------------
    # Covariate effects. Both are power functions of the covariate divided
    # by a centring constant that appears only inside the printed equation.
    # ------------------------------------------------------------------
    e_crcl_cl <- 0.456; label("Power exponent of CrCL on CL (unitless)")  # Table 2, "dCLdCrCL" 0.456 (RSE 8.6%, 95% CI 0.379-0.533; bootstrap median 0.466, 95% CI 0.39-0.59); Eq. 1 exponent
    e_wt_vc   <- 1.2;   label("Power exponent of WT on V1 (unitless)")    # Table 2, "dVdWT" 1.2 (RSE 22.2%, 95% CI 0.679-1.721; bootstrap median 1.28, 95% CI 0.60-1.75); Eq. 2 exponent

    # ------------------------------------------------------------------
    # Inter-individual variability, exponential on all four structural
    # parameters (Methods "PPK models of colistin sulfate": "Between-subject
    # variability (BSV) was assessed using an exponential function"; Eqs. 1-4
    # each carry an exp(eta) term).
    #
    # SCALE. The Table 2 legend states the scale explicitly and unambiguously
    # -- "omega^2 V1, VARIANCE of inter-individual variability for V1" (and
    # likewise for CL, Q and V2) -- so the tabulated numbers are variances and
    # are encoded here as-is, with no squaring. The implied CVs are 36% for
    # CL, 52% for V1, 190% for Q and 100% for V2. The two large ones are
    # consistent with their own reported precision and shrinkage: Q and V2 are
    # the parameters this sparse, trough-and-peak-dominated design constrains
    # worst (tvQ RSE 19.8%, tvV2 RSE 24.8%; omega^2 RSEs 49% and 29%;
    # eta-shrinkage 60.7% and 57.6%), and the Discussion concedes that with
    # only 14.29% first-dose peaks "there was insufficient information
    # regarding the initial rapid distribution phase of the drug in the body"
    # -- which is precisely the phase that identifies Q and V2. The vignette
    # cross-checks the variance reading against Figure 3 panel A, whose
    # hourly ARITHMETIC-mean curve sits above the typical-value curve by an
    # amount that depends on these variances.
    #
    # No off-diagonal covariances are reported, so the block is diagonal.
    # ------------------------------------------------------------------
    etalcl ~ 0.125  # Table 2, omega^2 CL 0.125 (RSE 22%, eta-shrinkage 38.5%; bootstrap median 0.117, 95% CI 0.05-0.19)
    etalvc ~ 0.244  # Table 2, omega^2 V1 0.244 (RSE 19%, eta-shrinkage 21.9%; bootstrap median 0.215, 95% CI 0.16-0.34)
    etalq  ~ 1.69   # Table 2, omega^2 Q  1.69  (RSE 49%, eta-shrinkage 60.7%; bootstrap median 1.657, 95% CI 0.12-6.00)
    etalvp ~ 0.993  # Table 2, omega^2 V2 0.993 (RSE 29%, eta-shrinkage 57.6%; bootstrap median 0.892, 95% CI 0.10-4.00)

    # ------------------------------------------------------------------
    # Residual error: COMBINED additive plus proportional in FORM, but with
    # both MAGNITUDES unreported, so both are fixed at zero.
    #
    # The form is stated twice. Results "PPK analysis": "Residual random
    # effects were assessed using additive plus proportional error models."
    # Methods: "Within-subject variability (WSV) was assessed by additive,
    # proportional, or combined (additive plus proportional) models".
    #
    # The magnitudes are not. Table 2's entire "Residual variability (sigma)"
    # block is a single row, "stdev0", with the value 1 and -- alone among
    # every row of the table -- no RSE and no confidence interval, only an
    # epsilon-shrinkage of 40.8%. That is the signature of the standard
    # NONMEM combined-error parameterisation
    #     W = SQRT(THETA(a)**2 + (THETA(p)*IPRED)**2)
    #     Y = IPRED + W*EPS(1)   with   $SIGMA 1 FIX
    # in which the reported "1" is the FIXED variance scale of EPS(1) and the
    # actual additive and proportional coefficients live in $ERROR as THETAs
    # that this paper does not tabulate. The same parameterisation is
    # documented for the sibling colistin sulfate model
    # `modellib('Ma_2026_colistinSulfate')`, whose paper DOES tabulate its
    # proportional coefficient.
    #
    # Reading the 1 as a real residual SD is refuted either way it is taken:
    # as an additive SD of 1 mg/L it would exceed the cohort's mean
    # steady-state concentration (0.93 +/- 0.35 mg/L, Table 1) and most of
    # the 0.16-4.91 mg/L observation range; as a proportional SD it would be
    # 100% CV, against which the paper's own external validation reports a
    # relative RMSE of 7.90% (Supplementary Table S3).
    #
    # Per the standing policy for unreported variance components, nothing is
    # invented: both terms are fixed at 0 so the model simulates the
    # typical-value / IIV-only trajectories the paper plots. The vignette
    # Errata records an order-of-magnitude reconstruction (roughly 13%
    # proportional, from the Supplementary Table S3 relative RMSE of 7.90%
    # de-shrunk by the reported 40.8% epsilon-shrinkage) for users who need a
    # non-zero residual, and notes that the near-equality of the peak (7.06%)
    # and trough (7.80%) RMSEs implies the additive term is negligible over
    # the observed concentration range.
    # ------------------------------------------------------------------
    propSd <- fixed(0); label("Proportional residual SD (fraction; 0 -- magnitude not reported in the source)")  # Table 2 reports only the FIXED sigma scale (stdev0 = 1), not the $ERROR coefficients
    addSd  <- fixed(0); label("Additive residual SD (mg/L; 0 -- magnitude not reported in the source)")          # Table 2 reports only the FIXED sigma scale (stdev0 = 1), not the $ERROR coefficients
  })

  model({
    # Individual PK parameters. Creatinine clearance enters CL and body
    # weight enters V1, each as a power function centred on the constant
    # printed in the corresponding equation (Sun 2025 Eqs. 1-2). Q and V2
    # carry inter-individual variability but no covariate (Eqs. 3-4).
    cl <- exp(lcl + etalcl) * (CRCL / 71.40)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 67.89)^e_wt_vc
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    # Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Colistin sulfate is given as an intravenous infusion, so dose records
    # target `central` directly (with a rate or duration) and there is no
    # depot compartment. The paper's simulated regimens use 1 h infusions,
    # or 3 h for the four regimens labelled "3h" in Figs. 2 and 4.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # Concentration in mg/L (= ug/mL), the assay's reporting unit.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
