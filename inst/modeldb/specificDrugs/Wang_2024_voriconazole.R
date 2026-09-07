Wang_2024_voriconazole <- function() {
  description <- "Two-compartment population pharmacokinetic model with first-order absorption and linear elimination for voriconazole in critically ill adults in a respiratory intensive care unit (Wang 2024); clearance carries five covariates - quick C-reactive protein, creatinine clearance, continuous renal replacement therapy, platelet count and prothrombin time - with continuous renal replacement therapy raising clearance 1.46-fold, and the absorption rate constant fixed to a published literature value"
  reference <- "Wang Y, Ye Q, Li P, Huang L, Qi Z, Chen W, Zhan Q, Wang C. Renal Replacement Therapy as a New Indicator of Voriconazole Clearance in a Population Pharmacokinetic Analysis of Critically Ill Patients. Pharmaceuticals (Basel). 2024;17(6):665. doi:10.3390/ph17060665"
  vignette <- "Wang_2024_voriconazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Wang 2024 Section 4.2 (voriconazole
  # given by intravenous infusion, orally, or by nasogastric tube) and
  # Section 4.3 (plasma UPLC-MS/MS assay, LLOQ 0.097 mg/L).
  compartmentData <- list(
    depot       = list(analyte = "voriconazole", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "voriconazole", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "voriconazole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRP = list(
      description        = "Quick C-reactive protein (qCRP), a point-of-care C-reactive protein assay",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying: measured on the day of each blood collection (Wang 2024 Section 4.2). Enters CL as (CRP/73.6)^-0.142; the 73.6 mg/L reference is the cohort median reported in Wang 2024 Table 1 and is reproduced exactly as the divisor of the final-model equation in Section 2.2. The negative exponent means clearance falls as inflammation rises, which the paper attributes (Discussion Section 3.2) to inflammation-driven downregulation of CYP2C19/CYP3A4. The source calls the assay 'quick CRP' (qCRP); it is recorded here under the general-scope CRP canonical, which explicitly spans assay variants and time-varying use.",
      source_name        = "qCRP"
    ),
    CRCL = list(
      description        = "Creatinine clearance by the Cockcroft-Gault equation, raw and NOT body-surface-area normalized",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Wang 2024 Section 4.2 states explicitly that 'CLCR was calculated using the Cockcroft-Gault equation', so this column is a raw mL/min value and must NOT be supplied on the BSA-normalized mL/min/1.73m2 scale that is the CRCL canonical's default. Enters CL as (CRCL/71.8)^0.218. The 71.8 mL/min divisor appears ONLY inside the final-model equation in Section 2.2 and does not equal the Table 1 cohort median of 68.5 mL/min; the printed equation is used, per the register's established precedent for this situation. Time-varying, measured on the day of each blood collection. Wang 2024 Discussion Section 3.4 calls the positive exponent 'surprising' given that under 2% of a voriconazole dose is renally excreted, and speculates that CLCR is partly a proxy for the extracorporeal clearance contributed by CRRT.",
      source_name        = "CLCR"
    ),
    RRT_CRRT_ACTIVE = list(
      description        = "Continuous renal replacement therapy running at the time of the record (1) or not (0)",
      units              = "binary",
      type               = "binary",
      reference_category = "0 (no CRRT running)",
      notes              = "Record-level and time-varying: CRRT was recorded per concentration, 185 of the 746 concentrations and 122 of the 501 on-machine occasions being on CRRT (Wang 2024 Section 2.1 and Table 1). Every CRRT patient in this cohort underwent continuous veno-venous hemofiltration (CVVH) with blood flow 120-150 mL/min, replacement-fluid rate 25-30 mL/kg/h and predilution, on Fresenius machines (Discussion Section 3.3) - a continuous modality, which is why RRT_CRRT_ACTIVE is used rather than the intermittent-hemodialysis counterpart. Enters CL multiplicatively as 1.46^CRRT, encoded here as the log-additive shift e_crrt_cl = log(1.46). This is the paper's headline finding: it contradicts the conventional view that voriconazole clearance is unaffected by renal replacement.",
      source_name        = "CRRT"
    ),
    PLT = list(
      description        = "Platelet count",
      units              = "10^9 cells/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying, from the routine blood examination drawn on the day of each blood collection. Enters CL as (PLT/144)^0.166. The 144 x 10^9/L divisor appears only inside the final-model equation in Section 2.2 and does not equal the Table 1 cohort median of 150.5 x 10^9/L; the printed equation is used, matching the precedent already recorded under this canonical for Stitt 2026 (equation divisor 196 against a table median of 197). Wang 2024 Discussion Section 3.5 reads the positive exponent as a liver-function marker rather than a platelet-mediated mechanism, thrombocytopenia accompanying hepatic dysfunction in this cohort.",
      source_name        = "PLT"
    ),
    PT_SEC = list(
      description        = "Prothrombin time, raw laboratory value in seconds",
      units              = "seconds",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying, measured on the day of each blood collection. Enters CL as (PT_SEC/15)^-0.875, much the largest of the four continuous covariate exponents. The 15 s reference exists ONLY inside the final-model equation in Section 2.2: Wang 2024 does not tabulate prothrombin time anywhere, so no published cohort distribution is available for simulation and the vignette holds this covariate at the 15 s reference. The negative exponent is the expected direction, a longer prothrombin time marking worse hepatic synthetic function and hence slower clearance of a hepatically metabolised triazole; the paper supports the interpretation with a Supplementary Figure S2 PT-versus-AST scatter. The canonical name carries the _SEC unit suffix because the bare token PT is already recorded in the register with an unrelated meaning (patient-versus-healthy-volunteer indicator, under DIS_GERD) and because this raw-seconds value must not be confused with PTR (ratio to the subject's own baseline) or INR_BASE (unitless INR); naming ratified by operator decision (sidecar request 001).",
      source_name        = "PT"
    )
  )

  covariatesDataExcluded <- list(
    ECMO_STATUS = list(
      description = "Extracorporeal membrane oxygenation in use",
      units       = "binary",
      type        = "binary",
      notes       = "Screened as a categorical covariate (Wang 2024 Section 4.4.3) and analysed as a stratifying factor in Table 4, but NOT retained on any pharmacokinetic parameter in the final model. Wang 2024 found no significant difference in CL, Vc, Vp, AUC24, Cmin or the reported half-life between the ECMO and non-ECMO groups (all p > 0.6). This is a substantive negative result for the paper - the title's contrast is that renal replacement therapy matters where extracorporeal membrane oxygenation does not - so the screen is recorded rather than dropped."
    ),
    AST = list(
      description = "Aspartate transaminase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Entered the model during forward selection (delta OFV 10.248) but was removed from the final model: Wang 2024 Section 2.2 states it 'had a poor relative standard error (RSE) (77%) and low estimate value (0.08)'. No usable point estimate is therefore published for it."
    ),
    ALT = list(
      description = "Alanine transaminase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Listed among the continuous covariates examined (Wang 2024 Section 4.4.3) and tabulated in Table 1, but not retained in the final model."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Listed among the continuous covariates examined (Wang 2024 Section 4.4.3) and tabulated in Table 1, but not retained in the final model."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Listed among the continuous covariates examined (Wang 2024 Section 4.4.3) and tabulated in Table 1 (median 34.0, under the header 'Albumin (mg/dL)' - a unit typo, since 34 g/L is the plausible value and 34 mg/dL is not), but not retained in the final model. Wang 2024 Discussion Section 3.3 nevertheless notes that 50.8% of the CRRT subgroup were hypoalbuminaemic."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Listed among the continuous covariates examined (Wang 2024 Section 4.4.3; cohort mean 65.3 kg), but not retained: this model carries NO allometric or other body-size term on any parameter, so clearance and both volumes are absolute population values rather than per-70-kg values."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Listed among the continuous covariates examined (Wang 2024 Section 4.4.3; cohort mean 64 years), but not retained in the final model."
    ),
    SEXF = list(
      description = "Female sex",
      units       = "binary",
      type        = "binary",
      notes       = "Screened as a categorical covariate (Wang 2024 Section 4.4.3; 287 of 408 participants, 70.3%, were men) but not retained in the final model."
    ),
    APACHE_II = list(
      description = "Acute Physiology and Chronic Health Evaluation II score",
      units       = "points",
      type        = "continuous",
      notes       = "Screened as a continuous covariate (Wang 2024 Section 4.4.3; median 19.0, IQR 14.0-25.0 on the day of blood collection) but not retained in the final model. The companion Sequential Organ Failure Assessment (SOFA) score, median 7.0, IQR 4.0-10.0, was screened in the same step and likewise not retained; it is recorded here rather than as its own entry because the register carries no SOFA canonical and this model uses neither score."
    ),
    CONMED_PPI = list(
      description = "Concomitant proton pump inhibitor",
      units       = "binary",
      type        = "binary",
      notes       = "Screened as a categorical covariate (Wang 2024 Section 4.4.3; 353 of 501 occasions, 70.5%) but not retained in the final model, despite the well-known CYP2C19 interaction between omeprazole-class agents and voriconazole."
    ),
    CONMED_STEROID = list(
      description = "Concomitant systemic glucocorticoid use",
      units       = "binary",
      type        = "binary",
      notes       = "Screened as a categorical covariate (Wang 2024 Section 4.4.3 names 'co-medications such as proton pump inhibitors and glucocorticoids'; Table 1 reports use in 197 of 501 occasions, 39.3%) but not retained in the final model. Recorded alongside CONMED_PPI because the paper screened the two co-medication classes in the same step, and because glucocorticoids are CYP3A4 inducers with a documented voriconazole interaction, so the negative screen is informative rather than incidental."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 408L,
    n_studies      = 1L,
    n_observations = 746L,
    age_median     = "64 years (mean)",
    weight_median  = "65.3 kg (mean)",
    sex_female_pct = 29.7,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Critically ill adults in a respiratory intensive care unit, every patient carrying either mild or severe lung infection, receiving voriconazole for suspected or documented invasive fungal infection. 104 patients (185 concentrations; 122 of 501 on-machine occasions) received continuous renal replacement therapy, uniformly as continuous veno-venous hemofiltration; 85 patients (154 concentrations) received extracorporeal membrane oxygenation, some concomitantly with CRRT.",
    dose_range     = "Voriconazole 200 mg every 12 h in 342 patients (83.8%), 150 mg q12h in 16 (3.9%), 100 mg q12h in 9 (2.2%), 200 mg every morning plus 100 mg every night in 7 (1.7%), other therapeutic-drug-monitoring-adjusted regimens in 34 (8.3%). Route on the day of pharmacokinetic sampling: intravenous infusion 68.1%, nasogastric 20.8%, oral 11.0%.",
    regions        = "Single center: China-Japan Friendship Hospital, Beijing, China.",
    renal_function = "Creatinine clearance (Cockcroft-Gault) median 68.5 mL/min, IQR 45.5-102.5; serum creatinine median 78.5 umol/L, IQR 54.4-126.0.",
    notes          = "Retrospective single-center study, 2017-2023. Concentrations measured by a validated UPLC-MS/MS assay, LLOQ 0.097 mg/L, calibration range 0.097-12.500 mg/L. NONMEM 7.2.0, first-order conditional estimation. Baseline demographics per Wang 2024 Table 1; final parameter estimates and 1000-sample nonparametric bootstrap per Table 2. Note that Table 1 reports statistics on a base of 501 on-machine occasions rather than 408 patients for everything except age, sex, weight, height, BMI, dosing method and dosage, so the two denominators are mixed within one table. Prothrombin time, although retained as a covariate on clearance, is not tabulated anywhere in the paper."
  )

  ini({
    # Structural parameters (Wang 2024 Table 2). These are ABSOLUTE
    # population values - the model carries no allometric or other
    # body-size term, weight having been screened and not retained.
    lka     <- fixed(log(1.20)); label("Absorption rate constant (1/h)")             # Wang 2024 Table 2, "Ka (/h) 1.20 (fixed)"; Section 4.4.1 states it "was fixed to a value of 1.2/h, as reported elsewhere [64,65]"
    lcl     <- log(3.55);        label("Clearance (L/h)")                            # Wang 2024 Table 2 (CL 3.55, RSE 3.5%, bootstrap median 3.55, 95% CI 3.33-3.77)
    lvc     <- log(33.5);        label("Central volume of distribution (L)")          # Wang 2024 Table 2 (Vc 33.50, RSE 19.1%, bootstrap median 33.11, 95% CI 22.70-43.38)
    lvp     <- log(138);         label("Peripheral volume of distribution (L)")       # Wang 2024 Table 2 (Vp 138.00, RSE 18.6%, bootstrap median 142.45, 95% CI 107.88-183.39)
    lq      <- log(52.8);        label("Intercompartmental clearance (L/h)")          # Wang 2024 Table 2 (Q 52.80, RSE 15.9%, bootstrap median 53.05, 95% CI 41.34-70.24)
    lfdepot <- log(0.835);       label("Bioavailability of the extravascular dose (unitless fraction)") # Wang 2024 Table 2 (F 0.835, RSE 5.8%, bootstrap median 0.83, 95% CI 0.75-0.93)

    # Covariate effects on clearance. All five come from the single
    # structural equation printed in Wang 2024 Section 2.2:
    #
    #   CL = CL_TV * (qCRP/73.6)^-0.142 * (CLCR/71.8)^0.218 * 1.46^CRRT
    #             * (PLT/144)^0.166 * (PT/15)^-0.875 * exp(eta_CL)
    #
    # That equation is the ONLY source of the four exponent SIGNS: Table 2
    # prints the same five coefficients as unsigned magnitudes (0.142,
    # 1.46, 0.166, 0.875, 0.218). Every sign cross-checks against the
    # Discussion, which reads clearance as rising with CRRT, CLCR and PLT
    # and falling as qCRP and PT rise.
    e_crp_cl    <- -0.142;      label("Exponent of (qCRP / 73.6 mg/L) on clearance (unitless)")               # Wang 2024 Section 2.2 equation; magnitude also in Table 2 (theta qCRP_CL 0.142, RSE 14.6%, bootstrap 95% CI 0.10-0.19)
    e_crcl_cl   <-  0.218;      label("Exponent of (CLCR / 71.8 mL/min) on clearance (unitless)")             # Wang 2024 Section 2.2 equation; magnitude also in Table 2 (theta CLCR_CL 0.218, RSE 15.6%, bootstrap 95% CI 0.14-0.30)
    e_plt_cl    <-  0.166;      label("Exponent of (PLT / 144 x 10^9/L) on clearance (unitless)")             # Wang 2024 Section 2.2 equation; magnitude also in Table 2 (theta PLT_CL 0.166, RSE 25.2%, bootstrap 95% CI 0.10-0.24)
    e_ptsec_cl  <- -0.875;      label("Exponent of (PT / 15 s) on clearance (unitless)")                      # Wang 2024 Section 2.2 equation; magnitude also in Table 2 (theta PT_CL 0.875, RSE 23.2%, bootstrap 95% CI 0.48-1.44)
    # CRRT is a SUPERSCRIPT in the printed equation (1.46^CRRT, confirmed
    # in a 400 dpi render of the published page), so it is a 1.46-fold
    # multiplier while CRRT is running and no effect otherwise. Stored on
    # the log scale so it enters clearance as a mu-referenced additive
    # shift, the form already used for categorical clearance effects in
    # Desai_2016_isavuconazole.R. A product reading (1.46 * CRRT) is
    # arithmetically impossible: it would send clearance to zero in the
    # 379 of 501 occasions that were off CRRT.
    e_crrt_cl   <- log(1.46);   label("Log fold-change in clearance while CRRT is running (unitless)")        # Wang 2024 Section 2.2 equation and Table 2 (theta CRRT_CL 1.46, RSE 5.9%, bootstrap median 1.46, 95% CI 1.29-1.65)

    # IIV. Wang 2024 Section 4.4.2 specifies exponential interindividual
    # variability, Pij = Ppop * exp(eta_ij) with eta of variance omega^2.
    # Table 2 reports the random effects as "% CV", which is read here as
    # sqrt(omega^2) * 100 rather than the exact log-normal
    # sqrt(exp(omega^2) - 1) * 100. The arbitration is internal to Table 2:
    # its residual block is headed "(%CV if proportional, SD if additive)"
    # and prints the proportional residual as 8.9, and for a proportional
    # error the only quantity that can mean is sigma * 100. The same
    # authors, in the same table, therefore convert a variance component to
    # a percentage by taking the square root. See the vignette Errata - the
    # alternative reading lowers omega^2 for Vp by about 30%.
    etalcl ~ 0.498^2 # Wang 2024 Table 2, IIV CL 49.80 %CV (RSE 4.4%, bootstrap median 49.29, 95% CI 45.05-53.31)
    etalvc ~ 0.667^2 # Wang 2024 Table 2, IIV Vc 66.70 %CV (RSE 26.6%, bootstrap median 66.65, 95% CI 47.20-90.25)
    etalvp ~ 0.817^2 # Wang 2024 Table 2, IIV Vp 81.70 %CV (RSE 21.8%, bootstrap median 78.39, 95% CI 46.57-115.06)
    # No eta on Q: Wang 2024 Table 2 reports its interindividual
    # variability as "0 (fixed)". A zero-variance eta is mechanically
    # identical to no eta and would make OMEGA singular, so it is omitted
    # rather than written as ~ fixed(0). No eta on F or Ka either - the
    # paper reports none.

    # Residual error. Wang 2024 Section 4.4.2 gives the combined model
    # explicitly as Cobs = Cpred * (1 + eps) + eps'.
    addSd  <- 0.192; label("Additive residual error (mg/L)")          # Wang 2024 Table 2, Additive 0.192 mg/L (RSE 28.1%, bootstrap median 0.20, 95% CI 0.04-0.33)
    propSd <- 0.089; label("Proportional residual error (fraction)")  # Wang 2024 Table 2, Proportional 8.9 under the header "(%CV if proportional, SD if additive)", i.e. sigma = 0.089 (RSE 9.5%, bootstrap median 8.8, 95% CI 3.78-12.19)
  })

  model({
    ka <- exp(lka)

    # Clearance: the Section 2.2 equation. The CRRT term is carried inside
    # exp() as a log-additive shift so the whole typical value stays
    # mu-referenced; the four continuous covariates keep the paper's
    # power form.
    cl <- exp(lcl + e_crrt_cl * RRT_CRRT_ACTIVE + etalcl) *
      (CRP / 73.6)^e_crp_cl *
      (CRCL / 71.8)^e_crcl_cl *
      (PLT / 144)^e_plt_cl *
      (PT_SEC / 15)^e_ptsec_cl

    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q  <- exp(lq)

    fdepot <- exp(lfdepot)

    # Flows and volumes are the stored parameterisation; the transfer
    # micro-constants are derived from them here and nowhere else. Writing
    # k12 / k21 as the primary parameters makes rxSolve()'s default
    # useLinCmt rewrite drop peripheral1 and silently solve a
    # one-compartment model; the vignette gates the simulated terminal
    # half-life against the analytic beta root to prove it did not.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Bioavailability applies to the extravascular (oral / nasogastric)
    # route only; intravenous doses go straight into central, which is
    # how F is identifiable at all in a cohort that received both.
    f(depot) <- fdepot

    # Dose in mg over volume in L gives mg/L, the unit of the reported
    # plasma concentrations (equivalently ug/mL).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
