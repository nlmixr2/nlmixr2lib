Sun_2025_maribavir <- function() {
  description <- paste0(
    "Updated two-compartment population PK model for oral maribavir in ",
    "healthy volunteers, phase I special populations, and hematopoietic ",
    "cell transplant (HCT) or solid organ transplant (SOT) recipients with ",
    "cytomegalovirus (CMV) infection (Sun 2025, n = 930, 7431 concentration ",
    "records pooled across phase 1, 2 and 3 studies including AURORA and ",
    "SOLSTICE). First-order absorption with an absorption lag time, ",
    "first-order elimination, estimated (not fixed) allometric body-weight ",
    "exponents on CL/F, Vc/F, Q/F and Vp/F, strong CYP3A4 inhibitor and ",
    "inducer effects and a CMV disease-state effect on CL/F, a dose effect ",
    "on Ka, and proton-pump-inhibitor effects on both Ka and relative ",
    "bioavailability. Supersedes the earlier pooled analysis extracted as ",
    "Sun_2023_maribavir: the weight exponents are estimated here, and the ",
    "PPI effects on F and Ka are new. The exposure-response analyses ",
    "reported alongside this PK model are not extracted; see the vignette ",
    "Assumptions and deviations section."
  )
  reference <- paste(
    "Sun K, Jomphe C, Gosselin NH, Pheng L, Durairaj C, Hang Y, Bhattacharya I.",
    "Population Pharmacokinetics and Exposure-Response Relationships of Maribavir",
    "in Transplant Recipients With First Episode or Refractory Cytomegalovirus.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(8):1346-1356.",
    "doi:10.1002/psp4.70054.",
    "Final NONMEM control stream in Supporting Information (file s001).",
    sep = " "
  )
  vignette <- "Sun_2025_maribavir"
  paper_specific_residual_sds <- c("propSdPhase1", "propSdPhase23")

  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE: checked against the source. The
  # NONMEM model is ADVAN4 TRANS4 with S2 = V2, so the scaled central
  # compartment is in mg and central/vc is mg/L = ug/mL, matching the
  # ug/mL and ug*h/mL units of Table 3.
  compartmentData <- list(
    depot       = list(analyte = "maribavir", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "maribavir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "maribavir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric power scaling with a 70 kg reference on all four",
        "disposition parameters. Unlike the earlier Sun 2023 analysis, which",
        "fixed the exponents at 0.75 and 1, this model ESTIMATED them, because",
        "the dataset now contains one individual under 18 years of age and the",
        "model was to support dose selection for an ongoing paediatric phase 3",
        "study. The estimated exponents are smaller than the usual theoretical",
        "values: 0.301 for the clearance terms and 0.536 for the volume terms.",
        "The supplement control stream shows that Q/F reuses the CL/F weight",
        "coefficient (MU_3 = THETA(3) + CLWT) and Vp/F reuses the Vc/F weight",
        "coefficient (MU_4 = THETA(4) + VCWT) -- THETA(9) and THETA(10) are",
        "present but commented 'Not used'. This is why Table 2 reports",
        "identical estimates, %RSE and 95% CI for the CL/F and Q/F weight rows",
        "and again for the Vc/F and Vp/F rows: they are ONE parameter each, not",
        "two that happen to agree. The model file therefore carries two weight",
        "exponents, not four. Overall median 73.0 kg (range 36.1-141), Table 1."
      ),
      source_name        = "WTBL"
    ),
    DOSE = list(
      description        = "Administered maribavir dose per administration",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Use case (a) of the DOSE canonical: the per-administration assigned",
        "dose level entering a power-form covariate effect on the first-order",
        "absorption rate, Ka = 0.697 * (DOSE/800)^-1.02, normalised at 800 mg.",
        "The exponent is negative, so Ka decreases as the dose increases. The",
        "control stream comment notes 'dose time changing for a few subjects',",
        "i.e. the column is per-record rather than strictly per-subject. Doses",
        "in the pooled dataset span single doses of 50-1600 mg (phase 1) and",
        "twice-daily regimens up to 1200 mg; the recommended clinical dose is",
        "400 mg twice daily."
      ),
      source_name        = "DOSE"
    ),
    DIS_CMV = list(
      description        = "Transplant recipient with cytomegalovirus infection/disease indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteer / non-CMV phase 1 participant)",
      notes              = paste(
        "1 = HCT or SOT recipient with CMV infection (the phase 2/3 and AURORA",
        "populations, n = 724); 0 = healthy volunteer or phase 1 participant",
        "without CMV, including the renal- and hepatic-impairment cohorts",
        "(n = 206). Derived in the control stream as HSCMV = 1 when the health",
        "status column HS equals 2, with the comment ';HV reference' marking the",
        "healthy volunteer as the reference level. Enters as a log-scale",
        "additive shift on CL/F: CL/F is 0.678x lower in transplant recipients",
        "with CMV, i.e. clearance is 32% lower, which the authors attribute to",
        "reduced liver and/or kidney function and concurrent medications.",
        "Note that transplant TYPE (HCT vs SOT, and organ within SOT) was tested",
        "and was NOT a significant predictor, so this covariate carries the",
        "whole patient-vs-healthy contrast."
      ),
      source_name        = "HSCMV"
    ),
    CONMED_CYP3A4_INH_STRONG = list(
      description        = "Concomitant strong CYP3A4 inhibitor coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no strong CYP3A4 inhibitor coadministration)",
      notes              = paste(
        "Time-varying per record (Table 1 footnote b: the same individual may",
        "appear as both No and Yes). Multiplicative power-form effect on CL/F:",
        "0.709^CONMED_CYP3A4_INH_STRONG, a 29% reduction in CL/F, which the",
        "authors note is consistent with the 35% reduction seen in the dedicated",
        "ketoconazole DDI study. 129 of 930 individuals (13%) had strong",
        "inhibitor exposure. The STRONG-specific canonical is used rather than",
        "the class-level CONMED_CYP3A4_INH because this analysis screened strong",
        "and moderate inhibitors separately (Table 1 lists both) and retained",
        "only the strong effect; the moderate-inhibitor coefficient THETA(20) is",
        "'(1) FIX' in the control stream, i.e. no effect. See",
        "covariatesDataExcluded$CONMED_CYP3A4_INH_MOD."
      ),
      source_name        = "CYP3AINH"
    ),
    CONMED_CYP3A4_IND = list(
      description        = "Concomitant strong CYP3A4 inducer coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no CYP3A4 inducer coadministration)",
      notes              = paste(
        "Time-varying per record. Multiplicative power-form effect on CL/F:",
        "2.27^CONMED_CYP3A4_IND, a 2.27-fold increase in CL/F, which the authors",
        "note is consistent with the 2.5-fold increase seen in the dedicated",
        "rifampin DDI study. Table 1 labels the covariate 'Strong CY3A4",
        "inducers' (18 of 930 individuals, 2%); the paper states explicitly that",
        "there were insufficient patients receiving moderate and weak inducers to",
        "evaluate their effect, so no strength-stratified inducer canonical is",
        "used and the class-level column carries the strong-inducer effect."
      ),
      source_name        = "CYP3AIND"
    ),
    CONMED_PPI = list(
      description        = "Concomitant proton-pump inhibitor use indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no proton-pump inhibitor coadministration)",
      notes              = paste(
        "Time-varying per record (Table 1 footnote b). New in this analysis",
        "relative to the earlier Sun 2023 model. Carries TWO effects, both",
        "multiplicative on the natural scale and both entered in the control",
        "stream as EXP(THETA) gates: relative bioavailability F = 0.905^PPI",
        "(F1 = PPIF) and absorption rate Ka = ... * 0.457^PPI (KA = ... *",
        "PPIKA). Because a lower F lowers exposure while a lower Ka flattens and",
        "delays the profile, the net effect is -9.5% on AUCss and -22.5% on",
        "Cmax,ss with essentially no change in Cmin,ss (Figure 2b) -- which the",
        "authors judge to be of little clinical significance given maribavir's",
        "flat exposure-response. 529 of 973 records-level classifications (54%)",
        "were PPI-exposed, the most prevalent co-medication in the analysis."
      ),
      source_name        = "PPI"
    ),
    STUDY_MARIBAVIR_PHASE1 = list(
      description        = "Phase 1 study cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (phase 2/3 study: SHP620-202, -203, -302, -303 and AURORA)",
      notes              = paste(
        "1 = the concentration record originates from a phase 1 study (n = 206",
        "individuals, 4231 records); 0 = a phase 2/3 study. Used ONLY to switch",
        "the proportional residual-error magnitude, exactly as in the control",
        "stream $ERROR block, which selects EPS(3) instead of EPS(1) when STUDY",
        "is 202, 203, 302 or 303. The authors attribute the difference to the",
        "different LC-MS assays, with different lower limits of quantification,",
        "used in the phase 1 versus the phase 2/3 studies. This column has no",
        "structural effect: it does not enter any PK parameter."
      ),
      source_name        = "STUDY"
    )
  )

  # Covariates that the source screened and carried in the final NONMEM control
  # stream but whose coefficients are FIXED (i.e. tested and not retained), plus
  # the covariates the paper reports as screened-and-rejected. Documented here
  # for provenance and deliberately NOT referenced in model(), because a
  # coefficient fixed at exactly 0 (or a multiplier fixed at exactly 1)
  # contributes nothing.
  covariatesDataExcluded <- list(
    CONMED_CYP3A4_INH_MOD = list(
      description        = "Concomitant moderate CYP3A4 inhibitor coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no moderate CYP3A4 inhibitor coadministration)",
      notes              = paste(
        "Screened as a multiplicative effect on CL/F via THETA(20), labelled",
        "'[CL~CYPMOD]', but the coefficient is '(1) FIX' -- a multiplier of",
        "exactly 1, i.e. no effect. The control-stream header line ';; 1. Based",
        "on: noCYPINHmCL' records that this run is the one built WITHOUT the",
        "moderate-CYP-inhibitor effect on CL. 103 of 930 individuals (11%) had",
        "moderate-inhibitor exposure (Table 1). Not reported in Table 2."
      ),
      source_name        = "CYP3AIHM"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Screened as an effect on both CL/F (THETA(14)) and Vc/F (THETA(15)) in",
        "the final control stream, but both coefficients are '(0) FIX'. The",
        "Discussion states there was no evidence that sex affected maribavir PK.",
        "Note that Figure 2a nonetheless shows ~24% higher steady-state exposure",
        "in females than males; that is a body-weight-mediated difference",
        "propagated through the allometric terms, not a separate sex effect.",
        "930 individuals were 41% female (Table 1)."
      ),
      source_name        = "SEXN"
    ),
    HEPIMP_MOD = list(
      description        = "Moderate hepatic impairment (Child-Pugh class B) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no moderate hepatic impairment)",
      notes              = paste(
        "Screened as an effect on Vc/F (THETA(13), labelled",
        "'[Vc~Child-Pugh Class B]') and derived in the control stream as",
        "HEPN2 = 1 when HEPN equals 2, but the coefficient is '(0) FIX'. 18 of",
        "930 individuals were Child-Pugh class B (Table S1). The Discussion",
        "lists hepatic impairment among the covariates with no evidence of an",
        "effect on maribavir PK."
      ),
      source_name        = "HEPN2"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 930L,
    n_studies      = "Not stated as a single count. The pooled dataset spans phase 1 (n = 206), phase 2/3 (n = 724, which includes SOLSTICE), and AURORA (n = 238, a subset of the phase 2/3 group), plus a phase 1 study in Japanese-descended and non-Hispanic Caucasian individuals (NCT05319353) newly added in this update.",
    n_observations = 7431L,
    age_range      = "12 to <18 years: 1 (<1%); 18 to <65 years: 761 (82%); 65 to <80 years: 168 (18%) (Table 1). The single individual under 18 is the reason the allometric exponents were estimated rather than fixed.",
    weight_range   = "36.1-141 kg; median 73.0 kg, mean 74.2 kg (SD 17.4) (Table 1). AURORA weights not reported.",
    sex_female_pct = 41,
    race_ethnicity = c(Caucasian = 77, Black = 13, Asian = 7, Other = 3),
    disease_state  = "Pooled analysis of healthy volunteers (157), phase 1 participants with hepatic impairment (10), renal impairment (19) or stable renal transplant (20), and HCT or SOT recipients with CMV infection (724). CMV category: no infection 206, asymptomatic infection 644, symptomatic infection 44, CMV organ disease 36. Transplant type: none 186, SOT 304, HCT 440.",
    dose_range     = "Single doses of 50-1600 mg and multiple doses up to 2400 mg/day across the pooled phase 1-3 dataset; the recommended and most-represented regimen is 400 mg twice daily.",
    regions        = "North America, Europe and Asia Pacific (region proportions reported only for the 238-patient AURORA exposure-response subset: North America 24.8%, Europe 58.0%, Asia Pacific 17.2%, Table S2).",
    co_medication  = "Proton-pump inhibitors 54%, strong CYP3A4 inhibitors 13%, moderate CYP3A4 inhibitors 11%, histamine H2 blockers 10%, antacids 8%, weak CYP3A4 inhibitors 6%, strong CYP3A4 inducers 2% (Table 1 and Table S1).",
    notes          = paste(
      "Below-limit-of-quantification data were handled by method M1 (all 297 BLQ",
      "records, 3.5% of post-dose values, excluded). Parameters were estimated in",
      "NONMEM 7.5.1 with IMPMAP; standard errors and 95% non-parametric",
      "confidence intervals came from bootstrap stratified by study. Structure",
      "was read from the final control stream in Supporting Information file",
      "s001; every parameter VALUE comes from the published Table 2, because the",
      "control stream's $THETA / $OMEGA / $SIGMA blocks are the run's INITIAL",
      "estimates (they are close to but not equal to the final values -- e.g.",
      "THETA(16) [KA~dose] is -1.17 initially against a final -1.02, and -1.17",
      "is in fact the lower bound of the published 95% CI)."
    )
  )

  ini({
    # ---- Structural parameters -------------------------------------------
    # NONMEM parameterises these as MU_n = THETA(n) with PARAM = EXP(MU_n +
    # ETA(n)), so each THETA is on the natural-log scale. The published Table 2
    # reports the BACK-TRANSFORMED typical values, so each log below is written
    # as log(<Table 2 estimate>) to keep the source value visible verbatim.
    # Reference subject: 70 kg, no CMV, no CYP3A4 perpetrator, no PPI, 800 mg.
    lcl   <- log(3.94);  label("Apparent clearance in the reference subject (L/h)")                        # Table 2 'CL/F (L/h)' 3.94 (%RSE 3.3, 95% CI 3.69-4.20)
    lvc   <- log(17.8);  label("Apparent central volume of distribution in the reference subject (L)")     # Table 2 'Vc/F (L)' 17.8 (%RSE 2.81, 95% CI 16.9-18.9)
    lq    <- log(1.24);  label("Apparent intercompartmental clearance in the reference subject (L/h)")      # Table 2 'Q/F (L/h)' 1.24 (%RSE 12.9, 95% CI 0.962-1.60)
    lvp   <- log(7.10);  label("Apparent peripheral volume of distribution in the reference subject (L)")  # Table 2 'Vp/F (L)' 7.10 (%RSE 4.2, 95% CI 6.04-8.35)
    lka   <- log(0.697); label("First-order absorption rate at the 800 mg reference dose (1/h)")             # Table 2 'Ka (1/h)' 0.697 (%RSE 8.21, 95% CI 0.594-0.819)
    ltlag <- log(0.212); label("Absorption lag time (h)")                                                       # Table 2 'Lag (h)' 0.212 (%RSE 5.16, 95% CI 0.192-0.235)

    # Relative bioavailability anchor. Table 2 reports 'F  1' with no %RSE and
    # no CI: F is fixed at 1 and only the PPI effect on it is estimated. The
    # control stream has no THETA for baseline F at all -- F1 is set directly to
    # the PPI gate -- so the anchor is structural, hence fixed().
    lfdepot <- fixed(log(1)); label("Relative bioavailability in the reference subject (fraction)")           # Table 2 'F' = 1, reported without uncertainty; control stream 'F1 = PPIF' with PPIF = 1 when PPI = 0

    # ---- Allometric body-weight exponents (ESTIMATED in this analysis) ----
    # Two parameters, not four: the control stream reuses CLWT for Q/F and VCWT
    # for Vp/F (THETA(9) and THETA(10) are present but commented 'Not used'),
    # which is why Table 2's CL/F and Q/F weight rows -- and its Vc/F and Vp/F
    # weight rows -- carry identical estimates, %RSE and 95% CI.
    e_wt_cl <- 0.301; label("Allometric (WT/70) exponent shared by CL/F and Q/F (unitless)")  # Table 2 'Effect of WT on CL/F' and 'Effect of WT on Q/F', both 0.301 (%RSE 24.1, 95% CI 0.159-0.443); control stream THETA(7) [CL~WT], used by MU_1 and MU_3
    e_wt_vc <- 0.536; label("Allometric (WT/70) exponent shared by Vc/F and Vp/F (unitless)") # Table 2 'Effect of WT on Vc/F' and 'Effect of WT on Vp/F', both 0.536 (%RSE 15.1, 95% CI 0.378-0.694); control stream THETA(8) [V~WT], used by MU_2 and MU_4

    # ---- Covariate effects on CL/F ---------------------------------------
    # The source writes these in power form as a multiplier raised to the
    # 0/1 indicator; stored here as logs so they enter the same log-scale sum
    # as lcl, which is exactly what the control stream does for the CMV term
    # (THETA(17)*HSCMV inside MU_1).
    e_cyp3a4_inh_cl <- log(0.709); label("Log-effect of concomitant strong CYP3A4 inhibitor on CL/F (unitless)") # Table 2 'Effect of CYP3AINH on CL/F' x0.709 (%RSE 2.0, 95% CI 0.681-0.737); control stream THETA(11)
    e_cyp3a4_ind_cl <- log(2.27);  label("Log-effect of concomitant CYP3A4 inducer on CL/F (unitless)")          # Table 2 'Effect of CYP3AIND on CL/F' x2.27 (%RSE 2.6, 95% CI 2.15-2.38); control stream THETA(12)
    e_dis_cmv_cl    <- log(0.678); label("Log-effect of transplant-recipient-with-CMV status on CL/F (unitless)") # Table 2 'Effect of CMV on CL/F' x0.678 (%RSE 4.04, 95% CI 0.626-0.733); control stream THETA(17) [CL~CMV]

    # ---- Covariate effects on absorption ---------------------------------
    e_dose_ka <- -1.02;         label("Power exponent of maribavir dose on Ka, normalised at 800 mg (unitless)") # Table 2 'Effect of dose on Ka' x(DOSE[mg]/800)^-1.02 (%RSE 7.3, 95% CI -1.17 to -0.875)
    e_ppi_ka  <- log(0.457);    label("Log-effect of concomitant proton-pump inhibitor on Ka (unitless)")        # Table 2 'Effect of PPI on Ka' x0.457 (%RSE 11.1, 95% CI 0.367-0.568); control stream THETA(18) [KA~PPI]
    e_ppi_f   <- log(0.905);    label("Log-effect of concomitant proton-pump inhibitor on F (unitless)")         # Table 2 'Effect of PPI on F' x0.905 (%RSE 32.2, 95% CI 0.849-0.964); control stream THETA(19) [F~PPI]

    # ---- Inter-individual variability ------------------------------------
    # STRUCTURE from the supplement control stream, which declares three
    # $OMEGA BLOCK(2) blocks over ETA(1..6) = CL, Vc, Q, Vp, Ka, ALAG1 (the
    # $PK assignment order), i.e. CL is correlated with Vc, Q with Vp, and Ka
    # with the lag time, with no correlation across blocks. The header line
    # ';; 2. Description: 3 OMEGA BLOCKs' confirms this is the final structure.
    #
    # DIAGONALS from the published Table 2 'IIV (%)' column, back-transformed
    # with the convention documented in the companion Sun 2023 maribavir paper
    # (Table S2 footnote c) by the same analysis group: CV = sqrt(omega^2)*100
    # when omega^2 <= 0.15, and CV = sqrt(exp(omega^2)-1)*100 otherwise. Only
    # Vc/F falls in the first branch, and the two branches differ there by less
    # than 2% in SD (0.064009 vs 0.062050), so the choice is immaterial. The
    # inverted values are self-consistent under the rule for all six entries:
    #   CL/F 49.5% -> log(1+0.495^2) = 0.2191556   Vc/F 25.3% -> 0.253^2 = 0.0640090
    #   Q/F   145% -> log(1+1.45^2)  = 1.1322082   Vp/F  110% -> log(1+1.10^2) = 0.7929925
    #   Ka   72.0% -> log(1+0.720^2) = 0.4176571   Lag  48.8% -> log(1+0.488^2) = 0.2136135
    #
    # OFF-DIAGONALS are NOT published: Table 2 has no correlation column and no
    # off-diagonal rows. The correlation COEFFICIENTS are therefore carried from
    # the control stream's own $OMEGA initial estimates (0.672580, 0.820357 and
    # -0.685656 for the three blocks) and rescaled onto the published final
    # variances, rather than invented. Because the initial variances are not the
    # final ones, these three covariances are the least well-sourced numbers in
    # this file -- see the vignette Assumptions and deviations section. Each
    # 2x2 block has |r| < 1 and so is positive definite.
    etalcl + etalvc ~ c(
      0.2191556,
      0.07966004, 0.0640090
    )
    etalq + etalvp ~ c(
      1.1322082,
      0.77732130, 0.7929925
    )
    etalka + etaltlag ~ c(
      0.4176571,
      -0.20480040, 0.2136135
    )

    # ---- Residual unexplained variability --------------------------------
    # The control stream $ERROR block is
    #   Y = IPRED + EPS(1)*IPRED + EPS(2)                          (phase 1)
    #   Y = IPRED + EPS(3)*IPRED + EPS(2)   for STUDY 202/203/302/303
    # i.e. a shared additive term and a study-phase-specific proportional term.
    # Table 2 reports all three as VARIANCES (sigma^2), so each is carried as
    # its square root here. The two proportional SDs are switched inside model()
    # by STUDY_MARIBAVIR_PHASE1, following the Cirincione_2017_exenatide_er
    # precedent for a study-specific residual magnitude.
    addSd         <- sqrt(0.000112); label("Additive residual error, all studies (ug/mL)")            # Table 2 'sigma^2 add' = 0.000112 (%RSE 33.8) -> SD 0.0105830; control stream $SIGMA EPS(2)
    propSdPhase1  <- sqrt(0.0705);   label("Proportional residual error, phase 1 studies (fraction)")  # Table 2 'sigma^2 prop Phase 1' = 0.0705 (%RSE 2.5) -> SD 0.265518; control stream $SIGMA EPS(1)
    propSdPhase23 <- sqrt(0.158);    label("Proportional residual error, phase 2/3 studies (fraction)") # Table 2 'sigma^2 prop Phase 2 & 3' = 0.158 (%RSE 3.2) -> SD 0.397492; control stream $SIGMA EPS(3)
  })

  model({
    # ---- Individual parameters -------------------------------------------
    # NONMEM writes the allometric terms inside the exponent as
    # LOG(WTBL/70)*THETA(n), which is identical to (WT/70)^THETA(n).
    # CL/F (L/h) = 3.94 * (WT/70)^0.301 * 0.709^CYP3A4INH * 2.27^CYP3A4IND
    #              * 0.678^CMV            (Table 2 Note)
    cl <- exp(lcl +
                e_cyp3a4_inh_cl * CONMED_CYP3A4_INH_STRONG +
                e_cyp3a4_ind_cl * CONMED_CYP3A4_IND +
                e_dis_cmv_cl    * DIS_CMV +
                etalcl) * (WT / 70)^e_wt_cl

    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q  <- exp(lq  + etalq)  * (WT / 70)^e_wt_cl
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc

    # Ka (1/h) = 0.697 * (DOSE/800)^-1.02 * 0.457^PPI  (Table 2 Note)
    ka   <- exp(lka + e_ppi_ka * CONMED_PPI + etalka) * (DOSE / 800)^e_dose_ka
    tlag <- exp(ltlag + etaltlag)

    # F = 1 * 0.905^PPI  (Table 2 Note; control stream F1 = PPIF)
    fdepot <- exp(lfdepot + e_ppi_f * CONMED_PPI)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    alag(depot) <- tlag
    f(depot)    <- fdepot

    # Dose in mg with vc in L gives mg/L = ug/mL, the units of Table 3.
    Cc <- central / vc

    # Study-phase-specific proportional term plus the shared additive term.
    propSd <- propSdPhase1 * STUDY_MARIBAVIR_PHASE1 +
      propSdPhase23 * (1 - STUDY_MARIBAVIR_PHASE1)
    Cc ~ prop(propSd) + add(addSd)
  })
}
