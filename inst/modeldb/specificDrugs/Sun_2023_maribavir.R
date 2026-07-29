Sun_2023_maribavir <- function() {
  description <- "Two-compartment population PK model for oral maribavir with first-order absorption, an absorption lag time, and dose-dependent absorption rate in adult transplant recipients with cytomegalovirus infection/disease (Sun 2023)"
  reference <- "Sun K, Hayes S, Farrell C, Song IH. Population pharmacokinetic modeling and simulation of maribavir to support dose selection and regulatory approval in adolescents with posttransplant refractory cytomegalovirus. CPT Pharmacometrics Syst Pharmacol. 2023;12(5):719-723. doi:10.1002/psp4.12943"
  vignette <- "Sun_2023_maribavir"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject (baseline weight). Allometric power scaling on CL/F, Vc/F, Q/F and Vp/F with a 70 kg reference; all four exponents were FIXED (0.75 for the clearance terms, 1 for the volume terms) rather than estimated, because the model with estimated exponents returned imprecise estimates for CL/F~weight and Vp/F~weight and was judged insufficiently robust for the adolescent simulations.",
      source_name        = "WTBL"
    ),
    DIS_CMV = list(
      description        = "Transplant recipient with cytomegalovirus infection/disease indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteer / non-CMV participant)",
      notes              = "1 = hematopoietic cell transplant (HCT) or solid organ transplant (SOT) recipient with CMV infection/disease (the phase II SHP620-202 / SHP620-203 and phase III SHP620-303 populations); 0 = healthy volunteer or phase I participant without CMV (including the renal- and hepatic-impairment cohorts). Enters as a log-scale additive shift on CL/F: CL/F is 0.756x lower in transplant recipients with CMV than in the healthy reference. Table S2 states the structural reference population explicitly as 'a 70-kg individual without CMV administered a 800 mg maribavir dose', so DIS_CMV = 0 is the reference.",
      source_name        = "HSCMV"
    ),
    DOSE = list(
      description        = "Administered maribavir dose per administration",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Use case (a) of the DOSE canonical: per-subject assigned dose level entering a power-form covariate effect on the first-order absorption rate, Ka = 0.336 * (DOSE / 800)^-1.94, normalised at 800 mg. The exponent is negative, so Ka decreases as the maribavir dose increases. Observed dose levels in the pooled analysis were single doses of 100, 200 and 400 mg and twice-daily regimens of 400, 800 and 1200 mg.",
      source_name        = "DOSE"
    ),
    CONMED_CYP3A4_INH = list(
      description        = "Concomitant CYP3A4 inhibitor coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no CYP3A4 inhibitor coadministration)",
      notes              = "Time-varying. Multiplicative power-form effect on CL/F: 0.700^CONMED_CYP3A4_INH (a 30 percent reduction in CL/F when 1), consistent with maribavir being cleared mainly by CYP3A4/CYP1A2 metabolism. The pooled dataset's inhibitor exposure comes principally from the dedicated ketoconazole drug-drug-interaction study 1263-102 (Table S1); the source does not enumerate which inhibitor strengths were pooled into the 1 category. This coefficient appears in the final NONMEM control stream but is not tabulated in Table S2 (which presents the model for the reference population, where the indicator is 0), so no RSE or confidence interval is reported for it.",
      source_name        = "CYP3AINH"
    ),
    CONMED_CYP3A4_IND = list(
      description        = "Concomitant CYP3A4 inducer coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no CYP3A4 inducer coadministration)",
      notes              = "Time-varying. Multiplicative power-form effect on CL/F: 2.242^CONMED_CYP3A4_IND (a 2.24-fold increase in CL/F when 1). The pooled dataset's inducer exposure comes principally from the dedicated rifampin drug-drug-interaction study 1263-110 (Table S1); the source does not enumerate which inducer strengths were pooled into the 1 category. This coefficient appears in the final NONMEM control stream but is not tabulated in Table S2 (which presents the model for the reference population, where the indicator is 0), so no RSE or confidence interval is reported for it.",
      source_name        = "CYP3AIND"
    )
  )

  # Covariates that the source screened and carried in the final NONMEM control
  # stream but whose coefficients were FIXED to 0 (i.e. tested and not retained).
  # They are documented here for provenance and are deliberately NOT referenced in
  # model(), because a coefficient fixed at exactly 0 contributes nothing.
  covariatesDataExcluded <- list(
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Screened as an effect on both CL/F (THETA(14)) and Vc/F (THETA(15)) in the final control stream, but both coefficients are '(0 FIX)' and the covariate is additionally hardcoded to SEXN = 0 in the published simulation stream. Not retained in the final model and not reported in Table S2.",
      source_name        = "SEXN"
    ),
    HEPIMP_MOD = list(
      description        = "Moderate hepatic impairment (Child-Pugh class B) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function)",
      notes              = "Screened as an effect on Vc/F (THETA(13), labelled '[Vc~Child-Pugh Class B]') in the final control stream, but the coefficient is '(0 FIX)' and the covariate is additionally hardcoded to HEPN2 = 0 in the published simulation stream. The hepatic-impairment cohort came from phase I study 1263-103 (10 participants with normal hepatic function and 10 with moderate impairment; Table S1). Not retained in the final model and not reported in Table S2.",
      source_name        = "HEPN2"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = "Not reported as a pooled total. Table S1 accounts for 182 participants across the 9 phase I studies and 235 across the 2 phase II studies (417 combined); the phase III SHP620-303 participant count is not given in Table S1. Table S3 reports model-estimated steady-state exposures for 253 phase III (400 mg b.i.d.) and 232 phase II (1200 mg b.i.d.) transplant recipients with CMV.",
    n_studies      = 12,
    age_range      = "Adults. Numeric age range not reported in the article or its supplement.",
    weight_range   = "Not reported numerically; the baseline body weight distribution is shown graphically only, in Figure S1. The model's allometric reference weight is 70 kg.",
    sex_female_pct = "Not reported.",
    race_ethnicity = "Not reported.",
    disease_state  = "Pooled analysis of healthy volunteers, phase I participants with renal or hepatic impairment, and hematopoietic cell transplant (HCT) or solid organ transplant (SOT) recipients with cytomegalovirus (CMV) infection refractory to (with or without resistance to) valganciclovir, ganciclovir, cidofovir or foscarnet.",
    dose_range     = "100-1200 mg orally: single doses of 100, 200 and 400 mg (phase I) and twice-daily regimens of 400, 800 and 1200 mg (phase II/III).",
    regions        = "Not reported.",
    notes          = "Study-by-study designs and per-study participant counts are in Table S1 of the supplement. The final parameter estimates are in Table S2; steady-state exposure summaries used as validation targets are in Table S3. The model was developed in NONMEM 7.4.3 / 7.5.0 (ADVAN4 TRANS4). This model file encodes the 'fixed weight effect exponents' model (control-stream run 171), which is the model the authors used for all reported simulations and exposure estimates; the alternative model with estimated weight exponents (run 174) was explicitly rejected as insufficiently robust and is not extracted."
  )

  ini({
    # Structural parameters. NONMEM parameterises these as MU_n = THETA(n) with
    # PARAM = EXP(MU_n + ETA(n)), so each THETA is already on the natural-log
    # scale and is carried here verbatim at full reported precision. The
    # back-transformed value in the trailing comment is the Table S2 estimate.
    # Reference subject: 70 kg, without CMV, 800 mg maribavir dose (Table S2 footnote).
    lcl   <-  1.32774;  label("Apparent clearance in the reference subject (CL/F, L/h)")                       # supplement s005.txt run 171 THETA(1) [ln_CL]; exp = 3.77 L/h (Table S2, 95% CI 3.50-4.06)
    lvc   <-  2.92052;  label("Apparent central volume of distribution in the reference subject (Vc/F, L)")    # supplement s005.txt run 171 THETA(2) [ln_V2]; exp = 18.6 L (Table S2, 95% CI 17.3-19.8)
    lq    <- -0.096712; label("Apparent intercompartmental clearance in the reference subject (Q/F, L/h)")     # supplement s005.txt run 171 THETA(3) [ln_Q]; exp = 0.908 L/h (Table S2, 95% CI 0.705-1.17)
    lvp   <-  2.15842;  label("Apparent peripheral volume of distribution in the reference subject (Vp/F, L)") # supplement s005.txt run 171 THETA(4) [ln_V3]; exp = 8.66 L (Table S2, 95% CI 7.05-10.6)
    lka   <- -1.092;    label("First-order absorption rate at the 800 mg reference dose (Ka, 1/h)")            # supplement s005.txt run 171 THETA(5) [ln_Ka]; exp = 0.336 1/h (Table S2, 95% CI 0.271-0.415)
    ltlag <- -1.30567;  label("Absorption lag time (h)")                                                      # supplement s005.txt run 171 THETA(6) [ln_ALAG1]; exp = 0.271 h (Table S2, 95% CI 0.241-0.304)

    # Allometric body-weight exponents, all FIXED rather than estimated: the model
    # with estimated exponents gave imprecise CL/F~weight and Vp/F~weight estimates
    # and was rejected for the adolescent simulations.
    e_wt_cl <- fixed(0.75); label("Allometric (WT) exponent on CL/F (unitless)") # supplement s005.txt run 171 THETA(7) '(0.75 FIX)'; Table S2 'CL/F~weight 0.75 fixed'
    e_wt_vc <- fixed(1);    label("Allometric (WT) exponent on Vc/F (unitless)") # supplement s005.txt run 171 THETA(8) '(1 FIX)'; Table S2 'Vc/F~weight 1 fixed'
    e_wt_q  <- fixed(0.75); label("Allometric (WT) exponent on Q/F (unitless)")  # supplement s005.txt run 171 THETA(9) '(0.75 FIX)'; Table S2 'Q/F~weight 0.75 fixed'
    e_wt_vp <- fixed(1);    label("Allometric (WT) exponent on Vp/F (unitless)") # supplement s005.txt run 171 THETA(10) '(1 FIX)'; Table S2 'Vp/F~weight 1 fixed'

    # Covariate effects
    e_dis_cmv_cl <- -0.280346; label("Log-effect of transplant-recipient-with-CMV status on CL/F (unitless)")  # supplement s005.txt run 171 THETA(17) [CL~CMV]; exp = 0.756 (Table S2 'CL/F~transplant patients with CMV' 0.756, 95% CI 0.690-0.827)
    e_dose_ka    <- -1.9439;   label("Power exponent of maribavir dose on Ka, normalised at 800 mg (unitless)") # supplement s005.txt run 171 THETA(16) [KA~dose]; Table S2 'Ka~dose' -1.94 (95% CI -2.19 to -1.70)

    # CYP3A4 perpetrator co-medication effects on CL/F. Stored as logs of the
    # multiplicative factors so they enter the same log-scale sum as lcl; the
    # source writes them in power form as THETA(11)^CYP3AINH * THETA(12)^CYP3AIND.
    # Not tabulated in Table S2 (see covariateData notes), so no RSE/CI available.
    e_cyp3a4_inh_cl <- log(0.700292); label("Log-effect of concomitant CYP3A4 inhibitor on CL/F (unitless)") # supplement s005.txt run 171 THETA(11) [CL~CYP3AINH] = 0.700292
    e_cyp3a4_ind_cl <- log(2.24181);  label("Log-effect of concomitant CYP3A4 inducer on CL/F (unitless)")   # supplement s005.txt run 171 THETA(12) [CL~CYP3AIND] = 2.24181

    # Inter-individual variability: NONMEM $OMEGA BLOCK(6) over
    # ETA(1..6) = CL, Vc, Q, Vp, Ka, ALAG1 (supplement s005.txt run 171).
    # The 21 values below are the lower triangle row-by-row, exactly as the
    # control stream lists them; the layout below puts one NONMEM row per
    # source line, so the trailing value on each line is that row's diagonal
    # variance and the leading values are its covariances:
    #   line 1 -> var(CL)                                        = 0.243536
    #   line 2 -> cov(CL,Vc),  var(Vc)                           = 0.11584
    #   line 3 -> cov(CL,Q),  cov(Vc,Q),  var(Q)                 = 0.600237
    #   line 4 -> cov(.,Vp) x3,           var(Vp)                = 0.722801
    #   line 5 -> cov(.,Ka) x4,           var(Ka)                = 1.19927
    #   line 6 -> cov(.,lag) x5,          var(lag time)          = 0.177434
    # Those six diagonals reproduce every Table S2 'IIV CV%' entry under that
    # table's footnote-c rule (CV = sqrt(omega^2)*100 when omega^2 <= 0.15,
    # else CV = sqrt(exp(omega^2)-1)*100): CL 52.5%, Vc 34.0%, Q 90.7%,
    # Vp 103%, Ka 152%, lag time 44.1%. The matrix is positive definite as
    # published, so no off-diagonal nudging was required.
    etalcl + etalvc + etalq + etalvp + etalka + etaltlag ~ c(
      0.243536,
      0.133832, 0.11584,
      -0.134868, -0.051648, 0.600237,
      0.0432159, 0.0911984, 0.429564, 0.722801,
      0.098576, 0.133826, -0.439322, -0.314392, 1.19927,
      0.00377063, -0.0215698, 0.150093, 0.114468, -0.32376, 0.177434
    )

    # Residual error. The source $ERROR block applies a study-specific
    # proportional term plus a shared additive term:
    #   Y = IPRED + EPS(1)*IPRED + EPS(2)                       (phase I studies)
    #   Y = IPRED + EPS(3)*IPRED + EPS(2)                       (studies 202/203/303)
    # $SIGMA is 0.0672821 [EPS(1)], '0 FIX' [EPS(2)], 0.136898 [EPS(3)].
    # This model file describes the transplant-recipient-with-CMV population
    # (phase II studies 202/203 and phase III study 303), so it carries EPS(3);
    # the published simulation control stream likewise hardcodes STUDY = 202 and
    # therefore selects EPS(3). The additive term EPS(2) is fixed at exactly zero
    # and so is omitted rather than written as add(fixed(0)) -- an additive SD of 0
    # contributes nothing. The phase I proportional SD is sqrt(0.0672821) = 0.259.
    propSd <- sqrt(0.136898); label("Proportional residual error, phase II/III transplant recipients (fraction)") # supplement s005.txt run 171 $SIGMA third element (EPS(3)) = 0.136898 -> SD 0.370
  })
  model({
    # Individual parameters. NONMEM writes the allometric terms inside the
    # exponent as LOG(WTBL/70)*THETA(n), which is identical to (WT/70)^THETA(n).
    # CYP3A4 perpetrator effects on CL/F: source form is
    # CL * THETA(11)^CYP3AINH * THETA(12)^CYP3AIND.
    cyp3a4_cl <- exp(e_cyp3a4_inh_cl * CONMED_CYP3A4_INH + e_cyp3a4_ind_cl * CONMED_CYP3A4_IND)

    cl   <- exp(lcl + e_dis_cmv_cl * DIS_CMV + etalcl) * (WT / 70)^e_wt_cl * cyp3a4_cl
    vc   <- exp(lvc + etalvc)                          * (WT / 70)^e_wt_vc
    q    <- exp(lq  + etalq)                           * (WT / 70)^e_wt_q
    vp   <- exp(lvp + etalvp)                          * (WT / 70)^e_wt_vp
    ka   <- exp(lka + etalka) * (DOSE / 800)^e_dose_ka
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # Dose in mg, volume in L -> mg/L = ug/mL, matching the Table S3 units
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
