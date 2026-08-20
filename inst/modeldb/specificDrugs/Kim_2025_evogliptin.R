Kim_2025_evogliptin <- function() {
  description <- "Population PK/PD model for evogliptin, a CYP3A4-metabolized dipeptidyl peptidase-4 (DPP-4) inhibitor, in adults spanning normal renal function through end-stage renal disease on hemodialysis (Kim 2025). Plasma evogliptin is described by a two-compartment model with first-order oral absorption, parameterized as apparent clearance CL/F, apparent central and peripheral volumes Vc/F and Vp/F, apparent inter-compartmental clearance Q/F, and absorption rate constant Ka. Body weight enters CL/F, Vc/F, Q/F, and Vp/F through fixed allometric exponents (0.75 on the clearances, 1 on the volumes) referenced to 65 kg. The only covariate effects retained after stepwise selection act on relative bioavailability F1, which is anchored to 1 at the healthy-subject median biochemistry (amylase 59.5 IU/L, triglyceride 112.5 mg/dL) and scaled by the power function F1 = (AMYL/59.5)^0.363 * (TRIG/112.5)^0.268. Both markers rise with worsening renal impairment, so the model expresses the paper's central finding that uremia inhibits the first-pass metabolism of a CYP3A4 substrate and thereby raises its bioavailability. A direct-link sigmoid Emax model maps the model-predicted plasma concentration onto percent inhibition of blood DPP-4 activity with no effect-compartment delay, reflecting the absence of hysteresis in the observed concentration-effect data."
  reference   <- "Kim B, Kim JE, Lee S, Oh J, Cho J-Y, Jang I-J, Lee S, Chung J-Y, Yoon S. Population pharmacokinetic and pharmacodynamic model of evogliptin: Severe uremia increases the bioavailability of evogliptin. CPT Pharmacometrics Syst Pharmacol. 2025;14(2):246-256. doi:10.1002/psp4.13263. Final NONMEM control streams for the PK and PD models are reproduced in Supporting Information (PSP4-14-246-s002.docx)."
  vignette    <- "Kim_2025_evogliptin"
  units       <- list(time = "h", dosing = "ug", concentration = "ug/L")

  compartmentData <- list(
    depot       = list(analyte = "evogliptin", units = "ug", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "evogliptin", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "evogliptin", units = "ug", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (single-dose studies). Applied as fixed allometric scaling to CL/F, Vc/F, Q/F, and Vp/F referenced to 65 kg, before any other covariate was tested (Kim 2025 Table 1 footnote d; supplement Table S4 model 5). Cohort mean weights ranged 60.8-69.8 kg across the five renal-function groups (supplement Table S3).",
      source_name        = "WT"
    ),
    AMYL = list(
      description        = "Blood amylase activity",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed; the arithmetic mean of the screening and Day 1 (pre-dose) clinical-laboratory values was used to damp transient variability (Kim 2025 Methods, 'Development of the covariate PK model'). Enters relative bioavailability as (AMYL / 59.5)^0.363, where 59.5 IU/L is the healthy-subject median. Amylase rises as renal function falls because it is renally cleared (cohort means 57.6 and 61.7 IU/L in the two healthy groups versus 124.3 in severe renal impairment and 145.4 in ESRD on hemodialysis; supplement Table S3), so in this model it acts as a surrogate for uremic burden rather than for pancreatic disease.",
      source_name        = "AMYL"
    ),
    TRIG = list(
      description        = "Blood triglyceride concentration",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed; the arithmetic mean of the screening and Day 1 (pre-dose) clinical-laboratory values was used (Kim 2025 Methods, 'Development of the covariate PK model'). Enters relative bioavailability as (TRIG / 112.5)^0.268, where 112.5 mg/dL is the healthy-subject median. Hypertriglyceridemia is the hallmark of uremic dyslipidemia, so triglyceride is interpreted as a uremia-severity surrogate; note it rises monotonically with renal impairment (means 119.3, 124.5, 163.1, 191.4 mg/dL from healthy to severe) but falls again in ESRD patients on hemodialysis (108.4 mg/dL), which is the mechanism by which the model predicts lower bioavailability in dialysed patients than in severe non-dialysed CKD (Kim 2025 Discussion). Units are load-bearing: the reference 112.5 is in mg/dL, not mmol/L.",
      source_name        = "TG"
    )
  )

  # Screened in the paper's stepwise covariate analysis but not retained in
  # the final model, so they are documented rather than implemented. Kim 2025
  # supplement Table S2 lists the full tested set; Table S4 shows that
  # backward elimination of TG from F1 raised the OFV by 13.0 while AST,
  # chloride, and the renal-function stratifiers never entered.
  covariatesDataExcluded <- list(
    AST = list(
      description = "Serum aspartate aminotransferase activity",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Tested on CL/F and on F1 in the stepwise analysis (supplement Table S2); not retained in the final model."
    ),
    CHLORIDE = list(
      description = "Blood chloride concentration",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Tested on F1 as a candidate uremic-toxin surrogate (supplement Table S2); not retained. Supplement Table S3 shows chloride rising with renal impairment (cohort means 102.9, 105.3, 106.3 mmol/L for mild, moderate, severe) but falling in ESRD patients on hemodialysis (95.6 mmol/L)."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (CKD-EPI 2021 equation)",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = "Used only to stratify subjects into CKD groups for cohort assignment, VPC stratification, and goodness-of-fit stratification; never retained as a covariate on any structural parameter. The paper's thesis is precisely that the biochemical consequences of uremia (amylase, triglyceride), rather than GFR itself, carry the bioavailability signal. The CKD-EPI 2021 equation is reproduced in the supplement."
    ),
    RRT_HEMODIAL_STATUS = list(
      description = "Hemodialysis status",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a categorical covariate (Kim 2025 Methods); not retained in the final model. Its effect is captured indirectly through the lower triglyceride levels observed in dialysed patients."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 46L,
    n_studies      = 2L,
    n_observations = "688 plasma evogliptin concentrations and 598 blood DPP-4 activity measurements (Kim 2025 Abstract and Results).",
    age_range      = "Cohort means 45.6-59.1 years (supplement Table S3)",
    weight_range   = "Cohort means 60.8-69.8 kg (supplement Table S3)",
    sex_female_pct = 47,
    race_ethnicity = c(Korean = 100),
    disease_state  = "Five strata by renal function: normal renal function (16 subjects pooled across the two studies), and mild (8), moderate (8), severe (6) renal impairment plus end-stage renal disease on hemodialysis (8 contributing PK data). Groups were assigned by MDRD-eGFR; EPI-GFR was recalculated for the population analysis. Cohort mean MDRD-eGFR was 100.2, 70.8, 50.5, 22.4, and 6.4 mL/min/1.73 m^2 respectively.",
    dose_range     = "Single oral 5 mg evogliptin",
    regions        = "Republic of Korea (Seoul National University Hospital)",
    renal_function = "Normal through end-stage renal disease on hemodialysis; the ESRD cohort was dosed both after (period 1) and before (period 2) a hemodialysis session, and only period 1 data entered the model.",
    notes          = "Pooled from two phase I studies, DA1229_RI_I (NCT02214693; renal-impairment cohorts and matched healthy subjects, sampled to 120 h) and DA1229_ESRD_I (NCT04195919; ESRD-on-hemodialysis cohort and matched healthy subjects, sampled to 48 h). Baseline demographics are tabulated in supplement Table S3. Sex percentage is the pooled female fraction across all six reported cohorts (22 of 47 enrolled = 46.8%); one ESRD subject dropped out before treatment, giving the 46 subjects who contributed data, and Table S3 does not state that subject's sex, so the fraction among the 46 is 45.7% or 47.8%. Assay LLOQ was 0.5 ug/L (DA1229_RI_I) and 0.2 ug/L (DA1229_ESRD_I); below-quantification values were treated as missing and not substituted."
  )

  ini({
    # ==================================================================
    # Structural PK -- Kim 2025 Table 1, final model column. Typical
    # values are quoted at the 65 kg allometric reference weight; the
    # printed equations in Results ('Final covariate model') give the
    # same numbers with the (WT/65) scaling made explicit.
    # ==================================================================
    lcl <- log(22.1);  label("Apparent clearance CL/F at 65 kg (L/h)")                           # Kim 2025 Table 1: CL/F = 22.1 L/h, RSE 7.3%, bootstrap 22.0 (19.2-25.9)
    lvc <- log(829);   label("Apparent central volume of distribution Vc/F at 65 kg (L)")        # Kim 2025 Table 1: Vc/F = 829.0 L, RSE 7.4%, bootstrap 834.7 (717.3-977.9)
    lka <- log(1.26);  label("First-order absorption rate constant Ka (1/h)")                    # Kim 2025 Table 1: Ka = 1.26 /h, RSE 12.2%, bootstrap 1.29 (1.01-1.67)
    lq  <- log(46.3);  label("Apparent inter-compartmental clearance Q/F at 65 kg (L/h)")        # Kim 2025 Table 1: Q/F = 46.3 L/h, RSE 27.0%, bootstrap 45.9 (17.9-67.8)
    lvp <- log(465);   label("Apparent peripheral volume of distribution Vp/F at 65 kg (L)")     # Kim 2025 Table 1: Vp/F = 465.0 L, RSE 11.4%, bootstrap 461.8 (364.0-567.3)

    # Relative bioavailability is an anchor, not an estimate: the paper
    # fixes the typical value to 1 at the healthy-subject median
    # biochemistry so that the amylase and triglyceride terms below read
    # as fold-changes relative to a healthy reference subject.
    lfdepot <- fixed(log(1)); label("Relative bioavailability F1 at the reference biochemistry (unitless)")  # Kim 2025 Table 1: F1 = 1, RSE NA, footnote c ("Typical value of relative bioavailability was set as 1 for a reference condition"); control stream THETA(6) = (1) FIX

    # Allometric exponents, all held at their canonical values and
    # applied before any other covariate was tested.
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F (unitless)")  # Kim 2025 Table 1: WT ~ CL/F = 0.75, RSE NA, footnote d ("Fixed allometric (WT) exponent")
    e_wt_vc <- fixed(1);    label("Allometric exponent of body weight on Vc/F (unitless)")  # Kim 2025 Table 1: WT ~ Vc/F = 1, RSE NA, footnote d
    e_wt_q  <- fixed(0.75); label("Allometric exponent of body weight on Q/F (unitless)")   # Kim 2025 Table 1: WT ~ Q/F = 0.75, RSE NA, footnote d
    e_wt_vp <- fixed(1);    label("Allometric exponent of body weight on Vp/F (unitless)")  # Kim 2025 Table 1: WT ~ Vp/F = 1, RSE NA, footnote d

    # Covariate effects on relative bioavailability. The exponents are
    # taken from the printed equation in Results ('Final covariate
    # model'), which carries one more significant figure than Table 1's
    # rounded 0.36 / 0.27; the control-stream THETA ordering confirms
    # THETA(10) is the amylase exponent and THETA(9) the triglyceride
    # exponent.
    e_amyl_fdepot <- 0.363; label("Power-model exponent of blood amylase on F1 (unitless)")        # Kim 2025 Results equation F1 = (Amyl/59.5)^0.363 * (TG/112.5)^0.268; Table 1 rounds to 0.36, RSE 34.2%, bootstrap 0.35 (0.17-0.64)
    e_trig_fdepot <- 0.268; label("Power-model exponent of blood triglyceride on F1 (unitless)")   # Kim 2025 Results equation F1 = (Amyl/59.5)^0.363 * (TG/112.5)^0.268; Table 1 rounds to 0.27, RSE 31.0%, bootstrap 0.26 (0.09-0.43)

    # ==================================================================
    # Interindividual variability. The paper's IIV model is log-normal,
    # theta_i = theta_TV * exp(eta_i) (Kim 2025 Methods), and Table 1
    # reports each omega as a bare percentage. Those percentages are the
    # log-scale standard deviation omega * 100, NOT a CV%: simulating
    # both readings against the paper's own Monte Carlo results (Table 2
    # scenarios A, B1, B2, C1) shows the omega-as-SD reading reproduces
    # the published Cmax and AUC0-120 medians to 0.8% and their 5th-95th
    # percentile spreads with roughly 30% less error than the CV%
    # reading on every metric. See the vignette 'IIV scale' section.
    # Only CL/F, Vc/F, Ka, and F1 carry IIV; Q/F and Vp/F are estimated
    # without random effects (Kim 2025 Results; control-stream $PK
    # assigns no ETA to Q or V3).
    # ==================================================================
    etalcl     ~ 0.202^2;  # Kim 2025 Table 1: IIV CL/F = 20.2%, RSE 35.4%, bootstrap 20.6 (7.1-30.2)
    etalvc     ~ 0.250^2;  # Kim 2025 Table 1: IIV Vc/F = 25.0%, RSE 36.5%, bootstrap 23.8 (9.2-36.5)
    etalka     ~ 0.768^2;  # Kim 2025 Table 1: IIV Ka = 76.8%, RSE 9.7%, bootstrap 76.7 (62.0-91.5)
    etalfdepot ~ 0.254^2;  # Kim 2025 Table 1: IIV F1 = 25.4%, RSE 32.2%, bootstrap 22.6 (12.0-36.5)

    # Residual error on plasma evogliptin. The control stream writes
    # W = SQRT(THETA(7)^2 * IPRED^2 + THETA(8)^2) with THETA(8) fixed at
    # 1e-5, so the model is proportional in all but name; THETA(7) is a
    # standard deviation (initial estimate 0.115, final 0.111).
    propSd <- 0.111;  label("Proportional residual error on plasma evogliptin (fraction)")  # Kim 2025 Table 1: proportional error = 11.1%, RSE 6.0%, bootstrap 11.1 (9.8-12.4)

    # ==================================================================
    # Direct-link sigmoid Emax PD model for inhibition of blood DPP-4
    # activity (Kim 2025 Methods, 'PD model of evogliptin'). No
    # effect-compartment delay: the paper justifies the direct link by
    # the alignment of maximum inhibition with Tmax and the absence of
    # hysteresis.
    # ==================================================================
    lemax <- log(88.9); label("Maximum DPP-4 activity inhibition Emax (% inhibition)")                       # Kim 2025 Table 1: Emax = 88.9 %inhibition, RSE 0.7%, bootstrap 89.0 (87.8-90.4)
    lec50 <- log(1.08); label("Evogliptin concentration producing 50% of Emax, EC50 (ug/L)")                 # Kim 2025 Table 1: EC50 = 1.08 ug/L, RSE 4.6%, bootstrap 1.08 (0.98-1.18)
    lhill <- log(2.14); label("Hill coefficient gamma of the concentration-effect curve (unitless)")         # Kim 2025 Table 1: gamma = 2.14, RSE 3.7%, bootstrap 2.14 (1.97-2.31)

    etalec50 ~ 0.259^2;  # Kim 2025 Table 1: IIV EC50 = 25.9%, RSE 14.7%, bootstrap 25.1 (18.7-33.4); same omega-as-SD reading as the PK IIVs

    # Additive residual error on percent DPP-4 inhibition. Table 1
    # reports 16.7 for this term. That number is the NONMEM $SIGMA
    # value, which is a VARIANCE (the PD control stream estimates
    # $SIGMA directly with no scaling THETA, unlike the PK model where
    # the reported 11.1% is a THETA standard deviation), so the residual
    # standard deviation is sqrt(16.7) = 4.09 %inhibition. Kim 2025
    # Figure 3 corroborates this: at saturating concentrations the
    # observed inhibition scatters roughly +/- 8 %inhibition about the
    # fitted curve, consistent with an SD near 4 and irreconcilable with
    # an SD of 16.7. See the vignette Errata.
    addSd_DPP4INH <- 4.087; label("Additive residual error on DPP-4 activity inhibition (% inhibition)")  # sqrt(16.7); Kim 2025 Table 1: PD additive error = 16.7, RSE 18.1%, bootstrap 16.2 (11.6-22.9), footnote e ("The unit of PD error is in DPP-4 activity %inhibition")
  })

  model({
    # 1. Individual PK parameters. Body weight scales the clearances and
    #    volumes allometrically about the 65 kg reference; relative
    #    bioavailability carries the two retained covariate effects.
    ka     <- exp(lka + etalka)
    cl     <- exp(lcl + etalcl) * (WT / 65)^e_wt_cl
    vc     <- exp(lvc + etalvc) * (WT / 65)^e_wt_vc
    q      <- exp(lq)           * (WT / 65)^e_wt_q
    vp     <- exp(lvp)          * (WT / 65)^e_wt_vp
    fdepot <- exp(lfdepot + etalfdepot) *
              (AMYL / 59.5)^e_amyl_fdepot *
              (TRIG / 112.5)^e_trig_fdepot

    # 2. Micro-constants (control stream: K12 = Q/V2, K21 = Q/V3, K = CL/V2)
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. Two-compartment disposition with first-order oral absorption
    #    (NONMEM ADVAN4 TRANS4 in the PK run, written out as $DES in the
    #    PD run).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - (k12 + kel) * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 4. Relative bioavailability applies to the absorbed oral dose.
    f(depot) <- fdepot

    # 5. Observations. Doses are in ug and vc in L, so central/vc is ug/L.
    Cc <- central / vc

    # Direct-link sigmoid Emax inhibition of blood DPP-4 activity,
    # driven by the model-predicted plasma concentration with no delay:
    #   DPP-4 inhibition (%) = Emax * Cpred^gamma / (EC50^gamma + Cpred^gamma)
    emax <- exp(lemax)
    ec50 <- exp(lec50 + etalec50)
    hill <- exp(lhill)
    DPP4INH <- emax * Cc^hill / (ec50^hill + Cc^hill)

    Cc      ~ prop(propSd)
    DPP4INH ~ add(addSd_DPP4INH)
  })
}
