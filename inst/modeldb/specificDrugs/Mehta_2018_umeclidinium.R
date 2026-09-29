Mehta_2018_umeclidinium <- function() {
  description <- "Two-compartment population PK model with first-order absorption for inhaled umeclidinium in patients with COPD receiving single-inhaler fluticasone furoate/umeclidinium/vilanterol triple therapy, with body weight, age and creatinine clearance on apparent clearance and body weight on apparent central volume, refit on the FULFIL study pooled with the historical umeclidinium/vilanterol program"
  reference <- "Mehta R, Pefani E, Beerahee M, Brealey N, Barnacle H, Birk R, Zhu CQ, Lipson DA. Population Pharmacokinetic Analysis of Fluticasone Furoate/Umeclidinium/Vilanterol via a Single Inhaler in Patients with COPD. J Clin Pharmacol. 2018;58(11):1461-1467. doi:10.1002/jcph.1253"
  vignette <- "Mehta_2018_fluticasoneFuroate_umeclidinium_vilanterol"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")
  # Unit note: doses are entered in ug and volumes are in L, so `Cc` is in
  # ug/L == ng/mL. Mehta 2018 reports umeclidinium concentrations and exposures
  # in pg/mL and pg*h/mL (assay LLQ 10 pg/mL, upper limit 2000 pg/mL, Sect.
  # 'Pharmacokinetic Assessments'); multiply `Cc` by 1000 to compare against
  # the published values. The source control stream (Supplementary Table 2)
  # carries the same factor as `S2 = V2/1000`; it is NOT reproduced inside the
  # model here so that dose / volume / clearance stay mutually consistent.

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters CL/F and V2/F as separate power terms normalised to 70 kg, per Supplementary Table 2 '$PK MU_1=LOG(THETA(1)) + WTEX1*LOG(WT/70) + ...' and 'MU_2=LOG(THETA(2)) + WTEX2*LOG(WT/70)'. Note the Methods narrative lists body weight, age and creatinine clearance on CL/F only; the deposited control stream additionally carries body weight on V2/F, and the stream is the authoritative statement of the fitted model. Mean weight in the FULFIL PK population was 81 kg (Table 3). Neither exponent is re-estimated or reported for the combined dataset, so both are held at a structural zero -- see the ini() comments.",
      source_name = "WT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power-model effect on CL/F normalised to 60 years, per Supplementary Table 2 '$PK ... + AGEX*LOG(AGE/60) ...'. Mean age in the FULFIL PK population was 64 years (Table 3); FULFIL required patients to be at least 40 years old. The exponent is not re-estimated or reported for the combined dataset and is held at a structural zero.",
      source_name = "AGE"
    ),
    CRCL = list(
      description = "Creatinine clearance (raw Cockcroft-Gault-style estimate, NOT BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Power-model effect on CL/F normalised to 110 mL/min, per Supplementary Table 2 '$PK ... + CRCLEX*LOG(CRCL/110)'. The source dataset column is `CRCL`, carried alongside a banded `CRCLCAT`; the control stream does not state a BSA normalisation, and the 110 mL/min reference is in the range expected for a raw creatinine clearance rather than a BSA-normalized eGFR in this population. Mehta 2018 reports no renal-function summary for the PK population. The exponent is not re-estimated or reported for the combined dataset and is held at a structural zero.",
      source_name = "CRCL"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "umeclidinium", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "umeclidinium", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "umeclidinium", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 74,
    n_studies = 1,
    age_mean = "64 years",
    weight_mean = "81 kg",
    sex_female_pct = 26,
    race_ethnicity = c(`White` = 100),
    disease_state = "symptomatic chronic obstructive pulmonary disease (mean 45% predicted FEV1)",
    dose_range = "umeclidinium 62.5 ug once daily by oral inhalation, as the fluticasone furoate/umeclidinium/vilanterol 100/62.5/25 ug single-inhaler triple combination (Ellipta)",
    regions = "global (162 centers in 15 countries: Russian Federation, Ukraine, Mexico, Germany, Greece, Czech Republic, Romania, Bulgaria, China, Estonia, Hungary, Italy, Poland, Republic of Korea, Slovakia)",
    notes = "Demographics are the FULFIL (CTT116853, NCT02345161) PK population of 74 patients randomized to fluticasone furoate/umeclidinium/vilanterol who provided serial (n = 10) or sparse (n = 64) samples at weeks 12 and 24 (Table 3); mean BMI 28 kg/m2, mean height 171 cm. The PARAMETER ESTIMATES in this file were obtained on a COMBINED dataset that pools these FULFIL data with the historical umeclidinium program data (source file 'Final_Anoro_FulFill_UMEC_AllDoses.csv') used to build the Goyal 2014 umeclidinium/vilanterol model; the size of the historical half is not restated in Mehta 2018, so n_subjects records the FULFIL contribution only. Data below the 10 pg/mL quantification limit were treated as censored using the NONMEM M3 full-likelihood approach (Ahn 2008)."
  )

  ini({
    # -----------------------------------------------------------------------
    # Structural parameters. Mehta 2018 Table 2 reports the umeclidinium
    # THETAs only on the natural scale, in its 'Model Parameter Estimates With
    # Combined Dataset (RSE%)' column; the 'Historical Model Parameter
    # Estimates' column of the same table is NOT used here. Supplementary
    # Table 2 parameterises them as MU_1 = LOG(THETA(1)) + <covariates>, i.e.
    # the THETAs are natural-scale, so each is wrapped in log() below.
    #
    # The deposited $THETA block holds round-number INITIAL estimates
    # ((0, 145) for CL, (0, 1000) for V2, (0, 500) for Q, (0, 1000) for V3,
    # (1, 5) for KA) that do not match either column of Table 2 -- they are
    # search starting points, not results, and are not used.
    # -----------------------------------------------------------------------
    lka <- log(40.3); label("Apparent first-order absorption rate constant after oral inhalation (1/h)") # Table 2 umeclidinium 'KA (h-1)' combined-dataset estimate 40.3 (RSE 300%); very imprecisely determined, absorption is effectively instantaneous relative to disposition
    lcl <- log(210); label("Apparent inhaled clearance CL/F at the covariate reference (L/h)") # Table 2 umeclidinium 'CL/F (L/h)' combined-dataset estimate 210 (RSE 2.9%)
    lvc <- log(1170); label("Apparent central volume of distribution V2/F at the covariate reference (L)") # Table 2 umeclidinium 'V2/F (L)' combined-dataset estimate 1170 (RSE 1.12%)
    lq <- log(854); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2 umeclidinium 'Q/F (L/h)' combined-dataset estimate 854 (RSE 5.4%)
    lvp <- log(16200); label("Apparent peripheral volume of distribution V3/F (L)") # Table 2 umeclidinium 'V3 /F (L)' combined-dataset estimate 16200 (RSE 7.28%)

    # -----------------------------------------------------------------------
    # Covariate effects. Supplementary Table 2 gives the exact functional form
    #   MU_1 = LOG(THETA(1)) + WTEX1*LOG(WT/70) + AGEX*LOG(AGE/60)
    #                        + CRCLEX*LOG(CRCL/110)
    #   MU_2 = LOG(THETA(2)) + WTEX2*LOG(WT/70)
    # i.e. plain power models on WT/70, AGE/60 and CRCL/110.
    #
    # The EXPONENTS are not reported for the combined dataset: Table 2 lists
    # only the five structural rows, and Methods states 'a covariate analysis
    # was not planned. The same covariate relationship was assumed'. The
    # stream's $THETA entries for them (0.4, 0.5, 0.4, 0.5) are round-number
    # initial estimates -- the same block whose structural entries (145, 1000,
    # 500, 1000, 5) are demonstrably initials rather than either column of
    # Table 2. They are therefore held at a structural zero rather than
    # transcribed from an initial estimate; the functional form is preserved so
    # a downstream user can supply exponents from the historical analysis.
    # See the vignette 'Assumptions and deviations' section.
    # -----------------------------------------------------------------------
    e_wt_cl <- fixed(0); label("Power exponent on (WT/70) for CL/F (unitless)") # Supplementary Table 2 '$PK WTEX1= THETA(8)' in 'MU_1=LOG(THETA(1)) + WTEX1*LOG(WT/70)'; no combined-dataset estimate reported in Table 2
    e_age_cl <- fixed(0); label("Power exponent on (AGE/60) for CL/F (unitless)") # Supplementary Table 2 '$PK AGEX= THETA(9)' in 'MU_1=... + AGEX*LOG(AGE/60)'; no combined-dataset estimate reported in Table 2
    e_crcl_cl <- fixed(0); label("Power exponent on (CRCL/110) for CL/F (unitless)") # Supplementary Table 2 '$PK CRCLEX=THETA(11)' in 'MU_1=... + CRCLEX*LOG(CRCL/110)'; no combined-dataset estimate reported in Table 2
    e_wt_vc <- fixed(0); label("Power exponent on (WT/70) for V2/F (unitless)") # Supplementary Table 2 '$PK WTEX2= THETA(10)' in 'MU_2=LOG(THETA(2)) + WTEX2*LOG(WT/70)'; no combined-dataset estimate reported in Table 2

    # -----------------------------------------------------------------------
    # Inter-individual variability. Supplementary Table 2 places an exponential
    # ETA on all five structural parameters, and Results confirms individual
    # MAP Bayes estimates were obtained for all 74 patients. The magnitudes are
    # not reported for the combined dataset (Table 2 has no OMEGA rows) and the
    # stream's $OMEGA block (0.5, 0.7, 0.3, 0.4, 0.3) is round-number initials,
    # so each is declared at fixed(0).
    #
    # The stream also carries INTER-OCCASION variability -- three $OMEGA
    # BLOCK(1) SAME pairs on CL and on V2 across OCC = 1, 2, 3, plus an eta on
    # the proportional residual term (ETA(12)). Per nlmixr2lib convention IOV
    # is not encoded in library models; it is recorded here and in the vignette
    # so the omission is auditable.
    # -----------------------------------------------------------------------
    etalka ~ fixed(0) # Supplementary Table 2 'KA = EXP(MU_5+ETA(5))'; combined-dataset variance not reported
    etalcl ~ fixed(0) # Supplementary Table 2 'CLB = EXP(MU_1+ETA(1))'; combined-dataset variance not reported
    etalvc ~ fixed(0) # Supplementary Table 2 'V2B = EXP(MU_2+ETA(2))'; combined-dataset variance not reported
    etalq ~ fixed(0) # Supplementary Table 2 'Q = EXP(MU_3+ETA(3))'; combined-dataset variance not reported
    etalvp ~ fixed(0) # Supplementary Table 2 'V3 = EXP(MU_4+ETA(4))'; combined-dataset variance not reported

    # -----------------------------------------------------------------------
    # Residual error. Supplementary Table 2 builds a COMBINED additive plus
    # proportional standard deviation on the untransformed concentration scale:
    #   SIG  = THETA(6)                      (additive SD, pg/mL)
    #   SIG2 = F*THETA(7)*EXP(ETA(12))       (proportional SD)
    #   SD   = SQRT(SIG*SIG + SIG2*SIG2)
    #   Y    = F + ERR(1)*SD                 with $SIGMA 1 FIXED
    # Neither magnitude is reported for the combined dataset and the stream's
    # initials (15 and 0.2) are round numbers, so both are declared at fixed(0).
    # The eta on the proportional term is not representable in this file and is
    # documented in the vignette instead.
    # -----------------------------------------------------------------------
    addSd <- fixed(0); label("Additive residual error (ng/mL)") # Supplementary Table 2 '$ERROR SIG=THETA(6)'; combined-dataset magnitude not reported
    propSd <- fixed(0); label("Proportional residual error (fraction)") # Supplementary Table 2 '$ERROR SIG2=F*THETA(7)*EXP(ETA(12))'; combined-dataset magnitude not reported
  })

  model({
    # Power covariate models, exactly as written in Supplementary Table 2's $PK
    # block: the log-additive MU_ formulation is algebraically identical to the
    # multiplicative power form used here.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (AGE / 60)^e_age_cl * (CRCL / 110)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    ka <- exp(lka + etalka)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
