Ly_2023_cabozantinib <- function() {
  description <- "Two-compartment population PK model for oral cabozantinib tablet (tyrosine kinase inhibitor) in healthy volunteers and patients with differentiated thyroid cancer, renal cell carcinoma, castration-resistant prostate cancer, or hepatocellular carcinoma (Ly 2023, n=1745 across the Phase 3 COSMIC-311 trial and 6 other studies). Absorption is described by two parallel processes sharing a single first-order rate constant Ka: a fraction F1 of the dose enters a depot feeding a chain of 4 transit compartments (the primary process, producing the observed peak near 3 h), and the remaining (1-F1) enters a second depot with a 19.1 h absorption lag (the delayed process, producing the second absorption phase near 24 h). Two-compartment disposition (central + peripheral1) with first-order elimination from central. Covariates are baseline body weight (power on 70 kg) on CL/F and Vc/F and female sex (fractional change) on CL/F; body weight has minimal impact on exposure but a marked impact on Vc/F. Residual error is proportional with separate magnitudes for healthy volunteers and for pooled cancer patients."
  reference <- paste(
    "Ly NS, Li J, Faggioni R, Roskos LK, Brose MS.",
    "Population pharmacokinetics and exposure-response analysis for the",
    "Phase 3 COSMIC-311 trial of cabozantinib for radioiodine-refractory",
    "differentiated thyroid cancer.",
    "Clin Pharmacokinet. 2023;62(4):629-639.",
    "doi:10.1007/s40262-023-01210-0.",
    sep = " "
  )
  vignette <- "Ly_2023_cabozantinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Both parallel absorption routes receive the dose: `depot` feeds the transit
  # chain (fraction F1) and `depot2` is the delayed route (fraction 1 - F1).
  # Declared explicitly because buildModelDb()'s auto-detection recognises only
  # `depot` and `central`, and `central` is NOT dosed in this model.
  dosing <- c("depot", "depot2")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Power effect on CL/F (exponent 0.144) and Vc/F (exponent 2.03), each centered on 70 kg. The 70 kg reference is stated explicitly in the printed covariate equations (Supplementary Material Eq. 5 and Eq. 6) and in Results ('The reference patient was a 70-kg male with DTC receiving 60 mg cabozantinib QD'). Weight had minimal impact on exposure (< 6% change in AUC0-24,ss, < 14% in Cmax,ss, < 4% in Cmin,ss) but a marked impact on Vc/F.",
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male, the typical-value reference)",
      notes              = "Time-fixed. Multiplicative fractional effect on CL/F: CL/F = theta_CL * (WT/70)^theta_WT,CL * (1 + SEX * theta_SEX) with theta_SEX = -0.214 (Supplementary Material Eq. 5). The supplement defines SEX as 'an indicator variable (1 for females and 0 for males)' and theta_SEX as 'the fractional change in cabozantinib oral clearance for females', so females have 21.4% lower CL/F; Results rounds this to 'approximately 20% lower'.",
      source_name        = "SEX"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-volunteer indicator, used only to select the residual-error magnitude",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (cancer patient; the pooled DTC / RCC / CRPC / HCC group)",
      notes              = "Time-fixed. Ly 2023 Results section 3.1.1: 'the structure of the residual error (RE) component of the PK model was best described by a RE term for healthy volunteers and a separate RE term for patients with various cancer types.' Selects propSd_hv = 0.266 (Table 2 'Healthy subjects' 26.6%) when 1 and propSd_pt = 0.363 (Table 2 'Patients' 36.3%) when 0. It has no effect on the structural PK parameters. Healthy volunteers are the 63 subjects of Study XL184-020; the remaining 1682 subjects are cancer patients.",
      source_name        = "POP"
    )
  )

  # Covariates that Ly 2023 tabulated or evaluated but did NOT retain in the
  # final (full) model. Race and patient population were assessed only from
  # post hoc individual estimates of the final model, with no coefficient
  # reported; age and the four laboratory covariates are tabulated as
  # "covariate information" (Table 1, Supplementary Table S2) and are named in
  # the Discussion's no-dose-adjustment conclusion, but no effect is estimated
  # for any of them. Documentation only -- none of these is referenced in
  # model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline",
      units       = "years",
      type        = "continuous",
      notes       = "Tabulated in Table 1 (COSMIC-311 median 66 y, range 32-85) and Supplementary Table S2 (pooled median 65 y, range 19-90). No age effect is estimated; Discussion concludes 'there is no need for dose adjustment for cabozantinib based on body weight, sex, race, age, ALT, AST, total bilirubin, or creatinine clearance'."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Tabulated in Table 1 (COSMIC-311 median 16 U/L) and Supplementary Table S2 (pooled mean 27.57 U/L). No effect estimated; see the Discussion no-dose-adjustment conclusion."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Tabulated in Table 1 (COSMIC-311 median 20 U/L) and Supplementary Table S2 (pooled mean 33.92 U/L). No effect estimated; see the Discussion no-dose-adjustment conclusion."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Tabulated in Table 1 (COSMIC-311 median 7 umol/L, range 3-21) and Supplementary Table S2 (pooled mean 9.05 umol/L). No effect estimated; see the Discussion no-dose-adjustment conclusion."
    ),
    CRCL = list(
      description = "Creatinine clearance estimated by the Cockcroft-Gault equation",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Reported as an absolute (NOT BSA-normalized) Cockcroft-Gault estimate per Supplementary Table S2 footnote b. Tabulated in Table 1 (COSMIC-311 median 83.4 mL/min, range 26.6-182.1) and Supplementary Table S2 (pooled median 85.95 mL/min). No effect estimated; see the Discussion no-dose-adjustment conclusion."
    ),
    RACE_WHITE = list(
      description = "White race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Race (White vs Black vs Asian) was 'evaluated using individual post hoc estimates from the final (full) model' (Methods 2.3), i.e. by inspection of post hoc CL/F and Vc/F rather than as an estimated covariate effect. White is the comparator group (1287/1745 = 73.8%, Supplementary Table S2). No coefficient is reported, so no effect can be encoded."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Assessed post hoc only (Methods 2.3). Results: 'There was a lower (not clinically significant) CL/F and Vc/F in the Asian population compared to the White population', which the Discussion attributes to lower body weight in the Asian subgroup and to its small size (216/1745 = 12.4%). No coefficient is reported, so no effect can be encoded."
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Assessed post hoc only (Methods 2.3). Only 28/1745 (1.6%) subjects overall and 0/101 in COSMIC-311 were Black (Table 1, Supplementary Table S2). No coefficient is reported."
    ),
    TUMTP_DTC = list(
      description = "Differentiated thyroid cancer indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Patient population (healthy volunteer vs CRPC / HCC / RCC / DTC) was evaluated post hoc from final-model individual estimates (Methods 2.3). Results: 'The CL/F and Vc/F range was overlapping in DTC patients and other cancer types', and Discussion: 'There was no marked difference between cabozantinib PK of patients with DTC and that of patients with other cancer types.' DTC n = 101 (5.8%). No coefficient is reported. Note that the separate DIS_HEALTHY covariate IS used, but only to select the residual-error magnitude, not to scale any structural parameter."
    ),
    TUMTP_RCC = list(
      description = "Renal cell carcinoma indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Assessed post hoc only (Methods 2.3). RCC n = 590 (33.8%) from Studies XL184-308 and CheckMate 9ER. No coefficient is reported."
    ),
    TUMTP_HRPC = list(
      description = "Castration-resistant prostate cancer indicator (paper writes CRPC; canonical column TUMTP_HRPC covers both HRPC and CRPC wordings)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Assessed post hoc only (Methods 2.3). CRPC n = 539 (30.9%) from Studies XL184-306 and XL184-307. No coefficient is reported."
    ),
    TUMTP_HCC = list(
      description = "Hepatocellular carcinoma indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Assessed post hoc only (Methods 2.3). HCC n = 452 (25.9%) from Study XL184-309. No coefficient is reported."
    )
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Ly 2023: the analyte is the parent
  # cabozantinib tablet (no metabolite is modelled) and the measured matrix is
  # plasma (Table 2 reports apparent oral CL/F and Vc/F; Table S3 reports plasma
  # concentrations in ng/mL).
  compartmentData <- list(
    depot       = list(analyte = "cabozantinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "cabozantinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit2    = list(analyte = "cabozantinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit3    = list(analyte = "cabozantinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit4    = list(analyte = "cabozantinib", units = "mg", specimen = "administration site", verified = TRUE),
    depot2      = list(analyte = "cabozantinib", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "cabozantinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cabozantinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 1745L,
    n_observations = 4746L,
    n_studies      = 7L,
    age_range      = "19-90 years",
    age_median     = "65 years (pooled across the seven studies; 66 years in COSMIC-311)",
    weight_range   = "35-190.7 kg",
    weight_median  = "77.55 kg (pooled across the seven studies; 70.0 kg in COSMIC-311)",
    sex_female_pct = 17.4,
    race_ethnicity = "White 73.8%, Asian 12.4%, Black 1.6%, Native American/Alaskan 0.3%, Other 3.2%, Unknown 8.8% (Supplementary Table S2)",
    disease_state  = "Pooled population: healthy volunteers (63, 3.6%); renal cell carcinoma (590, 33.8%); castration-resistant prostate cancer (539, 30.9%); hepatocellular carcinoma (452, 25.9%); radioiodine-refractory differentiated thyroid cancer (101, 5.8%, the COSMIC-311 cabozantinib arm).",
    dose_range     = "Single 20, 40 or 60 mg oral tablet (healthy volunteers, Study XL184-020); 60 mg once daily (Studies XL184-306, -307, -308, -309, -311); 40 mg once daily in combination with nivolumab (CheckMate 9ER). Dose interruptions and reductions (60 -> 40 -> 20 mg QD) were permitted in COSMIC-311 to manage adverse events.",
    regions        = "Multinational (Supplementary Table S1 enumerates the seven studies; no regional breakdown is reported)",
    formulations   = "Cabozantinib tablet (Cabometyx) only. The capsule formulation and the medullary-thyroid-cancer population of the earlier Lacy 2018 analysis were deliberately excluded; see Discussion.",
    studies        = c("XL184-020 (Phase 1, healthy volunteers, single 20 / 40 / 60 mg tablet, serial PK, n=63)",
                       "XL184-306 (Phase 3, CRPC, 60 mg QD, n=41)",
                       "XL184-307 (Phase 3, CRPC, 60 mg QD, n=498)",
                       "XL184-308 (Phase 3, RCC, 60 mg QD, n=282)",
                       "XL184-309 (Phase 3, HCC, 60 mg QD, n=452)",
                       "XL184-311 / COSMIC-311 (Phase 3, radioiodine-refractory DTC, 60 mg QD, n=101)",
                       "CheckMate 9ER (Phase 3, RCC, 40 mg QD + nivolumab, n=308)"),
    notes          = "Baseline demographics from Table 1 (COSMIC-311) and Supplementary Table S2 (all seven studies). 4746 quantifiable PK samples from 1745 subjects, of which 205 samples from 101 COSMIC-311 patients. Fourteen subjects were excluded from the PopPK analysis for missing information. XL184-020 was the only study with morning dosing and serial PK sampling; the other six used evening dosing and sparse sampling. The exposure-response component of Ly 2023 used only log-rank tests of Kaplan-Meier curves across exposure tertiles (Table 3) and fitted no parametric exposure-response model, so no exposure-response parameters are extractable from this paper; see the vignette."
  )

  ini({
    # ---- Structural population parameters (Ly 2023 Table 2, final (full) model) ----
    # Reference covariate set: male, 70 kg (Supplementary Material Eq. 5 and
    # Eq. 6 center weight on 70 kg; Results names "a 70-kg male with DTC
    # receiving 60 mg cabozantinib QD" as the reference patient). Apparent
    # (oral) parameters CL/F, Vc/F, Q/F and V3/F are reported in L/h and L,
    # Ka in 1/h, ALAG4 in h.
    lcl <- log(2.05)  ; label("Apparent oral clearance at the reference covariate set (CL/F, L/h)")                             # Ly 2023 Table 2 CL/F = 2.05 (ASE 0.0323, RSE 1.6%, 95% CI 1.98-2.11)
    lvc <- log(98.8)  ; label("Apparent central volume of distribution at the reference covariate set (Vc/F, L)")               # Ly 2023 Table 2 Vc/F = 98.8 (ASE 7.69, RSE 7.8%, 95% CI 83.8-114)
    lq  <- log(15.5)  ; label("Apparent inter-compartmental clearance (Q/F, L/h)")                                              # Ly 2023 Table 2 Q/F = 15.5 (ASE 0.976, RSE 6.3%, 95% CI 13.6-17.4)
    lvp <- log(178)   ; label("Apparent peripheral volume of distribution (V3/F, L)")                                           # Ly 2023 Table 2 V3/F = 178 (ASE 4.32, RSE 2.4%, 95% CI 170-187)

    # Single absorption rate constant. Ly 2023 reports exactly one absorption
    # rate (Table 2 "Ka"), and it serves BOTH parallel absorption processes:
    # it is the transfer rate constant of every step of the transit chain and
    # the first-order rate of the delayed depot2 process. The base model
    # (Supplementary Material Eq. 3) writes Ka = 3.38 * exp(eta_Ka); the final
    # (full) model value in Table 2 is 3.39.
    lka <- log(3.39)  ; label("First-order absorption rate constant, shared by the transit chain and the delayed depot (Ka, 1/h)")  # Ly 2023 Table 2 Ka = 3.39 (ASE 0.175, RSE 5.2%, 95% CI 3.04-3.73)

    # Absorption lag on the delayed (second) process only. NONMEM name ALAG4
    # -- the lag is on compartment 4, the depot of the delayed process, which
    # is also the compartment carrying F4 (Supplementary Material Eq. 4).
    ltlag <- log(19.1) ; label("Absorption lag time of the delayed depot2 process (ALAG4, h)")                                  # Ly 2023 Table 2 ALAG4 = 19.1 (ASE 0.0404, RSE 0.2%, 95% CI 19.1-19.2)

    # Split of the bioavailable dose between the two parallel absorption
    # processes. F1 = 0.735 is the fraction absorbed via the transit-compartment
    # process ("fraction of dose absorbed from first absorption depot", Table 2
    # abbreviation list); the remainder undergoes delayed absorption, which
    # Supplementary Material Eq. 4 writes as F4 = 1 - 0.736 (the base-model
    # value; the final-model value is 0.735).
    # logit(0.735) = log(0.735 / 0.265) = 1.02014.
    logitfrel <- 1.02014 ; label("Logit of the fraction of the bioavailable dose entering the transit-absorption process (F1 = expit(logitfrel))")  # Ly 2023 Table 2 F1 = 0.735 (ASE 0.017, RSE 2.3%, 95% CI 0.702-0.769); logit(0.735) = 1.02014

    # Overall oral bioavailability anchor. Supplementary Material Eq. 4 states
    # the model assumes "an oral bioavailability of 1", so F is a structural
    # assumption rather than an estimate; all disposition parameters are
    # therefore apparent (X/F) quantities.
    lfdepot <- fixed(log(1)) ; label("Overall oral bioavailability (F)")                                          # Ly 2023 Supplementary Material Eq. 4 note: "assuming an oral bioavailability of 1 in the model"

    # ---- Covariate effects (Ly 2023 Table 2; equations in Supplementary Material) ----
    # Continuous: power exponents on (WT / 70). Categorical: multiplicative
    # fractional change, (1 + SEX * theta_SEX).
    e_wt_cl   <-  0.144 ; label("Body-weight power exponent on CL/F (WT / 70 kg; unitless)")                                    # Ly 2023 Table 2 Weight on CL/F = 0.144 (ASE 0.0614, RSE 42.5%, 95% CI 0.0241-0.265); Supplementary Material Eq. 5
    e_wt_vc   <-  2.03  ; label("Body-weight power exponent on Vc/F (WT / 70 kg; unitless)")                                    # Ly 2023 Table 2 Weight on Vc/F = 2.03 (ASE 0.270, RSE 13.3%, 95% CI 1.50-2.55); Supplementary Material Eq. 6
    e_sexf_cl <- -0.214 ; label("Female (vs male) fractional change on CL/F (unitless)")                                        # Ly 2023 Table 2 Female on CL/F = -0.214 (ASE 0.0268, 95% CI -0.266 to -0.161); Supplementary Material Eq. 5

    # ---- Inter-individual variability ----
    # Ly 2023 Table 2 reports IIV as %CV (43.1, 100 and 39.3 for CL/F, Vc/F and
    # Ka), and Supplementary Material Eq. 1-3 give exponential (log-normal)
    # etas. For a log-normally distributed parameter the coefficient of
    # variation is sqrt(exp(omega^2) - 1), so omega^2 = log(1 + CV^2):
    #   CL/F : log(1 + 0.431^2) = 0.17038
    #   Vc/F : log(1 + 1.000^2) = 0.69315
    #   Ka   : log(1 + 0.393^2) = 0.14362
    # No IIV is reported on F1, ALAG4, Q/F or V3/F, and no omega off-diagonals
    # are reported, so Omega is diagonal over these three etas.
    etalcl ~ 0.17038  # Ly 2023 Table 2 IIV CL/F = 43.1 %CV (95% CI 41.0-45.0), shrinkage 19.6%; omega^2 = log(1 + 0.431^2)
    etalvc ~ 0.69315  # Ly 2023 Table 2 IIV Vc/F = 100 %CV (95% CI 88.7-111), shrinkage 66.1%; omega^2 = log(1 + 1.000^2)
    etalka ~ 0.14362  # Ly 2023 Table 2 IIV Ka = 39.3 %CV (95% CI 31.0-46.1), shrinkage 82.0%; omega^2 = log(1 + 0.393^2)

    # ---- Residual error ----
    # Two proportional residual-error magnitudes, selected by subject type.
    # Ly 2023 Results 3.1.1: "the structure of the residual error (RE)
    # component of the PK model was best described by a RE term for healthy
    # volunteers and a separate RE term for patients with various cancer
    # types." Table 2 reports both under "Residual variability" in units of %.
    propSd_hv <- 0.266 ; label("Proportional residual error SD, healthy volunteers (fraction)")                                 # Ly 2023 Table 2 Residual variability, Healthy subjects = 26.6% (ASE 0.598, RSE 2.3%, 95% CI 25.4-27.7)
    propSd_pt <- 0.363 ; label("Proportional residual error SD, pooled cancer patients (fraction)")                             # Ly 2023 Table 2 Residual variability, Patients = 36.3% (ASE 0.615, RSE 1.7%, 95% CI 35.1-37.5), footnote a "Pooled subjects with various cancer types"
  })

  model({
    # ---- 1. Reference covariate value ----
    # Supplementary Material Eq. 5 and Eq. 6 both center body weight on 70 kg.
    ref_wt <- 70  # kg

    # ---- 2. Individual PK parameters ----
    # Eq. 5: CL/F = theta_CL * (WT/70)^theta_WT,CL * (1 + SEX * theta_SEX)
    # Eq. 6: Vc/F = theta_Vc * (WT/70)^theta_WT,Vc
    cl <- exp(lcl + etalcl) *
          (WT / ref_wt)^e_wt_cl *
          (1 + e_sexf_cl * SEXF)
    vc <- exp(lvc + etalvc) *
          (WT / ref_wt)^e_wt_vc
    ka <- exp(lka + etalka)
    q  <- exp(lq)
    vp <- exp(lvp)

    tlag   <- exp(ltlag)
    frel   <- expit(logitfrel)
    fdepot <- exp(lfdepot)

    # ---- 3. Micro-constants ----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---- 4. ODE system: two parallel absorption processes + 2-compartment disposition ----
    # Primary (transit) process: a fraction F1 of the bioavailable dose enters
    # `depot` and passes through a chain of 4 transit compartments before
    # reaching central. Every step uses the single reported rate constant Ka,
    # so the mean transit time is 5/Ka = 1.47 h and the peak concentration
    # falls near 2.5 h, matching the "observed peak concentration around 3 h"
    # that Ly 2023 Results 3.1.1 and Discussion attribute to this process.
    # The chain length follows the paper's own compartment numbering: Table 2
    # names V3/F the PERIPHERAL volume (compartment 3) and puts the lag on
    # compartment 4, so compartments 2 and 3 are central and peripheral and
    # compartment 4 is the delayed depot -- leaving compartment 1 as a depot
    # distinct from the 4 transit compartments. See vignette Errata for the
    # quantitative scoring of this reading against the alternative.
    #
    # Delayed process: the remaining (1 - F1) enters `depot2`, which is held
    # for ALAG4 = 19.1 h and then absorbed first-order at Ka. This reproduces
    # "the increase in cabozantinib exposure at 24 h post-dose relative to the
    # cabozantinib exposure observed at 14 h post-dose in XL184-020".
    #
    # Dosing events must target BOTH `depot` and `depot2` with the same total
    # amount; the f() multipliers below split the dose into the F1 and (1 - F1)
    # fractions. See the vignette make_cohort() for the event-table pattern.
    d/dt(depot)       <- -ka * depot
    d/dt(transit1)    <-  ka * depot    - ka * transit1
    d/dt(transit2)    <-  ka * transit1 - ka * transit2
    d/dt(transit3)    <-  ka * transit2 - ka * transit3
    d/dt(transit4)    <-  ka * transit3 - ka * transit4
    d/dt(depot2)      <- -ka * depot2
    d/dt(central)     <-  ka * transit4 + ka * depot2 -
                          kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ---- 5. Bioavailability and lag ----
    f(depot)     <- fdepot * frel
    f(depot2)    <- fdepot * (1 - frel)
    alag(depot2) <- tlag

    # ---- 6. Observation and error ----
    # central is in mg and vc in L, giving mg/L; multiply by 1000 to report
    # ng/mL, the unit used throughout Ly 2023 (Table 2, Supplementary Table S3).
    Cc <- central / vc * 1000

    # Healthy volunteers carry propSd_hv (26.6%), cancer patients propSd_pt
    # (36.3%). Same switched-proportional-error pattern as
    # Marathe_2023_belzutifan.R.
    propSdEff <- propSd_hv * DIS_HEALTHY + propSd_pt * (1 - DIS_HEALTHY)
    Cc ~ prop(propSdEff)
  })
}
