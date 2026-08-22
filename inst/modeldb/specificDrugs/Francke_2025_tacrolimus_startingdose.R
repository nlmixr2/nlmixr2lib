Francke_2025_tacrolimus_startingdose <- function() {
  description <- "Two-compartment population PK model with first-order absorption and an absorption lag time for twice-daily oral immediate-release tacrolimus in adult living- and deceased-donor kidney transplant recipients (Francke 2025 STARTING-DOSE model). This is the companion reduced-covariate model to modellib('Francke_2025_tacrolimus'), refit on the same 13,427 concentrations from 1180 adults but restricted to covariates that are known BEFORE transplantation and are not expected to change much afterwards. Hematocrit and serum creatinine are therefore dropped, leaving four covariate effects on apparent oral clearance CL/F: CYP3A5 genotype (1.64 multiplier for *1/*3 heterozygotes, 1.9 for *1/*1 homozygotes, both relative to the *3/*3 non-expresser reference), CYP3A4*22 carriage (0.842 multiplier vs CYP3A4*1/*1), age (power exponent -0.332 centred at 57.5 years), and height (power exponent 0.97 centred at 170 cm). These explained 33.5% of the inter-individual variance in CL/F. Apparent peripheral volume V2/F is FIXED at 7670 L. Inter-individual variability is a full 3x3 block on CL/F, V1/F, and Q/F; ka and the lag time carry no IIV. Residual variability is a log-scale (exponential) error model with five separate magnitudes selected by the combination of transplant centre and bioanalytical method. The paper pairs this model with a closed-form starting-dose algorithm targeting an AUC0-12h of 234 ng*h/mL (a pre-dose concentration of 10 ng/mL); see the validation vignette."
  reference <- "Francke MI, Sassen SDT, Lloberas N, Colom H, Elens L, Moudio S, de Vries APJ, Moes DJAR, van Schaik RHN, Hesselink DA, de Winter BCM. A Population Pharmacokinetic Model and Dosing Algorithm to Guide the Tacrolimus Starting and Follow-Up Dose in Living and Deceased Donor Kidney Transplant Recipients. Clin Pharmacokinet. 2025;64:1379-1394. doi:10.1007/s40262-025-01533-0"
  vignette <- "Francke_2025_tacrolimus"
  paper_specific_residual_sds <- c(
    "expSdRotterdamImmuno", "expSdRotterdamLeidenLcms", "expSdBarcelonaImmuno",
    "expSdBarcelonaLcms", "expSdBrusselsImmuno"
  )
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Francke 2025 Supplementary Data S3
  # ($SUBROUTINE ADVAN4 TRANS4; S2 = V2/1000 scales the central compartment
  # from mg/L to ng/mL) and Section 2.1.2 (whole-blood tacrolimus).
  compartmentData <- list(
    depot       = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    CYP3A5_STAR1_HET = list(
      description        = "CYP3A5 *1/*3 heterozygote indicator: 1 if the recipient carries exactly one functional CYP3A5*1 allele, 0 otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 with CYP3A5_STAR1_HOM = 0 (CYP3A5 *3/*3 non-expresser, the model reference)",
      notes              = "Time-fixed germline genotype determined from rs776746 (CYP3A5 6986A>G). Paired with CYP3A5_STAR1_HOM so the three-level genotype is encoded by two binary indicators with *3/*3 as the implicit reference (both indicators = 0). Francke 2025 Supplementary Data S3 encodes the source column CYP3A5 as 1 = *3/*3, 2 = *1/*3, 3 = *1/*1, 4 = other, -999 = unknown; the control stream assigns the SAME multiplier to categories 2 AND 4, so the three Leiden recipients with an 'Other' CYP3A5 genotype (Table 1) are pooled with the *1/*3 heterozygotes and must be coded CYP3A5_STAR1_HET = 1. Unknown genotype (-999) receives a multiplier of 1 and is therefore coded as the *3/*3 reference (both indicators = 0). Cohort distribution (Table 1, n = 1180): *3/*3 884 (74.9%), *3/*1 224 (19.0%), *1/*1 39 (3.3%), Other 3 (0.3%), Unknown 30 (2.5%). Genotype is available before transplantation and so is retained in the starting-dose model.",
      source_name        = "CYP3A5"
    ),
    CYP3A5_STAR1_HOM = list(
      description        = "CYP3A5 *1/*1 homozygote indicator: 1 if the recipient carries two functional CYP3A5*1 alleles, 0 otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 with CYP3A5_STAR1_HET = 0 (CYP3A5 *3/*3 non-expresser, the model reference)",
      notes              = "Time-fixed germline genotype determined from rs776746 (CYP3A5 6986A>G). Paired with CYP3A5_STAR1_HET (see that entry for the source-column coding and the pooling of the 'Other' genotype stratum). Cohort n = 39 (3.3%; Table 1). The starting-dose model estimates this multiplier as 1.9, marginally below the full model's 1.93.",
      source_name        = "CYP3A5"
    ),
    SNP_CYP3A4_RS35599367 = list(
      description        = "CYP3A4*22 (rs35599367) reduced-function allele carrier indicator: 1 if the recipient carries at least one CYP3A4*22 allele, 0 if CYP3A4*1/*1.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A4*1/*1 wild-type, the model reference)",
      notes              = "Time-fixed germline genotype. Francke 2025 Supplementary Data S3 codes the source column CYP3A4 as 0 = *1/*1, 1 = *22 carrier, -999 = unknown; unknown receives a multiplier of 1 and is therefore coded as the wild-type reference (0). Dominant (carrier) genetic model -- the paper does not resolve *22 heterozygotes from *22 homozygotes. Cohort distribution (Table 1, n = 1180): *1 1035 (87.7%), *22 112 (9.5%), Unknown 33 (2.8%). Genotype is available before transplantation and so is retained in the starting-dose model.",
      source_name        = "CYP3A4"
    ),
    AGE = list(
      description        = "Recipient age at transplantation",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Enters CL/F as (AGE/57.5)^(-0.332); the centring constant 57.5 years is the dataset median used in Francke 2025 Supplementary Data S3 (Table 1 reports a per-recipient baseline median of 57.0 years, IQR 45.0-65.0). The starting-dose exponent (-0.332) is slightly steeper than the full model's (-0.309), which is expected because age partly absorbs the clearance signal carried by the dropped hematocrit and creatinine covariates. Missing age is coded to a multiplier of 1 per the control stream's IF(AGE.EQ.-999) branch; unknown in 40 of 1180 recipients (Table 1).",
      source_name        = "AGE"
    ),
    HT = list(
      description        = "Recipient height",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject, in CENTIMETRES. Enters CL/F as (HT/170)^0.97; the centring constant 170 is the dataset median and matches the Table 1 cohort median height of 170 cm (IQR 162-177) exactly, which is what establishes the unit. Francke 2025 Table 2 labels this row 'Height (m 2)', which is a typographical error: the Supplementary Data S3 control stream reads IF(HGT.GE.0) CLHGT = ((HGT/170)**THETA(17)) and the Discussion states 'height (in cm)'. Missing height is coded to a multiplier of 1 per the control stream's IF(HGT.EQ.-999) branch; unknown in 171 of 1180 recipients (Table 1).",
      source_name        = "HGT"
    ),
    STUDY_TACRO_FRANCKE = list(
      description        = "Integer (1-7) per-observation code for the combination of transplant centre and tacrolimus bioanalytical method, used to select the residual-error magnitude.",
      units              = "(integer 1-7)",
      type               = "categorical",
      reference_category = "None -- every observation belongs to exactly one of the seven strata and each stratum selects one of five estimated residual-error magnitudes.",
      notes              = "Per-observation (record-level) indicator. Francke 2025 Section 2.2.1 states that residual variability was modelled for the combination of transplant centre (Erasmus MC Rotterdam, LUMC Leiden, Bellvitge Barcelona, St. Luc Brussels) and measurement method (immunoassay vs LC-MS/MS), giving seven combinations, and that 'residual errors were combined where possible based on analytical method and fit of the model' -- collapsing the seven strata onto five estimated magnitudes. The Supplementary Data S3 $ERROR block gives the mapping exactly: codes 1 and 2 share ERR1 (Rotterdam immunoassays), codes 3 and 4 share ERR2 (Rotterdam and Leiden LC-MS/MS), code 5 is ERR3 (Barcelona immunoassay), code 6 is ERR4 (Barcelona LC-MS/MS), and code 7 is ERR5 (Brussels immunoassay). The source paper does NOT state what separates code 1 from code 2 or code 3 from code 4; because each pair shares a single estimated magnitude, the distinction is immaterial to the model's predictions. For a new single-centre LC-MS/MS dataset, code 3 (or 4) is the closest analogue.",
      source_name        = "DATA (referred to as CENTER in $INPUT)"
    )
  )

  covariatesDataExcluded <- list(
    HCT = list(
      description = "Hematocrit, expressed as a fraction of total blood volume (0-1, L/L)",
      units       = "fraction",
      type        = "continuous",
      notes       = "Retained in the Francke 2025 full model (exponent -0.51 centred at 0.33 L/L) but DELIBERATELY EXCLUDED from the starting-dose model. Section 3.4: 'The covariates that change substantially after transplantation (hematocrit and serum creatinine) were not included in the starting dose algorithm, as the values prior to transplantation cannot be used as predictors for the tacrolimus pharmacokinetics in the post-transplant period.' Use modellib('Francke_2025_tacrolimus') when post-transplant hematocrit is available."
    ),
    CREAT = list(
      description = "Serum creatinine concentration",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Retained in the Francke 2025 full model (exponent -0.0905 centred at 147 umol/L) but DELIBERATELY EXCLUDED from the starting-dose model, for the same reason as hematocrit (Section 3.4). Use modellib('Francke_2025_tacrolimus') when post-transplant serum creatinine is available."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1180L,
    n_studies      = 4L,
    age_range      = "IQR 45.0-65.0 years (all recipients at least 18 years old by protocol; range not reported)",
    age_median     = "57.0 years",
    weight_range   = "IQR 65.0-86.0 kg (range not reported)",
    weight_median  = "75.6 kg",
    height_median  = "170 cm (IQR 162-177)",
    sex_female_pct = 37.5,
    race_ethnicity = "Not reported. The four contributing centres are in the Netherlands (Rotterdam, Leiden), Spain (Barcelona), and Belgium (Brussels).",
    disease_state  = "Adult kidney transplant recipients (living-donor and deceased-donor grafts) receiving twice-daily oral immediate-release tacrolimus (Prograf, Astellas Pharma, or Adport, Sandoz Pharma) with mycophenolic acid and a tapering glucocorticoid course as maintenance immunosuppression. Blood-group-ABO-incompatible and HLA-incompatible transplant recipients were excluded. Donor type: living 536 (45.4%), deceased 342 (29.0%), unknown 302 (25.6%). Median time after transplantation at sampling 31 days (IQR 10-91).",
    dose_range     = "Median 4 mg per administration (IQR 2.5-6.5), twice daily; doses were set by local practice and adjusted by therapeutic drug monitoring. Median whole-blood tacrolimus concentration 9.2 ng/mL (IQR 6.9-12.6); median pre-dose concentration 8.8 ng/mL (IQR 6.7-11.8).",
    regions        = "Netherlands (Erasmus MC, University Medical Center Rotterdam, n = 547; Leiden University Medical Center, n = 100), Spain (Bellvitge University Hospital, Barcelona, n = 444), Belgium (Cliniques Universitaires St Luc, Brussels, n = 89).",
    genotypes      = "CYP3A5: *3/*3 884 (74.9%), *3/*1 224 (19.0%), *1/*1 39 (3.3%), Other 3 (0.3%), Unknown 30 (2.5%). CYP3A4: *1 1035 (87.7%), *22 112 (9.5%), Unknown 33 (2.8%).",
    notes          = "Baseline characteristics from Francke 2025 Table 1 (Total column, n = 1180). Same model-building dataset as the full model; Section 3.4 notes that because time after transplantation did not improve the population PK model, ALL tacrolimus concentrations (not only early post-transplant ones) were used to develop the starting-dose model. The starting-dose covariates explained 33.5% of the inter-individual variance in CL/F versus 35.0% for the full model (Discussion paragraph 6)."
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural parameters -- Francke 2025 Table 2, "Starting dose model"
    # column. Reference subject for the typical values: CYP3A5 *3/*3
    # non-expresser, CYP3A4*1/*1, age 57.5 years, height 170 cm (i.e. every
    # covariate factor equals 1). Apparent clearances (CL/F, Q/F) in L/h,
    # apparent volumes (V1/F, V2/F) in L, ka in 1/h, tlag in h. Values in
    # parentheses after each estimate are the %RSE reported by Table 2.
    # ---------------------------------------------------------------------
    ltlag <- log(0.375); label("Absorption lag time tlag (h)")                                     # Table 2 starting dose model: t_lag = 0.375 h (RSE 5.7%)
    lka   <- log(6.25) ; label("First-order absorption rate constant ka (1/h)")                    # Table 2 starting dose model: Ka = 6.25 (RSE 29.6%). The Table 2 column header reads "K a (L/h)" which is a typo; the parameter is a first-order rate constant with units 1/h.
    lcl   <- log(19.9) ; label("Apparent oral clearance CL/F at the reference covariate set (L/h)") # Table 2 starting dose model: CL/F = 19.9 L/h (RSE 1.6%)
    lvc   <- log(653)  ; label("Apparent central volume V1/F (L)")                                 # Table 2 starting dose model: V1/F = 653 L (RSE 4.9%)
    lq    <- log(9.14) ; label("Apparent inter-compartmental clearance Q/F (L/h)")                 # Table 2 starting dose model: Q/F = 9.14 L/h (RSE 3.8%)

    # V2/F is FIXED in every Table 2 column, at the value estimated from the
    # densely sampled subset (Section 3.2.1); Supplementary Data S3 $THETA
    # reads "(7670 FIX) ; 5 V2 - 7670".
    lvp   <- fixed(log(7670)); label("Apparent peripheral volume V2/F (L)")                 # Table 2 starting dose model: V2/F = 7670 L (f); Supplementary Data S3 $THETA(5) = 7670 FIX

    # ---------------------------------------------------------------------
    # Covariate effects on CL/F -- Francke 2025 Table 2 "Covariate effect on
    # CL" rows, "Starting dose model" column. Hematocrit and serum creatinine
    # are absent by design (Section 3.4). The structure is otherwise identical
    # to the full model: Eq 1 median-normalised power effects for the two
    # continuous covariates and Eq 2 multiplicative factors for the genotypes.
    #
    #   CL/F = 19.9
    #          * [1.0 (CYP3A5 *3/*3) | 1.64 (*1/*3) | 1.9 (*1/*1)]
    #          * [1.0 (CYP3A4*1)     | 0.842 (CYP3A4*22 carrier)]
    #          * (Age/57.5)^-0.332
    #          * (Height/170)^0.97
    # ---------------------------------------------------------------------
    e_cyp3a5_het_cl <- 1.64  ; label("CYP3A5 *1/*3 multiplier on CL/F (vs the *3/*3 non-expresser reference)") # Table 2 starting dose model: CYP3A5*1/*3 = 1.64 (RSE 2.9%)
    e_cyp3a5_hom_cl <- 1.9   ; label("CYP3A5 *1/*1 multiplier on CL/F (vs the *3/*3 non-expresser reference)") # Table 2 starting dose model: CYP3A5*1/*1 = 1.9 (RSE 4.7%)
    e_cyp3a4s22_cl  <- 0.842 ; label("CYP3A4*22 carrier multiplier on CL/F (vs the CYP3A4*1/*1 reference)")    # Table 2 starting dose model: CYP3A4*22 = 0.842 (RSE 3.9%)
    e_age_cl        <- -0.332; label("Age power exponent on CL/F, centred at 57.5 years (unitless)")           # Table 2 starting dose model: Age (years) = -0.332 (RSE 12.1%)
    e_ht_cl         <- 0.97  ; label("Height power exponent on CL/F, centred at 170 cm (unitless)")            # Table 2 starting dose model: Height = 0.97 (RSE 22.9%). The Table 2 row label "(m 2)" is a typo -- the control stream centres at 170 and the Discussion states "height (in cm)".

    # ---------------------------------------------------------------------
    # Inter-individual variability. Same $OMEGA BLOCK(3) structure as the full
    # model: exponential IIV on CL/F, V1/F, and Q/F; none on ka or tlag.
    # Table 2 reports the diagonal as %CV (see the full model file for the
    # falsification of the %CV reading against the paper's own 19.3% and 35%
    # summary statistics) and a separate correlation matrix.
    #
    # Log-scale variances via omega^2 = log(1 + CV^2):
    #   CL/F CV 41.8% -> log(1 + 0.418^2) = 0.161033  (omega 0.401289)
    #   V1/F CV 77.8% -> log(1 + 0.778^2) = 0.473301  (omega 0.687969)
    #   Q/F  CV 81.0% -> log(1 + 0.810^2) = 0.504465  (omega 0.710257)
    # Off-diagonals are cov_ij = r_ij * omega_i * omega_j:
    #   CL x V1  r = 0.584 -> 0.584 * 0.401289 * 0.687969 = 0.161228
    #   CL x Q   r = 0.079 -> 0.079 * 0.401289 * 0.710257 = 0.022516
    #   V1 x Q   r = 0.140 -> 0.140 * 0.687969 * 0.710257 = 0.068409
    # The resulting matrix is positive-definite (eigenvalues 0.5972, 0.4490,
    # 0.0926).
    # ---------------------------------------------------------------------
    # Block entries in lower-triangular order, all from Table 2 "Starting dose
    # model":
    #   [1,1] 0.161033 -- IIV CL/F 41.8% CV
    #   [2,1] 0.161228 -- correlation CL/F x V1 = 58.4%
    #   [2,2] 0.473301 -- IIV V1/F 77.8% CV
    #   [3,1] 0.022516 -- correlation CL/F x Q = 7.9%
    #   [3,2] 0.068409 -- correlation V1 x Q = 14%
    #   [3,3] 0.504465 -- IIV Q/F 81% CV
    etalcl + etalvc + etalq ~ c(
      0.161033,
      0.161228, 0.473301,
      0.022516, 0.068409, 0.504465
    )

    # ---------------------------------------------------------------------
    # Residual variability -- log-scale (exponential) error model, five
    # magnitudes selected by transplant centre x bioanalytical method. See the
    # full model file for the Supplementary Data S3 evidence that these THETAs
    # are standard deviations on the natural-log scale ($SIGMA 1 FIX with
    # Y = LOG(F) + THETA*EPS).
    # ---------------------------------------------------------------------
    expSdRotterdamImmuno     <- 0.245; label("Log-scale residual SD, Rotterdam immunoassays")        # Table 2 starting dose model: Rotterdam-immunoassays = 0.245 (RSE 3.3%)
    expSdRotterdamLeidenLcms <- 0.284; label("Log-scale residual SD, Rotterdam and Leiden LC-MS/MS") # Table 2 starting dose model: Rotterdam and Leiden-LC-MS/MS = 0.284 (RSE 2.3%)
    expSdBarcelonaImmuno     <- 0.381; label("Log-scale residual SD, Barcelona immunoassay")         # Table 2 starting dose model: Barcelona-immunoassay = 0.381 (RSE 4.3%)
    expSdBarcelonaLcms       <- 0.419; label("Log-scale residual SD, Barcelona LC-MS/MS")            # Table 2 starting dose model: Barcelona-LC-MS/MS = 0.419 (RSE 6.2%)
    expSdBrusselsImmuno      <- 0.323; label("Log-scale residual SD, Brussels immunoassay")          # Table 2 starting dose model: Brussels-immunoassay = 0.323 (RSE 4.3%)
  })

  model({
    # 1. Covariate factors on CL/F, each equal to 1 at the reference covariate
    #    set. Hematocrit and serum creatinine are absent by design.
    f_cyp3a5 <- 1 + (e_cyp3a5_het_cl - 1) * CYP3A5_STAR1_HET +
                    (e_cyp3a5_hom_cl - 1) * CYP3A5_STAR1_HOM

    f_cyp3a4 <- 1 + (e_cyp3a4s22_cl - 1) * SNP_CYP3A4_RS35599367

    f_age <- (AGE / 57.5) ^ e_age_cl
    f_ht  <- (HT  / 170)  ^ e_ht_cl

    # 2. Individual parameters. Bioavailability F was FIXED to 1 (Section
    #    2.2.1), so CL, Q, V1, and V2 are apparent (/F) quantities throughout.
    tlag <- exp(ltlag)
    ka   <- exp(lka)
    cl   <- exp(lcl + etalcl) * f_cyp3a5 * f_cyp3a4 * f_age * f_ht
    vc   <- exp(lvc + etalvc)
    q    <- exp(lq  + etalq)
    vp   <- exp(lvp)

    # 3. Micro-constants, matching ADVAN4 TRANS4.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 4. Two-compartment disposition with first-order oral absorption.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. Absorption lag time on the depot (ALAG1 in the control stream).
    alag(depot) <- tlag

    # 6. Observation. Whole-blood tacrolimus in ng/mL; central/vc is mg/L =
    #    ug/mL, so multiply by 1000. Reproduces S2 = V2/1000.
    Cc <- central / vc * 1000

    w_indiv <-
      expSdRotterdamImmuno     * ((STUDY_TACRO_FRANCKE == 1) + (STUDY_TACRO_FRANCKE == 2)) +
      expSdRotterdamLeidenLcms * ((STUDY_TACRO_FRANCKE == 3) + (STUDY_TACRO_FRANCKE == 4)) +
      expSdBarcelonaImmuno     * (STUDY_TACRO_FRANCKE == 5) +
      expSdBarcelonaLcms       * (STUDY_TACRO_FRANCKE == 6) +
      expSdBrusselsImmuno      * (STUDY_TACRO_FRANCKE == 7)

    Cc ~ lnorm(w_indiv)
  })
}
