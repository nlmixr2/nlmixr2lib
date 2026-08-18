Francke_2025_tacrolimus <- function() {
  description <- "Two-compartment population PK model with first-order absorption and an absorption lag time for twice-daily oral immediate-release tacrolimus (Prograf and Adport) in adult living- and deceased-donor kidney transplant recipients (Francke 2025 full model). Developed on 13,427 whole-blood concentrations from 1180 adults across four European transplant centres (Rotterdam, Leiden, Barcelona, Brussels). Apparent oral clearance CL/F carries six covariate effects, all multiplicative on a median-normalised power or categorical scale: CYP3A5 genotype (1.64 multiplier for *1/*3 heterozygotes, 1.93 for *1/*1 homozygotes, both relative to the *3/*3 non-expresser reference), CYP3A4*22 carriage (0.836 multiplier vs CYP3A4*1/*1), hematocrit (power exponent -0.51 centred at 0.33 L/L), age (power exponent -0.309 centred at 57.5 years), serum creatinine (power exponent -0.0905 centred at 147 umol/L), and height (power exponent 1.17 centred at 170 cm). Together these explained 35% of the inter-individual variance in CL/F (equivalently a 19.3% reduction on the CV scale). Apparent peripheral volume V2/F is FIXED at 7670 L, estimated from the subset of patients with at least one densely sampled concentration-time curve. Inter-individual variability is a full 3x3 block on CL/F, V1/F, and Q/F; ka and the lag time carry no IIV. Residual variability is a log-scale (exponential) error model with five separate magnitudes selected by the combination of transplant centre and bioanalytical method. The companion starting-dose model, which drops the two post-transplant laboratory covariates, is modellib('Francke_2025_tacrolimus_startingdose')."
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
      notes              = "Time-fixed germline genotype determined from rs776746 (CYP3A5 6986A>G). Paired with CYP3A5_STAR1_HOM so the three-level genotype is encoded by two binary indicators with *3/*3 as the implicit reference (both indicators = 0). Francke 2025 Supplementary Data S3 encodes the source column CYP3A5 as 1 = *3/*3, 2 = *1/*3, 3 = *1/*1, 4 = other, -999 = unknown; the control stream assigns the SAME multiplier THETA(12) = 1.64 to categories 2 AND 4, so the three Leiden recipients with an 'Other' CYP3A5 genotype (Table 1) are pooled with the *1/*3 heterozygotes and must be coded CYP3A5_STAR1_HET = 1. Unknown genotype (-999) receives a multiplier of 1 and is therefore coded as the *3/*3 reference (both indicators = 0). Cohort distribution (Table 1, n = 1180): *3/*3 884 (74.9%), *3/*1 224 (19.0%), *1/*1 39 (3.3%), Other 3 (0.3%), Unknown 30 (2.5%). Unlike the Andrews 2017 predecessor model, which pooled all CYP3A5*1 carriers into a single expresser stratum, Francke 2025 resolves heterozygotes from homozygotes separately (Discussion paragraph 3).",
      source_name        = "CYP3A5"
    ),
    CYP3A5_STAR1_HOM = list(
      description        = "CYP3A5 *1/*1 homozygote indicator: 1 if the recipient carries two functional CYP3A5*1 alleles, 0 otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 with CYP3A5_STAR1_HET = 0 (CYP3A5 *3/*3 non-expresser, the model reference)",
      notes              = "Time-fixed germline genotype determined from rs776746 (CYP3A5 6986A>G). Paired with CYP3A5_STAR1_HET (see that entry for the source-column coding and the pooling of the 'Other' genotype stratum). Corresponds to CYP3A5 = 3 in the Francke 2025 Supplementary Data S3 control stream, which assigns THETA(13) = 1.93. Cohort n = 39 (3.3%; Table 1).",
      source_name        = "CYP3A5"
    ),
    SNP_CYP3A4_RS35599367 = list(
      description        = "CYP3A4*22 (rs35599367) reduced-function allele carrier indicator: 1 if the recipient carries at least one CYP3A4*22 allele, 0 if CYP3A4*1/*1.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A4*1/*1 wild-type, the model reference)",
      notes              = "Time-fixed germline genotype. Francke 2025 Supplementary Data S3 codes the source column CYP3A4 as 0 = *1/*1, 1 = *22 carrier, -999 = unknown; unknown receives a multiplier of 1 and is therefore coded as the wild-type reference (0). Dominant (carrier) genetic model -- the paper does not resolve *22 heterozygotes from *22 homozygotes. Cohort distribution (Table 1, n = 1180): *1 1035 (87.7%), *22 112 (9.5%), Unknown 33 (2.8%).",
      source_name        = "CYP3A4"
    ),
    HCT = list(
      description        = "Hematocrit, expressed as a fraction of total blood volume (0-1, L/L)",
      units              = "fraction",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying post-transplant laboratory covariate. Francke 2025 reports hematocrit as a fraction (L/L), NOT as percent -- Table 1 gives a cohort median of 0.332 (IQR 0.300-0.370) and the Supplementary Data S3 control stream centres the effect at 0.33 L/L. The canonical-register HCT entry's units (%) are explicitly overridden here so the centring value 0.33 L/L reproduces the paper's equation directly; this follows the Andrews_2017_tacrolimus.R precedent, which makes the same override for the same reason. To use a dataset that records HCT in percent, multiply the column by 0.01 before passing it to this model. Enters CL/F as (HCT/0.33)^(-0.51): higher hematocrit lowers apparent whole-blood clearance because roughly 70-80% of whole-blood tacrolimus is bound to erythrocytes, so a larger erythrocyte mass reduces the unbound fraction available for clearance (Francke 2025 Discussion paragraph 4). Missing hematocrit is coded to a multiplier of 1 (i.e. treated as the 0.33 L/L median) per the control stream's IF(HCT.EQ.-999) branch. Unknown in 3 of 1180 recipients (Table 1). NOT a covariate in the companion starting-dose model, because pre-transplant hematocrit does not predict post-transplant hematocrit.",
      source_name        = "HCT"
    ),
    AGE = list(
      description        = "Recipient age at transplantation",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Enters CL/F as (AGE/57.5)^(-0.309); the centring constant 57.5 years is the dataset median used in Francke 2025 Supplementary Data S3 (Table 1 reports a per-recipient baseline median of 57.0 years, IQR 45.0-65.0 -- the small difference is because the control-stream constant is the median over concentration records rather than over recipients). Older recipients have lower apparent clearance. Missing age is coded to a multiplier of 1 per the control stream's IF(AGE.EQ.-999) branch; unknown in 40 of 1180 recipients (Table 1). The starting-dose model retains age with a slightly different exponent (-0.332).",
      source_name        = "AGE"
    ),
    CREAT = list(
      description        = "Serum creatinine concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying post-transplant laboratory covariate, reported in umol/L (NOT mg/dL) -- Table 1 gives a cohort median of 140 umol/L (IQR 110-191) and the Supplementary Data S3 control stream centres the effect at 147 umol/L (the median over concentration records). Enters CL/F as (CREAT/147)^(-0.0905). Francke 2025 Discussion paragraph 4 notes that less than 1% of tacrolimus is renally cleared, so serum creatinine is interpreted as a surrogate for general health status and body composition rather than a direct clearance mechanism. Missing creatinine is coded to a multiplier of 1 per the control stream's IF(CREAT.EQ.-999) branch; unknown in 22 of 1180 recipients (Table 1). NOT a covariate in the companion starting-dose model, because pre-transplant serum creatinine does not predict post-transplant serum creatinine.",
      source_name        = "CREAT"
    ),
    HT = list(
      description        = "Recipient height",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject, in CENTIMETRES. Enters CL/F as (HT/170)^1.17; the centring constant 170 is the dataset median and matches the Table 1 cohort median height of 170 cm (IQR 162-177) exactly, which is what establishes the unit. Francke 2025 Table 2 labels this row 'Height (m 2)', which is a typographical error: the Supplementary Data S3 control stream reads IF(HGT.GE.0) CLHGT = ((HGT/170)**THETA(17)) and the Discussion states 'height (in cm)'. Taller recipients have higher apparent clearance; Francke 2025 Discussion paragraph 5 argues height is a better dose predictor than body weight because it tracks ideal body weight and is not inflated by adiposity. Missing height is coded to a multiplier of 1 per the control stream's IF(HGT.EQ.-999) branch; unknown in 171 of 1180 recipients (Table 1).",
      source_name        = "HGT"
    ),
    STUDY_TACRO_FRANCKE = list(
      description        = "Integer (1-7) per-observation code for the combination of transplant centre and tacrolimus bioanalytical method, used to select the residual-error magnitude.",
      units              = "(integer 1-7)",
      type               = "categorical",
      reference_category = "None -- every observation belongs to exactly one of the seven strata and each stratum selects one of five estimated residual-error magnitudes.",
      notes              = "Per-observation (record-level) indicator. Francke 2025 Section 2.2.1 states that residual variability was modelled for the combination of transplant centre (Erasmus MC Rotterdam, LUMC Leiden, Bellvitge Barcelona, St. Luc Brussels) and measurement method (immunoassay vs LC-MS/MS), giving seven combinations, and that 'residual errors were combined where possible based on analytical method and fit of the model' -- collapsing the seven strata onto five estimated magnitudes. The Supplementary Data S3 $ERROR block gives the mapping exactly: codes 1 and 2 share ERR1 (Rotterdam immunoassays), codes 3 and 4 share ERR2 (Rotterdam and Leiden LC-MS/MS), code 5 is ERR3 (Barcelona immunoassay), code 6 is ERR4 (Barcelona LC-MS/MS), and code 7 is ERR5 (Brussels immunoassay). The source paper does NOT state what separates code 1 from code 2 or code 3 from code 4 (most plausibly two immunoassay platforms within Rotterdam, and the Rotterdam vs Leiden split of the shared LC-MS/MS magnitude); because each pair shares a single estimated magnitude, the distinction is immaterial to the model's predictions, and a user may code either member of a pair. The control stream refers to this column as DATA in $ERROR while $INPUT lists CENTER; the published stream does not reconcile the two names. For a new single-centre LC-MS/MS dataset, code 3 (or 4) is the closest analogue.",
      source_name        = "DATA (referred to as CENTER in $INPUT)"
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
    renal_function = "Serum creatinine median 140 umol/L (IQR 110-191); GFR median 44 mL/min (IQR 31-56, unknown in 390 recipients).",
    genotypes      = "CYP3A5: *3/*3 884 (74.9%), *3/*1 224 (19.0%), *1/*1 39 (3.3%), Other 3 (0.3%), Unknown 30 (2.5%). CYP3A4: *1 1035 (87.7%), *22 112 (9.5%), Unknown 33 (2.8%).",
    notes          = "Baseline characteristics from Francke 2025 Table 1 (Total column, n = 1180). International, multicentre, retrospective analysis pooling five previously reported cohorts across the four centres; 13,427 tacrolimus concentrations in total, of which 11,670 were pre-dose concentrations and 1757 came from 208 densely sampled AUC profiles in 185 recipients. Median 10 samples per recipient (IQR 6-16). Because the dataset is dominated by pre-dose concentrations, ka, V1/F, and V2/F were additionally estimated using only the densely sampled subset; V2/F was FIXED at that subset's estimate in the reported models (Section 3.2.1)."
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural parameters -- Francke 2025 Table 2, "Final model" column.
    # Reference subject for the typical values: CYP3A5 *3/*3 non-expresser,
    # CYP3A4*1/*1, hematocrit 0.33 L/L, age 57.5 years, serum creatinine
    # 147 umol/L, height 170 cm (i.e. every covariate factor equals 1).
    # Apparent clearances (CL/F, Q/F) in L/h, apparent volumes (V1/F, V2/F)
    # in L, ka in 1/h, tlag in h. Values in parentheses after each estimate
    # are the %RSE reported by Table 2.
    # ---------------------------------------------------------------------
    ltlag <- log(0.375); label("Absorption lag time tlag (h)")                                     # Table 2 final model: t_lag = 0.375 h (RSE 5.3%; bootstrap 0.345, 95% CI 0.254-0.496)
    lka   <- log(6.59) ; label("First-order absorption rate constant ka (1/h)")                    # Table 2 final model: Ka = 6.59 (RSE 30%; bootstrap 6.27, 95% CI 2.16-11.01). The Table 2 column header reads "K a (L/h)" which is a typo; the Abstract states "the mean absorption rate was 6.59/h" and the parameter is a first-order rate constant with units 1/h.
    lcl   <- log(20.7) ; label("Apparent oral clearance CL/F at the reference covariate set (L/h)") # Table 2 final model: CL/F = 20.7 L/h (RSE 1.6%; bootstrap 20.6, 95% CI 20.0-21.4)
    lvc   <- log(705)  ; label("Apparent central volume V1/F (L)")                                 # Table 2 final model: V1/F = 705 L (RSE 4.7%; bootstrap 701.4, 95% CI 635.6-774.7)
    lq    <- log(8.54) ; label("Apparent inter-compartmental clearance Q/F (L/h)")                 # Table 2 final model: Q/F = 8.54 L/h (RSE 4.3%; bootstrap 8.55, 95% CI 7.75-9.32)

    # V2/F was NOT estimated on the full dataset. Section 3.2.1: "V2 was fixed
    # on the value that was estimated using the data of patients with at least
    # one more densely sampled curve available." Table 2 marks it "(f)" in every
    # column, and Supplementary Data S3 $THETA reads "(7670 FIX) ; 5 V2 - 7670".
    lvp   <- fixed(log(7670)); label("Apparent peripheral volume V2/F (L)")                 # Table 2 final model: V2/F = 7670 L (f); Supplementary Data S3 $THETA(5) = 7670 FIX

    # ---------------------------------------------------------------------
    # Covariate effects on CL/F -- Francke 2025 Table 2 "Covariate effect on
    # CL" rows, and the published final-model equation (Section 3.2.2):
    #
    #   CL/F = 20.7
    #          * [1.0 (CYP3A5 *3/*3) | 1.64 (*1/*3) | 1.93 (*1/*1)]
    #          * [1.0 (CYP3A4*1)     | 0.836 (CYP3A4*22 carrier)]
    #          * (Hematocrit/0.33)^-0.51
    #          * (Age/57.5)^-0.309
    #          * (Creatinine/147)^-0.0905
    #          * (Height/170)^1.17
    #
    # Continuous effects follow Eq 1 (median-normalised power model); the two
    # genotype effects follow Eq 2 (multiplicative categorical factor). The
    # centring constants are hard-coded in model() below with the same source
    # citations, matching Supplementary Data S3 $PK.
    # ---------------------------------------------------------------------
    e_cyp3a5_het_cl   <- 1.64   ; label("CYP3A5 *1/*3 multiplier on CL/F (vs the *3/*3 non-expresser reference)") # Table 2 final model: CYP3A5*1/*3 = 1.64 (RSE 3%; bootstrap 1.64, 95% CI 1.55-1.73); Supplementary Data S3 THETA(12)
    e_cyp3a5_hom_cl   <- 1.93   ; label("CYP3A5 *1/*1 multiplier on CL/F (vs the *3/*3 non-expresser reference)") # Table 2 final model: CYP3A5*1/*1 = 1.93 (RSE 5%; bootstrap 1.94, 95% CI 1.74-2.12); Supplementary Data S3 THETA(13)
    e_cyp3a4s22_cl    <- 0.836  ; label("CYP3A4*22 carrier multiplier on CL/F (vs the CYP3A4*1/*1 reference)")    # Table 2 final model: CYP3A4*22 = 0.836 (RSE 4%; bootstrap 0.838, 95% CI 0.77-0.90); Supplementary Data S3 THETA(18)
    e_hct_cl          <- -0.51  ; label("Hematocrit power exponent on CL/F, centred at 0.33 L/L (unitless)")      # Table 2 final model: Hematocrit (L/L) = -0.51 (RSE 10%; bootstrap -0.51, 95% CI -0.61 to -0.41); Supplementary Data S3 THETA(14)
    e_age_cl          <- -0.309 ; label("Age power exponent on CL/F, centred at 57.5 years (unitless)")           # Table 2 final model: Age (years) = -0.309 (RSE 13%; bootstrap -0.308, 95% CI -0.385 to -0.323); Supplementary Data S3 THETA(15)
    e_creat_cl        <- -0.0905; label("Serum creatinine power exponent on CL/F, centred at 147 umol/L (unitless)") # Table 2 final model: Creatinine (umol/L) = -0.0905 (RSE 21%; bootstrap -0.0915, 95% CI -0.1285 to -0.0525); Supplementary Data S3 THETA(16)
    e_ht_cl           <- 1.17   ; label("Height power exponent on CL/F, centred at 170 cm (unitless)")            # Table 2 final model: Height = 1.17 (RSE 18%; bootstrap 1.19, 95% CI 0.75-1.60); Supplementary Data S3 THETA(17). The Table 2 row label "(m 2)" is a typo -- the control stream centres at 170 and the Discussion states "height (in cm)".

    # ---------------------------------------------------------------------
    # Inter-individual variability. Supplementary Data S3 uses $OMEGA BLOCK(3)
    # with exponential IIV on CL (ETA(1)), V1 (ETA(2)), and Q (ETA(3)); ka and
    # tlag carry no IIV. Table 2 reports the diagonal as %CV and gives a
    # separate correlation matrix.
    #
    # Scale check: Table 2 reports base-model IIV CL/F 51.3% and final-model
    # 41.4%. Reading the column as %CV reproduces BOTH published summary
    # numbers: (51.3 - 41.4)/51.3 = 19.3% (the Abstract's "explained 19.3% of
    # the inter-individual variability") and (51.3^2 - 41.4^2)/51.3^2 = 34.9%
    # (the Key Points / Results / Conclusions "explained 35% of the
    # variability"). The column is therefore %CV, not a variance.
    #
    # Log-scale variances via omega^2 = log(1 + CV^2):
    #   CL/F CV 41.4% -> log(1 + 0.414^2) = 0.158196  (omega 0.397739)
    #   V1/F CV 80.5% -> log(1 + 0.805^2) = 0.499578  (omega 0.706808)
    #   Q/F  CV 82.8% -> log(1 + 0.828^2) = 0.522112  (omega 0.722573)
    # Off-diagonals are cov_ij = r_ij * omega_i * omega_j (eta correlations are
    # scale-free, so Table 2's correlation matrix applies directly on the log
    # scale):
    #   CL x V1  r = 0.580 -> 0.580 * 0.397739 * 0.706808 = 0.163053
    #   CL x Q   r = 0.065 -> 0.065 * 0.397739 * 0.722573 = 0.018681
    #   V1 x Q   r = 0.125 -> 0.125 * 0.706808 * 0.722573 = 0.063840
    # The resulting matrix is positive-definite (eigenvalues 0.6131, 0.4740,
    # 0.0927). The omega^2 = log(1 + CV^2) convention follows the
    # Andrews_2017_tacrolimus.R precedent from the same research group.
    # ---------------------------------------------------------------------
    # Block entries in lower-triangular order, all from Table 2 "Final model":
    #   [1,1] 0.158196 -- IIV CL/F 41.4% CV (RSE 3.1%; bootstrap 41.0%)
    #   [2,1] 0.163053 -- correlation CL/F x V1 = 58%
    #   [2,2] 0.499578 -- IIV V1/F 80.5% CV (RSE 3.7%; bootstrap 79.7%)
    #   [3,1] 0.018681 -- correlation CL/F x Q = 6.5%
    #   [3,2] 0.063840 -- correlation V1 x Q = 12.5%
    #   [3,3] 0.522112 -- IIV Q/F 82.8% CV (RSE 6%; bootstrap 83.5%)
    etalcl + etalvc + etalq ~ c(
      0.158196,
      0.163053, 0.499578,
      0.018681, 0.063840, 0.522112
    )

    # ---------------------------------------------------------------------
    # Residual variability. Section 2.2.1: a logarithmic error model was
    # selected. Supplementary Data S3 confirms the scale unambiguously:
    #   $ERROR   ERR1 = SQRT(THETA(7)**2)
    #            IF(F.GT.0) IPRED = LOG(F)
    #            Y = IPRED + ERR1*EPS(1)
    #   $SIGMA   1 FIX   (all five)
    # With EPS ~ N(0, 1) held FIXED, each THETA IS the residual standard
    # deviation on the natural-log scale -- these are SDs, not variances. That
    # maps onto nlmixr2's lnorm() residual family directly.
    #
    # Five magnitudes, selected by transplant centre x bioanalytical method
    # (STUDY_TACRO_FRANCKE, see covariateData).
    # ---------------------------------------------------------------------
    expSdRotterdamImmuno     <- 0.242; label("Log-scale residual SD, Rotterdam immunoassays")               # Table 2 final model: Rotterdam-immunoassays = 0.242 (RSE 3.4%; bootstrap 0.242, 95% CI 0.226-0.258); Supplementary Data S3 THETA(7)
    expSdRotterdamLeidenLcms <- 0.281; label("Log-scale residual SD, Rotterdam and Leiden LC-MS/MS")        # Table 2 final model: Rotterdam and Leiden-LC-MS/MS = 0.281 (RSE 2.3%; bootstrap 0.281, 95% CI 0.268-0.293); Supplementary Data S3 THETA(8)
    expSdBarcelonaImmuno     <- 0.376; label("Log-scale residual SD, Barcelona immunoassay")                # Table 2 final model: Barcelona-immunoassay = 0.376 (RSE 4.3%; bootstrap 0.374, 95% CI 0.344-0.407); Supplementary Data S3 THETA(9)
    expSdBarcelonaLcms       <- 0.419; label("Log-scale residual SD, Barcelona LC-MS/MS")                   # Table 2 final model: Barcelona-LC-MS/MS = 0.419 (RSE 5.9%; bootstrap 0.418, 95% CI 0.371-0.467); Supplementary Data S3 THETA(10)
    expSdBrusselsImmuno      <- 0.313; label("Log-scale residual SD, Brussels immunoassay")                 # Table 2 final model: Brussels-immunoassay = 0.313 (RSE 4%; bootstrap 0.312, 95% CI 0.288-0.338); Supplementary Data S3 THETA(11)
  })

  model({
    # 1. Covariate factors on CL/F. Written so that every factor equals 1 at
    #    the reference covariate set, matching Supplementary Data S3 $PK
    #    (CLCOV = CLCYP3A5*CLHCT*CLAGE*CLCREAT*CLHGT*CLCYP3A4).

    # CYP3A5 genotype, three levels encoded by two paired binary indicators
    # with *3/*3 as the implicit reference. Recipients with an unknown or
    # "Other" genotype are handled at the data-coding layer (see
    # covariateData$CYP3A5_STAR1_HET$notes), not here.
    f_cyp3a5 <- 1 + (e_cyp3a5_het_cl - 1) * CYP3A5_STAR1_HET +
                    (e_cyp3a5_hom_cl - 1) * CYP3A5_STAR1_HOM

    # CYP3A4*22 carriage, dominant model against the CYP3A4*1/*1 reference.
    f_cyp3a4 <- 1 + (e_cyp3a4s22_cl - 1) * SNP_CYP3A4_RS35599367

    # Median-normalised power effects (Eq 1). Centring constants are the
    # dataset medians hard-coded in Supplementary Data S3 $PK: 0.33 L/L,
    # 57.5 years, 147 umol/L, and 170 cm.
    f_hct   <- (HCT   / 0.33) ^ e_hct_cl
    f_age   <- (AGE   / 57.5) ^ e_age_cl
    f_creat <- (CREAT / 147)  ^ e_creat_cl
    f_ht    <- (HT    / 170)  ^ e_ht_cl

    # 2. Individual parameters. Bioavailability F was FIXED to 1 (Section
    #    2.2.1), so CL, Q, V1, and V2 are apparent (/F) quantities throughout.
    #    IIV is exponential on CL/F, V1/F, and Q/F only.
    tlag <- exp(ltlag)
    ka   <- exp(lka)
    cl   <- exp(lcl + etalcl) * f_cyp3a5 * f_cyp3a4 * f_hct * f_age * f_creat * f_ht
    vc   <- exp(lvc + etalvc)
    q    <- exp(lq  + etalq)
    vp   <- exp(lvp)

    # 3. Micro-constants, matching ADVAN4 TRANS4 (K = CL/V2, K23 = Q/V2,
    #    K32 = Q/V3 in the control stream's NONMEM compartment numbering,
    #    where NONMEM's V2 is the central and V3 the peripheral volume).
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 4. Two-compartment disposition with first-order oral absorption.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. Absorption lag time on the depot (ALAG1 in the control stream).
    alag(depot) <- tlag

    # 6. Observation. Whole-blood tacrolimus in ng/mL. Dose is in mg and vc in
    #    L, so central/vc is mg/L = ug/mL; multiply by 1000 for ng/mL. This
    #    reproduces the control stream's S2 = V2/1000 scaling exactly.
    Cc <- central / vc * 1000

    # Log-scale (exponential) residual error, magnitude selected per
    # observation by the centre x assay stratum. Codes 1-2 share one estimate
    # and codes 3-4 share another, per the $ERROR block.
    w_indiv <-
      expSdRotterdamImmuno     * ((STUDY_TACRO_FRANCKE == 1) + (STUDY_TACRO_FRANCKE == 2)) +
      expSdRotterdamLeidenLcms * ((STUDY_TACRO_FRANCKE == 3) + (STUDY_TACRO_FRANCKE == 4)) +
      expSdBarcelonaImmuno     * (STUDY_TACRO_FRANCKE == 5) +
      expSdBarcelonaLcms       * (STUDY_TACRO_FRANCKE == 6) +
      expSdBrusselsImmuno      * (STUDY_TACRO_FRANCKE == 7)

    Cc ~ lnorm(w_indiv)
  })
}
