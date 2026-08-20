Aruldhas_2021_S_methadone <- function() {
  description <- paste(
    "Population pharmacokinetic model for S-methadone and its S-EDDP metabolite",
    "in post-operative paediatric surgical patients (Aruldhas 2021; 61 children",
    "aged 11-17 years given 0.1 mg/kg IV racemic methadone intra-operatively",
    "followed by 0.1 mg/kg oral racemic methadone every 12 h as either a tablet",
    "or an oral suspension). S-methadone disposition is a two-compartment model",
    "with first-order absorption (separate tablet and suspension Ka values), an",
    "estimated absolute bioavailability F common to both oral formulations, and",
    "a one-compartment S-EDDP metabolite whose central volume is set equal to",
    "the S-methadone central volume (scaling factor VF fixed to 1 due to",
    "metabolite unidentifiability). Total S-methadone clearance CL scales",
    "allometrically with body weight (fixed exponent 0.75, reference 70 kg).",
    "The fractional S-methadone clearance to S-EDDP (CLF) is a linear function",
    "of the CYP2B6 activity score and the number of active alleles of the",
    "intronic CYP3A4 SNP rs2246709 (Table 2 final model chosen in preference",
    "to the alternative rs11882424 model because rs2246709 is more commonly",
    "genotyped in clinical laboratories). The S-methadone central volume V2",
    "is a linear function of the plasma alpha-1 acid glycoprotein (AAG)",
    "concentration and the number of active alleles of the ORM1 SNP rs17650.",
    "Between-subject variability is estimated on CL, V2, peripheral volume,",
    "CLF, and S-EDDP clearance; residual error is proportional for both",
    "S-methadone and S-EDDP."
  )
  reference <- paste(
    "Aruldhas BW, Quinney SK, Overholser BR, Heathman MA, Masters AR, Ly RC,",
    "Gao H, Packiasabapathy S, Sadhasivam S (2021). Pharmacokinetic modeling",
    "of R and S-Methadone and their metabolites to study the effects of",
    "various covariates in post-operative children.",
    "CPT Pharmacometrics Syst Pharmacol 10(10):1183-1194.",
    "doi:10.1002/psp4.12687",
    sep = " "
  )
  vignette <- "Aruldhas_2021_methadone"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    depot        = list(analyte = "methadone", units = "mg", specimen = "administration site", verified = FALSE),
    central      = list(analyte = "methadone", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1  = list(analyte = "methadone", units = "mg", specimen = "plasma", verified = FALSE),
    central_eddp = list(analyte = "S-EDDP", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight; allometrically scales S-methadone total clearance CL",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline. Reference 70 kg with fixed allometric exponent",
        "0.75 on S-methadone CL per Aruldhas 2021 Results 'Covariate modeling on",
        "S methadone' paragraph 1 (same wording as for R-methadone): the",
        "exponent was estimated to be 0.67 but was fixed to 0.75 to prevent",
        "small-cohort bias. Allometric scaling of the volume parameters was not",
        "retained; only CL carries the WT covariate. Cohort weight median 53.60",
        "kg (IQR 47.90-60.10)."
      ),
      source_name        = "WT"
    ),
    AAG = list(
      description        = "Plasma alpha-1 acid glycoprotein (AAG) concentration; time-varying covariate on S-methadone central volume V2",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying: measured at every methadone-concentration time point per",
        "Aruldhas 2021 Methods 'Covariate models' paragraph 1. Centered on the",
        "cohort average 94.76 ug/mL per Aruldhas 2021 Results 'Covariate",
        "modeling on S methadone' paragraph 3, so the reported V2 typical value",
        "applies at AAG = 94.76 ug/mL. Baseline cohort AAG median 84 ug/mL (IQR",
        "62-109). AAG concentrations increase post-surgery (peaking at 3-4",
        "days, returning to baseline by 2-4 weeks; Aruldhas 2021 Discussion",
        "paragraph 3); the median AAG at the end of the sampling window was",
        "115.75 ug/mL. Register canonical `AAG` documents units g/L; Aruldhas",
        "2021 reports AAG in ug/mL and uses ug/mL as the model unit throughout",
        "(1 g/L = 1000 ug/mL)."
      ),
      source_name        = "AAG"
    ),
    SNP_ORM1_RS17650 = list(
      description        = "Number of active alleles of the ORM1 SNP rs17650 (0, 1, or 2); time-fixed genotype covariate on S-methadone central volume V2",
      units              = "(count, 0/1/2 alleles per subject)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Centered on n_active = 1 (heterozygous), so the reported V2 typical",
        "value applies at a heterozygous rs17650 genotype. rs17650 distinguishes",
        "the F (fast-migrating; variant) and S (slow-migrating; reference)",
        "allozymes of AAG encoded by ORM1; the F allozyme binds methadone with",
        "lower affinity than the S allozyme (Aruldhas 2021 Discussion",
        "paragraph 3). Effect is independent of AAG concentration. Cohort",
        "distribution: 7 wild-type / 30 heterozygous / 15 homozygous per Table",
        "S3."
      ),
      source_name        = "rs17650"
    ),
    CYP2B6 = list(
      description        = "CYP2B6 individual metabolic-activity score (five-level: 0, 0.5, 1, 1.5, 2 for PM / IM / NM / RM / UM); time-fixed genotype-derived phenotype covariate on the fractional S-methadone clearance to S-EDDP",
      units              = "(activity score, 0 to 2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Derived from Aldy-imputed *-haplotype diplotypes using the CPIC",
        "efavirenz guideline mapping (poor metabolizer = 0, intermediate = 0.5,",
        "normal / extensive = 1, rapid = 1.5, ultra-rapid = 2). Per-diplotype",
        "assignments for the Aruldhas 2021 cohort are in Supplemental Table S1.",
        "Centered on activity score 1 (normal metabolizer), so the reported CLF",
        "typical value applies at CYP2B6 = 1. Cohort distribution (from Table",
        "S1, n = 52 with genotype available): NM 29 (55.8%), IM 18 (34.6%), RM",
        "2 (3.8%), PM 3 (5.8%). Nine subjects (14.8%) with missing genotype had",
        "the activity score imputed as the cohort mode (NM = 1) per Aruldhas",
        "2021 Results 'Covariate modeling on S methadone' paragraph 4."
      ),
      source_name        = "CYP2B6 activity score"
    ),
    SNP_CYP3A4_RS2246709 = list(
      description        = "Number of active alleles of the intronic CYP3A4 SNP rs2246709 (0, 1, or 2); time-fixed genotype covariate on the fractional S-methadone clearance to S-EDDP",
      units              = "(count, 0/1/2 alleles per subject)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Centered on n_active = 1 (heterozygous), so the reported CLF typical",
        "value applies at heterozygous rs2246709 with CYP2B6 activity score 1.",
        "Aruldhas 2021 tested both an intronic CYP2B6 rs11882424 covariate",
        "(dOFV = -13.2) and this CYP3A4 rs2246709 covariate (dOFV = -9.0) on",
        "S-methadone CLF; the rs2246709 model was chosen as the final S-",
        "methadone model in preference to rs11882424 because rs2246709 is more",
        "commonly genotyped in clinical laboratories (Aruldhas 2021 Results",
        "'Covariate modeling on S methadone' paragraph 3 and Discussion",
        "paragraph 5). Cohort distribution: 23 wild-type / 24 heterozygous / 5",
        "homozygous per Table S3. Note: the very-large S-methadone coefficient",
        "(1.68 with RSE 57.8%) implies that the paper's model can produce a",
        "physically-invalid (negative) CLF for S-methadone at n_active = 0",
        "(wild-type homozygous) subjects: 1 + 1.68 * (0 - 1) = -0.68. This is",
        "a limitation of the paper's parameter estimates at the extreme of the",
        "covariate range; users simulating a wild-type-homozygous rs2246709",
        "subject should either restrict simulation to n_active >= 1 or accept",
        "the paper's noted limitation (documented in the vignette Assumptions",
        "and deviations)."
      ),
      source_name        = "rs2246709"
    ),
    FORM_TABLET = list(
      description        = "Oral formulation indicator (1 = tablet, 0 = oral suspension); per-dose covariate that selects the S-methadone absorption rate constant Ka",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (oral suspension; the reference formulation for the canonical Ka value)",
      notes              = paste(
        "Per-oral-dose indicator: 1 = tablet, 0 = suspension. Selects between",
        "the two independently-estimated absorption rate constants Ka_tab and",
        "Ka_susp reported in Aruldhas 2021 Table 2 (Ka_tab = 0.257 h^-1,",
        "Ka_susp = 0.432 h^-1 for S-methadone). Estimating separate Ka values",
        "for the two oral formulations improved the fit vs a single Ka (dOFV =",
        "-3.9; Results 'Structural model of S methadone' paragraph 1). Register",
        "canonical `FORM_TABLET` was originally introduced for tablet vs",
        "solution studies; per its Notes, per-model documentation should",
        "identify the non-tablet comparator here (suspension) as the reference",
        "(FORM_TABLET = 0). The paper's bioavailability F is common to both",
        "formulations (individual bioavailabilities were not identifiable and",
        "were assumed equivalent). Does not affect IV doses (Ka only applies to",
        "oral doses through the depot compartment)."
      ),
      source_name        = "formulation code (paper narrative; tablet vs suspension per-dose)"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 61L,
    n_studies        = 1L,
    age_range        = "11-17 years (paper Results paragraph 1 reports 11-17 years even though the eligibility range was 8-17.9 years)",
    age_median       = "14.74 years (IQR 13.62-15.66)",
    weight_range     = "cohort median 53.60 kg (IQR 47.90-60.10)",
    weight_median    = "53.60 kg",
    height_median    = "164.50 cm (IQR 158.00-171.50)",
    bmi_median       = "19.40 (IQR 17.61-22.50)",
    aag_baseline     = "median 84 ug/mL (IQR 62-109)",
    sex_female_pct   = 50.8,
    race_ethnicity   = c(White = 80.3, `African American` = 11.5, Hispanic = 4.9, `Native American` = 1.6, Unknown = 1.6),
    disease_state    = paste(
      "Children and adolescents undergoing pectus excavatum repair (n = 25) or",
      "posterior spinal fusion surgery for idiopathic scoliosis (n = 36).",
      "Exclusions: methadone allergy, developmental delay, neurologic disorder,",
      "liver or renal disease, pre-operative pain requiring analgesics."
    ),
    dose_range       = paste(
      "0.1 mg/kg racemic methadone IV intra-operatively (first dose), then 0.1",
      "mg/kg oral racemic methadone (either tablet or suspension per patient",
      "convenience) every 12 h post-operatively for four to six additional",
      "doses. Cohort median actual dose 0.087 mg/kg (IQR 0.069-0.094). The",
      "PK measurements below (5-8 samples per patient across the first three",
      "interdose intervals) inform the individual S-methadone and S-EDDP",
      "profiles."
    ),
    regions          = "United States (Indiana University School of Medicine, IRB #1707525204, ClinicalTrials.gov NCT03495388)",
    notes            = paste(
      "Racemic methadone dose is split 50/50 between the R and S enantiomers.",
      "This model file describes only the S-methadone enantiomer + its S-EDDP",
      "metabolite; the R-methadone + R-EDDP is in the sibling file",
      "Aruldhas_2021_R_methadone.R. To reproduce a racemic-dose simulation,",
      "dose both files with half of the racemic amount (mg of S-methadone in",
      "this file). A total of 430 methadone and 381 EDDP plasma concentration",
      "measurements were available for analysis after excluding the 0.46%",
      "methadone and 9.4% EDDP samples below the LLOQ (0.015 ng/mL; M1 method).",
      "Nine of the 61 children (14.8%) had missing genotype information",
      "imputed as the cohort mode. Analytical assay: liquid chromatography",
      "tandem mass spectrometry (linear 0.015-150 ng/mL; CV < 15%). AAG",
      "quantified by HPLC-UV (linear 20-1500 ug/mL; CV < 10%). Estimation:",
      "NONMEM 7.4 first-order conditional with eta-epsilon interaction",
      "(FOCE-I). Robustness assessed by bootstrap (n = 1000) and prediction-",
      "corrected VPC (n = 1000)."
    )
  )

  # Implementation notes (see vignette 'Assumptions and deviations' for the
  # full justification of each item):
  #
  # * Structure. Two-compartment S-methadone parent + one-compartment S-EDDP
  #   metabolite. Parent compartments are `central` (V2 = 98.3 L) and
  #   `peripheral1` (V4). The S-EDDP central compartment `central_eddp` uses
  #   the registered `eddp` metabolite suffix; its volume V3 = V2 * VF with
  #   VF fixed to 1 due to metabolite unidentifiability (paper Figure 1
  #   legend and Results 'Structural model of S methadone' paragraph 2).
  #
  # * V4 (parent peripheral) volume. As with the R-methadone sibling model,
  #   the paper's Table 2 row labelled "V3 (L) = 139 (RSE 19.3%)" with
  #   116% BSV for S-methadone is interpreted as V4 (parent peripheral
  #   volume) rather than V3 (metabolite central volume); V3 = V2 * VF with
  #   VF = 1 cannot carry independent typical value or BSV.
  #
  # * V2rs17650 coefficient sign. Aruldhas 2021 Table 2 lists the S-methadone
  #   V2rs17650 coefficient as "0.526 (28.4)" with 95% CI "-0.709 to -0.173"
  #   -- the CI is entirely negative but the point estimate is printed as
  #   positive. The R-methadone sibling reports -0.443 with a matching
  #   negative CI, and the paper text uniformly describes rs17650 as
  #   decreasing V2 (F-allozyme = lower binding affinity). This model uses
  #   -0.526 for the S-methadone V2rs17650 coefficient, treating the Table
  #   2 positive sign as a table transcription typo consistent with the
  #   paper's negative confidence interval.
  #
  # * S-methadone CLF final-model selection. Aruldhas 2021 evaluated two
  #   candidate models for the CYP-covariate effect on S-methadone CLF:
  #   (1) rs11882424 (dOFV = -13.2) and (2) rs2246709 (dOFV = -9.0). Model
  #   (2) was chosen as the "final model for the subsequent analyses"
  #   because rs2246709 is more commonly genotyped in clinical laboratories
  #   (paper Results 'Covariate modeling on S methadone' paragraph 3). Table
  #   2 reports the model-(2) parameter estimates and this extraction uses
  #   the model-(2) form.
  #
  # * Absorption. Two independently-estimated first-order Ka values, one
  #   per oral formulation, selected via the FORM_TABLET binary covariate.
  #   Suspension is the canonical Ka reference (0.432 h^-1) and the tablet
  #   is a log-additive multiplicative shift (log(0.257 / 0.432) = -0.520;
  #   giving 0.257 h^-1 when FORM_TABLET = 1). Paper text: dOFV = -3.9 for
  #   the two-Ka fit vs a shared Ka (Results 'Structural model of S
  #   methadone' paragraph 1).
  #
  # * Fractional metabolite formation clearance CLF. Same dimensionless-
  #   fraction interpretation as the R-methadone sibling model (see the
  #   Aruldhas_2021_R_methadone.R implementation notes for the full
  #   derivation). With CLF_ref = 0.135 for S-methadone (at CYP2B6 = 1 and
  #   rs2246709 = 1) and CL = 13.0 L/h, the typical parent-to-EDDP
  #   formation clearance is 1.76 L/h and the other-routes clearance is
  #   11.24 L/h.
  #
  # * S-methadone CLF at wild-type rs2246709. The very-large rs2246709
  #   coefficient (1.68) implies that CLF becomes negative at n_active = 0:
  #     CLF = 0.135 * (1 + 0.636 * (CYP2B6 - 1)) * (1 + 1.68 * (0 - 1))
  #         = 0.135 * (CYP2B6 term) * (-0.68)
  #   which is physically invalid. Users simulating a wild-type-homozygous
  #   rs2246709 subject should either restrict simulation to n_active >= 1
  #   or use the sibling R-methadone model as a proxy. This is a limitation
  #   of the paper's parameter estimates at the extreme of the covariate
  #   range, not of the extraction; documented in the vignette Assumptions
  #   and deviations.
  #
  # * Between-subject variability. NONMEM reports OMEGA on the internal
  #   log scale; Aruldhas 2021 Table 2 reports BSV as CV%. The packaged
  #   model uses the exact log-normal conversion omega^2 = log(CV^2 + 1)
  #   so the simulated cross-subject geometric coefficient of variation
  #   matches the reported value.
  #
  # * Residual error. Aruldhas 2021 reports proportional-only residual
  #   error (RUV drug = 0.165, RUV metabolite = 0.194; Table 2 for
  #   S-methadone), which maps directly to nlmixr2's `prop(propSd)`.
  ini({
    # ----- Absorption (Aruldhas 2021 Table 2, S-methadone) -----
    lka              <- log(0.432);          label("Suspension first-order absorption rate constant (1/h; canonical Ka, reference formulation)")  # Aruldhas 2021 Table 2 'Ka susp (h^-1) = 0.432 (22.8)' for S-methadone
    e_form_tablet_ka <- log(0.257 / 0.432);  label("Log-additive multiplicative shift on Ka for the tablet formulation (unitless)")               # Aruldhas 2021 Table 2 'Ka tablet (h^-1) = 0.257 (42.4)' / 'Ka susp (h^-1) = 0.432' for S-methadone; log(0.257/0.432) = -0.520 reproduces the 41% slower tablet absorption

    # ----- Parent S-methadone disposition (Aruldhas 2021 Table 2) -----
    lcl              <- log(13.0);           label("Typical S-methadone total apparent clearance CL (L/h) at 70 kg, CYP2B6 activity score 1, and heterozygous rs2246709")  # Aruldhas 2021 Table 2 'CL (L*h^-1) = 13.0 (16.7)' for S-methadone
    lvc              <- log(98.3);           label("Typical S-methadone central-compartment apparent volume V2 (L) at heterozygous rs17650 and cohort-average AAG 94.76 ug/mL")  # Aruldhas 2021 Table 2 'V2 (L) = 98.3 (12.9)' for S-methadone
    lq               <- log(105);            label("Typical S-methadone intercompartmental clearance Q (L/h)")                                    # Aruldhas 2021 Table 2 'Q (L*h^-1) = 105 (18.9)' for S-methadone
    lvp              <- log(139);            label("Typical S-methadone peripheral-compartment apparent volume V4 (L)")                            # Aruldhas 2021 Table 2 row labelled 'V3 (L) = 139 (19.3)' with 116% BSV for S-methadone; reassigned to V4 (parent peripheral) per the paper's Figure 1 legend and Results 'Structural model of S methadone' paragraph 1 (V3 = V2 * VF cannot carry an independent typical value or BSV when VF = 1 fixed; see implementation notes above)
    lfdepot          <- log(0.606);          label("Oral bioavailability F (fraction; common to tablet and suspension formulations)")             # Aruldhas 2021 Table 2 'F = 0.606 (14.8)' for S-methadone

    # ----- S-methadone -> S-EDDP fractional metabolite-formation clearance CLF (dimensionless) -----
    lclf             <- log(0.135);          label("Typical fraction of S-methadone total apparent clearance CL that becomes S-EDDP metabolite (dimensionless) at CYP2B6 activity score 1 and heterozygous rs2246709")  # Aruldhas 2021 Table 2 'CLF = 0.135 (23.7)' for S-methadone

    # ----- Metabolite (S-EDDP) disposition (Aruldhas 2021 Table 2) -----
    lcl_eddp         <- log(7.97);           label("Typical S-EDDP apparent clearance CL3 (L/h)")                                                  # Aruldhas 2021 Table 2 'CL3 (L*h^-1) = 7.97 (23.7)' for S-methadone
    lvf              <- fixed(log(1));       label("VF scaling factor between S-methadone V2 and S-EDDP V3 (V3 = V2 * VF; dimensionless, 1 due to metabolite unidentifiability)")  # Aruldhas 2021 Table 2 'VF = 1 FIX' and Figure 1 legend

    # ----- Allometric exponent (fixed) -----
    e_wt_cl                      <- fixed(0.75);  label("Allometric WT exponent on S-methadone CL (unitless)")                                    # Aruldhas 2021 Results 'Covariate modeling on S methadone' paragraph 1: fixed to 0.75 to avoid small-cohort bias despite the estimated exponent (0.67)

    # ----- Covariate effects on V2 (linear centered-additive multiplicative form) -----
    e_snp_orm1_rs17650_vc        <- -0.526;       label("Linear covariate coefficient of ORM1 rs17650 active-allele count on S-methadone V2 (unitless per active allele)")  # Aruldhas 2021 Table 2 'V2rs17650 = 0.526 (28.4)' with 95% CI (-0.709 to -0.173); the sign of the point estimate is treated as a Table 2 transcription typo (see implementation notes above); we use -0.526 consistent with the negative CI and the paper text
    e_aag_vc                     <- -0.00192;     label("Linear covariate coefficient of AAG concentration (centered at 94.76 ug/mL) on S-methadone V2 (unitless per ug/mL)")  # Aruldhas 2021 Table 2 'V2AAG = -0.00192 (51.8)' for S-methadone; V2 = V2_typ * (1 - 0.00192 * (AAG - 94.76))

    # ----- Covariate effects on CLF (linear centered-additive multiplicative form) -----
    e_cyp2b6_clf                 <- 0.636;        label("Linear covariate coefficient of CYP2B6 activity score on S-methadone CLF (unitless per activity-score unit)")  # Aruldhas 2021 Table 2 'CLFCYP2B6 = 0.636 (25.2)' for S-methadone; CLF = CLF_typ * (1 + 0.636 * (CYP2B6 - 1))
    e_snp_cyp3a4_rs2246709_clf   <- 1.68;         label("Linear covariate coefficient of CYP3A4 rs2246709 active-allele count on S-methadone CLF (unitless per active allele)")  # Aruldhas 2021 Table 2 'CLFrs2246709 = 1.68 (57.8)' for S-methadone; CLF = CLF_typ * (1 + 1.68 * (rs2246709 - 1)); can produce negative CLF at n_active = 0 (see implementation notes above)

    # ----- Between-subject variability (Aruldhas 2021 Table 2) -----
    # Exact log-normal conversion omega^2 = log(CV^2 + 1).
    etalcl           ~ 0.15361   # Aruldhas 2021 Table 2 'CL BSV %CV = 40.9' for S-methadone; log(1 + 0.409^2) = 0.15361
    etalvc           ~ 0.33170   # Aruldhas 2021 Table 2 'V2 BSV %CV = 63.6' for S-methadone; log(1 + 0.636^2) = 0.33170
    etalvp           ~ 0.79102   # Aruldhas 2021 Table 2 'V4 BSV %CV = 116' (labelled 'V3' in Table 2, see implementation notes) for S-methadone; log(1 + 1.16^2) = 0.79102
    etalclf          ~ 0.20505   # Aruldhas 2021 Table 2 'CLF BSV %CV = 47.9' for S-methadone; log(1 + 0.479^2) = 0.20505
    etalcl_eddp      ~ 0.10689   # Aruldhas 2021 Table 2 'CL3 BSV %CV = 33.7' for S-methadone; log(1 + 0.337^2) = 0.10689

    # ----- Residual error (Aruldhas 2021 Table 2) -----
    propSd           <- 0.165;   label("Proportional residual SD for S-methadone plasma concentration (fraction)")  # Aruldhas 2021 Table 2 'RUV drug = 0.165 (6.15)' for S-methadone
    propSd_eddp      <- 0.194;   label("Proportional residual SD for S-EDDP plasma concentration (fraction)")       # Aruldhas 2021 Table 2 'RUV metabolite = 0.194 (5.47)' for S-methadone
  })
  model({
    # ----- Reference values (Aruldhas 2021 covariate-model equations) -----
    ref_wt  <- 70      # kg, standard allometric reference
    ref_aag <- 94.76   # ug/mL, cohort average per Aruldhas 2021 Results 'Covariate modeling on S methadone' paragraph 3

    # ----- Individual parameters -----
    ka      <- exp(lka + e_form_tablet_ka * FORM_TABLET)
    cl      <- exp(lcl + etalcl) * (WT / ref_wt)^e_wt_cl
    vc      <- exp(lvc + etalvc) *
               (1 + e_snp_orm1_rs17650_vc * (SNP_ORM1_RS17650 - 1)) *
               (1 + e_aag_vc              * (AAG - ref_aag))
    q       <- exp(lq)
    vp      <- exp(lvp + etalvp)
    fdepot  <- exp(lfdepot)
    clf     <- exp(lclf + etalclf) *
               (1 + e_cyp2b6_clf               * (CYP2B6 - 1)) *
               (1 + e_snp_cyp3a4_rs2246709_clf * (SNP_CYP3A4_RS2246709 - 1))
    cl_eddp <- exp(lcl_eddp + etalcl_eddp)
    vf      <- exp(lvf)
    vc_eddp <- vc * vf
    cl_met_form <- cl * clf

    # ----- ODE system -----
    d/dt(depot)        <- -ka * depot
    d/dt(central)      <-  ka * depot -
                           cl * (central / vc) -
                           q  * (central / vc - peripheral1 / vp)
    d/dt(peripheral1)  <-  q  * (central / vc - peripheral1 / vp)
    d/dt(central_eddp) <-  cl_met_form * (central / vc) -
                           cl_eddp     * (central_eddp / vc_eddp)

    f(depot) <- fdepot

    # ----- Observations -----
    # Dose in mg / volumes in L => concentrations in mg/L (= ug/mL).
    Cc       <- central      / vc
    Cc_eddp  <- central_eddp / vc_eddp

    Cc       ~ prop(propSd)
    Cc_eddp  ~ prop(propSd_eddp)
  })
}
