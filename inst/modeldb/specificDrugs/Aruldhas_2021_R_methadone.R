Aruldhas_2021_R_methadone <- function() {
  description <- paste(
    "Population pharmacokinetic model for R-methadone and its R-EDDP metabolite",
    "in post-operative paediatric surgical patients (Aruldhas 2021; 61 children",
    "aged 11-17 years given 0.1 mg/kg IV racemic methadone intra-operatively",
    "followed by 0.1 mg/kg oral racemic methadone every 12 h as either a tablet",
    "or an oral suspension). R-methadone disposition is a two-compartment model",
    "with first-order absorption (separate tablet and suspension Ka values), an",
    "estimated absolute bioavailability F common to both oral formulations, and",
    "a one-compartment R-EDDP metabolite whose central volume is set equal to",
    "the R-methadone central volume (scaling factor VF fixed to 1 due to",
    "metabolite unidentifiability). Total R-methadone clearance CL scales",
    "allometrically with body weight (fixed exponent 0.75, reference 70 kg).",
    "The fractional R-methadone clearance to R-EDDP (CLF) is a linear function",
    "of the CYP2B6 activity score and the number of active alleles of the",
    "intronic CYP3A4 SNP rs2246709. The R-methadone central volume V2 is a",
    "linear function of the plasma alpha-1 acid glycoprotein (AAG)",
    "concentration and the number of active alleles of the ORM1 SNP rs17650.",
    "Between-subject variability is estimated on CL, V2, peripheral volume,",
    "CLF, and R-EDDP clearance; residual error is proportional for both",
    "R-methadone and R-EDDP."
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
    depot        = list(analyte = "R-methadone", units = "mg", specimen = "administration site", verified = FALSE),
    central      = list(analyte = "R-methadone", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1  = list(analyte = "R-methadone", units = "mg", specimen = "plasma", verified = FALSE),
    central_eddp = list(analyte = "R-EDDP", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight; allometrically scales R-methadone total clearance CL",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline. Reference 70 kg with fixed allometric exponent",
        "0.75 on R-methadone CL per Aruldhas 2021 Results 'Covariate modeling on",
        "R methadone' paragraph 1: 'Bodyweight was added allometrically to the",
        "clearance of the drug with a fixed exponent of 0.75.' Estimating the",
        "exponent improved the fit (1.2) but was fixed to 0.75 to prevent bias",
        "in a small-cohort estimation. Allometric scaling of the volume",
        "parameters worsened the fit (dOFV = 1.2) and was not retained; only CL",
        "carries the WT covariate. Cohort weight median 53.60 kg (IQR 47.90-",
        "60.10)."
      ),
      source_name        = "WT"
    ),
    AAG = list(
      description        = "Plasma alpha-1 acid glycoprotein (AAG) concentration; time-varying covariate on R-methadone central volume V2",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying: measured at every methadone-concentration time point per",
        "Aruldhas 2021 Methods 'Covariate models' paragraph 1. Centered on the",
        "cohort average 94.76 ug/mL per Aruldhas 2021 Results 'Covariate",
        "modeling on R methadone' paragraph 3, so the reported CV V2 typical",
        "value applies at AAG = 94.76 ug/mL. Baseline cohort AAG median 84 ug/",
        "mL (IQR 62-109). AAG concentrations increase post-surgery (peaking at",
        "3-4 days, returning to baseline by 2-4 weeks; Aruldhas 2021 Discussion",
        "paragraph 3); the median AAG at the end of the sampling window was",
        "115.75 ug/mL. Register canonical `AAG` documents units g/L; Aruldhas",
        "2021 reports AAG in ug/mL and uses ug/mL as the model unit throughout",
        "(1 g/L = 1000 ug/mL)."
      ),
      source_name        = "AAG"
    ),
    SNP_ORM1_RS17650 = list(
      description        = "Number of active alleles of the ORM1 SNP rs17650 (0, 1, or 2); time-fixed genotype covariate on R-methadone central volume V2",
      units              = "(count, 0/1/2 alleles per subject)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Centered on n_active = 1 (heterozygous), so the reported V2 typical",
        "value applies at a heterozygous rs17650 genotype. rs17650 distinguishes",
        "the F (fast-migrating; variant) and S (slow-migrating; reference)",
        "allozymes of AAG encoded by ORM1; the F allozyme binds methadone with",
        "lower affinity than the S allozyme (Aruldhas 2021 Discussion",
        "paragraph 3, citing Herve et al. 1998 [ref 10]). Effect is independent",
        "of AAG concentration (Aruldhas 2021 Results 'Covariate modeling on R",
        "methadone' paragraph 2). Cohort distribution: 7 wild-type / 30",
        "heterozygous / 15 homozygous per Table S3."
      ),
      source_name        = "rs17650"
    ),
    CYP2B6 = list(
      description        = "CYP2B6 individual metabolic-activity score (five-level: 0, 0.5, 1, 1.5, 2 for PM / IM / NM / RM / UM); time-fixed genotype-derived phenotype covariate on the fractional R-methadone clearance to R-EDDP",
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
        "2021 Results 'Covariate modeling on R methadone' paragraph 4."
      ),
      source_name        = "CYP2B6 activity score"
    ),
    SNP_CYP3A4_RS2246709 = list(
      description        = "Number of active alleles of the intronic CYP3A4 SNP rs2246709 (0, 1, or 2); time-fixed genotype covariate on the fractional R-methadone clearance to R-EDDP",
      units              = "(count, 0/1/2 alleles per subject)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Centered on n_active = 1 (heterozygous), so the reported CLF typical",
        "value applies at heterozygous rs2246709 with CYP2B6 activity score 1.",
        "Aruldhas 2021 identified rs2246709 as significantly associated with",
        "reduced R-methadone fractional clearance to R-EDDP after backward",
        "elimination (dOFV = -7.1; Table S2). Cohort distribution: 23 wild-type",
        "/ 24 heterozygous / 5 homozygous per Table S3."
      ),
      source_name        = "rs2246709"
    ),
    FORM_TABLET = list(
      description        = "Oral formulation indicator (1 = tablet, 0 = oral suspension); per-dose covariate that selects the R-methadone absorption rate constant Ka",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (oral suspension; the reference formulation for the canonical Ka value)",
      notes              = paste(
        "Per-oral-dose indicator: 1 = tablet, 0 = suspension. Selects between",
        "the two independently-estimated absorption rate constants Ka_tab and",
        "Ka_susp reported in Aruldhas 2021 Table 2 (Ka_tab = 0.123 h^-1,",
        "Ka_susp = 0.318 h^-1 for R-methadone). Estimating separate Ka values",
        "for the two oral formulations improved the fit vs a single Ka (dOFV =",
        "-8.3; Results 'Structural model of R methadone' paragraph 1). Register",
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
      "interdose intervals) inform the individual R-methadone and R-EDDP",
      "profiles."
    ),
    regions          = "United States (Indiana University School of Medicine, IRB #1707525204, ClinicalTrials.gov NCT03495388)",
    notes            = paste(
      "Racemic methadone dose is split 50/50 between the R and S enantiomers.",
      "This model file describes only the R-methadone enantiomer + its R-EDDP",
      "metabolite; the S-methadone + S-EDDP is in the sibling file",
      "Aruldhas_2021_S_methadone.R. To reproduce a racemic-dose simulation,",
      "dose both files with half of the racemic amount (mg of R-methadone in",
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
  # * Structure. Two-compartment R-methadone parent + one-compartment R-EDDP
  #   metabolite. Parent compartments are `central` (V2 = 176 L) and
  #   `peripheral1` (V4). The R-EDDP central compartment `central_eddp` uses
  #   the registered `eddp` metabolite suffix; its volume V3 = V2 * VF with
  #   VF fixed to 1 due to metabolite unidentifiability (paper Figure 1
  #   legend and Results 'Structural model of R methadone' paragraph 2).
  #
  # * V4 (parent peripheral) volume. The paper's Figure 1 legend explicitly
  #   names V4 as the parent peripheral-compartment volume, and the model
  #   description names V4 as one of the three parameters with estimated
  #   BSV. Paper Table 2 lists a row labelled "V3 (L) = 335 (RSE 23.6%)"
  #   with 62.23% BSV for R-methadone. This value is interpreted as V4
  #   (parent peripheral volume) rather than V3 (metabolite central volume)
  #   because (a) the paper text explicitly names V4 as the peripheral
  #   volume with BSV, (b) V3 = V2 * VF with VF = 1 fixed cannot carry an
  #   independent typical value or independent BSV, and (c) the metabolite
  #   BSV would flow through the metabolite clearance (CL3 already carries
  #   its own BSV in Table 2). The 335 L "V3" therefore is Vp for the
  #   R-methadone peripheral compartment.
  #
  # * Absorption. Two independently-estimated first-order Ka values, one
  #   per oral formulation, selected via the FORM_TABLET binary covariate.
  #   Suspension is the canonical Ka reference (0.318 h^-1) and the tablet
  #   is a log-additive multiplicative shift (log(0.123 / 0.318) = -0.951;
  #   giving 0.123 h^-1 when FORM_TABLET = 1). Paper text: dOFV = -8.3 for
  #   the two-Ka fit vs a shared Ka (Results 'Structural model of R
  #   methadone' paragraph 1). Bioavailability F common to both
  #   formulations because individual formulation-specific F values could
  #   not be identified.
  #
  # * Fractional metabolite formation clearance CLF. The paper's CLF is
  #   a dimensionless fraction of the total parent CL that becomes EDDP,
  #   despite being reported in Table 2 with an apparent "L/h" unit (paper
  #   Figure 1 caption: "CLF is the fraction of CL to EDDP (CL2)"; also
  #   Table 2 note: "CLF was the fraction of CL that contributes to CL2,
  #   the clearance toward the formation of the metabolite"). In model()
  #   the parent-to-metabolite formation clearance is `cl * clf` (mass per
  #   time) and the other-routes elimination is `cl * (1 - clf)`; the
  #   parent's total elimination is `cl * Cp` regardless of the CLF
  #   partition, and only the metabolite formation flux depends on CLF.
  #   With CLF_ref = 0.217 for R-methadone (at CYP2B6 = 1 and rs2246709 =
  #   1) and CL = 15.7 L/h, the typical parent-to-EDDP formation clearance
  #   is 3.4 L/h and the other-routes clearance is 12.3 L/h; the EDDP
  #   contribution matches the paper's stated CYP-mediated metabolism
  #   qualitative attribution while the "other routes" fraction absorbs
  #   both non-EDDP metabolism and renal / other elimination pathways.
  #
  # * Covariate parameterisations. All covariate effects use the standard
  #   centered-additive multiplicative form
  #     P_i = P_ref * (1 + Theta_cov * (cov_i - cov_ref))
  #   with the reference values chosen so P_ref is the paper's reported
  #   typical value at the reference covariate combination (see the CL,
  #   V2, and CLF equations in Results 'Covariate modeling on R
  #   methadone' paragraphs 2-3). The AAG reference is 94.76 ug/mL (the
  #   cohort average per the paper's equation). The rs17650 and rs2246709
  #   references are n_active = 1 (heterozygous). The CYP2B6 reference is
  #   activity score 1 (normal metabolizer).
  #
  # * Between-subject variability. NONMEM reports OMEGA on the internal
  #   log scale; Aruldhas 2021 Table 2 reports BSV as CV%. The packaged
  #   model uses the exact log-normal conversion omega^2 = log(CV^2 + 1)
  #   so the simulated cross-subject geometric coefficient of variation
  #   matches the reported value.
  #
  # * Residual error. Aruldhas 2021 reports proportional-only residual
  #   error (RUV drug = 0.165, RUV metabolite = 0.207; Table 2 for
  #   R-methadone), which maps directly to nlmixr2's `prop(propSd)`.
  ini({
    # ----- Absorption (Aruldhas 2021 Table 2, R-methadone) -----
    lka              <- log(0.318);          label("Suspension first-order absorption rate constant (1/h; canonical Ka, reference formulation)")  # Aruldhas 2021 Table 2 'Ka susp (h^-1) = 0.318 (28.3)' for R-methadone
    e_form_tablet_ka <- log(0.123 / 0.318);  label("Log-additive multiplicative shift on Ka for the tablet formulation (unitless)")               # Aruldhas 2021 Table 2 'Ka tablet (h^-1) = 0.123 (47.2)' / 'Ka susp (h^-1) = 0.318' for R-methadone; log(0.123/0.318) = -0.951 reproduces the 61% slower tablet absorption

    # ----- Parent R-methadone disposition (Aruldhas 2021 Table 2) -----
    lcl              <- log(15.7);           label("Typical R-methadone total apparent clearance CL (L/h) at 70 kg, CYP2B6 activity score 1, and heterozygous rs2246709")  # Aruldhas 2021 Table 2 'CL (L*h^-1) = 15.7 (27.1)' for R-methadone
    lvc              <- log(176);            label("Typical R-methadone central-compartment apparent volume V2 (L) at heterozygous rs17650 and cohort-average AAG 94.76 ug/mL")  # Aruldhas 2021 Table 2 'V2 (L) = 176 (16.6)' for R-methadone
    lq               <- log(69.2);           label("Typical R-methadone intercompartmental clearance Q (L/h)")                                    # Aruldhas 2021 Table 2 'Q (L*h^-1) = 69.2 (39.3)' for R-methadone
    lvp              <- log(335);            label("Typical R-methadone peripheral-compartment apparent volume V4 (L)")                            # Aruldhas 2021 Table 2 row labelled 'V3 (L) = 335 (23.6)' with 62.23% BSV for R-methadone; reassigned to V4 (parent peripheral) per the paper's Figure 1 legend and Results 'Structural model of R methadone' paragraph 1 (V3 = V2 * VF cannot carry an independent typical value or BSV when VF = 1 fixed; see implementation notes above)
    lfdepot          <- log(0.718);          label("Oral bioavailability F (fraction; common to tablet and suspension formulations)")             # Aruldhas 2021 Table 2 'F = 0.718 (13.2)' for R-methadone

    # ----- R-methadone -> R-EDDP fractional metabolite-formation clearance CLF (dimensionless) -----
    lclf             <- log(0.217);          label("Typical fraction of R-methadone total apparent clearance CL that becomes R-EDDP metabolite (dimensionless) at CYP2B6 activity score 1 and heterozygous rs2246709")  # Aruldhas 2021 Table 2 'CLF = 0.217 (31.9)' for R-methadone; dimensionless fraction of CL per Figure 1 caption and Table 2 note (see implementation notes above)

    # ----- Metabolite (R-EDDP) disposition (Aruldhas 2021 Table 2) -----
    lcl_eddp         <- log(25.7);           label("Typical R-EDDP apparent clearance CL3 (L/h)")                                                  # Aruldhas 2021 Table 2 'CL3 (L*h^-1) = 25.7 (37.3)' for R-methadone
    lvf              <- fixed(log(1));       label("VF scaling factor between R-methadone V2 and R-EDDP V3 (V3 = V2 * VF; dimensionless, 1 due to metabolite unidentifiability)")  # Aruldhas 2021 Table 2 'VF = 1 FIX' and Figure 1 legend

    # ----- Allometric exponent (fixed) -----
    e_wt_cl                      <- fixed(0.75);  label("Allometric WT exponent on R-methadone CL (unitless)")                                    # Aruldhas 2021 Results 'Covariate modeling on R methadone' paragraph 1: "Bodyweight was added allometrically to the clearance of the drug with a fixed exponent of 0.75"; estimating the exponent gave 1.2 but was fixed to 0.75 to avoid small-cohort bias

    # ----- Covariate effects on V2 (linear centered-additive multiplicative form) -----
    e_snp_orm1_rs17650_vc        <- -0.443;       label("Linear covariate coefficient of ORM1 rs17650 active-allele count on R-methadone V2 (unitless per active allele)")  # Aruldhas 2021 Table 2 'V2rs17650 = -0.443 (35.4)' for R-methadone; V2 = V2_typ * (1 + coef * (n_active - 1))
    e_aag_vc                     <- -0.00291;     label("Linear covariate coefficient of AAG concentration (centered at 94.76 ug/mL) on R-methadone V2 (unitless per ug/mL)")  # Aruldhas 2021 Table 2 'V2AAG = -0.00291 (31.3)' for R-methadone; V2 = V2_typ * (1 - 0.00291 * (AAG - 94.76))

    # ----- Covariate effects on CLF (linear centered-additive multiplicative form) -----
    e_cyp2b6_clf                 <- 0.745;        label("Linear covariate coefficient of CYP2B6 activity score on R-methadone CLF (unitless per activity-score unit)")  # Aruldhas 2021 Table 2 'CLFCYP2B6 = 0.745 (17.4)' for R-methadone; CLF = CLF_typ * (1 + 0.745 * (CYP2B6 - 1))
    e_snp_cyp3a4_rs2246709_clf   <- 0.450;        label("Linear covariate coefficient of CYP3A4 rs2246709 active-allele count on R-methadone CLF (unitless per active allele)")  # Aruldhas 2021 Table 2 'CLFrs2246709 = 0.450 (33.9)' for R-methadone; CLF = CLF_typ * (1 + 0.450 * (rs2246709 - 1))

    # ----- Between-subject variability (Aruldhas 2021 Table 2) -----
    # Exact log-normal conversion omega^2 = log(CV^2 + 1).
    etalcl           ~ 0.42902   # Aruldhas 2021 Table 2 'CL BSV %CV = 72.1' for R-methadone; log(1 + 0.721^2) = 0.42902
    etalvc           ~ 0.48178   # Aruldhas 2021 Table 2 'V2 BSV %CV = 79.4' for R-methadone; log(1 + 0.794^2) = 0.48178
    etalvp           ~ 0.31735   # Aruldhas 2021 Table 2 'V4 BSV %CV = 62.23' (labelled 'V3' in Table 2, see implementation notes) for R-methadone; log(1 + 0.6223^2) = 0.31735
    etalclf          ~ 0.34353   # Aruldhas 2021 Table 2 'CLF BSV %CV = 65.0' for R-methadone; log(1 + 0.65^2) = 0.34353
    etalcl_eddp      ~ 0.22095   # Aruldhas 2021 Table 2 'CL3 BSV %CV = 49.8' for R-methadone; log(1 + 0.498^2) = 0.22095

    # ----- Residual error (Aruldhas 2021 Table 2) -----
    propSd           <- 0.165;   label("Proportional residual SD for R-methadone plasma concentration (fraction)")  # Aruldhas 2021 Table 2 'RUV drug = 0.165 (5.17)' for R-methadone
    propSd_eddp      <- 0.207;   label("Proportional residual SD for R-EDDP plasma concentration (fraction)")       # Aruldhas 2021 Table 2 'RUV metabolite = 0.207 (5.38)' for R-methadone
  })
  model({
    # ----- Reference values (Aruldhas 2021 covariate-model equations) -----
    ref_wt  <- 70      # kg, standard allometric reference
    ref_aag <- 94.76   # ug/mL, cohort average per Aruldhas 2021 Results 'Covariate modeling on R methadone' paragraph 3

    # ----- Individual parameters -----
    # Absorption: single depot Ka selected per oral dose by the FORM_TABLET
    # covariate. Suspension is the canonical reference (FORM_TABLET = 0),
    # tablet is a log-additive multiplicative shift.
    ka      <- exp(lka + e_form_tablet_ka * FORM_TABLET)

    # Parent CL: allometric WT effect (fixed exponent 0.75) on the typical
    # value; log-normal BSV.
    cl      <- exp(lcl + etalcl) * (WT / ref_wt)^e_wt_cl

    # Parent central volume V2: linear centered-additive effects of the ORM1
    # rs17650 active-allele count and the (time-varying) AAG concentration.
    vc      <- exp(lvc + etalvc) *
               (1 + e_snp_orm1_rs17650_vc * (SNP_ORM1_RS17650 - 1)) *
               (1 + e_aag_vc              * (AAG - ref_aag))

    # Intercompartmental clearance Q (no covariates, no reported BSV).
    q       <- exp(lq)

    # Parent peripheral volume V4 (paper's "V3" row in Table 2; see notes).
    vp      <- exp(lvp + etalvp)

    # Oral bioavailability (common to both oral formulations).
    fdepot  <- exp(lfdepot)

    # Fractional metabolite-formation clearance CLF: linear centered-additive
    # effects of the CYP2B6 activity score and the CYP3A4 rs2246709 active-
    # allele count; log-normal BSV.
    clf     <- exp(lclf + etalclf) *
               (1 + e_cyp2b6_clf               * (CYP2B6 - 1)) *
               (1 + e_snp_cyp3a4_rs2246709_clf * (SNP_CYP3A4_RS2246709 - 1))

    # Metabolite clearance CL3 (no covariates, log-normal BSV).
    cl_eddp <- exp(lcl_eddp + etalcl_eddp)

    # Metabolite central volume V3 = V2 * VF (VF = 1 fixed).
    vf      <- exp(lvf)
    vc_eddp <- vc * vf

    # Parent-to-metabolite formation clearance (mass / time), derived from
    # the total parent CL and the fractional-formation clearance CLF.
    cl_met_form <- cl * clf

    # ----- ODE system -----
    # First-order oral absorption into the parent central compartment.
    d/dt(depot)        <- -ka * depot

    # Parent two-compartment disposition with total-clearance elimination.
    d/dt(central)      <-  ka * depot -
                           cl * (central / vc) -
                           q  * (central / vc - peripheral1 / vp)
    d/dt(peripheral1)  <-  q  * (central / vc - peripheral1 / vp)

    # Metabolite one-compartment disposition. Only the CLF-partitioned
    # parent-CL flux enters the metabolite compartment; the metabolite is
    # eliminated with its own apparent clearance CL3.
    d/dt(central_eddp) <-  cl_met_form * (central / vc) -
                           cl_eddp     * (central_eddp / vc_eddp)

    # Oral bioavailability applied to the depot compartment (IV doses target
    # `central` directly and are not affected).
    f(depot) <- fdepot

    # ----- Observations -----
    # Dose in mg / volumes in L => concentrations in mg/L (= ug/mL). Paper
    # assay units are ng/mL; 1 mg/L = 1000 ng/mL. Vignette converts to
    # ng/mL for comparison against published narrative and figures.
    Cc       <- central      / vc
    Cc_eddp  <- central_eddp / vc_eddp

    Cc       ~ prop(propSd)
    Cc_eddp  ~ prop(propSd_eddp)
  })
}
