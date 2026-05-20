Roberts_2016_topotecan <- function() {
  description <- "One-compartment population pharmacokinetic model for oral topotecan lactone in infants and very young children with primary central nervous system tumours (Roberts 2016). First-order absorption into a depot compartment is followed by first-order elimination from a central compartment. Apparent volume of distribution (V/F) and apparent clearance (CL/F) are scaled by body surface area as power functions centred on the cohort median (0.57 m^2); the ABCG2 rs4148157 G>A variant (heterozygous AG or homozygous AA carriers pooled vs the GG reference) carries an exponential covariate effect on the absorption rate constant Ka, yielding an approximately 2-fold higher Ka in carriers than in GG homozygotes."
  reference <- "Roberts JK, Birg AV, Lin T, Daryani VM, Panetta JC, Broniscer A, Robinson GW, Gajjar AJ, Stewart CF. Population pharmacokinetics of oral topotecan in infants and very young children with brain tumors demonstrates a role of ABCG2 rs4148157 on the absorption rate constant. Drug Metab Dispos. 2016;44(7):1116-1122. doi:10.1124/dmd.115.068676"
  vignette <- "Roberts_2016_topotecan"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    BSA = list(
      description        = "Body surface area at the time of the pharmacokinetic study.",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effects on V/F (exponent 0.78, theta_2) and CL/F (exponent 1.25, theta_3) centred on the cohort median BSA = 0.57 m^2 (Roberts 2016 Table 1; equations in the Results 'Population Pharmacokinetic Analysis' section). The paper does not state which BSA computation formula was used; the validation vignette therefore uses Mosteller (sqrt(height_cm * weight_kg / 3600)) as the default. The same BSA value is used both to compute the body-surface-area-normalised oral topotecan dose (0.8 mg/m^2 nominal) and as the covariate in the model.",
      source_name        = "BSA"
    ),
    SNP_ABCG2_RS4148157 = list(
      description        = "Binary indicator for the ABCG2 rs4148157 G>A intronic variant (1 = heterozygous AG or homozygous AA carrier, 0 = homozygous wild-type GG).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (homozygous wild-type GG)",
      notes              = "Time-fixed per subject (germline genotype). Heterozygous (AG) and homozygous (AA) variant carriers were pooled in Roberts 2016 because only one AA homozygote was present in the cohort (Table 2). The pooled indicator is encoded as 0 (GG) or 1 (AG/AA) per the GENECAT term in the Roberts 2016 individual-Ka equation, and enters as an exponential factor on Ka: Ka_i = Ka_pop * exp(theta_1 * GENECAT) * exp(eta_Ka). The Roberts 2016 cohort distribution (Table 2) was GG = 42, AG = 9, AA = 1, missing = 9 (carrier rate 19% among the 52 genotyped patients).",
      source_name        = "GENECAT (rs4148157 AG/AA carrier indicator)"
    )
  )

  population <- list(
    species        = "human (paediatric)",
    n_subjects     = 61L,
    n_studies      = 1L,
    n_observations = 182L,
    age_range      = "0.48-4.59 years",
    age_median     = "2.37 years",
    weight_range   = "5.20-17.50 kg",
    weight_median  = "12.60 kg",
    height_range   = "59.5-102 cm",
    height_median  = "89.1 cm",
    bsa_range      = "0.31-0.72 m^2",
    bsa_median     = "0.57 m^2",
    sex_female_pct = 37.7,
    disease_state  = "Newly diagnosed primary central nervous system tumours treated on the SJYC07 protocol (Risk-Adapted Therapy for Infants and Young Children with Embryonal Brain Tumors, High Grade Glioma, Choroid Plexus Carcinoma or Ependymoma; NCT00602667). Diagnoses include embryonal tumours, high-grade glioma, choroid plexus carcinoma, and ependymoma. Patients had adequate organ function (Lansky performance status >= 30) and normal renal function at enrolment.",
    dose_range     = "Oral topotecan 0.8 mg/m^2 once daily for 10 days on a 28-day cycle (maintenance phase of SJYC07); co-administered with oral cyclophosphamide 30 mg/m^2 for 21 days of the same cycle (cyclophosphamide is not part of the topotecan PK model). Topotecan was administered as a liquid mixed in a flavoured vehicle.",
    regions        = "Single-centre cohort enrolled at St. Jude Children's Research Hospital, Memphis, Tennessee (the multicentre SJYC07 protocol is conducted across multiple sites but only St. Jude patients contributed to this oral topotecan PK analysis).",
    pharmacogenomics = "Germline genotyping for ABCG2 (rs4148157, rs2622628, rs2725252) and ABCB1 (rs1045642, rs2032582, rs1128503) variants by Illumina Infinium Omni2.5 Exome-8 BeadChip in 52 of 61 patients; missing genotype handled as a missing covariate within Monolix.",
    notes          = "Sampling per limited-sampling model (Turner 2006): pre-dose, 15 min, 90 min, and 6 h after the first observed dose. 39% of concentrations were below the limit of quantitation (1 ng/mL) and were treated as left-censored using the Monolix M3-equivalent likelihood. Bioanalytical: isocratic HPLC with fluorescence detection (370 nm excitation, 520 nm emission). Cohort demographics per Roberts 2016 Table 1; pharmacogenomic frequencies per Table 2; final population parameter estimates per Table 3."
  )

  ini({
    # Structural parameters - typical-value reference is a GG homozygous
    # rs4148157 patient at the cohort median BSA = 0.57 m^2 (Roberts 2016
    # Table 3 final population estimates).
    lka <- log(0.61); label("First-order oral absorption rate constant Ka at the rs4148157 GG reference (1/h)")    # Roberts 2016 Table 3: Ka(GG) = 0.61 1/h (SE 0.11)
    lvc <- log(40.2); label("Apparent volume of distribution V/F at the BSA = 0.57 m^2 reference (L)")            # Roberts 2016 Table 3: V/F = 40.2 L (SE 7.0)
    lcl <- log(40.0); label("Apparent clearance CL/F at the BSA = 0.57 m^2 reference (L/h)")                       # Roberts 2016 Table 3: CL/F = 40.0 L/h (SE 2.9)

    # Covariate effects. The rs4148157 effect on Ka is exponential
    # (Ka_i = Ka_pop * exp(theta_1 * GENECAT) with GENECAT in {0,1}); the
    # BSA effects on V/F and CL/F are power-form (P_i = P_pop * (BSA_i /
    # BSA_pop)^theta) centred on the cohort median 0.57 m^2. All three
    # coefficients are estimated, not fixed.
    e_abcg2_ka <- 1.06;  label("Exponential coefficient for ABCG2 rs4148157 variant on Ka (unitless); Ka_variant / Ka_GG = exp(1.06) = 2.89") # Roberts 2016 Table 3 theta_1 (SE 0.25)
    e_bsa_vc   <- 0.78;  label("Power exponent for BSA on V/F (unitless); V/F * (BSA/0.57)^0.78")                   # Roberts 2016 Table 3 theta_2 (SE 0.48)
    e_bsa_cl   <- 1.25;  label("Power exponent for BSA on CL/F (unitless); CL/F * (BSA/0.57)^1.25")                 # Roberts 2016 Table 3 theta_3 (SE 0.39)

    # IIV - exponential log-normal between-subject variability per
    # Roberts 2016 Methods (P_i = mu_pop * exp(eta_i), eta ~ N(0, omega^2)).
    # Monolix 4.3.2 reports omega as the standard deviation (square root of
    # the log-scale variance) by default; variances used by ini() below are
    # the squares of the table-3 reported omegas. No off-diagonal
    # covariances were reported.
    etalka ~ 0.2672    # Roberts 2016 Table 3 omega_Ka = 0.517 (SE 0.12); variance = 0.517^2 = 0.2672 (approx CV 56% on Ka)
    etalvc ~ 0.1142    # Roberts 2016 Table 3 omega_V  = 0.338 (SE 0.17); variance = 0.338^2 = 0.1142 (approx CV 35% on V/F)
    etalcl ~ 0.2162    # Roberts 2016 Table 3 omega_CL = 0.465 (SE 0.11); variance = 0.465^2 = 0.2162 (approx CV 49% on CL/F)

    # Residual error. Roberts 2016 Methods states that the final model
    # used an additive error on the linear-concentration scale; Table 3
    # reports the additive SD in ng/mL.
    addSd <- 0.592; label("Additive residual SD (ng/mL)")    # Roberts 2016 Table 3 sigma_add = 0.592 (SE 0.44)
  })

  model({
    # Reference covariate values. The BSA reference is the cohort median
    # (Roberts 2016 Table 1; median 0.57 m^2 across 61 patients) and is
    # the value the source paper centres the BSA power model on.
    ref_bsa <- 0.57

    # Individual parameters. The rs4148157 covariate enters as
    # exp(e_abcg2_ka * SNP_ABCG2_RS4148157), where the binary indicator is
    # 1 for AG/AA carriers and 0 for the GG reference; this matches the
    # Monolix individual-Ka equation in Roberts 2016 Results.
    ka <- exp(lka + etalka + e_abcg2_ka * SNP_ABCG2_RS4148157)
    vc <- exp(lvc + etalvc) * (BSA / ref_bsa)^e_bsa_vc
    cl <- exp(lcl + etalcl) * (BSA / ref_bsa)^e_bsa_cl

    # One-compartment first-order absorption / first-order elimination.
    kel <- cl / vc
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Topotecan dose is provided in mg (oral, body-surface-area normalised
    # at 0.8 mg/m^2); compartmental amounts are therefore in mg, and
    # central / vc has units of mg/L. The source paper reports plasma
    # topotecan lactone in ng/mL, so the observation is scaled by 1000
    # (1 mg/L = 1000 ng/mL = 1 ug/mL).
    Cc <- (central / vc) * 1000

    Cc ~ add(addSd)
  })
}
