Sloan_2017_rifampicin <- function() {
  description <- "One-compartment population PK model for oral rifampin in Malawian adults with smear-positive pulmonary tuberculosis (Sloan 2017), developed using a two-stage NONMEM workflow: stage 1 fit a one-compartment + Savic 2007 transit-compartment absorption chain (NN, MTT, Ka) to 47 intensively-sampled patients, then stage 2 fit CL/F and V/F (plus IIVs and a multiplicative sex effect on CL) to 174 sparsely-sampled patients with absorption parameters fixed at the stage 1 estimates; F is fixed at 1, between-subject variability is on CL/F, V/F, and (fixed from stage 1) MTT, and an allometric weight model with fixed exponents 0.75 / 1.0 is referenced to 70 kg."
  reference <- paste(
    "Sloan DJ, McCallum AD, Schipani A, Egan D, Mwandumba HC, Ward SA,",
    "Waterhouse D, Banda G, Allain TJ, Owen A, Khoo SH, Davies GR. (2017).",
    "Genetic determinants of the pharmacokinetic variability of rifampin",
    "in Malawian adults with pulmonary tuberculosis.",
    "Antimicrob Agents Chemother 61(7):e00210-17.",
    "doi:10.1128/AAC.00210-17.",
    sep = " "
  )
  vignette <- "Sloan_2017_rifampicin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Drives allometric scaling of CL/F (fixed exponent 0.75) and V/F (fixed exponent 1.0), both standardized to a 70-kg patient (Sloan 2017 Methods 'Population pharmacokinetic analysis' paragraph 4: 'An allometric weight model was applied to standardize the pharmacokinetic parameters using a standard weight (wt_std) of 70 kg'). Sparse cohort: median 52 kg, range 34-74 kg (Table 1).",
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female) -- the typical value CL/F = 19.6 L/h applies to females; males have a multiplicative effect of 1.2x (95% CI 1.0-1.3) on CL/F per Sloan 2017 Table 3.",
      notes              = "Sloan 2017 Table 3 reports THETA_sex_male = 1.2 with the published narrative (Results paragraph 4) describing 'an increase of clearance of 17% in male patients'. The table value 1.2 is used canonically; the slight prose/table discrepancy is rounding (1.17 vs 1.2). Storage uses canonical SEXF (1 = female, 0 = male); the model applies the effect as 1.2^(1 - SEXF) so females (SEXF = 1) give factor 1 and males (SEXF = 0) give factor 1.2.",
      source_name        = "SEX"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age in years.",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a continuous covariate in the stage 2 model (Sloan 2017 Methods 'Population pharmacokinetic analysis' paragraph 3: 'The following covariates were explored: body weight, age, gender, HIV status, and SNP genotypes'). Did not significantly decrease the OFV and was not retained in the final model."
    ),
    HIV = list(
      description = "HIV infection status indicator (1 = HIV-infected, 0 = HIV-uninfected).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened in the stage 2 covariate search; not retained in the final model. 56.3% of the sparse cohort (98 of 174) were HIV-infected, of whom 28 (28.6%) were on antiretroviral therapy at recruitment."
    ),
    SLCO1B1_rs11045819 = list(
      description = "SLCO1B1 rs11045819 SNP genotype (CC = wild, AC = heterozygous, AA = variant); SLCO1B1 encodes the OATP1B1 hepatic uptake transporter.",
      units       = "(genotype)",
      type        = "categorical",
      notes       = "Tested as a covariate on CL/F and on relative bioavailability F (Sloan 2017 Results paragraph 7: 'Additive, dominant, and recessive models of effect were tested for each SNP'). Did not significantly improve model fit; variant allele frequency in the Malawian cohort was 0.07 (Sloan 2017 Discussion paragraph 4)."
    ),
    SLCO1B1_rs4149032 = list(
      description = "SLCO1B1 rs4149032 SNP genotype (TT = wild, CT = heterozygous, CC = variant).",
      units       = "(genotype)",
      type        = "categorical",
      notes       = "Tested on CL/F and F; not retained. Minor allele frequency 0.32 in the Malawian cohort; previously associated with reduced rifampin exposure in South African and Ugandan cohorts (Sloan 2017 Discussion paragraph 3)."
    ),
    AADAC_rs1803155 = list(
      description = "AADAC rs1803155 SNP genotype (GG = wild, CG = heterozygous, CC = variant); AADAC encodes a serine deacetylase responsible for 25-deacetylation of rifamycins.",
      units       = "(genotype)",
      type        = "categorical",
      notes       = "Tested on CL/F and F; not retained. Minor allele frequency 0.25 in the Malawian cohort."
    ),
    AADAC_rs61733692 = list(
      description = "AADAC rs61733692 SNP genotype.",
      units       = "(genotype)",
      type        = "categorical",
      notes       = "All 174 participants were wild-type (TT) for this SNP; not testable as a covariate."
    ),
    CES1_rs12149368 = list(
      description = "CES-1 rs12149368 SNP genotype; CES-1 encodes carboxylesterase-1, an alternative hepatic ester-bond-cleaving enzyme.",
      units       = "(genotype)",
      type        = "categorical",
      notes       = "Exploratory SNP screened as a covariate; only one of 174 participants was heterozygous, rest wild-type; not retained."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 174L,
    n_studies      = 1L,
    age_range      = "17-61 years (median 30)",
    age_median     = "30 years",
    weight_range   = "34-74 kg (median 52)",
    weight_median  = "52 kg",
    sex_female_pct = 30.5,
    race_ethnicity = "Black African (all participants were black Africans newly diagnosed with TB at Queen Elizabeth Central Hospital, Blantyre, Malawi).",
    disease_state  = "Smear-positive pulmonary tuberculosis (TB); HIV co-infection in 98 of 174 (56.3%), of whom 28 (28.6%) were on antiretroviral therapy at recruitment (median CD4 174 cells/uL, range 6-783).",
    dose_range     = "Oral rifampin / isoniazid fixed-dose combination tablets per WHO weight-adjusted bands: 300/150 mg (n=2, 1.1%), 450/225 mg (n=113, 64.9%), or 600/300 mg (n=59, 33.9%) once daily; corresponding to 8-12 mg/kg rifampin.",
    regions        = "Malawi (Blantyre; Queen Elizabeth Central Hospital, sparse cohort recruited 2010-2012).",
    co_medication  = "Standard first-line antitubercular regimen (rifampin + isoniazid +/- pyrazinamide + ethambutol per WHO guidelines). Antiretroviral therapy in 28 of 174 (16.1%).",
    notes          = "Sparse PK sampling at 0 (predose), 2, and 6 h post-dose on day 14 or 21 of TB treatment, drawn after an overnight fast at 7:30 a.m. Drug assay LLQ 0.5 ug/mL (LC-MS/MS); ~5% of samples below LLQ were imputed at LLQ/2. The stage 1 dataset comprised 47 intensively-sampled patients from a prior 2007-2008 cohort at the same hospital (Sloan 2017 Methods 'The intensively sampled pharmacokinetic data'); their absorption parameters Ka, NN, MMT, and IIV_MMT were carried forward as fixed values in the final stage 2 model. Pharmacogenetic SNPs in SLCO1B1, AADAC, and CES-1 were screened on CL/F and F but none significantly improved the model fit (Sloan 2017 Results paragraph 7); they are documented in covariatesDataExcluded."
  )

  ini({
    # ============================================================================
    # Structural PK parameters - final stage 2 typical values (Sloan 2017 Table 3
    # 'Parameter value estimates for the final model stage 2'). The typical-value
    # CL/F = 19.6 L/h is the female-reference value at WT = 70 kg; males get a
    # multiplicative factor of 1.2 (e_sex_cl below) and weight scales allometrically.
    # ============================================================================
    lcl <- log(19.6); label("Apparent oral clearance CL/F at WT = 70 kg, female reference (L/h)")        # Sloan 2017 Table 3: CL/F = 19.6 L/h (RSE 12%, 95% CI 16.7-22.5)
    lvc <- log(23.6); label("Apparent central volume of distribution V/F at WT = 70 kg (L)")            # Sloan 2017 Table 3: V/F = 23.6 L (RSE 9%, 95% CI 17.1-30.1)

    # ============================================================================
    # Absorption parameters - FIXED in the stage 2 model from the stage 1
    # estimates (Sloan 2017 Table 2 = Table 3 'FIX' rows). The two-stage workflow
    # used the 47 intensively-sampled patients in stage 1 to characterize the
    # absorption phase, then carried those values forward into the sparse-sampling
    # stage 2 because the sparse data (3 timepoints) carry little absorption
    # information (Sloan 2017 Methods 'Population pharmacokinetic analysis'
    # paragraph 1).
    # ============================================================================
    lka     <- fixed(log(0.277)); label("First-order absorption rate constant Ka (1/h)")                                                  # Sloan 2017 Table 3: Ka = 0.277/h (FIX); stage 1 estimate RSE 9.8%, 95% CI 0.23-0.32
    lmtt    <- fixed(log(0.326)); label("Mean transit time MMT (h)")                                                                      # Sloan 2017 Table 3: MMT = 0.326 h (FIX); stage 1 estimate RSE 35%, 95% CI 0.05-1.0
    lnn     <- fixed(log(1.5));   label("Number of absorption transit compartments NN (continuous, dimensionless)")                       # Sloan 2017 Table 3: NN = 1.5 (FIX); stage 1 estimate RSE 41%, 95% CI 0.09-3.9
    lfdepot <- fixed(log(1));     label("Oral bioavailability F (fixed at 1; CL/F and V/F are F-relative)")                               # Sloan 2017 Methods 'Population pharmacokinetic analysis' paragraph 4: 'F is the oral bioavailability, which was fixed to 1'

    # ============================================================================
    # Allometric weight model - canonical fixed exponents at reference weight 70 kg
    # (Sloan 2017 Methods 'Population pharmacokinetic analysis' paragraph 4:
    # 'An allometric weight model for clearance parameters is given by
    # CL/F_wt = (wt/wt_std)^(3/4) and for volume parameters is given by
    # V/F_wt = (wt/wt_std)^1').
    # ============================================================================
    allo_cl <- fixed(0.75); label("Allometric exponent on CL/F (unitless)")  # Sloan 2017 Methods: CL exponent fixed at 3/4
    allo_vc <- fixed(1.0);  label("Allometric exponent on V/F (unitless)")   # Sloan 2017 Methods: V exponent fixed at 1

    # ============================================================================
    # Sex effect on CL/F. Sloan 2017 Table 3 reports the multiplicative ratio
    # THETA_sex_male = 1.2 (95% CI 1.0-1.3): male CL = 1.2 x female CL. The
    # dichotomous covariate was 'introduced as a power model' (Methods
    # 'Population pharmacokinetic analysis' paragraph 6, equation 2), so it
    # enters as CL_i = CL_ref * theta_cov^X where X is the male indicator.
    # Converted to the canonical exponential-coefficient form e_sex_cl =
    # log(1.2) = 0.1823 and applied as (1 - SEXF) to preserve the male = 1
    # encoding while storing under canonical SEXF (1 = female).
    # ============================================================================
    e_sex_cl <- log(1.2); label("Exponential coefficient of male sex on CL (unitless; applied as (1 - SEXF))")  # Sloan 2017 Table 3: THETA_sex_male = 1.2 (RSE 13%, 95% CI 1.0-1.3); log(1.2) = 0.1823

    # ============================================================================
    # Inter-individual variability. Sloan 2017 reports IIVs as the variance
    # omega^2 directly (Methods 'Population pharmacokinetic analysis' paragraph 5:
    # 'eta_xi is the log interindividual variability for parameter x drawn from
    # a normal distribution with a mean of zero and variance omega^2'). IIV_MMT
    # is fixed from stage 1 because the sparse data carry little information
    # about MMT variability.
    # ============================================================================
    etalcl  ~ 0.076          # Sloan 2017 Table 3: IIVCL  = 0.076 (RSE 29%, 95% CI 0.033-0.11; shrinkage 22%); variance on log-CL
    etalvc  ~ 0.397          # Sloan 2017 Table 3: IIVV   = 0.397 (RSE 29%, 95% CI 0.17-0.63; shrinkage 28%); variance on log-V
    etalmtt ~ fixed(0.0706)  # Sloan 2017 Table 3: IIVMMT = 0.0706 (FIX from stage 1, RSE 75%, 95% CI 0.01-1.79); variance on log-MMT

    # ============================================================================
    # Residual error - proportional only (no additive or combined model retained;
    # Sloan 2017 Methods 'Population pharmacokinetic analysis' paragraph 2:
    # 'Proportional, additive, and combined proportional and additive error
    # models were considered'). Table 3 reports 'Proportional error (%) 0.22',
    # which is the SD on the linear concentration scale (~22% CV).
    # ============================================================================
    propSd <- 0.22; label("Proportional residual error (fraction)")  # Sloan 2017 Table 3: proportional error = 0.22 (RSE 12%, 95% CI 0.19-0.26)
  })

  model({
    # 1. Individual PK parameters with allometric weight scaling and sex effect on CL.
    #    Female (SEXF = 1) is the reference: CL_female = exp(lcl + etalcl) at 70 kg.
    #    Male (SEXF = 0): CL_male = CL_female * exp(e_sex_cl) = CL_female * 1.2.
    cl <- exp(lcl + etalcl + e_sex_cl * (1 - SEXF)) * (WT / 70) ^ allo_cl
    vc <- exp(lvc + etalvc)                          * (WT / 70) ^ allo_vc

    # 2. Absorption parameters (all fixed from stage 1).
    ka     <- exp(lka)
    mtt    <- exp(lmtt + etalmtt)
    nn     <- exp(lnn)
    fdepot <- exp(lfdepot)

    # 3. Elimination microconstant.
    kel <- cl / vc

    # 4. PK ODE system. Savic 2007 transit-compartment absorption (rxode2
    #    transit() implements the analytical gamma-PDF input flux for
    #    non-integer NN) feeds a virtual depot at rate ktr = (NN + 1) / MTT;
    #    the depot then absorbs into central via first-order ka. f(depot) <- 0
    #    suppresses the bolus arrival at depot; the transit() function reads
    #    the raw dose amount from podo(depot) regardless of f(depot) so the
    #    full dose enters the absorption process through the analytical chain.
    #    Bioavailability fdepot enters via the bio argument of transit().
    d/dt(depot)   <- transit(nn, mtt, fdepot) - ka * depot
    d/dt(central) <- ka * depot - kel * central

    f(depot) <- 0  # suppress bolus into depot; transit() drives the input rate

    # 5. Plasma rifampin concentration. Dose mg / volume L = mg/L (= ug/mL).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
