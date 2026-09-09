Paiboonvong_2025_sitafloxacin <- function() {
  description <- paste0(
    "Population PK model for oral sitafloxacin in plasma and pulmonary epithelial lining fluid ",
    "(ELF) of 12 critically ill Thai patients with pneumonia, developed by Paiboonvong 2025 to ",
    "drive Monte Carlo probability-of-target-attainment simulation against Streptococcus ",
    "pneumoniae. Plasma disposition is one-compartment; absorption is a Savic transit chain of ",
    "one transit compartment in which the absorption rate constant is constrained to equal the ",
    "transit rate constant (ka = ktr = (ntr + 1) / mtt with ntr fixed to 1). An ELF compartment ",
    "hangs off the central compartment as a peripheral compartment with a fixed physiological ",
    "volume of 0.025 L and an inter-compartmental clearance q_elf. IMPORTANT: in the authors' ",
    "own control stream the central-to-ELF micro-rate constant carries the plasma unbound ",
    "fraction, k_central_elf = q_elf * fu * pcelf / vc while k_elf_central = q_elf / v_elf, so ",
    "the estimated partition coefficient pcelf = 0.772 is the ELF-to-UNBOUND-plasma ratio and ",
    "the ELF-to-total-plasma concentration (and AUC) ratio the model actually produces is ",
    "fu * pcelf = 0.63 * 0.772 = 0.486. Body weight is an a-priori allometric covariate on CL/F ",
    "(exponent 0.75 fixed), V/F (exponent 1 fixed) and Q/F_ELF (exponent 0.75 fixed), all ",
    "normalized to the cohort median 52 kg. Age is a LINEAR covariate on relative ",
    "bioavailability, F = 1 + 0.0258 * (AGE - 57), which is only usable over the adult age range ",
    "the authors simulated (30-70 years): the expression reaches zero at 18.2 years and is ",
    "negative below it. All disposition parameters are apparent (/F): only oral data were ",
    "analysed and F was fixed to 1 with its inter-individual variability estimated. Residual ",
    "variability was additive on the natural-log concentration scale for both plasma and ELF, ",
    "encoded here as lnorm()."
  )
  reference <- paste(
    "Paiboonvong T, Montakantikul P, Panjasawatwong N, Singkham N, Punyawudho B.",
    "Population pharmacokinetics and pharmacodynamics of sitafloxacin in plasma and alveolar",
    "epithelial lining fluid of critically ill Thai patients with pneumonia.",
    "Pharmacol Res Perspect. 2025;13(2):e70081. doi:10.1002/prp2.70081. PMCID: PMC11930543.",
    "Final parameter estimates: Table 2. Model structure, the transit-absorption chain with",
    "KA = KTR, the ELF compartment and the age-on-F relationship: Results section 3.3 and",
    "Figure 1. Structural encoding including the placement of the unbound fraction inside the",
    "central-to-ELF micro-rate constant and the omega / sigma values on their estimation scales:",
    "Data S1 (Supporting Information), NONMEM control streams Run56 (estimation) and Run59",
    "(simulation). Underlying clinical study, sampling design and bioanalysis:",
    "reference 7 of the paper.",
    "See also modellib('Wu_2025_sitafloxacin') and modellib('Rodjun_2023_sitafloxacin')",
    "for independent sitafloxacin popPK models in non-critically-ill populations.",
    sep = " "
  )
  vignette <- "Paiboonvong_2025_sitafloxacin"

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. The control stream numbers these 1 = gut, 3 = transit,
  # 2 = blood, 4 = ELF; the canonical names below carry the same roles.
  compartmentData <- list(
    depot    = list(analyte = "sitafloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "sitafloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    central  = list(analyte = "sitafloxacin", units = "mg", specimen = "plasma", verified = TRUE),
    elf      = list(analyte = "sitafloxacin", units = "mg", specimen = "epithelial lining fluid", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Added a priori as an allometric function on CL/F (exponent fixed to 0.75), V/F (exponent fixed to 1)",
        "and Q/F_ELF (exponent fixed to 0.75), normalized to the cohort MEDIAN weight of 52 kg rather than the",
        "usual 70 kg reference (Methods 2.2 'Body weight was added as an allometric function a priori'; Table 2",
        "reports CL/F and V/F per 52 kg explicitly, and Data S1 Run59 codes ((WT/52)**0.75) and (WT/52)).",
        "The allometric term on Q/F_ELF appears only in the control stream, which flags it in its own header",
        "comment ('Add allometric function on QELF'); Table 2 reports Q/F_ELF per 52 kg consistently with it.",
        "Observed range in the 12-patient cohort: median 52 kg, IQR 44-68 kg (Results 3.1)."
      ),
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate retained in the final model. Enters relative bioavailability LINEARLY and centered",
        "on the cohort median age of 57 years: F = 1 + 0.0258 * (AGE - 57), i.e. every 1-year increase in age",
        "raises F by 2.58 percentage points of the reference value (Results 3.3; Data S1 Run59",
        "F1AGE = (1 + THETA(8)*(AGE - 57))).",
        "CAUTION: because the form is linear rather than exponential or power, F is exactly zero at",
        "AGE = 57 - 1/0.0258 = 18.2 years and NEGATIVE below it. The authors simulated only ages 30, 40, 50, 60",
        "and 70 years (Methods 2.4) and the observed cohort spanned an IQR of 40-65 years, so the relationship",
        "must not be extrapolated to children or young adults. NONMEM bounded THETA(8) to (-0.056, 0.032).",
        "Observed range in the cohort: median 57 years, IQR 40-65 years (Results 3.1)."
      ),
      source_name        = "AGE"
    )
  )

  # Screened in the stepwise covariate procedure (Methods 2.2) but NOT retained
  # in the final model, so they are documented rather than carried in
  # covariateData. All three are present as columns in the analysis dataset
  # ($INPUT of Data S1 Run56).
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "categorical",
      notes       = paste(
        "Screened for an effect on the PK parameters by forward selection / backward elimination and not",
        "retained (Methods 2.2; Results 3.3 reports only age as significant). Six of the 12 patients were male",
        "(Results 3.1). The dataset column is SEX with 0 = male and 1 = female (Data S1 Run56 $INPUT),",
        "which is already the canonical SEXF polarity."
      ),
      source_name = "SEX"
    ),
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault equation (raw, NOT BSA-normalized)",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Screened and not retained. The authors flag this explicitly as a surprise, since creatinine clearance",
        "is a known covariate on sitafloxacin CL/F in non-critically-ill patients, and attribute the negative",
        "finding to the narrow renal-function range of the 12-patient cohort (66.7% mild-to-moderate impairment)",
        "and to a single augmented-renal-clearance outlier at 235 mL/min (Discussion).",
        "Median 68 mL/min, IQR 30-96 mL/min (Results 3.1)."
      ),
      source_name = "CLCR"
    ),
    APACHE_II = list(
      description = "Acute Physiology and Chronic Health Evaluation II score",
      units       = "(score)",
      type        = "continuous",
      notes       = paste(
        "Screened as a continuous illness-severity covariate and not retained (Methods 2.2).",
        "Median 21, IQR 18-33 (Results 3.1)."
      ),
      source_name = "APACHE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 12L,
    n_studies      = 1L,
    age_median     = "57 years (IQR 40-65)",
    weight_median  = "52 kg (IQR 44-68)",
    sex_female_pct = 50,
    disease_state  = paste(
      "Critically ill patients with pneumonia admitted to an intensive care unit.",
      "Median APACHE II score 21 (IQR 18-33), median serum albumin 2.0 g/dL (IQR 1.8-2.2),",
      "median creatinine clearance 68 mL/min (IQR 30-96); 66.7% had mild-to-moderate renal",
      "impairment and one patient had augmented renal clearance (CLcr 235 mL/min)."
    ),
    dose_range     = "Sitafloxacin 200 mg orally as a single dose under fasting conditions",
    regions        = "Thailand",
    notes          = paste(
      "Baseline demographics: Results 3.1. The model was built on 83 plasma and 12 ELF concentrations;",
      "one plasma concentration below the 0.025 mg/L LLOQ was excluded (Results 3.3).",
      "Plasma sampling at pre-dose, 0.5, 1, 2, 3, 8 and 12 h; one bronchoalveolar-lavage sample per",
      "patient, randomly assigned to the 0.5-2, 3-4, 5-6 or 7-9 h window (Methods 2.1). Concentrations",
      "were assayed by LC-MS/MS and analysed on the natural-log scale.",
      "The Monte Carlo simulations that generated the paper's target-attainment results used 5000 virtual",
      "adults at each of the ages 30, 40, 50, 60 and 70 years, dosed 50 or 100 mg every 12 h to steady",
      "state, with the PK/PD target fAUC/MIC > 30 and the MIC distribution of S. pneumoniae isolated from",
      "Thai patients (Methods 2.4, Table 1)."
    )
  )

  ini({
    # ---- Plasma disposition (Table 2; Data S1 Run59 $THETA) ---------------
    # Reference subject: 52 kg (cohort median body weight).
    lcl <- log(7.03); label("Log apparent oral clearance CL/F at 52 kg (L/h)")                                   # Table 2: CL/F = 7.03 L/h/52kg, RSE 24.9%, SIR 95% CI 4.36-11.0; Run59 THETA(1)
    lvc <- log(116); label("Log apparent central volume of distribution V/F at 52 kg (L)")                       # Table 2: V/F = 116 L/52kg, RSE 13.0%, SIR 95% CI 93.9-153; Run59 THETA(2)

    # ---- Transit absorption (Table 2; Data S1 Run59) ----------------------
    # One fixed transit compartment with the absorption rate constant
    # constrained to the transit rate constant, KA = KTR = (NN+1)/MTT with
    # NN = 1 (Results 3.3; Run59 "NN = 1", "KTR = (NN+1)/MTT", "K13 = KTR",
    # "K32 = KTR"). Only MTT is estimated.
    lmtt <- log(1.48); label("Log mean transit absorption time MTT (h)")                                         # Table 2: MTT = 1.48 h, RSE 27.1%, SIR 95% CI 0.858-2.44; Run59 THETA(4)

    # Relative bioavailability was FIXED to 1 and only its IIV estimated
    # (Methods 2.2: "Relative oral bioavailability (F) was fixed to 1, and its
    # inter-individual variability (IIV) was estimated"); Run59 THETA(3) is
    # "(1) FIX". The age covariate multiplies this reference value.
    lfdepot <- fixed(log(1)); label("Log relative oral bioavailability of the depot at the reference age (unitless)") # Table 2: F = 1 fix; Run59 THETA(3) (1) FIX

    # ---- ELF limb (Table 2; Data S1 Run59) --------------------------------
    lq_elf <- log(0.0441); label("Log apparent central-to-ELF inter-compartmental clearance Q/F_ELF at 52 kg (L/h)") # Table 2: Q/F_ELF = 0.0441 L/h/52kg, RSE 36.5%, SIR 95% CI 0.0194-0.0812; Run59 THETA(5)
    # V/F_ELF "could not be reliably estimated" and was fixed to a lung
    # physiological value taken from reference 11 of the paper (Results 3.3).
    lv_elf <- fixed(log(0.025)); label("Log apparent ELF volume of distribution V/F_ELF (L)")                    # Table 2: V/F_ELF = 0.025 L fix; Run59 THETA(6) (0.025) FIX
    lpcelf <- log(0.772); label("Log ELF-to-unbound-plasma partition coefficient (unitless)")                    # Table 2: PC = 0.772, RSE 18.2%, SIR 95% CI 0.556-1.08; Run59 THETA(7)
    # The unbound fraction is a fixed literature value that the authors put
    # INSIDE the central-to-ELF micro-rate constant, not only in their
    # post-hoc free-AUC calculation: Run59 codes "K24 = QELF*0.63*PC/VBD"
    # against "K42 = QELF/VELF". See the description and the vignette Errata.
    fu <- fixed(0.63); label("Fraction of sitafloxacin unbound in plasma (unitless)")                            # Table 2: fu = 0.63 fix; Methods 2.4 cites reference 7; Run59 K24 = QELF*0.63*PC/VBD

    # ---- Allometric exponents (a priori, all fixed) -----------------------
    e_wt_cl    <- fixed(0.75); label("Allometric exponent of body weight on CL/F (unitless)")                    # Methods 2.2: "allometric function a priori on CL/F with the exponent fixed to 0.75"; Run59 TVCLBD = THETA(1)*((WT/52)**0.75)
    e_wt_vc    <- fixed(1);    label("Allometric exponent of body weight on V/F (unitless)")                     # Methods 2.2: "volume of distribution (V/F) with the exponent fixed to 1"; Run59 TVVBD = THETA(2)*(WT/52)
    e_wt_q_elf <- fixed(0.75); label("Allometric exponent of body weight on Q/F_ELF (unitless)")                 # Run59 TVQELF = THETA(5)*((WT/52)**0.75) and header comment "Add allometric function on QELF"

    # ---- Covariate effect (Table 2; Data S1 Run59) ------------------------
    # LINEAR, centered on the median age of 57 years. See covariateData$AGE
    # for the extrapolation caution: F reaches zero at 18.2 years.
    e_age_fdepot <- 0.0258; label("Linear effect of age, centered at 57 years, on relative bioavailability (1/year)") # Table 2: AGE on F = 2.58%, RSE 7.64%, SIR 95% CI 2.09-2.84; Results 3.3 "F = 1 + (Age - 57) * 0.0258"; Run59 THETA(8)

    # ---- Inter-individual variability (Table 2; Data S1 Run59 $OMEGA) ----
    # Table 2 reports IIV as %CV; Run59 $OMEGA carries the log-scale
    # VARIANCES, which is what nlmixr2 wants. The two agree exactly under
    # CV = sqrt(exp(omega^2) - 1): sqrt(exp(0.566)-1) = 87.2%,
    # sqrt(exp(0.0285)-1) = 17.0%, sqrt(exp(0.626)-1) = 93.3%,
    # sqrt(exp(0.0528)-1) = 23.3%.
    #
    # No eta is carried on V/F, Q/F_ELF or V/F_ELF: the paper states the data
    # did not support estimating them (Results 3.3, Discussion limitations)
    # and Run59 fixes OMEGA(2,2), OMEGA(5,5) and OMEGA(6,6) to 0.
    etalcl    ~ 0.566   # Run59 $OMEGA 1; Table 2: IIV CL/F 87.2 %CV, RSE 17.8%, SIR 95% CI 45.9-125
    etalfdepot ~ 0.0285 # Run59 $OMEGA 3; Table 2: IIV F 17.0 %CV, RSE 17.1%, SIR 95% CI 9.87-22.1
    etalmtt   ~ 0.626   # Run59 $OMEGA 4; Table 2: IIV MTT 93.3 %CV, RSE 27.1%, SIR 95% CI 64.7-177
    etalpcelf ~ 0.0528  # Run59 $OMEGA 7; Table 2: IIV PC 23.3 %CV, RSE 18.4%, SIR 95% CI 17.2-33.0

    # ---- Residual variability (Table 2; Data S1 Run59 $SIGMA) ------------
    # "An additive residual error model on the logarithmic scale was used to
    # describe the residual variability (RUV) for both plasma and ELF"
    # (Results 3.3), i.e. Y = LOG(IPRED) + EPS, which is a log-normal
    # (exponential) residual on the linear concentration scale -> lnorm().
    # Table 2's sigma rows are that log-scale variance re-expressed as %CV:
    # sqrt(exp(0.122)-1) = 36.0% and sqrt(exp(0.249)-1) = 53.2%.
    expSd      <- sqrt(0.122); label("Log-scale residual SD for plasma concentrations (natural-log units)")      # Run59 $SIGMA 1 = 0.122; Table 2: sigma plasma 36.0, RSE 10.5%, SIR 95% CI 29.4-44.6
    expSd_Celf <- sqrt(0.249); label("Log-scale residual SD for ELF concentrations (natural-log units)")         # Run59 $SIGMA 2 = 0.249; Table 2: sigma ELF 53.2, RSE 14.4%, SIR 95% CI 36.3-70.8
  })

  model({
    # 1. Individual parameters. Body weight is allometric on CL/F, V/F and
    #    Q/F_ELF, normalized to the cohort median of 52 kg.
    cl    <- exp(lcl + etalcl) * (WT / 52)^e_wt_cl
    vc    <- exp(lvc) * (WT / 52)^e_wt_vc
    q_elf <- exp(lq_elf) * (WT / 52)^e_wt_q_elf
    v_elf <- exp(lv_elf)
    mtt   <- exp(lmtt + etalmtt)
    pcelf <- exp(lpcelf + etalpcelf)

    # Relative bioavailability: fixed reference value of 1, multiplied by the
    # linear age effect centered at 57 years, with exponential IIV.
    fdepot <- exp(lfdepot + etalfdepot) * (1 + e_age_fdepot * (AGE - 57))

    # 2. Micro-constants. One transit compartment, with the absorption rate
    #    constant constrained to equal the transit rate constant.
    ntr <- 1
    ktr <- (ntr + 1) / mtt
    kel <- cl / vc
    # The unbound fraction sits in the central-to-ELF direction only, exactly
    # as the authors coded it, so ELF equilibrates to fu * pcelf times the
    # TOTAL plasma concentration.
    k_central_elf <- q_elf * fu * pcelf / vc
    k_elf_central <- q_elf / v_elf

    # 3. ODE system (Data S1 Run59 $DES, states renamed to canonical roles:
    #    A(1) = depot, A(3) = transit1, A(2) = central, A(4) = elf).
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(central)  <- ktr * transit1 - kel * central -
      k_central_elf * central + k_elf_central * elf
    d/dt(elf)      <- k_central_elf * central - k_elf_central * elf

    # 4. Bioavailability
    f(depot) <- fdepot

    # 5. Observations
    Cc   <- central / vc
    Celf <- elf / v_elf

    Cc   ~ lnorm(expSd)
    Celf ~ lnorm(expSd_Celf)
  })
}
