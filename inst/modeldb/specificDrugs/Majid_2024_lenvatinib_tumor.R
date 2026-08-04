Majid_2024_lenvatinib_tumor <- function() {
  description <- "Integrated tumor-growth-inhibition and serum-biomarker PK/PD model for lenvatinib in radioiodine-refractory differentiated thyroid cancer (Majid 2024). First-order tumor growth (rate KG) is opposed by three shrinkage terms acting on the RECIST 1.1 sum of longest diameters: an Emax function of lenvatinib average daily AUC whose potency decays exponentially with an acquired-resistance rate lambda, and two linear terms proportional to the model-predicted relative change from baseline of Tie-2 and Ang-2. The Tie-2 and Ang-2 time courses are re-integrated inside this model from the companion four-biomarker turnover model with their parameters fixed. Lenvatinib exposure enters as the AUC_LEN covariate; there is no PK ODE. Time is in weeks. Companion models: modellib('Majid_2024_lenvatinib') and modellib('Majid_2024_lenvatinib_biomarkers')."
  reference <- paste(
    "Majid O, Hayato S, Sreerama Reddy SH, Lalovic B, Hihara T, Hoshi T,",
    "Funahashi Y, Aluri J, Takenaka O, Yasuda S, Hussein Z.",
    "Population pharmacokinetic-pharmacodynamic modeling of serum biomarkers as",
    "predictors of tumor dynamics following lenvatinib treatment in patients with",
    "radioiodine-refractory differentiated thyroid cancer (RR-DTC).",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(6):954-969.",
    "doi:10.1002/psp4.13130.",
    "Equation 4 and Table 3; NONMEM control stream in Supplementary Text S3.",
    "Tumor-dynamics form follows Claret et al. J Clin Oncol. 2009;27(25):4103-4108.",
    sep = " "
  )
  vignette <- "Majid_2024_lenvatinib"
  paper_specific_compartments <- c("tie2", "ang2")
  paper_specific_etas <- c("etaibase_prop", "etaibase_add")
  units <- list(
    time          = "week",
    dosing        = "n/a (no drug-dosing events; lenvatinib exposure enters as the AUC_LEN covariate, not via a PK ODE)",
    concentration = "mm (RECIST 1.1 sum of longest diameters of target lesions; not a drug concentration)"
  )

  covariateData <- list(
    AUC_LEN = list(
      description        = "Lenvatinib average steady-state daily AUC over the interval between two tumor assessments, driving the Emax tumor-shrinkage term and the Tie-2 / Ang-2 turnover sub-models.",
      units              = "ng*h/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying, step-wise between tumor assessments. Majid 2024 Methods: 'Exposure is lenvatinib AUC calculated based on the average dose between two tumor assessments.' In the source sequential fit it is the NONMEM data column LEVAVAUC, derived from the upstream population PK run (Text S1: AUC = 1000 * F1 * DGRP / CL). Text S3 rescales it to ug*h/mL with 'AUC = LEVAVAUC/1000' and correspondingly divides the carried EC50 values by 1000; this model keeps everything in ng*h/mL so AUC_LEN is interchangeable with the value used by Majid_2024_lenvatinib_biomarkers.R. Set to 0 for placebo subjects and off-treatment periods, which makes the drug-shrinkage term vanish exactly.",
      source_name        = "LEVAVAUC"
    ),
    TUM_SLD = list(
      description        = "Observed baseline RECIST 1.1 sum of longest diameters of target lesions at study entry; per-subject, time-fixed. Used as the deterministic component of the tumor ODE initial condition.",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Text S3 sets BAS = BTM (the baseline tumor measurement data column) at TIME = 0 and forms the individual initial condition as IBASE = BAS + ETA(5)*THETA(5)*BAS + ETA(6)*THETA(6) with OMEGA(5,5) = OMEGA(6,6) = 1 FIX -- i.e. the observed baseline perturbed by the model's own combined residual error, the IPP-style baseline-residual construction also used in Hansson_2013b_sunitinib.R. Tumor-size PK/PD population baseline: mean 70.2 mm, median 59.5 mm, range 10.2-331.0 mm (Table S2).",
      source_name        = "BTM"
    )
  )

  compartmentData <- list(
    tie2 = list(
      analyte  = "soluble TEK tyrosine kinase 2 (Tie-2)",
      units    = "ng/mL",
      specimen = "serum",
      verified = TRUE
    ),
    ang2 = list(
      analyte  = "angiopoietin 2 (Ang-2)",
      units    = "ng/mL",
      specimen = "serum",
      verified = TRUE
    ),
    tumor = list(
      analyte  = "tumor burden (RECIST 1.1 sum of longest diameters of target lesions)",
      units    = "mm",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 558L,
    n_studies      = 2L,
    age_range      = "not reported as a range; 299 of 558 subjects were < 65 years and 258 were >= 65 years (Table S2)",
    weight_range   = "31.0-190 kg",
    weight_median  = "75 kg",
    sex_female_pct = 48.6,
    race_ethnicity = c(
      White                             = 74.0,
      Asian_other_than_Chinese_Japanese = 9.3,
      Japanese                          = 7.3,
      Missing                           = 5.0,
      Black_African_American            = 2.2,
      Others                            = 2.0,
      Chinese                           = 0.2
    ),
    disease_state  = "Radioiodine-refractory differentiated thyroid cancer (RR-DTC), with measurable target lesions by RECIST 1.1.",
    dose_range     = "Lenvatinib 18 or 24 mg orally once daily starting dose with protocol-driven reductions, or placebo.",
    regions        = "Multicenter: phase 3 study E7080-G000-303 and phase 2 post-marketing study E7080-G000-211.",
    n_observations = "3413 tumor-size observations from 558 RR-DTC patients.",
    tumor_burden   = "Baseline tumor size (sum of longest diameters of target lesions): mean 70.2 mm (SD 44.1), median 59.5 mm, range 10.2-331.0 mm (Table S2).",
    notes          = "Demographics from Majid 2024 Table S2 (tumor-size PK/PD population, N = 558). Sex 287 male / 271 female gives 48.6 percent female. Tumor assessments were performed every 8 weeks during the randomization phase (28-day cycles) and every 12 weeks in the extension phase, by independent reviewer for study 303 and by investigator for study 211. Two subjects (IDs 3266 and 3261) were excluded in Text S3 $DATA IGNORE."
  )

  ini({
    # ---- Tumor-dynamics parameters (Majid 2024 Table 3, Equation 4).
    # Time base is weeks. The inter-individual variability on KG, Emax and
    # lambda is ADDITIVE in the source (Text S3: KG = TVKG + ETA(1),
    # EMAX = THETA(2) + ETA(2), TMEFF = THETA(4) + ETA(4)) and Table 3
    # reports those three IIV rows as '(SD)', so these three typical values
    # are NOT log-transformed here. EC50 and KAng-2 carry exponential IIV
    # and Table 3 reports them as 'CV%'.
    kgrow_pop <- 0.00252;   label("Tumor growth rate KG (1/week)")                                             # Majid 2024 Table 3: KG = 0.00252 /week (%RSE 12.5; bootstrap 0.00248 [0.00167-0.00340])
    emax_pop  <- 0.0755;    label("Maximum lenvatinib effect on the tumor shrinkage rate Emax (1/week)")       # Majid 2024 Table 3: E max = 0.0755 /week (%RSE 3.80)
    lec50     <- log(1420); label("Lenvatinib average daily AUC giving half of Emax, EC50 (ng*h/mL)")          # Majid 2024 Table 3: EC50 = 1420 ng*h/mL (%RSE 8.31); Text S3 carries it as THETA(3) = 1.42 on the ug*h/mL scale
    lam_pop   <- 0.259;     label("Resistance appearance rate constant lambda (1/week)")                       # Majid 2024 Table 3: lambda = 0.259 /week (%RSE 1.79)
    ktie2_pop <- -0.0112;   label("Tumor size reduction rate constant for relative Tie-2 change (1/week)")     # Majid 2024 Table 3: KTie-2 = -0.0112 /week (bootstrap -0.0109 [-0.0208, -0.0005])
    kang2_pop <- -0.0144;   label("Tumor size reduction rate constant for relative Ang-2 change (1/week)")     # Majid 2024 Table 3: KAng-2 = -0.0144 /week (%RSE 3.68)

    # ---- Tie-2 and Ang-2 turnover parameters carried in from the companion
    # four-biomarker model. In the source these enter Text S3 as per-subject
    # empirical Bayes estimates in the input data ($INPUT BM02 MRT2 IMAX2
    # IC502 HILL2 DPSLO2 BM03 MRT3 IMAX3 IC503 HILL3 DPSLO3), i.e. the
    # tumor model was fitted sequentially with the biomarker parameters
    # already estimated. They are therefore fixed() here at the biomarker
    # model's typical values (Majid 2024 Table 2), converted from the
    # biomarker model's hours to this model's weeks where applicable.
    lrbase_tie2 <- fixed(log(14.6));            label("Tie-2 baseline BM0 (ng/mL)")                          # Majid 2024 Table 2: Tie-2 Baseline = 14.6
    lmrt_tie2   <- fixed(log(354 / 24 / 7));    label("Tie-2 mean residence time (week)")                    # Majid 2024 Table 2: Tie-2 MRT = 354 h; Text S3 'MRTT = MRT2/24/7'
    hill_tie2   <- fixed(0.313);                label("Hill coefficient for the Tie-2 drug effect")          # Majid 2024 Table 2: Tie-2 Hill coefficient = 0.313
    lrbase_ang2 <- fixed(log(3.36));            label("Ang-2 baseline BM0 (ng/mL)")                          # Majid 2024 Table 2: Ang-2 Baseline = 3.36
    lmrt_ang2   <- fixed(log(173 / 24 / 7));    label("Ang-2 mean residence time (week)")                    # Majid 2024 Table 2: Ang-2 MRT = 173 h; Text S3 'MRTA = MRT3/24/7'
    hill_ang2   <- fixed(4.27);                 label("Hill coefficient for the Ang-2 drug effect")          # Majid 2024 Table 2: Ang-2 Hill coefficient = 4.27
    lemax_bm    <- fixed(log(0.344));           label("Common biomarker maximum drug effect Emax (fraction)")# Majid 2024 Table 2: All common E max = 0.344
    lec50_bm    <- fixed(log(930));             label("Common biomarker EC50 on lenvatinib daily AUC (ng*h/mL)") # Majid 2024 Table 2: All common EC50 = 930 ng*h/mL
    ldp_slope   <- fixed(log(2.93e-6 * 24 * 7));label("Common biomarker linear disease-progression slope (1/week)") # Majid 2024 Table 2: DP slope 2.93e-6 /h; Text S3 'DPSLOT = DPSLO2*24*7'

    # ---- Inter-individual variability (Majid 2024 Table 3). The three
    # additive-eta rows are reported as standard deviations, so the
    # variances below are SD^2; the two exponential-eta rows are reported as
    # 'CV%' which the table footnote defines as 'square root of variance x
    # 100', so those variances are (CV%/100)^2. Both readings are
    # corroborated by the Text S3 $OMEGA initial estimates (variances):
    # KG 5.5e-5 vs 5.07e-5, Emax 0.003 vs 0.00212, EC50 0.55 vs 0.7225,
    # lambda 0.0125 vs 0.0106, KAng-2 0.75 vs 0.5535.
    etakgrow ~ 5.06944e-05; label("IIV variance on KG (additive eta)")                # Majid 2024 Table 3: KG (SD) = 0.00712 -> (0.00712)^2; Text S3 $OMEGA init .000055 [A]
    etaemax  ~ 0.00211600;  label("IIV variance on Emax (additive eta)")              # Majid 2024 Table 3: E max (SD) = 0.0460 -> (0.0460)^2; Text S3 $OMEGA init .003 [A]
    etalec50 ~ 0.72250000;  label("IIV variance on log EC50")                         # Majid 2024 Table 3: EC50 85.0 CV% -> (0.850)^2; Text S3 $OMEGA init .55 [P]
    etalam   ~ 0.01060900;  label("IIV variance on lambda (additive eta)")            # Majid 2024 Table 3: lambda (SD) = 0.103 -> (0.103)^2; Text S3 $OMEGA init .0125 [A]
    etakang2 ~ 0.55353600;  label("IIV variance on log KAng-2")                       # Majid 2024 Table 3: KAng-2 74.4 CV% -> (0.744)^2; Text S3 $OMEGA init 0.75 [P]

    # KTie-2 carries no IIV: Text S3 $OMEGA slot 7 is '0 FIX', and the paper
    # states 'Since IIV for KTie2 was very high (189 percent) despite
    # bootstrapped confidence intervals tightly encompassing the NONMEM
    # value, for the final model this parameter was not estimated.'

    # ---- Unit-variance etas that carry the baseline measurement error into
    # the tumor initial condition (Text S3 OMEGA slots 5 and 6, both
    # '1 FIX'). They are consumed only inside tumor(0).
    etaibase_prop ~ fixed(1)   # Majid 2024 Text S3 $OMEGA 5: '1 FIX'; used as ETA(5)*THETA(5)*BAS
    etaibase_add  ~ fixed(1)   # Majid 2024 Text S3 $OMEGA 6: '1 FIX'; used as ETA(6)*THETA(6)

    # ---- Residual error. Text S3 $ERROR builds the standard deviation as
    # W = IPRE*THETA(5) + THETA(6), a LINEAR sum of the proportional and
    # additive components, which is nlmixr2's combined1() form (not the
    # default combined2 sum-of-squares). Table 3's 'Proportional (CV%)' row
    # prints 0.0321, which is the fraction (3.21 percent), not a percentage;
    # the Text S3 $THETA initial estimate of 0.03 confirms the scale.
    propSd <- 0.0321; label("Proportional residual error on tumor size (fraction)")  # Majid 2024 Table 3: Proportional = 0.0321 (%RSE 2.29); Text S3 $THETA 5 init (0,.03) [Prop_RES]
    addSd  <- 1.21;   label("Additive residual error on tumor size (mm)")            # Majid 2024 Table 3: Additive SD = 1.21 mm (%RSE 1.98); Text S3 $THETA 6 init (0,1.25) [Add_RES]
  })

  model({
    # ---- Individual tumor-dynamics parameters. Additive etas on KG, Emax
    # and lambda; exponential etas on EC50 and KAng-2 (Text S3 $PK).
    kgrow <- kgrow_pop + etakgrow
    emax  <- emax_pop  + etaemax
    ec50  <- exp(lec50 + etalec50)
    lam   <- lam_pop   + etalam
    ktie2 <- ktie2_pop
    kang2 <- kang2_pop * exp(etakang2)

    # ---- Carried biomarker parameters (all fixed).
    rbase_tie2 <- exp(lrbase_tie2)
    mrt_tie2   <- exp(lmrt_tie2)
    rbase_ang2 <- exp(lrbase_ang2)
    mrt_ang2   <- exp(lmrt_ang2)
    emax_bm    <- exp(lemax_bm)
    ec50_bm    <- exp(lec50_bm)
    dp_slope   <- exp(ldp_slope)

    # ---- Tie-2 and Ang-2 turnover, re-integrated here exactly as in the
    # companion biomarker model (Majid 2024 Equation 2) but on a weekly time
    # base. Lenvatinib inhibits Kin, so both pools deplete under treatment.
    #
    # NOTE ON A DEPOSITED-CODE SLIP. Text S3 defines IMAX1T, HILLT and
    # IC50T for the Tie-2 arm and then never uses them: the EFFT line reads
    # 'IMAX1*(AUC+.1)**HILL/(IC50**HILL+(AUC+.1)**HILL)', i.e. it picks up
    # the VEGF-indexed columns instead. That is provably a copy-paste slip
    # rather than model intent, on three independent grounds: (i) the
    # Tie-2-specific variables are defined immediately above and left
    # unused; (ii) IC50 is the VEGF column on the ng*h/mL scale while AUC in
    # Text S3 is on the ug*h/mL scale, a 1000-fold mismatch that would drive
    # the Tie-2 drug effect to essentially zero and contradict the paper's
    # headline finding that Tie-2 is a significant predictor; and (iii) the
    # parallel Ang-2 arm three lines below correctly uses IMAX1A, HILLA and
    # IC50A. This model follows the paper's printed Equations 2 and 4 and
    # uses the Tie-2 parameters for the Tie-2 arm. See the vignette Errata.
    kout_tie2 <- 1 / mrt_tie2
    dp_tie2   <- rbase_tie2 * (1 + dp_slope * t)
    kin_tie2  <- dp_tie2 * kout_tie2
    eff_tie2  <- emax_bm * AUC_LEN^hill_tie2 / (ec50_bm^hill_tie2 + AUC_LEN^hill_tie2)

    kout_ang2 <- 1 / mrt_ang2
    dp_ang2   <- rbase_ang2 * (1 + dp_slope * t)
    kin_ang2  <- dp_ang2 * kout_ang2
    eff_ang2  <- emax_bm * AUC_LEN^hill_ang2 / (ec50_bm^hill_ang2 + AUC_LEN^hill_ang2)

    # ---- Drug-driven tumor shrinkage with exponentially decaying potency
    # (Majid 2024 Equation 4: DRV = Emax * AUC / (AUC + EC50 * e^(lambda*T))).
    drv <- emax * AUC_LEN / (AUC_LEN + ec50 * exp(lam * t))

    # ---- Biomarker-driven shrinkage terms. Text S3 forms each as the
    # relative change from baseline of the treated biomarker pool times its
    # rate constant. Both rate constants are negative, so a fall in Tie-2 or
    # Ang-2 contributes a positive shrinkage rate.
    rel_tie2 <- (tie2 - rbase_tie2) / rbase_tie2
    rel_ang2 <- (ang2 - rbase_ang2) / rbase_ang2

    # ---- Initial conditions.
    tie2(0)  <- rbase_tie2
    ang2(0)  <- rbase_ang2
    tumor(0) <- TUM_SLD * (1 + etaibase_prop * propSd) + etaibase_add * addSd

    d/dt(tie2)  <- kin_tie2 * (1 - eff_tie2) - kout_tie2 * tie2
    d/dt(ang2)  <- kin_ang2 * (1 - eff_ang2) - kout_ang2 * ang2
    d/dt(tumor) <- kgrow * tumor -
      (drv + ktie2 * rel_tie2 + kang2 * rel_ang2) * tumor

    # Text S3 also carries two 'untreated' reference states (DADT(2) and
    # DADT(4), the Tie-2 and Ang-2 pools with no drug effect). They are
    # integrated but never read by the tumor ODE or the $ERROR block, so
    # they are omitted here.

    tumor ~ add(addSd) + prop(propSd) + combined1()
  })
}
