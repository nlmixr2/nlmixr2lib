Courlet_2023_cabamiquine <- function() {
  description <- paste(
    "Semi-mechanistic joint blood- and liver-stage population PK/PD model",
    "for the Plasmodium elongation factor 2 (eEF2) inhibitor cabamiquine",
    "(formerly DDD107498 / M5717) in healthy adult men challenged with",
    "Plasmodium falciparum (Courlet 2023, Table 1 and Fig. 2). The PK layer",
    "is the three-compartment transit-absorption plus enterohepatic-",
    "recycling model of modellib('Courlet_2023_cabamiquine_pk'), whose",
    "individual parameters were carried forward. The PD layer is a parasite",
    "turnover model at two life-cycle stages: a liver compartment seeded",
    "with Finc * inoculum infected hepatocytes that grows at a fixed",
    "kgr,l = ln(30000)/(6*24) = 0.072 /h and releases merozoites into blood",
    "through a logistic transfer function Tr_lb(t) centred on T50 = 6 days,",
    "and a blood compartment seeded at baseline parasitemia P0 that grows at",
    "kgr,b. Cabamiquine concentration in the central compartment drives",
    "parasite killing at both stages through a sigmoid Emax term with a",
    "shared maximum kill rate kki and a shared inverse-exponential onset",
    "delay (1 - exp(-kt * t)) that lets the killing rate maximise over the",
    "first ~48 h after dosing. Blood-stage potency differs between the two",
    "challenge populations: EC50,b,SpzCh = EC50,b,IBSM * (1 - 0.83 *",
    "STUDY_SPZCH). Liver-stage potency (EC50,l = 0.66 ng/mL) is higher than",
    "blood-stage potency, reproducing the paper's finding of superior",
    "efficacy against the liver stage. Set STUDY_SPZCH = 0 to simulate the",
    "induced blood stage malaria (IBSM) design (blood inoculation only) and",
    "STUDY_SPZCH = 1 for the sporozoite challenge (SpzCh) design (liver",
    "inoculation only).",
    sep = " "
  )
  reference <- paste(
    "Courlet P, Wilkins JJ, Oeuvray C, Gao W, Khandelwal A.",
    "Semi-mechanistic population pharmacokinetic/pharmacodynamic modeling",
    "of a Plasmodium elongation factor 2 inhibitor cabamiquine for",
    "prevention and cure of malaria.",
    "Antimicrob Agents Chemother. 2023;67(12):e00891-23.",
    "doi:10.1128/aac.00891-23.",
    "PK layer carried forward from the same paper's population PK model;",
    "see modellib('Courlet_2023_cabamiquine_pk').",
    sep = " "
  )
  vignette <- "Courlet_2023_cabamiquine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Parasite life-cycle states and the killing-onset bookkeeping state are
  # paper-mechanistic (Courlet 2023 Fig. 2, right panel: "Population
  # parasitemia model") and have no canonical counterpart.
  paper_specific_compartments <- c(
    "parasite_liver",
    "parasite_blood",
    "kill_onset"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Courlet 2023 Fig. 2.
  compartmentData <- list(
    depot          = list(analyte = "cabamiquine", units = "mg", specimen = "administration site", verified = TRUE),
    central        = list(analyte = "cabamiquine", units = "mg", specimen = "whole blood",         verified = TRUE),
    peripheral1    = list(analyte = "cabamiquine", units = "mg", specimen = "whole blood",         verified = TRUE),
    peripheral2    = list(analyte = "cabamiquine", units = "mg", specimen = "whole blood",         verified = TRUE),
    gallbladder    = list(analyte = "cabamiquine", units = "mg", specimen = "bile",                verified = TRUE),
    parasite_liver = list(analyte = "Plasmodium falciparum liver-stage parasites", units = "parasites", specimen = "tissue", verified = TRUE),
    parasite_blood = list(analyte = "Plasmodium falciparum blood-stage parasites", units = "parasites/mL", specimen = "blood cell", verified = TRUE),
    kill_onset     = list(analyte = "killing-onset fraction (dimensionless)", units = "fraction", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling of the PK layer with exponents fixed at 0.75 on",
        "the apparent clearances and 1 on the apparent volumes (Courlet 2023",
        "Supplementary Material 1). The paper does not state the reference",
        "weight; the packaged model uses the conventional 70 kg. See",
        "vignette Assumptions and deviations."
      ),
      source_name        = "WT"
    ),
    DOSE_CABAMIQUINE_MG = list(
      description        = "Administered cabamiquine single oral dose",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-subject single oral dose, entering the PK layer as the",
        "empirical power effect on the apparent central volume",
        "V2/F = 2363 * (WT/70)^1 * DOSE_CABAMIQUINE_MG^(-0.50). The",
        "reference (centering) dose is not reported in the paper; the",
        "uncentred dose in mg is used. Doses are FREE BASE milligrams, as is",
        "the `amt` on the dose record: the IBSM cohorts were dosed as the",
        "succinate salt and must be converted at 1 mg salt = 0.797 mg free",
        "base (Methods 'Data set'), so the 150 / 400 / 800 mg salt cohorts",
        "are 119.6 / 318.8 / 637.6 mg free base. See vignette Assumptions",
        "and deviations and modellib('Courlet_2023_cabamiquine_pk')."
      ),
      source_name        = "DOSE"
    ),
    STUDY_SPZCH = list(
      description        = "Sporozoite challenge (SpzCh) study indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (induced blood stage malaria, IBSM, challenge)",
      notes              = paste(
        "1 = participant enrolled in the phase Ib sporozoite challenge study",
        "(inoculated with approximately 3,200 P. falciparum NF54",
        "sporozoites, so infection starts at the liver stage); 0 = subject",
        "enrolled in the phase I IBSM study (inoculated with 3D7-infected",
        "erythrocytes, so infection starts directly at the blood stage).",
        "Serves two roles. (1) It scales the blood-stage potency:",
        "EC50,b,SpzCh = EC50,b,IBSM * (1 - theta * STUDY_SPZCH) with",
        "theta = 0.83 (Courlet 2023 Table 1 footnote c). (2) It selects",
        "which life-cycle stage is inoculated, so the same model file",
        "reproduces both trial designs: the liver state is seeded with",
        "Finc * inoculum when STUDY_SPZCH = 1 and the blood state is seeded",
        "with P0 when STUDY_SPZCH = 0."
      ),
      source_name        = "SpzCh"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 61L,
    n_studies      = 2L,
    age_range      = "18-55 years",
    sex_female_pct = 0,
    race_ethnicity = c(White = 100),
    disease_state  = paste(
      "Healthy malaria-naive adult men experimentally challenged with",
      "Plasmodium falciparum. IBSM cohort (22 participants, study 1 part 2):",
      "inoculated with 3D7-infected erythrocytes 8 days before dosing.",
      "SpzCh cohort (39 participants, study 2; 9 placebo and 30 active):",
      "inoculated with approximately 3,200 NF54 sporozoites, with drug given",
      "2 h (early liver stage) or 96 h (late liver stage) after inoculation."
    ),
    dose_range     = paste(
      "IBSM: single oral doses of 150, 400, and 800 mg cabamiquine succinate",
      "salt (1 mg salt = 0.797 mg free base). SpzCh: single oral doses of",
      "30, 60, 80, 100, 150, 200, 400, and 800 mg free base."
    ),
    regions        = "Australia (QIMR Berghofer) and the Netherlands",
    trial_registration = "ClinicalTrials.gov NCT03261401 (study 1)",
    notes          = paste(
      "Parasitemia was measured by qPCR with a lower limit of quantification",
      "of 1 parasite/mL; 60% of parasitemia measures were below the LLOQ",
      "(72% of those from the SpzCh study), handled with the M3 method",
      "during estimation. The blood-stage layer used 501 parasitemia",
      "measures from 22 IBSM participants and the liver-stage layer 50",
      "evaluable measures from 39 SpzCh participants (Results 'PK/PD data').",
      "The Discussion notes the pooled parasitemia population was 61 healthy",
      "malaria-naive adult Caucasian men. Blood volume was assumed to be 5 L",
      "for all subjects (Methods 'Joint blood- and liver-stage modeling')."
    )
  )

  ini({
    # =========================================================================
    # PK layer. Structural values are Courlet 2023 Supplementary Material 1
    # ("PK model parameter estimates"); the PK model was fit first and its
    # individual parameters were carried forward into the PK/PD model
    # (Methods 'Blood-stage modeling': "Individual PK parameter estimates
    # were carried forward into the PKPD model"). See
    # modellib('Courlet_2023_cabamiquine_pk') for the standalone PK model.
    # =========================================================================
    lcl    <- log(17.8);   label("Apparent cabamiquine clearance CL/F for a 70 kg subject (L/h)")                 # Suppl. Material 1 row 'Clearance (CL/F, L/h) = 17.8 (RSE 5.64%)'
    lvc    <- log(2363);   label("Apparent central volume V2/F at a 1 mg dose for a 70 kg subject (L)")           # Suppl. Material 1 row 'Central volume of distribution (V2/F, L) = 2363 (RSE 3.23%)'
    lvp    <- log(2051);   label("Apparent first peripheral volume V3/F for a 70 kg subject (L)")                 # Suppl. Material 1 row 'Peripheral volume of distribution 1 (V3/F, L) = 2051 (RSE 7.76%)'
    lq     <- log(6.37);   label("Apparent first intercompartmental clearance Q2/F for a 70 kg subject (L/h)")    # Suppl. Material 1 row 'Intercompartmental clearance 1 (Q2/F, L/h) = 6.37 (RSE 0.0589%)'
    lvp2   <- log(2548);   label("Apparent second peripheral volume V4/F for a 70 kg subject (L)")                # Suppl. Material 1 row 'Peripheral volume of distribution 2 (V4/F, L) = 2548 (RSE 3.23%)'
    lq2    <- log(58.5);   label("Apparent second intercompartmental clearance Q3/F for a 70 kg subject (L/h)")   # Suppl. Material 1 row 'Intercompartmental clearance 2 (Q3/F, L/h) = 58.5 (RSE 3.08%)'
    lka    <- log(8.27);   label("First-order absorption rate constant ka, depot to central (1/h)")               # Suppl. Material 1 row 'Absorption rate constant (ka /h) = 8.27 (RSE 21.7%)'
    lktr   <- log(13.7);   label("Transit rate between absorption compartments ktr (1/h)")                        # Suppl. Material 1 row 'Transit rate between compartments (ktr /h) = 13.7 (RSE 5.66%)'
    lmtt   <- log(0.21);   label("Mean transit time MTT through the transit chain (h)")                           # Suppl. Material 1 row 'Mean transit time (MTT, h) = 0.21 (RSE 11.6%)'
    lkbm   <- fixed(log(0.0039)); label("Central to recycling (biliary) transfer rate constant k2g (1/h)")        # Suppl. Material 1 row 'Central to depot transfer rate constant (k2g /h) = 0.0039 Fixed'
    lkehc  <- fixed(log(15));     label("Recycling to depot release rate constant kg1 (1/h), gated at MTIME")     # Suppl. Material 1 row 'Depot to absorption transfer rate constant (kg1 /h) = 15 Fixed*' (footnote: 'Fixed to a high number; transfer is near-instant.')
    lmtime <- log(25.05);         label("Depot emptying time MTIME, onset of the recycling release (h)")          # Suppl. Material 1 row 'Depot emptying time (h) = 25.05 (RSE 13.0%)'
    lfdepot <- fixed(log(1));     label("Reference relative oral bioavailability F (fraction)")                   # Apparent-parameter anchor; Suppl. Material 1 reports IIV on F but no typical value

    e_wt_cl   <- fixed(0.75); label("Allometric exponent on CL/F, Q2/F and Q3/F with body weight (unitless)")     # Suppl. Material 1 rows 'Weight on CL/F = 0.75 Fixed', 'Weight on Q2/F = 0.75 Fixed', 'Weight on Q3/F = 0.75 Fixed'
    e_wt_vc   <- fixed(1);    label("Allometric exponent on V2/F, V3/F and V4/F with body weight (unitless)")     # Suppl. Material 1 rows 'Weight on V2/F = 1 Fixed', 'Weight on V3/F = 1 Fixed', 'Weight on V4/F = 1 Fixed'
    e_dose_vc <- -0.50;       label("Power exponent of administered dose on V2/F (unitless)")                     # Suppl. Material 1 row 'Dose on V2/F = -0.50 (RSE 5.18%)' (the main text quotes -0.530; see vignette Errata)

    # =========================================================================
    # PD layer -- blood stage. Values are Courlet 2023 Table 1 ("Parameter
    # estimates for the joint blood- and liver-stage model"). Table 1
    # footnote a marks the parameters that were fixed to the values estimated
    # during the earlier blood-stage-only modeling step; those carry fixed().
    # =========================================================================
    lp0     <- fixed(log(0.03));  label("Baseline blood parasitemia at inoculation P0 (parasites/mL)")            # Table 1 row 'Baseline parasitemia (P 0 , parasites/mL) a = 0.03' (footnote a: fixed to the blood-stage modeling estimate)
    lemax   <- fixed(log(0.21));  label("Maximum parasite killing rate kki (1/h), shared by blood and liver")     # Table 1 row 'Parasite killing rate (k ki /h) a = 0.21' (footnote a; Methods: 'Killing and delay rate constants (kki and kt, respectively) were assumed to be similar for blood and liver')
    lkgrow  <- fixed(log(0.064)); label("Blood-stage parasite growth rate constant kgr,b (1/h)")                  # Table 1 row 'Parasite growth rate (k gr,b /h) a = 0.064' (footnote a)
    lhill   <- fixed(log(19));    label("Hill coefficient gamma of the killing sigmoid (unitless)")               # Table 1 row 'Hill coefficient (gamma) a = 19' (footnote a; Methods: 'The Hill coefficient was fixed to a value obtained during previous analyses (12)')
    lkt     <- fixed(log(0.030)); label("Killing-onset delay rate constant kt (1/h); ln(2)/kt = 23.1 h")          # Table 1 row 'Delay rate constant (k t /h) a = 0.030' (footnote a; Results: 'an inverse-exponential model with a typical half-life of 23.1 h before maximal killing')
    lec50   <- log(7.60);         label("Blood-stage EC50 in the IBSM population, EC50,b,IBSM (ng/mL)")           # Table 1 row 'EC 50,b,IBSM (ng/mL) = 7.60 (RSE 10.2%)'
    e_spzch_ec50 <- 0.83;         label("Fractional reduction of the blood-stage EC50 in the SpzCh population (unitless)") # Table 1 footnote c: 'EC50,b,SpzCh = EC50,b,IBSM x (1 - theta x SpzCh) ... The effect parameter (i.e., theta) was estimated to be 0.83 with RSE 0.4%.' Check: 7.60 * (1 - 0.83) = 1.29 = the tabulated EC50,b,SpzCh

    # =========================================================================
    # PD layer -- liver stage. Table 1 footnote b marks the parameters that
    # were fixed based on literature data (reference 13); those carry fixed().
    # =========================================================================
    lec50_liver <- log(0.66);         label("Liver-stage EC50, EC50,l (ng/mL)")                                   # Table 1 row 'EC 50,l (ng/mL) = 0.66 (RSE 19.9%)'
    lkgrow_liver <- fixed(log(0.072)); label("Liver-stage parasite growth rate constant kgr,l (1/h)")             # Table 1 row 'Parasite growth rate in liver (k gr,l ) b = 0.072' (footnote b). Fig. 2: 'k gr,l fixed to ln(30000)/(6x24) = 0.072 h^-1'; Methods: '30,000 merozoites are released by each infected hepatocyte after 6 days'
    lklb    <- fixed(log(6));         label("Rate of merozoites invading erythrocytes on leaving hepatocytes klb (1/h)") # Table 1 row 'Rate of merozoites invading erythrocytes (k lb /h) b = 6' (footnote b); ln(2)/6 = 6.9 min, matching Methods 'a quick invasion of erythrocytes by merozoites (10 min)'
    lt50    <- fixed(log(144));       label("Time by which half the infected hepatocytes have burst T50 (h)")     # Table 1 row 'Time by which half the hepatocytes burst after infection (T 50 , h) b = 144' (footnote b); Fig. 2: 'with T50 fixed to 6 days'
    lsigma_lb <- fixed(log(0.1));     label("Spread of the hepatocyte-burst release around T50, sigma_lb (h)")    # Table 1 row 'Distribution of release around time (sigma lb ) b = 0.1' (footnote b); the logistic 5-95% span is +/- ln(19)*0.1 = +/- 17.7 min, matching Methods 'a release of the major part of merozoites at day 6 +/- 30 min'
    lfinc   <- fixed(log(0.0012));    label("Fraction of inoculated sporozoites that invade hepatocytes Finc (fraction)") # Table 1 row 'Fraction of infectious inoculated sporozoites invading hepatocytes (F inc ) b = 0.0012' (footnote b); 0.0012 * 3200 = 3.84, matching Methods 'a typical value of four infected hepatocytes with IIV'
    linoculum <- fixed(log(3200));    label("Number of sporozoites inoculated in the SpzCh challenge (sporozoites)") # Methods 'Data set': 'the drug administered 2 h or 96 h after inoculation of approximately 3,200 sporozoites (strain NF54)'
    lvblood <- fixed(log(5));         label("Blood volume used to convert liver-released merozoites to parasites/mL (L)") # Methods 'Joint blood- and liver-stage modeling': 'The blood volume was assumed to be 5 L for all subjects (14)'

    # =========================================================================
    # Inter-individual variability. Table 1 reports IIV as a coefficient of
    # variation in percent, converted to log-normal variance via
    # omega^2 = log(1 + CV^2). Every PD IIV carries Table 1 footnote a or b
    # (fixed from the blood-stage step or from literature) except the IIV on
    # EC50,b,IBSM, which was estimated (RSE 16.8%).
    #
    # The PK IIVs are Supplementary Material 1 values, carried forward with
    # the PK layer.
    # =========================================================================
    etalcl     ~ 0.07059  # Suppl. Material 1 'IIV on CL/F = 27%'  -> omega^2 = log(1 + 0.27^2) = 0.07059
    etalvc     ~ 0.06062  # Suppl. Material 1 'IIV on V2/F = 25%'  -> omega^2 = log(1 + 0.25^2) = 0.06062
    etalka     ~ 2.51231  # Suppl. Material 1 'IIV on ka = 340%'   -> omega^2 = log(1 + 3.40^2) = 2.51231
    etalkbm    ~ 0.50974  # Suppl. Material 1 'IIV on k2g = 82%'   -> omega^2 = log(1 + 0.82^2) = 0.50974
    etalktr    ~ 0.06062  # Suppl. Material 1 'IIV on ktr = 25%'   -> omega^2 = log(1 + 0.25^2) = 0.06062
    etalmtt    ~ 0.79470  # Suppl. Material 1 'IIV on MTT = 110%'  -> omega^2 = log(1 + 1.10^2) = 0.79470
    etalfdepot ~ 0.05161  # Suppl. Material 1 'IIV on F = 23%'     -> omega^2 = log(1 + 0.23^2) = 0.05161

    etalp0     ~ fixed(1.5853296)  # Table 1 'IIV on P0 a = 197'            -> omega^2 = log(1 + 1.97^2) = 1.5853296 (footnote a)
    etalemax   ~ fixed(0.0354637)  # Table 1 'IIV on kki a = 19'            -> omega^2 = log(1 + 0.19^2) = 0.0354637 (footnote a)
    etalec50   ~ 0.1844027         # Table 1 'IIV on EC50,b,IBSM = 45 (RSE 16.8%)' -> omega^2 = log(1 + 0.45^2) = 0.1844027 (estimated)
    etalkgrow  ~ fixed(0.0142973)  # Table 1 'IIV on kgr,b a = 12'          -> omega^2 = log(1 + 0.12^2) = 0.0142973 (footnote a)
    etalfinc   ~ fixed(0.0194104)  # Table 1 'IIV on Finc b = 14'           -> omega^2 = log(1 + 0.14^2) = 0.0194104 (footnote b)

    # Table 1 fixes the correlation between the blood and liver growth-rate
    # random effects at exactly 1 ('Correlation between kgr,b and kgr,l = 1
    # Fixed'; Results: 'A 100% correlation was assumed between the growth
    # rate constants in blood and liver'). A perfectly correlated 2x2 omega
    # block is exactly singular and makes rxSolve fail in the Cholesky
    # decomposition, so the liver growth rate reuses the SINGLE blood
    # random effect, rescaled by the ratio of the two reported IIV standard
    # deviations. This is algebraically identical to a correlation-1 block
    # and is applied in model() as
    #   kgrow_liver = exp(lkgrow_liver + sd_ratio_kgrow * etalkgrow).
    # Table 1 'IIV on kgr,l a = 11' -> omega^2 = log(1 + 0.11^2) = 0.0120274,
    # so sd_ratio = sqrt(0.0120274 / 0.0142973) = 0.9171881.
    sd_ratio_kgrow <- fixed(0.9171881); label("Ratio of the liver to blood growth-rate IIV standard deviations (unitless)") # Derived from Table 1 'IIV on kgr,l a = 11' and 'IIV on kgr,b a = 12' with the correlation fixed at 1

    # =========================================================================
    # Residual error.
    # =========================================================================
    expSd_parasitemia <- 1.68; label("Additive residual error on log parasitemia (log parasites/mL)")             # Table 1 row 'Residual error [parasitemia, log(/mL)] = 1.68 (RSE 3.30%)'
    propSd            <- 0.17; label("Proportional residual error on blood cabamiquine concentration (fraction)") # Suppl. Material 1 row 'Residual error (%) = 17 (RSE 1.88%)' (carried forward with the PK layer)
  })

  model({
    # =========================================================================
    # 1. Individual PK parameters (see Courlet_2023_cabamiquine_pk.R for the
    #    standalone PK model and its documentation).
    # =========================================================================
    allcl <- (WT / 70)^e_wt_cl
    allv  <- (WT / 70)^e_wt_vc

    cl  <- exp(lcl + etalcl) * allcl
    q   <- exp(lq)           * allcl
    q2  <- exp(lq2)          * allcl
    vc  <- exp(lvc + etalvc) * allv * DOSE_CABAMIQUINE_MG^e_dose_vc
    vp  <- exp(lvp)          * allv
    vp2 <- exp(lvp2)         * allv

    ka     <- exp(lka  + etalka)
    ktr    <- exp(lktr + etalktr)
    mtt    <- exp(lmtt + etalmtt)
    kbm    <- exp(lkbm + etalkbm)
    kehc   <- exp(lkehc)
    # Bare name is tmtime, not mtime: `mtime` is a reserved rxode2 keyword.
    tmtime <- exp(lmtime)
    fdepot <- exp(lfdepot + etalfdepot)

    # Savic transit-chain shape nn = mtt * ktr - 1 (typical 1.877), floored
    # at 1 -- the smallest value for which the transit input is continuously
    # differentiable at the dose time. Roughly a third of random draws fall
    # below that, and they are both unphysical (a chain shorter than one
    # compartment) and unsolvable. See the PK model file and the vignette
    # Assumptions and deviations.
    nn <- mtt * ktr - 1
    nn <- max(nn, 1)

    # =========================================================================
    # 2. Individual PD parameters.
    #    The blood-stage EC50 is scaled down in the sporozoite-challenge
    #    population per Table 1 footnote c.
    #    The liver growth rate reuses the blood growth-rate random effect,
    #    rescaled, to encode the correlation fixed at 1.
    # =========================================================================
    p0          <- exp(lp0 + etalp0)
    emax        <- exp(lemax + etalemax)
    kgrow       <- exp(lkgrow + etalkgrow)
    kgrow_liver <- exp(lkgrow_liver + sd_ratio_kgrow * etalkgrow)
    hill        <- exp(lhill)
    kt          <- exp(lkt)
    ec50_liver  <- exp(lec50_liver)
    klb         <- exp(lklb)
    t50         <- exp(lt50)
    sigma_lb    <- exp(lsigma_lb)
    finc        <- exp(lfinc + etalfinc)
    inoculum    <- exp(linoculum)
    vblood      <- exp(lvblood)

    ec50 <- exp(lec50 + etalec50) * (1 - e_spzch_ec50 * STUDY_SPZCH)

    # =========================================================================
    # 3. PK input. Explicit Savic 2007 transit density (the rxode2 transit()
    #    macro rescales its dose lookup by bioavailability and therefore
    #    delivers zero dose under f(depot) <- 0 in UI-form models), plus the
    #    enterohepatic recycling release gated at MTIME.
    # =========================================================================
    #    tad() and podo() both return NA until the first dose record, so the
    #    input and the recycling gate are evaluated only after a dose has been
    #    given; otherwise the NA would propagate into the parasite states. The
    #    `else` branch also covers the SpzCh placebo arm, which has no dose.
    tdose <- tad(depot)
    if (tdose > 0) {
      transit_in <- exp(log(podo(depot)) + log(fdepot) + log(ktr) +
                          nn * log(ktr * tdose) - ktr * tdose - lgamma(nn + 1))
      gate <- (tdose >= tmtime)
    } else {
      transit_in <- 0
      gate <- 0
    }

    d/dt(depot)       <-  transit_in + gate * kehc * gallbladder - ka * depot
    d/dt(central)     <-  ka * depot -
                            cl * central / vc -
                            q  * central / vc + q  * peripheral1 / vp -
                            q2 * central / vc + q2 * peripheral2 / vp2 -
                            kbm * central
    d/dt(peripheral1) <-  q  * central / vc - q  * peripheral1 / vp
    d/dt(peripheral2) <-  q2 * central / vc - q2 * peripheral2 / vp2
    d/dt(gallbladder) <-  kbm * central - gate * kehc * gallbladder

    f(depot) <- 0

    # Cabamiquine concentration in the central compartment, in ng/mL: this is
    # the Ccab that drives parasite killing at both stages (Fig. 2).
    Cc <- 1000 * central / vc

    # =========================================================================
    # 4. Killing-onset delay. Fig. 2 gives the killing rates as
    #      k_kill = kki * (1 - exp(-kt * t)) * Ccab^g / (Ccab^g + EC50^g)
    #    where t is time since dosing (Results: the delay is 'in the first
    #    48 h' after cabamiquine dosing). tad() / tafd() return NA before the
    #    first dose, which would poison the pre-dose parasite ODEs, so the
    #    factor is integrated instead as a state:
    #      d/dt(kill_onset) = kt * (1 - kill_onset) * dosed
    #    with kill_onset(0) = 0 and `dosed` switching on when drug first
    #    appears in the body. That yields exactly 1 - exp(-kt * (t - tdose))
    #    for t >= tdose and 0 before, with no NA.
    # =========================================================================
    drug_total <- depot + central + peripheral1 + peripheral2 + gallbladder
    dosed      <- (drug_total > 0)
    d/dt(kill_onset) <- kt * (1 - kill_onset) * dosed
    kill_onset(0)    <- 0

    # Sigmoid Emax killing terms, shared kki and shared onset delay,
    # stage-specific EC50 (Fig. 2, right panel).
    kkill_blood <- emax * kill_onset * Cc^hill / (Cc^hill + ec50^hill)
    kkill_liver <- emax * kill_onset * Cc^hill / (Cc^hill + ec50_liver^hill)

    # =========================================================================
    # 5. Parasite life cycle (Fig. 2, right panel).
    #    Tr_lb(t) is the logistic cumulative fraction of infected hepatocytes
    #    that have burst, centred on T50 with spread sigma_lb; merozoites
    #    leave the liver at klb * Tr_lb(t) and enter the blood, where the
    #    count is divided by the assumed 5 L (5000 mL) blood volume to give
    #    parasites/mL.
    # =========================================================================
    #    sigma_lb is only 0.1 h, so the raw logistic argument reaches
    #    (0 - 144) / 0.1 = -1440 at the start of the simulation and
    #    exp(1440) overflows to Inf. The overflow makes the solver's error
    #    control fail for a substantial minority of simulated subjects, so
    #    the exponent is clamped to a range that is still far beyond
    #    saturation (exp(500) ~ 1e217, i.e. tr_lb ~ 1e-217 versus 0).
    z_lb  <- max(min((t - t50) / sigma_lb, 500), -500)
    tr_lb <- 1 / (1 + exp(-z_lb))

    d/dt(parasite_liver) <- kgrow_liver * parasite_liver -
                              klb * tr_lb * parasite_liver -
                              kkill_liver * parasite_liver
    d/dt(parasite_blood) <- kgrow * parasite_blood -
                              kkill_blood * parasite_blood +
                              klb * tr_lb * parasite_liver / (1000 * vblood)

    # Inoculation. The SpzCh design starts the infection at the liver stage
    # (Finc * inoculum infected hepatocytes) and the IBSM design starts it
    # directly at the blood stage (baseline parasitemia P0); STUDY_SPZCH
    # selects between them so one model file reproduces both trial designs.
    parasite_liver(0) <- finc * inoculum * STUDY_SPZCH
    parasite_blood(0) <- p0 * (1 - STUDY_SPZCH)

    # =========================================================================
    # 6. Observations. Parasitemia is measured by qPCR in parasites/mL and
    #    Table 1 reports its residual error on the natural-log scale, i.e. an
    #    exponential (log-normal) residual on the linear scale.
    # =========================================================================
    parasitemia <- parasite_blood

    Cc          ~ prop(propSd)
    parasitemia ~ lnorm(expSd_parasitemia)
  })
}
