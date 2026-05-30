Kawamura_2018_eribulin <- function() {
  description <- "Three-compartment IV PK driver coupled with a Friberg-style semi-mechanistic PD model for eribulin-induced neutropenia in Japanese patients with recurrent or metastatic breast cancer (Kawamura 2018). Plasma eribulin concentrations are produced by a 3-compartment model with linear elimination from the central compartment whose parameters are FIXED from the Majid 2014 popPK analysis (reproduced verbatim in Kawamura 2018 section 2.3): CL depends on body weight (allometric 0.75), serum albumin, alkaline phosphatase, and total bilirubin; V1, V2, V3 scale linearly with body weight; Q2 and Q3 scale allometrically with body weight. The PD layer (proliferation + three transit compartments + circulating neutrophils + feedback) is estimated on 401 patients / 5199 ANC measurements (Table 2): MTT = 104.5 h, Kprol = 0.0377 /h, Kout = 0.0295 /h, Gamma = 0.203, Slope = 0.0413 mL/ng (linear drug effect). Serum albumin influences Kprol (negative exponent), MTT (positive exponent), and Kout (positive exponent); a binary low-baseline-ANC indicator (BNEU3 = 1 when baseline ANC < 3000/uL) multiplies Kprol. IIV is reported on Kprol, Kout, and Slope (no IIV on MTT or Gamma). Additive residual error on circulating ANC (sigma = 1.15 cells/nL = 1150 cells/uL). Eribulin doses must be supplied in milligrams of eribulin-FREE-BASE equivalent (1.4 mg/m^2 mesilate = 1.23 mg/m^2 free base, conversion factor 1.23/1.4)."
  reference <- paste(
    "Kawamura T, Kasai H, Fermanelli V, Takahashi T, Sakata Y, Matsuoka T, Ishii M, Tanigawara Y. (2018).",
    "Pharmacodynamic analysis of eribulin safety in breast cancer patients using real-world postmarketing surveillance data.",
    "Cancer Sci 109(9):2822-2829. doi:10.1111/cas.13708.",
    "PK structure fixed from Majid O, Gupta A, Reyderman L, Olivo M, Hussein Z. (2014).",
    "Population pharmacometric analyses of eribulin in patients with locally advanced or metastatic breast cancer previously treated with anthracyclines and taxanes.",
    "J Clin Pharmacol 54(10):1134-1143. doi:10.1002/jcph.315.",
    "PD structure based on Friberg LE, Henningsson A, Maas H, Nguyen L, Karlsson MO. (2002).",
    "Model of chemotherapy-induced myelosuppression with parameter consistency across drugs.",
    "J Clin Oncol 20(24):4713-4721. doi:10.1200/JCO.2002.02.140;",
    "see modellib('Friberg_2002_paclitaxel') and modellib('Ozawa_2007_docetaxel') for the Friberg-family template.",
    sep = " "
  )
  vignette <- "Kawamura_2018_eribulin"
  units <- list(
    time          = "hour",
    dosing        = "mg",
    concentration = "mg/L",
    anc           = "cells/uL"
  )
  # Dosing unit note: dose is expected in milligrams of eribulin FREE BASE.
  # Convert from the clinical eribulin mesilate dose via D_free_base =
  # D_mesilate * 1.23/1.4 (Kawamura 2018 section 2.3 conversion factor;
  # 1.4 mg/m^2 mesilate = 1.23 mg/m^2 free base). Concentration unit note:
  # central/vc gives mg/L; this is 1000 * the ng/mL unit used by Majid 2014
  # and Kawamura 2018 -- the factor-of-1000 conversion is applied inside the
  # model() body where slope (paper unit mL/ng) multiplies Cc.

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on the Majid 2014 PK parameters with reference 68.7 kg (Kawamura 2018 section 2.3 equations): CL, Q2, Q3 scale as (WT/68.7)^0.75; V1, V2, V3 scale linearly as (WT/68.7). Kawamura 2018 Table 1 reports the dosing schedule per mg/m^2 (median 1.4 mg/m^2 eribulin mesilate, range 0.7-1.4) but does not tabulate body weight directly; users must supply WT as an explicit covariate per the underlying Majid 2014 model.",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Serum albumin at baseline.",
      units              = "g/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 4.0 g/dL on both the Majid 2014 PK CL term (positive exponent 0.946, Kawamura 2018 section 2.3) and the Kawamura 2018 PD Kprol/MTT/Kout terms (Kawamura 2018 Table 2 footnote: MTT = tvMTT * (ALB/4)^thetaALBMTT; Kprol = tvKprol * (ALB/4)^thetaALBKprol * (1 + BNEU3 * thetaBNEU3Kprol); Kout = tvKout * (ALB/4)^thetaALBKout). Cohort median 3.9 g/dL, range 1.3-5.1 g/dL (Table 1).",
      source_name        = "ALB"
    ),
    ALP = list(
      description        = "Serum alkaline phosphatase at baseline.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 132 U/L on the Majid 2014 PK CL term (negative exponent -0.209, Kawamura 2018 section 2.3). The Kawamura paper does not tabulate the cohort distribution for ALP but reports it was collected as part of postmarketing surveillance; 182 patients with missing ALB / ALP / BILI were excluded from the PD analysis (Kawamura 2018 Figure 2).",
      source_name        = "ALP"
    ),
    TBILI = list(
      description        = "Total serum bilirubin at baseline.",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 0.5 mg/dL on the Majid 2014 PK CL term (negative exponent -0.180, Kawamura 2018 section 2.3). Source column name in Kawamura 2018 is BILI; the canonical register uses TBILI. The Kawamura paper does not tabulate the cohort distribution for total bilirubin.",
      source_name        = "BILI"
    ),
    NEUT = list(
      description        = "Per-subject baseline absolute neutrophil count (observed value before drug administration), used as the initial condition of the Friberg proliferation / transit / circulating compartments AND as the BNEU input to the homeostatic-feedback term (NEUT/circ)^gamma AND as the threshold input to the binary BNEU3 indicator (1 if NEUT < 3000 cells/uL, else 0) that multiplies Kprol.",
      units              = "cells/uL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column name in Kawamura 2018 is BNEU (baseline neutrophils). Cohort median 3200 cells/uL, range 943-15 000 cells/uL (Table 1). The BNEU3 derived indicator in Kawamura 2018 Table 2 footnote uses threshold 3000 cells/uL (paper notation: 'If baseline neutrophil counts < 3000, then BNEU3 = 1; if >= 3000, then BNEU3 = 0'). The canonical NEUT register units are cells/mm^3 (= cells/uL numerically); this model documents cells/uL explicitly to match the paper's reporting convention.",
      source_name        = "BNEU"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 401L,
    n_studies       = 1L,
    age_range       = "26-84 years (median 58)",
    weight_range    = "not tabulated in Kawamura 2018; the underlying Majid 2014 popPK reference weight is 68.7 kg",
    sex_female_pct  = 100,
    race_ethnicity  = c(Asian = 100),
    disease_state   = "Recurrent or metastatic breast cancer (RBC/MBC) in Japanese patients receiving eribulin mesilate as a first eribulin treatment, with no concomitant granulocyte colony-stimulating factor (G-CSF). 401 of 608 surveyed patients met inclusion criteria (Kawamura 2018 Figure 2); 207 were excluded for missing ALB/ALP/BILI (n=182) or missing baseline neutrophil count (n=25).",
    dose_range      = "Eribulin mesilate 0.7-1.4 mg/m^2 IV (median 1.4 mg/m^2), delivered on per-patient schedules including standard (day 1 + day 8 q21d, n=275), biweekly (day 1 + day 15 q28d, n=64), triweekly (day 1 q21d, n=50), and other schedules (n=12 excluded from the schedule simulation). Approximately 30% of patients had dosing-frequency reductions in cycle 1 due to neutropenic toxicity (Kawamura 2018 Discussion). Doses must be supplied to this model in milligrams of eribulin FREE BASE; conversion from the mesilate dose is D_free_base = D_mesilate * 1.23/1.4 (Kawamura 2018 section 2.3).",
    regions         = "Japan (325 centers; postmarketing surveillance July-December 2011)",
    ecog_distribution = "ECOG performance status 0-1 / 2 / >=3 = 192 / 172 / 37 patients (Kawamura 2018 Table 1)",
    prior_chemotherapy = "Median 4 (range 0-13) prior chemotherapy regimens including taxanes; distribution 0 / 1 / 2-4 / >=5 = 11 / 31 / 194 / 165 patients (Kawamura 2018 Table 1)",
    notes           = "Postmarketing-surveillance cohort with retrospective PD analysis. Patients with high myelosuppression, known hypersensitivity to eribulin mesilate, or pregnancy were excluded from the surveillance per Good Post-Marketing Study Practice. 5199 neutrophil count measurements across 401 patients fed the PD analysis. Bootstrap validation: 182/200 runs converged, mean/final-estimate ratios 97.5-142.8% (Kawamura 2018 Table 2 right columns). Estimation was performed in Phoenix NLME 7.0 (FO method). Model parameter estimates were 'similar to those previously reported using premarketing clinical trial data' (the upstream van Hasselt 2013 eribulin-neutropenia popPK/PD analysis cited as ref 13)."
  )

  ini({
    # =========================================================================
    # PK layer -- Majid 2014 three-compartment model with linear elimination,
    # reproduced verbatim in Kawamura 2018 section 2.3. All PK parameters are
    # FIXED here because Kawamura 2018 did not estimate them on the
    # postmarketing-surveillance cohort -- they used the Majid 2014 equations as
    # a deterministic input to simulate plasma eribulin concentrations C(t).
    # Reference values: WT 68.7 kg, ALB 4.0 g/dL, ALP 132 U/L, BILI 0.5 mg/dL.
    # =========================================================================
    lcl  <- fixed(log(3.11)); label("CL: eribulin clearance (L/h)")                                  # Kawamura 2018 section 2.3 (Majid 2014 reproduction): CL[L/h] = 3.11 * (WT/68.7)^0.75 * (ALB/4.0)^0.946 * (ALP/132)^-0.209 * (BILI/0.5)^-0.180
    lvc  <- fixed(log(4.06)); label("V1: central volume of distribution (L)")                         # Kawamura 2018 section 2.3: V1[L] = 4.06 * (WT/68.7)
    lq   <- fixed(log(2.64)); label("Q2: intercompartmental clearance central <-> peripheral1 (L/h)") # Kawamura 2018 section 2.3: Q2[L/h] = 2.64 * (WT/68.7)^0.75
    lvp  <- fixed(log(2.42)); label("V2: peripheral1 volume of distribution (L)")                     # Kawamura 2018 section 2.3: V2[L] = 2.42 * (WT/68.7)
    lq2  <- fixed(log(6.60)); label("Q3: intercompartmental clearance central <-> peripheral2 (L/h)") # Kawamura 2018 section 2.3: Q3[L/h] = 6.60 * (WT/68.7)^0.75
    lvp2 <- fixed(log(121));  label("V3: peripheral2 volume of distribution (L)")                     # Kawamura 2018 section 2.3: V3[L] = 121 * (WT/68.7)

    # PK covariate exponents -- all FIXED from Majid 2014 (Kawamura 2018 section 2.3 equations).
    e_wt_cl    <- fixed(0.75);   label("Allometric exponent of WT on CL (fixed)")                        # Kawamura 2018 section 2.3 CL equation: (WT/68.7)^0.75
    e_alb_cl   <- fixed(0.946);  label("Power exponent of ALB on CL (fixed)")                            # Kawamura 2018 section 2.3 CL equation: (ALB/4.0)^0.946
    e_alp_cl   <- fixed(-0.209); label("Power exponent of ALP on CL (fixed)")                            # Kawamura 2018 section 2.3 CL equation: (ALP/132)^-0.209
    e_tbili_cl <- fixed(-0.180); label("Power exponent of BILI on CL (fixed)")                           # Kawamura 2018 section 2.3 CL equation: (BILI/0.5)^-0.180
    e_wt_q     <- fixed(0.75);   label("Allometric exponent of WT on Q2 (fixed)")                        # Kawamura 2018 section 2.3 Q2 equation: (WT/68.7)^0.75
    e_wt_q2    <- fixed(0.75);   label("Allometric exponent of WT on Q3 (fixed)")                        # Kawamura 2018 section 2.3 Q3 equation: (WT/68.7)^0.75
    e_wt_vc    <- fixed(1);      label("Linear exponent of WT on V1 (fixed)")                            # Kawamura 2018 section 2.3 V1 equation: linear (WT/68.7), i.e. exponent 1
    e_wt_vp    <- fixed(1);      label("Linear exponent of WT on V2 (fixed)")                            # Kawamura 2018 section 2.3 V2 equation: linear (WT/68.7), i.e. exponent 1
    e_wt_vp2   <- fixed(1);      label("Linear exponent of WT on V3 (fixed)")                            # Kawamura 2018 section 2.3 V3 equation: linear (WT/68.7), i.e. exponent 1

    # =========================================================================
    # PD layer -- Friberg-style myelosuppression model estimated on the 401-
    # patient cohort by Kawamura 2018. Final parameter estimates from Table 2.
    # =========================================================================
    lkprol <- log(0.0377); label("Kprol: neutrophil proliferation rate constant (1/h)")           # Kawamura 2018 Table 2: tvKprol = 0.0377 (bootstrap 95% CI 0.0303-0.0472)
    lmtt   <- log(104.5);  label("MTT: mean transit time through neutrophil maturation chain (h)") # Kawamura 2018 Table 2: tvMTT   = 104.5  (bootstrap 95% CI 82.1-124.1)
    lkout  <- log(0.0295); label("Kout: circulating-neutrophil elimination rate constant (1/h)")  # Kawamura 2018 Table 2: tvKout  = 0.0295 (bootstrap 95% CI 0.0142-0.0489)
    lgamma <- log(0.203);  label("Gamma: feedback exponent on (BNEU/circ) (unitless)")            # Kawamura 2018 Table 2: tvGamma = 0.203  (bootstrap 95% CI 0.157-0.239)
    lslope <- log(0.0413); label("Slope: linear drug-effect coefficient on eribulin concentration (mL/ng)") # Kawamura 2018 Table 2: tvSlope = 0.0413 (bootstrap 95% CI 0.0320-0.0496)

    # Covariate effects on PD parameters -- Kawamura 2018 Table 2 (point estimates) and
    # Table 2 footnote (functional forms):
    #   MTT   = tvMTT   * (ALB/4)^thetaALBMTT
    #   Kprol = tvKprol * (ALB/4)^thetaALBKprol * (1 + BNEU3 * thetaBNEU3Kprol)
    #   Kout  = tvKout  * (ALB/4)^thetaALBKout
    # where BNEU3 = 1 if baseline NEUT < 3000 cells/uL, else 0.
    e_alb_kprol   <- -0.759; label("Power exponent of ALB on Kprol (unitless; negative -> low albumin raises Kprol)") # Kawamura 2018 Table 2: thetaALBKprol  = -0.759 (bootstrap 95% CI -1.110 to -0.427)
    e_alb_mtt     <-  0.605; label("Power exponent of ALB on MTT (unitless; positive -> low albumin shortens MTT)")    # Kawamura 2018 Table 2: thetaALBMTT    =  0.605 (bootstrap 95% CI 0.278 to 0.973)
    e_alb_kout    <-  0.357; label("Power exponent of ALB on Kout (unitless; positive -> low albumin lowers Kout)")    # Kawamura 2018 Table 2: thetaALBKout   =  0.357 (bootstrap 95% CI -0.144 to 0.950; CI crosses zero)
    e_bneu3_kprol <-  0.0704; label("Additive multiplier of (NEUT<3000) on Kprol via (1 + BNEU3*theta) (unitless)")    # Kawamura 2018 Table 2: thetaBNEU3Kprol =  0.0704 (bootstrap 95% CI 0.0432 to 0.0953)

    # IIV -- Kawamura 2018 Table 2 reports omega^2 directly as variances on the
    # log-eta scale (parameters are log-transformed structural thetas, IIV is
    # exp(eta) multiplicative). No IIV reported on MTT or Gamma.
    etalkprol ~ 0.00417  # Kawamura 2018 Table 2: omega2_Kprol = 0.00417
    etalkout  ~ 0.374    # Kawamura 2018 Table 2: omega2_Kout  = 0.374
    etalslope ~ 0.163    # Kawamura 2018 Table 2: omega2_slope = 0.163

    # Residual error -- Kawamura 2018 Table 2: sigma = 1.15 cells/nL (additive
    # residual SD on circulating neutrophil count). Converting to the cells/uL
    # unit declared in units$anc: 1 nL = 10^-3 uL, so 1 cell/nL = 10^3 cells/uL,
    # i.e. 1.15 cells/nL = 1150 cells/uL. Sanity check: cohort median baseline
    # NEUT = 3200 cells/uL = 3.2 cells/nL, so sigma/baseline = 1.15/3.2 ~= 36%
    # of baseline ANC -- a large but plausible additive RUV for routine
    # clinical-laboratory neutrophil counts in this postmarketing cohort.
    addSd_ANC <- 1150; label("Additive residual SD on circulating ANC (cells/uL; paper sigma = 1.15 cells/nL = 1150 cells/uL)") # Kawamura 2018 Table 2: sigma = 1.15 / nL (bootstrap 95% CI 1.03-1.23)
  })

  model({
    # ---- PK individual parameters (typical-value only; all theta FIXED from Majid 2014) ----
    cl   <- exp(lcl)  * (WT/68.7)^e_wt_cl * (ALB/4.0)^e_alb_cl * (ALP/132)^e_alp_cl * (TBILI/0.5)^e_tbili_cl
    vc   <- exp(lvc)  * (WT/68.7)^e_wt_vc
    q    <- exp(lq)   * (WT/68.7)^e_wt_q
    vp   <- exp(lvp)  * (WT/68.7)^e_wt_vp
    q2   <- exp(lq2)  * (WT/68.7)^e_wt_q2
    vp2  <- exp(lvp2) * (WT/68.7)^e_wt_vp2

    kel <- cl  / vc
    k12 <- q   / vc
    k21 <- q   / vp
    k13 <- q2  / vc
    k31 <- q2  / vp2

    # ---- PD individual parameters with covariate effects ----
    # BNEU3 binary indicator: 1 if baseline NEUT below the 3000 cells/uL threshold, else 0.
    bneu3  <- (NEUT < 3000)
    kprol  <- exp(lkprol + etalkprol) * (ALB/4)^e_alb_kprol * (1 + bneu3 * e_bneu3_kprol)
    mtt    <- exp(lmtt)               * (ALB/4)^e_alb_mtt
    kout   <- exp(lkout  + etalkout)  * (ALB/4)^e_alb_kout
    gamma  <- exp(lgamma)
    slope  <- exp(lslope + etalslope)
    ktr    <- 4 / mtt

    # ---- Eribulin plasma concentration ----
    # Dose is supplied in mg of eribulin free base; central holds the amount in
    # mg; vc is in L; so Cc = central/vc is in mg/L (the nlmixr2 default for
    # small-molecule popPK output). Note 1 mg/L = 1000 ng/mL -- the unit used
    # by Majid 2014 and Kawamura 2018; the conversion is applied where the
    # slope parameter is used (see edrug below).
    Cc <- central / vc

    # ---- Drug effect (linear in concentration) ----
    # Kawamura 2018 reports the linear drug-effect coefficient as Slope = 0.0413
    # mL/ng (Table 2). To use this paper-literal slope value with Cc in mg/L,
    # convert Cc to ng/mL via the factor 1000 (1 mg/L = 1000 ng/mL):
    #   edrug[unitless] = slope[mL/ng] * Cc[ng/mL]
    #                   = slope[mL/ng] * Cc[mg/L] * 1000
    edrug <- slope * Cc * 1000

    # ---- Three-compartment IV PK (Kawamura 2018 Figure 1 lower panel; Majid 2014) ----
    d/dt(central)     <- -(kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # ---- Friberg myelosuppression chain (Kawamura 2018 section 2.3 equations) ----
    # Compartment naming follows the Friberg-family precedent in
    # nlmixr2lib (Friberg_2002_paclitaxel, Ozawa_2007_docetaxel):
    #   precursor1 = Prol (proliferating progenitor pool)
    #   precursor2..4 = Transit1..3 (maturation chain)
    #   circ = Neu (circulating neutrophils, observed)
    # Feedback term uses the per-subject baseline neutrophil count NEUT (paper's
    # BNEU) in the ratio (BNEU/Neu)^Gamma.
    d/dt(precursor1) <- kprol * precursor1 * (1 - edrug) * (NEUT / circ)^gamma - ktr * precursor1
    d/dt(precursor2) <- ktr   * precursor1 - ktr * precursor2
    d/dt(precursor3) <- ktr   * precursor2 - ktr * precursor3
    d/dt(precursor4) <- ktr   * precursor3 - ktr * precursor4
    d/dt(circ)       <- ktr   * precursor4 - kout * circ

    # ---- Initial conditions ----
    # All Friberg compartments start at the per-subject baseline neutrophil
    # count NEUT. Kawamura 2018 estimated Kprol = 0.0377, MTT = 104.5 (so
    # Ktr = 4/MTT = 0.0383), and Kout = 0.0295 INDEPENDENTLY -- not under the
    # Friberg 2002 / Ozawa 2007 steady-state constraint Kprol = Ktr = Kout.
    # The published estimates satisfy Kprol ~= Ktr (~1.5% mismatch) but
    # Kout < Ktr by ~23%, so a no-drug simulation starting from
    # precursor1..4 = circ = NEUT is NOT at exact steady state: ANC drifts up
    # to ~30% above NEUT around t = 100 h, undershoots to ~10-15% below NEUT
    # around t = 480 h, and settles to ~95% of NEUT at the asymptotic
    # equilibrium (~5% steady-state mismatch). This drift is an inherent
    # property of the published parameterisation -- not a transcription bug --
    # and is small compared to drug-driven perturbations (cycle-1 nadirs of
    # 500-1000 cells/uL at 21 days under the standard dosing schedule).
    # Initialising all five PD compartments at NEUT matches Kawamura 2018's
    # own simulation convention (Figure 3 / Figure 4) and the broader
    # Friberg-family convention in nlmixr2lib (Friberg_2002_paclitaxel,
    # Ozawa_2007_docetaxel) but the magnitude of the drift here is
    # appreciable enough to call out in the validation vignette.
    precursor1(0) <- NEUT
    precursor2(0) <- NEUT
    precursor3(0) <- NEUT
    precursor4(0) <- NEUT
    circ(0)       <- NEUT

    # ---- Observations ----
    # ANC is the observed circulating neutrophil count (cells/uL). Additive
    # residual error per Kawamura 2018 Table 2 (sigma = 1.15 cells/nL).
    ANC <- circ
    ANC ~ add(addSd_ANC)
  })
}
