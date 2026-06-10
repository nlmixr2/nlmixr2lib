Stroh_2013_rolofylline <- function() {
  description <- paste(
    "Simultaneous three-output population PK model for IV rolofylline",
    "(adenosine A1 receptor antagonist) and both M1-trans and M1-cis",
    "active hydroxyl metabolites in 36 healthy adult male volunteers",
    "after single 1-60 mg IV infusions over 1-2 h (Stroh 2013, study",
    "KW-3902 IV-EU01). Parent rolofylline disposition is two-compartment",
    "with linear total clearance CL1 directing the entire parent loss to",
    "metabolite formation; the fraction FM of formed material is converted",
    "directly to M1-cis while (1 - FM) is converted to M1-trans. M1-trans",
    "disposition is two-compartment with distributional clearance CL4 and",
    "an additional unidirectional stereochemical interconversion clearance",
    "CL3 from M1-trans to M1-cis. M1-cis disposition is one-compartment",
    "with first-order clearance CL5. Random effects were not estimable for",
    "parent and M1-trans distributional clearances (CL2, CL4) or for the",
    "M1-cis central volume (V5) and were fixed to zero per the source.",
    "No covariates were retained: the cohort was a single Phase 1",
    "dose-escalation in white male volunteers and the paper screened no",
    "demographic effects. Structural identifiability of the final model",
    "was confirmed via the DAISY software tool (Bellu et al.).",
    sep = " "
  )
  reference   <- paste(
    "Stroh M, Hutmacher MM, Pang J, Lutz R, Magara H, Stone J.",
    "Simultaneous Pharmacokinetic Model for Rolofylline and both M1-trans",
    "and M1-cis Metabolites. AAPS J. 2013;15(2):498-504.",
    "doi:10.1208/s12248-012-9443-5",
    sep = " "
  )
  vignette    <- "Stroh_2013_rolofylline"
  units       <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  population <- list(
    species         = "human",
    n_subjects      = 36L,
    n_studies       = 1L,
    age_range       = "18-42 years (mean 29 years)",
    weight_range    = "60.2-90.0 kg (mean 73.3 kg)",
    height_range    = "161-191 cm (mean 175 cm)",
    sex_female_pct  = 0,
    race_ethnicity  = c(White = 100),
    disease_state   = "Healthy adult male volunteers",
    dose_range      = paste(
      "Single-dose IV infusion of rolofylline at 1, 2.5, 5, 10, 20, 30,",
      "40, 50, or 60 mg. Infusion duration 1 h for doses up to 30 mg and",
      "2 h for 40-60 mg; volumetric rate fixed at 1 mL/min.",
      sep = " "
    ),
    regions         = "Europe (study KW-3902 IV-EU01, Hoechst Marion Roussel bioanalysis)",
    n_observations  = paste(
      "1,914 post-dose plasma concentrations across rolofylline, M1-trans,",
      "and M1-cis combined; 223 BLQ observations (rolofylline 65,",
      "M1-trans 52, M1-cis 106) treated as missing.",
      sep = " "
    ),
    notes           = paste(
      "37 white male subjects were enrolled; 36 received study medication",
      "across the dose-escalation cohort. Estimation in NONMEM Version VI",
      "level 2.0 using ADVAN 7 (general linear model) and FOCE-I. Below-LOQ",
      "(LOQ = 0.5 ng/mL for all three analytes) treated as missing.",
      sep = " "
    )
  )

  ini({
    # ---------------- Parent (rolofylline) structural parameters --------------
    # Two-compartment disposition; CL1 (the "parent to metabolite clearance")
    # represents the total parent elimination, with all material routed to the
    # metabolite pools. Vss = V1 + V2 = 37.8 + 201 = 238.8 ~= 239 L matches
    # the paper's text-reported Vss of 239 L (Stroh 2013 Results).
    lvc          <- log(37.8);  label("V1: rolofylline central volume (L)")                                # Stroh 2013 Table II: V1 = 37.8 L, RSE 3.15%
    lcl          <- log(24.4);  label("CL1: rolofylline-to-metabolite total clearance (L/h)")              # Stroh 2013 Table II: CL1 = 24.4 L/h, RSE 4.39%
    lvp          <- log(201);   label("V2: rolofylline peripheral volume (L)")                             # Stroh 2013 Table II: V2 = 201 L, RSE 5.08%
    lq           <- log(13.2);  label("CL2: rolofylline distributional clearance (L/h)")                   # Stroh 2013 Table II: CL2 = 13.2 L/h, RSE 3.47%

    # ---------------- M1-trans metabolite structural parameters ---------------
    # CL3 is the unidirectional M1-trans -> M1-cis stereochemical interconversion
    # clearance (the only elimination route for M1-trans in the final model);
    # CL4 is the M1-trans distributional clearance.
    lvc_m1trans  <- log(26.1);  label("V3: M1-trans central volume (L)")                                   # Stroh 2013 Table II: V3 = 26.1 L, RSE 9.16%
    lcl_m1trans  <- log(19.6);  label("CL3: M1-trans -> M1-cis interconversion clearance (L/h)")           # Stroh 2013 Table II: CL3 = 19.6 L/h, RSE 6.89%
    lvp_m1trans  <- log(41.7);  label("V4: M1-trans peripheral volume (L)")                                # Stroh 2013 Table II: V4 = 41.7 L, RSE 7.29%
    lq_m1trans   <- log(28.4);  label("CL4: M1-trans distributional clearance (L/h)")                      # Stroh 2013 Table II: CL4 = 28.4 L/h, RSE 11.7%

    # ---------------- M1-cis metabolite structural parameters -----------------
    # M1-cis disposition is one-compartment; CL5 is its first-order elimination.
    lvc_m1cis    <- log(3.78);  label("V5: M1-cis central volume (L)")                                     # Stroh 2013 Table II: V5 = 3.78 L, RSE 34.13%
    lcl_m1cis    <- log(91.6);  label("CL5: M1-cis elimination clearance (L/h)")                           # Stroh 2013 Table II: CL5 = 91.6 L/h, RSE 6.77%

    # ---------------- Branching fraction --------------------------------------
    # FM is the apparent fraction of rolofylline metabolised directly to M1-cis
    # (the complementary 1 - FM goes to M1-trans). Logit / bounded forms were
    # not used in the source; exp() IIV is applied per the paper's "exponential
    # random-effects" wording.
    lfm          <- log(0.194); label("FM: fraction of parent metabolised directly to M1-cis (unitless)")  # Stroh 2013 Table II: FM = 0.194, RSE 11.8%

    # ---------------- Inter-individual variability (log-scale variances) ------
    # All etas are exponential per "The inter-individual variability (IIV) was
    # captured using exponential random-effects terms" (Stroh 2013 Methods).
    # The paper reports BSV as CV%; log-scale variance is omega^2 = log(CV^2 + 1)
    # with CV on the decimal scale (e.g., 12.3% -> CV = 0.123).
    # Random effects fixed to zero in the source: CL2 (parent distributional),
    # CL4 (M1-trans distributional), V5 (M1-cis central) - omitted here.
    etalvc         ~ 0.015016   # 12.3% CV; log(0.123^2 + 1); Stroh 2013 Table II BSV(V1) = 12.3%, RSE 29.8%
    etalcl         ~ 0.044780   # 21.4% CV; log(0.214^2 + 1); Stroh 2013 Table II BSV(CL1) = 21.4%, RSE 22.7%
    etalvp         ~ 0.047268   # 22.0% CV; log(0.220^2 + 1); Stroh 2013 Table II BSV(V2) = 22.0%, RSE 51.0%
    etalvc_m1trans ~ 0.255027   # 53.9% CV; log(0.539^2 + 1); Stroh 2013 Table II BSV(V3) = 53.9%, RSE 24.1%
    etalcl_m1trans ~ 0.164633   # 42.3% CV; log(0.423^2 + 1); Stroh 2013 Table II BSV(CL3) = 42.3%, RSE 21.1%
    etalvp_m1trans ~ 0.102786   # 32.9% CV; log(0.329^2 + 1); Stroh 2013 Table II BSV(V4) = 32.9%, RSE 30.7%
    etalcl_m1cis   ~ 0.142239   # 39.1% CV; log(0.391^2 + 1); Stroh 2013 Table II BSV(CL5) = 39.1%, RSE 27.8%
    etalfm         ~ 0.156055   # 41.1% CV; log(0.411^2 + 1); Stroh 2013 Table II BSV(FM) = 41.1%, RSE 33.4%

    # ---------------- Residual error (per output) -----------------------------
    # Parent: proportional only (additive omitted in Stroh 2013 Table II).
    # M1-trans and M1-cis: combined additive + proportional. CV% on the
    # proportional term is taken as the proportional SD on the linear scale.
    propSd         <- 0.261;    label("Rolofylline (parent) proportional residual SD (fraction)")          # Stroh 2013 Table II: 26.1% CV, RSE 10.6%
    propSd_m1trans <- 0.174;    label("M1-trans proportional residual SD (fraction)")                      # Stroh 2013 Table II: 17.4% CV, RSE 17.3%
    addSd_m1trans  <- 0.217;    label("M1-trans additive residual SD (ng/mL)")                             # Stroh 2013 Table II: 0.217 ng/mL, RSE 177% (essentially at the noise floor; retained as reported)
    propSd_m1cis   <- 0.150;    label("M1-cis proportional residual SD (fraction)")                        # Stroh 2013 Table II: 15.0% CV, RSE 20.4%
    addSd_m1cis    <- 0.614;    label("M1-cis additive residual SD (ng/mL)")                               # Stroh 2013 Table II: 0.614 ng/mL, RSE 39.0%
  })

  model({
    # ---- Individual structural parameters ------------------------------------
    # Distributional clearances (CL2, CL4) and the M1-cis central volume (V5)
    # carry no eta because the paper fixed those random effects to zero.
    vc          <- exp(lvc         + etalvc)
    cl          <- exp(lcl         + etalcl)
    vp          <- exp(lvp         + etalvp)
    q           <- exp(lq)
    vc_m1trans  <- exp(lvc_m1trans + etalvc_m1trans)
    cl_m1trans  <- exp(lcl_m1trans + etalcl_m1trans)
    vp_m1trans  <- exp(lvp_m1trans + etalvp_m1trans)
    q_m1trans   <- exp(lq_m1trans)
    vc_m1cis    <- exp(lvc_m1cis)
    cl_m1cis    <- exp(lcl_m1cis   + etalcl_m1cis)
    fm          <- exp(lfm         + etalfm)

    # ---- ODE system ----------------------------------------------------------
    # IV infusion enters `central`. The total parent loss is cl / vc * central;
    # a fraction fm of that loss appears as M1-cis directly while (1 - fm)
    # appears as M1-trans. M1-trans further interconverts unidirectionally to
    # M1-cis at clearance cl_m1trans. Distribution between central /
    # peripheral compartments uses inter-compartmental clearance q (parent) and
    # q_m1trans (M1-trans). M1-cis is one-compartment with first-order CL5.
    # No molecular-weight correction is applied between species; the paper's
    # NONMEM ADVAN 7 mass-balance treats amounts in compatible units.
    d/dt(central)             <- -(cl + q) / vc * central + q / vp * peripheral1
    d/dt(peripheral1)         <- q / vc * central - q / vp * peripheral1
    d/dt(central_m1trans)     <- (1 - fm) * cl / vc * central -
                                 (cl_m1trans + q_m1trans) / vc_m1trans * central_m1trans +
                                 q_m1trans / vp_m1trans * peripheral1_m1trans
    d/dt(peripheral1_m1trans) <- q_m1trans / vc_m1trans * central_m1trans -
                                 q_m1trans / vp_m1trans * peripheral1_m1trans
    d/dt(central_m1cis)       <- fm * cl / vc * central +
                                 cl_m1trans / vc_m1trans * central_m1trans -
                                 cl_m1cis / vc_m1cis * central_m1cis

    # ---- Observations --------------------------------------------------------
    # Dosing is in mg and volumes are in L; central / vc gives mg/L. Multiply
    # by 1000 to express in ng/mL to match the paper's bioanalysis units.
    Cc         <- (central          / vc        ) * 1000
    Cc_m1trans <- (central_m1trans  / vc_m1trans) * 1000
    Cc_m1cis   <- (central_m1cis    / vc_m1cis  ) * 1000

    Cc         ~ prop(propSd)
    Cc_m1trans ~ add(addSd_m1trans) + prop(propSd_m1trans)
    Cc_m1cis   ~ add(addSd_m1cis)   + prop(propSd_m1cis)
  })
}
