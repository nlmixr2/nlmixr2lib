Patel_2015_oseltamivir <- function() {
  description <- paste(
    "Joint parent + metabolite population PK model for oral oseltamivir and its",
    "active metabolite oseltamivir carboxylate in adults with end-stage renal disease",
    "receiving automated peritoneal dialysis (Patel 2015). Oseltamivir follows a",
    "two-compartment disposition with lagged first-order absorption; the fraction",
    "F_met = 0.964 of the absorbed dose is converted to carboxylate by hepatic",
    "first-pass metabolism and enters an empirical first-pass compartment that",
    "releases carboxylate into the systemic circulation at ka_OC = 0.109 1/h, while",
    "the remaining (1 - F_met) is absorbed as unchanged prodrug into the parent",
    "central compartment and converted systemically at CL_pm = 9.63 L/h/70 kg.",
    "Oseltamivir carboxylate follows a one-compartment disposition (Vc fixed at",
    "32.8 L/70 kg from Rayner et al.) with FOUR parallel elimination arms, three of",
    "which the study measured directly: peritoneal dialysate clearance split by",
    "exchange modality (CCPD 0.319 vs CAPD 0.170 L/h/70 kg, gated by the",
    "time-varying RRT_CCPD_ACTIVE / RRT_CAPD_ACTIVE session indicators so the two",
    "alternate as the dialysis schedule cycles), renal clearance (0.736 L/h/70 kg,",
    "gated off in anuric subjects by ANURIA), and an unidentified 'other' route",
    "(0.319 L/h/70 kg). Dialysate and urine are carried as cumulative-amount",
    "collection compartments. All clearances and volumes are standardized to a 70 kg",
    "adult by fixed allometric exponents. Doses and amounts are in micrograms and",
    "concentrations in ug/L (= ng/mL), matching the units the source model was",
    "fitted in; a 75 mg oral dose is 75000 ug."
  )
  reference <- paste(
    "Patel K, Rayner CR, Giraudon M, Kamal MA, Morcos PN, Robson R, Kirkpatrick CM.",
    "Pharmacokinetics and safety of oseltamivir in patients with end-stage renal",
    "disease treated with automated peritoneal dialysis.",
    "Br J Clin Pharmacol. 2015;79(4):624-635. doi:10.1111/bcp.12526"
  )
  vignette <- "Patel_2015_oseltamivir"
  units    <- list(time = "h", dosing = "ug", concentration = "ug/L")

  # Patel 2015 Figure 1 compartment 2 ("First-pass (OC)") is an empirical
  # compartment that delays the appearance of first-pass-generated carboxylate in
  # plasma; the authors introduced it to account for the ~28 h delay between the
  # parent and metabolite peaks. It is not a Savic absorption-chain transit.
  # Named to match the same structural feature in Standing_2012_oseltamivir.R,
  # which declares it the same way.
  paper_specific_compartments <- "transit_oselcarb"

  compartmentData <- list(
    depot = list(
      analyte = "oseltamivir", units = "ug",
      specimen = "administration site", verified = TRUE
    ),
    transit_oselcarb = list(
      analyte = "oseltamivir carboxylate", units = "ug",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "oseltamivir", units = "ug",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "oseltamivir", units = "ug",
      specimen = "plasma", verified = TRUE
    ),
    central_oselcarb = list(
      analyte = "oseltamivir carboxylate", units = "ug",
      specimen = "plasma", verified = TRUE
    ),
    dialysate_oselcarb = list(
      analyte = "oseltamivir carboxylate", units = "ug",
      specimen = "dialysate", verified = TRUE
    ),
    urine_oselcarb = list(
      analyte = "oseltamivir carboxylate", units = "ug",
      specimen = "urine", verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight (kg). Drives allometric size scaling of every clearance and volume, referenced to 70 kg: (WT/70)^0.75 on CL_pm, CL_D, and all four oseltamivir carboxylate clearance arms; (WT/70)^1.0 on Vc_OP, Vp_OP, and Vc_OC.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Patel 2015 reports every disposition parameter in Table 2 with units 'l h-1 (70 kg)-1' or 'l (70 kg)-1' and states in Results (Pharmacokinetic modelling) that 'the estimated clearances and volumes were standardized to a 70 kg adult human to allow for future comparison with other subpopulations, including extrapolation to infants and children'. The paper does NOT print the exponents; 0.75 / 1.0 are the theory-based values that this idiom denotes and are carried here as fixed() -- see the vignette Assumptions and deviations section. Cohort weight range 60-92 kg (Table 1), so the terms are interpolating over a narrow band within this study and extrapolate beyond it only under the fixed-exponent assumption.",
      source_name        = "Bodyweight (kg)"
    ),
    RRT_CCPD_ACTIVE = list(
      description        = "Time-varying indicator that a continuous cycler-assisted peritoneal dialysis (CCPD) exchange is running (1) or not (0). Gates the CCPD arm of the oseltamivir carboxylate dialysate clearance.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no cycler-assisted exchange running)",
      notes              = "The study regimen (Methods, Treatments; confirmed by Table 3 footnote (section symbol)) ran three CCPD exchanges of 2.5 L over an 8 h daytime block, so RRT_CCPD_ACTIVE = 1 for 0-8 h of each 24 h cycle and 0 thereafter, with oseltamivir given immediately before the block started. Genuinely time-varying within every subject. The paper's Monte Carlo simulations re-schedule this gate across three different automated-peritoneal-dialysis prescriptions (CAPD-only, intermediate, and the intensive regimen studied), which is why the CCPD / CAPD split is load-bearing rather than cosmetic. Both gates are 0 during a dry period, leaving only the renal and other routes active.",
      source_name        = "CCPD"
    ),
    RRT_CAPD_ACTIVE = list(
      description        = "Time-varying indicator that a continuous ambulatory peritoneal dialysis (CAPD) exchange is dwelling (1) or not (0). Gates the CAPD arm of the oseltamivir carboxylate dialysate clearance.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no ambulatory exchange dwelling)",
      notes              = "Two CAPD exchanges of 2.0 L over the 16 h overnight block, so RRT_CAPD_ACTIVE = 1 for 8-24 h of each 24 h cycle in the studied regimen. The osmotic agent differed between the two overnight exchanges (dextrose 1.5-4.25% w/v for one, icodextrin 7.5% w/v for the other), but Patel 2015 Results found the two comparable for oseltamivir carboxylate -- geometric mean CAPD clearance 0.195 L/h with icodextrin versus 0.177 L/h with dextrose, from which the paper concludes 'icodextrin had minimal impact on oseltamivir pharmacokinetics' -- so a single CAPD gate carries both and no osmotic-agent covariate is needed.",
      source_name        = "CAPD"
    ),
    ANURIA = list(
      description        = "Per-subject binary indicator of anuria: 1 = anuric (no residual renal elimination), 0 = residual urine production. Switches the renal arm of oseltamivir carboxylate elimination off entirely via the (1 - ANURIA) gate.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (residual urine production preserved)",
      notes              = "Ascertained per patient in Patel 2015 Table 1 as a Yes / No 'Anuric' column; 5 of the 9 analysed patients were anuric and 4 produced urine. Note the orientation -- the model gates by (1 - ANURIA), because the covariate is 1 for the subgroup in which the renal arm is ABSENT. The paper reports urine output only as prose ranges ('generally < 1000 ml day-1', with one patient at 1463-2776 ml day-1) and never as a per-subject volume, which is why the binary rather than URINE_VOL_24H is the right column here. Residual renal clearance was NOT predictable from serum creatinine or estimated creatinine clearance: Results states that 'inclusion of creatinine clearance as a covariate effect did not significantly improve the model', and the Discussion is explicit that 'this effect was not predicted by serum creatinine measurement or estimated creatinine clearance'. So ANURIA is the only renal-function covariate the model can use, and among urine producers the 117% between-subject variability on the renal arm is genuinely unexplained.",
      source_name        = "Anuric"
    )
  )

  covariatesDataExcluded <- list(
    CREAT = list(
      description = "Baseline serum creatinine (umol/L), reported per patient in Patel 2015 Table 1 (range 320-1461 umol/L).",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened but NOT retained. Results (Pharmacokinetic modelling): 'inclusion of creatinine clearance as a covariate effect did not significantly improve the model and did not explain the random variability associated with clearance by this route (see Supplementary Figure S2)'. Documented here to preserve the provenance of the paper's covariate screen; carries no coefficient and is not referenced in model()."
    ),
    BSA = list(
      description = "Body surface area (m^2), reported per patient in Patel 2015 Table 1 (range 1.66-2.09 m^2, an inclusion criterion of 1.7-2.3 m^2).",
      units       = "m^2",
      type        = "continuous",
      notes       = "Screened but NOT retained. Results: 'Covariate effects were not explored owing to the small sample size and because visual inspection of diagnostic plots showed no significant parameter-covariate relationships.' Body size enters the model through WT allometry only."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 9,
    n_studies      = 1,
    age_range      = "24-70 years",
    age_median     = "mean 52.8 years (Results, Study population)",
    weight_range   = "60-92 kg",
    weight_median  = "79 kg (median of the 10 enrolled patients, Table 1)",
    sex_female_pct = 50,
    race_ethnicity = c(
      Caucasian = 20, `Pacific Islander` = 40, `Indian/African` = 30, Maori = 10
    ),
    disease_state  = "End-stage renal disease on stable peritoneal dialysis for at least 3 months, with total Kt/V > 1.7 and a peritoneal equilibration test not indicating low transporter status. Five patients were anuric and five had residual urine production.",
    renal_function = "Anuric (n = 5) or residual urine production (n = 5; n = 4 in the PK analysis set). Baseline serum creatinine 320-1461 umol/L. Baseline creatinine clearance in the three urine producers with reported values was < 6, 13, and 14 mL/min.",
    dose_range     = "Single 75 mg oral dose of oseltamivir (75000 ug) on day 1, administered immediately before the start of the automated peritoneal dialysis block.",
    regions        = "New Zealand (two specialist clinical study facilities; ClinicalTrials.gov NCT01556633).",
    n_observations = "Plasma oseltamivir and oseltamivir carboxylate predose and at 0.5, 1.33, 2, 2.5, 3, 4, 5, 6.67, 8, 10, 12, 14, 16, 20, 24, 28, 32, 48, 72, 96, 120, 144 and 168 h postdose. Dialysate collected at the beginning and end of each peritoneal-dialysis interval for the first 48 h and every 24 h thereafter; urine collected predose and over 0-24, 24-48, 48-72, 72-96, 96-120, 120-144 and 144-168 h. Values below the limit of quantification were handled by the Beal M3 method.",
    dialysis       = "Standardized aggressive automated peritoneal dialysis: three continuous cycler-assisted (CCPD) exchanges of 2.5 L over an 8 h daytime block, then two continuous ambulatory (CAPD) exchanges of 2.0 L over 16 h overnight. Dextrose 1.5-4.25% w/v for all exchanges except one 8 h CAPD exchange using icodextrin 7.5% w/v.",
    notes          = "Of 27 patients screened, 10 were enrolled and all 10 were included in the safety evaluation. Patient 2 was excluded from the pharmacokinetic analysis for a protocol violation (QTcF outside the permitted range after dosing), so the model was fitted to n = 9; the demographic percentages above are over the 10 enrolled patients of Table 1. All subjects were taking concomitant medication. The regimen was deliberately chosen as a worst-case high-clearance scenario rather than as typical practice."
  )

  ini({
    # ----------------------------------------------------------------------
    # Absorption and first-pass conversion (Patel 2015 Table 2, "Absorption").
    # ----------------------------------------------------------------------
    # F_met is carried on the LOGIT scale, not as an exponential-IIV fraction.
    # Table 2 gives F_met = 0.964 with BSV 16.9%, and Methods says BSV was
    # "estimated using an exponential variability model"; but an exponential eta
    # of that size on 0.964 puts F_met above 1 for ~41% of subjects, and Figure 1
    # sends ka_OP*(1 - F_met) into the parent central compartment, so those
    # subjects would get a negative parent flux. The logit encoding is the
    # canonical fix for exactly this leak (see the logitfm register entry) and
    # reproduces the same F_met spread as reading the 16.9% on the complementary
    # (1 - F_met) scale: 0.958-0.969 versus 0.957-0.970 at +/- 1 SD.
    logitfm <- 3.287572; label("Logit of the fraction of absorbed oseltamivir converted to carboxylate by first pass (unitless)")  # Patel 2015 Table 2: F_met = 0.964; logit(0.964) = 3.287572
    lka     <- log(1.36); label("Absorption rate constant of oseltamivir, ka_OP (1/h)")                                            # Patel 2015 Table 2: ka_OP = 1.36 1/h
    ltlag   <- log(0.485); label("Absorption lag time of oseltamivir, Alag_OP (h)")                                                # Patel 2015 Table 2: Alag_OP = 0.485 h
    lka_oselcarb <- log(0.109); label("Appearance rate constant of oseltamivir carboxylate from the first-pass compartment, ka_OC (1/h)")  # Patel 2015 Table 2: ka_OC = 0.109 1/h

    # ----------------------------------------------------------------------
    # Oseltamivir (parent) disposition (Patel 2015 Table 2).
    # ----------------------------------------------------------------------
    lcl_met <- log(9.63); label("Systemic clearance of oseltamivir to carboxylate in plasma, CL_pm (L/h per 70 kg)")  # Patel 2015 Table 2: CL_pm = 9.63 l/h/(70 kg)
    lvc     <- log(16.7); label("Central volume of distribution of oseltamivir, Vc_OP (L per 70 kg)")                 # Patel 2015 Table 2: Vc_OP = 16.7 l/(70 kg)
    lq      <- log(6.80); label("Intercompartmental clearance of oseltamivir, CL_D (L/h per 70 kg)")                  # Patel 2015 Table 2: CL_D = 6.80 l/h/(70 kg)
    lvp     <- log(307);  label("Peripheral volume of distribution of oseltamivir, Vp_OP (L per 70 kg)")              # Patel 2015 Table 2: Vp_OP = 307 l/(70 kg)

    # ----------------------------------------------------------------------
    # Oseltamivir carboxylate (metabolite) disposition (Patel 2015 Table 2).
    # Vc_OC was FIXED to the value established by Rayner et al. (reference 18)
    # to avoid a parameter-identifiability problem; a BSV term on it was still
    # estimated. Methods: "For oseltamivir carboxylate, the central volume of
    # distribution (Vc_OC) was fixed to that previously established
    # {32.8 l (70 kg)-1; Rayner et al. [18]} to avoid issues with parameter
    # identification."
    # ----------------------------------------------------------------------
    lvc_oselcarb <- fixed(log(32.8)); label("Central volume of distribution of oseltamivir carboxylate, Vc_OC (L per 70 kg; from Rayner et al.)")  # Patel 2015 Table 2: Vc_OC = 32.8 l/(70 kg) FIX
    lcl_ccpd_oselcarb  <- log(0.319); label("Clearance of oseltamivir carboxylate by a CCPD exchange, CL_OCCCPD (L/h per 70 kg)")   # Patel 2015 Table 2: CL_OCCCPD = 0.319 l/h/(70 kg)
    lcl_capd_oselcarb  <- log(0.170); label("Clearance of oseltamivir carboxylate by a CAPD exchange, CL_OCCAPD (L/h per 70 kg)")   # Patel 2015 Table 2: CL_OCCAPD = 0.170 l/h/(70 kg)
    lcl_renal_oselcarb <- log(0.736); label("Renal clearance of oseltamivir carboxylate in urine producers, CL_OCURINE (L/h per 70 kg)")  # Patel 2015 Table 2: CL_OCURINE = 0.736 l/h/(70 kg)
    lcl_other_oselcarb <- log(0.319); label("Clearance of oseltamivir carboxylate by other unidentified routes, CL_OCOTH (L/h per 70 kg)")  # Patel 2015 Table 2: CL_OCOTH = 0.319 l/h/(70 kg)

    # ----------------------------------------------------------------------
    # Allometric exponents. Patel 2015 standardizes every clearance and volume
    # to 70 kg (Table 2 units and the Results statement quoted in
    # covariateData$WT$notes) but never prints the exponents; 0.75 / 1.0 are
    # the theory-based values that the "standardized to a 70 kg adult" idiom
    # denotes, and are the same pair used by Standing_2012_oseltamivir.R for
    # this drug. Carried as fixed() because they are structural assumptions,
    # not estimates from this paper.
    # ----------------------------------------------------------------------
    e_wt_cl <- fixed(0.75); label("Allometric exponent on all clearances (unitless; not printed in the source)")
    e_wt_vc <- fixed(1.0);  label("Allometric exponent on all volumes (unitless; not printed in the source)")

    # ----------------------------------------------------------------------
    # Between-subject variability. Patel 2015 Table 2 reports BSV as CV%;
    # converted to log-scale variance by omega^2 = log(1 + CV^2), except for
    # logitfm whose 16.9% is carried as the SD on the logit scale (see the
    # logitfm comment above and the Jung_2024_clopidogrel precedent):
    # omega^2 = 0.169^2. Parameters shown as "-" in the BSV column
    # (ka_OC, CL_pm, CL_D, Vp_OP) carry no eta.
    # ----------------------------------------------------------------------
    etalogitfm ~ 0.028561            # Patel 2015 Table 2: BSV(F_met)   = 16.9%; carried on the logit scale, omega^2 = 0.169^2
    etalka     ~ 1.150857            # Patel 2015 Table 2: BSV(ka_OP)   = 147%;  omega^2 = log(1 + 1.47^2)
    etaltlag   ~ 0.254211            # Patel 2015 Table 2: BSV(Alag_OP) = 53.8%; omega^2 = log(1 + 0.538^2)
    etalvc     ~ 0.092325            # Patel 2015 Table 2: BSV(Vc_OP)   = 31.1%; omega^2 = log(1 + 0.311^2)
    etalvc_oselcarb ~ 0.139561       # Patel 2015 Table 2: BSV(Vc_OC)   = 38.7%; omega^2 = log(1 + 0.387^2). Estimated even though the typical value is FIXED.
    # Table 2 reports the SAME BSV of 7.70% for both dialysate arms, which most
    # likely reflects a single shared OMEGA in the original NONMEM run. They are
    # carried here as two independent etas -- the literal transcription of the
    # table -- because the two arms are never active at the same time (CCPD and
    # CAPD alternate), so their between-subject correlation cannot affect any
    # simulated profile. See the vignette Assumptions and deviations section.
    etalcl_ccpd_oselcarb  ~ 0.005911 # Patel 2015 Table 2: BSV(CL_OCCCPD)  = 7.70%; omega^2 = log(1 + 0.077^2)
    etalcl_capd_oselcarb  ~ 0.005911 # Patel 2015 Table 2: BSV(CL_OCCAPD)  = 7.70%; omega^2 = log(1 + 0.077^2)
    etalcl_renal_oselcarb ~ 0.862426 # Patel 2015 Table 2: BSV(CL_OCURINE) = 117%;  omega^2 = log(1 + 1.17^2)
    etalcl_other_oselcarb ~ 0.127655 # Patel 2015 Table 2: BSV(CL_OCOTH)   = 36.9%; omega^2 = log(1 + 0.369^2)

    # ----------------------------------------------------------------------
    # Residual unexplained variability (Patel 2015 Table 2, "Residual error").
    # Table 2 footnote: "Data were modelled with units of concentration (in
    # micrograms per litre) for plasma and units of amount (in micrograms) for
    # dialysate and urine." The dialysate observation carries BOTH a
    # proportional and an additive term; the other three carry one each.
    # ----------------------------------------------------------------------
    propSd <- 0.359; label("Proportional residual error on plasma oseltamivir (fraction)")                                     # Patel 2015 Table 2: CV_OP_plasma  = 35.9 CV%
    addSd_oselcarb <- 73.1; label("Additive residual error on plasma oseltamivir carboxylate (ug/L)")                          # Patel 2015 Table 2: SD_OC_plasma  = 73.1 ug/l
    propSd_dialysateOselcarb <- 0.140; label("Proportional residual error on dialysate oseltamivir carboxylate (fraction)")    # Patel 2015 Table 2: CV_OC_dial    = 14.0 CV%
    addSd_dialysateOselcarb  <- 11.6;  label("Additive residual error on dialysate oseltamivir carboxylate (ug)")              # Patel 2015 Table 2: SD_OC_dial    = 11.6 ug
    propSd_urineOselcarb <- 0.355; label("Proportional residual error on urine oseltamivir carboxylate (fraction)")            # Patel 2015 Table 2: CV_OC_urine   = 35.5 CV%
  })

  model({
    # Reference body weight for the allometric terms (Patel 2015 Table 2 units).
    ref_wt <- 70

    # ------------------------------------------------------------------
    # 1. Fraction converted to carboxylate by first pass, back-transformed
    #    from the logit scale so it is bounded in (0, 1) for every subject.
    # ------------------------------------------------------------------
    logitfm_ind <- logitfm + etalogitfm
    fm          <- 1 / (1 + exp(-logitfm_ind))

    # ------------------------------------------------------------------
    # 2. Individual parameters. Clearances scale as (WT/70)^0.75 and volumes
    #    as (WT/70)^1.0; ka_OP, Alag_OP and ka_OC are rate constants and are
    #    not size-scaled.
    # ------------------------------------------------------------------
    ka          <- exp(lka + etalka)
    tlag        <- exp(ltlag + etaltlag)
    ka_oselcarb <- exp(lka_oselcarb)

    cl_met <- exp(lcl_met) * (WT / ref_wt)^e_wt_cl
    vc     <- exp(lvc + etalvc) * (WT / ref_wt)^e_wt_vc
    q      <- exp(lq)  * (WT / ref_wt)^e_wt_cl
    vp     <- exp(lvp) * (WT / ref_wt)^e_wt_vc

    vc_oselcarb <- exp(lvc_oselcarb + etalvc_oselcarb) * (WT / ref_wt)^e_wt_vc

    cl_ccpd_oselcarb  <- exp(lcl_ccpd_oselcarb  + etalcl_ccpd_oselcarb)  * (WT / ref_wt)^e_wt_cl
    cl_capd_oselcarb  <- exp(lcl_capd_oselcarb  + etalcl_capd_oselcarb)  * (WT / ref_wt)^e_wt_cl
    cl_other_oselcarb <- exp(lcl_other_oselcarb + etalcl_other_oselcarb) * (WT / ref_wt)^e_wt_cl

    # Renal arm is gated OFF entirely in anuric subjects. Patel 2015 Results:
    # "In patients with residual urine production (n = 4), renal elimination was
    # the major route of metabolite clearance."
    cl_renal_oselcarb <- (1 - ANURIA) *
      exp(lcl_renal_oselcarb + etalcl_renal_oselcarb) * (WT / ref_wt)^e_wt_cl

    # Dialysate clearance in force at this time point. The two exchange types
    # alternate as the automated-peritoneal-dialysis schedule cycles, so exactly
    # one gate is 1 during dialysis and both are 0 during a dry period.
    cl_dialysate_oselcarb <-
      RRT_CCPD_ACTIVE * cl_ccpd_oselcarb + RRT_CAPD_ACTIVE * cl_capd_oselcarb

    # ------------------------------------------------------------------
    # 3. Micro-rate constants (1/h).
    # ------------------------------------------------------------------
    kel_met <- cl_met / vc          # parent -> carboxylate systemic conversion
    k12     <- q / vc               # parent central -> peripheral1
    k21     <- q / vp               # parent peripheral1 -> central
    kel_dialysate_oselcarb <- cl_dialysate_oselcarb / vc_oselcarb
    kel_renal_oselcarb     <- cl_renal_oselcarb     / vc_oselcarb
    kel_other_oselcarb     <- cl_other_oselcarb     / vc_oselcarb

    # ------------------------------------------------------------------
    # 4. ODE system (Patel 2015 Figure 1). States hold amounts in ug.
    #      depot              1: Gut               -- oseltamivir
    #      transit_oselcarb   2: First-pass (OC)   -- carboxylate
    #      central            3: Central (OP)      -- oseltamivir
    #      peripheral1        4: Peripheral (OP)   -- oseltamivir
    #      central_oselcarb   5: Central (OC)      -- carboxylate
    #      dialysate_oselcarb 6: Dialysate (OC)    -- carboxylate, cumulative
    #      urine_oselcarb     7: Urine (OC)        -- carboxylate, cumulative
    #    Mass leaving the gut is conserved: the two absorption arms sum to
    #    ka*fm + ka*(1 - fm) = ka.
    # ------------------------------------------------------------------
    d/dt(depot)            <- -ka * depot
    d/dt(transit_oselcarb) <-  ka * fm * depot - ka_oselcarb * transit_oselcarb
    d/dt(central)          <-  ka * (1 - fm) * depot -
                               kel_met * central -
                               k12 * central +
                               k21 * peripheral1
    d/dt(peripheral1)      <-  k12 * central - k21 * peripheral1
    d/dt(central_oselcarb) <-  ka_oselcarb * transit_oselcarb +
                               kel_met * central -
                               (kel_dialysate_oselcarb +
                                kel_renal_oselcarb +
                                kel_other_oselcarb) * central_oselcarb
    d/dt(dialysate_oselcarb) <- kel_dialysate_oselcarb * central_oselcarb
    d/dt(urine_oselcarb)     <- kel_renal_oselcarb     * central_oselcarb

    # ------------------------------------------------------------------
    # 5. Absorption lag on the oral dose (Patel 2015 Table 2: Alag_OP).
    # ------------------------------------------------------------------
    alag(depot) <- tlag

    # ------------------------------------------------------------------
    # 6. Observations. Plasma outputs are concentrations in ug/L (= ng/mL);
    #    dialysate and urine outputs are cumulative amounts in ug, matching
    #    the Table 2 footnote. The dialysate and urine states are cumulative,
    #    so a user reproducing per-interval collected amounts must reset them
    #    at each interval boundary (see the vignette).
    # ------------------------------------------------------------------
    Cc                <- central          / vc
    Cc_oselcarb       <- central_oselcarb / vc_oselcarb
    dialysateOselcarb <- dialysate_oselcarb
    urineOselcarb     <- urine_oselcarb

    Cc                ~ prop(propSd)
    Cc_oselcarb       ~ add(addSd_oselcarb)
    dialysateOselcarb ~ add(addSd_dialysateOselcarb) + prop(propSd_dialysateOselcarb)
    urineOselcarb     ~ prop(propSd_urineOselcarb)
  })
}
