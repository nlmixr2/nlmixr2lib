Jiao_2008_mycophenolic_acid <- function() {
  description <- "Population PK model with enterohepatic circulation (EHC) for mycophenolic acid (MPA) and its 7-O-glucuronide metabolite (MPAG) in healthy Chinese male volunteers after a single 500 mg oral dose of mycophenolate mofetil (MMF, Cellcept). Five-compartment chain model (Figure 2 of Jiao 2008): a gastrointestinal depot, a two-compartment MPA disposition (central + peripheral), a one-compartment MPAG disposition (central_mpag), and a gallbladder accumulation compartment (gallbladder_mpag). First-order absorption with an absorption-lag time. Complete (fm = 1, fixed) one-pass conversion of MPA to MPAG by glucuronidation; MPAG is renally cleared in parallel with biliary excretion into the gallbladder. EHC is encoded as time-gated bolus releases of the gallbladder pool back into the GI depot at two postprandial meal times (4 and 10 h post-dose, study-1 design), with rate constant k51 acting over a 0.01 h window; the recycled MPAG is reabsorbed via the same first-order ka as the oral dose. The fraction of MPAG biliary-routed at the branch is encoded as EHCP = k45 / (k40 + k45). Body-weight scaling: paper Eq 5 (linear-proportional 'slope without intercept') with reference 65.5 kg applied to CL_MPA/F, Q/F, and V_3/F via fixed allometric exponent 1. Cross-parameter IIV linkage: eta(CL_MPAG/F) = psi_q_cl_mpag * eta(Q/F) reproduces the paper's joint eta structure where psi_q_cl_mpag is the paper's 'q' parameter. UGT1A9 polymorphisms were screened but not retained in the final model (no significant effect)."
  reference <- paste(
    "Jiao Z, Ding JJ, Shen J, Liang HQ, Zhong LJ, Wang Y, Zhong MK, Lu WY.",
    "Population pharmacokinetic modelling for enterohepatic circulation",
    "of mycophenolic acid in healthy Chinese and the influence of",
    "polymorphisms in UGT1A9.",
    "Br J Clin Pharmacol. 2008;65(6):893-907.",
    "doi:10.1111/j.1365-2125.2008.03109.x.",
    sep = " "
  )
  vignette <- "Jiao_2008_mycophenolic_acid"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight (kg).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used as a linear-proportional ('slope without intercept') scalar on CL_MPA/F, Q/F, and V_3/F per Jiao 2008 Methods Eq 5: P_i = theta1 * (WT_i / WT_m), with reference WT_m = 65.5 kg (cohort median; Table 1). Encoded with a fixed allometric exponent of 1, i.e. exponent fixed at the value implied by Eq 5. The paper screened a linear-with-intercept form (Eq 6) and a power form (Eq 7); Results: 'the slope linear model without intercept was selected according to the OFV value'. Inclusion of WT reduced IIV on Q/F, V_3/F and CL_MPA/F by 26.5%, 34.8% and 15.9% respectively (Results).",
      source_name        = "WT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age (years).",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in the GAM and stepwise covariate analysis (Methods) but not retained in the final model."
    ),
    HT = list(
      description = "Height (m).",
      units       = "m",
      type        = "continuous",
      notes       = "Screened in the GAM and stepwise covariate analysis (Methods) but not retained in the final model."
    ),
    CREAT = list(
      description = "Serum creatinine (umol/L).",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened in the GAM analysis; only CRCL (derived from CREAT) was carried into the stepwise procedure, and CRCL was rejected in the final model (effect on CLMPAG/F was not retained)."
    ),
    ALB = list(
      description = "Serum albumin (g/L).",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened in the GAM and stepwise covariate analysis but not retained in the final model."
    ),
    HGB = list(
      description = "Haemoglobin (g/L).",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened in the GAM and stepwise covariate analysis but not retained in the final model."
    ),
    CRCL = list(
      description = "Creatinine clearance by Cockcroft-Gault (mL/min).",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened in the GAM (identified as influential on CLMPAG/F) but not retained in the final stepwise / backwards-elimination model; weight on CL_MPA/F, Q/F and V_3/F dominated the OFV reduction (Results)."
    ),
    UGT1A9 = list(
      description = "UGT1A9 promoter / coding-region single-nucleotide polymorphism panel: -2208C/T, -2152C/T, -2141C/T, -1887T/G, -1818T/C, -665C/T, -440T/C, -331C/T, -275T/A, -109_-98 T(n) insertion, -87G/A, +98T/C (UGT1A9*3).",
      units       = "categorical (per-SNP allele)",
      type        = "categorical",
      notes       = "Genotyped per Methods and tested for an effect on CLMPA/F via Eq 8 (paper notation). Polymorphisms of UGT1A9 did not show any influence on CLMPA/F (Results: '-275T/A and -2152C/T polymorphisms were not found in this study cohort and only one individual was identified as a carrier of a functional UGT1A9*3 allele'). Not retained in the final model."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 42L,
    n_studies      = 2L,
    n_observations = "590 MPA + 589 MPAG plasma concentrations (Results).",
    age_range      = "19-26 years",
    age_median     = "21 years (Table 1)",
    weight_range   = "56.5-89.0 kg",
    weight_median  = "65.5 kg (Table 1)",
    sex_female_pct = 0,
    race_ethnicity = "Chinese (100%)",
    disease_state  = "Healthy adult male volunteers.",
    dose_range     = "Single 0.5 g (500 mg) oral dose of mycophenolate mofetil (MMF, Cellcept 0.25 g capsules, Shanghai Roche) given as two 0.25 g test or two 0.25 g reference formulations after an overnight (>= 10 h) fast.",
    regions        = "China (Huashan Hospital, Fudan University, Shanghai).",
    sampling_window= "Predose and 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12, 24, 36, 48 h postdose (study 1); predose and 0.17, 0.33, 0.5, 0.75, 1, 2, 3, 4, 5, 6, 8, 10, 12, 24, 36, 48 h postdose (study 2). Standardized lunch at 4 h post-dose in both studies; standardized dinner at 10 h (study 1) or 9.5 h (study 2); next-day breakfast at 24.25 h.",
    notes          = "Pharmacokinetic data pooled across two open-label, single-dose, randomized crossover bioequivalence studies (20 + 22 healthy volunteers). Only the Cellcept reference-formulation concentration-time data were used for the population EHC model. UGT1A9 genotype distribution given in Table 2: -2208, -2152, -2141, -665, -275 were monomorphic; -1818 T/T 23.8%, T/C 47.6%, C/C 28.6%; -1887 T/T 83.3%, T/G 14.3%, G/G 2.4%; -440/-331 in linkage with -440 T/C 4.8%, C/C 95.2% / -331 C/T 4.8%, T/T 95.2%; -109_-98 T(n) 9/9 47.6%, 9/10 4.8%, 10/10 47.6%; -87 G/G 97.6%, G/A 2.4%; +98 (UGT1A9*3) T/T 97.6%, T/C 2.4%. Estimation: first-order conditional with interaction (FOCE-I) for the final-model estimates reported in Table 3."
  )

  # Implementation notes (see vignette Assumptions and deviations for the
  # full justification):
  # * Compartment naming. depot replaces paper's compartment 1 (GI tract);
  #   central / peripheral1 are MPA central / peripheral (paper 2 / 3);
  #   central_mpag is MPAG central (paper 4); gallbladder_mpag is the
  #   gallbladder pool of MPAG (paper 5). The _mpag suffix matches the
  #   registered metabolite-suffix list (R/conventions.R::registeredMetabolites);
  #   the bare 'gallbladder' name is the canonical EHC compartment, and the
  #   '_mpag' suffix on the gallbladder state documents that the held species
  #   is MPAG (the EHC re-entry to the GI lumen is faithfully MPAG, with the
  #   paper's fm = 1 reconversion to MPA on re-absorption captured by routing
  #   the released amount through depot).
  # * Cross-parameter eta linkage. Table 3 reports an additional parameter
  #   'q' (= 1.33, RSE 27.2%) defined as eta(CLMPAG/F) = q * eta(Q/F),
  #   reducing the joint IIV structure to a single random effect on Q/F that
  #   scales onto CL_MPAG/F. Encoded inside model() as
  #   cl_mpag <- exp(lcl_mpag + psi_q_cl_mpag * etalq), with psi_q_cl_mpag
  #   estimated. No separate eta on lcl_mpag is declared. This deterministic
  #   linear linkage is the paper-faithful encoding of a singular perfect
  #   correlation between the two random effects.
  # * Gallbladder emptying. Two EHC events were modelled at 4 and 10 h
  #   post-dose with bolus emptying duration 0.01 h (study-1 design; study 2
  #   used 4 and 9.5 h). The packaged model fixes the emptying times to 4 and
  #   10 h and the bolus duration to 0.01 h; users wanting to simulate a
  #   study-2 design can override ltet2 to log(9.5). The emptying-window
  #   gating uses an inclusive-low / exclusive-high indicator on absolute
  #   simulation time t (rxode2's lowercase time variable). For multi-dose
  #   simulation beyond 24 h the meal-time pattern is hard-coded relative to
  #   the first dose; the packaged model is best suited to single-dose
  #   simulation per the paper's design.
  # * EHC re-routing. Paper Methods state 'all MPAG secreted from GB to
  #   intestine were completely converted to MPA and was followed by
  #   reabsorption into the system'. The gallbladder pool empties into the
  #   GI depot, and absorption via the same first-order ka delivers the
  #   recycled amount to MPA central. No additional MPAG -> MPA conversion
  #   compartment is introduced; the conversion is captured by routing the
  #   gallbladder release through depot rather than directly to central_mpag.
  # * MPA elimination. Per identifiability assumption fm = 1 (paper Methods,
  #   p. 5: 'the conversion ratio from MPA to MPAG (fm) was fixed at 100%'),
  #   so paper's k20 = 0 and k24 = k20*; all MPA elimination from central is
  #   one-pass conversion to MPAG. CL_MPA/F = k20* * V_2/F is therefore the
  #   apparent rate of MPA-to-MPAG conversion. The corresponding flux enters
  #   central_mpag 1:1 (paper Methods 'completely converted', and fm = 1).
  # * Weight scaling. Paper Eq 5 ('slope linear without intercept') maps to
  #   a fixed allometric exponent of 1, encoded as e_wt_cl / e_wt_q /
  #   e_wt_vp = fixed(1) so that the source-trace is explicit. The
  #   reference weight is 65.5 kg (Table 1 median).
  # * Residual error. Paper reports an exponential residual model for both
  #   MPA and MPAG (Methods Eq 3: Y = IPRED * exp(epsilon)). This maps
  #   exactly to nlmixr2's lnorm() form with expSd on the log scale, where
  #   the reported '%' value is the log-scale residual SD.
  # * Concentration units. Bioanalytical reporting was in mg/L (Methods).
  #   The plasma concentrations Cc (MPA) and Cc_mpag (MPAG) inherit
  #   central / vc and central_mpag / vc_mpag in mg/L given dose in mg.
  # * Parameter names. Paper-named parameters not in the canonical PK
  #   register (ltet1, ltet2, ldur_gb, lk51, lehcp, psi_q_cl_mpag) are
  #   documented in the labels and source-trace comments; this follows the
  #   precedent of Ide_2009_pravastatin (ltg) for paper-specific EHC
  #   timing constants.
  ini({
    # ---- Absorption ----
    # Paper k_12 (absorption rate) and t_lag (Table 3).
    lka     <- log(3.53);    label("First-order absorption rate constant ka (1/h)")  # Table 3 Final Model k_12 = 3.53 1/h (RSE 12.4%)
    ltlag   <- log(0.0956);  label("Absorption lag time tlag (h)")                   # Table 3 Final Model t_lag = 0.0956 h (RSE 15.8%)

    # ---- MPA disposition (2-compartment central + peripheral1) ----
    lcl     <- log(10.2);    label("Apparent MPA clearance CL_MPA/F (L/h)")          # Table 3 Final Model CL_MPA/F = 10.2 L/h (RSE 5.7%)
    lvc     <- log(12.5);    label("Apparent MPA central volume V_2/F (L)")          # Table 3 Final Model V_2/F = 12.5 L (RSE 8.3%)
    lq      <- log(16.1);    label("Apparent MPA inter-compartmental clearance Q/F (L/h)")  # Table 3 Final Model Q/F = 16.1 L/h (RSE 5.1%)
    lvp     <- log(213);     label("Apparent MPA peripheral volume V_3/F (L)")       # Table 3 Final Model V_3/F = 213 L (RSE 9.1%)

    # ---- MPAG disposition (1-compartment central_mpag) ----
    lcl_mpag <- log(1.38);   label("Apparent MPAG renal clearance CL_MPAG/F (L/h)")  # Table 3 Final Model CL_MPAG/F = 1.38 L/h (RSE 6.9%)
    lvc_mpag <- log(4.40);   label("Apparent MPAG central volume V_4/F (L)")         # Table 3 Final Model V_4/F = 4.40 L (RSE 6.4%)

    # ---- Enterohepatic recirculation ----
    # EHCP = k45 / (k40 + k45) is the fraction of MPAG flux routed into
    # the gallbladder (vs. renal); the packaged model derives k45 from
    # EHCP and the renal-clearance branch inside model().
    lehcp    <- log(0.291);  label("EHCP: fraction of MPAG biliary-routed at the branch (unitless)")  # Table 3 Final Model EHCP = 29.1% (RSE 10.4%)
    lk51     <- log(67.5);   label("Gallbladder-to-GI rate constant k51 (1/h, active only in the bolus emptying window)")  # Table 3 Final Model k51 = 67.5 1/h (RSE 12.7%)

    # Gallbladder emptying timing (paper Methods, p. 5: 'GB emptying is
    # postulated to be at mealtime, i.e. at 4 and 10 h postdose in the
    # first study and 4 and 9.5 h postdose in the second study').
    # Packaged model fixes the study-1 schedule; users may override.
    ltet1    <- fixed(log(4));    label("First gallbladder emptying time ET1 (h post-dose, study-1 schedule)")    # Methods p. 5
    ltet2    <- fixed(log(10));   label("Second gallbladder emptying time ET2 (h post-dose, study-1 schedule)")   # Methods p. 5
    ldur_gb  <- fixed(log(0.01)); label("Gallbladder emptying bolus duration (h)")                                # Methods p. 5: 'assuming emptying duration was 0.01 h'

    # ---- Body-weight scaling on CL_MPA/F, Q/F, V_3/F ----
    # Paper Eq 5 ('slope linear without intercept') with WT_m = 65.5 kg.
    # Implemented as (WT/65.5)^exponent with exponent fixed at 1.
    e_wt_cl  <- fixed(1);     label("Allometric exponent on CL_MPA/F (unitless, fixed; paper Eq 5)")  # Methods Eq 5
    e_wt_q   <- fixed(1);     label("Allometric exponent on Q/F (unitless, fixed; paper Eq 5)")       # Methods Eq 5
    e_wt_vp  <- fixed(1);     label("Allometric exponent on V_3/F (unitless, fixed; paper Eq 5)")     # Methods Eq 5

    # ---- Cross-parameter eta linkage (paper 'q') ----
    # Paper Table 3 reports an estimated scalar q such that
    # eta(CL_MPAG/F) = q * eta(Q/F). Encoded as a paper-named parameter;
    # not in the canonical PK-parameter register (analogous precedent:
    # Ide 2009 'ltg').
    psi_q_cl_mpag <- 1.33;   label("Cross-parameter eta scalar: eta(CL_MPAG/F) = psi_q_cl_mpag * eta(Q/F) (paper q)")  # Table 3 Final Model q = 1.33 (RSE 27.2%)

    # ---- Inter-individual variability (exponential IIV) ----
    # Paper Methods Eq 1 (Pi = P' * exp(eta_i)) with paper-reported CV%.
    # Variance on the log scale via omega^2 = log(1 + CV^2).
    # Note: no separate eta on lcl_mpag is declared; the IIV on CL_MPAG/F
    # is conveyed via psi_q_cl_mpag * etalq (see model() body).
    etaltlag    ~ 0.281     # Table 3 IIV t_lag    = 57.3% CV -> log(1 + 0.573^2) = 0.281
    etalka      ~ 0.308     # Table 3 IIV k_12     = 60.3% CV -> log(1 + 0.603^2) = 0.308
    etalq       ~ 0.0186    # Table 3 IIV Q/F      = 13.7% CV -> log(1 + 0.137^2) = 0.0186
    etalcl      ~ 0.0351    # Table 3 IIV CL_MPA/F = 18.9% CV -> log(1 + 0.189^2) = 0.0351
    etalvc      ~ 0.1125    # Table 3 IIV V_2/F    = 34.5% CV -> log(1 + 0.345^2) = 0.1125
    etalvp      ~ 0.0502    # Table 3 IIV V_3/F    = 22.7% CV -> log(1 + 0.227^2) = 0.0502
    etalvc_mpag ~ 0.0520    # Table 3 IIV V_4/F    = 23.1% CV -> log(1 + 0.231^2) = 0.0520
    etalehcp    ~ 0.0807    # Table 3 IIV EHCP     = 29.0% CV -> log(1 + 0.290^2) = 0.0807

    # ---- Residual error (exponential / log-normal per Methods Eq 3) ----
    # Y = IPRED * exp(epsilon) with epsilon ~ N(0, expSd^2); the paper-
    # reported value is the log-scale residual SD.
    expSd        <- 0.453;   label("Log-normal residual SD for MPA Cc (log-scale SD)")              # Table 3 Final Model e for MPA  = 45.3% (RSE 9.3%)
    expSd_mpag   <- 0.208;   label("Log-normal residual SD for MPAG Cc_mpag (log-scale SD)")        # Table 3 Final Model e for MPAG = 20.8% (RSE 16.3%)
  })

  model({
    # ---- 1. Individual structural parameters ----
    ka       <- exp(lka       + etalka)
    tlag     <- exp(ltlag     + etaltlag)
    vc       <- exp(lvc       + etalvc)
    cl       <- exp(lcl       + etalcl)      * (WT / 65.5)^e_wt_cl
    q        <- exp(lq        + etalq)       * (WT / 65.5)^e_wt_q
    vp       <- exp(lvp       + etalvp)      * (WT / 65.5)^e_wt_vp
    vc_mpag  <- exp(lvc_mpag  + etalvc_mpag)
    # IIV on CL_MPAG/F is linked to etalq via the paper's 'q' scalar
    # (Table 3): eta(CL_MPAG/F) = psi_q_cl_mpag * eta(Q/F).
    cl_mpag  <- exp(lcl_mpag  + psi_q_cl_mpag * etalq)
    ehcp     <- exp(lehcp     + etalehcp)
    k51      <- exp(lk51)
    et1      <- exp(ltet1)
    et2      <- exp(ltet2)
    dur_gb   <- exp(ldur_gb)

    # ---- 2. Micro-constant derivations ----
    # MPA central elimination is one-pass conversion to MPAG (paper fm = 1
    # identifiability assumption -> k20 = 0, k24 = k20* = CL_MPA/F / V_2/F).
    k_central_to_mpag <- cl / vc                          # = k24 (paper)
    k12_central_to_p  <- q  / vc                          # = k23 (paper)
    k21_p_to_central  <- q  / vp                          # = k32 (paper)
    # MPAG branch: total elimination from central_mpag splits into a renal
    # arm (CL_MPAG/F at rate k40 = CL_MPAG/F / V_4/F) and a biliary arm
    # (rate k45 = k40 * EHCP / (1 - EHCP) from EHCP = k45 / (k40 + k45)).
    k40_mpag_renal    <- cl_mpag / vc_mpag                # = k40 (paper)
    k45_mpag_to_gb    <- k40_mpag_renal * ehcp / (1 - ehcp)  # = k45 (paper)

    # ---- 3. Gallbladder emptying window ----
    # Paper Methods: bolus emptying of duration dur_gb (= 0.01 h) at meal
    # times et1 and et2; rate k51 active only inside the window.
    # rxode2 exposes absolute simulation time as the lowercase variable t.
    in_emptying  <- ((t >= et1) * (t < (et1 + dur_gb))) +
                    ((t >= et2) * (t < (et2 + dur_gb)))
    k51_active   <- k51 * in_emptying

    # ---- 4. ODE system ----
    # State (amount, mg-equivalent):
    #   depot              -- GI lumen (oral input)
    #   central            -- MPA central
    #   peripheral1        -- MPA peripheral
    #   central_mpag       -- MPAG central
    #   gallbladder_mpag   -- gallbladder MPAG pool (released to depot during
    #                         the meal-time bolus windows and reabsorbed via ka)
    # The gallbladder->depot bolus encodes the paper's assumption that
    # 'all MPAG secreted from GB to intestine were completely converted to
    # MPA and was followed by reabsorption into the system' (Methods p. 4).
    d/dt(depot)            <- -ka * depot + k51_active * gallbladder_mpag
    d/dt(central)          <-  ka * depot - k_central_to_mpag * central -
                               k12_central_to_p * central +
                               k21_p_to_central * peripheral1
    d/dt(peripheral1)      <-  k12_central_to_p * central -
                               k21_p_to_central * peripheral1
    d/dt(central_mpag)     <-  k_central_to_mpag * central -
                               k40_mpag_renal * central_mpag -
                               k45_mpag_to_gb * central_mpag
    d/dt(gallbladder_mpag) <-  k45_mpag_to_gb * central_mpag -
                               k51_active * gallbladder_mpag

    # Absorption lag on the oral depot.
    alag(depot) <- tlag

    # ---- 5. Observation and residual error ----
    # Concentration in mg/L given dose in mg and volumes in L.
    Cc       <- central      / vc
    Cc_mpag  <- central_mpag / vc_mpag

    Cc       ~ lnorm(expSd)
    Cc_mpag  ~ lnorm(expSd_mpag)
  })
}
