AlvarezJimenez_2016_scopolamine <- function() {
  description <- "Two-compartment linear intravenous PK model for scopolamine hydrobromide in healthy adults (18-78 years), coupled with 13 indirect-response EMAX PD models describing neurophysiological, psychomotor, and cognitive tests: 0-back reaction time, saccadic inaccuracy, saccadic peak velocity, adaptive tracker performance, EEG alpha/delta/theta band power in Fz-Cz and Pz-Oz leads, and N-back 0/1/2 correct-answer percentages (log-odds). Age and body weight modify clearance; body weight modifies peripheral volume; age modifies kIN for 0-back RT and delta/theta EEG baselines and EC50 for saccadic peak velocity."
  reference <- paste(
    "Alvarez-Jimenez R, Groeneveld GJ, van Gerven JMA, Goulooze SC,",
    "Baakman AC, Hay JL, Stevens J.",
    "Model-based exposure-response analysis to quantify age related differences",
    "in the response to scopolamine in healthy subjects.",
    "Br J Clin Pharmacol. 2016;82(4):1011-1021. doi:10.1111/bcp.13031.",
    sep = " "
  )
  vignette <- "AlvarezJimenez_2016_scopolamine"
  units <- list(time = "min", dosing = "mg", concentration = "pg/mL")

  # PD response compartments -- turnover states for each of the 13 PD
  # endpoints (no canonical compartment role fits) plus the drug-reduction
  # compartments for the three N-back logit models.
  paper_specific_compartments <- c(
    "nbRT", "sacInacc", "sacPV", "adTrack",
    "eegAfc", "eegApo", "eegDfc", "eegDpo", "eegTfc", "eegTpo",
    "red_nb0", "red_nb1", "red_nb2"
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference weight 78.5 kg (cohort mean). Enters as (WT/78.5)^CWC on CL and (WT/78.5)^CWV on Vp, per Supplemental Equations A.1-A.2.",
      source_name        = "WGT"
    ),
    AGE = list(
      description        = "Subject age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference age 28 years. Enters as (AGE/28)^CAC on CL (Supplemental Equation A.1), exp(-AGE/CAE) on saccadic peak velocity EC50 (Eq B.1), and exp(AGE/CAB) on kIN for 0-back RT, delta Fz-Cz, delta Pz-Oz, and theta Fz-Cz EEG endpoints (Eq B.2).",
      source_name        = "AGE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 135L,
    n_studies      = 4L,
    age_range      = "18-78 years",
    age_median     = "39 years",
    weight_range   = "not reported explicitly; cohort mean 78.6 kg (SD 9.08)",
    weight_median  = "78.6 kg",
    sex_female_pct = 9.6,
    race_ethnicity = "not reported",
    disease_state  = "healthy volunteers",
    dose_range     = "0.3 mg (subjects >=65 years) or 0.5 mg (subjects <65 years) intravenous scopolamine hydrobromide infused over 15 min",
    regions        = "Netherlands (Centre for Human Drug Research, Leiden)",
    notes          = "Pooled analysis of four clinical studies. Three earlier studies enrolled only male subjects (drug teratogenicity constraints for parallel study arms). Study 4 (older-subject cohort, 65-78 y, n=36) enrolled 13 females (36% female in that cohort). See Table 1 of the source paper. Older subjects (>=65 y) received the 0.3 mg dose; younger subjects the 0.5 mg dose. Body weight and age summarised in Table 1."
  )

  ini({
    # ------------------------------------------------------------------
    # Pharmacokinetics -- 2-cmt linear model (Supplemental Eq A.1-A.5)
    # Reference age 28 y, reference weight 78.5 kg per Eq A.1-A.2.
    # ------------------------------------------------------------------
    lcl  <- log(1.09)     ; label("Clearance alpha (L/min)")                                # Table 2 CL[alpha] = 1.09
    lvc  <- log(2.66)     ; label("Central volume of distribution (L)")                     # Table 2 Vc = 2.66
    lvp  <- log(62.10)    ; label("Peripheral volume of distribution beta (L)")             # Table 2 Vp[beta] = 62.10
    lq   <- log(1.01)     ; label("Inter-compartmental clearance (L/min)")                  # Table 2 Q = 1.01

    # Covariate coefficients (CAC, CWC on CL; CWV on Vp; supplemental Eq A.1-A.2)
    cac  <- -0.12         ; label("Age exponent on CL (unitless, base of AGE/28)")          # Table 2 CAC = -0.12
    cwc  <- 0.56          ; label("Weight exponent on CL (unitless, base of WT/78.5)")      # Table 2 CWC = 0.56
    cwv  <- 0.38          ; label("Weight exponent on Vp (unitless, base of WT/78.5)")      # Table 2 CWV = 0.38

    # PK IIV -- diagonal; paper reports no PK omega block. CV converted via
    # variance = log(1 + CV^2).
    etalcl ~ 0.0106       ; label("IIV var on lcl (~10.3% CV)")                             # Table 2 IIV(CL) = 10.3%
    etalvc ~ 0.437        ; label("IIV var on lvc (~74.1% CV, shrinkage 49.3%)")            # Table 2 IIV(Vc) = 74.1%
    etalvp ~ 0.00899      ; label("IIV var on lvp (~9.5% CV)")                              # Table 2 IIV(Vp) = 9.5%

    # PK residual (proportional; sigma^2 = 0.045 -> propSd = sqrt(0.045))
    propSd <- 0.21213     ; label("Proportional residual error on Cc (fraction)")           # Table 2 Error sigma^2 = 0.045

    # ------------------------------------------------------------------
    # PD endpoints -- indirect-response EMAX (Supplemental Eq B.1-B.4).
    # E0 = kIN/kOUT (steady-state baseline); drug modulates production via
    # (1 + effect) for endpoints where scopolamine INCREASES the response
    # (0-back RT, saccadic inaccuracy, delta/theta EEG bands), and
    # (1 - effect) for endpoints where scopolamine DECREASES the response
    # (saccadic peak velocity, adaptive tracker, alpha EEG bands, N-back
    # correct answers).
    # ------------------------------------------------------------------

    # ---------- 0-back reaction time (msec) ----------
    lec50_nbRT  <- log(1300)   ; label("EC50 for 0-back RT (pg/mL)")                        # Table 2 EC50 = 1300
    lemax_nbRT  <- log(0.837)  ; label("EMAX for 0-back RT (fraction)")                     # Table 2 EMAX = 0.837
    lkin_nbRT   <- log(3.42)   ; label("kIN reference for 0-back RT (msec/min)")            # Table 2 kIN[lambda] = 3.42
    lkout_nbRT  <- log(0.00974); label("kOUT for 0-back RT (1/min)")                        # Table 2 kOUT = 0.00974
    cab_nbRT    <- 229         ; label("Age scale (CAB) on kIN for 0-back RT (years)")      # Table 2 CAB = 229
    etalec50_nbRT ~ 0.632      ; label("IIV var on lec50_nbRT (~99.8% CV, shrinkage 33.2%)")# Table 2 IIV(EC50) = 99.8%
    etalkin_nbRT  ~ 0.0106     ; label("IIV var on lkin_nbRT (~10.5% CV)")                  # Table 2 IIV(kIN) = 10.5%
    propSd_nbRT   <- 0.09995   ; label("Proportional residual error for 0-back RT (fraction)") # Table 2 sigma^2 = 0.00999

    # ---------- Saccadic inaccuracy (%) ----------
    lec50_sacInacc <- log(55.1)   ; label("EC50 for saccadic inaccuracy (pg/mL)")           # Table 2 EC50 = 55.1
    lemax_sacInacc <- log(0.388)  ; label("EMAX for saccadic inaccuracy (fraction)")        # Table 2 EMAX = 0.388
    lkin_sacInacc  <- log(22.8)   ; label("kIN for saccadic inaccuracy (%/min)")            # Table 2 kIN = 22.8
    lkout_sacInacc <- log(3.33)   ; label("kOUT for saccadic inaccuracy (1/min)")           # Table 2 kOUT = 3.33
    # Omega block between EC50 and EMAX (OC (omega correlation) = 0.06
    # interpreted as covariance;
    # var_ec50 = log(1+0.137^2) = 0.0186, var_emax = log(1+0.499^2) = 0.2216).
    etalec50_sacInacc + etalemax_sacInacc ~ c(0.0186, 0.06, 0.2216)                          # Table 2 IIV 13.7% / 49.9%, OC (omega correlation) = 0.06
    etalkin_sacInacc  ~ 0.0175     ; label("IIV var on lkin_sacInacc (~13.2% CV)")           # Table 2 IIV(kIN) = 13.2%
    propSd_sacInacc   <- 0.22204   ; label("Proportional residual error for saccadic inaccuracy") # Table 2 sigma^2 = 0.0493

    # ---------- Saccadic peak velocity (deg/sec) ----------
    lec50_sacPV  <- log(2530)  ; label("EC50 reference (epsilon) for saccadic peak velocity (pg/mL)") # Table 2 EC50 = 2530
    lemax_sacPV  <- log(0.232) ; label("EMAX for saccadic peak velocity (fraction)")         # Table 2 EMAX = 0.232
    lkin_sacPV   <- log(1280)  ; label("kIN for saccadic peak velocity (deg/sec/min)")       # Table 2 kIN = 1280
    lkout_sacPV  <- log(2.73)  ; label("kOUT for saccadic peak velocity (1/min)")            # Table 2 kOUT = 2.73
    cae_sacPV    <- 34.0       ; label("Age scale (CAE) on EC50 for saccadic peak velocity (years)") # Table 2 CAE = 34.0
    etalec50_sacPV ~ 1.322     ; label("IIV var on lec50_sacPV (~277.0% CV)")                # Table 2 IIV(EC50) = 277.0%
    etalkin_sacPV  ~ 0.00975   ; label("IIV var on lkin_sacPV (~9.9% CV)")                   # Table 2 IIV(kIN) = 9.9%
    # Paper labels residual as (prop) with sigma^2 = 768.0 which is physically
    # implausible for a proportional error (2770% CV). Interpreted here as
    # additive with SD = sqrt(768) = 27.71 deg/sec (~6% of typical peak velocity
    # ~470 deg/sec). See vignette Errata for this deviation.
    addSd_sacPV  <- 27.71      ; label("Additive residual error for saccadic peak velocity (deg/sec)") # Table 2 sigma^2 = 768.0

    # ---------- Adaptive tracker performance (%) ----------
    lec50_adTrack <- log(386)    ; label("EC50 for adaptive tracker (pg/mL)")                # Table 2 EC50 = 386
    lgamma_adTrack <- log(1.10)  ; label("Hill gamma for adaptive tracker (unitless)")       # Table 2 gamma = 1.10
    lkin_adTrack   <- log(0.636) ; label("kIN for adaptive tracker (%/min)")                 # Table 2 kIN = 0.636
    lkout_adTrack  <- log(0.0294); label("kOUT for adaptive tracker (1/min)")                # Table 2 kOUT = 0.0294
    # EMAX fixed to 1 per paper (see Model development section).
    etalec50_adTrack ~ 0.556   ; label("IIV var on lec50_adTrack (~85.5% CV)")               # Table 2 IIV(EC50) = 85.5%
    etalkin_adTrack  ~ 0.0491  ; label("IIV var on lkin_adTrack (~22.4% CV)")                # Table 2 IIV(kIN) = 22.4%
    addSd_adTrack    <- 2.939  ; label("Additive residual error for adaptive tracker (%)")   # Table 2 sigma^2 = 8.64

    # ---------- EEG alpha Fz-Cz (uV) ----------
    lec50_eegAfc <- log(2.55)     ; label("EC50 for EEG alpha Fz-Cz (pg/mL)")                # Table 2 EC50 = 2.55
    lemax_eegAfc <- log(0.232)    ; label("EMAX for EEG alpha Fz-Cz (fraction)")             # Table 2 EMAX = 0.232
    lkin_eegAfc  <- log(0.206)    ; label("kIN for EEG alpha Fz-Cz (uV/min)")                # Table 2 kIN = 0.206
    lkout_eegAfc <- log(0.0699)   ; label("kOUT for EEG alpha Fz-Cz (1/min)")                # Table 2 kOUT = 0.0699
    etalec50_eegAfc ~ 2.882      ; label("IIV var on lec50_eegAfc (~1530.0% CV)")            # Table 2 IIV(EC50) = 1530.0%
    etalemax_eegAfc ~ 0.652      ; label("IIV var on lemax_eegAfc (~95.4% CV)")              # Table 2 IIV(EMAX) = 95.4%
    etalkin_eegAfc  ~ 0.190      ; label("IIV var on lkin_eegAfc (~46.6% CV)")               # Table 2 IIV(kIN) = 46.6%
    propSd_eegAfc   <- 0.160     ; label("Proportional residual error for EEG alpha Fz-Cz")  # Table 2 sigma^2 = 0.0256

    # ---------- EEG alpha Pz-Oz (uV) ----------
    lec50_eegApo <- log(14.7)     ; label("EC50 for EEG alpha Pz-Oz (pg/mL)")                # Table 2 EC50 = 14.7
    lemax_eegApo <- log(0.443)    ; label("EMAX for EEG alpha Pz-Oz (fraction)")             # Table 2 EMAX = 0.443
    lkin_eegApo  <- log(3.05)     ; label("kIN for EEG alpha Pz-Oz (uV/min)")                # Table 2 kIN = 3.05
    lkout_eegApo <- log(0.553)    ; label("kOUT for EEG alpha Pz-Oz (1/min)")                # Table 2 kOUT = 0.553
    etalec50_eegApo ~ 2.647      ; label("IIV var on lec50_eegApo (~1309.3% CV)")            # Table 2 IIV(EC50) = 1309.3%
    etalemax_eegApo ~ 0.129      ; label("IIV var on lemax_eegApo (~37.6% CV)")              # Table 2 IIV(EMAX) = 37.6%
    etalkin_eegApo  ~ 0.297      ; label("IIV var on lkin_eegApo (~57.9% CV)")               # Table 2 IIV(kIN) = 57.9%
    propSd_eegApo   <- 0.207     ; label("Proportional residual error for EEG alpha Pz-Oz")  # Table 2 sigma^2 = 0.0429

    # ---------- EEG delta Fz-Cz (uV) ----------
    lec50_eegDfc <- log(469)      ; label("EC50 for EEG delta Fz-Cz (pg/mL)")                # Table 2 EC50 = 469
    lemax_eegDfc <- log(0.419)    ; label("EMAX for EEG delta Fz-Cz (fraction)")             # Table 2 EMAX = 0.419
    lkin_eegDfc  <- log(0.0354)   ; label("kIN reference (lambda) for EEG delta Fz-Cz (uV/min)") # Table 2 kIN = 0.0354
    lkout_eegDfc <- log(0.0166)   ; label("kOUT for EEG delta Fz-Cz (1/min)")                # Table 2 kOUT = 0.0166
    cab_eegDfc   <- 161           ; label("Age scale (CAB) on kIN for EEG delta Fz-Cz (years)") # Table 2 CAB = 161
    # Two omega blocks reported (EC50 vs EMAX; EMAX vs kIN). Interpreted as
    # covariances per paper caption's imprecise labeling of OC (omega
    # correlation). Encoded as a
    # 3-parameter block; the EC50-vs-kIN off-diagonal (not reported) is set
    # to zero via fixed().
    etalec50_eegDfc + etalemax_eegDfc + etalkin_eegDfc ~
      c(1.209,
        1.410, 2.418,
        fixed(0), 0.153, 0.0800)                                                             # Table 2 IIV 154.5%/325.4%/28.9%, OC (omega correlation) = 1.410 (EC50 vs EMAX), 0.153 (EMAX vs kIN)
    propSd_eegDfc   <- 0.184     ; label("Proportional residual error for EEG delta Fz-Cz")  # Table 2 sigma^2 = 0.0338

    # ---------- EEG delta Pz-Oz (uV) ----------
    lec50_eegDpo <- log(1230)     ; label("EC50 for EEG delta Pz-Oz (pg/mL)")                # Table 2 EC50 = 1230
    lemax_eegDpo <- log(0.537)    ; label("EMAX for EEG delta Pz-Oz (fraction)")             # Table 2 EMAX = 0.537
    lkin_eegDpo  <- log(0.726)    ; label("kIN reference (lambda) for EEG delta Pz-Oz (uV/min)") # Table 2 kIN = 0.726
    lkout_eegDpo <- log(0.348)    ; label("kOUT for EEG delta Pz-Oz (1/min)")                # Table 2 kOUT = 0.348
    cab_eegDpo   <- 210           ; label("Age scale (CAB) on kIN for EEG delta Pz-Oz (years)") # Table 2 CAB = 210
    etalec50_eegDpo ~ 2.541      ; label("IIV var on lec50_eegDpo (~1166.2% CV)")            # Table 2 IIV(EC50) = 1166.2%
    etalemax_eegDpo ~ 0.264      ; label("IIV var on lemax_eegDpo (~55.8% CV, shrinkage 48.0%)") # Table 2 IIV(EMAX) = 55.8%
    etalkin_eegDpo  ~ 0.0862     ; label("IIV var on lkin_eegDpo (~30.3% CV)")               # Table 2 IIV(kIN) = 30.3%
    propSd_eegDpo   <- 0.174     ; label("Proportional residual error for EEG delta Pz-Oz")  # Table 2 sigma^2 = 0.0304

    # ---------- EEG theta Fz-Cz (uV) ----------
    lec50_eegTfc <- log(4110)     ; label("EC50 for EEG theta Fz-Cz (pg/mL)")                # Table 2 EC50 = 4110
    lemax_eegTfc <- log(0.594)    ; label("EMAX for EEG theta Fz-Cz (fraction)")             # Table 2 EMAX = 0.594
    lkin_eegTfc  <- log(0.863)    ; label("kIN reference (lambda) for EEG theta Fz-Cz (uV/min)") # Table 2 kIN = 0.863
    lkout_eegTfc <- log(0.594)    ; label("kOUT for EEG theta Fz-Cz (1/min)")                # Table 2 kOUT = 0.594
    cab_eegTfc   <- 142           ; label("Age scale (CAB) on kIN for EEG theta Fz-Cz (years)") # Table 2 CAB = 142
    etalec50_eegTfc ~ 3.166      ; label("IIV var on lec50_eegTfc (~2162.7% CV)")            # Table 2 IIV(EC50) = 2162.7%
    etalkin_eegTfc  ~ 0.0798     ; label("IIV var on lkin_eegTfc (~28.4% CV)")               # Table 2 IIV(kIN) = 28.4%
    propSd_eegTfc   <- 0.172     ; label("Proportional residual error for EEG theta Fz-Cz")  # Table 2 sigma^2 = 0.0297

    # ---------- EEG theta Pz-Oz (uV) ----------
    lec50_eegTpo <- log(8240)     ; label("EC50 for EEG theta Pz-Oz (pg/mL)")                # Table 2 EC50 = 8240
    lemax_eegTpo <- log(1.15)     ; label("EMAX for EEG theta Pz-Oz (fraction)")             # Table 2 EMAX = 1.15
    lkin_eegTpo  <- log(13)       ; label("kIN for EEG theta Pz-Oz (uV/min)")                # Table 2 kIN = 13
    lkout_eegTpo <- log(6.08)     ; label("kOUT for EEG theta Pz-Oz (1/min)")                # Table 2 kOUT = 6.08
    etalec50_eegTpo ~ 2.987      ; label("IIV var on lec50_eegTpo (~1778.6% CV)")            # Table 2 IIV(EC50) = 1778.6%
    etalkin_eegTpo  ~ 0.107      ; label("IIV var on lkin_eegTpo (~33.4% CV)")               # Table 2 IIV(kIN) = 33.4%
    propSd_eegTpo   <- 0.195     ; label("Proportional residual error for EEG theta Pz-Oz")  # Table 2 sigma^2 = 0.0379

    # ---------- N-back 0 (% correct, log-odds observation) ----------
    # Modelled as direct effect on log-odds baseline. Paper's four buffer
    # compartments for the delay are approximated by a single indirect-
    # response state with rate kIN_nb0. EMAX fixed to 1 per paper narrative.
    lec50_nb0 <- log(2140)        ; label("EC50 for 0-back correct-answer ratio (pg/mL)")     # Table 2 EC50 = 2140
    lemax_nb0 <- fixed(log(1))    ; label("EMAX for 0-back correct-answer ratio (fixed at 1)") # Table 2 EMAX = 1 (fixed per paper)
    lkin_nb0  <- log(0.0354)      ; label("kIN for 0-back ratio (1/min); approximated single-compartment for 4-transit chain") # Table 2 kIN = 0.0354
    lkout_nb0 <- log(0.036)       ; label("kOUT for 0-back ratio (1/min)")                    # Table 2 kOUT = 0.036
    e0_nb0    <- 3.77             ; label("Baseline log-odds of 0-back correct answers")      # Table 2 E0 = 3.77 (expit(3.77) = 97.7%)
    # Omega block: EMAX and E0 correlation (OC (omega correlation) = 0.40
    # interpreted as covariance).
    # etalemax_nb0 (log-scale IIV on EMAX; EMAX fixed at 1 so IIV represents
    # subject-specific scale factor exp(eta) applied to 1) and etae0_nb0
    # (additive IIV on log-odds baseline).
    etalemax_nb0 + etae0_nb0 ~ c(3.760, 0.40, 0.132)                                          # Table 2 IIV 631.3% / 37.8%, OC (omega correlation) = 0.40
    addSd_nb0 <- 0.0624           ; label("Additive residual error for 0-back log-odds (log-odds units)") # Table 2 sigma^2 = 0.0039

    # ---------- N-back 1 (% correct, log-odds observation) ----------
    lec50_nb1 <- log(0.624)       ; label("EC50 for 1-back correct-answer ratio (pg/mL)")     # Table 2 EC50 = 0.624
    lemax_nb1 <- log(0.961)       ; label("EMAX for 1-back correct-answer ratio (fraction)")  # Table 2 EMAX = 0.961
    lkin_nb1  <- log(1.10)        ; label("kIN for 1-back ratio (1/min)")                     # Table 2 kIN = 1.10
    lkout_nb1 <- log(221)         ; label("kOUT for 1-back ratio (1/min)")                    # Table 2 kOUT = 221
    e0_nb1    <- 3.33             ; label("Baseline log-odds of 1-back correct answers")      # Table 2 E0 = 3.33 (expit(3.33) = 96.5%)
    etalkin_nb1  ~ 3.535         ; label("IIV var on lkin_nb1 (~578.5% CV)")                 # Table 2 IIV(kIN) = 578.5%
    etalemax_nb1 ~ 0.898         ; label("IIV var on lemax_nb1 (~121.6% CV)")                # Table 2 IIV(EMAX) = 121.6%
    addSd_nb1 <- 0.0768           ; label("Additive residual error for 1-back log-odds (log-odds units)") # Table 2 sigma^2 = 0.0059

    # ---------- N-back 2 (% correct, log-odds observation) ----------
    lec50_nb2 <- log(212)         ; label("EC50 for 2-back correct-answer ratio (pg/mL)")     # Table 2 EC50 = 212
    lemax_nb2 <- log(1.09)        ; label("EMAX for 2-back correct-answer ratio (fraction)")  # Table 2 EMAX = 1.09
    lkin_nb2  <- log(0.54)        ; label("kIN for 2-back ratio (1/min)")                     # Table 2 kIN = 0.54
    lkout_nb2 <- log(0.317)       ; label("kOUT for 2-back ratio (1/min)")                    # Table 2 kOUT = 0.317
    e0_nb2    <- 2.84             ; label("Baseline log-odds of 2-back correct answers")      # Table 2 E0 = 2.84 (expit(2.84) = 94.5%)
    # Omega block: EMAX and E0 (OC = 0.23 interpreted as covariance).
    etalemax_nb2 + etae0_nb2 ~ c(0.868, 0.23, 0.115)                                          # Table 2 IIV 109.8% / 35.9%, OC (omega correlation) = 0.23
    addSd_nb2 <- 0.0922           ; label("Additive residual error for 2-back log-odds (log-odds units)") # Table 2 sigma^2 = 0.0085
  })

  model({
    # ---------------- PK ----------------
    # Individual PK parameters. Age and weight scale CL and Vp per Eq A.1-A.2.
    cl <- exp(lcl + etalcl) * (AGE / 28)^cac * (WT / 78.5)^cwc
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp) * (WT / 78.5)^cwv
    q  <- exp(lq)

    # Plasma scopolamine amount / concentration. linCmt() returns amount / vc
    # (concentration) in the same units as dose / vc; dose is in mg and vc in
    # L, so the raw linCmt() value is mg/L = ug/mL. Convert to pg/mL to match
    # the published EC50 units (Table 2). See vignette Errata for the unit-
    # conversion rationale.
    Cc <- 1e6 * linCmt()

    # ---------------- PD helpers ----------------
    # Individual PD parameters and drug-effect intensities for each endpoint.

    # 0-back RT (drug increases response)
    kin_nbRT_i  <- exp(lkin_nbRT + etalkin_nbRT) * exp(AGE / cab_nbRT)
    kout_nbRT_i <- exp(lkout_nbRT)
    ec50_nbRT_i <- exp(lec50_nbRT + etalec50_nbRT)
    emax_nbRT_i <- exp(lemax_nbRT)
    eff_nbRT    <- emax_nbRT_i * Cc / (ec50_nbRT_i + Cc)
    nbRT(0) <- kin_nbRT_i / kout_nbRT_i
    d/dt(nbRT) <- kin_nbRT_i * (1 + eff_nbRT) - kout_nbRT_i * nbRT

    # Saccadic inaccuracy (drug increases response)
    kin_sacInacc_i  <- exp(lkin_sacInacc + etalkin_sacInacc)
    kout_sacInacc_i <- exp(lkout_sacInacc)
    ec50_sacInacc_i <- exp(lec50_sacInacc + etalec50_sacInacc)
    emax_sacInacc_i <- exp(lemax_sacInacc + etalemax_sacInacc)
    eff_sacInacc    <- emax_sacInacc_i * Cc / (ec50_sacInacc_i + Cc)
    sacInacc(0) <- kin_sacInacc_i / kout_sacInacc_i
    d/dt(sacInacc) <- kin_sacInacc_i * (1 + eff_sacInacc) - kout_sacInacc_i * sacInacc

    # Saccadic peak velocity (drug decreases response; age modifies EC50 per Eq B.1)
    kin_sacPV_i  <- exp(lkin_sacPV + etalkin_sacPV)
    kout_sacPV_i <- exp(lkout_sacPV)
    ec50_sacPV_i <- exp(lec50_sacPV + etalec50_sacPV) * exp(-AGE / cae_sacPV)
    emax_sacPV_i <- exp(lemax_sacPV)
    eff_sacPV    <- emax_sacPV_i * Cc / (ec50_sacPV_i + Cc)
    sacPV(0) <- kin_sacPV_i / kout_sacPV_i
    d/dt(sacPV) <- kin_sacPV_i * (1 - eff_sacPV) - kout_sacPV_i * sacPV

    # Adaptive tracker (drug decreases response; EMAX fixed at 1; Hill gamma estimated)
    kin_adTrack_i  <- exp(lkin_adTrack + etalkin_adTrack)
    kout_adTrack_i <- exp(lkout_adTrack)
    ec50_adTrack_i <- exp(lec50_adTrack + etalec50_adTrack)
    gamma_adTrack_i <- exp(lgamma_adTrack)
    eff_adTrack    <- Cc^gamma_adTrack_i / (ec50_adTrack_i^gamma_adTrack_i + Cc^gamma_adTrack_i)
    adTrack(0) <- kin_adTrack_i / kout_adTrack_i
    d/dt(adTrack) <- kin_adTrack_i * (1 - eff_adTrack) - kout_adTrack_i * adTrack

    # EEG alpha Fz-Cz (drug decreases response)
    kin_eegAfc_i  <- exp(lkin_eegAfc + etalkin_eegAfc)
    kout_eegAfc_i <- exp(lkout_eegAfc)
    ec50_eegAfc_i <- exp(lec50_eegAfc + etalec50_eegAfc)
    emax_eegAfc_i <- exp(lemax_eegAfc + etalemax_eegAfc)
    eff_eegAfc    <- emax_eegAfc_i * Cc / (ec50_eegAfc_i + Cc)
    eegAfc(0) <- kin_eegAfc_i / kout_eegAfc_i
    d/dt(eegAfc) <- kin_eegAfc_i * (1 - eff_eegAfc) - kout_eegAfc_i * eegAfc

    # EEG alpha Pz-Oz (drug decreases response)
    kin_eegApo_i  <- exp(lkin_eegApo + etalkin_eegApo)
    kout_eegApo_i <- exp(lkout_eegApo)
    ec50_eegApo_i <- exp(lec50_eegApo + etalec50_eegApo)
    emax_eegApo_i <- exp(lemax_eegApo + etalemax_eegApo)
    eff_eegApo    <- emax_eegApo_i * Cc / (ec50_eegApo_i + Cc)
    eegApo(0) <- kin_eegApo_i / kout_eegApo_i
    d/dt(eegApo) <- kin_eegApo_i * (1 - eff_eegApo) - kout_eegApo_i * eegApo

    # EEG delta Fz-Cz (drug increases response; age modifies kIN per Eq B.2)
    kin_eegDfc_i  <- exp(lkin_eegDfc + etalkin_eegDfc) * exp(AGE / cab_eegDfc)
    kout_eegDfc_i <- exp(lkout_eegDfc)
    ec50_eegDfc_i <- exp(lec50_eegDfc + etalec50_eegDfc)
    emax_eegDfc_i <- exp(lemax_eegDfc + etalemax_eegDfc)
    eff_eegDfc    <- emax_eegDfc_i * Cc / (ec50_eegDfc_i + Cc)
    eegDfc(0) <- kin_eegDfc_i / kout_eegDfc_i
    d/dt(eegDfc) <- kin_eegDfc_i * (1 + eff_eegDfc) - kout_eegDfc_i * eegDfc

    # EEG delta Pz-Oz (drug increases response; age modifies kIN per Eq B.2)
    kin_eegDpo_i  <- exp(lkin_eegDpo + etalkin_eegDpo) * exp(AGE / cab_eegDpo)
    kout_eegDpo_i <- exp(lkout_eegDpo)
    ec50_eegDpo_i <- exp(lec50_eegDpo + etalec50_eegDpo)
    emax_eegDpo_i <- exp(lemax_eegDpo + etalemax_eegDpo)
    eff_eegDpo    <- emax_eegDpo_i * Cc / (ec50_eegDpo_i + Cc)
    eegDpo(0) <- kin_eegDpo_i / kout_eegDpo_i
    d/dt(eegDpo) <- kin_eegDpo_i * (1 + eff_eegDpo) - kout_eegDpo_i * eegDpo

    # EEG theta Fz-Cz (drug increases response; age modifies kIN per Eq B.2)
    kin_eegTfc_i  <- exp(lkin_eegTfc + etalkin_eegTfc) * exp(AGE / cab_eegTfc)
    kout_eegTfc_i <- exp(lkout_eegTfc)
    ec50_eegTfc_i <- exp(lec50_eegTfc + etalec50_eegTfc)
    emax_eegTfc_i <- exp(lemax_eegTfc)
    eff_eegTfc    <- emax_eegTfc_i * Cc / (ec50_eegTfc_i + Cc)
    eegTfc(0) <- kin_eegTfc_i / kout_eegTfc_i
    d/dt(eegTfc) <- kin_eegTfc_i * (1 + eff_eegTfc) - kout_eegTfc_i * eegTfc

    # EEG theta Pz-Oz (drug increases response)
    kin_eegTpo_i  <- exp(lkin_eegTpo + etalkin_eegTpo)
    kout_eegTpo_i <- exp(lkout_eegTpo)
    ec50_eegTpo_i <- exp(lec50_eegTpo + etalec50_eegTpo)
    emax_eegTpo_i <- exp(lemax_eegTpo)
    eff_eegTpo    <- emax_eegTpo_i * Cc / (ec50_eegTpo_i + Cc)
    eegTpo(0) <- kin_eegTpo_i / kout_eegTpo_i
    d/dt(eegTpo) <- kin_eegTpo_i * (1 + eff_eegTpo) - kout_eegTpo_i * eegTpo

    # N-back 0 (log-odds observation; direct-effect approximation to paper's
    # 4-transit-compartment structure). Response = E0 - drug_reduction, where
    # drug_reduction follows first-order kinetics with rate kout_nb0. EMAX
    # fixed at 1 with a subject-specific scale via etalemax_nb0.
    kin_nb0_i  <- exp(lkin_nb0)
    kout_nb0_i <- exp(lkout_nb0)
    ec50_nb0_i <- exp(lec50_nb0)
    emax_nb0_i <- exp(lemax_nb0 + etalemax_nb0)
    e0_nb0_i   <- e0_nb0 + etae0_nb0
    eff_nb0    <- emax_nb0_i * Cc / (ec50_nb0_i + Cc)
    d/dt(red_nb0) <- kin_nb0_i * eff_nb0 - kout_nb0_i * red_nb0
    nb0 <- e0_nb0_i - red_nb0

    # N-back 1
    kin_nb1_i  <- exp(lkin_nb1 + etalkin_nb1)
    kout_nb1_i <- exp(lkout_nb1)
    ec50_nb1_i <- exp(lec50_nb1)
    emax_nb1_i <- exp(lemax_nb1 + etalemax_nb1)
    eff_nb1    <- emax_nb1_i * Cc / (ec50_nb1_i + Cc)
    d/dt(red_nb1) <- kin_nb1_i * eff_nb1 - kout_nb1_i * red_nb1
    nb1 <- e0_nb1 - red_nb1

    # N-back 2
    kin_nb2_i  <- exp(lkin_nb2)
    kout_nb2_i <- exp(lkout_nb2)
    ec50_nb2_i <- exp(lec50_nb2)
    emax_nb2_i <- exp(lemax_nb2 + etalemax_nb2)
    e0_nb2_i   <- e0_nb2 + etae0_nb2
    eff_nb2    <- emax_nb2_i * Cc / (ec50_nb2_i + Cc)
    d/dt(red_nb2) <- kin_nb2_i * eff_nb2 - kout_nb2_i * red_nb2
    nb2 <- e0_nb2_i - red_nb2

    # ---------------- Observations ----------------
    # Multi-output residual-error assignments. Simulation callers can
    # request every output by observing on time only (no cmt column);
    # for fitting, source data must carry a dvid= column with numeric
    # indices matching the observation order below (Cc = 1, nbRT = 2,
    # sacInacc = 3, ...).
    Cc            ~ prop(propSd)
    nbRT     ~ prop(propSd_nbRT)
    sacInacc ~ prop(propSd_sacInacc)
    sacPV    ~ add(addSd_sacPV)
    adTrack  ~ add(addSd_adTrack)
    eegAfc   ~ prop(propSd_eegAfc)
    eegApo   ~ prop(propSd_eegApo)
    eegDfc   ~ prop(propSd_eegDfc)
    eegDpo   ~ prop(propSd_eegDpo)
    eegTfc   ~ prop(propSd_eegTfc)
    eegTpo   ~ prop(propSd_eegTpo)
    nb0   ~ add(addSd_nb0)
    nb1   ~ add(addSd_nb1)
    nb2   ~ add(addSd_nb2)
  })
}
