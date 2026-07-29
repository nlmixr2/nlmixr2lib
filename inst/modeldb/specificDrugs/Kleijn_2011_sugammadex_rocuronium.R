Kleijn_2011_sugammadex_rocuronium <- function() {
  description <- paste(
    "Integrated population PK-PD model for sugammadex-mediated reversal",
    "of rocuronium-induced neuromuscular blockade (Kleijn 2011). Both",
    "sugammadex and rocuronium have two-compartment PK from IV bolus",
    "dosing into the central compartment; the sugammadex-rocuronium",
    "inclusion complex has its own two-compartment PK with parameters",
    "set equal to free sugammadex (Bom 2002 framework). Complex",
    "formation is dynamic with fixed equilibrium dissociation constant",
    "kd = 0.0559 uM and estimated dissociation rate k2 = 0.034 1/min",
    "(association k1 = k2 / kd = 0.61 1/(min*uM)). The rocuronium",
    "central concentration drives an effect compartment via ke0 =",
    "0.134 1/min; neuromuscular blockade (T4/T1 twitch ratio x 100)",
    "follows a sigmoid Emax form with Emax set equal to E0 so the",
    "readout decreases monotonically from baseline E0 ~ 104 toward 0",
    "as the effect-compartment rocuronium concentration rises.",
    "Sugammadex-mediated reversal enters as an additional first-order",
    "elimination of rocuronium from the effect compartment driven by",
    "the central free-sugammadex concentration (ks = 0.033 1/(min*uM)).",
    "Both plasma assays measured total drug (free + complex), so the",
    "Cc and Cc_roc outputs return total sugammadex and total rocuronium.",
    "Allometric scaling on all volumes (exponent 1), flows (0.75), and",
    "rate constants (-0.25) at reference WT = 70 kg; sugammadex CL is",
    "NOT allometrically scaled (creatinine-clearance covariate replaces",
    "size scaling). All units are molar inside the model (dose in",
    "umol, concentrations in uM); see vignette for mg-to-umol",
    "conversion with rocuronium MW = 529.78 g/mol and sugammadex MW =",
    "2178.01 g/mol (octasodium salt)."
  )
  reference <- paste(
    "Kleijn HJ, Zollinger DP, van den Heuvel MW, Kerbusch T.",
    "Population pharmacokinetic-pharmacodynamic analysis for",
    "sugammadex-mediated reversal of rocuronium-induced neuromuscular",
    "blockade. Br J Clin Pharmacol. 2011 Sep;72(3):415-433.",
    "doi:10.1111/j.1365-2125.2011.04000.x. PMID: 21501216.",
    "Complex-formation framework from Bom AHJ, Bradley M, Cameron K,",
    "Clark JK, Van Egmond J, Feilden H, MacLean EJ, Muir AW, Palin R,",
    "Rees DC, Zhang MQ. A novel concept of reversing neuromuscular",
    "block: chemical encapsulation of rocuronium bromide by a",
    "cyclodextrin-based synthetic host. Angew Chem Int Ed Engl. 2002.",
    sep = " "
  )
  vignette <- "Kleijn_2011_sugammadex_rocuronium"
  units <- list(
    time          = "min",
    dosing        = "umol",
    concentration = "umol/L"
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric reference 70 kg applied to all volumes (exponent 1),",
        "all flows (exponent 0.75) and all rate constants (exponent -0.25)",
        "for the rocuronium PK, the complex PK, the rocuronium effect",
        "compartment equilibration rate ke0, and the sugammadex-mediated",
        "effect-compartment elimination rate ks. Sugammadex CL is NOT",
        "allometrically scaled because the paper-fitted creatinine-clearance",
        "(REN) effect on CL absorbs the size effect (Methods 'Renal",
        "function'). Sugammadex CL and V1 carry an additional linear",
        "BW shift centred at the dataset population mean 74.5 kg, on top",
        "of the allometric scaling for V1 (Table 4 'CL_BW = 1 + theta *",
        "(BW - 74.5)' and 'V1_BW = 1 + theta * (BW - 74.5)').",
        sep = " "
      ),
      source_name        = "BW"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Affects rocuronium CL (linear shift centred at 43 years; Kleijn",
        "2011 Table 4 'CL_AGE = 1 + theta * (AGE - 43.0)') and rocuronium",
        "V2 (exponential shift centred at 43 years; Table 4 'V2_AGE =",
        "Exp[theta * (AGE - 43.0)]'). Paediatric subjects: paediatric",
        "creatinine clearance was derived from Schwartz 1976 with",
        "denormalisation by Dubois & Dubois BSA, then used directly as",
        "the CRCL covariate in mL/min (see Methods 'Renal function').",
        "Allometric scaling on volumes / flows handles the bulk of the",
        "paediatric size effect; the AGE covariate captures residual",
        "non-size paediatric effects.",
        sep = " "
      ),
      source_name        = "AGE"
    ),
    CRCL = list(
      description        = "Creatinine clearance (uncorrected, Cockcroft-Gault for adults; Schwartz formula denormalised by Dubois BSA for paediatric subjects)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference 119 mL/min (Kleijn 2011 Table 4 covariate-centring",
        "value for CR). Drives the sugammadex CL renal-function effect",
        "REN = (2 * CRCL / (CRCL + 119))^1.29 (Methods), the sugammadex",
        "V2 exponential shift 'V2_CR = Exp[theta * (CR - 119)]' (Table 4),",
        "the rocuronium V1 exponential shift 'V1_CR = Exp[theta * (CR -",
        "119)]' (Table 4), and is used for typical-subject simulations",
        "in Kleijn 2011 Table 3 for the four renal-function strata. The",
        "REN saturation form drives sugammadex CL from 0 at CRCL = 0 up",
        "to 2 * 0.093 = 0.186 L/min asymptotically; at the reference 119",
        "mL/min REN = 1 and the typical sugammadex CL recovers the",
        "structural 0.093 L/min. Adult CRCL = Cockcroft-Gault, no BSA",
        "normalisation; paediatric CRCL = Schwartz / BSA_paed *",
        "BSA_adult, also raw mL/min.",
        sep = " "
      ),
      source_name        = "CRCL"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian; Caucasian + Afro-American + Hispanic; the dataset reference category)",
      notes              = paste(
        "1 = Asian, 0 = non-Asian. Multiplicative effect on rocuronium",
        "inter-compartmental clearance Q2_RAC = 1 + theta_race_q_roc =",
        "0.788 (i.e., -21.2% in Asians, Table 4) and on sugammadex",
        "central volume V1_RAC = 1 + theta_race_vc = 0.84 (i.e., -16%",
        "in Asians, Table 4). The paper's Asian cohort came primarily",
        "from study 19.4.208 Part A (Japanese patients). No race effect",
        "on rocuronium CL / V1 / V2 or on sugammadex CL / Q2 / V2.",
        sep = " "
      ),
      source_name        = "RAC"
    ),
    CONMED_SEVO = list(
      description        = "Sevoflurane volatile-anaesthesia maintenance indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no sevoflurane; propofol-only TIVA or another agent)",
      notes              = paste(
        "1 = subject received sevoflurane gas maintenance during the",
        "perioperative observation window; 0 = subject received propofol",
        "without sevoflurane. Kleijn 2011 Methods: 'Typically, patients",
        "received sevoflurane anaesthesia prior to rocuronium",
        "administration until reversal of NMB'. Multiplicative effect on",
        "the rocuronium effect-compartment equilibration rate keo_SEV =",
        "1 + theta_keo_sevo = 0.433 (i.e., -56.7% in sevoflurane patients,",
        "Table 4; equilibration slows from t1/2 5.2 min to 12 min) and on",
        "the rocuronium NMB EC50_SEV = 1 + theta_ec50_sevo = 0.605 (i.e.,",
        "-39.5% in sevoflurane patients, Table 4; EC50 drops from 1.62",
        "uM to 0.98 uM), reproducing the well-known sevoflurane",
        "potentiation of aminosteroid NMBAs. No sevoflurane effect was",
        "retained on the sugammadex-mediated reversal rate ks (Results).",
        sep = " "
      ),
      source_name        = "SEV"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 446L,
    n_studies      = 8L,
    age_range      = "1-91 years (paediatric, adolescent, adult, elderly, and old-elderly cohorts)",
    age_median     = "43 years",
    weight_range   = "9.60-139 kg (median 74.5 kg)",
    sex_female_pct = round(157 / 446 * 100, 1),
    race_ethnicity = c(
      "Non-Asian" = round(393 / 446 * 100, 1),
      "Asian"     = round(53 / 446 * 100, 1)
    ),
    disease_state  = paste(
      "Patients under general anaesthesia for elective surgery requiring",
      "rocuronium-induced neuromuscular blockade and sugammadex reversal",
      "(adults, paediatrics, and elderly), pooled with healthy male",
      "volunteers (study 19.4.101) and a renal-impairment cohort",
      "(CrCl < 30 mL/min; n = 14; study 19.4.304).",
      sep = " "
    ),
    dose_range     = paste(
      "Rocuronium IV bolus 0.6-1.2 mg/kg (maintenance doses 0.1-0.2",
      "mg/kg in some studies). Sugammadex IV bolus 0.1-16 mg/kg",
      "(placebo in some arms). Sugammadex given either at a fixed",
      "time (3 min, 5 min, 15 min) after rocuronium or at reappearance",
      "of T2 (moderate-blockade scenario).",
      sep = " "
    ),
    regions        = "Multi-regional (Japanese sub-cohort, Caucasian, Afro-American, Hispanic)",
    notes          = paste(
      "Pooled across eight model-building clinical trials (Table 1):",
      "19.4.101 part II (healthy male volunteers), 19.4.201 / 19.4.202",
      "/ 19.4.205 (adult patients), 19.4.208 part A (Japanese) and",
      "part B (Caucasian), 19.4.304 (renal impairment), 19.4.305",
      "(elderly), and 19.4.306 (paediatric, infant to adolescent).",
      "Two further trials (19.4.207 and 19.4.210) were reserved for",
      "external validation. The PK model fitted free sugammadex, free",
      "rocuronium, and complex jointly under the assumption that the",
      "two assays measured total drug; the PD model linking effect-",
      "compartment rocuronium to T4/T1 twitch ratio was fitted in two",
      "stages (placebo patients first; sugammadex-treated patients in",
      "a post hoc step) because simultaneous fitting did not",
      "converge. The complex PK parameters were SET equal to free",
      "sugammadex PK parameters (not separately estimated). Sevoflurane",
      "was tested as a covariate on keo, EC50, and ks; effects were",
      "retained only on keo and EC50. Allometric scaling exponents",
      "are FIXED at 0.75 / 1 / -0.25 (Methods 'Size effects'). The",
      "in-vitro equilibrium dissociation constant kd = 0.0559 uM",
      "(Methods, citing Bom 2002) is FIXED in the model and only the",
      "dissociation rate k2 is estimated; the association rate is",
      "derived as k1 = k2 / kd."
    )
  )

  ini({
    # ============================================================
    # SUGAMMADEX PK (Table 4 'Pharmacokinetic model sugammadex'
    # column). Parent / unsuffixed; canonical compartment names
    # `central` and `peripheral1`.
    # ============================================================
    # CL = CL_BW * REN * theta * (NOT BW-scaled) where REN absorbs
    # the size-effect on CL; see covariateData[[WT]]$notes.
    lcl  <- log(0.093);  label("Sugammadex CL at CRCL = 119 mL/min, BW = 74.5 kg (L/min)")  # Table 4 sug: CL = 0.093 L/min, RSE 1.47%
    lvc  <- log(4.42);   label("Sugammadex V1 at BW = 70 kg, non-Asian (L)")                # Table 4 sug: V1 = 4.42 L, RSE 2.35%
    lq   <- log(0.206);  label("Sugammadex Q2 at BW = 70 kg (L/min)")                       # Table 4 sug: Q2 = 0.206 L/min, RSE 4.69%
    lvp  <- log(6.35);   label("Sugammadex V2 at BW = 70 kg, CRCL = 119 mL/min (L)")        # Table 4 sug: V2 = 6.35 L, RSE 2.82%

    # Sugammadex CL covariate effects (Table 4 'Pharmacokinetic model
    # sugammadex' rows). The REN saturation form is CL = exp(lcl) *
    # (2 * CRCL / (CRCL + 119))^e_crcl_cl; at CRCL = 119 the bracket
    # is 1 and structural CL is recovered.
    e_crcl_cl <- 1.29;     label("Sugammadex CL renal-function power exponent (unitless)")  # Table 4 sug: REN = [2*CR/(CR+119)]^theta = 1.29, RSE 4.77%
    e_wt_cl   <- 0.00378;  label("Sugammadex CL linear BW shift per kg from BW = 74.5 kg")  # Table 4 sug: CL_BW = 1 + theta*(BW-74.5) = 0.00378, RSE 19.1%

    # Sugammadex V1 covariate effects.
    e_wt_vc       <- -0.00354;  label("Sugammadex V1 linear BW shift per kg from BW = 74.5 kg")  # Table 4 sug: V1_BW = 1 + theta*(BW-74.5) = -0.00354, RSE 23.1%
    e_race_asian_vc <- -0.16;   label("Sugammadex V1 multiplicative Asian-vs-non-Asian fractional change")  # Table 4 sug: V1_RAC (Asian) = 1 + theta = -0.16, RSE 22.9%

    # Sugammadex V2 covariate effects.
    e_crcl_vp <- -0.00305;  label("Sugammadex V2 exponential CRCL shift per mL/min from 119 mL/min")  # Table 4 sug: V2_CR = Exp[theta*(CR-119)] = -0.00305, RSE 18.3%

    # ============================================================
    # ROCURONIUM PK (Table 4 'Pharmacokinetic model rocuronium'
    # column). `_roc` sibling-drug suffix; compartments central_roc,
    # peripheral1_roc.
    # ============================================================
    lcl_roc <- log(0.269); label("Rocuronium CL at AGE = 43 y, BW = 70 kg (L/min)")                # Table 4 roc: CL = 0.269 L/min, RSE 1.79%
    lvc_roc <- log(4.73);  label("Rocuronium V1 at BW = 70 kg, CRCL = 119 mL/min (L)")             # Table 4 roc: V1 = 4.73 L, RSE 2.68%
    lq_roc  <- log(0.279); label("Rocuronium Q2 at BW = 70 kg, non-Asian (L/min)")                 # Table 4 roc: Q2 = 0.279 L/min, RSE 5.27%
    lvp_roc <- log(6.76);  label("Rocuronium V2 at BW = 70 kg, AGE = 43 y (L)")                    # Table 4 roc: V2 = 6.76 L, RSE 2.22%

    # Rocuronium CL covariate effects.
    e_age_cl_roc <- -0.00678; label("Rocuronium CL linear AGE shift per year from 43 years")     # Table 4 roc: CL_AGE = 1 + theta*(AGE-43) = -0.00678, RSE 16.1%

    # Rocuronium V1 covariate effects.
    e_crcl_vc_roc <- -0.00143; label("Rocuronium V1 exponential CRCL shift per mL/min from 119 mL/min")  # Table 4 roc: V1_CR = Exp[theta*(CR-119)] = -0.00143, RSE 26.4%

    # Rocuronium Q2 covariate effects.
    e_race_asian_q_roc <- -0.212; label("Rocuronium Q2 multiplicative Asian-vs-non-Asian fractional change")  # Table 4 roc: Q2_RAC (Asian) = 1 + theta = -0.212, RSE 38.2%

    # Rocuronium V2 covariate effects.
    e_age_vp_roc <- 0.00613; label("Rocuronium V2 exponential AGE shift per year from 43 years")  # Table 4 roc: V2_AGE = Exp[theta*(AGE-43)] = 0.00613, RSE 20.6%

    # ============================================================
    # COMPLEX-FORMATION KINETICS (Table 4 'Pharmacokinetic model
    # complex' column). The complex has its own central /
    # peripheral1 compartments with PK SET equal to sugammadex
    # (Methods 'Pharmacokinetic model development'); no separate
    # ini() parameters for complex PK. Only the dissociation
    # kinetics are estimated; kd is fixed from in-vitro.
    # ============================================================
    kd  <- fixed(0.0559);  label("Equilibrium dissociation constant of the sugammadex-rocuronium complex (uM, fixed from in-vitro per Bom 2002)")  # Table 4 complex: kd = 0.0559 uM, FIXED
    lk2 <- -3.38;          label("Log of complex dissociation rate constant k2 (1/min)")                                                          # Table 4 complex: log_e(k2) = -3.38, RSE 16.5%; k2 = exp(-3.38) = 0.034 1/min, association t1/2 = 1.1 min, dissociation t1/2 = 20.4 min

    # ============================================================
    # ROCURONIUM PD (Table 4 'Pharmacodynamic model rocuronium'
    # column). Effect-compartment hysteresis driving a sigmoid Emax
    # function on T4/T1 twitch ratio (x 100). Emax is SET equal to
    # E0 per the paper; the readout decreases monotonically from
    # baseline E0 toward 0 as effect-compartment rocuronium rises.
    # ============================================================
    lke0  <- log(0.134); label("Rocuronium effect-compartment equilibration rate at BW = 70 kg, no sevoflurane (1/min)")  # Table 4 roc PD: keo = 0.134 1/min, RSE 6.49%
    lec50 <- log(1.62);  label("Rocuronium NMB EC50 in effect compartment, no sevoflurane (uM)")                          # Table 4 roc PD: EC50 = 1.62 uM, RSE 3.68%
    lhill <- log(7.52);  label("Hill exponent of the NMB sigmoid Emax (unitless)")                                        # Table 4 roc PD: Hill = 7.52, RSE 5.84%
    le0   <- log(104);   label("Typical baseline T4/T1 twitch ratio (paper reports E0 = theta * 100 with theta = 1.04 so E0 = 104 on T4/T1 x 100 scale)")  # Table 4 roc PD: E0 = theta * 100 = 1.04 * 100 = 104, RSE 1.51%

    # Sevoflurane covariate effects on keo and EC50 (Table 4 roc PD).
    e_conmed_sevo_ke0  <- -0.567; label("Sevoflurane multiplicative fractional change on keo (unitless)")   # Table 4 roc PD: keo_SEV = 1 + theta = -0.567, RSE 8.94%; sevoflurane slows equilibration t1/2 from 5.2 min to 12 min
    e_conmed_sevo_ec50 <- -0.395; label("Sevoflurane multiplicative fractional change on EC50 (unitless)")  # Table 4 roc PD: EC50_SEV = 1 + theta = -0.395, RSE 12.0%; sevoflurane drops EC50 from 1.62 uM to 0.98 uM (40% potentiation)

    # ============================================================
    # SUGAMMADEX PD (Table 4 'Pharmacodynamic model sugammadex'
    # column). Sugammadex modulates rocuronium elimination from the
    # rocuronium effect compartment via a second-order rate
    # constant; no separate sugammadex effect compartment is
    # required (Bom 2002 framework). ks is allometrically scaled
    # with exponent -0.25, same as keo.
    # ============================================================
    lks <- -3.43;  label("Log of sugammadex-mediated rocuronium effect-compartment elimination rate ks at BW = 70 kg (1/(min*uM))")  # Table 4 sug PD: log_e(ks) = -3.43, RSE 0.222%; ks = exp(-3.43) = 0.0324 1/(min*uM) (paper text rounds to 0.033); reversal t1/2 at sug 2 mg/kg = 1.43 min, sug 4 mg/kg = 0.71 min

    # ============================================================
    # IIV (Table 4 'IIV ... (shrinkage XX%)' rows reported in CV%).
    # omega^2 = log(1 + CV^2). The paper retains only one off-
    # diagonal correlation in the final model (keo - EC50 r = 0.37);
    # all other parameters have independent log-normal IIV.
    # ============================================================
    etalcl ~ log(1 + 0.224^2)  # Table 4 sug: IIV CL = 22.4% CV; shrinkage 25%

    etalcl_roc ~ log(1 + 0.324^2)  # Table 4 roc: IIV CL = 32.4% CV; shrinkage 12%
    etalvc_roc ~ log(1 + 0.238^2)  # Table 4 roc: IIV V1 = 23.8% CV; shrinkage 28%
    etalvp_roc ~ log(1 + 0.322^2)  # Table 4 roc: IIV V2 = 32.2% CV; shrinkage 37%

    # Correlated keo + EC50 IIV (Table 4 roc PD).
    # - IIV keo  = 41.7% CV (shrinkage 0%) -> var_keo  = log(1 + 0.417^2)
    # - IIV EC50 = 24.9% CV (shrinkage 0%) -> var_ec50 = log(1 + 0.249^2)
    # - r(keo, EC50) = 0.37 (95% CI 0.191-0.455, RSE 32.9%)
    # - cov = r * sqrt(var_keo * var_ec50)
    #       = 0.37 * sqrt(0.15976 * 0.06024)
    #       = 0.0363
    etalke0 + etalec50 ~ c(
      log(1 + 0.417^2),
      0.37 * sqrt(log(1 + 0.417^2) * log(1 + 0.249^2)),
      log(1 + 0.249^2)
    )
    etalhill ~ log(1 + 0.411^2)  # Table 4 roc PD: IIV Hill = 41.1% CV; shrinkage 2%
    etale0   ~ log(1 + 0.111^2)  # Table 4 roc PD: IIV E0 = 11.1% CV; shrinkage 2%

    etalks ~ log(1 + 1.14^2)  # Table 4 sug PD: IIV ks = 114% CV; shrinkage 3%

    # ============================================================
    # Residual error (Table 4 'Residual error' rows). The paper used
    # log-additive residual on plasma concentrations -- in nlmixr2
    # the proportional form ~ prop(propSd) reproduces this for the
    # CV range reported. Cc (parent / sugammadex) and Cc_roc are
    # both total drug (free + complex; the assays did not
    # discriminate). The NMB readout uses an additive residual on
    # the T4/T1 x 100 scale; the paper reports two PD residuals,
    # 2.70 (rocuronium-PD fit on placebo patients) and 3.24
    # (sugammadex-PD fit on combined data), reflecting the
    # sequential two-stage estimation. Using 3.24 here because the
    # primary use of this model is forward simulation of reversal
    # scenarios under both drugs.
    # ============================================================
    propSd     <- 0.363;  label("Sugammadex (parent) total-plasma proportional residual SD (fraction)")          # Table 4 sug: residual error 36.3% (log-additive), shrinkage 3%
    propSd_roc <- 0.200;  label("Rocuronium total-plasma proportional residual SD (fraction)")                   # Table 4 roc: residual error 20.0% (log-additive), shrinkage 14%
    addSd_NMB  <- 3.24;   label("NMB T4/T1 x 100 additive residual SD; sugammadex-PD-fit residual per Table 4")  # Table 4 sug PD: residual error 3.24 (T4/T1 x 100 scale), shrinkage 0%; rocuronium-PD-fit residual 2.70 (Table 4) applies to the placebo subset that fitted the rocuronium spontaneous-recovery PD stage
  })

  model({
    # ---------------------------------------------------------------
    # Allometric ratio (reference WT = 70 kg; Kleijn 2011 Table 4
    # footnote and Methods 'Size effects'). Exponents are FIXED at
    # 0.75 / 1 / -0.25 (per Methods); not estimated.
    # ---------------------------------------------------------------
    wt_ratio <- WT / 70

    # ---------------------------------------------------------------
    # SUGAMMADEX PK individual parameters
    # ---------------------------------------------------------------
    # REN saturates from 0 (CRCL = 0) to 2 at CRCL >> 119; at the
    # reference CRCL = 119 mL/min, REN = 1 and structural CL is
    # recovered.
    ren_factor <- (2 * CRCL / (CRCL + 119))^e_crcl_cl
    cl_bw_sug  <- 1 + e_wt_cl * (WT - 74.5)
    cl <- exp(lcl + etalcl) * ren_factor * cl_bw_sug
    # Note: NO allometric (BW/70)^0.75 on sugammadex CL per paper.

    vc_bw_sug   <- 1 + e_wt_vc * (WT - 74.5)
    vc_race_sug <- 1 + e_race_asian_vc * RACE_ASIAN
    vc <- exp(lvc) * vc_bw_sug * vc_race_sug * wt_ratio

    q  <- exp(lq) * wt_ratio^0.75

    vp <- exp(lvp) * exp(e_crcl_vp * (CRCL - 119)) * wt_ratio

    # Micro-rates for the sugammadex PK subsystem (also reused for
    # the complex subsystem since complex PK = sugammadex PK).
    kel_sug <- cl / vc
    k12_sug <- q  / vc
    k21_sug <- q  / vp

    # ---------------------------------------------------------------
    # ROCURONIUM PK individual parameters
    # ---------------------------------------------------------------
    cl_age_roc <- 1 + e_age_cl_roc * (AGE - 43)
    cl_roc <- exp(lcl_roc + etalcl_roc) * cl_age_roc * wt_ratio^0.75

    vc_roc <- exp(lvc_roc + etalvc_roc) * exp(e_crcl_vc_roc * (CRCL - 119)) * wt_ratio

    q_race_roc <- 1 + e_race_asian_q_roc * RACE_ASIAN
    q_roc <- exp(lq_roc) * q_race_roc * wt_ratio^0.75

    vp_roc <- exp(lvp_roc + etalvp_roc) * exp(e_age_vp_roc * (AGE - 43)) * wt_ratio

    kel_roc <- cl_roc / vc_roc
    k12_roc <- q_roc  / vc_roc
    k21_roc <- q_roc  / vp_roc

    # ---------------------------------------------------------------
    # COMPLEX-FORMATION kinetics. kd is fixed from in-vitro; k2 is
    # estimated on the log-scale; k1 is derived to maintain the
    # equilibrium constraint kd = k2 / k1.
    # ---------------------------------------------------------------
    k2 <- exp(lk2)
    k1 <- k2 / kd

    # ---------------------------------------------------------------
    # ROCURONIUM PD individual parameters
    # ---------------------------------------------------------------
    ke0_sev_factor <- 1 + e_conmed_sevo_ke0 * CONMED_SEVO
    ke0 <- exp(lke0 + etalke0) * ke0_sev_factor * wt_ratio^(-0.25)

    ec50_sev_factor <- 1 + e_conmed_sevo_ec50 * CONMED_SEVO
    ec50 <- exp(lec50 + etalec50) * ec50_sev_factor

    hill <- exp(lhill + etalhill)
    e0   <- exp(le0   + etale0)

    # ---------------------------------------------------------------
    # SUGAMMADEX PD individual parameter
    # ---------------------------------------------------------------
    ks <- exp(lks + etalks) * wt_ratio^(-0.25)

    # ---------------------------------------------------------------
    # Plasma concentrations (uM = umol / L). Free / complex are
    # tracked separately in distinct compartments; assay reports
    # total = free + complex.
    # ---------------------------------------------------------------
    csug_free <- central          / vc      # uM free sugammadex
    croc_free <- central_roc      / vc_roc  # uM free rocuronium
    ccomplex  <- central_complex  / vc      # uM complex (in sug-volume)

    # ---------------------------------------------------------------
    # Complex-formation flux (umol/min) in the central compartment.
    # Per Bom 2002 / Kleijn 2011 Methods, complexation happens only
    # in the central compartments (not in peripheral). The forward
    # flux is k1 * Csug_free * Croc_free with units 1/(min*uM) *
    # uM * uM = uM/min, scaled to amount by vc (the complex's
    # central volume, which equals the sugammadex central volume).
    # ---------------------------------------------------------------
    flux_form <- k1 * csug_free * croc_free * vc
    flux_diss <- k2 * central_complex

    # ---------------------------------------------------------------
    # ODEs -- amounts in umol. Both drugs are IV bolus; dose lands
    # directly in the central compartment (no depot for either).
    # ---------------------------------------------------------------
    # Sugammadex free (central + peripheral1)
    d/dt(central) <- -kel_sug * central -
                      k12_sug * central + k21_sug * peripheral1 -
                      flux_form + flux_diss
    d/dt(peripheral1) <- k12_sug * central - k21_sug * peripheral1

    # Rocuronium free (central_roc + peripheral1_roc); no
    # sugammadex-mediated elimination of free rocuronium from the
    # central compartment -- the sugammadex action is on the effect
    # compartment, per Methods 'Pharmacokinetic-pharmacodynamic model
    # development'. Mass balance of free rocuronium in central is
    # standard 2-compartment plus complex formation.
    d/dt(central_roc) <- -kel_roc * central_roc -
                          k12_roc * central_roc + k21_roc * peripheral1_roc -
                          flux_form + flux_diss
    d/dt(peripheral1_roc) <- k12_roc * central_roc - k21_roc * peripheral1_roc

    # Complex (central_complex + peripheral1_complex), PK = sugammadex
    d/dt(central_complex) <- -kel_sug * central_complex -
                              k12_sug * central_complex + k21_sug * peripheral1_complex +
                              flux_form - flux_diss
    d/dt(peripheral1_complex) <- k12_sug * central_complex - k21_sug * peripheral1_complex

    # Rocuronium effect compartment. Track concentration directly
    # (standard nlmixr2 effect-compartment idiom: the effect
    # compartment has nominal unit volume so the ODE state IS the
    # effect-site concentration in uM). Sugammadex-mediated reversal
    # appears as an extra first-order elimination of the effect-site
    # rocuronium driven by central free-sugammadex concentration.
    d/dt(effect_roc) <- ke0 * (croc_free - effect_roc) -
                        ks * csug_free * effect_roc

    # ---------------------------------------------------------------
    # Outputs. Plasma assays measured total drug (free + complex);
    # Methods 'Rocuronium and sugammadex plasma concentrations'.
    # ---------------------------------------------------------------
    Cc     <- csug_free + ccomplex   # Total sugammadex in central (uM)
    Cc_roc <- croc_free + ccomplex   # Total rocuronium in central (uM)

    # ---------------------------------------------------------------
    # NMB readout: T4/T1 twitch ratio on the x 100 scale. Sigmoid
    # Emax form with Emax = E0 forced (Methods 'Pharmacokinetic-
    # pharmacodynamic model development'). When effect-site
    # rocuronium = 0, NMB = E0 (baseline). When effect-site
    # rocuronium >> EC50, NMB approaches 0 (full neuromuscular
    # block). Reversal-time endpoints TOF90 / TOF70 are read off
    # this trace at NMB = 90 and NMB = 70.
    # ---------------------------------------------------------------
    NMB <- e0 * (1 - effect_roc^hill / (ec50^hill + effect_roc^hill))

    # Residual error.
    Cc     ~ prop(propSd)
    Cc_roc ~ prop(propSd_roc)
    NMB    ~ add(addSd_NMB)
  })
}
