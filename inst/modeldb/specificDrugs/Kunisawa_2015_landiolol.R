Kunisawa_2015_landiolol <- function() {
  description <- "Two-compartment intravenous population PK model with lag time for landiolol hydrochloride (an ultra-short-acting cardioselective beta1-adrenergic receptor blocker) in adult patients with peripheral arterial disease undergoing peripheral arterial surgery, with linear body-weight normalization on CL, Vc, Q and Vp (Kunisawa 2015)"
  reference <- paste(
    "Kunisawa T, Yamagishi A, Suno M, Nakade S, Honda N, Kurosawa A,",
    "Sugawara A, Tasaki Y, Iwasaki H.",
    "Target-controlled infusion and population pharmacokinetics of landiolol",
    "hydrochloride in patients with peripheral arterial disease.",
    "Ther Clin Risk Manag. 2015;11:107-114.",
    "doi:10.2147/TCRM.S74867."
  )
  vignette <- "Kunisawa_2015_landiolol"
  units <- list(time = "hr", dosing = "ug", concentration = "ng/mL") # Methods: TCI delivers landiolol hydrochloride at ug/kg/min infusion rates, plasma assayed by HPLC-fluorescence (Suno et al. 2009); time converted from minutes-as-reported to hours for nlmixr2lib convention

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear body-weight normalization (per-kg reporting in Table 2). Parameters reported in mL/min/kg (CL, Q) and mL/kg (V1, V2); the underlying NONMEM model scales individual CL/Vc/Q/Vp linearly by WT and the published typical values are the per-kg coefficients. Implemented as (WT/70)^1 with reference weight 70 kg; cohort mean weight 58.9 kg (range 40.4-71.8 kg, Table 1). Body weight was also tested as an additional covariate via forward selection but no significant effect was retained beyond the built-in linear normalization (Results).",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",                              # Methods: adult patients undergoing peripheral arterial surgery
    n_subjects     = 8L,                                   # Methods + Table 1: eight patients with peripheral arterial disease
    n_studies      = 1L,                                   # Single-centre prospective TCI study (UMIN000015077)
    age_range      = "64-84 years",                        # Table 1
    age_median     = "73 years",                           # Table 1 (mean 73, SD 7)
    weight_range   = "40.4-71.8 kg",                       # Table 1
    weight_median  = "58.9 kg",                            # Table 1 (mean 58.9, SD 10.9)
    sex_female_pct = 25,                                   # Table 1: 6 male / 2 female (8 total)
    race_ethnicity = c(Asian = 100),                       # Asahikawa Medical University Hospital, Japan
    disease_state  = "Adult patients scheduled for peripheral arterial surgery (peripheral arterial disease, PAD); ASA physical status 2 or 3; excluded if pre-existing arrhythmia (atrial fibrillation, conduction-system disturbance) or recent treatment with alpha-methyldopa, clonidine, or beta-blockers (Methods).",
    dose_range     = "Target-controlled IV infusion (Harvard pump under STANPUMP control using Honda et al's two-compartment parameters) of landiolol hydrochloride at target plasma concentrations of 500 ng/mL and 1,000 ng/mL for 30 min each (~50% and 100% of the highest clinical dose range 10-40 ug/kg/min from the package insert).",
    regions        = "Japan (Asahikawa Medical University, Hokkaido)",
    sampling       = "112 plasma concentrations across 8 subjects (rich). Samples drawn at 1, 2, 5, and 25 min after starting each TCI segment; at the target-concentration change; and at 1, 2, 5, 10, 15, and 20 min after the end of infusion (Methods + Figure 1).",
    co_medication  = "General anaesthesia maintained with TCI propofol (Diprifusor, BIS-titrated 40-60) and remifentanil (Minto TCI to effect-site 2 ng/mL pre-incision, raised to 8 ng/mL pre-surgical-stimulus). Rocuronium 0.6 mg/kg for intubation. Dopamine 3 ug/kg/min for haemodynamic stability starting 20 min before skin incision. Rescue: atropine 0.5 mg IV for HR <= 45 bpm (none required); ephedrine 5 mg IV for hypotension with HR <= 60 bpm; phenylephrine 0.05 mg IV for hypotension without bradycardia. No co-administration of beta-blockers, alpha-methyldopa, or clonidine.",
    assay          = "HPLC with fluorescence detection per Suno et al. (J Chromatogr B 2008/2009); samples collected in chilled ethanol/EDTA-2Na/neostigmine to prevent pseudocholinesterase-mediated hydrolysis of the ester.",
    notes          = "Baseline laboratory values (Table 1) were mostly within normal range; mildly abnormal albumin, cholinesterase, BUN, and serum creatinine were not retained as covariates because of the limited frequency and extent of abnormalities. Body weight, lean body mass, and age were the demographic covariates tested via forward selection; none were retained. PAD-cohort V1 and CL were approximately 64% and 84% of the matched healthy-volunteer values (Conclusion). Estimation method: FOCE-I in NONMEM VII level 1.2 (ICON Development Solutions). Final model selected by Akaike information criterion (AIC of 2-cmt + lag = 1,247.171); 1- and 3-compartment models did not converge."
  )

  ini({
    # Structural PK parameters -- reference weight 70 kg (Kunisawa 2015 Table 2 reports mL/min/kg
    # and mL/kg per-kg values; converted here to L/h and L at 70 kg for nlmixr2lib units).
    lcl  <- log(128.94);  label("Total clearance CL at 70 kg (L/h)")                              # Table 2 (PAD column): TVCL = 30.7 mL/min/kg (SE 2.08) = 30.7 * 70 * 60 / 1000 = 128.94 L/h at 70 kg
    lvc  <- log(4.55);    label("Central volume of distribution Vc at 70 kg (L)")                 # Table 2 (PAD column): TVV1 = 65.0 mL/kg (SE 4.97) = 65.0 * 70 / 1000 = 4.55 L at 70 kg
    lq   <- log(202.86);  label("Intercompartmental clearance Q at 70 kg (L/h)")                  # Table 2 (PAD column): TVQ = 48.3 mL/min/kg (SE 16.4) = 48.3 * 70 * 60 / 1000 = 202.86 L/h at 70 kg
    lvp  <- log(3.808);   label("Peripheral volume of distribution Vp at 70 kg (L)")              # Table 2 (PAD column): TVV2 = 54.4 mL/kg (SE 4.54) = 54.4 * 70 / 1000 = 3.808 L at 70 kg
    ltlag <- log(0.01055); label("Lag time on central (h)")                                        # Table 2 (PAD column): TVALAG = 0.633 min (SE 0.000173) = 0.633 / 60 = 0.01055 h

    # Body-weight scaling exponents -- fixed at 1 because the source paper reports parameters per-kg (linear scaling), not estimated allometric exponents
    e_wt_cl <- fixed(1); label("Body-weight scaling exponent on CL (linear normalization)")       # Table 2 reports CL in mL/min/kg (linear per-kg scaling)
    e_wt_vc <- fixed(1); label("Body-weight scaling exponent on Vc (linear normalization)")       # Table 2 reports V1 in mL/kg (linear per-kg scaling)
    e_wt_q  <- fixed(1); label("Body-weight scaling exponent on Q (linear normalization)")        # Table 2 reports Q in mL/min/kg (linear per-kg scaling)
    e_wt_vp <- fixed(1); label("Body-weight scaling exponent on Vp (linear normalization)")       # Table 2 reports V2 in mL/kg (linear per-kg scaling)

    # IIV -- exponential error model on CL only (Results + Table 2). V1, Q, V2, and ALAG were
    # tested with random effects but not retained in the final model.
    # NONMEM exponential IIV CL_i = TVCL * exp(eta_CL); the omega^2 estimate 0.0183 is the
    # variance on the log scale of CL and maps directly to etalcl.
    etalcl ~ 0.0183 # Table 2 (PAD column): omega_CL^2 = 0.0183 (SE 0.00768); reported %CV = 13.5

    # Residual error -- "exponential error model" (Methods + Results). NONMEM Y = F * EXP(EPS)
    # is exactly a log-normal residual: log(Y) = log(F) + EPS with EPS ~ N(0, sigma^2).
    # This maps to nlmixr2 lnorm() with expSd = sqrt(sigma^2) on the log scale.
    expSd <- 0.278; label("Lognormal residual SD on log(Cc) (unitless on log scale)")             # Table 2 (PAD column): sigma^2 = 0.0773 (SE 0.0234); reported %CV = 27.8; expSd = sqrt(0.0773) = 0.278
  })

  model({
    # Individual structural parameters with linear body-weight scaling to a 70 kg reference.
    cl   <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc   <- exp(lvc)          * (WT / 70)^e_wt_vc
    q    <- exp(lq)           * (WT / 70)^e_wt_q
    vp   <- exp(lvp)          * (WT / 70)^e_wt_vp
    lagt <- exp(ltlag)

    # Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment IV with lag on central (TCI infusion is recorded into `central`,
    # ALAG shifts both the start and end of each infusion segment by lagt). The
    # original NONMEM model used the standard ADVAN3 / TRANS4 ALAG mechanism.
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    lag(central)      <- lagt

    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
