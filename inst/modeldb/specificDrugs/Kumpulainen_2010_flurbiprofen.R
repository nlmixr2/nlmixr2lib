Kumpulainen_2010_flurbiprofen <- function() {
  description <- "Three-compartment population PK model with a separate cerebrospinal-fluid (CSF) compartment for flurbiprofen in 64 healthy children aged 3 months to 13 years (Kumpulainen 2010). Two parallel absorption routes: oral syrup via an absorption compartment with lag time (K12) and a single first-order ka, and IV flurbiprofen axetil prodrug via a separate dosing compartment that converts to flurbiprofen with a first-order rate constant (K42). Plasma kinetics scaled allometrically by weight (exponents fixed at 0.75 for CL and 1 for all volumes, including the CSF volume held fixed at 0.15 L/70 kg per literature). The paper's QCSF + UPTK parameterisation is encoded as canonical influx / efflux clearances clin = QCSF * UPTK and clef = QCSF, with fraction unbound (fu) gating only the central-to-CSF flux."
  reference <- paste(
    "Kumpulainen E, Valitalo P, Kokki M, Lehtonen M, Hooker A, Ranta V-P, Kokki H.",
    "Plasma and cerebrospinal fluid pharmacokinetics of flurbiprofen in children.",
    "Br J Clin Pharmacol. 2010;70(4):557-566.",
    "doi:10.1111/j.1365-2125.2010.03720.x.",
    sep = " "
  )
  vignette <- "Kumpulainen_2010_flurbiprofen"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Drives allometric scaling on CL (exponent fixed at 0.75) and on all volumes V_central, V_shallow, V_deep, V_CSF (shared exponent fixed at 1.0). Reference weight 70 kg per the paper's (WT/70) normalisation. Median weight in the study was 20 kg (range 7-76 kg).",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 64L,
    n_studies      = 1L,
    age_range      = "3 months to 13 years",
    age_median     = "5.2 years",
    weight_range   = "7-76 kg",
    weight_median  = "20 kg",
    height_range   = "60-167 cm",
    height_median  = "110 cm",
    sex_female_pct = 23.4,
    race_ethnicity = "Not reported (single-centre Finnish cohort; 64 healthy children scheduled for elective lower-body surgery with spinal anaesthesia)",
    disease_state  = "Healthy children scheduled for elective lower-body surgery (herniotomy, orthopaedic, genitourinary, other) under spinal anaesthesia. EudraCT 2006-000310-20; Research Ethics Committee of the Hospital District of Northern Savo (no. 12/2006).",
    dose_range     = "Single preoperative dose: 37 children received 1 mg/kg oral flurbiprofen syrup (Froben, 5 mg/mL, Abbott); 27 children received a 10-min IV injection of 0.9 mg/kg flurbiprofen axetil (Ropion 10 mg/mL, Kaken Pharmaceutical) which corresponds to approximately 0.65 mg/kg flurbiprofen equivalent. The model receives the flurbiprofen-equivalent mass on the IV depot2 compartment.",
    regions        = "Finland (Kuopio University Hospital, single centre).",
    notes          = "304 total plasma + 62 protein-free plasma + 60 total CSF flurbiprofen concentrations from 64 children. NONMEM VI with first-order conditional estimation with interaction (FOCE-I). Premedication: transmucosal midazolam + ketamine; intraoperative sedation: propofol + thiopental; postoperative analgesia: paracetamol + ketoprofen with fentanyl rescue. Concomitant medications were not modelled."
  )

  ini({
    # ============================================================
    # Structural PK parameters - Kumpulainen 2010 Table 2 final
    # model. Each value is the back-transformed typical-value
    # estimate from NONMEM; the ini() carries log(value) per the
    # nlmixr2lib log-transformed-PK convention. Reference weight
    # 70 kg.
    # ============================================================

    # Absorption (parallel two-route)
    lka    <- log(5.5);  label("Oral absorption rate constant K12 from depot to central (1/h)")          # Table 2: K12 = 5.5 1/h (RSE 0.24)
    lka2   <- log(29);   label("IV flurbiprofen axetil conversion rate constant K42 from depot2 to central (1/h)") # Table 2: K42 = 29 1/h (RSE 0.32); half-life 1.4 min
    ltlag  <- log(0.11); label("Lag time for oral absorption (h)")                                       # Table 2: lag time = 0.11 h (RSE 0.16) = 6 min

    # Disposition
    lcl    <- log(0.96); label("Flurbiprofen total clearance CL (L/h, 70 kg reference)")                 # Table 2: CL = 0.96 L/h/70kg (RSE 0.057)
    lvc    <- log(3.6);  label("Central volume of distribution V_central (L, 70 kg reference)")          # Table 2: V_central = 3.6 L/70kg (RSE 0.11)
    lq     <- log(1.5);  label("Intercompartmental clearance to shallow peripheral Q2 (L/h)")            # Table 2: Q2 = 1.5 L/h (RSE 0.39); paper's Q (shallow peripheral) maps to nlmixr2 lq (peripheral1)
    lvp    <- log(1.8);  label("Shallow peripheral volume V_shallow (L, 70 kg reference)")               # Table 2: V_shallow = 1.8 L/70kg (RSE 0.20)
    lq2    <- log(0.18); label("Intercompartmental clearance to deep peripheral Q (L/h)")                # Table 2: Q (deep peripheral) = 0.18 L/h (RSE 0.30); paper's Q (deep peripheral) maps to nlmixr2 lq2 (peripheral2)
    lvp2   <- log(2.7);  label("Deep peripheral volume V_deep (L, 70 kg reference)")                     # Table 2: V_deep = 2.7 L/70kg (RSE 0.18)
    lfdepot <- log(0.81); label("Oral bioavailability F (unitless)")                                     # Table 2: F = 0.81 (RSE 0.055)

    # CSF distribution - reparameterised from the paper's QCSF + UPTK
    # form to the canonical clin / clef (Campagne 2019 precedent).
    # Paper Table 2 reports QCSF = 0.12 L/h (RSE 0.27) and the active
    # uptake multiplier UPTK = 6.8 (RSE 0.070). The paper's transfer
    # rates are:
    #   K_central_to_CSF = QCSF * UPTK * fu / V_central
    #   K_CSF_to_central = QCSF / V_CSF
    # Equivalently:
    #   clin = QCSF * UPTK    -> influx clearance applied to fu * Cc
    #   clef = QCSF           -> efflux clearance applied to Ccsf
    # The implied uptake ratio UPTK = clin / clef = 6.8 is preserved by
    # construction; the vignette's "Assumptions and deviations"
    # documents the reparameterisation.
    lclin  <- log(0.12 * 6.8); label("CSF influx clearance clin = QCSF * UPTK (L/h)")                    # Table 2: QCSF * UPTK = 0.12 * 6.8 = 0.816 L/h (uptake ratio UPTK = 6.8)
    lclef  <- log(0.12);       label("CSF efflux clearance clef = QCSF (L/h)")                           # Table 2: QCSF = 0.12 L/h (RSE 0.27)
    lvcsf  <- fixed(log(0.15)); label("CSF volume of distribution V_CSF (L, 70 kg reference; fixed)")    # Text (Results, popPK model): "The volume of the CSF compartment (0.15 l) was retrieved from the literature [13] and also scaled by weight." Held fixed at 0.15 L/70 kg.

    # Plasma protein binding
    lfu    <- log(0.00031); label("Plasma fraction unbound fu (unitless)")                                # Table 2: protein-free fraction = 0.00031 (= 0.031%, RSE 0.043); also Results "fraction unbound found to be constant... protein-free fraction of 0.031%"

    # ============================================================
    # Allometric scaling exponents - fixed per Kumpulainen 2010
    # Discussion and Table 2 footnote: estimated exponents (1.01 for
    # V, 0.774 for CL) were "in close agreement with" the canonical
    # allometric values 1 and 0.75 and were "fixed to 1 and 0.75 in
    # this study". The 1.0 exponent applies to all volumes
    # (V_central, V_shallow, V_deep, V_CSF) per the (WT/70)
    # annotation on each volume row in Table 2 and the text "the
    # volume of the CSF compartment... also scaled by weight". The
    # 0.75 exponent applies only to CL (elimination); the
    # intercompartmental clearances Q2, Q and QCSF carry no WT
    # scaling in Table 2 and are encoded the same way here.
    # ============================================================
    e_wt_cl    <- fixed(0.75); label("Allometric exponent on CL (fixed)")                                # Table 2: CL (l/h) * (WT/70)^0.75; Discussion: "the exponents were fixed to 1 and 0.75 in this study"
    e_wt_vc_vp <- fixed(1);    label("Allometric exponent shared across all volumes V_central / V_shallow / V_deep / V_CSF (fixed)") # Table 2: V_central / V_shallow / V_deep all * (WT/70); Results: V_CSF "also scaled by weight"

    # ============================================================
    # Inter-individual variability - Kumpulainen 2010 Table 2.
    # The omega values reported in the table are SDs of eta on the
    # log scale (consistent with the approximately 30% CV
    # description in Results: omega_CL = 0.28 ~ 28% CV, omega_Vd =
    # 0.28 ~ 28%, omega_K12 = 0.81 ~ 80%, omega_fu = 0.30 ~ 30%).
    # Internal nlmixr2 variance = (paper omega)^2. The exponential
    # model P_i = TVP * exp(eta_i) is the standard NONMEM
    # log-normal form (Methods: "BSV was assigned to volume of
    # distribution, clearance, oral absorption rate and fraction
    # unbound with an exponential model"). BSV is uncorrelated
    # (Methods: "Between-subject variability was best described as
    # uncorrelated"). BSV is applied only to V_central (the paper's
    # singular "volume of distribution"); the peripheral and CSF
    # volumes carry no IIV in Table 2.
    # ============================================================
    etalcl ~ 0.0784  # omega_CL = 0.28 (RSE 0.20; bootstrap 0.22-0.34) -> variance 0.28^2 = 0.0784
    etalvc ~ 0.0784  # omega_Vd = 0.28 (RSE 0.25; bootstrap 0.19-0.35) -> variance 0.28^2 = 0.0784
    etalka ~ 0.6561  # omega_K12 = 0.81 (RSE 0.40; bootstrap 0.42-1.4) -> variance 0.81^2 = 0.6561
    etalfu ~ 0.09    # omega_fu = 0.30 (RSE 0.28; bootstrap 0.19-0.36) -> variance 0.30^2 = 0.09

    # ============================================================
    # Residual error - Kumpulainen 2010 Table 2. Proportional model
    # with separate variances for plasma and CSF observations
    # (Methods: "a proportional error model with separate variances
    # for plasma and CSF observations proved adequate"). The sigma
    # values are SDs of the proportional residual (linear scale,
    # i.e. fractional CV). The same plasma sigma is applied to both
    # total and unbound plasma observations (Methods: "Unbound
    # observations were included in the same compartment as total
    # plasma observations. They were defined as special cases of
    # total plasma observations"); nlmixr2 requires per-output
    # residual parameters, so propSd (total Cc) and propSd_Cu
    # (unbound Cu) carry the same source value.
    # ============================================================
    propSd      <- 0.13; label("Proportional residual SD on total plasma Cc (fraction)")            # Table 2: sigma_plasma = 0.13 (RSE 0.065; bootstrap 0.11-0.14); shared between total and unbound plasma observations in the source
    propSd_Cu   <- 0.13; label("Proportional residual SD on unbound plasma Cu (fraction)")          # Table 2: sigma_plasma = 0.13 applied to unbound plasma observations per Methods modelling strategy
    propSd_Ccsf <- 0.50; label("Proportional residual SD on CSF concentration Ccsf (fraction)")     # Table 2: sigma_CSF = 0.50 (RSE 0.091; bootstrap 0.43-0.60)
  })

  model({
    # ---- Individual PK parameters (back-transform + allometric) ----
    # Reference weight 70 kg per Table 2's (WT/70) annotation. The
    # 0.75 exponent on CL and 1.0 exponent on all volumes are fixed
    # (Table 2 footnote, Discussion). The intercompartmental
    # clearances Q2, Q and QCSF are not allometrically scaled by
    # weight in Table 2; the (WT/70) factor appears only on the CL
    # and V rows in the source table, so it is applied here only
    # to cl and all volumes.
    cl    <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc    <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    vp    <- exp(lvp)          * (WT / 70)^e_wt_vc_vp
    vp2   <- exp(lvp2)         * (WT / 70)^e_wt_vc_vp
    vcsf  <- exp(lvcsf)        * (WT / 70)^e_wt_vc_vp
    q     <- exp(lq)
    q2    <- exp(lq2)
    ka    <- exp(lka + etalka)
    ka2   <- exp(lka2)
    tlag  <- exp(ltlag)
    clin  <- exp(lclin)
    clef  <- exp(lclef)
    fu    <- exp(lfu + etalfu)

    # ---- Concentrations -------------------------------------------
    Cc    <- central / vc       # total plasma concentration (mg/L)
    Cu    <- fu * Cc            # unbound plasma concentration (mg/L)
    Ccsf  <- csf / vcsf         # CSF concentration (mg/L)

    # ---- ODE system (Figure 1 + Methods 'Modelling strategy') -----
    # Two parallel absorption depots feed central:
    #   depot  (oral syrup, absorption rate ka = K12, lag tlag)
    #   depot2 (IV flurbiprofen axetil prodrug, conversion rate ka2 = K42)
    # Central exchanges with two peripheral compartments (shallow
    # peripheral1 via q; deep peripheral2 via q2) and a CSF
    # compartment via active influx clin = QCSF * UPTK applied to
    # fu * Cc and passive efflux clef = QCSF applied to Ccsf.
    # Elimination is first-order from central (cl).
    d/dt(depot)       <- -ka * depot
    d/dt(depot2)      <- -ka2 * depot2
    d/dt(central)     <-  ka * depot + ka2 * depot2 -
                          cl * Cc -
                          q * (Cc - peripheral1 / vp) -
                          q2 * (Cc - peripheral2 / vp2) -
                          clin * fu * Cc + clef * Ccsf
    d/dt(peripheral1) <-  q  * (Cc - peripheral1 / vp)
    d/dt(peripheral2) <-  q2 * (Cc - peripheral2 / vp2)
    d/dt(csf)         <-  clin * fu * Cc - clef * Ccsf

    # ---- Bioavailability and lag time ------------------------------
    # Oral syrup: estimated bioavailability F on the depot
    # compartment. IV axetil: structural F = 1 anchor (the depot2
    # dose record carries the flurbiprofen-equivalent mass set by
    # the user, so no additional scaling is needed).
    f(depot)    <- exp(lfdepot)
    alag(depot) <- tlag

    # ---- Residual error -------------------------------------------
    # Plasma sigma applied to both total and unbound plasma
    # observations (shared source value); CSF sigma applied to CSF.
    Cc   ~ prop(propSd)
    Cu   ~ prop(propSd_Cu)
    Ccsf ~ prop(propSd_Ccsf)
  })
}
