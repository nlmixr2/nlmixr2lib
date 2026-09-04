Edwards_2016_obeticholicAcid_pbpk <- function() {
  description <- "PBPK (semi-mechanistic bile-acid enterohepatic-recirculation model). Obeticholic acid (OCA) and its glycine (glyco-OCA) and taurine (tauro-OCA) conjugates in healthy adults and subjects with Child-Pugh A/B/C cirrhosis. 19 ODEs spanning systemic, portal, sinusoidal, liver, bile duct, gallbladder and gut spaces per analyte, plus an oral depot; blood/bile flows act on concentration, transport and biotransformation on amount. Meal-gated gallbladder emptying; four hepatic-impairment mechanisms (reduced hepatic uptake, portal-systemic shunting with arterial buffer response, reduced functional liver volume, preferential tauro-conjugation)."
  reference <- "Edwards JE, LaCerte C, Peyret T, Gosselin NH, Marier JF, Hofmann AF, Shapiro D. Modeling and Experimental Studies of Obeticholic Acid Exposure and the Impact of Cirrhosis Stage. Clin Transl Sci. 2016 Dec;9(6):328-336. doi:10.1111/cts.12421"
  vignette <- "Edwards_2016_obeticholicAcid"

  # All 19 anatomical-space x chemical-form states are paper-mechanistic: the
  # model is an OCA-specific adaptation of the Molino 1986 CDCA physiological
  # PK model, whose seven "spaces" (systemic, portal, sinusoidal, liver, bile
  # duct, gallbladder, gut) do not map onto the canonical central/peripheral
  # names. Same treatment as the sibling bile-acid model Zuo_2016_UDCA.R.
  # OCA itself has no bile duct or gallbladder state: Figure 1 shows no
  # liver -> bile duct arrow in the OCA column (there is no "t20"), i.e. only
  # the conjugates are secreted into bile.
  paper_specific_compartments <- c(
    "systemic_oca", "portal_oca", "sinusoidal_oca", "liver_oca", "gut_oca",
    "systemic_goca", "portal_goca", "sinusoidal_goca", "liver_goca",
    "bileduct_goca", "gallbladder_goca", "gut_goca",
    "systemic_toca", "portal_toca", "sinusoidal_toca", "liver_toca",
    "bileduct_toca", "gallbladder_toca", "gut_toca"
  )

  units <- list(time = "h", dosing = "nmol", concentration = "nM")

  # Issue #482. Every state holds an amount in nmol; the model is solved on a
  # molar basis because the paper's residual error and LLOQ are reported in nM
  # ("Bioanalytical": OCA 420.6, glyco-OCA 477.7, tauro-OCA 527.8 g/mol).
  # verified = TRUE: the seven spaces and their contents are read off the
  # Figure 1 diagram and cross-checked against the Table 1 "Mean Simulated OCA
  # Distribution (% Nanomoles Total OCA)" block, which is itemised by exactly
  # these spaces.
  compartmentData <- list(
    depot            = list(analyte = "OCA", units = "nmol", specimen = "administration site", verified = TRUE),
    systemic_oca     = list(analyte = "OCA", units = "nmol", specimen = "plasma", verified = TRUE),
    portal_oca       = list(analyte = "OCA", units = "nmol", specimen = "whole blood", verified = TRUE),
    sinusoidal_oca   = list(analyte = "OCA", units = "nmol", specimen = "whole blood", verified = TRUE),
    liver_oca        = list(analyte = "OCA", units = "nmol", specimen = "tissue", verified = TRUE),
    gut_oca          = list(analyte = "OCA", units = "nmol", specimen = "administration site", verified = TRUE),
    systemic_goca    = list(analyte = "glyco-OCA", units = "nmol", specimen = "plasma", verified = TRUE),
    portal_goca      = list(analyte = "glyco-OCA", units = "nmol", specimen = "whole blood", verified = TRUE),
    sinusoidal_goca  = list(analyte = "glyco-OCA", units = "nmol", specimen = "whole blood", verified = TRUE),
    liver_goca       = list(analyte = "glyco-OCA", units = "nmol", specimen = "tissue", verified = TRUE),
    bileduct_goca    = list(analyte = "glyco-OCA", units = "nmol", specimen = "bile", verified = TRUE),
    gallbladder_goca = list(analyte = "glyco-OCA", units = "nmol", specimen = "bile", verified = TRUE),
    gut_goca         = list(analyte = "glyco-OCA", units = "nmol", specimen = "administration site", verified = TRUE),
    systemic_toca    = list(analyte = "tauro-OCA", units = "nmol", specimen = "plasma", verified = TRUE),
    portal_toca      = list(analyte = "tauro-OCA", units = "nmol", specimen = "whole blood", verified = TRUE),
    sinusoidal_toca  = list(analyte = "tauro-OCA", units = "nmol", specimen = "whole blood", verified = TRUE),
    liver_toca       = list(analyte = "tauro-OCA", units = "nmol", specimen = "tissue", verified = TRUE),
    bileduct_toca    = list(analyte = "tauro-OCA", units = "nmol", specimen = "bile", verified = TRUE),
    gallbladder_toca = list(analyte = "tauro-OCA", units = "nmol", specimen = "bile", verified = TRUE),
    gut_toca         = list(analyte = "tauro-OCA", units = "nmol", specimen = "administration site", verified = TRUE)
  )

  covariateData <- list(
    HEPIMP_MILD = list(
      description       = "Mild hepatic impairment (Child-Pugh Class A, score 5-6).",
      units             = "(binary)",
      type              = "binary",
      reference_category = "0 (normal hepatic function, or a non-mild impairment category)",
      notes             = "Child-Pugh classification, not NCI ODWG (Edwards 2016 Results, 'Hepatic impairment model development': 'Child-Pugh Score: Class A/Mild 5-6 points, Class B/Moderate 7-9 points, Class C/Severe 10-15 points'). Mutually exclusive with HEPIMP_MOD and HEPIMP_SEV; all three 0 selects the healthy-volunteer physiology. Selects the mild column of all four hepatic-impairment mechanisms in Supplementary Table S2.",
      source_name       = "Child-Pugh Class A"
    ),
    HEPIMP_MOD = list(
      description       = "Moderate hepatic impairment (Child-Pugh Class B, score 7-9).",
      units             = "(binary)",
      type              = "binary",
      reference_category = "0 (normal hepatic function, or a non-moderate impairment category)",
      notes             = "Child-Pugh classification, not NCI ODWG. Mutually exclusive with HEPIMP_MILD and HEPIMP_SEV. Selects the moderate column of Supplementary Table S2.",
      source_name       = "Child-Pugh Class B"
    ),
    HEPIMP_SEV = list(
      description       = "Severe hepatic impairment (Child-Pugh Class C, score 10-15).",
      units             = "(binary)",
      type              = "binary",
      reference_category = "0 (normal hepatic function, or a non-severe impairment category)",
      notes             = "Child-Pugh classification, not NCI ODWG. Mutually exclusive with HEPIMP_MILD and HEPIMP_MOD. Selects the severe column of Supplementary Table S2.",
      source_name       = "Child-Pugh Class C"
    ),
    MEAL_FLAG = list(
      description       = "Indicator that the current time falls within a post-prandial gallbladder-contraction window.",
      units             = "(binary)",
      type              = "binary",
      reference_category = "0 (no gallbladder contraction active)",
      notes             = "Time-varying covariate; must be supplied on every observation row. Set to 1 over [t_meal, t_meal + 1.5 h] for each standardized meal and 0 otherwise, per Edwards 2016 'Source data': 'Gallbladder contraction was assumed to last 90 min after the start of a meal.' Gates the gallbladder -> gut output rate f23; outside the window the gallbladder only fills (via f22) and does not empty.",
      source_name       = "meal consumption information"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 399L,
    n_studies       = 5L,
    age_range       = "adults >= 18 years; healthy-volunteer pool mean (SD) age 37.0 (9.8) years, hepatic-impairment cohort 55.0 (5.6) years",
    weight_range    = "healthy-volunteer pool mean (SD) 76.4 (11.8) kg; hepatic-impairment cohort 81.7 (16.9) kg",
    sex_female_pct  = 41,
    race_ethnicity  = "Healthy-volunteer pool (Study 1): 65.6% white, 32.5% black, 0.6% Asian, 1.3% other. Hepatic-impairment cohort (Study 2): 90.6% white, 3.1% black, 3.1% Asian, 3.1% other.",
    disease_state   = "Model development: healthy volunteers with normal hepatic function (Study 1, n = 160; 8,248 plasma samples) then subjects with Child-Pugh A/B/C cirrhosis plus normal-function controls (Study 2, n = 32; 928 samples). External validation: healthy volunteers (Studies 3 and 4, n = 24 and n = 160) and cirrhotic subjects with portal hypertension (Study 5 / PESTO, n = 23).",
    dose_range      = "5, 10 and 25 mg oral OCA; single dose and once-daily multiple dosing to steady state",
    notes           = "Sex percentage is derived from the two model-development cohorts, which were 59% and 72% male. The 22 structural parameters were estimated on the healthy-volunteer data and then held fixed while only the hepatic-impairment parameters were estimated. Population fit in Phoenix NLME v1.3 with Lindstrom-Bates FOCE. BLQ samples (38.2% / 9.4% / 24.4% for OCA / glyco-OCA / tauro-OCA in Study 1) were imputed to LLOQ/2."
  )

  ini({
    # ---- Blood, bile and gastrointestinal flows (L/h) --------------------
    # Supplemental Table 3 reports f3 and f4 as "Fixed"; they are inherited
    # physiological values from the base CDCA model (Molino 1986, reference 18).
    lf_por_sin  <- fixed(log(39.6));  label("Hepatic portal flow, portal -> sinusoidal, f3 (L/h)")     # Suppl Table 3
    lf_sys_sin  <- fixed(log(14.4));  label("Hepatic arterial flow, systemic -> sinusoidal, f4 (L/h)") # Suppl Table 3
    lf_bd_gb    <- log(0.856);        label("Flow from bile duct to gallbladder, f22 (L/h)")           # Suppl Table 3
    lf_bd_gut   <- log(7.29);         label("Flow from bile duct to gut, f24 (L/h)")                   # Suppl Table 3
    lk_gb_gut   <- fixed(log(1.2));   label("Rate of output from gallbladder to gut, f23 (1/h)")       # Suppl Table 3
    lkout       <- log(0.612);        label("Rate of fecal elimination of OCA, Kout (L/h)")            # Suppl Table 3

    # ---- Hepatic uptake and efflux (1/h) --------------------------------
    lt_sin_liv_oca  <- log(1698);        label("OCA transport rate, sinusoidal -> liver, t10 (1/h)")        # Suppl Table 3
    lt_sin_liv_goca <- log(1210);        label("Glyco-OCA transport rate, sinusoidal -> liver, t9 (1/h)")   # Suppl Table 3
    lt_sin_liv_toca <- log(1615);        label("Tauro-OCA transport rate, sinusoidal -> liver, t11 (1/h)")  # Suppl Table 3
    lt_liv_sin_oca  <- fixed(log(1.62)); label("OCA transport rate, liver -> sinusoidal, t13 (1/h)")        # Suppl Table 3
    lt_liv_sin_conj <- fixed(log(1.62)); label("Glyco- and tauro-OCA transport rate, liver -> sinusoidal, t12 (1/h)") # Suppl Table 3

    # ---- Biliary secretion of the conjugates (1/h) ----------------------
    # There is no OCA counterpart: Figure 1 has no liver -> bile duct arrow in
    # the OCA column, so unconjugated OCA is not secreted into bile.
    lt_liv_bd_goca <- log(7.44); label("Glyco-OCA transport rate, liver -> bile duct, t19 (1/h)") # Suppl Table 3
    lt_liv_bd_toca <- log(9.28); label("Tauro-OCA transport rate, liver -> bile duct, t21 (1/h)") # Suppl Table 3

    # ---- Intestinal absorption into the portal space (1/h) --------------
    lt_gut_por_oca  <- log(0.857); label("OCA rate of absorption, gut -> portal, t34 (1/h)")       # Suppl Table 3
    lt_gut_por_goca <- log(0.904); label("Glyco-OCA rate of absorption, gut -> portal, t33 (1/h)") # Suppl Table 3
    lt_gut_por_toca <- log(1.62);  label("Tauro-OCA rate of absorption, gut -> portal, t35 (1/h)") # Suppl Table 3

    # ---- Biotransformation (1/h) ----------------------------------------
    # Conjugation occurs in the liver; deconjugation by gut bacteria.
    lb_conj_gly    <- log(1.44);   label("OCA rate of conjugation with glycine, b15 (1/h)")   # Suppl Table 3
    lb_conj_tau    <- log(0.312);  label("OCA rate of conjugation with taurine, b16 (1/h)")   # Suppl Table 3
    lb_deconj_goca <- log(0.0431); label("Glyco-OCA rate of deconjugation to OCA, b36 (1/h)") # Suppl Table 3
    lb_deconj_toca <- log(0.0200); label("Tauro-OCA rate of deconjugation to OCA, b37 (1/h)") # Suppl Table 3

    # ---- Oral absorption -------------------------------------------------
    lka <- log(5.32); label("OCA first-order rate constant of oral absorption, Ka (1/h)") # Suppl Table 3

    # ---- Hepatic-impairment effects (Supplemental Table 4) ---------------
    # Log-scale multiplicative deviations, applied as exp(effect) per the
    # Supplemental Table 2 equations. Estimated on Study 2 while the 22
    # healthy-volunteer structural parameters were held fixed.
    e_hepimp_mild_uptake <- -0.132;  label("Effect of mild hepatic impairment on hepatic uptake of OCA and conjugates (log scale)")     # Suppl Table 4
    e_hepimp_mod_uptake  <- -1.86;   label("Effect of moderate hepatic impairment on hepatic uptake of OCA and conjugates (log scale)") # Suppl Table 4
    e_hepimp_sev_uptake  <- -2.37;   label("Effect of severe hepatic impairment on hepatic uptake of OCA and conjugates (log scale)")   # Suppl Table 4
    e_hepimp_mild_tauro  <- 0.00481; label("Effect of mild hepatic impairment on OCA tauro-conjugation (log scale)")     # Suppl Table 4
    e_hepimp_mod_tauro   <- 1.05;    label("Effect of moderate hepatic impairment on OCA tauro-conjugation (log scale)") # Suppl Table 4
    e_hepimp_sev_tauro   <- 1.56;    label("Effect of severe hepatic impairment on OCA tauro-conjugation (log scale)")   # Suppl Table 4

    # ---- Between-subject variability ------------------------------------
    # Supplemental Table 3 reports BSV on only two parameters, as %CV from
    # Phoenix NLME. Converted to log-scale variance by omega^2 = log(CV^2 + 1):
    #   Ka  195%  -> log(1.95^2 + 1)  = 1.5692
    #   f22  78.1% -> log(0.781^2 + 1) = 0.4762
    etalka    ~ 1.5692  # Suppl Table 3, BSV 195% CV on Ka
    etalf_bd_gb ~ 0.4762 # Suppl Table 3, BSV 78.1% CV on flow from bile duct to gallbladder

    # ---- Residual error ---------------------------------------------------
    # Healthy-volunteer values (Supplemental Table 3). The hepatic-impairment
    # model re-estimated a separate residual-error set (Supplemental Table 4:
    # proportional 122 / 112 / 123%, additive 0.993 / 0.273 / 0.532 nM); see the
    # vignette's Assumptions and deviations section.
    propSd         <- 0.880; label("Proportional residual error, OCA")       # Suppl Table 3, 88.0%
    addSd          <- 0.546; label("Additive residual error, OCA (nM)")      # Suppl Table 3
    propSd_Cc_goca <- 0.626; label("Proportional residual error, glyco-OCA") # Suppl Table 3, 62.6%
    addSd_Cc_goca  <- 0.675; label("Additive residual error, glyco-OCA (nM)") # Suppl Table 3
    propSd_Cc_toca <- 0.716; label("Proportional residual error, tauro-OCA") # Suppl Table 3, 71.6%
    addSd_Cc_toca  <- 0.469; label("Additive residual error, tauro-OCA (nM)") # Suppl Table 3
  })

  model({
    # ================= Physiological compartment volumes (L) =============
    # Printed in the Volumes legend inside the Figure 1 panel (p. 330). These
    # are the values the paper describes in Methods as "fixed to physiological
    # values from the base CDCA model" (Molino 1986, reference 18); the figure
    # panel reproduces them, so no upstream lookup is required.
    v_sys <- 2.46 # Systemic circulation  -- Figure 1 legend
    v_por <- 0.42 # Portal circulation    -- Figure 1 legend
    v_sin <- 0.12 # Sinusoidal circulation -- Figure 1 legend
    v_liv <- 0.95 # Liver                 -- Figure 1 legend
    v_bd  <- 0.02 # Bile duct             -- Figure 1 legend
    v_gb  <- 0.02 # Gallbladder           -- Figure 1 legend
    v_gut <- 0.90 # Intestine (gut)       -- Figure 1 legend

    # ============ Hepatic-impairment physiology (Suppl Table 2) ==========
    # Fixed multipliers taken from the Simcyp cirrhotic library via Johnson
    # 2010 (Edwards 2016 reference 11), applied for Child-Pugh A / B / C.
    # The three indicators are mutually exclusive, so each multiplier
    # collapses to 1 (or 0 for the shunt) in healthy subjects.
    f4_mult <- 1 + HEPIMP_MILD * (1.408 - 1) + HEPIMP_MOD * (1.625 - 1) + HEPIMP_SEV * (1.915 - 1)
    f3_mult <- 1 + HEPIMP_MILD * (0.910 - 1) + HEPIMP_MOD * (0.635 - 1) + HEPIMP_SEV * (0.554 - 1)
    vl_mult <- 1 + HEPIMP_MILD * (0.891 - 1) + HEPIMP_MOD * (0.710 - 1) + HEPIMP_SEV * (0.610 - 1)

    uptake_eff <- HEPIMP_MILD * e_hepimp_mild_uptake +
      HEPIMP_MOD * e_hepimp_mod_uptake +
      HEPIMP_SEV * e_hepimp_sev_uptake
    tauro_eff <- HEPIMP_MILD * e_hepimp_mild_tauro +
      HEPIMP_MOD * e_hepimp_mod_tauro +
      HEPIMP_SEV * e_hepimp_sev_tauro

    # ===================== Individual parameters =========================
    f3   <- exp(lf_por_sin) * f3_mult   # portal -> sinusoidal
    f4   <- exp(lf_sys_sin) * f4_mult   # systemic -> sinusoidal (hepatic artery)
    # Portal-systemic shunt f2: "the coefficient for portal to sinusoidal flow
    # was progressively decreased ... and was matched by a progressive increase
    # in flow from the portal to systemic circulation of equal magnitude. The
    # latter flow does not occur in healthy individuals." (Results, Hepatic
    # impairment). So f2 = f3_healthy - f3_impaired, and is 0 when healthy.
    f2   <- exp(lf_por_sin) - f3
    # f1 (systemic -> portal) and f5 (sinusoidal -> systemic) are shown in
    # Figure 1 but are not tabulated; they are fixed by blood-flow continuity
    # on the portal and sinusoidal spaces respectively:
    #   portal:     f1 = f2 + f3  = 39.6 L/h in every impairment stratum
    #   sinusoidal: f5 = f3 + f4
    f1   <- f2 + f3
    f5   <- f3 + f4
    f22  <- exp(lf_bd_gb + etalf_bd_gb)
    f24  <- exp(lf_bd_gut)
    f23  <- exp(lk_gb_gut)
    kout <- exp(lkout)

    # Hepatic uptake is reduced in cirrhosis; sinusoidal efflux is not.
    t10 <- exp(lt_sin_liv_oca + uptake_eff)
    t9  <- exp(lt_sin_liv_goca + uptake_eff)
    t11 <- exp(lt_sin_liv_toca + uptake_eff)
    t13 <- exp(lt_liv_sin_oca)
    t12 <- exp(lt_liv_sin_conj)

    t19 <- exp(lt_liv_bd_goca)
    t21 <- exp(lt_liv_bd_toca)
    t34 <- exp(lt_gut_por_oca)
    t33 <- exp(lt_gut_por_goca)
    t35 <- exp(lt_gut_por_toca)

    b15 <- exp(lb_conj_gly)
    b16 <- exp(lb_conj_tau + tauro_eff) # preferential tauro-conjugation in cirrhosis
    b36 <- exp(lb_deconj_goca)
    b37 <- exp(lb_deconj_toca)

    ka <- exp(lka + etalka)

    # Functional liver volume shrinks with cirrhosis. It does not enter the
    # ODEs (every liver flux is a first-order rate constant on amount), only
    # the reported liver concentration - which is why Table 1 shows liver
    # exposure rising far less than systemic exposure.
    v_liv_i <- v_liv * vl_mult

    # Gallbladder empties only during the 90-minute post-prandial window.
    gb_out <- f23 * MEAL_FLAG

    # ===================== Concentrations (nM) ===========================
    c_sys_oca  <- systemic_oca / v_sys
    c_por_oca  <- portal_oca / v_por
    c_sin_oca  <- sinusoidal_oca / v_sin
    c_gut_oca  <- gut_oca / v_gut

    c_sys_goca <- systemic_goca / v_sys
    c_por_goca <- portal_goca / v_por
    c_sin_goca <- sinusoidal_goca / v_sin
    c_bd_goca  <- bileduct_goca / v_bd

    c_sys_toca <- systemic_toca / v_sys
    c_por_toca <- portal_toca / v_por
    c_sin_toca <- sinusoidal_toca / v_sin
    c_bd_toca  <- bileduct_toca / v_bd

    # ========================= ODE system =================================
    # Amounts in nmol. Flows f (L/h) act on the source-space concentration;
    # transport t, biotransformation b, Ka and f23 (1/h) act on amount.
    # Topology transcribed from the Figure 1 diagram (p. 330).
    d/dt(depot) <- -ka * depot

    # ---- OCA (no bile duct / gallbladder state; see Figure 1) ------------
    d/dt(systemic_oca) <- f5 * c_sin_oca + f2 * c_por_oca -
      f1 * c_sys_oca - f4 * c_sys_oca
    d/dt(portal_oca) <- f1 * c_sys_oca - f3 * c_por_oca - f2 * c_por_oca +
      t34 * gut_oca
    d/dt(sinusoidal_oca) <- f3 * c_por_oca + f4 * c_sys_oca - f5 * c_sin_oca -
      t10 * sinusoidal_oca + t13 * liver_oca
    d/dt(liver_oca) <- t10 * sinusoidal_oca - t13 * liver_oca -
      b15 * liver_oca - b16 * liver_oca
    d/dt(gut_oca) <- ka * depot + b36 * gut_goca + b37 * gut_toca -
      t34 * gut_oca - kout * c_gut_oca

    # ---- Glyco-OCA -------------------------------------------------------
    d/dt(systemic_goca) <- f5 * c_sin_goca + f2 * c_por_goca -
      f1 * c_sys_goca - f4 * c_sys_goca
    d/dt(portal_goca) <- f1 * c_sys_goca - f3 * c_por_goca - f2 * c_por_goca +
      t33 * gut_goca
    d/dt(sinusoidal_goca) <- f3 * c_por_goca + f4 * c_sys_goca - f5 * c_sin_goca -
      t9 * sinusoidal_goca + t12 * liver_goca
    d/dt(liver_goca) <- t9 * sinusoidal_goca - t12 * liver_goca +
      b15 * liver_oca - t19 * liver_goca
    d/dt(bileduct_goca) <- t19 * liver_goca - f22 * c_bd_goca - f24 * c_bd_goca
    d/dt(gallbladder_goca) <- f22 * c_bd_goca - gb_out * gallbladder_goca
    d/dt(gut_goca) <- f24 * c_bd_goca + gb_out * gallbladder_goca -
      t33 * gut_goca - b36 * gut_goca

    # ---- Tauro-OCA -------------------------------------------------------
    d/dt(systemic_toca) <- f5 * c_sin_toca + f2 * c_por_toca -
      f1 * c_sys_toca - f4 * c_sys_toca
    d/dt(portal_toca) <- f1 * c_sys_toca - f3 * c_por_toca - f2 * c_por_toca +
      t35 * gut_toca
    d/dt(sinusoidal_toca) <- f3 * c_por_toca + f4 * c_sys_toca - f5 * c_sin_toca -
      t11 * sinusoidal_toca + t12 * liver_toca
    d/dt(liver_toca) <- t11 * sinusoidal_toca - t12 * liver_toca +
      b16 * liver_oca - t21 * liver_toca
    d/dt(bileduct_toca) <- t21 * liver_toca - f22 * c_bd_toca - f24 * c_bd_toca
    d/dt(gallbladder_toca) <- f22 * c_bd_toca - gb_out * gallbladder_toca
    d/dt(gut_toca) <- f24 * c_bd_toca + gb_out * gallbladder_toca -
      t35 * gut_toca - b37 * gut_toca

    # ======================== Observations ================================
    # Plasma concentrations are the systemic-circulation concentrations.
    Cc      <- c_sys_oca
    Cc_goca <- c_sys_goca
    Cc_toca <- c_sys_toca

    # Total OCA (Results, Source data): "total OCA concentrations were
    # calculated as the sum of OCA, glyco-OCA, and tauro-OCA" on a molar
    # basis. Liver total OCA drives the Table 1 hepatic-exposure column.
    Cc_total    <- c_sys_oca + c_sys_goca + c_sys_toca
    Cliver      <- (liver_oca + liver_goca + liver_toca) / v_liv_i

    Cc      ~ add(addSd) + prop(propSd)
    Cc_goca ~ add(addSd_Cc_goca) + prop(propSd_Cc_goca)
    Cc_toca ~ add(addSd_Cc_toca) + prop(propSd_Cc_toca)
  })
}
