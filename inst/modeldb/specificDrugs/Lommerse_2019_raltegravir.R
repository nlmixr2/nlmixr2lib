Lommerse_2019_raltegravir <- function() {
  description <- "Integrated maternal-neonatal two-compartment first-order-absorption population PK model of oral raltegravir, coupling the maternal and fetal central compartments via a very fast intercompartmental clearance (1000 L/h) during pregnancy to enforce instantaneous placental equilibrium, and decoupling at birth (time t=0 in the model). Neonatal apparent clearance rises from nil at birth to CL_max (9.44 L/h at 25 kg) with a first-order maturation rate constant CL_tau (11.3 1/year, 90% mature by ~11 weeks); neonatal absorption rate constant rises from KA_base (0.0915 1/h) to KA_max (0.43 1/h) with a first-order maturation rate constant KA_tau (63.2 1/year, 90% mature by ~12 days). Neonate CL, Q, and volumes are allometrically scaled with fixed exponents 0.75 and 1.0 to a reference weight of 25 kg. Maternal disposition parameters (V2 3.52 L, V3 27 L, CL 9.73 L/h, Q 0.866 L/h; all at 25 kg reference) are fixed from the Rizk 2015 pediatric popPK (ref [12]); maternal KA (0.175 1/h) and bioavailability (F 0.517) are estimated. IIV is on neonate CL and KA and maternal F; residual error is combined additive (11.9 nM) and proportional (54%) (Lommerse 2019)."
  reference <- "Lommerse J, Clarke D, Kerbusch T, Merdjan H, Witjes H, Teppler H, Mirochnick M, Acosta EP, Wenning L, Nachman S, Chain A. Maternal-Neonatal Raltegravir Population Pharmacokinetics Modeling: Implications for Initial Neonatal Dosing. CPT Pharmacometrics Syst Pharmacol. 2019;8(9):643-653. doi:10.1002/psp4.12443"
  vignette <- "Lommerse_2019_raltegravir"
  units <- list(time = "hour", dosing = "mg", concentration = "nM")

  paper_specific_compartments <- c(
    "depot_mother", "central_mother", "peripheral_mother",
    "depot_neonate", "central_neonate", "peripheral_neonate"
  )

  covariateData <- list(
    WT = list(
      description        = "Neonate body weight (time-varying, kg)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying per observation. Drives the allometric scaling of neonate CL, Q (exponent 0.75) and central/peripheral V (exponent 1.0) to a reference weight of 25 kg (Table 2 caption). Postnatal weight growth in the paper's simulations follows the empirical fit `BW(kg) = 2.935 + 8.909 * (1 - exp(-1.103 * PNA_years))` (Methods). Neonate birth weight in the analysis dataset ranged 2.2-3.4 kg for raltegravir-exposed and 2.2-5.3 kg for raltegravir-unexposed neonates (Table 1). Maternal body weight is not modelled as a covariate (Table 1 lists it as 'Unknown'); the maternal fixed disposition is applied at a typical maternal weight of 60 kg baked into the model with the same allometric exponents.",
      source_name        = "BW"
    ),
    PNA = list(
      description        = "Postnatal age (chronological time since birth, months; time-varying)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Drives the neonate CL and KA first-order maturation kernels. Paper uses PNA in YEARS (`CL_tau = 11.3 1/year`, `KA_tau = 63.2 1/year`); canonical PNA units are months, so the maturation equations are reparameterised inline as `1 - exp(-cl_tau * PNA/12)` and `1 - exp(-ka_tau * PNA/12)`. During pregnancy the neonate is a fetus and the CL/KA maturation is inactive; the paper reports 'no elimination capacity at birth' (CL_base = 0 fixed), so PNA >= 0 by construction and `pmax(PNA, 0)` clamps any negative in-utero time to zero.",
      source_name        = "PNA"
    )
  )

  population <- list(
    species              = "human",
    n_subjects           = 104L,
    n_subjects_mothers   = 19L,
    n_subjects_neonates  = 61L,
    n_subjects_infants   = 24L,
    n_studies            = 3L,
    n_observations       = 759L,
    age_range            = "Mothers: adult, gravid at 3rd trimester or delivery. Raltegravir-unexposed neonates: 0-11 days at PK sampling (P1110 Cohort 1) or 0-6 weeks (P1110 Cohort 2). Raltegravir-exposed neonates: 0-11 days (P1110) or 0-2 days (P1097). Infants: 5 weeks to <2.4 years (P1066).",
    weight_range         = "Raltegravir-unexposed neonates 2.2-5.3 kg; raltegravir-exposed neonates 2.2-4.1 kg; infants 3.7-14 kg; mothers unknown (Table 1 reports 'Unknown').",
    sex_female_pct       = 40,
    race_ethnicity       = "Predominantly African American / Black (majority per Methods 'Covariate analysis', not enumerated).",
    disease_state        = "HIV-1-exposed neonates (with and without in utero raltegravir exposure), HIV-1-infected infants and toddlers, and HIV-1-infected pregnant women.",
    dose_range           = "Mothers: 400 mg twice-daily raltegravir during pregnancy (IMPAACT P1097; last dose within 48 h of delivery). Raltegravir-unexposed neonates: 2 x 3 mg/kg single doses one week apart (P1110 Cohort 1) or the ramped 6-week regimen 1.5 mg/kg QD (days 1-7) -> 3 mg/kg BID (days 8-28) -> 6 mg/kg BID (days 28-42) (P1110 Cohort 2). Raltegravir-exposed neonates: 1.5 mg/kg single dose (P1110). Infants: 6 mg/kg BID granules for oral suspension (P1066).",
    regions              = "United States, Brazil, South Africa, Thailand (IMPAACT network sites).",
    notes                = "Data pooled from IMPAACT P1110 (raltegravir-unexposed and exposed term neonates 0-6 weeks), P1066 (infants and toddlers 4 weeks to <2 years), and P1097 (19 mother-infant pairs; mothers received raltegravir during pregnancy). Cohort composition and demographics from Table 1. Model estimated in NONMEM (ICON Inc.), SAS 9.4 for data assembly, R 3.1.3 for post-processing. Bootstrap N=1000 stratified by study, cohort, and in-utero-exposure status; 237 minimisation-terminated and 3 boundary runs omitted (Table 2 bootstrap columns)."
  )

  ini({
    # =====================================================================
    # NEONATE PK -- reference weight 25 kg (Table 2 caption); PNA in years for
    # the maturation kernels (paper's original parameterisation), converted
    # inline to months (canonical PNA units) in model() via PNA/12.
    # =====================================================================
    lvc_neonate  <- log(7.04);   label("Neonate central volume of distribution at 25 kg (V2, L)")            # Table 2, V2 = 7.04 L
    lvp_neonate  <- log(10.3);   label("Neonate peripheral volume of distribution at 25 kg (V3, L)")         # Table 2, V3 = 10.3 L
    lcl_max      <- log(9.44);   label("Neonate maximum apparent clearance at full maturation, 25 kg (CL_max, L/h)")   # Table 2, CL_max = 9.44 L/h
    lq_neonate   <- log(0.786);  label("Neonate intercompartmental clearance at 25 kg (Q, L/h)")             # Table 2, Q = 0.786 L/h
    lka_max      <- log(0.43);   label("Neonate maximum absorption rate constant at full maturation (KA_max, 1/h)")    # Table 2, KA_max = 0.43 1/h
    lka_base     <- log(0.0915); label("Neonate absorption rate constant at birth (KA_base, 1/h)")           # Table 2, KA_base = 0.0915 1/h
    cl_tau       <- 11.3;        label("Neonate CL maturation rate constant (1/year); 90% mature by ~11 weeks")        # Table 2, CL_tau = 11.3 1/year
    ka_tau       <- 63.2;        label("Neonate KA maturation rate constant (1/year); 90% mature by ~12 days")         # Table 2, KA_tau = 63.2 1/year
    # CL_base (neonate CL at birth) is fixed at 0 in the paper (Table 2, CL_base = 0);
    # this is baked into the maturation kernel `cl_neonate = cl_max * (1 - exp(-cl_tau * PNA_years))`,
    # which evaluates to 0 at PNA=0 by construction. No explicit CL_base parameter is needed.
    # F4 (neonate bioavailability after birth) fixed to 1 as the reference anchor for the
    # relative-bioavailability parameterisation.
    lfdepot_neonate <- fixed(log(1)); label("Neonate oral bioavailability (F4, fraction)")       # Table 2, F4 fixed = 1

    # =====================================================================
    # MATERNAL PK -- Rizk 2015 pediatric popPK (ref [12]) values at the same
    # 25 kg reference. All maternal disposition parameters are fixed; only
    # KA and F are estimated to inform the mother-to-fetus transfer.
    # =====================================================================
    lcl_mother    <- fixed(log(9.73));  label("Maternal apparent clearance at 25 kg, from Rizk 2015 (CL, L/h)")     # Table 2, CL fixed = 9.73 L/h
    lvc_mother    <- fixed(log(3.52));  label("Maternal central volume at 25 kg, from Rizk 2015 (V2, L)")           # Table 2, V2 fixed = 3.52 L
    lvp_mother    <- fixed(log(27));    label("Maternal peripheral volume at 25 kg, from Rizk 2015 (V3, L)")        # Table 2, V3 fixed = 27 L
    lq_mother     <- fixed(log(0.866)); label("Maternal intercompartmental clearance at 25 kg, from Rizk 2015 (Q, L/h)") # Table 2, Q fixed = 0.866 L/h
    lka_mother    <- log(0.175);        label("Maternal first-order absorption rate constant (KA, 1/h)")                  # Table 2, KA = 0.175 1/h
    lfdepot_mother <- log(0.517);       label("Maternal oral bioavailability relative to neonate F=1 anchor (F, fraction)") # Table 2, F = 0.517

    # =====================================================================
    # PLACENTAL COUPLING (in utero) -- set to 1000 L/h to enforce
    # instantaneous maternal-fetal equilibrium during pregnancy; the
    # coupling is gated OFF at and after birth (t >= 0).
    # =====================================================================
    q_link_pregnant <- fixed(1000); label("Placental intercompartmental clearance during pregnancy (L/h; instant equilibrium)") # Methods: "the intercompartmental CL linking the maternal and fetal central compartments was set to 1,000 L/hour"

    # =====================================================================
    # ALLOMETRIC EXPONENTS (all fixed per Methods)
    # =====================================================================
    e_wt_cl_q  <- fixed(0.75); label("Shared allometric exponent on CL and Q (unitless; Holford 1996)")  # Methods: fixed to 0.75
    e_wt_vc_vp <- fixed(1.0);  label("Shared allometric exponent on Vc and Vp (unitless; isometric)")    # Methods: fixed to 1.0

    # =====================================================================
    # IIV -- paper reports the magnitudes as fractions in Table 2 (matching
    # the "IIV on CL and KA was 33% and 20%" text). Encoded as log-normal
    # variance = log(1 + CV^2) so the omega on the log scale matches the
    # paper's "33%" and "20%" and "31%" verbal magnitudes.
    # IIV is applied to the ASYMPTOTIC anchors lcl_max and lka_max (the
    # fully-matured typical values); the corresponding maturation kernel
    # remains a shared population-level function of PNA.
    # =====================================================================
    etalcl_max        ~ 0.10336  # log(1 + 0.33^2);  Table 2 IIV on CL = 0.33  (paper text 33% CV)
    etalka_max        ~ 0.03772  # log(1 + 0.196^2); Table 2 IIV on KA = 0.196 (paper text 20% CV)
    etalfdepot_mother ~ 0.09239  # log(1 + 0.311^2); Table 2 IIV on F  = 0.311 (maternal-F)

    # =====================================================================
    # RESIDUAL ERROR -- combined additive (nM) + proportional (fraction).
    # The paper reports ONE shared residual pair applied to both maternal
    # and neonatal plasma observations (Table 2 residual-error block).
    # rxode2/nlmixr2 requires each endpoint to have its own residual-error
    # parameters, so `propSd_Cmother` / `propSd_Cneonate` and
    # `addSd_Cmother` / `addSd_Cneonate` are declared with IDENTICAL
    # numeric values from Table 2 to preserve the shared-residual fit.
    # =====================================================================
    propSd_Cmother  <- 0.54;  label("Proportional residual error on maternal plasma (fraction)")  # Table 2, RUV-prop = 0.54 (shared)
    addSd_Cmother   <- 11.9;  label("Additive residual error on maternal plasma (nM)")            # Table 2, RUV-add = 11.9 nM (shared)
    propSd_Cneonate <- 0.54;  label("Proportional residual error on neonate plasma (fraction)")   # Table 2, RUV-prop = 0.54 (shared)
    addSd_Cneonate  <- 11.9;  label("Additive residual error on neonate plasma (nM)")             # Table 2, RUV-add = 11.9 nM (shared)
  })

  model({
    # ---------------------------------------------------------------
    # 1. Time / age transforms
    # ---------------------------------------------------------------
    # Paper anchors t = 0 at birth. Pre-birth (in utero): t < 0.
    # Post-birth: t >= 0. PNA is postnatal age in months (canonical); the
    # `(PNA >= 0) * PNA` clamp forces the maturation kernels to 0 pre-birth
    # (consistent with "no elimination capacity at birth"). `pmax()` is not
    # available inside rxode2 model bodies (variadic-argument restriction);
    # the multiplicative clamp is the standard rxode2 idiom.
    pna_months_clamped <- (PNA >= 0) * PNA
    pna_years          <- pna_months_clamped / 12  # convert canonical months -> years for the paper's kernels

    # ---------------------------------------------------------------
    # 2. Neonate individual PK parameters (allometric to WT, matured with PNA)
    # ---------------------------------------------------------------
    matur_cl_neonate <- 1 - exp(-cl_tau * pna_years)   # 0 at PNA=0; asymptote 1
    matur_ka_neonate <- 1 - exp(-ka_tau * pna_years)   # 0 at PNA=0; asymptote 1

    cl_neonate <- exp(lcl_max + etalcl_max) * (WT / 25)^e_wt_cl_q  * matur_cl_neonate
    vc_neonate <- exp(lvc_neonate)          * (WT / 25)^e_wt_vc_vp
    vp_neonate <- exp(lvp_neonate)          * (WT / 25)^e_wt_vc_vp
    q_neonate  <- exp(lq_neonate)           * (WT / 25)^e_wt_cl_q
    # KA blends KA_base -> KA_max with the IIV applied to the asymptotic
    # (KA_max) anchor; KA_base is a shared population-level constant.
    ka_neo_base_i <- exp(lka_base)
    ka_neo_max_i  <- exp(lka_max  + etalka_max)
    ka_neonate    <- ka_neo_base_i + (ka_neo_max_i - ka_neo_base_i) * matur_ka_neonate

    # ---------------------------------------------------------------
    # 3. Maternal individual PK parameters
    # ---------------------------------------------------------------
    # Table 1 does not report maternal weight; fix at a typical 60 kg
    # maternal weight for the same fixed allometric exponents.
    wt_mother  <- 60.0
    cl_mother  <- exp(lcl_mother) * (wt_mother / 25)^e_wt_cl_q
    vc_mother  <- exp(lvc_mother) * (wt_mother / 25)^e_wt_vc_vp
    vp_mother  <- exp(lvp_mother) * (wt_mother / 25)^e_wt_vc_vp
    q_mother   <- exp(lq_mother)  * (wt_mother / 25)^e_wt_cl_q
    ka_mother  <- exp(lka_mother)
    f_mother   <- exp(lfdepot_mother + etalfdepot_mother)

    # ---------------------------------------------------------------
    # 4. Placental coupling gate: 1000 L/h in utero (t < 0), 0 at/after birth.
    # ---------------------------------------------------------------
    q_link                <- q_link_pregnant * (t < 0)
    k_mother_to_neonate   <- q_link / vc_mother
    k_neonate_to_mother   <- q_link / vc_neonate

    # ---------------------------------------------------------------
    # 5. Micro-constants
    # ---------------------------------------------------------------
    kel_neonate <- cl_neonate / vc_neonate
    k12_neonate <- q_neonate  / vc_neonate
    k21_neonate <- q_neonate  / vp_neonate
    kel_mother  <- cl_mother  / vc_mother
    k12_mother  <- q_mother   / vc_mother
    k21_mother  <- q_mother   / vp_mother

    # ---------------------------------------------------------------
    # 6. ODE system
    # ---------------------------------------------------------------
    d/dt(depot_mother)       <- -ka_mother * depot_mother
    d/dt(central_mother)     <-  ka_mother * depot_mother -
                                  kel_mother * central_mother -
                                  k12_mother * central_mother + k21_mother * peripheral_mother -
                                  k_mother_to_neonate * central_mother +
                                  k_neonate_to_mother * central_neonate
    d/dt(peripheral_mother)  <-  k12_mother * central_mother - k21_mother * peripheral_mother
    d/dt(depot_neonate)      <- -ka_neonate * depot_neonate
    d/dt(central_neonate)    <-  ka_neonate * depot_neonate -
                                  kel_neonate * central_neonate -
                                  k12_neonate * central_neonate + k21_neonate * peripheral_neonate +
                                  k_mother_to_neonate * central_mother -
                                  k_neonate_to_mother * central_neonate
    d/dt(peripheral_neonate) <-  k12_neonate * central_neonate - k21_neonate * peripheral_neonate

    # ---------------------------------------------------------------
    # 7. Bioavailability
    # ---------------------------------------------------------------
    f(depot_mother)  <- f_mother
    f(depot_neonate) <- exp(lfdepot_neonate)  # fixed = 1

    # ---------------------------------------------------------------
    # 8. Observations -- plasma concentrations in nM. Doses are in mg and
    # volumes in L, so central / vc is in mg/L = ug/mL; multiplied by
    # 1000/0.4444 gives nM using the paper's conversion factor
    # (0.4444 ng/mL per nM = raltegravir MW 444.4 g/mol / 1000).
    # Output names Cmother / Cneonate are paper-named (single-token,
    # analogous to Ccsf) so they are exempt from the Cc_<metab> pattern.
    # ---------------------------------------------------------------
    Cmother  <- (central_mother  / vc_mother)  * 1000 / 0.4444
    Cneonate <- (central_neonate / vc_neonate) * 1000 / 0.4444

    Cmother  ~ add(addSd_Cmother)  + prop(propSd_Cmother)
    Cneonate ~ add(addSd_Cneonate) + prop(propSd_Cneonate)
  })
}
