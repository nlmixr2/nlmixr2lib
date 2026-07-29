Li_2014_penicillinG_cattle <- function() {
  description <- "Preclinical (cattle). Three-compartment population pharmacokinetic model for penicillin G in cattle, with four parallel first-order absorption depots covering intramuscular penicillin sodium, intramuscular procaine penicillin, subcutaneous procaine penicillin, and oral procaine penicillin (the oral depot feeds the liver compartment directly), plus separate liver and kidney tissue compartments connected to the central compartment by inter-compartmental clearance; pooled meta-analysis of 100 cattle from 30 published studies and FARAD records (Li 2014)."
  reference   <- "Li M, Gehring R, Tell L, Baynes R, Huang Q, Riviere JE. Interspecies mixed-effect pharmacokinetic modeling of penicillin G in cattle and swine. Antimicrob Agents Chemother. 2014;58(8):4495-4503. doi:10.1128/AAC.02806-14"
  vignette    <- "Li_2014_penicillinG"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used as a power-model covariate on Vp2 (second peripheral volume) and on the liver inter-compartmental clearance, per Li 2014 Table 2 (covariate factors). Form: P_i = P_pop * (WT/300)^theta. Reference WT (300 kg) is a rounded midrange of the cattle dataset (Li 2014 Table 1 weights 43.5-633 kg); the paper did not report the exact normalisation weight. Effect sizes (theta_1 = -0.005 on Vp2, theta_2 = 0.008 on liver clearance) are very small, so the exact reference is numerically inconsequential.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used as a power-model covariate on Vc (central volume) per Li 2014 Table 2 covariate factor theta_3 = -0.011. Form: Vc_i = Vc_pop * (AGE/1)^(-0.011) * exp(eta). Reference AGE (1 year) is a rounded midrange of Li 2014 Table 1 ages (0.01-5.8 years); the paper did not report the exact normalisation age. The effect coefficient is very small so the reference value is numerically inconsequential.",
      source_name        = "AGE"
    )
  )

  population <- list(
    species        = "cattle (Bos taurus); mixed steers, heifers, cows, and calves pooled from 30 published studies",
    n_subjects     = 100L,
    n_studies      = 30L,
    age_range      = "0.01-5.8 years (calves through adult cattle, pooled across studies)",
    weight_range   = "43.5-633 kg (pooled across studies; not every study reported weights)",
    sex_female_pct = NA_real_,
    disease_state  = "Healthy animals only (animals with various diseased conditions were excluded; Li 2014 Methods)",
    dose_range     = "Penicillin sodium 1-30 mg/kg IV or IM; procaine penicillin 1.5-66 mg/kg IM/SC/PO (Li 2014 Table 1)",
    regions        = "Published cattle PK studies (US and international); FARAD records",
    n_observations = 368L,
    notes          = "Pooled plasma/serum data (368 concentrations from 100 cattle) plus 26 liver and 13 kidney tissue residue concentrations. Cattle muscle data were too sparse to model. Routes carried in the dosing dataset via the rxode2 cmt column: 'central' for IV; 'depot1' for intramuscular penicillin sodium; 'depot2' for intramuscular procaine penicillin; 'depot3' for subcutaneous procaine penicillin; 'depot4' for oral procaine penicillin (which Li 2014 Methods assumed is absorbed directly into the liver compartment, so the depot4 mass flux feeds liver rather than central). The companion swine model from this paper is modellib('Li_2014_penicillinG_swine')."
  )

  ini({
    # Structural parameters from Li 2014 Table 2 (cattle, final PK model).
    # Volumes in litres; clearances in L/h; absorption rate constants in 1/h;
    # bioavailabilities as fractions (paper Table 2 reports %, converted here).
    lvc       <- log(3.45);  label("Central volume of distribution (Vc, L)")                                   # Li 2014 Table 2 V1
    lvp       <- log(20.8);  label("First peripheral volume of distribution (Vp, L)")                          # Li 2014 Table 2 V2
    lvp2      <- log(30.3);  label("Second peripheral volume of distribution (Vp2, L)")                        # Li 2014 Table 2 V3
    lv_liver  <- log(11.2);  label("Liver compartment volume of distribution (V_liver, L)")                    # Li 2014 Table 2 VL
    lv_kidney <- log(12.5);  label("Kidney compartment volume of distribution (V_kidney, L)")                  # Li 2014 Table 2 VK
    lcl       <- log(105);   label("Central clearance (CL, L/h)")                                              # Li 2014 Table 2 CL1
    lq        <- log(24.8);  label("Inter-compartmental clearance central <-> peripheral1 (Q, L/h)")           # Li 2014 Table 2 CL2
    lq2       <- log(9.65);  label("Inter-compartmental clearance central <-> peripheral2 (Q2, L/h)")          # Li 2014 Table 2 CL3
    lq_liver  <- log(21.5);  label("Inter-compartmental clearance central <-> liver (Q_liver, L/h)")           # Li 2014 Table 2 CLL
    lq_kidney <- log(21.5);  label("Inter-compartmental clearance central <-> kidney (Q_kidney, L/h)")         # Li 2014 Table 2 CLK
    lka_im_na  <- log(0.31); label("First-order absorption rate constant after IM penicillin sodium (Ka, 1/h)")     # Li 2014 Table 2 Kim1
    lka_im_pro <- log(0.22); label("First-order absorption rate constant after IM procaine penicillin (Ka, 1/h)")   # Li 2014 Table 2 Kim2
    lka_sc     <- log(0.56); label("First-order absorption rate constant after SC procaine penicillin (Ka, 1/h)")   # Li 2014 Table 2 Ksc
    lka_oral   <- log(0.38); label("First-order absorption rate constant after oral procaine penicillin (Ka, 1/h)") # Li 2014 Table 2 Kpo
    lfdepot_im_na  <- log(0.801); label("Bioavailability of IM penicillin sodium (F, fraction)")        # Li 2014 Table 2 Fim1 (80.1%)
    lfdepot_im_pro <- log(0.913); label("Bioavailability of IM procaine penicillin (F, fraction)")      # Li 2014 Table 2 Fim2 (91.3%)
    lfdepot_sc     <- log(0.752); label("Bioavailability of SC procaine penicillin (F, fraction)")      # Li 2014 Table 2 Fsc  (75.2%)
    lfdepot_oral   <- log(0.424); label("Bioavailability of oral procaine penicillin (F, fraction)")    # Li 2014 Table 2 Fpo  (42.4%)

    # Covariate effects (Li 2014 Table 2 covariate factors). Power model:
    # P_i = P_pop * (cov/ref)^theta * exp(eta). The paper did not report the
    # normalisation value; rounded midrange values are used (WT/300 kg, AGE/1 yr).
    e_wt_vp2      <- fixed(-0.005); label("Effect of body weight on Vp2, power exponent (unitless)")        # Li 2014 Table 2 theta_1
    e_wt_q_liver  <- fixed( 0.008); label("Effect of body weight on Q_liver, power exponent (unitless)")    # Li 2014 Table 2 theta_2
    e_age_vc      <- fixed(-0.011); label("Effect of age on Vc, power exponent (unitless)")                 # Li 2014 Table 2 theta_3

    # IIV (Li 2014 Table 2 "IIV" column). Phoenix NLME default reports omega^2
    # (variance of the random effect on log scale). Parameters listed with
    # "NE (NA)" in the source had IIV not estimated and are omitted from the
    # ini() block (their model() expressions carry no eta term): Q2 (CL3),
    # Q_liver (CLL), Ka_sc (Ksc), Ka_oral (Kpo). The cattle IIVs for Vc
    # (9.59), V_liver (2.76), V_kidney (2.91), Q_kidney (4.04), and the
    # bioavailability terms (1.24-2.31) are unusually large; see the vignette
    # Errata for the interpretation choice (taken as printed).
    etalvc            ~ 9.59  # Li 2014 Table 2 IIV V1
    etalvp            ~ 0.89  # Li 2014 Table 2 IIV V2
    etalvp2           ~ 0.45  # Li 2014 Table 2 IIV V3
    etalv_liver       ~ 2.76  # Li 2014 Table 2 IIV VL
    etalv_kidney      ~ 2.91  # Li 2014 Table 2 IIV VK
    etalcl            ~ 0.89  # Li 2014 Table 2 IIV CL1
    etalq             ~ 0.37  # Li 2014 Table 2 IIV CL2
    etalq_kidney      ~ 4.04  # Li 2014 Table 2 IIV CLK
    etalka_im_na      ~ 0.05  # Li 2014 Table 2 IIV Kim1
    etalka_im_pro     ~ 0.35  # Li 2014 Table 2 IIV Kim2
    etalfdepot_im_na  ~ 1.24  # Li 2014 Table 2 IIV Fim1
    etalfdepot_im_pro ~ 1.84  # Li 2014 Table 2 IIV Fim2
    etalfdepot_sc     ~ 2.31  # Li 2014 Table 2 IIV Fsc
    etalfdepot_oral   ~ 2.01  # Li 2014 Table 2 IIV Fpo

    # Residual error: proportional error model with a single residual error
    # value 23% reported in Li 2014 Table 2; that single value is applied to
    # plasma and tissue outputs alike (nlmixr2 requires a distinct endpoint
    # parameter per output, so the same numeric is repeated three times below).
    propSd         <- 0.23;  label("Proportional residual error, plasma (fraction)")  # Li 2014 Table 2
    propSd_Cliver  <- 0.23;  label("Proportional residual error, liver (fraction)")   # Li 2014 Table 2 (same value as plasma)
    propSd_Ckidney <- 0.23;  label("Proportional residual error, kidney (fraction)")  # Li 2014 Table 2 (same value as plasma)
  })

  model({
    # Individual parameters with IIV (where estimated) and covariate effects.
    vc        <- exp(lvc + e_age_vc * log(AGE / 1) + etalvc)
    vp        <- exp(lvp + etalvp)
    vp2       <- exp(lvp2 + e_wt_vp2 * log(WT / 300) + etalvp2)
    v_liver   <- exp(lv_liver + etalv_liver)
    v_kidney  <- exp(lv_kidney + etalv_kidney)
    cl        <- exp(lcl + etalcl)
    q         <- exp(lq + etalq)
    q2        <- exp(lq2)
    q_liver   <- exp(lq_liver + e_wt_q_liver * log(WT / 300))
    q_kidney  <- exp(lq_kidney + etalq_kidney)
    ka_im_na  <- exp(lka_im_na + etalka_im_na)
    ka_im_pro <- exp(lka_im_pro + etalka_im_pro)
    ka_sc     <- exp(lka_sc)
    ka_oral   <- exp(lka_oral)
    f_im_na   <- exp(lfdepot_im_na + etalfdepot_im_na)
    f_im_pro  <- exp(lfdepot_im_pro + etalfdepot_im_pro)
    f_sc      <- exp(lfdepot_sc + etalfdepot_sc)
    f_oral    <- exp(lfdepot_oral + etalfdepot_oral)

    # ODE system (Li 2014 equation 1 and supporting equations).
    # The Methods state that the oral procaine penicillin depot feeds the
    # liver compartment directly. Each tissue compartment (liver, kidney) is
    # connected to the central compartment by a single inter-compartmental
    # clearance (Q_liver = CLL, Q_kidney = CLK); see vignette Errata for the
    # interpretation that ignores the printed "- CLL * C_L" tissue-elimination
    # term as a typesetting artefact (the table description "Clearance between
    # the central and the X compartment" denotes inter-compartmental clearance,
    # and the Methods state that CL of peripheral 1 ~ CL of tissue compartments).
    # depot1 = IM penicillin sodium; depot2 = IM procaine penicillin;
    # depot3 = SC procaine penicillin; depot4 = oral procaine penicillin
    # (depot4 feeds the liver compartment, not the central compartment).
    d/dt(depot1) <- -ka_im_na  * depot1
    d/dt(depot2) <- -ka_im_pro * depot2
    d/dt(depot3) <- -ka_sc     * depot3
    d/dt(depot4) <- -ka_oral   * depot4
    d/dt(central) <-
        ka_im_na  * depot1 +
        ka_im_pro * depot2 +
        ka_sc     * depot3 -
        cl       * central / vc -
        q        * (central / vc - peripheral1 / vp) -
        q2       * (central / vc - peripheral2 / vp2) -
        q_liver  * (central / vc - liver       / v_liver) -
        q_kidney * (central / vc - kidney      / v_kidney)
    d/dt(peripheral1) <- q  * (central / vc - peripheral1 / vp)
    d/dt(peripheral2) <- q2 * (central / vc - peripheral2 / vp2)
    d/dt(liver)       <- ka_oral * depot4 + q_liver  * (central / vc - liver  / v_liver)
    d/dt(kidney)      <- q_kidney * (central / vc - kidney / v_kidney)

    # Route-specific bioavailability on each absorption depot.
    f(depot1) <- f_im_na
    f(depot2) <- f_im_pro
    f(depot3) <- f_sc
    f(depot4) <- f_oral

    # Observation variables: plasma/serum, liver, and kidney concentrations
    # share the single proportional residual error reported in Li 2014 Table 2.
    Cc      <- central / vc
    Cliver  <- liver   / v_liver
    Ckidney <- kidney  / v_kidney
    Cc      ~ prop(propSd)
    Cliver  ~ prop(propSd_Cliver)
    Ckidney ~ prop(propSd_Ckidney)
  })
}
