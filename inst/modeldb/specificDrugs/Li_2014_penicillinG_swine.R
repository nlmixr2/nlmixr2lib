Li_2014_penicillinG_swine <- function() {
  description <- "Preclinical (swine). Three-compartment population pharmacokinetic model for penicillin G in swine, with two parallel first-order absorption depots covering intramuscular penicillin potassium and intramuscular procaine penicillin, plus separate kidney and muscle tissue compartments connected to the central compartment by inter-compartmental clearance; pooled meta-analysis of 89 pigs from 13 published studies and one unpublished FDA dataset (Li 2014)."
  reference   <- "Li M, Gehring R, Tell L, Baynes R, Huang Q, Riviere JE. Interspecies mixed-effect pharmacokinetic modeling of penicillin G in cattle and swine. Antimicrob Agents Chemother. 2014;58(8):4495-4503. doi:10.1128/AAC.02806-14"
  vignette    <- "Li_2014_penicillinG"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used as a power-model covariate on Vp (first peripheral volume) per Li 2014 Table 3 covariate factor theta_1 = 0.132. Form: Vp_i = Vp_pop * (WT/50)^0.132. Reference WT (50 kg) is a rounded midrange of the swine dataset (Li 2014 Table 1 weights 3.3-221 kg); the paper did not report the exact normalisation weight.",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "swine (Sus scrofa domesticus); piglets and adult pigs pooled from 13 published studies plus one FDA unpublished dataset",
    n_subjects     = 89L,
    n_studies      = 14L,
    age_range      = "0.02-0.259 years (piglet through grower pig, pooled across studies where reported)",
    weight_range   = "3.3-221.1 kg (piglet through adult, pooled across studies)",
    sex_female_pct = NA_real_,
    disease_state  = "Healthy animals only (animals with various diseased conditions were excluded; Li 2014 Methods)",
    dose_range     = "Penicillin potassium 7.64-52.6 mg/kg IV or IM; procaine penicillin 8.47-33 mg/kg IM (Li 2014 Table 1)",
    regions        = "Published swine PK studies (US and international); FARAD records; one unpublished FDA dataset (Shelver et al.)",
    n_observations = 443L,
    notes          = "Pooled plasma data (443 concentrations from 89 pigs) plus 97 kidney and 84 muscle tissue residue concentrations. Swine liver data were too sparse to model. Routes carried in the dosing dataset via the rxode2 cmt column: 'central' for IV; 'depot1' for intramuscular penicillin potassium; 'depot2' for intramuscular procaine penicillin. The companion cattle model from this paper is modellib('Li_2014_penicillinG_cattle')."
  )

  ini({
    # Structural parameters from Li 2014 Table 3 (swine, final PK model).
    lvc       <- log(3.05);  label("Central volume of distribution (Vc, L)")                                   # Li 2014 Table 3 V1
    lvp       <- log(1.65);  label("First peripheral volume of distribution (Vp, L)")                          # Li 2014 Table 3 V2
    lvp2      <- log(4.65);  label("Second peripheral volume of distribution (Vp2, L)")                        # Li 2014 Table 3 V3
    lv_kidney <- log(4.38);  label("Kidney compartment volume of distribution (V_kidney, L)")                  # Li 2014 Table 3 VK
    lv_muscle <- log(1.10);  label("Muscle compartment volume of distribution (V_muscle, L)")                  # Li 2014 Table 3 VM
    lcl       <- log(16.9);  label("Central clearance (CL, L/h)")                                              # Li 2014 Table 3 CL1
    lq        <- log(13.7);  label("Inter-compartmental clearance central <-> peripheral1 (Q, L/h)")           # Li 2014 Table 3 CL2
    lq2       <- log(0.52);  label("Inter-compartmental clearance central <-> peripheral2 (Q2, L/h)")          # Li 2014 Table 3 CL3
    lq_kidney <- log(12.1);  label("Inter-compartmental clearance central <-> kidney (Q_kidney, L/h)")         # Li 2014 Table 3 CLK
    lq_muscle <- log(14.8);  label("Inter-compartmental clearance central <-> muscle (Q_muscle, L/h)")         # Li 2014 Table 3 CLM
    lka_im_k   <- log(3.03); label("First-order absorption rate constant after IM penicillin potassium (Ka, 1/h)")  # Li 2014 Table 3 Kim1
    lka_im_pro <- log(0.48); label("First-order absorption rate constant after IM procaine penicillin (Ka, 1/h)")   # Li 2014 Table 3 Kim2
    lfdepot_im_k   <- log(0.733); label("Bioavailability of IM penicillin potassium (F, fraction)")    # Li 2014 Table 3 Fim1 (73.3%)
    lfdepot_im_pro <- log(0.642); label("Bioavailability of IM procaine penicillin (F, fraction)")     # Li 2014 Table 3 Fim2 (64.2%)

    # Covariate effect (Li 2014 Table 3). Power model:
    # Vp_i = Vp_pop * (WT/50)^0.132 * exp(eta).
    e_wt_vp  <- fixed(0.132); label("Effect of body weight on Vp, power exponent (unitless)")  # Li 2014 Table 3 theta_1

    # IIV (Li 2014 Table 3 "IIV" column). Phoenix NLME default reports omega^2
    # (variance of the random effect on log scale). Parameters listed with
    # "NE (NA)" had IIV not estimated and are omitted (their model() expressions
    # carry no eta term): Vp2 (V3), Q (CL2), Q2 (CL3).
    etalvc            ~ 0.13  # Li 2014 Table 3 IIV V1
    etalvp            ~ 0.03  # Li 2014 Table 3 IIV V2
    etalv_kidney      ~ 0.02  # Li 2014 Table 3 IIV VK
    etalv_muscle      ~ 0.03  # Li 2014 Table 3 IIV VM
    etalcl            ~ 0.12  # Li 2014 Table 3 IIV CL1
    etalq_kidney      ~ 0.06  # Li 2014 Table 3 IIV CLK
    etalq_muscle      ~ 0.06  # Li 2014 Table 3 IIV CLM
    etalka_im_k       ~ 0.04  # Li 2014 Table 3 IIV Kim1
    etalka_im_pro     ~ 0.13  # Li 2014 Table 3 IIV Kim2
    etalfdepot_im_k   ~ 0.02  # Li 2014 Table 3 IIV Fim1
    etalfdepot_im_pro ~ 0.14  # Li 2014 Table 3 IIV Fim2

    # Residual error: proportional error model with a single residual error
    # value 35% reported in Li 2014 Table 3; applied to plasma and tissue
    # outputs alike (nlmixr2 requires a distinct endpoint parameter per
    # output, so the same numeric is repeated three times below).
    propSd         <- 0.35;  label("Proportional residual error, plasma (fraction)")  # Li 2014 Table 3
    propSd_Ckidney <- 0.35;  label("Proportional residual error, kidney (fraction)")  # Li 2014 Table 3 (same value as plasma)
    propSd_Cmuscle <- 0.35;  label("Proportional residual error, muscle (fraction)")  # Li 2014 Table 3 (same value as plasma)
  })

  model({
    # Individual parameters with IIV (where estimated) and covariate effect.
    vc        <- exp(lvc + etalvc)
    vp        <- exp(lvp + e_wt_vp * log(WT / 50) + etalvp)
    vp2       <- exp(lvp2)
    v_kidney  <- exp(lv_kidney + etalv_kidney)
    v_muscle  <- exp(lv_muscle + etalv_muscle)
    cl        <- exp(lcl + etalcl)
    q         <- exp(lq)
    q2        <- exp(lq2)
    q_kidney  <- exp(lq_kidney + etalq_kidney)
    q_muscle  <- exp(lq_muscle + etalq_muscle)
    ka_im_k   <- exp(lka_im_k + etalka_im_k)
    ka_im_pro <- exp(lka_im_pro + etalka_im_pro)
    f_im_k    <- exp(lfdepot_im_k + etalfdepot_im_k)
    f_im_pro  <- exp(lfdepot_im_pro + etalfdepot_im_pro)

    # ODE system (Li 2014 equation 1, swine variant). Each tissue compartment
    # (kidney, muscle) is connected to the central compartment by a single
    # inter-compartmental clearance (Q_kidney = CLK, Q_muscle = CLM); see vignette
    # Errata for the interpretation that ignores the printed "- CLK * C_K" term
    # in the kidney equation as a typesetting artefact (the swine muscle
    # equation as printed already has only the inter-compartmental term).
    # depot1 = IM penicillin potassium; depot2 = IM procaine penicillin.
    d/dt(depot1) <- -ka_im_k   * depot1
    d/dt(depot2) <- -ka_im_pro * depot2
    d/dt(central) <-
        ka_im_k   * depot1 +
        ka_im_pro * depot2 -
        cl       * central / vc -
        q        * (central / vc - peripheral1 / vp) -
        q2       * (central / vc - peripheral2 / vp2) -
        q_kidney * (central / vc - kidney      / v_kidney) -
        q_muscle * (central / vc - muscle      / v_muscle)
    d/dt(peripheral1) <- q  * (central / vc - peripheral1 / vp)
    d/dt(peripheral2) <- q2 * (central / vc - peripheral2 / vp2)
    d/dt(kidney)      <- q_kidney * (central / vc - kidney / v_kidney)
    d/dt(muscle)      <- q_muscle * (central / vc - muscle / v_muscle)

    # Route-specific bioavailability on each absorption depot.
    f(depot1) <- f_im_k
    f(depot2) <- f_im_pro

    # Observation variables: plasma, kidney, and muscle concentrations share
    # the single proportional residual error reported in Li 2014 Table 3.
    Cc      <- central / vc
    Ckidney <- kidney  / v_kidney
    Cmuscle <- muscle  / v_muscle
    Cc      ~ prop(propSd)
    Ckidney ~ prop(propSd_Ckidney)
    Cmuscle ~ prop(propSd_Cmuscle)
  })
}
