Ogungbenro_2009_penciclovir <- function() {
  description <- "Two-compartment population PK model with first-order absorption and lag time for penciclovir in pooled adults and children (Ogungbenro 2009). Famciclovir is the oral prodrug of penciclovir; both oral famciclovir and intravenous penciclovir doses are described jointly (six clinical studies, 69 subjects of whom 23 are children, 160 occasions, 1676 plasma penciclovir observations). Allometric body-weight scaling with reference 70 kg (exponent 0.75 shared on CL and Q, exponent 1.0 shared on V1 and V2), an empirical piecewise age effect on CL with separate K parameters for AGE < 40 years (rising-with-youth limb) and AGE >= 40 years (declining-with-elderly limb), and a power function of creatinine clearance on CL with reference 100 mL/min (Cockcroft-Gault, raw mL/min). Inter-individual variability on ka, CL, V1 (fixed at omega^2 = 0.003), V2, and Q; combined proportional plus additive residual error (additive variance fixed at 0.01 mg^2/L^2)."
  reference <- paste(
    "Ogungbenro K, Matthews I, Looby M, Kaiser G, Graham G, Aarons L. (2009).",
    "Population pharmacokinetics and optimal design of paediatric studies for famciclovir.",
    "British Journal of Clinical Pharmacology 68(4):546-560.",
    "doi:10.1111/j.1365-2125.2009.03479.x",
    sep = " "
  )
  vignette <- "Ogungbenro_2009_penciclovir"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling with reference 70 kg: exponent 0.75 (shared) on CL and Q, exponent 1.0 (shared) on V1 and V2 (Ogungbenro 2009 Methods, allometric-model paragraph: 'F_WT_CL = (WT/WTSTD)^(3/4)' and 'F_WT_V = (WT/WTSTD)^1' with WTSTD = 70 kg). Cohort range 13.9-94.6 kg (Table 2 Combined column: mean 59.3, SD 23.7; children mean 29.5, SD 12.2; adults mean 74.1, SD 9.7).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Empirical piecewise age effect on CL with reference age 40 years (Ogungbenro 2009 Methods, age-model paragraph): F_AGE_CL = (K_AGE - AGE)/(K_AGE - 40), with K_AGE = 159 for AGE < 40 years and K_AGE = 113 for AGE >= 40 years. Both limbs equal 1 at AGE = 40 so the covariate is continuous at the boundary. Cohort range 2-63 years (Table 2 Combined column: mean 26.5, SD 15.8; children 8.1 +/- 3.4 years, adults 35.8 +/- 10.6 years).",
      source_name        = "AGE"
    ),
    CRCL = list(
      description        = "Creatinine clearance estimated by the Cockcroft-Gault equation (raw mL/min, NOT BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Ogungbenro 2009 Methods computes individual creatinine clearance using the Cockcroft-Gault equation; values are raw mL/min and NOT BSA-normalized to mL/min/1.73 m^2. The covariate-columns register permits raw Cockcroft-Gault under the canonical CRCL name when the per-model description records the assay form (Delattre 2010 amikacin precedent). The covariate effect is a power function `(CRCL/100)^e_crcl_cl` with reference 100 mL/min (Methods, creatinine-clearance paragraph). Cohort range 27.6-175.6 mL/min (Table 2 Combined column: mean 87.9, SD 34.5; children 58.2 +/- 19.9, adults 102.8 +/- 30.4).",
      source_name        = "CLCR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 69L,
    n_studies      = 6L,
    age_range      = "2-63 years",
    age_median     = "mean 26.5 years (SD 15.8); children 8.1 +/- 3.4 years, adults 35.8 +/- 10.6 years",
    weight_range   = "13.9-94.6 kg",
    weight_median  = "mean 59.3 kg (SD 23.7); children 29.5 +/- 12.2 kg, adults 74.1 +/- 9.7 kg",
    sex_female_pct = 10.1,
    race_ethnicity = "Not reported in source paper",
    disease_state  = "Pooled adult and paediatric cohorts across six clinical trials supplied by Novartis: (1) single ascending oral dose in healthy adults, (2) IV-oral bioavailability crossover in healthy adults, (3) single oral dose in adults with renal impairment, (4) single IV infusion in immunocompromised paediatric patients (2-17 years) with two subjects re-dosed orally after washout, (5) multiple oral doses in paediatric chronic hepatitis B patients (6-11 years), and (6) multiple IV infusion in adults.",
    dose_range     = "Oral famciclovir (prodrug; single or multiple doses) and intravenous penciclovir (single or multiple-infusion doses). Adult reference dose was 500 mg famciclovir orally; paediatric simulations targeted a 10 mg/kg oral dose. F = 0.598 represents the bioavailability of systemic penciclovir delivered via oral famciclovir (combining prodrug conversion plus absorption efficiency).",
    regions        = "Not specified (data supplied by Novartis AG, Basel, Switzerland)",
    notes          = "Baseline demographics from Ogungbenro 2009 Table 2; study details from Table 1. 23 children (39 occasions, 322 concentrations) and 46 adults (121 occasions, 1354 concentrations) contributed 1676 plasma penciclovir observations in total. Sex distribution 62 M / 7 F across the combined cohort (Table 2). The covariates Age, Weight, Sex, Serum creatinine, and Creatinine clearance were available; Sex was not retained in the final model. Final model selection used a 3.84-point reduction in NONMEM VI objective function value per added parameter, supplemented by visual diagnostics. Bootstrap analysis (1000 replicates, 923 successful) confirmed parameter estimates and standard errors."
  )

  ini({
    # Structural parameters -- Ogungbenro 2009 Table 4 'Original data' column.
    # Reference subject is a 70 kg adult at AGE = 40 years and CRCL = 100 mL/min
    # (Cockcroft-Gault, raw mL/min). Seven typical-value point estimates from
    # Table 4: ka = 1.86 1/h, CL = 31.2 L/h, V1 = 28.6 L, V2 = 54.5 L,
    # Q = 60.2 L/h, F = 0.598, Tlag = 0.206 h.
    lka     <- log(1.86);  label("Absorption rate constant of penciclovir released from oral famciclovir prodrug (1/h)")  # Table 4: ka = 1.86 (CV% 10.3)
    lcl     <- log(31.2);  label("Clearance at 70 kg, AGE = 40 y, CRCL = 100 mL/min (L/h)")                                # Table 4: CL = 31.2 (CV% 6.0)
    lvc     <- log(28.6);  label("Central volume V1 at 70 kg (L)")                                                          # Table 4: V1 = 28.6 (CV% 6.0)
    lvp     <- log(54.5);  label("Peripheral volume V2 at 70 kg (L)")                                                       # Table 4: V2 = 54.5 (CV% 4.9)
    lq      <- log(60.2);  label("Intercompartmental clearance Q at 70 kg (L/h)")                                           # Table 4: Q  = 60.2 (CV% 7.1)
    lfdepot <- log(0.598); label("Relative bioavailability of penciclovir delivered as oral famciclovir (unitless)")        # Table 4: F  = 0.598 (CV% 2.9)
    ltlag   <- log(0.206); label("Absorption lag time after oral famciclovir (h)")                                          # Table 4: Tlag = 0.206 (CV% 2.2)

    # Allometric body-weight exponents -- fixed a priori per Anderson-Holford
    # size scaling. Ogungbenro 2009 Methods writes the allometric model as
    # F_WT_CL = (WT/WTSTD)^(3/4) (applied to CL AND Q -- both clearance terms)
    # and F_WT_V = (WT/WTSTD)^1 (applied to V1 AND V2 -- both volume terms).
    # The shared-exponent form e_<cov>_<param1>_<param2> matches the
    # three-token convention in references/parameter-names.md.
    e_wt_cl_q  <- fixed(0.75); label("Shared allometric WT exponent on CL and Q (unitless)")  # Methods, allometric paragraph (clearance terms, 0.75)
    e_wt_vc_vp <- fixed(1.00); label("Shared allometric WT exponent on V1 and V2 (unitless)") # Methods, allometric paragraph (volume terms, 1.0)

    # Piecewise empirical age effect on CL -- Ogungbenro 2009 Methods,
    # age-model paragraph; Table 4. The empirical fractional effect is
    #   F_AGE_CL = (K_AGE - AGE) / (K_AGE - 40)
    # with K_AGE = 159 for AGE < 40 (rising-with-youth limb; F_AGE_CL > 1 for
    # children, peaking at F = 1.32 in a 2-year-old) and K_AGE = 113 for
    # AGE >= 40 (declining-with-elderly limb; F_AGE_CL = 0.73 at AGE = 60).
    # Both limbs equal 1 at the reference AGE = 40, so the covariate is
    # continuous at the boundary.
    e_age_young_cl <- 159;  label("K_AGE for AGE < 40 in F_AGE = (K - AGE)/(K - 40) (years)")    # Table 4: K_AGE<40 = 159 (CV% 37.4)
    e_age_old_cl   <- 113;  label("K_AGE for AGE >= 40 in F_AGE = (K - AGE)/(K - 40) (years)")    # Table 4: K_AGE>=40 = 113 (CV% 24.4)

    # Power-function effect of creatinine clearance on CL -- Ogungbenro 2009
    # Methods, creatinine-clearance paragraph; Table 4 'Exponent of FCL_CR'.
    # Reference CRCL = 100 mL/min (raw Cockcroft-Gault).
    e_crcl_cl <- 0.28;  label("Power exponent of CRCL on CL via (CRCL/100)^e_crcl_cl (unitless)")  # Table 4: F_CL_CR exponent = 0.28 (CV% 45.7)

    # Inter-individual variability. Ogungbenro 2009 fits BSV terms in
    # NONMEM VI with the standard log-normal individual-PK form
    #   theta_i = theta_TV * exp(eta_i),  eta_i ~ N(0, omega^2),
    # and the Table 4 'BSV' column reports the omega^2 variance estimates
    # directly (NONMEM $OMEGA convention). BSV on V1 was held fixed at the
    # small value 0.003 because attempts to estimate a full ka/CL/V1/V2/Q
    # variance-covariance block led to stability issues (Results paragraph:
    # 'the final model only has diagonal elements'); BSV on F and Tlag was
    # excluded altogether ('BSV_F was low ... removed from the final model';
    # 'BSV_Tlag could not be estimated').
    etalka  ~ 0.640         # Table 4: BSV(ka) = 0.640 (CV% 25.9)
    etalcl  ~ 0.230         # Table 4: BSV(CL) = 0.23  (CV% 22.3)
    etalvc  ~ fixed(0.003)  # Table 4: BSV(V1) = 0.003 fix
    etalvp  ~ 0.255         # Table 4: BSV(V2) = 0.255 (CV% 29.3)
    etalq   ~ 0.342         # Table 4: BSV(Q)  = 0.342 (CV% 59.4)

    # Residual error -- Ogungbenro 2009 Methods reports a combined additive
    # plus proportional model fit in NONMEM VI. The Table 4 'Proportional
    # error' and 'Additive error' entries are NONMEM $SIGMA estimates, i.e.
    # variances on each epsilon. nlmixr2's propSd / addSd take standard
    # deviations, so we convert via sqrt(.). The additive variance was held
    # fixed at 0.01 mg^2/L^2 by the authors (Table 4 footnote 'fix').
    propSd <- sqrt(0.221);        label("Proportional residual SD (fraction)") # Table 4: proportional error variance = 0.221 (CV% 9.6); SD = sqrt(0.221) ~ 0.470
    addSd  <- fixed(sqrt(0.01));  label("Additive residual SD (mg/L)")         # Table 4: additive error variance = 0.01 fix; SD = sqrt(0.01) = 0.1 mg/L
  })

  model({
    # Piecewise empirical age effect on CL (Ogungbenro 2009 Methods,
    # age-model paragraph; Table 5 final-model equations). The K parameters
    # in ini() are e_age_young_cl (= 159 in the rising-with-youth limb) and
    # e_age_old_cl (= 113 in the declining-with-elderly limb).
    f_age_cl <- ifelse(
      AGE < 40,
      (e_age_young_cl - AGE) / (e_age_young_cl - 40),
      (e_age_old_cl   - AGE) / (e_age_old_cl   - 40)
    )

    # Power-function creatinine-clearance effect on CL (Methods,
    # creatinine-clearance paragraph). Reference 100 mL/min Cockcroft-Gault.
    f_crcl_cl <- (CRCL / 100)^e_crcl_cl

    # Individual PK parameters. Allometric body-weight scaling with
    # reference 70 kg; shared 0.75 exponent on the two clearance terms
    # (CL and Q) and shared 1.0 exponent on the two volume terms (V1 and
    # V2). Age and CRCL covariates multiply CL only (Table 5 row for CL).
    ka     <- exp(lka  + etalka)
    cl     <- exp(lcl  + etalcl) * (WT / 70)^e_wt_cl_q  * f_age_cl * f_crcl_cl
    vc     <- exp(lvc  + etalvc) * (WT / 70)^e_wt_vc_vp
    vp     <- exp(lvp  + etalvp) * (WT / 70)^e_wt_vc_vp
    q      <- exp(lq   + etalq)  * (WT / 70)^e_wt_cl_q
    fdepot <- exp(lfdepot)
    tlag   <- exp(ltlag)

    # Two-compartment disposition with first-order absorption from depot.
    # Oral famciclovir doses enter `depot` (cmt = depot) and are absorbed
    # into central with rate ka, bioavailability fdepot, and lag tlag.
    # Intravenous penciclovir doses enter `central` directly (cmt = central);
    # for those events fdepot and tlag have no effect because they apply
    # only to the depot compartment.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                               k12 * central - k21 * peripheral1

    f(depot)    <- fdepot
    alag(depot) <- tlag

    # Plasma penciclovir concentration; dose mg / V1 L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
