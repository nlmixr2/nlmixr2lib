Frymoyer_2013_mycophenolic_acid <- function() {
  description <- "Population PK model for unbound mycophenolic acid (MPA, the active moiety of mycophenolate mofetil MMF) in adult allogeneic haematopoietic cell transplantation (alloHCT) recipients (Frymoyer 2013). Two-compartment disposition with first-order absorption and linear elimination; oral bioavailability F = 0.560. Absorption lag time follows a two-class mixture: Group 1 (no delay, ALAG = 0, 91 % of subjects) and Group 2 (delayed absorption, ALAG = 1.96 h, 9 % of subjects), gated by the latent MIX_LAGGED_ABS class indicator. Creatinine clearance (Cockcroft-Gault with ideal body weight, NOT BSA-normalized) is the only retained covariate, entering CL as a power scaling (CRCL / 86 mL/min)^0.207. Inter-individual variability is log-normal on all structural parameters and on F; residual error is proportional. Doses are MPA-equivalent (mg) -- MMF mass must be converted externally via F_MW = 0.739 (oral) or 0.682 (IV) per Frymoyer 2013 Methods."
  reference <- paste(
    "Frymoyer A, Verotta D, Jacobson P, Long-Boyle J.",
    "Population pharmacokinetics of unbound mycophenolic acid in adult",
    "allogeneic haematopoietic cell transplantation: effect of",
    "pharmacogenetic factors.",
    "Br J Clin Pharmacol. 2013;75(2):463-475.",
    "doi:10.1111/j.1365-2125.2012.04372.x.",
    sep = " "
  )
  vignette <- "Frymoyer_2013_mycophenolic_acid"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance estimated by the Cockcroft-Gault equation using ideal body weight (NOT BSA-normalized).",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Cockcroft-Gault raw mL/min, per Frymoyer 2013 Methods 'Patients' subsection: 'Creatinine clearance (CLcr) was estimated by the Cockcroft-Gault equation using ideal body weight'. Cohort median 86 mL/min (range 30-197 mL/min; Table 1). Enters CL as a power scaling CL = 1610 * (CRCL/86)^0.207 (Table 3; Results 'Final population PK model' equation).",
      source_name        = "CLcr"
    ),
    MIX_LAGGED_ABS = list(
      description        = paste(
        "Per-subject latent mixture-model class indicator from the",
        "Frymoyer 2013 absorption model. 1 = subject in the",
        "delayed-absorption Group 2 subpopulation (lag time fixed at 1.96",
        "h); 0 = subject in the no-delay Group 1 subpopulation (lag time",
        "fixed at 0). Not a measured clinical covariate -- the mixture",
        "assignment is the per-subject latent-class index from a NONMEM",
        "$MIXTURE block (Frymoyer 2013 Methods / Population PK model",
        "building paragraph 2: 'a mixture model was implemented that",
        "allowed for two sub-populations with different absorption lag",
        "time')."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Group 1 = no-delay majority class; 91 % of the Frymoyer 2013 source cohort)",
      notes              = paste(
        "The population probability of Group 2 (delayed absorption) is the",
        "estimated mixture fraction P2 = 0.09 (Table 3; identified 8 of",
        "132 subjects in the source cohort). Group 1 (no-delay) is the",
        "complement with P1 = 0.91 (Table 3, RSE 4.0 %). For",
        "typical-value simulation set MIX_LAGGED_ABS = 0 to reproduce the",
        "dominant 91 % no-delay phenotype (ka = 0.400 /h, lag = 0); set",
        "MIX_LAGGED_ABS = 1 to reproduce the rarer 9 % delayed-absorption",
        "phenotype (ka = 0.400 /h, lag = 1.96 h). For population",
        "simulation, draw MIX_LAGGED_ABS ~ Bernoulli(0.09) per subject.",
        "Frymoyer 2013 also tested allowing ka and F to vary between",
        "subpopulations (an extended mixture model) but reported no",
        "improvement; only ALAG differs by class."
      ),
      source_name        = "MIXTURE (NONMEM $MIXTURE assignment; the binary indicator is MIX_LAGGED_ABS = as.integer(MIXTURE == 2))"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 132L,
    n_studies       = 3L,
    n_observations  = "1171 quantifiable unbound MPA concentrations across 179 intensive PK profiles (80 intravenous, 99 oral); 47 subjects had paired oral and IV profiles; 8 plasma unbound MPA concentrations were below the level of quantification (1 ng/mL) and included in the analysis.",
    age_range       = "19-69 years",
    age_median      = "52 years (Table 1)",
    weight_range    = "50-149 kg (actual body weight)",
    weight_median   = "80 kg (Table 1; actual)",
    bsa_range       = "1.5-2.8 m^2",
    bsa_median      = "2.0 m^2 (Table 1)",
    crcl_range      = "30-197 mL/min",
    crcl_median     = "86 mL/min (Cockcroft-Gault using ideal body weight; Table 1)",
    albumin_range   = "2.0-4.3 g/dL",
    albumin_median  = "2.9 g/dL (Table 1)",
    bilirubin_range = "0.1-20 mg/dL (10 subjects with total bilirubin > 3 mg/dL)",
    sex_female_pct  = 38,
    race_ethnicity  = "Not reported (single-centre US adult HCT cohort).",
    disease_state   = "Adult recipients of nonmyeloablative fludarabine-based allogeneic haematopoietic cell transplantation receiving mycophenolate mofetil (MMF) + ciclosporin for acute GVHD immunosuppression. Underlying diseases include non-Hodgkin's lymphoma (23 %), Hodgkin's lymphoma (11 %), chronic myelogenous leukaemia (14 %), acute leukaemia (33 %), myelodysplastic syndrome (14 %) and other (5 %). Stem-cell source: umbilical cord blood (62 %), peripheral blood stem cells (32 %), bone marrow (6 %). Renal impairment per FDA definitions: 33 % mild (CrCl > 50-80 mL/min), 5 % moderate (30-50 mL/min), 0 % severe.",
    dose_range      = "MMF (Cellcept) 1 g every 12 h (n = 98), 1 g every 8 h (n = 19), or 1.5 g every 12 h (n = 15), administered oral or intravenously. MMF doses are converted to MPA-equivalent mass for the model via F_MW = 0.739 (oral MMF, MW 433.5; MPA MW 320.3) and F_MW = 0.682 (IV MMF, MW 469.96 for the hydrochloride salt; MPA MW 320.3) -- Frymoyer 2013 Methods / Population PK analysis paragraph 1.",
    regions         = "USA (University of Minnesota / University of California San Francisco retrospective meta-analysis pooling three previously published intensive-PK studies).",
    co_medication   = "100 % ciclosporin (concomitant by inclusion criterion); 95 % fluoroquinolones; 72 % fluconazole; 58 % proton pump inhibitors; 27 % glucocorticoids; 19 % ursodiol; 12 % hormonal contraceptives; 2 % seizure prophylaxis. 38 % were on at least one potential UGT inducer (corticosteroid / hormonal contraceptive / anticonvulsant) at the time of PK sampling.",
    notes           = "Retrospective pooling of three previously published intensive-PK studies (Frymoyer 2013 references [4, 30, 31]). PK sampling performed one or two times within the first 15 days post-transplant at steady state. Pharmacogenetic variants tested but not retained in the final PK model: UGT1A8*2 (rs1042597), UGT1A8*3 (rs17863762), UGT1A9 98T>C, UGT1A9 -2152C>T (rs17868320), UGT1A9 -275T>A (rs6714486), UGT1A10*2 (rs10187694), UGT2B7 802C>T (rs7439366), MRP2 -24C>T (rs717620), MRP2 3972C>T (rs3740066), MRP2 1249G>A (rs2273697); for UGT1A9*2 and UGT1A10*2 no variant alleles were observed. Continuous covariates screened (forward addition) but not retained after backward elimination: body weight, body surface area, albumin, CSA trough concentration, day of PK sampling relative to stem cell infusion. Categorical covariates screened but not retained: gender, total bilirubin, concomitant medications (CSA, fluoroquinolones, fluconazole, proton pump inhibitors, glucocorticoids, ursodiol, hormonal contraceptives, seizure prophylaxis), genotype. Only CRCL survived backward elimination (P < 0.01)."
  )

  ini({
    # ------------------------------------------------------------------
    # Final population PK model parameter estimates (Frymoyer 2013 Table 3,
    # "Final model" column). Median typical values; %CV are reported on
    # the linear scale and converted to log-normal omega^2 via
    # omega^2 = log(1 + CV^2).
    # ------------------------------------------------------------------

    # Structural disposition -- two-compartment first-order absorption with
    # linear elimination. CL is parameterized at the reference CRCL of 86
    # mL/min (the cohort median; Table 1).
    lka     <- log(0.400); label("Absorption rate constant Ka (1/h)")                                  # Table 3 ka = 0.400 1/h (RSE 9.3 %)
    lcl     <- log(1610);  label("Clearance CL at reference CRCL = 86 mL/min (L/h)")                   # Table 3 CL = 1610 L/h (RSE 5.8 %)
    lvc     <- log(1230);  label("Central volume of distribution Vc (L)")                              # Table 3 V_central = 1230 L (RSE 12.5 %)
    lq      <- log(541);   label("Intercompartmental clearance Q (L/h)")                               # Table 3 Q = 541 L/h (RSE 13.1 %)
    lvp     <- log(6140);  label("Peripheral volume of distribution Vp (L)")                           # Table 3 V_peripheral = 6140 L (RSE 21.0 %)
    lfdepot <- log(0.560); label("Oral bioavailability F on depot dose (unitless)")                    # Table 3 F = 0.560 (RSE 7.4 %)

    # CRCL covariate on CL (power form). Reference 86 mL/min = cohort
    # median (Table 1).
    e_crcl_cl <- 0.207; label("Power-form exponent of CRCL on CL (unitless)")                          # Table 3 d (CRCL effect on CL) = 0.207 (RSE 39.3 %); Results "Final population PK model" CL = 1610 * (CLcr/86)^0.207

    # Mixture-model lag time for the delayed-absorption Group 2 subpopulation.
    # Group 1 lag is fixed at 0 (Table 3 "ALAG_P1, lag time for population
    # without delayed absorption (h) -- 0, fixed"); Group 2 lag is
    # estimated at 1.96 h. The mixture proportions (P1 = 0.91, P2 = 0.09)
    # are encoded structurally via the MIX_LAGGED_ABS covariate (see
    # covariateData[[MIX_LAGGED_ABS]]$notes) and are not estimated here --
    # the typical-value library does not estimate mixture proportions.
    ltlag <- log(1.96); label("Absorption lag time for the delayed-absorption Group 2 subpopulation (h)") # Table 3 ALAG_P2 = 1.96 h (RSE 1.4 %)

    # Inter-individual variability (log-normal; omega^2 = log(1 + CV^2)).
    # The paper "Allowing for inter-patient variability in lag time did
    # not improve the model and therefore was not included in the final
    # model" (Results paragraph 4 of Population PK model building) -- so
    # no etaltlag.
    etalka     ~ 0.4272 # Table 3 IIV ka             = 73.0 % CV -> log(1 + 0.730^2)
    etalcl     ~ 0.1311 # Table 3 IIV CL             = 37.4 % CV -> log(1 + 0.374^2)
    etalvc     ~ 0.1318 # Table 3 IIV V_central      = 37.5 % CV -> log(1 + 0.375^2)
    etalq      ~ 0.4524 # Table 3 IIV Q              = 75.6 % CV -> log(1 + 0.756^2)
    etalvp     ~ 0.8920 # Table 3 IIV V_peripheral   = 120  % CV -> log(1 + 1.20^2)
    etalfdepot ~ 0.1820 # Table 3 IIV F              = 44.7 % CV -> log(1 + 0.447^2)

    # Residual error -- proportional only (Methods: "Residual variability
    # was assumed to follow a proportional error model").
    propSd <- 0.423; label("Proportional residual error (fraction)")                                   # Table 3 Residual intra-individual variability = 42.3 % CV (RSE 6.9 %)
  })

  model({
    # Individual parameters (log-normal IIV).
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (CRCL / 86)^e_crcl_cl
    vc <- exp(lvc + etalvc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    # Mixture-class-gated absorption lag. Group 1 (MIX_LAGGED_ABS = 0):
    # lag = 0 (structurally fixed in Table 3). Group 2 (MIX_LAGGED_ABS =
    # 1): lag = exp(ltlag) = 1.96 h.
    tlag <- exp(ltlag) * MIX_LAGGED_ABS

    # Micro-constants for two-compartment disposition.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system: depot -> central <-> peripheral1.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Oral bioavailability is applied to depot doses only; IV doses
    # delivered directly to the central compartment use the default F = 1.
    # Doses must be MPA-equivalent (mg) -- the user converts MMF mass
    # externally via F_MW = 0.739 (oral) or F_MW = 0.682 (IV) per
    # Frymoyer 2013 Methods (see population$dose_range).
    f(depot)    <- exp(lfdepot + etalfdepot)
    alag(depot) <- tlag

    # Observation: dose in mg, vc in L -> central/vc in mg/L = ug/mL.
    # Multiply by 1000 to express Cc in ng/mL (the unit used for unbound
    # MPA in Frymoyer 2013).
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
