Han_2010_voriconazole <- function() {
  description <- "Two-compartment population pharmacokinetic model with first-order absorption and first-order elimination for intravenous and oral voriconazole in adult lung transplant recipients during the early postoperative period (Han 2010). Bioavailability is estimated for the oral route. The base structural model is reported as the primary result; three separate single-covariate sub-models -- cystic fibrosis (CF) and postoperative time (POT) on bioavailability, and body weight (WT) on peripheral volume -- are reported in the paper but were not combined into a final model; the base-model typical-value parameter estimates are encoded here, and the three covariate sub-models are reproduced in the validation vignette."
  reference <- "Han K, Capitano B, Bies R, Potoski BA, Husain S, Gilbert S, Paterson DL, McCurry K, Venkataramanan R. Bioavailability and Population Pharmacokinetics of Voriconazole in Lung Transplant Recipients. Antimicrob Agents Chemother. 2010. doi:10.1128/AAC.00504-10"
  vignette <- "Han_2010_voriconazole"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 13L,
    n_studies      = 1L,
    age_range      = "19-70 years",
    age_mean       = "50.9 +/- 16.1 years",
    weight_range   = "46-91 kg",
    weight_mean    = "68.0 +/- 15.2 kg",
    ideal_bodyweight_range = "45.5-75.3 kg",
    ideal_bodyweight_mean  = "59.6 +/- 8.2 kg",
    sex_female_pct = 46.2,
    race_ethnicity = c(Caucasian = 92.3, Other = 7.7),
    primary_diagnosis = c(
      CysticFibrosis_n          = 3,
      Emphysema_n               = 5,
      IdiopathicPulmonaryFibrosis_n = 4,
      Scleroderma_n             = 1
    ),
    days_posttransplant_oral_study = "8.5 +/- 4.4 days (range 3-19)",
    disease_state  = "Adult lung transplant recipients during the early postoperative period; all received tacrolimus as primary immunosuppression. 3 patients (23.1%) had cystic fibrosis; the remaining 10 had emphysema, idiopathic pulmonary fibrosis, or scleroderma. One patient did not complete the oral study.",
    dose_range     = "Loading two intravenous voriconazole doses 6 mg/kg given as 2-hour infusion every 12 h immediately post-transplant, followed by oral doses 200 mg every 12 h for 3 months posttransplant. Blood sampling at pre-dose and at 0.5, 1, 1.5, 2, 4, 6, 8, and 12 h following the second intravenous dose and following an oral dose (the 5th to 37th dose, mean 15th).",
    regions        = "Single center: University of Pittsburgh Medical Center, Pittsburgh, PA, USA.",
    notes          = "Prospective single-center observational study. Voriconazole plasma concentrations measured by validated HPLC (LLOQ implicit in linearity range 0.2-9 ug/mL, R^2 = 0.9998; assay precision 1.3-9.0%; assay bias 0.7-3.1%). NONMEM 6.2.0 (GloboMax) with FOCE-I. 9 samples per interval over one intravenous and one oral interval per patient; steady state confirmed by paired-t test on trough concentrations C0 vs C12 (P = 0.82, mean difference -2.7%). Baseline demographics per Han 2010 Table 1 (the 'Characteristics of patients' table); base-model population PK parameter estimates per Han 2010 Results section 'Population pharmacokinetic analysis'."
  )

  ini({
    # Structural parameters: Han 2010 base model (no covariates), reported
    # in the Results section paragraph beginning "A two-compartment model
    # with first-order absorption and elimination adequately described the
    # data. The population estimates (interindividual) of bioavailability,
    # clearance, volume of distribution of the central compartment (Vc) and
    # peripheral compartment (Vp), intercompartment clearance (Q), and
    # absorption rate constant (ka) were 45.9% (82.9%), 3.45 liters/h
    # (107%), 54.7 liters (78.4%), 143 liters (88.3%), 22.6 liters/h
    # (50.1%) and 0.591 h-1 (115.2%)."
    lka     <- log(0.591); label("Absorption rate constant (1/h)")          # Han 2010 Results: ka = 0.591 h-1
    lcl     <- log(3.45);  label("Clearance (L/h)")                         # Han 2010 Results: CL = 3.45 L/h
    lvc     <- log(54.7);  label("Central volume of distribution (L)")      # Han 2010 Results: Vc = 54.7 L
    lvp     <- log(143);   label("Peripheral volume of distribution (L)")   # Han 2010 Results: Vp = 143 L
    lq      <- log(22.6);  label("Inter-compartmental clearance (L/h)")     # Han 2010 Results: Q = 22.6 L/h
    lfdepot <- log(0.459); label("Oral bioavailability (fraction)")         # Han 2010 Results: F = 45.9%

    # IIV. Han 2010 Methods states "Correlations between pharmacokinetic
    # parameters were always incorporated and estimated" -- the model
    # estimates an OMEGA BLOCK covariance, but the off-diagonal correlations
    # are not reported numerically anywhere in the paper or table. The
    # diagonal CV% values are reported in the Results section parenthetical
    # adjacent to each typical-value parameter (above). A diagonal IIV with
    # the reported CV% values is used here as the most faithful extraction
    # achievable without fabricating correlation magnitudes; the deviation
    # is documented in the vignette's Assumptions and deviations section.
    # Han 2010 reports CV% (approximate-CV convention) so the internal
    # variance is omega^2 = (CV/100)^2:
    #   ka  CV 115.2% -> 1.152^2 = 1.327
    #   CL  CV 107%   -> 1.07^2  = 1.1449
    #   Vc  CV  78.4% -> 0.784^2 = 0.6147
    #   Vp  CV  88.3% -> 0.883^2 = 0.7797
    #   Q   CV  50.1% -> 0.501^2 = 0.2510
    #   F   CV  82.9% -> 0.829^2 = 0.6872
    etalka     ~ 1.327    # Han 2010 Results: IIV ka  = 115.2% CV
    etalcl     ~ 1.1449   # Han 2010 Results: IIV CL  = 107% CV
    etalvc     ~ 0.6147   # Han 2010 Results: IIV Vc  = 78.4% CV
    etalvp     ~ 0.7797   # Han 2010 Results: IIV Vp  = 88.3% CV
    etalq      ~ 0.2510   # Han 2010 Results: IIV Q   = 50.1% CV
    etalfdepot ~ 0.6872   # Han 2010 Results: IIV F   = 82.9% CV

    # Combined residual error: Han 2010 Methods 'Population pharmacokinetic
    # analysis' specifies the form C_obs = C_pred * (1 + epsilon) + epsilon'
    # and the Results section reports "The proportional and additive
    # residual variability was 0.31 and 0.49 ug/mL, respectively."
    # The Methods text labels the epsilon variances symbolically as
    # sigma^2 and sigma'^2, but the explicit ug/mL unit on the additive
    # value (concentration scale rather than concentration-squared) and
    # the conventional NONMEM popPK reporting practice indicate the
    # printed values are standard deviations rather than variances. The
    # SD interpretation is used here; downstream users who prefer the
    # alternative variance interpretation should set propSd = sqrt(0.31)
    # = 0.557 and addSd = sqrt(0.49) = 0.700 ug/mL. The SD-versus-variance
    # ambiguity is flagged in the vignette's Assumptions and deviations.
    propSd <- 0.31; label("Proportional residual error (fraction)")   # Han 2010 Results: sigma_prop = 0.31
    addSd  <- 0.49; label("Additive residual error (ug/mL)")          # Han 2010 Results: sigma_add  = 0.49 ug/mL
  })

  model({
    # Individual PK parameters: lognormal IIV on each structural parameter.
    ka     <- exp(lka     + etalka)
    cl     <- exp(lcl     + etalcl)
    vc     <- exp(lvc     + etalvc)
    vp     <- exp(lvp     + etalvp)
    q      <- exp(lq      + etalq)
    fdepot <- exp(lfdepot + etalfdepot)

    # Micro-constants for the two-compartment system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment ODE system with first-order oral absorption. The
    # intravenous infusions in the source study were delivered as 2 h
    # infusions directly into the central compartment; the oral doses
    # enter via the depot compartment subject to bioavailability fdepot.
    # Downstream dosing datasets specify the route via the cmt column
    # (cmt = "depot" for oral; cmt = "central" with RATE set for the
    # 2 h intravenous infusion).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Oral bioavailability applies only when the dose enters via the
    # depot compartment; intravenous doses bypass the depot and so are
    # unaffected by f(depot).
    f(depot) <- fdepot

    # Observation. Concentration units mg / L = ug/mL match the paper's
    # plasma-concentration reporting scale (linearity range 0.2-9 ug/mL).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
