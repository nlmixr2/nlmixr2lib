Brekkan_2016_factorIX <- function() {
  description <- "Three-compartment population PK model for plasma-derived factor IX (FIX) activity in patients with moderate or severe haemophilia B, developed by Brekkan et al. 2016 to support pharmacokinetic dose individualisation. Disposition is described by linear three-compartment kinetics with intravenous input and first-order elimination from the central compartment; allometric body-weight scaling on CL/Q (0.75) and V1/V2/V3 (1.0) is fixed with a reference weight of 70 kg, and an endogenous baseline FIX activity is estimated as a structural parameter."
  reference <- "Brekkan A, Berntorp E, Jensen K, Nielsen EI, Jonsson S. Population pharmacokinetics of plasma-derived factor IX: procedures for dose individualization. J Thromb Haemost. 2016;14(4):724-732. doi:10.1111/jth.13271"
  vignette <- "Brekkan_2016_factorIX"
  units <- list(time = "hour", dosing = "U", concentration = "U/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on CL, Q2, Q3 (exponent fixed at 0.75) and on V1, V2, V3 (exponent fixed at 1.0) with reference weight 70 kg (Brekkan 2016 Table 2 footnote and Equation 1). Among three approaches tested for the CL exponent (estimated, fixed at 0.75, fixed to the original published 1.26), 0.75 was retained because all three gave similar fits.",
      source_name        = "WT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age in years",
      units       = "years",
      type        = "continuous",
      notes       = "Age was tested as a covariate on CL using a fractional change per year deviation from median age (Brekkan 2016 Equation 2) but was not statistically significant (P < 0.05) and was not retained in the final model. Documented here only to record that it was screened."
    ),
    PRODUCT = list(
      description = "Factor IX product identity",
      units       = NA_character_,
      type        = "categorical",
      notes       = "Seven plasma-derived FIX products were tested as a categorical covariate on CL: AlphaNine (reference, most common), Factor IX Grifols, Immunine, Octanine, Nanotiv, Preconativ, Mononine. The likelihood-ratio test was statistically significant (P < 0.01, df = 6) but the largest deviation from AlphaNine was within +/- 20% of the typical CL (range 252-378 mL/h vs. 315 mL/h reference), so the clinical-significance criterion was not met and the product effect was dropped from the final model (Brekkan 2016 Results, p. 727). Documented here to record the screen."
    ),
    OCC = list(
      description = "Occasion identifier",
      units       = NA_character_,
      type        = "count",
      notes       = "Inter-occasion variability (IOV) on CL and V1 with a correlation between the two random effects was retained in the published model (CV 21.4% on CL, 20.1% on V1, correlation 0.902; Brekkan 2016 Table 2). The static library model omits IOV because there is no occasion variable on the standard library event grid; this is a deliberate simplification documented in the vignette's Assumptions and deviations section."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 34L,
    n_studies      = 5L,
    age_range      = "Not stated as range; mean 27.5 (SD 10.9) years across all studies; per-study means 23.1-42.8 years",
    age_median     = "Not reported",
    weight_range   = "Not stated as range; mean 66.8 (SD 13.5) kg across all studies; per-study means 61.0-69.7 kg",
    weight_median  = "Not reported; reference weight in the model is 70 kg",
    sex_female_pct = 0,
    race_ethnicity = "Not reported",
    disease_state  = "Moderate or severe haemophilia B; 35 of 1794 samples had FIX activity below 0.01 U/mL (the common lower limit of quantification) and were retained.",
    dose_range     = "Single-dose PK following intravenous administration of plasma-derived factor IX concentrate; specific dose ranges per study are not summarised in the paper.",
    regions        = "Pooled across five previously published studies (Lissitchkov et al., Aznar et al., Berntorp et al., Bjorkman et al., Carlsson et al.); the latter three plus two additional unpublished studies fed the original model of Berntorp et al. [reference 15 in Brekkan 2016].",
    n_samples      = "1,794 FIX activity samples across 34 unique patients; 7-17 samples per patient per occasion. Several patients contributed to more than one study.",
    products       = "Seven plasma-derived FIX products: AlphaNine (most common; Grifols), Factor IX Grifols, Immunine (Baxter/Immuno), Octanine (Octapharma), Nanotiv (Pharmacia), Preconativ (Pharmacia), Mononine (Armour Pharmaceutical).",
    notes          = "Patient demographics summarised in Brekkan 2016 Table 1. The model is intended for plasma-derived FIX products only; the Discussion notes the final model should be used with caution to describe FIX activity following administration of products that are not plasma-derived. Haemophilia B is X-linked, so the cohort is essentially all male."
  )

  ini({
    # Structural parameters - typical values for a 70 kg reference subject
    # (Brekkan 2016 Table 2). CL and Q in mL/h; V1, V2, V3 in mL; baseline
    # FIX activity in U/mL.
    lcl    <- log(319.8)  ; label("Clearance for the reference 70 kg subject (CL, mL/h)")                     # Brekkan 2016 Table 2: CL = 319.8 mL/h
    lvc    <- log(5922)   ; label("Central volume of distribution for the reference 70 kg subject (V1, mL)") # Brekkan 2016 Table 2: V1 = 5922 mL
    lvp    <- log(828.9)  ; label("Peripheral volume of distribution 1 for the reference 70 kg subject (V2, mL)") # Brekkan 2016 Table 2: V2 = 828.9 mL
    lq     <- log(1049)   ; label("Intercompartmental clearance to peripheral 1 (Q2, mL/h)")                 # Brekkan 2016 Table 2: Q2 = 1049 mL/h
    lvp2   <- log(2234)   ; label("Peripheral volume of distribution 2 for the reference 70 kg subject (V3, mL)") # Brekkan 2016 Table 2: V3 = 2234 mL
    lq2    <- log(160.4)  ; label("Intercompartmental clearance to peripheral 2 (Q3, mL/h)")                 # Brekkan 2016 Table 2: Q3 = 160.4 mL/h

    # Endogenous baseline FIX activity. Brekkan 2016 estimates this as a
    # structural parameter rather than subtracting an empirical baseline before
    # fitting. The simulation study set baseline = 0 to represent severe
    # haemophilia B; library users can do the same by setting rbase = 0 in
    # rxode2.
    lrbase <- log(0.01588); label("Endogenous baseline factor IX activity (U/mL)")                          # Brekkan 2016 Table 2: baseline = 0.01588 U/mL

    # Allometric body-weight exponents fixed at canonical values (Brekkan 2016
    # Methods, p. 725: "The exponent was fixed to 1 and 0.75 for volume and
    # intercompartmental CL parameters, respectively" -- and Results, p. 727:
    # "The allometric exponent for body weight on CL was fixed to 0.75 since
    # the three approaches evaluated [estimated, 1.26, 0.75] resulted in
    # similar model fits"). Reference weight 70 kg.
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL and Q (unitless)")              # Brekkan 2016 Results, p. 727: CL exponent fixed at 0.75
    e_wt_vc <- fixed(1.00); label("Allometric exponent of body weight on V1, V2, V3 (unitless)")            # Brekkan 2016 Methods, p. 725: volume exponent fixed at 1

    # Inter-individual variability. Brekkan 2016 Table 2 footnote describes
    # the IIV column as "coefficient of variation of interindividual
    # variability of clearance, volumes and baseline", with lognormal
    # parameters (Methods, p. 725). Per the convention used in sister factor
    # IX papers (Diao 2014 explicit footnote: "sqrt(variance) * 100"), the
    # reported value is the SD of the log-scale eta directly, so
    # omega^2 = (reported_value)^2. Covariances are derived as
    # corr * sqrt(var1) * sqrt(var2).
    #
    # CL/V1 correlated block (Brekkan 2016 Table 2):
    #   omega^2_CL = 0.127^2 = 0.016129
    #   cov_CL_V1  = 0.705 * 0.127 * 0.157 = 0.014053
    #   omega^2_V1 = 0.157^2 = 0.024649
    etalcl + etalvc ~ c(0.016129,
                        0.014053, 0.024649)        # Brekkan 2016 Table 2: IIV CL = 0.127, IIV V1 = 0.157, corr(CL,V1) = 0.705

    # V2/V3 correlated block (Brekkan 2016 Table 2):
    #   omega^2_V2 = 0.667^2 = 0.444889
    #   cov_V2_V3  = 0.814 * 0.667 * 1.020 = 0.553828
    #   omega^2_V3 = 1.020^2 = 1.040400
    etalvp + etalvp2 ~ c(0.444889,
                         0.553828, 1.040400)       # Brekkan 2016 Table 2: IIV V2 = 0.667, IIV V3 = 1.020, corr(V2,V3) = 0.814

    # Independent IIV on baseline (Brekkan 2016 Table 2):
    #   omega^2_baseline = 0.207^2 = 0.042849
    # Eta shrinkage on baseline was 34.3% per Table 2 footnote; eta shrinkage
    # on all other IIV parameters was < 22%.
    etalrbase ~ 0.042849                            # Brekkan 2016 Table 2: IIV baseline = 0.207

    # Residual error - combined additive + proportional (Brekkan 2016 Table 2,
    # footnotes section dagger/double-dagger): additive SD reported as 0.0067
    # U/mL, proportional reported as CV 0.0695 (used directly as propSd in
    # nlmixr2's linear-space proportional form). Epsilon shrinkage was < 10%.
    addSd  <- 0.0067 ; label("Additive residual error on FIX activity (U/mL)")    # Brekkan 2016 Table 2: additive error SD = 0.0067 U/mL
    propSd <- 0.0695 ; label("Proportional residual error (fraction)")            # Brekkan 2016 Table 2: proportional error CV = 0.0695
  })

  model({
    # Individual PK parameters with allometric scaling (Brekkan 2016 Eq. 1
    # and Table 2 footnote). Reference weight 70 kg; exponents fixed at 0.75
    # for CL/Q2/Q3 and 1.0 for V1/V2/V3. IIV is on CL, V1, V2, V3, and
    # baseline; the paper retained no IIV on Q2 or Q3 after model
    # simplification (Brekkan 2016 Results, p. 727: "the original model could
    # be reduced by four parameters related to the IIV part").
    cl    <- exp(lcl  + etalcl)  * (WT / 70)^e_wt_cl
    vc    <- exp(lvc  + etalvc)  * (WT / 70)^e_wt_vc
    vp    <- exp(lvp  + etalvp)  * (WT / 70)^e_wt_vc
    q     <- exp(lq)             * (WT / 70)^e_wt_cl
    vp2   <- exp(lvp2 + etalvp2) * (WT / 70)^e_wt_vc
    q2    <- exp(lq2)            * (WT / 70)^e_wt_cl
    rbase <- exp(lrbase + etalrbase)

    # Micro-constants for the three-compartment model with first-order
    # elimination from the central compartment.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # ODE system. Dose enters central (FIX is given IV by short infusion);
    # peripheral1 maps to the paper's V2/Q2 compartment and peripheral2 maps
    # to V3/Q3 (Brekkan 2016 Material and methods, p. 725).
    d/dt(central)     <- -(kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # Observed FIX activity in U/mL. Brekkan 2016 estimated the baseline FIX
    # activity (endogenous production) as a structural parameter added to the
    # model-predicted concentration. To represent a patient with severe
    # haemophilia B, set rbase = 0 in the simulation (the paper does this in
    # the dose-individualisation simulations, p. 726).
    Cc <- central / vc + rbase
    Cc ~ add(addSd) + prop(propSd)
  })
}
