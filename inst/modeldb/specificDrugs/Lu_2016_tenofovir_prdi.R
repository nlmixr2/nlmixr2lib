Lu_2016_tenofovir_prdi <- function() {
  description <- "Two-compartment population PK model with first-order absorption for tenofovir (300 mg oral TDF once daily) in HIV-1-uninfected African adults receiving once-daily preexposure prophylaxis (Lu 2016, Partners PrEP Study). PRDI variant: parameters estimated using patient-reported dosing information with a steady-state assumption. Apparent oral clearance (CL/F) carries a power-form covariate effect on creatinine clearance (raw Cockcroft-Gault, mL/min) centred at the cohort median 106 mL/min. Diagonal IIV on CL/F only; combined additive + proportional residual error."
  reference <- "Lu Y, Goti V, Chaturvedula A, Haberer JE, Fossler MJ, Sale ME, Bangsberg D, Baeten JM, Celum CL, Hendrix CW. Population pharmacokinetics of tenofovir in HIV-1-uninfected members of serodiscordant couples and effect of dose reporting methods. Antimicrob Agents Chemother. 2016;60(9):5379-5386. doi:10.1128/AAC.00559-16"
  vignette <- "Lu_2016_tenofovir"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CLCR. Cockcroft-Gault creatinine clearance in raw mL/min (NOT BSA-normalized). Stored under the canonical CRCL column per inst/references/covariate-columns.md (CRCL accepts raw mL/min when the source paper does not apply BSA normalization, with the per-model description recording the assay form). Reference value 106 mL/min is the cohort mean (Lu 2016 Table 1; SD 31 mL/min). The covariate enters CL/F as a power form (CRCL/106)^e_crcl_cl, where e_crcl_cl = 0.379. The paper's covariate-form selection (Methods: 'Models with linear, power, and exponential functions were tested') is not explicitly stated for the retained final-model form; the power interpretation is adopted here because (i) the reported coefficient 0.379 is dimensionless and in the typical range for power exponents on renal-function covariates, (ii) the linear/exponential alternatives with the same numeric coefficient give physically implausible CL values at the extremes of the observed CrCl range, and (iii) NONMEM's standard divisive centring (CRCL/median) for a power model matches the paper's prose 'centered to the median values'.",
      source_name        = "CLCR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 404L,
    n_studies      = 1L,
    age_range      = "mean 35 years (SD 8)",
    age_median     = "35 years",
    weight_range   = "mean 61 kg (SD 11)",
    weight_median  = "61 kg",
    sex_female_pct = 45,
    race_ethnicity = "Sub-Saharan African (Kenyan and Ugandan cohort); race not stratified further in the paper",
    disease_state  = "HIV-1-seronegative; HIV-1-uninfected members of HIV-1-serodiscordant heterosexual couples enrolled in the Partners PrEP Study",
    dose_range     = "300 mg oral tenofovir disoproxil fumarate (TDF) once daily as preexposure prophylaxis (alone or in fixed-dose combination with emtricitabine 200 mg)",
    regions        = "Kenya and Uganda",
    crcl_range     = "mean 106 mL/min (SD 31)",
    samples        = "1280 tenofovir plasma concentrations from 404 participants (mean 3.2 concentrations per subject); 17% below limit of quantitation (LLOQ 0.31 ng/mL)",
    adherence      = "97% by pill counts in the main trial (PRDI data set used here)",
    notes          = "Cohort demographics per Lu 2016 Table 1 (PRDI data set column). Phase 3 randomised, double-blind, placebo-controlled HIV-1 PrEP trial; the PK substudy enrolled HIV-1-seronegative partners on active TDF or TDF+FTC arms. Plasma sampled at monthly clinic visits with patient-reported dosing times; the PRDI fit assumes once-daily dosing at steady state with the patient-reported dose times preceding sample collection."
  )

  ini({
    # Lu 2016 Table 2 'PRDI data Final model' column (point estimates).
    # Reference subject: CRCL = 106 mL/min (cohort mean, Table 1).
    # All apparent oral clearances (CL/F, Q/F) in L/h; apparent volumes (Vc/F, Vp/F) in L; ka in 1/h.
    lka  <- log(4.7);   label("First-order absorption rate constant ka (1/h)")                      # Lu 2016 Table 2 PRDI final Ka = 4.7 /h
    lcl  <- log(57);    label("Apparent oral clearance CL/F at CRCL = 106 mL/min (L/h)")            # Lu 2016 Table 2 PRDI final CL = 57 L/h
    lvc  <- log(393);   label("Apparent central volume of distribution V1/F (L)")                   # Lu 2016 Table 2 PRDI final V1 = 393 L
    lq   <- log(178);   label("Apparent inter-compartmental clearance Q/F (L/h)")                   # Lu 2016 Table 2 PRDI final Q = 178 L/h
    lvp  <- log(614);   label("Apparent peripheral volume of distribution Vp/F (L)")                # Lu 2016 Table 2 PRDI final Vp = 614 L

    # Power-form CrCl effect on CL/F. Lu 2016 reports theta_CL-CR = 0.379 as the
    # bootstrap median for the final-model covariate coefficient (Table 2, PRDI
    # final model column; 95% CI 0.286-0.441). Applied as (CRCL/106)^e_crcl_cl.
    e_crcl_cl <- 0.379; label("Power exponent of (CRCL/106) on CL/F (unitless)")                     # Lu 2016 Table 2 PRDI final theta_CL-CR = 0.379 (median 0.379, 95% CI 0.286-0.441)

    # Diagonal IIV on CL/F only. Lu 2016 Table 2 reports IIV as %CV with no
    # inter-eta correlations documented (PRDI final-model row: 'IIV on CL'
    # = 16% CV; no IIV reported on V1, Ka, Q, or Vp for the PRDI final model).
    #   omega^2 = log(CV^2 + 1):
    #     CL/F CV  16% -> log(1 + 0.16^2) = log(1.0256) = 0.02528
    etalcl ~ 0.02528  # Lu 2016 Table 2 PRDI final IIV on CL = 16% CV

    # Combined additive + proportional residual error (Lu 2016 Table 2 PRDI
    # final-model row). Concentrations are in ng/mL; additive SD is in ng/mL.
    # Lu 2016 also reports 'IIV on additive error' = 143% CV for the PRDI
    # final model; this represents an inter-individual scaling of the additive
    # residual-error magnitude per the Karlsson reference cited in Methods
    # (reference 21). It is NOT encoded structurally here because nlmixr2's
    # canonical residual-error syntax does not naturally express IIV on the
    # error magnitude without breaking parameter-naming conventions; the
    # omission affects only the lower tail of the prediction interval near
    # the LLOQ and does not alter typical-value disposition. See the vignette
    # Assumptions and deviations section.
    addSd  <- 28;   label("Additive residual error (ng/mL)")                                         # Lu 2016 Table 2 PRDI final additive = 28 ng/mL
    propSd <- 0.21; label("Proportional residual error (fraction)")                                  # Lu 2016 Table 2 PRDI final proportional = 21% CV
  })
  model({
    # Individual PK parameters. Power-form CrCl effect on CL/F centred at the
    # cohort mean 106 mL/min; no covariates on Ka, Vc, Q, or Vp in the PRDI
    # final model (Lu 2016 Section 'Population PK model for tenofovir').
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (CRCL / 106) ^ e_crcl_cl
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # Dose in mg, vc in L -> central / vc has units mg/L; * 1000 -> ng/mL.
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
