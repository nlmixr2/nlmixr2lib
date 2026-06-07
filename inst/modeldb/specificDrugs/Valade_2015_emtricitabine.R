Valade_2015_emtricitabine <- function() {
  description <- "Two-compartment oral population PK model for emtricitabine (FTC) in HIV-1-infected men on combined antiretroviral therapy, with an asymmetric effect compartment of negligible volume describing seminal plasma distribution via distinct blood-plasma-to-seminal-plasma transfer rate (k1e) and seminal-plasma elimination rate (ke1) constants (Valade 2015, EVARIST ANRS-EP 49 study)"
  reference <- "Valade E, Treluyer JM, Illamola SM, Bouazza N, Foissac F, De Sousa Mendes M, Lui G, Chenevier-Gobeaux C, Suzan-Monti M, Rouzioux C, Assoumou L, Viard JP, Hirt D, Urien S, Ghosn J, for the Evarist ANRS-EP 49 Study Group. Emtricitabine seminal plasma and blood plasma population pharmacokinetics in HIV-infected men in the EVARIST ANRS-EP 49 study. Antimicrob Agents Chemother. 2015;59(11):6800-6806. doi:10.1128/AAC.01517-15"
  vignette <- "Valade_2015_emtricitabine"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CLCR. Computed by the Cockcroft-Gault equation in raw mL/min (NOT BSA-normalized to mL/min/1.73 m^2). Stored under the canonical CRCL column per inst/references/covariate-columns.md (CRCL accepts raw mL/min when the source paper does not apply BSA normalization, with the per-model description recording the assay form). Reference value 113 mL/min (population median, Valade 2015 Table 1). Applied as a power-scaling covariate on CL/F: CL/F = 14.8 * (CRCL / 113)^0.178.",
      source_name        = "CLCR"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = "Tested as a continuous covariate on CL/F, Vc/F, Q/F, Vp/F, k1e, and ke1 (Valade 2015 Methods; equation `param = theta * (CO / median_CO)^beta_CO`). Not retained in the final model: after inclusion of CRCL on CL/F, no further covariate decreased OFV by >= 3.84 (Results, paragraph following Table 2)."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Tested as a continuous covariate on CL/F, Vc/F, Q/F, Vp/F, k1e, and ke1 (Valade 2015 Methods). Not retained in the final model (Results)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Tested as a continuous covariate on FTC PK parameters (Valade 2015 Methods). Not retained in the final model; CRCL (Cockcroft-Gault) was retained instead on CL/F (Results)."
    ),
    CONMED_ARV = list(
      description = "Associated antiretroviral drug regimen (binary indicators for EFV, NVP, ETR, RAL, DRV/r, ATZ/r, LPV/r, IP/r, other combination)",
      units       = "(binary set)",
      type        = "categorical",
      notes       = "Binary co-medication covariates tested per the equation `param = theta * beta_CO^CO` (Valade 2015 Methods). None retained in the final model; in particular, no co-medication explained the IIV on k1e (Results paragraph on FTC distribution variability)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 122L,
    n_studies      = 1L,
    age_range      = "27-63 years",
    age_median     = "43 years",
    weight_range   = "46-108 kg",
    weight_median  = "73 kg",
    sex_female_pct = 0,
    race_ethnicity = "Not reported (Paris-region MSM cohort recruited at 6 clinical centres)",
    disease_state  = "HIV-1-infected men who have sex with men, stable combined antiretroviral therapy for >= 3 months, suppressed blood-plasma HIV RNA load (< 50 copies/mL) for >= 6 months",
    dose_range     = "Daily 200 mg oral emtricitabine combined with tenofovir disoproxil fumarate (TDF) plus one of: EFV, NVP, ETR, RAL, DRV/r, ATZ/r, LPV/r, IP/r, or other combination",
    regions        = "France (6 clinical centres in Paris and nearby suburbs: Hopital Hotel Dieu, Hopital Bichat, Hopital Bicetre, Hopital Foch, Hopital Saint-Louis, Hopital Lariboisiere)",
    creatinine_range = "Serum creatinine median 78 umol/L (range 28-113); Cockcroft-Gault CrCl median 113 mL/min (range 67-368)",
    n_observations   = "236 blood-plasma FTC concentrations + 209 seminal-plasma FTC concentrations at steady state; one BP concentration < LLOQ (0.00625 mg/L) was set to half-LLOQ (Valade 2015 Results).",
    notes          = "EVARIST ANRS EP 49 study (Ghosn et al.). All-male cohort (HIV-1-infected MSM), so sex_female_pct is 0. Baseline demographics per Valade 2015 Table 1. Two visits per patient (inclusion and 1 month after) at varied time intervals between FTC intake and sampling."
  )

  ini({
    # Structural PK parameters (Valade 2015 Table 2 final BP + SP model).
    # ka was not estimable in this dataset; fixed to 0.53 1/h from a previously
    # reported adult FTC popPK estimate (Valade 2015 Results paragraph 2 and
    # reference 21).
    lka  <- fixed(log(0.53)); label("Absorption rate constant ka (1/h, FIXED from previously reported adult FTC popPK estimate, Valade 2015 reference 21)") # Valade 2015 Table 2: ka = 0.53 1/h (asterisk = fixed)
    lcl  <- log(14.8);        label("Apparent oral clearance CL/F (L/h) at median CRCL = 113 mL/min") # Valade 2015 Table 2: CL/F = 14.8 L/h (RSE 4%)
    lvc  <- log(51.6);        label("Apparent central volume of distribution Vc/F (L)")              # Valade 2015 Table 2: Vc/F = 51.6 L (RSE 11%)
    lq   <- log(8.19);        label("Apparent inter-compartmental clearance Q/F (L/h)")              # Valade 2015 Table 2: Q/F = 8.19 L/h (RSE 26%)
    lvp  <- log(106);         label("Apparent peripheral volume of distribution Vp/F (L)")           # Valade 2015 Table 2: Vp/F = 106 L (RSE 44%)

    # Asymmetric effect-compartment rate constants for seminal-plasma distribution.
    # The effect compartment has negligible volume and does not modify the BP
    # disposition. Valade 2015 distinguishes the input rate (k1e) from the output
    # rate (ke1) -- in contrast to the Holford-Sheiner symmetric form (k1e = ke0)
    # used by existing nlmixr2lib effect-compartment models (Park 2001, Berges
    # 2015, Przybylowski 2015). Paper notation preserved (lk1e, lke1) to keep
    # the model readable side-by-side with Valade 2015 Table 2 and Eq. 4.
    lk1e <- log(0.341);       label("Blood-plasma-to-seminal-plasma transfer rate constant k1e (1/h)") # Valade 2015 Table 2: k1e = 0.341 1/h (RSE 28%)
    lke1 <- log(0.113);       label("Seminal-plasma elimination rate constant ke1 (1/h)")             # Valade 2015 Table 2: ke1 = 0.113 1/h (RSE 27%)

    # Covariate effect: power-scaling of CL/F on CRCL with reference 113 mL/min
    # (Valade 2015 final covariate submodel for BP PK).
    e_crcl_cl <- 0.178; label("Power-scaling exponent of CRCL on CL/F (unitless)") # Valade 2015 Table 2: beta_CLCR on CL/F = 0.178 (RSE 35%); CL/F = 14.8 * (CRCL/113)^0.178

    # Inter-individual variability (Valade 2015 Table 2; exponential IIV in
    # Monolix 4.1.4 with reported omega = SD on the log scale). The text
    # describes the IIV on k1e as 53.3 percent (Results "FTC distribution in
    # the male genital tract was variable"), consistent with omega = 0.533
    # under the small-omega approximation CV ~ omega (exact log-normal
    # CV = sqrt(exp(omega^2) - 1) gives 57.5 percent at omega = 0.533).
    etalcl  ~ 0.065025  # Valade 2015 Table 2: omega_CL/F = 0.255 -> variance = 0.255^2 = 0.065025; RSE 12%
    etalk1e ~ 0.284089  # Valade 2015 Table 2: omega_k1e  = 0.533 -> variance = 0.533^2 = 0.284089; RSE 10%

    # Residual error: proportional error model on BP (Cc) and SP (Cseminal)
    # observations (Valade 2015 Results: "Residual variabilities for both blood
    # plasma and seminal plasma were best described by proportional error models").
    # In Monolix the proportional residual y_obs = y_pred * (1 + eps) with
    # eps ~ N(0, sigma^2) maps directly to nlmixr2 `prop(propSd)` with
    # propSd = sigma.
    propSd          <- 0.339; label("Proportional residual error on blood plasma Cc (fraction)")       # Valade 2015 Table 2: sigma_plasma = 0.339 (RSE 6%)
    propSd_Cseminal <- 0.357; label("Proportional residual error on seminal plasma Cseminal (fraction)") # Valade 2015 Table 2: sigma_seminal_plasma = 0.357 (RSE 7%)
  })
  model({
    # Individual PK parameters
    ka  <- exp(lka)
    cl  <- exp(lcl + etalcl) * (CRCL / 113)^e_crcl_cl
    vc  <- exp(lvc)
    q   <- exp(lq)
    vp  <- exp(lvp)
    k1e <- exp(lk1e + etalk1e)
    ke1 <- exp(lke1)

    # Micro-rate constants (Valade 2015 Eq. notation: k10 = CL/Vc, k12 = Q/Vc,
    # k21 = Q/Vp)
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment oral PK (Valade 2015 ODEs, paragraph following Fig. 1).
    # A(1) = depot (gut); A(2) = central; A(3) = peripheral.
    d/dt(depot)       <- -ka  * depot
    d/dt(central)     <-  ka  * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central                                - k21 * peripheral1

    # Effect compartment for seminal plasma. Valade 2015 explicitly defines Ae
    # (the effect-compartment state) as the FTC concentration in seminal
    # plasma (negligible volume), so the state is parameterised in concentration
    # units directly: dAe/dt = k1e * (A(2)/Vc) - ke1 * Ae.
    d/dt(effect)      <-  k1e * (central / vc) - ke1 * effect

    # Observations. Cc is the blood-plasma concentration (mg/L: dose mg / Vc L);
    # Cseminal is the seminal-plasma concentration (already in mg/L because
    # the effect state is a concentration).
    Cc       <- central / vc
    Cseminal <- effect

    Cc       ~ prop(propSd)
    Cseminal ~ prop(propSd_Cseminal)
  })
}
