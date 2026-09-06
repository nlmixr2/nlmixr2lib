Deng_2024_magnesiumSulfate <- function() {
  description <- "One-compartment population PK model of magnesium sulfate (MgSO4-7H2O) given by intravenous infusion, with creatinine clearance, body mass index and concomitant furosemide effects on clearance and a concomitant furosemide effect on volume, in Chinese women with preeclampsia (Deng 2024)."
  reference <- "Deng J, Peng L, Wang Y, Li J, Tang L, Yu Y. Population pharmacokinetics and dose optimization of magnesium sulfate in Chinese preeclampsia population. BMC Pregnancy Childbirth 2024;24:424. doi:10.1186/s12884-024-06620-x"
  vignette  <- "Deng_2024_magnesiumSulfate"
  units     <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central = list(analyte = "magnesium", units = "mg", specimen = "serum", verified = FALSE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Maternal creatinine clearance. Raw (NOT BSA-normalised) mL/min; Deng 2024 reports it only as 'creatinine clearance (CCR)' and never states the estimating equation, so the assay form is unknown. The cohort is renally hyperfiltrating (mean 182.18 +/- 67.15 mL/min), which the paper attributes to the 40-65% pregnancy-associated rise in glomerular filtration rate.",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL, encoded here as (CRCL/175)^e_crcl_cl with exponent +0.39 from Deng 2024 Table 2 (dCLdCCR). IMPORTANT: the reference value 175 is NOT the constant printed in the paper's own equation. Deng 2024 prints 'CL(L/h) = 2.98 x (CCR/51)^0.39 x (BMI/51)^-0.54 x (1 - 0.16 x furosemide) x exp(etaCL)' -- normalising BOTH covariates by 51, which is the number of subjects (n = 51), not a covariate reference value for either. 175 mL/min is the study's own median CCR, stated verbatim in the Monte Carlo simulation section ('when CCR and BMI are taken at the median (CCR:175 ml/min, BMI:29 kg/m2)'). See the vignette's Assumptions and deviations section for the full argument.",
      source_name        = "CCR"
    ),
    BMI = list(
      description        = "Maternal body mass index at baseline.",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL, encoded here as (BMI/29)^e_bmi_cl with exponent -0.54 from Deng 2024 Table 2 (dCLdBMI). The reference value 29 kg/m^2 is the study median stated in the Monte Carlo simulation section; Deng 2024 Table 1 gives the median as 29.13 (IQR 27.16-33.30), a 0.4% difference that changes the covariate multiplier by 0.24%. As with CRCL, this replaces the '51' printed in the paper's equation. The negative exponent means higher BMI lowers magnesium clearance; the paper flags this as opposite in direction to body-weight effects in other magnesium popPK models and attributes it to preeclampsia-associated water retention.",
      source_name        = "BMI"
    ),
    CONMED_FUROSEMIDE = list(
      description        = "Concomitant furosemide (loop diuretic) administration indicator; 1 = the woman received furosemide during the treatment window, 0 = she did not.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant furosemide)",
      notes              = "16 of 51 women (31.37%) received furosemide (Deng 2024 Table 1). Enters BOTH clearance and volume as a fractional-change term: cl is multiplied by (1 + e_conmed_furosemide_cl * CONMED_FUROSEMIDE) = 0.84 and vc by (1 + e_conmed_furosemide_vc * CONMED_FUROSEMIDE) = 0.75 when furosemide is present, per Deng 2024 Table 2 (dCLdfurosemide = -0.16, dVdfurosemide = -0.25) and the printed final-model equations. Both effects RAISE serum magnesium, which the paper attributes to furosemide being prescribed for the volume overload and renal impairment of progressing preeclampsia rather than to a direct drug interaction; the paper describes this as the first report of a furosemide effect on magnesium popPK. Because the indicator is a marker of disease progression as much as of coadministration, do not transport it to a population in which furosemide is used for another reason.",
      source_name        = "furosemide"
    )
  )

  # Covariates screened by the stepwise covariate model but NOT retained in the
  # final model (Deng 2024 Results: "Age, calcium, albumin, labetalol,
  # nifedipine and other clinical variables did not have statistically
  # significant relationships with PK parameters"). Documented for provenance
  # only; no point estimates are reported for any of them.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Maternal age at enrolment",
      units       = "years",
      type        = "continuous",
      notes       = "Median 31 (IQR 28-35), mean 31.86 +/- 5.22 years (Deng 2024 Table 1). Screened in the stepwise covariate model on CL and V; not retained and no point estimate reported."
    ),
    ALB = list(
      description = "Maternal serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Mean 30.45 +/- 4.66 g/L (Deng 2024 Table 1). Screened in the stepwise covariate model; not retained and no point estimate reported."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 51L,
    n_studies      = 1L,
    age_range      = "18-45 years by inclusion criteria; observed median 31 (IQR 28-35), mean 31.86 +/- 5.22 years",
    age_median     = "31 years",
    weight_range   = "not reported; body size is characterised by BMI only (median 29.13 kg/m^2, IQR 27.16-33.30)",
    sex_female_pct = 100,
    race_ethnicity = "Chinese (women enrolled at the Affiliated Suzhou Hospital of Nanjing Medical University, Suzhou, Jiangsu)",
    disease_state  = "Preeclampsia requiring intravenous MgSO4 for seizure prophylaxis; gestational age 32.31 +/- 3.92 weeks at treatment. Exclusions: myasthenia gravis or other neuromuscular disorder, severe renal insufficiency, hypermagnesaemia, hypocalcaemia, hypokalaemia, heart block, stillbirth, fetal malformation, diabetes mellitus, thyroid disease, intrahepatic cholestasis.",
    dose_range     = "Day 1: 5 g MgSO4-7H2O intravenous loading dose over 30-120 min followed by a 10 g maintenance dose over 6-8 h by infusion pump. Days 2-5: 10 g maintenance dose only. Maternal blood was sampled at 0, 4, 5 and 12 h after the day-2 maintenance dose, so the model was fit to maintenance-dose data with no loading dose in the observation window.",
    regions        = "China (Suzhou, Jiangsu)",
    notes          = "Prospective observational study, April 2021 - April 2023; 199 serum magnesium concentrations from 51 women (2-4 samples each). Estimation in Phoenix NLME 8.3 by FOCE-ELS; model evaluated by 1000-replicate nonparametric bootstrap, goodness-of-fit plots and a visual predictive check. Doses are administered as MgSO4-7H2O (heptahydrate, MW 246.47) and serum concentrations are reported as elemental magnesium; the PK parameters in this file use dose units of mg of elemental Mg and concentration units of mg/L of elemental Mg, matching the sibling models Salinger_2013_magnesiumSulfate.R and Easterling_2018_magnesium_sulfate.R. Convert administered MgSO4-7H2O grams to mg Mg by multiplying by 24.305/246.47 = 0.0986 (5 g = 493.1 mg Mg, 10 g = 986.2 mg Mg, 15 g = 1479.3 mg Mg). The mg/L concentration unit is fixed by the Deng 2024 Fig. 3 visual-predictive-check y-axis ('Prediction-corrected concentrations (mg/L)', observations spanning roughly 20-52 mg/L); the additive residual SD of 3.65 is on that same mg/L scale. Additional covariates screened and rejected: serum calcium (2.17 +/- 0.15, tabulated as g/L but almost certainly mmol/L), labetalol (50/51, 98.03%), nifedipine (22/51, 43.14%), adverse-reaction occurrence (21/51, 41.18%), gestational age, diagnosis and mode of delivery (cesarean 40, vaginal 4, rivanol-induced abortion 7). IMPORTANT -- this model carries NO endogenous magnesium baseline term, because Deng 2024 Table 2 estimates none; it therefore predicts administered-magnesium concentrations only, whereas the data it was fit to are total serum magnesium. The paper reports an observed baseline of 0.76 mmol/L (18.5 mg/L, IQR 0.71-0.86) and its own Fig. 4 simulations visibly start every curve at 0.80 mmol/L, so the authors added a baseline post hoc that is not a parameter of the published model. Both sibling magnesium models DO estimate a baseline (Salinger 2013 BL = 20.8 mg/L, Easterling 2018 BL = 22.48 mg/L). Add roughly 18.5 mg/L to Cc before comparing this model against measured total serum magnesium. See the vignette's Assumptions and deviations section."
  )

  ini({
    # Structural parameters from Deng 2024 Table 2 (final model estimates).
    # Reference subject: CCR = 175 mL/min, BMI = 29 kg/m^2, no furosemide.
    # Dose units used here: mg of elemental Mg (multiply g of MgSO4-7H2O by 98.6).
    # Concentration units: mg/L of elemental Mg (administered drug only; see
    # the population notes on the absent endogenous baseline).
    lcl <- log(2.98);  label("Clearance for the reference subject (CL, L/h)")             # Deng 2024 Table 2 tvCL, 95% CI 1.29-4.62
    lvc <- log(25.07); label("Volume of distribution for the reference subject (V, L)")   # Deng 2024 Table 2 tvV, 95% CI 23.31-26.82

    # Covariate effects from Deng 2024 Table 2 and the final-model equations
    # printed in Results ("The equations for V and CL were as follows").
    # NOTE on the two reference values: the printed equation normalises BOTH
    # CCR and BMI by 51 -- the subject count, not a covariate reference. Taken
    # literally it puts clearance at the study median covariates at
    # 2.98 * (175/51)^0.39 * (29/51)^-0.54 = 6.54 L/h, which contradicts
    # Table 2's own "tv" (typical value) label, the Abstract ("CL was
    # 2.98 L/h"), and the Discussion's comparison of 2.98 against a literature
    # range of 1.38-5.00 L/h ("lower than the values typically described").
    # The study medians 175 mL/min and 29 kg/m^2 are used instead; the paper
    # states them verbatim in the Monte Carlo section. This choice is
    # corroborated by the SHAPE of the paper's own Fig. 4: dividing the
    # modelled by the published drug-attributable concentration at t = 1..7 h
    # gives a near-constant ratio (0.659, CV 1.0% without furosemide; 0.633,
    # CV 1.6% with), i.e. the profiles differ by amplitude alone. Under the
    # literal /51 reading the same ratio drifts monotonically from 0.61 to
    # 0.40 (CV 15%), so the shape is wrong. Fig. 4's amplitude is NOT
    # reproducible under either reading; see the vignette Errata.
    e_crcl_cl              <-  0.39; label("Power exponent of (CRCL/175) on CL (unitless)")                       # Deng 2024 Table 2 dCLdCCR, 95% CI 0.31-0.46
    e_bmi_cl               <- -0.54; label("Power exponent of (BMI/29) on CL (unitless)")                         # Deng 2024 Table 2 dCLdBMI, 95% CI (-0.70)-(-0.38)
    e_conmed_furosemide_cl <- -0.16; label("Fractional change in CL with concomitant furosemide (unitless)")      # Deng 2024 Table 2 dCLdfurosemide, 95% CI (-0.22)-(-0.096)
    e_conmed_furosemide_vc <- -0.25; label("Fractional change in V with concomitant furosemide (unitless)")       # Deng 2024 Table 2 dVdfurosemide, 95% CI (-0.34)-(-0.16)

    # Inter-individual variability from Deng 2024 Table 2. The table labels
    # these rows "omega^2 V" and "omega^2 CL", i.e. VARIANCES on the
    # log scale, giving CV 15.3% on V and 29.2% on CL. The adjacent "CV (%)"
    # column (26.89 for V, 17.38 for CL) is the relative standard error of the
    # variance estimate, not the IIV CV: were it the IIV CV, the smaller
    # variance (V) would have to carry the smaller CV, and it does not.
    # Reported shrinkage 13.04% (V) and 13.05% (CL). No eta correlation is
    # reported, so the two are coded independently.
    etalvc ~ 0.023; label("IIV on volume of distribution (variance, log scale)")  # Deng 2024 Table 2 omega^2 V
    etalcl ~ 0.082; label("IIV on clearance (variance, log scale)")               # Deng 2024 Table 2 omega^2 CL

    # Residual error from Deng 2024 Table 2 ("stdev0", 95% CI 3.20-4.09). The
    # paper selected an additive error model over proportional and combined
    # alternatives on AIC/BIC/OFV grounds; the SD is on the mg/L scale of the
    # Fig. 3 visual predictive check.
    addSd <- 3.65; label("Additive residual error (SD, mg/L)")  # Deng 2024 Table 2 stdev0
  })
  model({
    # Individual PK parameters, replicating the final-model equations printed
    # in Deng 2024 Results with the covariate reference values discussed in
    # ini() above:
    #   CL_i = 2.98 * (CCR_i/175)^0.39 * (BMI_i/29)^-0.54
    #          * (1 - 0.16 * furosemide) * exp(etaCL)
    #   V_i  = 25.07 * (1 - 0.25 * furosemide) * exp(etaV)
    cl <- exp(lcl + etalcl) *
      (CRCL / 175)^e_crcl_cl *
      (BMI / 29)^e_bmi_cl *
      (1 + e_conmed_furosemide_cl * CONMED_FUROSEMIDE)
    vc <- exp(lvc + etalvc) *
      (1 + e_conmed_furosemide_vc * CONMED_FUROSEMIDE)

    kel <- cl / vc

    # One-compartment disposition. Magnesium sulfate is given intravenously
    # only (loading infusion and maintenance infusion), so all doses go into
    # central via the events file with a rate or duration set per record.
    d/dt(central) <- -kel * central

    # Observation: administered magnesium only. Deng 2024 Table 2 estimates no
    # endogenous baseline, so add ~18.5 mg/L (the reported 0.76 mmol/L
    # baseline) before comparing against measured total serum magnesium.
    Cc <- central / vc

    Cc ~ add(addSd)
  })
}
