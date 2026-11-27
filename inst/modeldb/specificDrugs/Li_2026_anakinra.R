Li_2026_anakinra <- function() {
  description <- "One-compartment population PK model of intravenous anakinra (recombinant IL-1 receptor antagonist) in extremely preterm neonates (24-28 weeks gestational age) over the first 21 days of life (Li 2026, Anakinra Pilot trial). Linear elimination; clearance and volume scale allometrically with total body mass (fixed exponents 0.75 and 1, 70 kg reference), and clearance matures with postnatal age from 44.6% of its mature value at birth with a maturation half-time of 1.8 days. The assay does not separate anakinra from endogenous IL-1Ra, so the predicted concentration adds a static endogenous IL-1Ra baseline (6640 pg/mL) to the anakinra concentration. Log-normal IIV on CL and on the baseline; proportional residual error. Includes the assumed direct-effect IL-1 inhibition model (Imax = 1, IC50 = 0.04 mg/L) the authors used for dose-finding simulations."
  reference <- "Li J, Metz DK, Green E, et al. Model-Informed Dose Optimization of Anakinra in Preterm Neonates Implications for Safe and Effective Dosing. Clin Pharmacol Ther. 2026. doi:10.1002/cpt.70491"
  vignette <- "Li_2026_anakinra"

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "anakinra", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying (control stream WTKG; missing values carried forward). Allometric size on CL (exponent 0.75) and V (exponent 1), standardised to 70 kg (Eq. 6-7). Cohort baseline median 0.90 kg (0.54-1.30), 0.53-1.60 kg during therapy (Table 1).",
      source_name = "WTKG"
    ),
    PNA = list(
      description = "Postnatal age",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. Canonical units are months; the paper's maturation function is in DAYS (control stream PNAD), so model() converts with PNA * 30.4375 days per month. PNA = 0 at birth. The function describes only the first 21 days of life (cohort range 0-20 days); it reaches 98.8% of the mature CL by day 10 and is not intended for extrapolation beyond the neonatal period. The trial dataset records PNA as the day of life (whole days), and the authors' simulated exposures (Figure 4) are reproduced with PNA held at the day of life over each day; a continuously advancing PNA gives an AUC0-24h about 11% lower on the day of birth, under 1% different from day 6.",
      source_name = "PNAD"
    )
  )

  covariatesDataExcluded <- list(
    PAGE = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened as a sigmoid maturation function on CL (Table S1 models 8, 9 and 11), alone and combined with the PNA function; the PNA function alone fitted best (model 10) and PMA was not retained. Control stream column PMAW."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened as an effect on the endogenous IL-1Ra baseline (Table S1 model 12, dOFV -0.45); not retained."
    ),
    GA = list(
      description = "Gestational age at birth",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Carried in the dataset (control stream GA) and enters PMA, but has no direct effect in the final model. Cohort median 26.3 weeks (24.0-27.9)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 25L,
    n_studies = 1L,
    n_observations = 178L,
    age_range = "Postnatal age 0-20 days; gestational age at birth 24.0-27.9 weeks (median 26.3); postmenstrual age 24.0-30.9 weeks over the study",
    weight_range = "Baseline 0.54-1.30 kg (median 0.90); 0.53-1.60 kg (median 0.91) during therapy",
    sex_female_pct = 60,
    disease_state = "Extremely preterm neonates (born 24 weeks 0 days to 27 weeks 6 days gestation) receiving anakinra from birth to prevent inflammation-related complications of prematurity (bronchopulmonary dysplasia, its pulmonary hypertension, necrotising enterocolitis, white matter injury).",
    dose_range = "IV anakinra from birth for 21 days: stage 1, 1.0 mg/kg every 48 h (n = 9); stage 2, 0.8 mg/kg every 24 h (n = 16).",
    regions = "Australia (Monash Children's Hospital, Melbourne; 19 neonates) and New Zealand (Starship Children's Hospital, Auckland; 6 neonates).",
    baseline_il1ra = "Observed baseline endogenous IL-1Ra median 8754.9 pg/mL (294.5-87099.7).",
    notes = "Li 2026 Table 1 and Methods. Anakinra Pilot Phase I/IIa trial (NCT05280340). 178 samples (4-8 per neonate, median 8) taken pre-dose, at 12 and 21-24 h post-dose on days 1, 3, 7, 14 and 21, and 0.5-6 h post-dose on days 1 and 21. AlphaLISA immunoassay measuring endogenous plus exogenous IL-1Ra (LLOQ 1.3 pg/mL). NONMEM 7.5.1, FOCE-I; 1,000-replicate bootstrap (954 successful)."
  )

  ini({
    # Structural parameters standardised to a 70 kg adult with mature
    # postnatal-age clearance. The supplement control stream $THETA carries
    # the final estimates (they equal Table 2 and Eq. 6).
    lcl <- log(3.97); label("Clearance at full postnatal maturation for 70 kg (L/h)") # Table 2 'CL, L/h/70 kg' = 3.97 (RSE 14.3%); ctl THETA(1) = 3.97; Eq. 6
    lvc <- log(28.9); label("Volume of distribution for 70 kg (L)") # Table 2 'V, L/70 kg' = 28.9 (RSE 15.7%); ctl THETA(2) = 28.9; Eq. 7

    # Static endogenous IL-1Ra baseline added to the anakinra prediction
    # (ctl IPRED = F + BSL). 6640 pg/mL = 0.00664 mg/L (1 mg/L = 1e6 pg/mL).
    lbl_il1ra <- log(0.00664); label("Baseline endogenous IL-1Ra concentration (mg/L)") # Table 2 'BSL, pg/mL' = 6,640 (RSE 27.0%); ctl THETA(3) = 6640 pg/mL

    # Allometric exponents on total body mass, fixed a priori (Methods).
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL (unitless)") # Table 2 'WT on CL' = 0.75 (Fixed); ctl FWTCL = (WTKG/70)**0.75
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on V (unitless)") # Table 2 'WT on V' = 1 (Fixed); ctl FWTV = (WTKG/70)**1

    # Postnatal-age maturation of CL (Eq. 3 and 6):
    # Fmat = 1 - FCL * exp(-ln2 * PNA_days / PNAT50), written in the
    # <param>_pna_famp form as 1 + cl_pna_famp * exp(-cl_pna_kdes * PNA_days)
    # with cl_pna_famp = -FCL and cl_pna_kdes = ln2 / PNAT50.
    cl_pna_famp <- -0.554; label("Signed fractional difference of CL from its mature value at birth, -FCL (fraction)") # Table 2 'FCL' = 0.55 (RSE 13.4%); ctl THETA(4) PNAMAX = 0.554; Eq. 6 '1-0.554'
    lcl_pna_kdes <- log(log(2) / 1.8); label("Rate constant of postnatal CL maturation (1/day of PNA)") # Table 2 'PNAT50, day' = 1.80 (RSE 24.3%); ctl THETA(5) = 1.8; rate = ln(2) / 1.8

    # Assumed direct-effect IL-1 inhibition model used for the dose-finding
    # simulations (Eq. 5); not estimated from data in this study.
    limax <- fixed(log(1)); label("Maximum fractional IL-1 inhibition, Imax (fraction)") # Eq. 5 'Imax = 1' (assumed)
    lic50 <- fixed(log(0.04)); label("IL-1Ra concentration at 50% of maximum inhibition, IC50 (mg/L)") # Eq. 5 'IC50 = 0.04 mg/L' (literature mean)

    # IIV: Table 2 reports sqrt(omega); the control stream $OMEGA carries the
    # variances. IIV on V was fixed to 0 (unidentifiable) and is omitted.
    etalcl ~ 0.0798 # ctl $OMEGA PPVCL = 0.0798; Table 2 'IIV CL' = 0.28 = sqrt(0.0798)
    etalbl_il1ra ~ 0.898 # ctl $OMEGA PPBSL = 0.898; Table 2 'IIV BSL' = 0.95 = sqrt(0.898)

    # Proportional residual error on the total IL-1Ra concentration.
    propSd <- 0.765; label("Proportional residual error (fraction)") # ctl $SIGMA RUVPRO = 0.585 (variance) -> SD sqrt(0.585) = 0.765; Table 2 'Proportional error' = 0.77
  })
  model({
    # Postnatal-age maturation of CL over PNA in days.
    pna_days <- PNA * 30.4375
    cl_pna_kdes <- exp(lcl_pna_kdes)
    fpna_cl <- 1 + cl_pna_famp * exp(-cl_pna_kdes * pna_days)

    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * fpna_cl
    vc <- exp(lvc) * (WT / 70)^e_wt_vc
    bl_il1ra <- exp(lbl_il1ra + etalbl_il1ra)

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Cc is what the assay reports: anakinra plus the endogenous IL-1Ra
    # baseline. Subtract bl_il1ra for the anakinra-only concentration.
    Cc <- central / vc + bl_il1ra

    # Fractional IL-1 inhibition driven by total IL-1Ra (Eq. 5); at zero
    # anakinra the baseline alone gives 0.00664 / (0.04 + 0.00664) = 0.14,
    # matching the t = 0 value of Figure 3b and 3d.
    imax <- exp(limax)
    ic50 <- exp(lic50)
    inhib <- imax * Cc / (ic50 + Cc)

    Cc ~ prop(propSd)
  })
}
