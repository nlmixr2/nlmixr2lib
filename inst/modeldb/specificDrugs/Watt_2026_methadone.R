Watt_2026_methadone <- function() {
  description <- "One-compartment population PK model of methadone in hospitalized children from birth to 21 years given intravenous or enteral methadone per standard of care for pain or iatrogenic opiate withdrawal (Watt 2026). First-order enteral absorption with ka fixed at 2.72 1/h, enteral bioavailability 0.64, allometric total-body-weight scaling on CL (exponent 0.75) and linear total-body-weight scaling on V (exponent 1) with both exponents fixed and NO weight normalization (the published equations scale raw kg, so CL = 0.325 * WT^0.75 = 7.87 L/h and V = 5.22 * WT = 365 L at 70 kg), correlated interindividual variability on CL and V, and proportional residual error. Body weight is the only covariate: the univariate screen flagged total bilirubin, serum creatinine and obesity on CL and postnatal age and total bilirubin on V, none survived backward elimination, fat-free mass performed no better than total body weight, and postnatal- / postmenstrual-age maturation functions did not improve the fit, so the base model is also the final irreducible model."
  reference <- "Watt KM, Thompson EJ, Lam L, Zimmerman K, Hornik CP, Atz AM, Fernandez A, Hupp SR, Bhatt-Mehta V, Benjamin DK Jr, Anand R, Cohen-Wolkowiez M, Gonzalez D, Smith PB, Capparelli EV; Best Pharmaceuticals for Children Act - Pediatric Trials Network Steering Committee. Population Pharmacokinetics to Support Intravenous and Enteral Methadone Dosing in Children. J Clin Pharmacol. 2026;66(1):e70143. doi:10.1002/jcph.70143"
  vignette <- "Watt_2026_methadone"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/L"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Watt 2026 Methods (plasma methadone
  # assayed by LC-MS/MS; enteral doses given orally or via feeding tube).
  compartmentData <- list(
    depot   = list(analyte = "methadone", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "methadone", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "The only covariate retained in the final model. Watt 2026 Equations 4-5 scale RAW kilograms with no reference-weight normalization: V (L) = 5.22 * WT and CL (L/h) = 0.325 * WT^0.75, which the Results restate as 365 L/70 kg and 7.87 L/h/70 kg. Studied range 0.72-159 kg (median 13.0 kg, Table 1). Fat-free mass was tested in place of total body weight and gave similar clearance and interindividual variability (83% vs 88% on CL, 70% vs 74% on V), so total body weight was kept for dosing practicality.",
      source_name        = "TBW"
    )
  )

  # Covariates Watt 2026 screened but did NOT retain in the final model
  # (Covariate Analysis: none survived backward elimination at P < ~.005).
  # Documentation only -- checkModelConventions() does not require these to
  # appear in model().
  covariatesDataExcluded <- list(
    FFM = list(
      description = "Fat-free mass",
      units       = "kg",
      type        = "continuous",
      notes       = "Tested in place of total body weight as the size descriptor. Estimated by the sex-specific Al-Sallami equations for children >= 2 years (Watt 2026 Equations 1-2): FFM_female = [1.11 + (1 - 1.11) / (1 + (age / 7.1)^-1.1)] * 9270 * TBW / (8780 + 244 * BMI) and FFM_male = [0.88 + (1 - 0.88) / (1 + (age / 13.4)^-12.7)] * 9270 * TBW / (6680 + 216 * BMI); 90% of TBW when height was missing; TBW used unchanged for children < 2 years, for whom the equations are not validated. Not selected -- model performance was similar and dosing by FFM is impractical."
    ),
    PNA = list(
      description = "Postnatal age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened on V (significant univariately) and as the age driver of a sigmoidal CL maturation function Fage = age^HILL / (TM50^HILL + age^HILL) with TM50 in days (Watt 2026 Equation 3). Neither survived: the maturation function did not improve the fit and the ETA-versus-PNA plots (Figure S1) show no trend. Range 0-19.02 years."
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units       = "weeks",
      type        = "continuous",
      notes       = "Tested as the alternative age driver of the same sigmoidal CL maturation function (Watt 2026 Equation 3, TM50 in weeks for PMA). Did not improve the fit."
    ),
    HT = list(
      description = "Body height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened on CL and V; not retained. Range 30-189 cm. Also an input to the BMI used by the fat-free-mass equations."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL and V; not retained. 50/99 male (Table 1). Sex enters the fat-free-mass equations, which were themselves not selected."
    ),
    RACE_WHITE = list(
      description = "White race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL and V; not retained. 71/99 (72%) White."
    ),
    RACE_BLACK = list(
      description = "Black or African American race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL and V; not retained. 20/99 (20%) Black or African American."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL and V; not retained. 1/99 (1%) Asian."
    ),
    RACE_MULTI = list(
      description = "Multiracial indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL and V; not retained. 6/99 (6%) multiple races."
    ),
    RACE_HISPANIC = list(
      description = "Hispanic or Latino ethnicity indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL and V as 'ethnicity'; not retained. 13/99 (19% of those reporting) Hispanic or Latino."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Screened on CL and V; not retained. Collected in POPS only; median 3.3 g/dL (range 1.9-4.5)."
    ),
    AAG = list(
      description = "Alpha-1-acid glycoprotein",
      units       = "mg/mL",
      type        = "continuous",
      notes       = "Screened on CL and V; not retained. Methadone is highly bound to AAG. Collected in MTH01 only; median 2.34 mg/mL (range 0.66-6.41)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened on CL and V; not retained. Collected in POPS only; median 52 U/L (range 11-280)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened on CL and V; not retained. Collected in POPS only; median 48 U/L (range 16-178)."
    ),
    HCT = list(
      description = "Hematocrit",
      units       = "%",
      type        = "continuous",
      notes       = "Screened on CL and V; not retained. Collected in MTH01 only; median 30.9% (range 30.9-39.5)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened on CL and V; significant on CL in the univariate screen but eliminated in backward selection. Collected in POPS only; median 0.3 mg/dL (range 0.1-1.2)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened on CL and V; significant on BOTH CL and V in the univariate screen but eliminated in backward selection. Collected in POPS only; median 0.8 mg/dL (range 0.2-4.9). Reported in US units here, matching the source; the canonical register unit is umol/L."
    ),
    CONMED_AZOLE = list(
      description = "Concomitant azole antifungal therapy (itraconazole, ketoconazole, voriconazole or fluconazole)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL, V AND bioavailability; not retained. Only 5/99 (5%) exposed, so the screen had little power against a real CYP3A4-inhibition interaction."
    ),
    MEAL_PREDOSE_2H = list(
      description = "Food intake within 2 h before an enteral methadone dose",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on bioavailability only; not retained. Only 11/99 participants received both IV and enteral doses, so F covariate effects were poorly estimable (Discussion, limitations)."
    ),
    ROUTE_IV = list(
      description = "Intravenous administration indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a covariate on bioavailability and not retained as such. Route is handled STRUCTURALLY in this model rather than as a covariate: IV doses are given into `central` and enteral doses into `depot`, where f(depot) = 0.64 applies. 232/1798 (13%) of recorded doses were IV."
    ),
    OBESE_BMI95 = list(
      description = "Pediatric obesity indicator, body mass index >= 95th percentile for age and sex",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL and V; significant on CL in the univariate screen but eliminated in backward selection. Per CDC guidance, children < 2 years (44/99) were not classified; among the remaining 55 children 22 (22% of the cohort) were obese. NO canonical register entry exists for a pediatric BMI-percentile obesity flag -- none was proposed because the covariate was screened and dropped, so this name is descriptive only and is not a ratified canonical. The register's DIS_OBESE_MORBID is a different concept (morbidly obese cohort indicator) and BMIZ is the continuous z-score."
    ),
    CYP3A4_STAR1B = list(
      description = "CYP3A4*1B genotype (*1/*1 or *1/*1B versus *1B/*1B)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened on CL and V; not retained. Genotyped in MTH01 only: 19 *1/*1 or *1/*1B, 4 *1B/*1B, 76 unknown (Table 1) -- 77% missing, so the screen had almost no power. NO canonical register entry exists for the CYP3A4*1B (rs2740574) star-allele encoding; the registered SNP_CYP3A4_* entries are different variants and the registered CYP3A4 entry is a metabolic-activity score, not a genotype. This name is descriptive only and is not a ratified canonical."
    ),
    CYP2B6_STAR6 = list(
      description = "CYP2B6*6 genotype (*1/*1 or *1/*6 versus *6/*6)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened on CL and V; not retained. Genotyped in MTH01 only: 20 *1/*1 or *1/*6, 3 *6/*6, 76 unknown (Table 1). CYP2B6 is a principal methadone-metabolizing enzyme, so this is a mechanistically plausible covariate that these data could not test. NO canonical register entry exists for the *6 star-allele encoding; SNP_CYP2B6_RS3745274_T_COUNT is the 516G>T tag SNP rather than the star-allele diplotype, and equating them is a mapping claim this extraction does not make. This name is descriptive only and is not a ratified canonical."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 99L,
    n_studies      = 2L,
    age_range      = "0-19.02 years (postnatal age); eligibility birth to < 21 years",
    age_median     = "2.29 years",
    weight_range   = "0.72-159.0 kg",
    weight_median  = "13.0 kg",
    sex_female_pct = 50.5,
    race_ethnicity = c(White = 72.0, Black = 20.0, Asian = 1.0, Multiple = 6.0, UnknownOrNotReported = 1.0),
    disease_state  = "Hospitalized children prescribed methadone per standard of care for pain or iatrogenic opiate withdrawal",
    dose_range     = "IV 0.11 (0.01-0.39) mg/kg and enteral 0.10 (0.01-0.61) mg/kg per dose, median (range); 1798 recorded doses of which 232 (13%) were IV",
    regions        = "USA (23 enrolling children's hospitals)",
    notes          = "Pooled from two prospective, multi-center, open-label PK studies (Watt 2026 Table 1): MTH01 (NCT01945736, 5 sites, n = 26, multiple-dose enteral methadone for iatrogenic withdrawal in children >= 90 days to < 18 years, scheduled sampling windows per Table S1) and POPS (NCT01431326, 18 sites, n = 73, opportunistic standard-of-care sampling in children < 21 years, max 10 samples). 263 of 273 collected samples were analyzable (6 insufficient quantity, 4 dropped for high weighted residuals); none below the 0.1-100 ng/mL assay validation range. Median (range) methadone concentration 42.2 (0.9-729.2) ng/mL sampled a median 4 h (0-31.8) after the last dose. 10 participants received IV doses only, 78 enteral only, and 11 both. Laboratory covariates were split by protocol: alpha-1-acid glycoprotein, hematocrit and CYP genotypes in MTH01 only; albumin, ALT, AST, serum creatinine and total bilirubin in POPS only."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Watt 2026 Table 2 point estimates, with
    # the size model given by Equations 4-7. The published equations
    # scale RAW kilograms (there is no reference-weight divisor), so
    # exp(lcl) and exp(lvc) are the per-kg^0.75 and per-kg coefficients
    # and NOT clearances/volumes at a reference weight.
    # ------------------------------------------------------------------
    lcl     <- log(0.325); label("Clearance coefficient, CL = exp(lcl) * WT^0.75 (L/h/kg^0.75; 7.87 L/h at 70 kg)")  # Watt 2026 Table 2 (CL 0.325 L/h/kg, %RSE 21.6, bootstrap 0.171-0.516) with Equation 5
    lvc     <- log(5.22);  label("Central volume coefficient, V = exp(lvc) * WT (L/kg; 365 L at 70 kg)")              # Watt 2026 Table 2 (V 5.22 L/kg, %RSE 21.8, bootstrap 2.70-8.01) with Equation 4
    lka     <- fixed(log(2.72)); label("First-order enteral absorption rate constant (ka, 1/h)")                      # Watt 2026 Table 2 (KA 2.72 FIX) and Equation 6; fixed at the base-model initial estimate because the data could not characterize absorption, and within the published 0.33-5.93 1/h range
    lfdepot <- log(0.64);  label("Enteral bioavailability (F, fraction)")                                             # Watt 2026 Table 2 (F 0.64, %RSE 21.8, bootstrap 0.32-1.00) and Equation 7

    # Size exponents. Watt 2026 Methods: total body weight "was included
    # as a covariate for CL (allometric scaling - WT^0.75) and V (linear
    # scaling - WT^1.0) in the base model"; both are structural
    # assumptions carried into the final model with no RSE or CI, so
    # both are fixed.
    e_wt_cl <- fixed(0.75); label("Total-body-weight exponent on CL (unitless)")  # Watt 2026 Methods, Population Pharmacokinetic Analysis; Equation 5
    e_wt_vc <- fixed(1.0);  label("Total-body-weight exponent on V (unitless)")   # Watt 2026 Methods, Population Pharmacokinetic Analysis; Equation 4

    # ------------------------------------------------------------------
    # Correlated interindividual variability on CL and V.
    #
    # Watt 2026 Table 2 prints the two diagonals as 'interindividual
    # variability (CV%)' -- 88.1 for CL and 74.0 for V -- and the
    # off-diagonal as a bare 'CL ~ V interindividual variability
    # covariance' of 0.296. The off-diagonal is therefore on the raw
    # OMEGA scale, so the diagonals must be too, or the printed
    # covariance would be incommensurable with them and no correlation
    # could be recovered from the table. That fixes the reading as
    # omega SD = CV% / 100 (variances 0.881^2 and 0.740^2, correlation
    # 0.454) rather than the exact log-normal omega^2 = log(1 + CV^2)
    # (which would give correlation 0.591). The paper's own dosing
    # simulations corroborate it: reproducing Table S2's toxicity
    # attainment across all 54 dose-by-route-by-dose-number cells,
    # with the virtual-population weight scale calibrated to the
    # published median AUC0-tau, matches to within 2 percentage points
    # under this reading and only within 4 under the alternative. See
    # the vignette's Assumptions and deviations section.
    # ------------------------------------------------------------------
    etalcl + etalvc ~ c(0.776161, 0.296, 0.547600)  # Watt 2026 Table 2: CL IIV 88.1 CV% -> 0.881^2; 'CL ~ V interindividual variability covariance' 0.296; V IIV 74.0 CV% -> 0.740^2

    # ------------------------------------------------------------------
    # Residual error -- proportional only (Results: "A one-compartment
    # PK model with proportional residual error described the methadone
    # concentration versus time data well").
    # ------------------------------------------------------------------
    propSd  <- 0.240; label("Proportional residual error (fraction)")  # Watt 2026 Table 2 (residual proportional error 24.0 CV%, %RSE 11.6, bootstrap 18.9-29.0)
  })

  model({
    # 1. Individual parameters. Raw body weight in kg, no reference-weight
    # normalization (Watt 2026 Equations 4-5).
    cl <- exp(lcl + etalcl) * WT^e_wt_cl
    vc <- exp(lvc + etalvc) * WT^e_wt_vc
    ka <- exp(lka)
    fdepot <- exp(lfdepot)

    # 2. Micro-constant
    kel <- cl / vc

    # 3. One-compartment system with a first-order enteral depot. Enteral
    # doses enter `depot`; IV doses bypass it by dosing `central` directly
    # from the event table.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 4. Enteral bioavailability applies only to the depot input, so IV
    # doses are complete by construction (F_IV = 1).
    f(depot) <- fdepot

    # 5. Plasma methadone concentration in mg/L (dose mg, volume L).
    # Multiply by 1000 to compare with the paper's ng/mL targets.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
