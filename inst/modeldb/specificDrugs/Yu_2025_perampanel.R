Yu_2025_perampanel <- function() {
  description <- "One-compartment popPK model with first-order absorption of perampanel in Chinese pediatric epilepsy patients on routine therapeutic drug monitoring (Yu 2025). NONMEM 7.5.0 FOCE-I fit of 454 plasma concentrations from 151 patients aged 0.58-17.9 years, 120 (79.5%) of whom were younger than 12 years. The absorption rate constant was fixed to the adult value of 3.37 /h (Fujita 2023) because no pediatric estimate was available. Apparent clearance carries a linear-centralized age term ((AGE+10)/8.8)^1.31 plus multiplicative comedication effects: co-administered oxcarbazepine raises CL/F 1.51-fold and carbamazepine 1.88-fold (CYP3A4/5 induction), while sodium valproate lowers it to 0.745-fold (CYP3A4/5 inhibition). Apparent volume of distribution is proportional to the log of body weight. Inter-individual variability is exponential on clearance only; residual error is proportional."
  reference <- "Yu L, Mao F, Chen S, Liu J, Xiao J, Chen M, Luo H, Yu Z, Dai H. Development and Validation of a Population Pharmacokinetics Model of Perampanel for Pediatric Epilepsy Patients for Optimized Dosing. Drug Des Devel Ther. 2025;19:3119-3128. doi:10.2147/DDDT.S499085"
  vignette <- "Yu_2025_perampanel"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/L"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot   = list(analyte = "perampanel", units = "mg", specimen = "administration site", verified = FALSE),
    central = list(analyte = "perampanel", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    AGE = list(
      description        = "Age at the time of the therapeutic-drug-monitoring sample",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying in principle; Yu 2025 Table 1 reports baseline age. Enters apparent clearance as the linear-centralized power term ((AGE + 10)/8.8)^e_age_cl (Yu 2025 Equation 1). The normalisation constants 10 and 8.8 are printed in Equation 1 and are NOT the sample median (median age 9.00 years); they are used exactly as published. Note that ((AGE + 10)/8.8) = 1 at AGE = -1.2 years, so lcl is an extrapolated coefficient rather than a typical value at any observed age -- see the lcl label. Studied range 0.58-17.9 years (Table 1).",
      source_name        = "Age"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters the apparent volume of distribution as V = lvc * log10(WT) (Yu 2025 Equation 2, V(L) = 227 * LGBW, where LGBW is described in the text as 'the log value of body weight'). The published equation applies no normalisation, so lvc is a coefficient in L per base-10-log unit of body weight in kg, not a volume at a reference weight. Studied range 9.00-89.0 kg, median 28.1 kg (Table 1). See the vignette Errata for the log-base reading.",
      source_name        = "BW / LGBW"
    ),
    CONMED_OXC = list(
      description        = "Concomitant oxcarbazepine",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant oxcarbazepine)",
      notes              = "Yu 2025 Equation 1 defines the comedication covariates as taking the value 0 when absent and 1 when the drug is administered simultaneously with perampanel. 26 of 151 patients (17.2%) received oxcarbazepine (Table 1). Multiplicative effect on apparent clearance (cl *= e_oxc_cl^CONMED_OXC) with e_oxc_cl = 1.51 (Table 2 final model); oxcarbazepine is a CYP3A4 inducer and perampanel is cleared predominantly by CYP3A4/5.",
      source_name        = "OXC"
    ),
    CONMED_VPA = list(
      description        = "Concomitant sodium valproate",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant sodium valproate)",
      notes              = "Yu 2025 Equation 1; 0 when absent, 1 when co-administered. 53 of 151 patients (35.1%) received sodium valproate (Table 1); a further single patient (0.66%) received magnesium valproate, which the covariate definition in Equation 1 does not cover. Multiplicative effect on apparent clearance (cl *= e_vpa_cl^CONMED_VPA) with e_vpa_cl = 0.745 (Table 2 final model), i.e. valproate lowers CL/F by about 25% via CYP3A4/5 inhibition.",
      source_name        = "VPA"
    ),
    CONMED_CBZ = list(
      description        = "Concomitant carbamazepine",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant carbamazepine)",
      notes              = "Yu 2025 Equation 1; 0 when absent, 1 when co-administered. 8 of 151 patients (5.30%) received carbamazepine (Table 1). Multiplicative effect on apparent clearance (cl *= e_cbz_cl^CONMED_CBZ) with e_cbz_cl = 1.88 (Table 2 final model); carbamazepine is a strong CYP3A4 inducer.",
      source_name        = "CBZ"
    )
  )

  # Covariates that Yu 2025 screened but did NOT retain in the final model.
  # Documented here for provenance only; these names are not referenced in
  # model() and checkModelConventions() does not require them to be.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened in the covariate analysis (Yu 2025 Methods, 'PPK Modeling and Validation') but not retained in the final model."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened but not retained. Correlated with body weight; Yu 2025 dropped one of any pair of covariates with correlation coefficient > 0.3."
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = "Screened but not retained (correlated with body weight)."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened but not retained (correlated with body weight)."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened but not retained."
    ),
    TPROT = list(
      description = "Total protein",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened but not retained."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Entered on clearance at forward inclusion (Table S1 step 7, dOFV -5.03) but removed at backward elimination (step 8, dOFV +5.03 < 10.83)."
    ),
    CONMED_PB = list(
      description = "Concomitant phenobarbital",
      units       = "(binary)",
      type        = "binary",
      notes       = "Entered on clearance at forward inclusion (Table S1 step 6, dOFV -5.89) but removed at backward elimination (step 9, dOFV +5.89 < 10.83). Only 2 of 151 patients (1.32%) received phenobarbital."
    ),
    CONMED_PHT = list(
      description = "Concomitant phenytoin",
      units       = "(binary)",
      type        = "binary",
      notes       = "Listed among the enzyme-inducing antiseizure medications screened (Yu 2025 Methods) but not retained; no phenytoin use is tabulated in Table 1."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 151L,
    n_studies      = 1L,
    age_range      = "0.58-17.9 years",
    age_median     = "9.00 years (IQR 6.34-11.88)",
    weight_range   = "9.00-89.0 kg",
    weight_median  = "28.1 kg (IQR 21.2-41.0)",
    sex_female_pct = 36.4,
    disease_state  = "Pediatric patients with epilepsy (focal 51.0%, generalized 46.4%, focal evolving to generalized 2.65%) receiving oral perampanel; no patient had severe hepatic or renal impairment",
    dose_range     = "1.00-8.00 mg total daily oral perampanel; median daily dose 2.00 mg (IQR 2.00-4.00). Monte Carlo dosing simulations spanned 2-12 mg/d.",
    regions        = "China (Second Affiliated Hospital, Zhejiang University School of Medicine, Hangzhou)",
    notes          = "Retrospective single-centre study of routine TDM records collected February 2021 - September 2023 (Yu 2025 Methods). 454 plasma concentrations from 151 patients; 120 patients (79.5%) were younger than 12 years and 31 were aged 12-18 years. TDM was performed about 3 weeks after perampanel initiation, with samples drawn in the morning approximately 12 h after the previous dose, so essentially all observations are near-trough steady-state concentrations. Assay: validated HPLC with a range of 15-1500 ng/mL; observed perampanel concentrations 30.0-1082 ng/mL, median 242 ng/mL (IQR 144-340). Median height 134 cm, BMI 16.6 kg/m^2, BSA 1.06 m^2, CrCl 149 mL/min. Comedications: sodium valproate 35.1%, levetiracetam 27.8%, topiramate 18.5%, oxcarbazepine 17.2%, lamotrigine 8.61%, lacosamide 6.62%, clonazepam 5.96%, carbamazepine 5.30%, clobazam 4.64%, zonisamide 3.97%, vigabatrin 1.99%, phenobarbital 1.32%, magnesium valproate 0.66%, nitrazepam 5.83%/9.68% by age stratum."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural PK parameters. Yu 2025 Table 2 "Final Model / Estimate"
    # column, and Equations (1) and (2) on page 4.
    #
    # A one-compartment model with first-order absorption and first-order
    # elimination best described the data (Yu 2025 "PPK Model
    # Development"). All disposition parameters are apparent (CL/F, V/F)
    # because dosing was oral throughout; no bioavailability term was
    # estimated, so F is folded into cl and vc.
    # ------------------------------------------------------------------

    # KA was fixed, not estimated: Table 2 reports "3.37 FIXED" with no SE
    # or RSE, and footnote a states that KA was taken from Fujita Y et al.
    # Ther Drug Monit. 2023;45(5):653-659. The Discussion notes "There is
    # no absorption constant available for pediatric patients, and an
    # adult constant was applied in modeling."
    lka <- fixed(log(3.37)); label("Perampanel first-order absorption rate constant ka (1/h), from the adult value of Fujita 2023")  # Yu 2025 Table 2 KA = 3.37 FIXED

    # 0.177 is the CL/F coefficient in Equation (1), i.e. the value CL/F
    # takes when ((AGE + 10)/8.8)^1.31 = 1 and no comedication is present.
    # That happens at AGE = -1.2 years, so 0.177 L/h is an extrapolated
    # intercept rather than a typical value at any studied age. At the
    # sample median age of 9.00 years the model gives
    # 0.177 * ((9 + 10)/8.8)^1.31 = 0.485 L/h.
    lcl <- log(0.177); label("Perampanel apparent clearance CL/F coefficient (L/h) at (AGE + 10)/8.8 = 1 with no comedication")  # Yu 2025 Table 2 CL/F = 0.177; Equation (1)

    # 227 is the V coefficient in Equation (2), V(L) = 227 * LGBW. The
    # published equation applies no normalisation to LGBW, so 227 has
    # units of L per log10 unit of body weight in kg -- it is not a volume
    # at a reference weight. At the sample median weight of 28.1 kg the
    # model gives 227 * log10(28.1) = 329 L.
    lvc <- log(227); label("Perampanel apparent volume of distribution V/F coefficient (L per log10 unit of body weight in kg)")  # Yu 2025 Table 2 V = 227; Equation (2)

    # ------------------------------------------------------------------
    # Covariate effects on apparent clearance. Yu 2025 Equation (1):
    #   CL/F = 0.177 * ((Age + 10)/8.8)^1.31 * 1.51^OXC * 0.745^VPA * 1.88^CBZ
    # The three comedication covariates take the value 0 when absent and 1
    # when co-administered, so each multiplier applies once when present.
    # ------------------------------------------------------------------
    e_age_cl <- 1.31;  label("Power-function exponent of linear-centralized age (AGE + 10)/8.8 on CL/F (unitless)")  # Yu 2025 Table 2 AGE = 1.31; Equation (1)
    e_oxc_cl <- 1.51;  label("Multiplicative effect of concomitant oxcarbazepine on CL/F (unitless)")                # Yu 2025 Table 2 Oxcarbazepine = 1.51; Equation (1)
    e_vpa_cl <- 0.745; label("Multiplicative effect of concomitant sodium valproate on CL/F (unitless)")             # Yu 2025 Table 2 Sodium valproate = 0.745; Equation (1)
    e_cbz_cl <- 1.88;  label("Multiplicative effect of concomitant carbamazepine on CL/F (unitless)")                # Yu 2025 Table 2 Carbamazepine = 1.88; Equation (1)

    # ------------------------------------------------------------------
    # Inter-individual variability. Yu 2025 estimated IIV on clearance
    # only; the Discussion states "it is difficult to estimate the
    # interindividual viability of V", so no eta is carried on the volume.
    # Table 2 reports the clearance term as "omega CL = 0.0963" and the
    # table footnote defines it as "interindividual variance for CL", so
    # 0.0963 is entered directly as the log-normal variance
    # (equivalently CV = sqrt(exp(0.0963) - 1) = 31.8%).
    # ------------------------------------------------------------------
    etalcl ~ 0.0963  # Yu 2025 Table 2 omega CL = 0.0963, defined as a variance in the table footnote

    # ------------------------------------------------------------------
    # Residual error. Table 2 reports "sigma = 0.130" with the footnote
    # "sigma, residual variability for proportional error". As with omega,
    # NONMEM reports $SIGMA on the variance scale, so the proportional
    # residual standard deviation is sqrt(0.130) = 0.3606 (36.1%).
    # ------------------------------------------------------------------
    propSd <- 0.360555; label("Proportional residual error on perampanel plasma concentration (fraction)")  # Yu 2025 Table 2 sigma = 0.130 (variance); sqrt(0.130) = 0.360555
  })

  model({
    # 1. Individual PK parameters
    #    CL/F_i = 0.177 * ((AGE + 10)/8.8)^1.31
    #                   * 1.51^CONMED_OXC * 0.745^CONMED_VPA * 1.88^CONMED_CBZ
    #                   * exp(eta_lcl)                              Equation (1)
    #    V/F_i  = 227 * log10(WT)                                   Equation (2)
    #    ka     = 3.37 /h, fixed                                    Table 2
    ka <- exp(lka)
    cl <-
      exp(lcl + etalcl) *
      ((AGE + 10) / 8.8)^e_age_cl *
      e_oxc_cl^CONMED_OXC *
      e_vpa_cl^CONMED_VPA *
      e_cbz_cl^CONMED_CBZ
    vc <- exp(lvc) * log10(WT)

    # 2. Micro-constant
    kel <- cl / vc

    # 3. One-compartment ODEs with first-order oral absorption
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 4. Observation and error
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
