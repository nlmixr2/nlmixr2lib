Moes_2016_tacrolimus <- function() {
  description <- "Two-compartment population pharmacokinetic model for oral once-daily tacrolimus (Advagraf) in stable adult liver transplant recipients (Moes 2016), with first-order elimination from the central compartment and a delayed first-order absorption phase described by three sequential transit compartments sharing the absorption rate constant ka, a fixed oral bioavailability F = 0.23, a categorical donor + recipient CYP3A5*3 combination effect on apparent oral clearance (reference both nonexpressers; donor nonexpresser + recipient *1 carrier +33%; donor *1 carrier + recipient nonexpresser +33%; both *1 carriers +71%), independent log-normal IIV on CL, Vc, and ka, and proportional residual error on whole-blood concentration."
  reference   <- "Moes DJAR, van der Bent SAS, Swen JJ, van der Straaten T, Inderson A, Olofsen E, Verspaget HW, Guchelaar HJ, den Hartigh J, van Hoek B. Population pharmacokinetics and pharmacogenetics of once daily tacrolimus formulation in stable liver transplant recipients. Eur J Clin Pharmacol. 2016;72(2):163-174. doi:10.1007/s00228-015-1963-3"
  vignette    <- "Moes_2016_tacrolimus"
  units       <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    CYP3A5_EXPR = list(
      description        = "Recipient CYP3A5 expresser indicator: 1 if the liver-transplant recipient carries at least one functional CYP3A5*1 allele (genotype *1/*1 or *1/*3 at rs776746), 0 if homozygous *3/*3.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (recipient CYP3A5 *3/*3 nonexpresser)",
      notes              = "Time-fixed germline genotype derived from rs776746. In Moes 2016 the model-development cohort (n = 49) had recipient genotype A/A = 36, G/A = 10, G/G = 3 (Table 2) -- i.e., 13 of 49 (26.5%) recipients carried at least one *1 allele and are encoded as CYP3A5_EXPR = 1. The recipient genotype represents intestinal CYP3A5 expression and contributes to first-pass extraction independently of the donor liver genotype. Combined with CYP3A5_EXPR_DONOR to reconstruct the Moes 2016 four-level combination categories C1-C4 inside model().",
      source_name        = "Recipient CYP3A5*3"
    ),
    CYP3A5_EXPR_DONOR = list(
      description        = "Donor CYP3A5 expresser indicator: 1 if the transplanted liver was donated by a donor carrying at least one functional CYP3A5*1 allele (genotype *1/*1 or *1/*3 at rs776746), 0 if the donor was homozygous *3/*3.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (donor CYP3A5 *3/*3 nonexpresser graft)",
      notes              = "Time-fixed donor germline genotype (recovered from donor spleen or liver biopsy at transplantation). In Moes 2016 the model-development cohort (n = 49) had donor genotype A/A = 40, G/A = 9, G/G = 0 (Table 2) -- i.e., 9 of 49 (18.4%) donors carried at least one *1 allele and are encoded as CYP3A5_EXPR_DONOR = 1. The donor genotype represents the engrafted-liver hepatic CYP3A5 expression and is biologically distinct from the recipient's intestinal CYP3A5 contribution. Combined with CYP3A5_EXPR to reconstruct the Moes 2016 four-level combination categories C1-C4 inside model().",
      source_name        = "Donor CYP3A5*3"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 49L,
    n_observations   = 282L,
    age_range        = "29-69 years",
    age_median       = "55 years",
    age_mean_sd      = "54 +/- 11 years",
    weight_range     = "50-131 kg",
    weight_median    = "84 kg",
    weight_mean_sd   = "84 +/- 18 kg (Table 1 popPK column; the Results section narrative reports a different mean of 77.5 +/- 11.8 kg with range 50-121 kg -- see Errata in the validation vignette)",
    sex_female_pct   = 36.7,
    race_ethnicity   = "Caucasian 92% (45 of 49); remaining 8% not separately stratified.",
    disease_state    = "Stable adult liver transplant recipients (bilirubin and albumin within reference range, stable graft function for at least 3 months) converted from twice-daily tacrolimus (Prograf) to once-daily tacrolimus (Advagraf) with at least 2 weeks of stable Advagraf dosing prior to PK sampling.",
    dose_range       = "Once-daily oral Advagraf 0.5-14 mg/day; mean 3.6 +/- 2.2 mg/day, median 3 mg/day (Table 1 popPK column). Each PK profile was a single AUC0-6h profile after at least 2 weeks of dosing on the new once-daily formulation.",
    regions          = "The Netherlands (Leiden University Medical Center).",
    co_medication    = "Primary diagnoses included alcoholic liver disease (24.5%), primary sclerosing cholangitis (18%), hepatitis C (8%), cystic liver disease (8%), nonalcoholic steatohepatitis (6%), primary biliary cirrhosis (4%), cryptogenic liver disease (4%), and other (24%). Concomitant medications including prednisolone (low doses <=10 mg) were tested as covariates and not retained.",
    sampling_design  = "Single AUC0-6h PK profile per subject (predose plus 1, 2, 3, 4, and 6 hours postdose); 282 whole-blood tacrolimus concentrations across 49 subjects, all collected after the patient had been stable on once-daily Advagraf for at least 2 weeks.",
    cyp3a5_distribution_recipient = "A/A (rs776746) n = 36 (73.5%), G/A n = 10 (20.4%), G/G n = 3 (6.1%); 13 of 49 (26.5%) CYP3A5*1 carriers (Table 2).",
    cyp3a5_distribution_donor     = "A/A n = 40 (81.6%), G/A n = 9 (18.4%), G/G n = 0; 9 of 49 (18.4%) CYP3A5*1 carriers (Table 2).",
    cyp3a5_combination_distribution = "C1 (both nonexpressers) n = 32; C2 (recipient *1 carrier + donor nonexpresser) n = 8; C3 (recipient nonexpresser + donor *1 carrier) n = 4; C4 (both *1 carriers) n = 5 (Table 2).",
    notes            = "Baseline demographics from Table 1 popPK column. The model-development dataset was the subset of 49 of 66 enrolled subjects for whom both recipient and donor DNA were available. The remaining 17 subjects (donor DNA not available) contributed only to the limited-sampling-strategy development, not to the population PK / pharmacogenetic analysis whose parameters this file encodes."
  )

  ini({
    # --- Structural PK parameters ---
    # Moes 2016 Table 4 "Final model" column (mean values). With oral
    # bioavailability F fixed at 0.23 in NONMEM and applied to the dose, the
    # reported CL, Vc, Q, and Vp are actual (non-apparent) values; the
    # apparent oral parameters are obtained by dividing by F = 0.23.
    # Time in hours; rates in 1/h; clearances in L/h; volumes in L.
    lka <- log(3.76)  ; label("First-order absorption rate constant ka driving every transition through the three-transit-compartment chain (1/h)")  # Moes 2016 Table 4 Final model Ka = 3.76 1/h (RSE 10%)
    lcl <- log(4.21)  ; label("Plasma clearance CL at the reference CYP3A5 combination C1 (donor + recipient both *3/*3 nonexpressers) (L/h)")        # Moes 2016 Table 4 Final model CL = 4.21 L/h (RSE 8%)
    lvc <- log(88.3)  ; label("Central volume of distribution Vc (L)")                                                                                  # Moes 2016 Table 4 Final model Vc = 88.3 L (RSE 12%)
    lq  <- log(14.0)  ; label("Inter-compartmental clearance Q (L/h)")                                                                                  # Moes 2016 Table 4 Final model Q = 14 L/h (RSE 22%)
    lvp <- log(145.0) ; label("Peripheral volume of distribution Vp (L)")                                                                               # Moes 2016 Table 4 Final model Vp = 145 L (RSE 41%)

    # --- Fixed oral bioavailability ---
    # Moes 2016 Methods Base model: "The value for bioavailability was fixed
    # to 0.23 which was based on literature [22]." NONMEM applied F to the
    # dose; the table-4 CL, Vc, Q, and Vp are therefore reported as actual
    # (not apparent) parameters and the model file applies the same F = 0.23
    # via f(depot).
    lfdepot <- fixed(log(0.23)) ; label("Fixed oral bioavailability F (unitless)")  # Moes 2016 Methods Base model "F fixed to 0.23" and Table 4 row "F (fixed) = 0.23"

    # --- Covariate effects: CYP3A5*3 donor + recipient combination on CL ---
    # Moes 2016 Methods covariate-effect equation:
    #   TV(CL) = THETA_CL * (1 + theta_cov)
    # i.e., theta_cov is the fractional deviation of the comparator
    # combination's typical CL from the reference combination C1. Moes 2016
    # Table 4 reports the final-model values C1 = 0% (reference), C2 = 33%,
    # C3 = 33%, C4 = 71%. The four levels are reconstructed inside model()
    # from the two binary inputs CYP3A5_EXPR (recipient) and
    # CYP3A5_EXPR_DONOR (donor):
    #   C1 = both nonexpressers      (CYP3A5_EXPR = 0, DONOR = 0)
    #   C2 = recipient + donor nonex (CYP3A5_EXPR = 1, DONOR = 0)
    #   C3 = recipient nonex + donor (CYP3A5_EXPR = 0, DONOR = 1)
    #   C4 = both expressers         (CYP3A5_EXPR = 1, DONOR = 1)
    e_cyp3a5_c2_cl <- 0.33 ; label("Recipient *1 carrier + donor nonexpresser fractional CL deviation vs C1 (unitless)") # Moes 2016 Table 4 Final model C2 = +33%
    e_cyp3a5_c3_cl <- 0.33 ; label("Recipient nonexpresser + donor *1 carrier fractional CL deviation vs C1 (unitless)") # Moes 2016 Table 4 Final model C3 = +33%
    e_cyp3a5_c4_cl <- 0.71 ; label("Recipient + donor both *1 carriers fractional CL deviation vs C1 (unitless)")         # Moes 2016 Table 4 Final model C4 = +71%

    # --- Inter-individual variability ---
    # Moes 2016 Table 4 reports BSV as percent CV on CL, Vc, and Ka in the
    # final model. Conversion to log-scale variance via
    # omega^2 = log(CV^2 + 1):
    #   CL  CV 42.8% -> log(0.428^2 + 1) = 0.16839
    #   Vc  CV 86.3% -> log(0.863^2 + 1) = 0.55657
    #   Ka  CV 65.9% -> log(0.659^2 + 1) = 0.36067
    # No off-diagonal correlations were retained in the final model
    # (Methods Base model: "Random effect parameters for interindividual
    # variability in clearance (CL), volume of central compartment (Vc),
    # and rate of absorption (Ka) were identified.").
    etalcl ~ 0.16839  # Moes 2016 Table 4 Final model IIV CL = 42.8% CV (RSE 13%, shrinkage 0%)
    etalvc ~ 0.55657  # Moes 2016 Table 4 Final model IIV Vc = 86.3% CV (RSE 14%, shrinkage 9%)
    etalka ~ 0.36067  # Moes 2016 Table 4 Final model IIV Ka = 65.9% CV (RSE 14%, shrinkage 15%)

    # --- Residual unexplained variability ---
    # Moes 2016 Table 4 Final model sigma1 = "proportional error (%) = 13".
    propSd <- 0.13 ; label("Proportional residual error (fraction)")  # Moes 2016 Table 4 Final model sigma1 = 13% (RSE 8%, shrinkage 23%)
  })

  model({
    # --- Donor + recipient CYP3A5 combination effect on CL ---
    # Reconstruct the four-level combination from the two binary covariate
    # inputs (recipient + donor). Exactly one indicator is 1 for any given
    # subject; the others are 0. C1 (reference) contributes zero fractional
    # deviation, so the C1 indicator is implicit.
    ind_c2 <-      CYP3A5_EXPR  * (1 - CYP3A5_EXPR_DONOR)  # recipient *1, donor nonexpresser
    ind_c3 <- (1 - CYP3A5_EXPR) *      CYP3A5_EXPR_DONOR   # recipient nonexpresser, donor *1
    ind_c4 <-      CYP3A5_EXPR  *      CYP3A5_EXPR_DONOR   # both *1

    cyp3a5_factor <- e_cyp3a5_c2_cl * ind_c2 +
                     e_cyp3a5_c3_cl * ind_c3 +
                     e_cyp3a5_c4_cl * ind_c4

    # --- Individual PK parameters ---
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (1 + cyp3a5_factor)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # --- ODE system ---
    # Two-compartment disposition with a chain of three transit compartments
    # between the depot and the central compartment, all sharing the
    # absorption rate constant ka (Moes 2016 Results "Structural model
    # development" and Figure 1).
    d/dt(depot)       <- -ka * depot
    d/dt(transit1)    <-  ka * depot    - ka * transit1
    d/dt(transit2)    <-  ka * transit1 - ka * transit2
    d/dt(transit3)    <-  ka * transit2 - ka * transit3
    d/dt(central)     <-  ka * transit3 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # --- Bioavailability ---
    f(depot) <- exp(lfdepot)

    # --- Observation and error ---
    # Tacrolimus is measured by LC-MS/MS in whole blood and reported in
    # ug/L (= ng/mL). With dose in mg and Vc in L, central / vc gives mg/L,
    # which numerically equals ug/mL = 1000 ug/L; multiply by 1000 to
    # report ug/L directly.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
