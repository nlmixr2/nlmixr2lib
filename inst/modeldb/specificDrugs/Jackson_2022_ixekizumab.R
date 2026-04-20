Jackson_2022_ixekizumab <- function() {
  description <- "Two-compartment linear population PK model for subcutaneous ixekizumab in paediatric patients with moderate-to-severe plaque psoriasis (IXORA-PEDS; Jackson 2022)"
  reference <- "Jackson K, Chua L, Velez de Mendizabal N, et al. Population pharmacokinetic and exposure-efficacy analysis of ixekizumab in paediatric patients with moderate-to-severe plaque psoriasis (IXORA-PEDS). Br J Clin Pharmacol. 2022;88(3):1074-1086. doi:10.1111/bcp.15034"
  vignette <- "Jackson_2022_ixekizumab"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on CL, Q, V2, V3 with reference weight 58.6 kg (paediatric dataset reference; Jackson 2022 Table 2 footnotes). Baseline weight used; over-108-week weight change in 23% of subjects was not modelled (Discussion, page 1083).",
      source_name        = "WT"
    ),
    ADA_TITRE = list(
      description        = "Antidrug-antibody titre (continuous reciprocal dilution)",
      units              = "titre",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Log-linear multiplicative effect on CL per Jackson 2022 Table 2 footnote: CL_i = CL * (WT/58.6)^0.989 * (1 + 0.0292 * log_e[ADA titre]). ADA-negative samples (85.8% of the dataset) are encoded as ADA_TITRE = 1 so that log_e(1) = 0 cancels the effect (NONMEM convention; canonical encoding per inst/references/covariate-columns.md).",
      source_name        = "ADA titre"
    )
  )

  population <- list(
    n_subjects     = 184L,
    n_studies      = 1L,
    age_range      = "6-17 years",
    age_median     = "15 years (>50 kg), 10 years (25-50 kg), 7 years (<25 kg)",
    weight_range   = "21.5-136 kg",
    weight_median  = "65.2 kg (>50 kg group), 40 kg (25-50 kg group), 21.7 kg (<25 kg group); reference value used in the PopPK model is 58.6 kg",
    sex_female_pct = 56.5,
    race_ethnicity = c(White = 81.5, `Black or African American` = 3.26, Asian = 3.26, `American Indian or Alaska Native` = 1.63, Other = 7.61, Missing = 2.72),
    disease_state  = "Moderate-to-severe plaque psoriasis (PASI >=12, sPGA >=3, BSA >=10% at screening/baseline) in paediatric patients aged 6 to <18 years",
    dose_range     = "Weight-based Q4W SC: 20 mg (<25 kg), 40 mg (25-50 kg), 80 mg (>50 kg) after an initial 40/80/160 mg loading dose respectively",
    regions        = "US (37.5%), Europe (41.3%), Rest of World (21.2%)",
    weight_groups  = "4 patients <25 kg (2.2%), 45 patients 25-50 kg (24.5%), 135 patients >50 kg (73.4%)",
    ada_incidence  = "85.8% ADA-negative, 14.2% treatment-emergent-ADA positive (78 samples in 49 patients); maximum post-baseline titres 1:10 to 1:2560",
    injection_site = "Abdomen 28.6%, Arm 57.0%, Thigh 14.4% (injection site was not retained as a covariate in the paediatric model)",
    notes          = "IXORA-PEDS Phase 3 trial (NCT03073200). Baseline demographics per Jackson 2022 Table 1; 558 measurable IXE serum concentration measurements from 184 patients."
  )

  ini({
    # Structural parameters — reference body weight 58.6 kg (Jackson 2022 Table 2 footnotes)
    lka     <- log(0.00801);   label("Absorption rate constant (Ka, 1/hour)")                                     # Table 2: Ka = 0.00801 h^-1
    lcl     <- log(0.0120);    label("Apparent clearance at 58.6 kg reference weight, ADA-negative (CL, L/hour)") # Table 2: CL = 0.0120 L/h
    lvc     <- log(2.72);      label("Central volume of distribution at 58.6 kg reference weight (V2, L)")        # Table 2: V2 = 2.72 L
    lvp     <- log(2.11);      label("Peripheral volume of distribution at 58.6 kg reference weight (V3, L)")     # Table 2: V3 = 2.11 L
    lq      <- log(0.0119);    label("Intercompartmental clearance at 58.6 kg reference weight (Q, L/hour)")      # Table 2: Q = 0.0119 L/h
    lfdepot <- fixed(log(0.72)); label("Subcutaneous bioavailability (F1, fraction) -- fixed at adult model value") # Table 2 footnote c: F1 fixed at 0.72 (adult model)

    # Allometric weight exponents (reference weight 58.6 kg)
    allo_cl <- 0.989;          label("Allometric exponent on CL and Q (unitless)")                                # Table 2: weight on CL and Q = 0.989
    allo_v  <- 0.998;          label("Allometric exponent on V2 and V3 (unitless)")                               # Table 2: weight on V2 and V3 = 0.998

    # ADA-titre effect on CL (log-linear multiplicative)
    e_ada_cl <- 0.0292;        label("ADA-titre effect on CL, log-linear (unitless slope on log_e titre)")        # Table 2: ADA titre on CL = 0.0292

    # IIV on CL only (omega^2 = log(CV^2 + 1); paper reports 28.4% using %IIV = 100*sqrt(e^OMEGA - 1))
    etalcl ~ 0.07759           # Table 2: IIV on CL = 28.4%; log(1 + 0.284^2) = 0.07759

    # Residual error (proportional) -- Table 2: 27.7%
    propSd  <- 0.277;          label("Proportional residual error (fraction)")                                    # Table 2: residual error (proportional) = 27.7%
  })
  model({
    # ADA-titre effect on CL: multiplicative, log-linear. ADA-negative samples use ADA_TITRE = 1 so log_e(1) = 0.
    ada_cl <- 1 + e_ada_cl * log(ADA_TITRE)

    # Individual PK parameters with allometric weight scaling (reference 58.6 kg)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 58.6)^allo_cl * ada_cl
    vc <- exp(lvc)          * (WT / 58.6)^allo_v
    vp <- exp(lvp)          * (WT / 58.6)^allo_v
    q  <- exp(lq)           * (WT / 58.6)^allo_cl

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- exp(lfdepot)

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
