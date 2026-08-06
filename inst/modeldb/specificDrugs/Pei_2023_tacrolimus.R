Pei_2023_tacrolimus <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and",
    "elimination for oral tacrolimus in adult heart transplant recipients",
    "(Pei 2023 Pharmaceutics, Phoenix NLME 8.3). Apparent oral clearance",
    "CL/F carries four covariate effects: a power effect of total",
    "bilirubin, an exponential effect of concomitant voriconazole, an",
    "exponential three-level CYP3A5*3 (rs776746) genotype effect with the",
    "*3/*3 nonexpresser as reference, and an exponential IL-10 G-1082A",
    "(rs1800896) heterozygote effect. Absorption rate constant fixed at",
    "0.30 1/h; no covariate was retained on Vd/F. Exponential IIV on both",
    "CL/F and Vd/F, proportional residual error. This is the companion",
    "top-down model to the whole-body PBPK model of the same paper",
    "(Pei_2023_tacrolimus_pbpk.R); the two identify different covariates,",
    "which the paper's Discussion addresses explicitly."
  )
  reference <- paste(
    "Pei L, Li R, Zhou H, Du W, Gu Y, Jiang Y, Wang Y, Chen X, Sun J,",
    "Zhu J (2023). A Physiologically Based Pharmacokinetic Approach to",
    "Recommend an Individual Dose of Tacrolimus in Adult Heart Transplant",
    "Recipients. Pharmaceutics 15(11):2580.",
    "doi:10.3390/pharmaceutics15112580.",
    "Population PK final estimates from Supplementary Table S4."
  )
  vignette <- "Pei_2023_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Tacrolimus concentrations in Pei 2023 are whole-blood
  # (CMIA assay, Methods 2.3).
  compartmentData <- list(
    depot   = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    TBILI = list(
      description        = "Total bilirubin",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying laboratory value collected from post-transplant day",
        "1 to day 30 (Methods 2.2, 'biochemical parameters ... total",
        "bilirubin (TBIL)'). Table 2 gives the cohort median (IQR) as",
        "15.90 (12.13, 21.10) umol/L; the units column of Table 2 reads",
        "'umol/L'. Enters CL/F as the Phoenix NLME continuous-covariate",
        "power form (TBILI / 15.90)^-0.19, normalised at the cohort",
        "median. The paper does not print the covariate equation; the",
        "power form is the Phoenix NLME default for a continuous covariate",
        "on a log-normally distributed parameter, and it is the only",
        "reading consistent with the magnitude of the estimate: an",
        "exponential form exp(-0.19 * (TBILI - 15.90)) would place a",
        "63 percent drop in CL/F across the interquartile range, which is",
        "not compatible with the reported IIV. Results 3.1: 'The apparent",
        "clearance (CL/F) of tacrolimus showed a negative correlation with",
        "TBIL, voriconazole, CYP3A5*3 (rs776746), and IL-10 G-1082A",
        "(rs1800896).'"
      ),
      source_name        = "TBIL"
    ),
    CONMED_VORICONAZOLE = list(
      description        = "Concomitant voriconazole coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant voriconazole)",
      notes              = paste(
        "Time-varying. Supplement 2.1.1: 'Categorical covariates, such as",
        "SEX (male = 1 and female = 2) and Voriconazole",
        "(co-administration = 1 and none = 0), were included in the",
        "analysis using indicator variables.' Voriconazole is a strong",
        "CYP3A inhibitor and reduces tacrolimus CL/F by exp(-0.64) = 0.53,",
        "i.e. a 47 percent decrease. Methods 2.2 lists voriconazole first",
        "among the concomitant CYP3A-inhibiting medications recorded",
        "(voriconazole, fluconazole, calcium antagonists, proton pump",
        "inhibitors, amiodarone); only voriconazole was retained in the",
        "final model."
      ),
      source_name        = "Voriconazole"
    ),
    CYP3A5_STAR1_HET = list(
      description        = "CYP3A5*1/*3 heterozygote indicator (one functional CYP3A5*1 allele)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5*3/*3 nonexpresser, when paired with CYP3A5_STAR1_HOM = 0)",
      notes              = paste(
        "Time-fixed per subject (germline genotype at rs776746). Table S4",
        "reports two separate rs776746 coefficients ('rs776746-TC' 0.77",
        "and 'rs776746-TT' 1.18), so the three genotype strata carry",
        "distinct typical values and the paired HET / HOM decomposition",
        "applies (SLCO1B1_HAP15_HET / _HOM precedent). Genotype counts",
        "(Table 2, n = 86 genotyped): CC 38 (44.2 percent), CT 45 (52.3",
        "percent), TT 3 (3.5 percent). The C allele is CYP3A5*3: its",
        "frequency here is (2*38 + 45) / 172 = 0.703, matching the",
        "reported CYP3A5*3 frequency in Chinese Han, and the resulting",
        "clearance ordering *3/*3 < *1/*3 < *1/*1 is monotone and matches",
        "Figure 1. CT therefore maps to CYP3A5_STAR1_HET = 1 and TT to",
        "CYP3A5_STAR1_HOM = 1, with CC (*3/*3) as the reference. Cross",
        "check: exp-weighting the two coefficients over the observed CT / TT",
        "counts gives (45 * exp(0.77) + 3 * exp(1.18)) / 48 = 2.23, against",
        "the paper's own reported ratio of Bayesian CL/F estimates,",
        "27.3 / 12.2 = 2.24 (Results 3.1)."
      ),
      source_name        = "rs776746-TC"
    ),
    CYP3A5_STAR1_HOM = list(
      description        = "CYP3A5*1/*1 homozygote indicator (two functional CYP3A5*1 alleles)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5*3/*3 nonexpresser, when paired with CYP3A5_STAR1_HET = 0)",
      notes              = paste(
        "Time-fixed per subject (germline genotype at rs776746). Maps to",
        "the TT genotype (3 of 86 genotyped subjects, 3.5 percent) - see",
        "the CYP3A5_STAR1_HET notes for the C = CYP3A5*3 orientation and",
        "the allele-frequency check. Table S4 'rs776746-TT' coefficient",
        "1.18 (CV 13.94 percent, bootstrap median 1.17, 95 percent CI",
        "0.99-1.62)."
      ),
      source_name        = "rs776746-TT"
    ),
    SNP_IL10_RS1800896_HET = list(
      description        = "IL-10 G-1082A (rs1800896) heterozygote indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (rs1800896 T/T homozygote)",
      notes              = paste(
        "Time-fixed per subject (germline genotype). Table 2 reports only",
        "two genotype strata at rs1800896 in the 86 genotyped subjects: TT",
        "63 (73.3 percent) and CT 23 (26.7 percent); no CC subjects were",
        "observed, so the single Table S4 coefficient ('rs1800896-TC',",
        "-0.35) is a heterozygote-versus-TT contrast and is encoded as a",
        "heterozygote indicator rather than an allele count. Carriers have",
        "exp(-0.35) = 0.70 times the CL/F of TT subjects. Discussion:",
        "'IL-10 is a potent modulator of CYP3A enzyme activity, inhibiting",
        "CYP3A-associated drug metabolism'."
      ),
      source_name        = "rs1800896-TC"
    )
  )

  covariatesDataExcluded <- list(
    SNP_CYP3A4_RS2242480_VAR_COUNT = list(
      description = "CYP3A4 rs2242480 (CYP3A4*18B) variant-allele count",
      units       = "(count, 0/1/2 alleles per subject)",
      type        = "continuous",
      notes       = paste(
        "Screened and reported as significant in univariate analysis -",
        "Results 3.1: 'CYP3A4*18B genotypes were also significant",
        "covariates of CL/F. For CYP3A4*18B, there was a significant",
        "difference in CL/F of different genotypes: 25.4 +/- 12.0 L/h for",
        "CYP3A4*18B*18B/*1*18B, and 16.3 +/- 10.5 L/h for CYP3A4*1*1",
        "(p < 0.001)' - but NOT retained in the final popPK model (Table",
        "S4 has no rs2242480 row). The Discussion explains why: 'the",
        "CYP3A5*3 mutant was also typically associated with the",
        "CYP3A4*18B mutation, indicating a linkage disequilibrium between",
        "the two genotypes', and the stepwise procedure kept only the",
        "covariate producing the largest objective-function reduction",
        "(Supplement 2.1.2). CYP3A4*18B does enter the companion PBPK",
        "model of the same paper (Pei_2023_tacrolimus_pbpk.R) through the",
        "FACYP3A4 activity level."
      )
    ),
    SEXF = list(
      description = "Sex, female indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Tested as an indicator variable (Supplement 2.1.1: 'SEX",
        "(male = 1 and female = 2)') but not retained in the final model.",
        "Note the paper's coding is 1 = male / 2 = female, so",
        "SEXF = SEX - 1 on the canonical scale. Discussion: 'the mean",
        "tacrolimus concentration was higher in female patients",
        "(3.77 +/- 2.82 ng/mL) compared to male patients",
        "(2.75 +/- 3.31 ng/mL) ... the variability of tacrolimus does not",
        "include gender in many studies, which may be due to the large",
        "difference in the male to female ratio' (93 men, 22 women)."
      )
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Collected (Methods 2.2 demographic characteristics) and used by",
        "the companion PBPK model, but not retained on either CL/F or",
        "Vd/F in the final popPK model (Table S4). Results 3.1: 'No",
        "relevant influencing variables were found for Vd/F.'"
      )
    ),
    HCT = list(
      description = "Hematocrit",
      units       = "% (volume fraction times 100)",
      type        = "continuous",
      notes       = paste(
        "Collected (Methods 2.2 biochemical parameters) and the second",
        "most sensitive input of the companion PBPK model (Table 5), but",
        "not retained in the final popPK model. The Discussion contrasts",
        "the two covariate sets directly: 'The main influencing factors of",
        "tacrolimus PK that were determined using the bottom-up approach",
        "were different from the results of the population PK analysis.'"
      )
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 115L,
    n_studies        = 1L,
    n_concentrations = 443L,
    age_range        = "adults >= 18 years; median (IQR) 52.00 (46.00, 61.00) years (Table 2)",
    age_median       = "52.00 years",
    weight_range     = "median (IQR) 67.50 (57.50, 75.00) kg (Table 2)",
    weight_median    = "67.50 kg",
    sex_female_pct   = 19.1,
    disease_state    = paste(
      "Adult heart transplant recipients at Nanjing First Hospital",
      "(November 2012 - January 2023) on tacrolimus (Prograf) with",
      "mycophenolate and corticosteroids; multi-organ transplant",
      "recipients and patients with incomplete clinical data excluded.",
      "Retrospective data from post-transplant day 1 to day 30; median",
      "(IQR) postoperative day 23.00 (19.00, 29.00). Concomitant",
      "CYP3A-inhibiting medications recorded were voriconazole,",
      "fluconazole, calcium antagonists, proton pump inhibitors and",
      "amiodarone."
    ),
    dose_range       = paste(
      "Oral tacrolimus titrated to a whole-blood trough target of 10-15",
      "ng/mL during the early postoperative period; median (IQR) daily",
      "dose 5.00 (4.00, 6.00) mg."
    ),
    regions          = "China (Nanjing First Hospital)",
    sampling_design  = paste(
      "443 steady-state whole-blood trough concentrations from 115",
      "recipients collected during routine therapeutic drug monitoring.",
      "Whole-blood tacrolimus by CMIA; LLOQ 2 ng/mL, quantitative range",
      "2-30 ng/mL."
    ),
    genotyping       = paste(
      "20 SNPs genotyped in 86 of the 115 subjects. rs35599367, rs1135840,",
      "rs150461093, rs2229109 and rs4253728 failed Hardy-Weinberg",
      "equilibrium (p < 0.05) and were excluded from the analysis.",
      "rs776746 genotype counts CC 38 / CT 45 / TT 3; rs1800896 counts",
      "TT 63 / CT 23; rs2242480 counts CC 45 / CT 40 / TT 1 (Table 2)."
    ),
    notes            = paste(
      "Software: Phoenix NLME 8.3 (Certara) with the FOCE algorithm.",
      "Structural-model selection by OFV / AIC / BIC; covariates by",
      "stepwise forward inclusion (p < 0.05) and backward exclusion",
      "(p < 0.01), keeping the covariate with the largest objective",
      "reduction when several were significant. Evaluated by",
      "prediction-corrected VPC (Figure S2) and a 1000-replicate",
      "bootstrap (Table S4). Final estimates and bootstrap medians agreed",
      "closely and all point estimates fell inside the 2.5-97.5th",
      "bootstrap percentiles. Ethics: Nanjing Medical University College",
      "Ethics Committee KY20190404-03-KS-01."
    )
  )

  ini({
    # ---- Structural parameters (Pei 2023 Table S4, final model) ----------
    # Ka was FIXED at 0.30 1/h - Table S4 prints "0.30(fixed)" in both the
    # final-model and bootstrap columns, with no CV%. Trough-only sampling
    # cannot identify an absorption rate constant.
    lka <- fixed(log(0.30)); label("First-order absorption rate constant (1/h)")  # Pei 2023 Table S4, Ka "0.30(fixed)"

    # Apparent oral clearance for the reference subject: CYP3A5*3/*3
    # (rs776746 CC), IL-10 rs1800896 T/T, no concomitant voriconazole, and
    # total bilirubin at the cohort median of 15.90 umol/L. Cross-check:
    # the paper's own Bayesian per-genotype summary gives 12.2 +/- 6.0 L/h
    # for CYP3A5*3/*3 (Results 3.1), matching this typical value.
    lcl <- log(12.35);  label("Apparent oral clearance CL/F for the reference subject (L/h)")  # Pei 2023 Table S4, CL/F 12.35 (CV 7.17%, bootstrap median 12.00, 95% CI 11.15-13.00)
    lvc <- log(656.80); label("Apparent volume of distribution Vd/F (L)")                      # Pei 2023 Table S4, Vd/F 656.80 (CV 13.65%, bootstrap median 669.13, 95% CI 511.39-744.68)

    # ---- Covariate effects on CL/F (Pei 2023 Table S4) -------------------
    # Total bilirubin, power form on the median-normalised covariate
    # (Phoenix NLME continuous-covariate default):
    #   CL/F *= (TBILI / 15.90)^e_tbili_cl
    e_tbili_cl <- -0.19; label("Power exponent of (TBILI / 15.90) on CL/F (unitless)")  # Pei 2023 Table S4, TBIL -0.19 (CV -18.78%, 95% CI -0.26 to -0.14)

    # Categorical effects, exponential form on the indicator (Phoenix NLME
    # categorical-covariate default): CL/F *= exp(theta * indicator).
    # The exponential reading is confirmed by the rs776746 coefficients:
    # exp-weighting 0.77 (CT, n = 45) and 1.18 (TT, n = 3) over the
    # observed counts gives a CYP3A5*1-carrier-to-*3/*3 CL/F ratio of 2.23,
    # matching the paper's own reported 27.3 / 12.2 = 2.24 (Results 3.1);
    # a linear (1 + theta) reading would give 1.80 and is refuted.
    e_vori_cl        <- -0.64; label("Exponential coefficient of concomitant voriconazole on CL/F (unitless)")            # Pei 2023 Table S4, Voriconazole -0.64 (CV -32.19%, 95% CI -0.91 to -0.44)
    e_cyp3a5_het_cl  <-  0.77; label("Exponential coefficient of CYP3A5*1/*3 (rs776746 CT) on CL/F (unitless)")           # Pei 2023 Table S4, rs776746-TC 0.77 (CV 11.59%, 95% CI 0.69-0.86)
    e_cyp3a5_hom_cl  <-  1.18; label("Exponential coefficient of CYP3A5*1/*1 (rs776746 TT) on CL/F (unitless)")           # Pei 2023 Table S4, rs776746-TT 1.18 (CV 13.94%, 95% CI 0.99-1.62)
    e_il10_het_cl    <- -0.35; label("Exponential coefficient of IL-10 rs1800896 heterozygotes on CL/F (unitless)")       # Pei 2023 Table S4, rs1800896-TC -0.35 (CV -27.58%, 95% CI -0.44 to -0.19)

    # ---- Inter-individual variability (Pei 2023 Table S4) ----------------
    # Supplement Eq S1: Pi = Pt * exp(eta_i), eta ~ N(0, omega^2). Table S4
    # reports the omega^2 values directly (the column header is
    # "omega_Vd/F^2" and "omega_CL/F^2"), so no CV-to-variance conversion
    # is needed.
    etalcl ~ 0.14   # Pei 2023 Table S4, omega^2 CL/F = 0.14 (CV 14.3%, bootstrap median 0.13, 95% CI 0.09-0.17)
    etalvc ~ 0.70   # Pei 2023 Table S4, omega^2 Vd/F = 0.70 (CV 35.7%, bootstrap median 0.80, 95% CI 0.23-1.36)

    # ---- Residual error (Pei 2023 Table S4) ------------------------------
    # Supplement Eq S2: Ci = C * (1 + eps), eps ~ N(0, sigma^2). Table S4
    # reports "Proportional (%) 1.48", read literally as sigma = 0.0148.
    # The alternative reading (1.48 = 148 percent) is excluded by the
    # Figure S2 pcVPC, whose 5th-95th prediction band is far too narrow for
    # a 148 percent residual and whose 5th percentile stays well above
    # zero.
    propSd <- 0.0148; label("Proportional residual error (fraction)")  # Pei 2023 Table S4, Proportional (%) 1.48 (CV 5.7%, bootstrap median 1.49, 95% CI 1.36-1.58)
  })

  model({
    # ---- 1. Covariate factors on CL/F ------------------------------------
    # Total bilirubin, power form normalised at the Table 2 cohort median.
    tbili_factor <- (TBILI / 15.90)^e_tbili_cl

    # Concomitant voriconazole, exponential form on the 0/1 indicator.
    vori_factor <- exp(e_vori_cl * CONMED_VORICONAZOLE)

    # CYP3A5*3 (rs776746) three-level genotype, exponential form on the
    # paired mutually exclusive indicators. The *3/*3 reference has both
    # indicators = 0, so the factor collapses to 1.
    cyp3a5_factor <- exp(e_cyp3a5_het_cl * CYP3A5_STAR1_HET +
                         e_cyp3a5_hom_cl * CYP3A5_STAR1_HOM)

    # IL-10 G-1082A (rs1800896) heterozygote, exponential form; T/T is the
    # reference.
    il10_factor <- exp(e_il10_het_cl * SNP_IL10_RS1800896_HET)

    # ---- 2. Individual parameters ----------------------------------------
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * tbili_factor * vori_factor *
          cyp3a5_factor * il10_factor
    vc <- exp(lvc + etalvc)

    # ---- 3. ODE system (one compartment, first-order absorption) ---------
    kel <- cl / vc
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ---- 4. Observation and error ----------------------------------------
    # central is in mg and vc in L, so central / vc is mg/L = ug/mL; the
    # factor 1000 converts to the ng/mL reported by Pei 2023.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
