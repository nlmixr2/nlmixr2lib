Sanchez_2011_efavirenz <- function() {
  description <- "One-compartment population PK/pharmacogenetic model for oral efavirenz in Caucasian HIV-infected adults (Sanchez 2011), with GGT, CYP2B6*6 genotype (linked 516G>T + 785A>G), and ABCC4 (MRP4) 1497C>T carrier covariate effects on apparent oral clearance CL/F. Absorption rate ka fixed at 0.3 h^-1 (sparse TDM data could not estimate it); no covariate effect on V/F."
  reference <- paste(
    "Sanchez A, Cabrera S, Santos D, Valverde MP, Fuertes A,",
    "Dominguez-Gil A, Garcia MJ, and the Tormes Group.",
    "Population pharmacokinetic/pharmacogenetic model for optimization",
    "of efavirenz therapy in Caucasian HIV-infected patients.",
    "Antimicrob Agents Chemother. 2011;55(11):5616-5623.",
    "doi:10.1128/AAC.00194-11."
  )
  vignette <- "Sanchez_2011_efavirenz"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    GGT = list(
      description        = "Serum gamma-glutamyltransferase activity",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying biochemical liver-function covariate used as a linear-additive shift on apparent oral efavirenz clearance CL/F per Sanchez 2011 Table 4 / final-model equation: CL/F drops by 0.00279 L/h per +1 U/L of GGT. The effect is clinically relevant only at extreme values (Discussion paragraph 3: a 20% decrease in CL/F is reached around GGT 1,155 U/L). Cohort distribution: mean 121.21 +/- 156.79 U/L, range 8-1,612 U/L (Table 1).",
      source_name        = "GGT"
    ),
    SNP_CYP2B6_RS3745274_T_COUNT = list(
      description        = "Count of CYP2B6 c.516G>T (rs3745274, p.Q172H) T-alleles per subject (0/1/2). 0 = GG homozygous wild-type, 1 = GT heterozygous, 2 = TT homozygous variant. In the Sanchez 2011 Caucasian cohort the CYP2B6*6 haplotype was defined by the joint presence of 516G>T and 785A>G; the two SNPs are in tight linkage disequilibrium (Table 2 reports near-identical heterozygous/homozygous frequencies: 516G>T 60.80/32.80/6.40 % WT/het/hom and 785A>G 58.40/35.20/6.40 %), so *6 status maps one-to-one onto 516G>T genotype.",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (germline genotype). The source paper encodes two mutually-exclusive binary indicators 'CYP2B6*6 [G/T]' and 'CYP2B6*6 [T/T]' as multiplicative factors on CL/F (paper Table 4 / final-model equation). The canonical count column reconstructs them deterministically: *6 [G/T] = (count == 1), *6 [T/T] = (count == 2). Cohort distribution: 76 (60.80%) GG wild-type, 41 (32.80%) GT heterozygous, 8 (6.40%) TT homozygous (Table 2 row 516G>T; matching 785A>G frequencies in Table 2). The Sanchez 2011 paper combines these two SNPs into the *6 haplotype because both define the *6 allele; in any Caucasian cohort with the same tight 516/785 linkage disequilibrium, the canonical 516G>T count column is the right encoding (same precedent as Schipani_2011_nevirapine.R and Olagunju_2018_efavirenz.R).",
      source_name        = "CYP2B6*6 [G/T] and CYP2B6*6 [T/T] (paired indicators reconstructing the 516G>T T-allele count via *6 = joint 516G>T + 785A>G in tight LD)"
    ),
    SNP_ABCC4_1497CT_CARRIER = list(
      description        = "Binary carrier indicator for the ABCC4 (MRP4) c.1497C>T variant. 1 = subject carries at least one 1497T allele (heterozygous CT or homozygous TT); 0 = homozygous 1497CC wild-type. In the Sanchez 2011 cohort no 1497TT homozygotes were observed, so the indicator is effectively a heterozygous-vs-wild-type indicator in this dataset.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (homozygous 1497CC wild-type)",
      notes              = "Time-fixed per subject (germline genotype). Sanchez 2011 Table 2 reports 121 (96.80%) CC wild-type, 4 (3.20%) CT heterozygous, 0 (0.00%) TT homozygous out of n = 125 successfully genotyped patients. Multiplicative factor 0.793 on CL/F for 1497CT carriers (Discussion paragraph 8: 'EFV CL/F decreased by a factor of 0.79 for patients with a heterozygous genotype, possibly due to a decreased protein expression').",
      source_name        = "MRP4 1497C>T"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 128,
    n_studies      = 1,
    age_range      = "18-77 years",
    age_median     = "45 years (mean 45.06 +/- 9.16)",
    weight_range   = "39-113 kg",
    weight_median  = "65 kg (mean 64.98 +/- 12.20)",
    sex_female_pct = 32.82,
    race_ethnicity = c(Caucasian = 96.87, Other = 3.13),
    disease_state  = "HIV-positive adults on EFV-based antiretroviral therapy (600 mg EFV once daily plus two NRTIs) for at least 3 months at an unchanged dose for at least 1 month, with treatment adherence > 90% and no co-medication with known CYP inducers or inhibitors. 38.66% of analyzed concentrations were drawn from subjects with concomitant hepatitis C (HCV) co-infection.",
    dose_range     = "Initial dose 600 mg orally once daily; approximately 20% of subjects required TDM-driven dose adjustments in the range 200-1,600 mg/day (mean daily dose 608.75 +/- 104.36 mg/day; Table 1).",
    regions        = "Spain (University Hospital of Salamanca outpatient TDM clinic)",
    cyp2b6_freq    = "516G>T (rs3745274): GG 60.80%, GT 32.80%, TT 6.40% (n = 125). 785A>G: AA 58.40%, AG 35.20%, GG 6.40% (n = 125). The two SNPs are in tight linkage disequilibrium and jointly define the CYP2B6*6 haplotype that the paper uses as a paired heterozygous/homozygous indicator on CL/F (Sanchez 2011 Table 2).",
    abcc4_freq     = "ABCC4 (MRP4) 1497C>T: CC 96.80%, CT 3.20%, TT 0.00% (n = 125). Low minor-allele frequency in this Caucasian cohort, with no 1497TT homozygotes observed.",
    notes          = "869 EFV plasma concentrations from 128 patients (mean 4.59 +/- 2.84 samples per patient). Sparse therapeutic-drug-monitoring (TDM) data drawn at the midpoint of the dosing interval, 8-20 h post dose at steady state. Mean observed EFV concentration 3.18 +/- 1.61 ug/mL (range 0.84-15.16 ug/mL). 90 SNPs (CYP2A6, CYP2B6, CYP2C19, CYP2C8, CYP2C9, CYP2D6, CYP3A4, CYP3A5, MDR1, MRP1, MRP2, MRP4, UGT2B7, ABCA1, BCRP) plus 12 demographic / biochemical covariates were screened; only GGT, CYP2B6*6, and MRP4 1497C>T were retained in the final CL/F model (paper Results 'final model adopted for CL/F' paragraph; Table 4)."
  )

  ini({
    # ---- Structural PK parameters (Sanchez 2011 Table 4 final-model estimates) ----
    lka <- fixed(log(0.3)) ; label("First-order absorption rate constant ka (1/h) -- fixed")  # Sanchez 2011 Methods 'Population PK/PG model development' paragraph 2: ka 'could not be estimated and was fixed at 0.3 h^-1, a Ka value previously reported (13)' (Cabrera 2009)
    lcl <- log(12.2)       ; label("CL/F (L/h) -- intercept at GGT = 0, CYP2B6*6 wild-type, ABCC4 1497CC wild-type")  # Sanchez 2011 Table 4 theta_1 = 12.2 L/h (RSE 4.36%)
    lvc <- log(247)        ; label("V/F (L)")                                                  # Sanchez 2011 Table 4 theta_2 = 247 L (RSE 14.21%); no covariate effects retained in the final V/F model (Results 'Regarding V/F' paragraph)

    # ---- Covariate effects on CL/F (paper Table 4 final-model equation) ----
    # CL/F = (theta_1 + theta_3 * GGT) * theta_4^I(*6 het) * theta_5^I(*6 hom) * theta_6^I(ABCC4 carrier)
    # GGT enters linearly-additive on the intercept; the three polymorphism factors enter
    # multiplicatively (raised to a 0/1 indicator power so the factor is applied only when the
    # corresponding genotype indicator equals 1).
    e_ggt_cl          <- -0.00279 ; label("Linear-additive slope on CL/F per +1 U/L of GGT (L/h per U/L)")  # Sanchez 2011 Table 4 theta_3 = -0.00279 L^2/(U*h) (RSE 35.39%); negative slope -- higher GGT (cholestatic liver dysfunction) lowers CL/F
    e_cyp2b6_6het_cl  <-  0.602   ; label("Multiplicative factor on CL/F for CYP2B6*6 heterozygotes (unitless)") # Sanchez 2011 Table 4 theta_4 = 0.602 (RSE 7.83%); G/T heterozygotes have ~40% lower CL/F than wild-type
    e_cyp2b6_6hom_cl  <-  0.354   ; label("Multiplicative factor on CL/F for CYP2B6*6 homozygotes (unitless)")   # Sanchez 2011 Table 4 theta_5 = 0.354 (RSE 15.61%); T/T homozygotes have ~65% lower CL/F than wild-type
    e_abcc4_1497ct_cl <-  0.793   ; label("Multiplicative factor on CL/F for ABCC4 (MRP4) 1497C>T carriers (unitless)") # Sanchez 2011 Table 4 theta_6 = 0.793 (RSE 12.52%); carriers have ~21% lower CL/F than 1497CC wild-type

    # ---- IIV (diagonal omega; proportional / log-normal errors per Methods paragraph 1) ----
    # Sanchez 2011 Results 'final model adopted for CL/F' paragraph and Table 4:
    #   CV CL/F = 28.40% (RSE 18.29%); CV V/F = 86.91% (RSE 20.59%).
    # CV-to-variance conversion: omega^2 = log(CV^2 + 1)
    etalcl ~ 0.07758  # log(0.284^2 + 1) = 0.07758
    etalvc ~ 0.56253  # log(0.8691^2 + 1) = 0.56253

    # ---- Residual error (proportional; additive structure did not improve fit, Methods para 2) ----
    propSd <- 0.1682 ; label("Proportional residual error (fraction)")  # Sanchez 2011 Table 4 sigma = 16.82% (RSE 7.86%); proportional residual error model
  })

  model({
    # 1. Decompose the canonical 0/1/2 CYP2B6 516G>T allele count into the paired
    #    heterozygous and homozygous *6 indicators used by the paper. In the
    #    Sanchez 2011 Caucasian cohort, 516G>T genotype is concordant with *6
    #    haplotype because 516G>T and 785A>G are in tight linkage disequilibrium
    #    (see covariateData notes).
    het_6 <- (SNP_CYP2B6_RS3745274_T_COUNT == 1)
    hom_6 <- (SNP_CYP2B6_RS3745274_T_COUNT == 2)

    # 2. Typical CL/F: linear-additive intercept + GGT slope, with three
    #    multiplicative genotype factors raised to their 0/1 indicator powers
    #    (so each factor reduces to 1 when the indicator = 0).
    #    Reproduces Sanchez 2011 Table 4 / final-model equation literally.
    tvcl <- (exp(lcl) + e_ggt_cl * GGT) *
            e_cyp2b6_6het_cl^het_6 *
            e_cyp2b6_6hom_cl^hom_6 *
            e_abcc4_1497ct_cl^SNP_ABCC4_1497CT_CARRIER

    # 3. Individual PK parameters
    ka <- exp(lka)
    cl <- tvcl * exp(etalcl)
    vc <- exp(lvc + etalvc)

    # 4. Micro-constants
    kel <- cl / vc

    # 5. ODE system: one-compartment with first-order absorption (oral)
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 6. Observation and error
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
