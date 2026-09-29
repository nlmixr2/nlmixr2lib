Wang_2019_lamotrigine <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order absorption and first-order elimination for oral lamotrigine (LTG) in 89 Chinese patients with epilepsy (4-63 years) sampled sparsely at steady-state trough (Wang 2019 Table 4, final model). Ka was FIXED at 1.97 1/h from the literature because sampling carried no absorption-phase information. Apparent oral clearance (CL/F) carries three multiplicative fractional covariate factors: concomitant valproic acid (-38.6%), concomitant rifampicin (+64.7%) and the SLC22A1 1222G>A (rs628031) AA genotype (-52.5%). Apparent volume (V/F) carries two: the ABCG2 34G>A (rs2231137) AA genotype (-42.0%) and the combined MDR1/ABCB1 2677TT + 3435TT genotype pair (+139%). Interindividual variability is exponential on CL/F and V/F; residual error is combined additive plus proportional."
  reference <- paste(
    "Wang Z-z, Zhang Y-f, Huang W-c, Wang X-p, Ni X-j, Lu H-y, Hu J-q,",
    "Deng S-h, Zhu X-q, Xie H-s, Chen H-z, Zhang M, Qiu C, Wen Y-g, Shang D-w.",
    "Effects of comedication and genetic factors on the population",
    "pharmacokinetics of lamotrigine: a prospective analysis in Chinese",
    "patients with epilepsy. Front Pharmacol. 2019;10:832.",
    "doi:10.3389/fphar.2019.00832. PMCID PMC6669232.",
    "Parameters from Table 4 (final model); cohort demographics from Table 1;",
    "genotype frequencies from Tables 2 and 3.",
    sep = " "
  )
  vignette <- "Wang_2019_lamotrigine"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "lamotrigine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "lamotrigine", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    CONMED_VPA = list(
      description = "Concomitant valproic acid therapy",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant valproic acid)",
      notes = "Recorded per observation: 'Comedication status was recorded when the LTG concentration was determined during the study' (Subjects and Methods, 'Patients and Study Design'), so the indicator is time-varying within a subject. 157 of 419 observations (37.5%) were on LTG + VPA and a further 19 (4.5%) on LTG + VPA + RFP (Table 1). Enters CL/F as the fractional factor (1 + e_conmed_vpa_cl * CONMED_VPA) with e_conmed_vpa_cl = -0.386, i.e. a 38.6% reduction, consistent with valproate's inhibition of lamotrigine glucuronidation.",
      source_name = "VPA"
    ),
    CONMED_RIF = list(
      description = "Concomitant rifampicin therapy",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant rifampicin)",
      notes = "Recorded per observation at the LTG sampling time point, so time-varying within a subject; 57 of 419 observations (13.6%) were on LTG + RFP and a further 19 (4.5%) on LTG + VPA + RFP (Table 1). Wang 2019 defines no induction-onset lag: blood was drawn only at steady state after at least 7 days on the same dose (Methods, 'Patients and Study Design'), so the indicator represents chronic co-administration at post-induction equilibrium. Enters CL/F as the fractional factor (1 + e_conmed_rif_cl * CONMED_RIF) with e_conmed_rif_cl = 0.647, i.e. a 64.7% increase, attributed to induction of UGT1A4 / UGT2B7 glucuronidation.",
      source_name = "RFP"
    ),
    SNP_SLC22A1_RS628031_HOM = list(
      description = "SLC22A1 (OCT1) 1222G>A (rs628031) homozygous-variant indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (pooled 1222GG wild-type homozygotes and 1222GA heterozygotes)",
      notes = "1 = 1222AA genotype, 0 = GG or GA. Recessive encoding: Wang 2019 retained only the AA stratum. Cohort frequencies (Table 3): GG 41 (46.1%), GA 43 (48.3%), AA 5 (5.6%). Enters CL/F as (1 + e_snp_slc22a1_rs628031_hom_cl * SNP_SLC22A1_RS628031_HOM) with e_snp_slc22a1_rs628031_hom_cl = -0.525, a 52.5% reduction. OCT1 mediates active uptake of lamotrigine into hepatocytes, so a loss-of-function variant slows clearance (Discussion, citing Dickens 2012 and Shikata 2007).",
      source_name = "SLC22A1-1222AA"
    ),
    SNP_ABCG2_RS2231137_HOM = list(
      description = "ABCG2 (BCRP) 34G>A (rs2231137, V12M) homozygous-variant indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (pooled 34GG wild-type homozygotes and 34GA heterozygotes)",
      notes = "1 = 34AA genotype, 0 = GG or GA. Recessive encoding: Wang 2019 retained only the AA stratum. Cohort frequencies (Results, 'Frequencies of UGTs, MDR1, ABBC2, and SLC22A1 Genotyping Variants'): GG 36.0%, AA 13.5%. Enters V/F as (1 + e_snp_abcg2_rs2231137_hom_vc * SNP_ABCG2_RS2231137_HOM) with e_snp_abcg2_rs2231137_hom_vc = -0.420, a 42.0% reduction. Distinct SNP from the register's rs2231142 (Q141K, c.421C>A) entries.",
      source_name = "ABCG2-34AA"
    ),
    SNP_ABCB1_RS2032582_HOM = list(
      description = "ABCB1 / MDR1 2677G>T/A (rs2032582) homozygous 2677TT indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (any genotype other than 2677TT)",
      notes = "1 = 2677TT genotype, 0 = otherwise. 12 of 89 patients (13.5%) were TT homozygotes (Results). Used ONLY in combination with SNP_ABCB1_RS1045642_HOM: Wang 2019 combined the two SNPs into one covariate variable because '2677TT carriers were likely to have C3435T TT genotypes (91.67%)' (Discussion), so model() multiplies the two indicators to form the joint 2677TT + 3435TT stratum. rs2032582 is tri-allelic (2677G>T/A); the 2677A allele is pooled into the 0 reference here because Wang 2019 defines the covariate on the TT genotype alone.",
      source_name = "MDR1-2677TT"
    ),
    SNP_ABCB1_RS1045642_HOM = list(
      description = "ABCB1 / MDR1 3435C>T (rs1045642) homozygous-variant indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (pooled 3435CC wild-type homozygotes and 3435CT heterozygotes)",
      notes = "1 = 3435TT genotype, 0 = CC or CT. Cohort frequencies (Results): CC 47.2%, CT 38.2%, TT 14.6%. Used ONLY in combination with SNP_ABCB1_RS2032582_HOM (see that entry); the product of the two indicators enters V/F as (1 + e_snp_abcb1_2677tt_3435tt_vc * ...) with e_snp_abcb1_2677tt_3435tt_vc = 1.390, i.e. a 139% increase (2.39-fold).",
      source_name = "MDR1-C3435TT"
    )
  )

  # Covariates screened by Wang 2019 but NOT retained in the final model. Documented
  # here for provenance; they are deliberately absent from model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened as a linear mean-centred continuous covariate (Methods) and not retained. The Discussion explains why: van Dijkman 2018 found adult and paediatric clearance per kg to be close (0.0319 vs 0.0374 L/h/kg), so the 4-63 year span carried no age signal in this cohort."
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened as a linear mean-centred continuous covariate and not retained; Supplementary Figure 1 shows adolescent and adult body weights were similar, which the Discussion offers as a partial explanation. No allometric term appears in the final model."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened with the fractional discrete-covariate model (Methods coded COVdis = 1 for male, 0 for female) and not retained."
    ),
    CRCL = list(
      description = "Renal function (serum creatinine based)",
      units = "mL/min",
      type = "continuous",
      notes = "Serum creatinine was collected from patient charts and renal function was screened as a covariate (Methods); not retained."
    ),
    SMOKE = list(
      description = "Current-smoker indicator",
      units = "(binary)",
      type = "binary",
      notes = "Smoking status was recorded and screened (Methods); not retained."
    ),
    SNP_UGT1A4_RS2011425 = list(
      description = "UGT1A4 142T>G (rs2011425) genotype",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained. 'No obvious genetic effect of UGT enzymes was found relative to the concentrations of LTG in Chinese patients' (Abstract). Cohort TT 61.8% / TG 32.6% / GG 5.6% (Table 2)."
    ),
    SNP_UGT2B7_RS7668258 = list(
      description = "UGT2B7 -161C>T (rs7668258) genotype",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained (Table 2); representative of the four UGT2B7 SNPs (rs4356975, rs7668258, rs7662029, rs7438135) that were all screened and all rejected."
    ),
    SNP_ABCG2_RS2231142 = list(
      description = "ABCG2 421C>A (rs2231142, Q141K) genotype",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained: 'In our final model, ABCG2 421C>A had no effect on LTG concentrations' (Discussion). Distinct from the retained rs2231137 (34G>A) indicator."
    ),
    SNP_ABCC2_RS2273697 = list(
      description = "ABCC2 1249G>A (rs2273697) genotype",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained. Cohort GG 80.9% / GA 16.9% / AA 2.2% (Table 3)."
    ),
    SNP_SLC22A1_RS2282143 = list(
      description = "SLC22A1 1022C>T (rs2282143) genotype",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained; only 2 patients (2.2%) were TT carriers (Table 3)."
    ),
    SNP_ABCB1_RS1128503 = list(
      description = "ABCB1 / MDR1 1236C>T (rs1128503) genotype",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained as an independent covariate. The Results narrative lists MDR1-C1236T among the transporter polymorphisms identified for V/F, but Table 4 and the final-model equations carry only the combined 2677TT + 3435TT covariate; see vignette Errata."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 89,
    n_studies = 1,
    n_observations = 419,
    age_range = "4-63 years; mean 28 years (Table 1). The authors restrict the model's recommended applicability to 13-65 years because very young children were under-represented (Discussion).",
    weight_range = "15-94 kg; mean 59 kg (Table 1).",
    sex_female_pct = 52.8,
    race_ethnicity = c(Asian = 100),
    disease_state = "Out-patients with epilepsy receiving lamotrigine titration therapy, as monotherapy or in combination.",
    dose_range = "6.25-300 mg per administration; mean 118 mg (Table 1).",
    regions = "China (single centre: Guangzhou Huiai Hospital, The Affiliated Brain Hospital of Guangzhou Medical University).",
    co_medication = "Lamotrigine alone 187 obs (44.6%), lamotrigine + valproic acid 157 (37.5%), lamotrigine + rifampicin 57 (13.6%), lamotrigine + both 19 (4.5%) (Table 1). Both comedications are retained covariates on CL/F.",
    notes = "Prospective therapeutic-drug-monitoring study, July 2014 to January 2017. Sampling was sparse and at steady state (after at least 7 days on the same dose), 'mainly at 0.5 and 1 h prior to the next dosing', so the dataset is essentially trough-only and carries no absorption-phase information - the reason Ka is FIXED. Serum LTG was measured by HPLC-MS/MS, calibrated 0.20-25.00 ug/mL with an LLOQ of 0.20 ug/mL. Duration of therapy mean 42.8 weeks (range 4-204). The final model was evaluated by NPDE (mean 0.0106, variance 0.953, p = 0.504) and by external validation against an independent retrospective TDM cohort of 114 patients with 384 concentrations (mean prediction error 2.8 mg/L, RMSE 4.8 mg/L)."
  )

  ini({
    # Structural parameters - Wang 2019 Table 4, 'Final model / Estimates' column.

    # Ka FIXED. "The absorption rate constant (Ka) was fixed to the previously
    # reported value of 1.97 h-1 because of the lack of serum samples around the
    # absorption phase (Milosheska et al., 2016)." (Methods, 'Population
    # Pharmacokinetic Modeling'). Table 4 prints '1.97 FIX' with no %CV and no
    # bootstrap CI. A sensitivity analysis at half and twice 1.97 left CL/F and
    # V/F essentially unchanged (Results and Discussion).
    lka <- fixed(log(1.97))
    label("Absorption rate constant Ka (1/h)")
    lcl <- log(1.12)
    label("Apparent oral clearance CL/F (L/h)") # Table 4: 1.12 (%CV 14.6, bootstrap 95% CI 0.95-1.52)
    lvc <- log(12.7)
    label("Apparent central volume of distribution V/F (L)") # Table 4: 12.7 (%CV 28.4, bootstrap 95% CI 9.57-20.38)

    # Covariate effects. Wang 2019 uses the fractional discrete-covariate model
    # Pij = Ptv,j * (1 + theta_j * COVdis) (Methods), and the final-model text
    # restates each effect as a multiplicative factor E = 1 + theta that equals 1
    # when the covariate is absent. Effects are multiplicative ACROSS covariates:
    # the Discussion notes that concurrent VPA and RFP 'offset' one another, which
    # holds only for a product (0.614 * 1.647 = 1.011).
    #
    # Where the Results paragraph's restated E values disagree with the Table 4
    # theta column (E_VPA printed as 0.624 rather than 0.614; E_ABCG2 as 0.680
    # rather than 0.580), the Table 4 estimates are used: they carry the %CV and
    # the bootstrap CI, they agree exactly with the three other E values, and they
    # are the values corroborated by the Abstract's own percentages. See the
    # vignette Errata.
    e_conmed_vpa_cl <- -0.386
    label("Fractional effect of concomitant valproic acid on CL/F") # Table 4 'theta VPA on CL/F': -0.386 (%CV 19.1, bootstrap 95% CI -0.55 to -0.25)
    e_conmed_rif_cl <- 0.647
    label("Fractional effect of concomitant rifampicin on CL/F") # Table 4 'theta RFP on CL/F': 0.647 (%CV 15.4, bootstrap 95% CI 0.48-0.86)
    e_snp_slc22a1_rs628031_hom_cl <- -0.525
    label("Fractional effect of SLC22A1 1222AA on CL/F") # Table 4 'theta SLC22A1-1222AA on CL/F': -0.525 (%CV 29.5, bootstrap 95% CI -0.79 to -0.19)
    e_snp_abcg2_rs2231137_hom_vc <- -0.420
    label("Fractional effect of ABCG2 34AA on V/F") # Table 4 'theta ABCG2-34AA on V/F': -0.420 (%CV 35.5, bootstrap 95% CI -0.65 to -0.14)
    e_snp_abcb1_2677tt_3435tt_vc <- 1.390
    label("Fractional effect of MDR1 2677TT + 3435TT on V/F") # Table 4 'theta MDR1-2677TT + C3435TT on V/F': 1.390 (%CV 43.0, bootstrap 95% CI 0.21-3.25)

    # IIV. Exponential ('The IIV was modeled with an exponential error model
    # assuming a normal distribution with a mean of zero and variance of omega_p^2',
    # Methods). Table 4 reports the final-model IIV as a percentage: CL 15.2%,
    # V 20.1%. Converted to log-normal variances via omega^2 = log(1 + CV^2):
    #   CL: log(1 + 0.152^2) = 0.02284
    #   V : log(1 + 0.201^2) = 0.03961
    # No eta on Ka: the trough-only sampling 'could not estimate Ka and the IIV of
    # Ka' (Discussion, limitations), which is why Ka is fixed.
    etalcl ~ 0.02284 # Table 4 'CL INTER VAR, %' final = 15.2 (%CV 113.2, bootstrap 95% CI 2.47-47.89)
    etalvc ~ 0.03961 # Table 4 'V INTER VAR, %' final = 20.1 (%CV 172.6, bootstrap 95% CI 1.00-98.95)

    # Residual error. Table 4 reports BOTH an additive and a proportional term for
    # the final model, so the combined form is encoded. The proportional term is
    # printed in the '%' -labelled row as 0.167, i.e. a fraction of 16.7%, which is
    # what its own bootstrap CI of 10.24-20.77 (%) brackets; the base-model row in
    # the same column is printed as a percent (21.8). Concentration units are mg/L
    # (= ug/mL, the paper's unit).
    addSd <- 0.797
    label("Additive residual error (mg/L)") # Table 4 'Additive error, mg/L' final: 0.797 (%CV 31.0, bootstrap 95% CI 0.28-1.35)
    propSd <- 0.167
    label("Proportional residual error (fraction)") # Table 4 'Proportional error, %' final: 0.167 (%CV 24.4, bootstrap 95% CI 10.24-20.77%)
  })

  model({
    # Ka carries no IIV (see ini()).
    ka <- exp(lka)

    # CL/F: three multiplicative fractional covariate factors, each reducing to 1
    # when the covariate is 0 (Wang 2019 Results, final-model equations).
    cl <- exp(lcl + etalcl) *
      (1 + e_conmed_vpa_cl * CONMED_VPA) *
      (1 + e_conmed_rif_cl * CONMED_RIF) *
      (1 + e_snp_slc22a1_rs628031_hom_cl * SNP_SLC22A1_RS628031_HOM)

    # V/F: the ABCG2 34AA factor, plus the COMBINED MDR1 covariate. Wang 2019
    # merged the 2677G>T and 3435C>T SNPs into a single variable
    # (MDR1-2677+C3435T) because 91.67% of 2677TT carriers were also 3435TT
    # (Discussion), so the joint indicator is the product of the two homozygous
    # indicators and is 1 only for subjects who are TT at both loci.
    mdr1TTTT <- SNP_ABCB1_RS2032582_HOM * SNP_ABCB1_RS1045642_HOM
    vc <- exp(lvc + etalvc) *
      (1 + e_snp_abcg2_rs2231137_hom_vc * SNP_ABCG2_RS2231137_HOM) *
      (1 + e_snp_abcb1_2677tt_3435tt_vc * mdr1TTTT)

    kel <- cl / vc

    # One compartment, first-order oral absorption and first-order elimination.
    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - kel * central

    # Dose in mg, V/F in L -> mg/L, the unit of the paper's 3-15 ug/mL therapeutic
    # reference range.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
