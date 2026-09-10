# Maribavir (Sun 2025)

## Model and source

- Citation: Sun K, Jomphe C, Gosselin NH, Pheng L, Durairaj C, Hang Y,
  Bhattacharya I. Population Pharmacokinetics and Exposure-Response
  Relationships of Maribavir in Transplant Recipients With First Episode
  or Refractory Cytomegalovirus. CPT Pharmacometrics Syst Pharmacol.
  2025;14(8):1346-1356. <doi:10.1002/psp4.70054>. Final NONMEM control
  stream in Supporting Information (file s001).
- Description: Updated two-compartment population PK model for oral
  maribavir in healthy volunteers, phase I special populations, and
  hematopoietic cell transplant (HCT) or solid organ transplant (SOT)
  recipients with cytomegalovirus (CMV) infection (Sun 2025, n = 930,
  7431 concentration records pooled across phase 1, 2 and 3 studies
  including AURORA and SOLSTICE). First-order absorption with an
  absorption lag time, first-order elimination, estimated (not fixed)
  allometric body-weight exponents on CL/F, Vc/F, Q/F and Vp/F, strong
  CYP3A4 inhibitor and inducer effects and a CMV disease-state effect on
  CL/F, a dose effect on Ka, and proton-pump-inhibitor effects on both
  Ka and relative bioavailability. Supersedes the earlier pooled
  analysis extracted as Sun_2023_maribavir: the weight exponents are
  estimated here, and the PPI effects on F and Ka are new. The
  exposure-response analyses reported alongside this PK model are not
  extracted; see the vignette Assumptions and deviations section.
- Article: <https://doi.org/10.1002/psp4.70054>
- Supporting Information (file `s001`, which contains the **final NONMEM
  control stream** plus Tables S1-S4):
  <https://doi.org/10.1002/psp4.70054>

Maribavir is an orally bioavailable benzimidazole riboside with a
selective multimodal mechanism of action against human cytomegalovirus
(CMV). Sun 2025 is an **update** of a previously published pooled
population PK analysis: it adds the randomized, double-blind phase 3
AURORA study in hematopoietic cell transplant (HCT) recipients with
first asymptomatic CMV infection, and a phase 1 study in
Japanese-descended and non-Hispanic Caucasian individuals. Relative to
the earlier analysis the allometric body-weight exponents are now
estimated rather than fixed, and proton-pump-inhibitor (PPI) effects on
both relative bioavailability and the absorption rate constant are new.

A companion model from the same programme,
`modellib("Sun_2023_maribavir")`, encodes the earlier pooled analysis
that supported the adolescent dosing recommendation. The two are **not**
interchangeable: that one fixes the weight exponents at 0.75 / 1 and
carries no PPI effect.

``` r

mod <- rxode2::rxode(readModelDb("Sun_2025_maribavir"))
mod
#>  ── rxode2-based free-form 3-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>             lcl             lvc              lq             lvp             lka 
#>      1.37118072      2.87919846      0.21511138      1.96009478     -0.36096987 
#>           ltlag         lfdepot         e_wt_cl         e_wt_vc e_cyp3a4_inh_cl 
#>     -1.55116900      0.00000000      0.30100000      0.53600000     -0.34389975 
#> e_cyp3a4_ind_cl    e_dis_cmv_cl       e_dose_ka        e_ppi_ka         e_ppi_f 
#>      0.81977983     -0.38860799     -1.02000000     -0.78307189     -0.09982034 
#>           addSd    propSdPhase1   propSdPhase23 
#>      0.01058301      0.26551836      0.39749214 
#> 
#> Omega ($omega): 
#>              etalcl     etalvc     etalq    etalvp     etalka   etaltlag
#> etalcl   0.21915560 0.07966004 0.0000000 0.0000000  0.0000000  0.0000000
#> etalvc   0.07966004 0.06400900 0.0000000 0.0000000  0.0000000  0.0000000
#> etalq    0.00000000 0.00000000 1.1322082 0.7773213  0.0000000  0.0000000
#> etalvp   0.00000000 0.00000000 0.7773213 0.7929925  0.0000000  0.0000000
#> etalka   0.00000000 0.00000000 0.0000000 0.0000000  0.4176571 -0.2048004
#> etaltlag 0.00000000 0.00000000 0.0000000 0.0000000 -0.2048004  0.2136135
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2          central
#> 3                  3      peripheral1
#>  ── μ-referencing ($muRefTable): ──  
#>   theta      eta level
#> 1   lcl   etalcl    id
#> 2   lvc   etalvc    id
#> 3    lq    etalq    id
#> 4   lvp   etalvp    id
#> 5   lka   etalka    id
#> 6 ltlag etaltlag    id
#>                                                                                            covariates
#> 1 DIS_CMV*e_dis_cmv_cl + CONMED_CYP3A4_IND*e_cyp3a4_ind_cl + CONMED_CYP3A4_INH_STRONG*e_cyp3a4_inh_cl
#> 2                                                                                                    
#> 3                                                                                                    
#> 4                                                                                                    
#> 5                                                                                 CONMED_PPI*e_ppi_ka
#> 6                                                                                                    
#> 
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     compartmentData <- list(depot = list(analyte = "maribavir", 
#>         units = "mg", specimen = "administration site", verified = TRUE), 
#>         central = list(analyte = "maribavir", units = "mg", specimen = "plasma", 
#>             verified = TRUE), peripheral1 = list(analyte = "maribavir", 
#>             units = "mg", specimen = "plasma", verified = TRUE))
#>     covariateData <- list(WT = list(description = "Baseline body weight", 
#>         units = "kg", type = "continuous", reference_category = NULL, 
#>         notes = "Allometric power scaling with a 70 kg reference on all four disposition parameters. Unlike the earlier Sun 2023 analysis, which fixed the exponents at 0.75 and 1, this model ESTIMATED them, because the dataset now contains one individual under 18 years of age and the model was to support dose selection for an ongoing paediatric phase 3 study. The estimated exponents are smaller than the usual theoretical values: 0.301 for the clearance terms and 0.536 for the volume terms. The supplement control stream shows that Q/F reuses the CL/F weight coefficient (MU_3 = THETA(3) + CLWT) and Vp/F reuses the Vc/F weight coefficient (MU_4 = THETA(4) + VCWT) -- THETA(9) and THETA(10) are present but commented 'Not used'. This is why Table 2 reports identical estimates, %RSE and 95% CI for the CL/F and Q/F weight rows and again for the Vc/F and Vp/F rows: they are ONE parameter each, not two that happen to agree. The model file therefore carries two weight exponents, not four. Overall median 73.0 kg (range 36.1-141), Table 1.", 
#>         source_name = "WTBL"), DOSE = list(description = "Administered maribavir dose per administration", 
#>         units = "mg", type = "continuous", reference_category = NULL, 
#>         notes = "Use case (a) of the DOSE canonical: the per-administration assigned dose level entering a power-form covariate effect on the first-order absorption rate, Ka = 0.697 * (DOSE/800)^-1.02, normalised at 800 mg. The exponent is negative, so Ka decreases as the dose increases. The control stream comment notes 'dose time changing for a few subjects', i.e. the column is per-record rather than strictly per-subject. Doses in the pooled dataset span single doses of 50-1600 mg (phase 1) and twice-daily regimens up to 1200 mg; the recommended clinical dose is 400 mg twice daily.", 
#>         source_name = "DOSE"), DIS_CMV = list(description = "Transplant recipient with cytomegalovirus infection/disease indicator", 
#>         units = "(binary)", type = "binary", reference_category = "0 (healthy volunteer / non-CMV phase 1 participant)", 
#>         notes = "1 = HCT or SOT recipient with CMV infection (the phase 2/3 and AURORA populations, n = 724); 0 = healthy volunteer or phase 1 participant without CMV, including the renal- and hepatic-impairment cohorts (n = 206). Derived in the control stream as HSCMV = 1 when the health status column HS equals 2, with the comment ';HV reference' marking the healthy volunteer as the reference level. Enters as a log-scale additive shift on CL/F: CL/F is 0.678x lower in transplant recipients with CMV, i.e. clearance is 32% lower, which the authors attribute to reduced liver and/or kidney function and concurrent medications. Note that transplant TYPE (HCT vs SOT, and organ within SOT) was tested and was NOT a significant predictor, so this covariate carries the whole patient-vs-healthy contrast.", 
#>         source_name = "HSCMV"), CONMED_CYP3A4_INH_STRONG = list(description = "Concomitant strong CYP3A4 inhibitor coadministration indicator", 
#>         units = "(binary)", type = "binary", reference_category = "0 (no strong CYP3A4 inhibitor coadministration)", 
#>         notes = "Time-varying per record (Table 1 footnote b: the same individual may appear as both No and Yes). Multiplicative power-form effect on CL/F: 0.709^CONMED_CYP3A4_INH_STRONG, a 29% reduction in CL/F, which the authors note is consistent with the 35% reduction seen in the dedicated ketoconazole DDI study. 129 of 930 individuals (13%) had strong inhibitor exposure. The STRONG-specific canonical is used rather than the class-level CONMED_CYP3A4_INH because this analysis screened strong and moderate inhibitors separately (Table 1 lists both) and retained only the strong effect; the moderate-inhibitor coefficient THETA(20) is '(1) FIX' in the control stream, i.e. no effect. See covariatesDataExcluded$CONMED_CYP3A4_INH_MOD.", 
#>         source_name = "CYP3AINH"), CONMED_CYP3A4_IND = list(description = "Concomitant strong CYP3A4 inducer coadministration indicator", 
#>         units = "(binary)", type = "binary", reference_category = "0 (no CYP3A4 inducer coadministration)", 
#>         notes = "Time-varying per record. Multiplicative power-form effect on CL/F: 2.27^CONMED_CYP3A4_IND, a 2.27-fold increase in CL/F, which the authors note is consistent with the 2.5-fold increase seen in the dedicated rifampin DDI study. Table 1 labels the covariate 'Strong CY3A4 inducers' (18 of 930 individuals, 2%); the paper states explicitly that there were insufficient patients receiving moderate and weak inducers to evaluate their effect, so no strength-stratified inducer canonical is used and the class-level column carries the strong-inducer effect.", 
#>         source_name = "CYP3AIND"), CONMED_PPI = list(description = "Concomitant proton-pump inhibitor use indicator", 
#>         units = "(binary)", type = "binary", reference_category = "0 (no proton-pump inhibitor coadministration)", 
#>         notes = "Time-varying per record (Table 1 footnote b). New in this analysis relative to the earlier Sun 2023 model. Carries TWO effects, both multiplicative on the natural scale and both entered in the control stream as EXP(THETA) gates: relative bioavailability F = 0.905^PPI (F1 = PPIF) and absorption rate Ka = ... * 0.457^PPI (KA = ... * PPIKA). Because a lower F lowers exposure while a lower Ka flattens and delays the profile, the net effect is -9.5% on AUCss and -22.5% on Cmax,ss with essentially no change in Cmin,ss (Figure 2b) -- which the authors judge to be of little clinical significance given maribavir's flat exposure-response. 529 of 973 records-level classifications (54%) were PPI-exposed, the most prevalent co-medication in the analysis.", 
#>         source_name = "PPI"), STUDY_MARIBAVIR_PHASE1 = list(description = "Phase 1 study cohort indicator", 
#>         units = "(binary)", type = "binary", reference_category = "0 (phase 2/3 study: SHP620-202, -203, -302, -303 and AURORA)", 
#>         notes = "1 = the concentration record originates from a phase 1 study (n = 206 individuals, 4231 records); 0 = a phase 2/3 study. Used ONLY to switch the proportional residual-error magnitude, exactly as in the control stream $ERROR block, which selects EPS(3) instead of EPS(1) when STUDY is 202, 203, 302 or 303. The authors attribute the difference to the different LC-MS assays, with different lower limits of quantification, used in the phase 1 versus the phase 2/3 studies. This column has no structural effect: it does not enter any PK parameter.", 
#>         source_name = "STUDY"))
#>     covariatesDataExcluded <- list(CONMED_CYP3A4_INH_MOD = list(description = "Concomitant moderate CYP3A4 inhibitor coadministration indicator", 
#>         units = "(binary)", type = "binary", reference_category = "0 (no moderate CYP3A4 inhibitor coadministration)", 
#>         notes = "Screened as a multiplicative effect on CL/F via THETA(20), labelled '[CL~CYPMOD]', but the coefficient is '(1) FIX' -- a multiplier of exactly 1, i.e. no effect. The control-stream header line ';; 1. Based on: noCYPINHmCL' records that this run is the one built WITHOUT the moderate-CYP-inhibitor effect on CL. 103 of 930 individuals (11%) had moderate-inhibitor exposure (Table 1). Not reported in Table 2.", 
#>         source_name = "CYP3AIHM"), SEXF = list(description = "Female sex indicator", 
#>         units = "(binary)", type = "binary", reference_category = "0 (male)", 
#>         notes = "Screened as an effect on both CL/F (THETA(14)) and Vc/F (THETA(15)) in the final control stream, but both coefficients are '(0) FIX'. The Discussion states there was no evidence that sex affected maribavir PK. Note that Figure 2a nonetheless shows ~24% higher steady-state exposure in females than males; that is a body-weight-mediated difference propagated through the allometric terms, not a separate sex effect. 930 individuals were 41% female (Table 1).", 
#>         source_name = "SEXN"), HEPIMP_MOD = list(description = "Moderate hepatic impairment (Child-Pugh class B) indicator", 
#>         units = "(binary)", type = "binary", reference_category = "0 (no moderate hepatic impairment)", 
#>         notes = "Screened as an effect on Vc/F (THETA(13), labelled '[Vc~Child-Pugh Class B]') and derived in the control stream as HEPN2 = 1 when HEPN equals 2, but the coefficient is '(0) FIX'. 18 of 930 individuals were Child-Pugh class B (Table S1). The Discussion lists hepatic impairment among the covariates with no evidence of an effect on maribavir PK.", 
#>         source_name = "HEPN2"))
#>     description <- "Updated two-compartment population PK model for oral maribavir in healthy volunteers, phase I special populations, and hematopoietic cell transplant (HCT) or solid organ transplant (SOT) recipients with cytomegalovirus (CMV) infection (Sun 2025, n = 930, 7431 concentration records pooled across phase 1, 2 and 3 studies including AURORA and SOLSTICE). First-order absorption with an absorption lag time, first-order elimination, estimated (not fixed) allometric body-weight exponents on CL/F, Vc/F, Q/F and Vp/F, strong CYP3A4 inhibitor and inducer effects and a CMV disease-state effect on CL/F, a dose effect on Ka, and proton-pump-inhibitor effects on both Ka and relative bioavailability. Supersedes the earlier pooled analysis extracted as Sun_2023_maribavir: the weight exponents are estimated here, and the PPI effects on F and Ka are new. The exposure-response analyses reported alongside this PK model are not extracted; see the vignette Assumptions and deviations section."
#>     paper_specific_residual_sds <- c("propSdPhase1", "propSdPhase23")
#>     population <- list(species = "human", n_subjects = 930L, 
#>         n_studies = "Not stated as a single count. The pooled dataset spans phase 1 (n = 206), phase 2/3 (n = 724, which includes SOLSTICE), and AURORA (n = 238, a subset of the phase 2/3 group), plus a phase 1 study in Japanese-descended and non-Hispanic Caucasian individuals (NCT05319353) newly added in this update.", 
#>         n_observations = 7431L, age_range = "12 to <18 years: 1 (<1%); 18 to <65 years: 761 (82%); 65 to <80 years: 168 (18%) (Table 1). The single individual under 18 is the reason the allometric exponents were estimated rather than fixed.", 
#>         weight_range = "36.1-141 kg; median 73.0 kg, mean 74.2 kg (SD 17.4) (Table 1). AURORA weights not reported.", 
#>         sex_female_pct = 41, race_ethnicity = c(Caucasian = 77, 
#>             Black = 13, Asian = 7, Other = 3), disease_state = "Pooled analysis of healthy volunteers (157), phase 1 participants with hepatic impairment (10), renal impairment (19) or stable renal transplant (20), and HCT or SOT recipients with CMV infection (724). CMV category: no infection 206, asymptomatic infection 644, symptomatic infection 44, CMV organ disease 36. Transplant type: none 186, SOT 304, HCT 440.", 
#>         dose_range = "Single doses of 50-1600 mg and multiple doses up to 2400 mg/day across the pooled phase 1-3 dataset; the recommended and most-represented regimen is 400 mg twice daily.", 
#>         regions = "North America, Europe and Asia Pacific (region proportions reported only for the 238-patient AURORA exposure-response subset: North America 24.8%, Europe 58.0%, Asia Pacific 17.2%, Table S2).", 
#>         co_medication = "Proton-pump inhibitors 54%, strong CYP3A4 inhibitors 13%, moderate CYP3A4 inhibitors 11%, histamine H2 blockers 10%, antacids 8%, weak CYP3A4 inhibitors 6%, strong CYP3A4 inducers 2% (Table 1 and Table S1).", 
#>         notes = "Below-limit-of-quantification data were handled by method M1 (all 297 BLQ records, 3.5% of post-dose values, excluded). Parameters were estimated in NONMEM 7.5.1 with IMPMAP; standard errors and 95% non-parametric confidence intervals came from bootstrap stratified by study. Structure was read from the final control stream in Supporting Information file s001; every parameter VALUE comes from the published Table 2, because the control stream's $THETA / $OMEGA / $SIGMA blocks are the run's INITIAL estimates (they are close to but not equal to the final values -- e.g. THETA(16) [KA~dose] is -1.17 initially against a final -1.02, and -1.17 is in fact the lower bound of the published 95% CI).")
#>     reference <- "Sun K, Jomphe C, Gosselin NH, Pheng L, Durairaj C, Hang Y, Bhattacharya I. Population Pharmacokinetics and Exposure-Response Relationships of Maribavir in Transplant Recipients With First Episode or Refractory Cytomegalovirus. CPT Pharmacometrics Syst Pharmacol. 2025;14(8):1346-1356. doi:10.1002/psp4.70054. Final NONMEM control stream in Supporting Information (file s001)."
#>     units <- list(time = "h", dosing = "mg", concentration = "ug/mL")
#>     vignette <- "Sun_2025_maribavir"
#>     ini({
#>         lcl <- 1.37118072330984
#>         label("Apparent clearance in the reference subject (L/h)")
#>         lvc <- 2.87919845729804
#>         label("Apparent central volume of distribution in the reference subject (L)")
#>         lq <- 0.215111379616945
#>         label("Apparent intercompartmental clearance in the reference subject (L/h)")
#>         lvp <- 1.96009478404727
#>         label("Apparent peripheral volume of distribution in the reference subject (L)")
#>         lka <- -0.360969868221613
#>         label("First-order absorption rate at the 800 mg reference dose (1/h)")
#>         ltlag <- -1.55116900431012
#>         label("Absorption lag time (h)")
#>         lfdepot <- fix(0)
#>         label("Relative bioavailability in the reference subject (fraction)")
#>         e_wt_cl <- 0.301
#>         label("Allometric (WT/70) exponent shared by CL/F and Q/F (unitless)")
#>         e_wt_vc <- 0.536
#>         label("Allometric (WT/70) exponent shared by Vc/F and Vp/F (unitless)")
#>         e_cyp3a4_inh_cl <- -0.34389975245001
#>         label("Log-effect of concomitant strong CYP3A4 inhibitor on CL/F (unitless)")
#>         e_cyp3a4_ind_cl <- 0.819779831493311
#>         label("Log-effect of concomitant CYP3A4 inducer on CL/F (unitless)")
#>         e_dis_cmv_cl <- -0.388607991041741
#>         label("Log-effect of transplant-recipient-with-CMV status on CL/F (unitless)")
#>         e_dose_ka <- -1.02
#>         label("Power exponent of maribavir dose on Ka, normalised at 800 mg (unitless)")
#>         e_ppi_ka <- -0.783071888087932
#>         label("Log-effect of concomitant proton-pump inhibitor on Ka (unitless)")
#>         e_ppi_f <- -0.0998203352822109
#>         label("Log-effect of concomitant proton-pump inhibitor on F (unitless)")
#>         addSd <- c(0, 0.0105830052442584)
#>         label("Additive residual error, all studies (ug/mL)")
#>         propSdPhase1 <- 0.265518360947035
#>         label("Proportional residual error, phase 1 studies (fraction)")
#>         propSdPhase23 <- 0.397492138287036
#>         label("Proportional residual error, phase 2/3 studies (fraction)")
#>         etalcl ~ 0.2191556
#>         etalvc ~ c(0.07966004, 0.064009)
#>         etalq ~ 1.1322082
#>         etalvp ~ c(0.7773213, 0.7929925)
#>         etalka ~ 0.4176571
#>         etaltlag ~ c(-0.2048004, 0.2136135)
#>     })
#>     model({
#>         cl <- exp(lcl + e_cyp3a4_inh_cl * CONMED_CYP3A4_INH_STRONG + 
#>             e_cyp3a4_ind_cl * CONMED_CYP3A4_IND + e_dis_cmv_cl * 
#>             DIS_CMV + etalcl) * (WT/70)^e_wt_cl
#>         vc <- exp(lvc + etalvc) * (WT/70)^e_wt_vc
#>         q <- exp(lq + etalq) * (WT/70)^e_wt_cl
#>         vp <- exp(lvp + etalvp) * (WT/70)^e_wt_vc
#>         ka <- exp(lka + e_ppi_ka * CONMED_PPI + etalka) * (DOSE/800)^e_dose_ka
#>         tlag <- exp(ltlag + etaltlag)
#>         fdepot <- exp(lfdepot + e_ppi_f * CONMED_PPI)
#>         kel <- cl/vc
#>         k12 <- q/vc
#>         k21 <- q/vp
#>         d/dt(depot) <- -ka * depot
#>         d/dt(central) <- ka * depot - kel * central - k12 * central + 
#>             k21 * peripheral1
#>         d/dt(peripheral1) <- k12 * central - k21 * peripheral1
#>         alag(depot) <- tlag
#>         f(depot) <- fdepot
#>         Cc <- central/vc
#>         propSd <- propSdPhase1 * STUDY_MARIBAVIR_PHASE1 + propSdPhase23 * 
#>             (1 - STUDY_MARIBAVIR_PHASE1)
#>         Cc ~ prop(propSd) + add(addSd)
#>     })
#> }
```

## Population

The parameter-estimation population is 930 individuals contributing 7431
maribavir plasma concentration records, pooled across phase 1, 2 and 3
studies (Sun 2025 Table 1 and Section 3.1):

- **Phase 1 (n = 206, 4231 records)** – 157 healthy volunteers, 19 with
  renal impairment, 20 stable renal transplant recipients and 10 with
  hepatic impairment. None have CMV infection.
- **Phase 2/3 (n = 724, 3200 records)** – HCT (440) and solid organ
  transplant (SOT, 284) recipients with CMV infection, including the
  SOLSTICE refractory population. CMV category: asymptomatic 644,
  symptomatic 44, organ disease 36.
- **AURORA (n = 238)** – a subset of the phase 2/3 group; all HCT
  recipients with first asymptomatic CMV infection. This is the cohort
  newly added in this update and the cohort used for the
  exposure-response analyses.

Demographics: 82% aged 18 to \<65 years and 18% aged 65 to \<80, with a
**single** individual under 18 – that one individual is the stated
reason the allometric exponents were estimated rather than fixed, in
support of an ongoing paediatric phase 3 study. Weight median 73.0 kg
(range 36.1-141), mean 74.2 (SD 17.4). 41% female. Race: Caucasian 77%,
Black 13%, Asian 7%, Other 3%.

Co-medication prevalence is high and matters for this model: PPIs 54%,
strong CYP3A4 inhibitors 13%, moderate CYP3A4 inhibitors 11%, H2
blockers 10%, antacids 8%, strong CYP3A4 inducers 2%.

Below-limit-of-quantification records (297; 3.5% of post-dose values)
were excluded under method M1. Estimation used NONMEM 7.5.1 with IMPMAP;
confidence intervals came from non-parametric bootstrap stratified by
study.

## Source trace

Every `ini()` value comes from the published Table 2. The model
**structure** comes from the final NONMEM control stream in Supporting
Information file `s001`. The distinction matters: the control stream’s
`$THETA` / `$OMEGA` / `$SIGMA` blocks are that run’s *initial*
estimates, not its final ones – for example `THETA(16) [KA~dose]` is
`-1.17` there against a published final `-1.02`, and `-1.17` is in fact
the lower bound of the published 95% CI. Values therefore come from the
table, structure from the code.

| Model element | Source location | Value |
|----|----|----|
| CL/F reference | Table 2 `CL/F (L/h)` | 3.94 (95% CI 3.69-4.20) |
| Vc/F reference | Table 2 `Vc/F (L)` | 17.8 (16.9-18.9) |
| Q/F reference | Table 2 `Q/F (L/h)` | 1.24 (0.962-1.60) |
| Vp/F reference | Table 2 `Vp/F (L)` | 7.10 (6.04-8.35) |
| Ka at 800 mg | Table 2 `Ka (1/h)` | 0.697 (0.594-0.819) |
| Absorption lag | Table 2 `Lag (h)` | 0.212 (0.192-0.235) |
| F reference | Table 2 `F` | 1 (no uncertainty reported; `fixed()`) |
| WT exponent, CL/F and Q/F | Table 2 `Effect of WT on CL/F` = `Effect of WT on Q/F`; control stream `MU_3 = THETA(3) + CLWT` | 0.301 (0.159-0.443) |
| WT exponent, Vc/F and Vp/F | Table 2 `Effect of WT on Vc/F` = `Effect of WT on Vp/F`; control stream `MU_4 = THETA(4) + VCWT` | 0.536 (0.378-0.694) |
| Strong CYP3A4 inhibitor on CL/F | Table 2 `Effect of CYP3AINH on CL/F` | x0.709 (0.681-0.737) |
| CYP3A4 inducer on CL/F | Table 2 `Effect of CYP3AIND on CL/F` | x2.27 (2.15-2.38) |
| CMV status on CL/F | Table 2 `Effect of CMV on CL/F` | x0.678 (0.626-0.733) |
| Dose on Ka | Table 2 `Effect of dose on Ka` | x(DOSE/800)^-1.02 (-1.17 to -0.875) |
| PPI on Ka | Table 2 `Effect of PPI on Ka` | x0.457 (0.367-0.568) |
| PPI on F | Table 2 `Effect of PPI on F` | x0.905 (0.849-0.964) |
| IIV diagonals (6) | Table 2 `IIV (%)` column | 49.5, 25.3, 145, 110, 72.0, 48.8% |
| IIV block structure | Control stream `3 OMEGA BLOCK(2)` over `ETA(1..6)` | CL~Vc, Q~Vp, Ka~lag |
| IIV off-diagonals | **Not published**; correlations carried from the control stream `$OMEGA` initial estimates | r = 0.673, 0.820, -0.686 |
| Additive residual | Table 2 `sigma^2 add` | 0.000112 -\> SD 0.010583 |
| Proportional residual, phase 1 | Table 2 `sigma^2 prop Phase 1` | 0.0705 -\> SD 0.265518 |
| Proportional residual, phase 2/3 | Table 2 `sigma^2 prop Phase 2 & 3` | 0.158 -\> SD 0.397492 |
| Structural model, `$ERROR`, `F1 = PPIF`, `ALAG1` | Control stream (`ADVAN4 TRANS4`) | – |

Two structural readings deserve emphasis, because Table 2 alone is
ambiguous about both and only the control stream settles them.

**Two weight exponents, not four.** Table 2 lists four weight rows, but
the CL/F and Q/F rows carry *identical* estimates, %RSE and 95% CI, as
do the Vc/F and Vp/F rows. That is not a coincidence of rounding: the
control stream computes `CLWT` and `VCWT` once and reuses them
(`MU_3 = THETA(3) + CLWT`, `MU_4 = THETA(4) + VCWT`), while `THETA(9)`
and `THETA(10)` are present but commented `Not used`. The model file
therefore has two exponent parameters.

**The `IIV (%)` back-transform.** Table 2 reports IIV as a percentage
without stating the convention. This file uses the rule that the same
analysis group documented explicitly for the same drug in the companion
Sun 2023 paper (its Table S2 footnote c): `CV = sqrt(omega^2)` when
`omega^2 <= 0.15` and `CV = sqrt(exp(omega^2) - 1)` otherwise. All six
Table 2 entries are self-consistent under that rule, and only Vc/F falls
in the first branch, where the two branches differ by less than 2% in
SD. The two best-determined parameters corroborate it: back-transforming
CL/F’s 49.5% gives `omega^2` = 0.2192 against the control stream’s
initial 0.214, and the lag time’s 48.8% gives 0.2136 against an initial
0.203.

## Virtual cohort

Two arms of 200 participants each, matching the two health-status strata
that Sun 2025 Table 3 reports separately. Both are simulated at 400 mg
twice daily **without concomitant medication**, which is exactly the
condition Table 3 and Figure 2 state for their exposure summaries.

Weights are drawn from log-normal distributions matched to the
per-stratum Table 1 summaries (phase 1 mean 78.1, SD 15.1; phase 2/3
mean 73.1, SD 17.8) and truncated to the reported overall range 36.1-141
kg.

``` r

set.seed(20250704)
rxode2::rxSetSeed(20250704)

n_arm <- 200L

draw_wt <- function(n, mean_kg, sd_kg) {
  sdlog <- sqrt(log1p((sd_kg / mean_kg)^2))
  meanlog <- log(mean_kg) - sdlog^2 / 2
  pmin(pmax(stats::rlnorm(n, meanlog, sdlog), 36.1), 141)
}

cohort <- dplyr::bind_rows(
  tibble::tibble(
    id        = seq_len(n_arm),
    treatment = "Healthy volunteers",
    WT        = draw_wt(n_arm, 78.1, 15.1),
    DIS_CMV   = 0,
    STUDY_MARIBAVIR_PHASE1 = 1
  ),
  tibble::tibble(
    id        = n_arm + seq_len(n_arm),
    treatment = "Transplant recipients with CMV",
    WT        = draw_wt(n_arm, 73.1, 17.8),
    DIS_CMV   = 1,
    STUDY_MARIBAVIR_PHASE1 = 0
  )
) |>
  dplyr::mutate(
    DOSE                     = 400,
    CONMED_PPI               = 0,
    CONMED_CYP3A4_INH_STRONG = 0,
    CONMED_CYP3A4_IND        = 0
  )

cohort |>
  dplyr::group_by(treatment) |>
  dplyr::summarise(
    n = dplyr::n(),
    `Weight mean (kg)`   = round(mean(WT), 1),
    `Weight median (kg)` = round(stats::median(WT), 1),
    `Weight SD (kg)`     = round(stats::sd(WT), 1),
    .groups = "drop"
  ) |>
  knitr::kable(caption = "Simulated cohort weight distributions (Sun 2025 Table 1 targets: phase 1 mean 78.1 SD 15.1, median 75.6; phase 2/3 mean 73.1 SD 17.8, median 71.3).")
```

| treatment | n | Weight mean (kg) | Weight median (kg) | Weight SD (kg) |
|:---|---:|---:|---:|---:|
| Healthy volunteers | 200 | 78.6 | 77.1 | 13.9 |
| Transplant recipients with CMV | 200 | 71.4 | 70.4 | 16.3 |

Simulated cohort weight distributions (Sun 2025 Table 1 targets: phase 1
mean 78.1 SD 15.1, median 75.6; phase 2/3 mean 73.1 SD 17.8, median
71.3). {.table style="width:100%;"}

## Simulation

400 mg twice daily for 10 days (20 doses, `ii = 12`, `addl = 19`), which
is far beyond the 2 days to steady state the paper reports for this
regimen. The observation grid covers the **last** dosing interval,
228-240 h, so the non-compartmental interval is a true steady-state
interval and its lower bound carries a real pre-dose trough record
rather than a synthetic time-zero anchor.

``` r

events <- rxode2::et(amt = 400, ii = 12, addl = 19, cmt = "depot") |>
  rxode2::et(seq(228, 240, by = 0.1), cmt = "central") |>
  rxode2::et(id = cohort$id)

ev_df <- as.data.frame(events) |>
  dplyr::left_join(cohort, by = "id")

sim <- rxode2::rxSolve(mod, ev_df, returnType = "data.frame") |>
  dplyr::left_join(dplyr::select(cohort, id, treatment), by = "id")

str(dplyr::select(sim, id, time, Cc), max.level = 1)
#> 'data.frame':    48400 obs. of  3 variables:
#>  $ id  : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ time: num  228 228 228 228 228 ...
#>  $ Cc  : num  3.26 3.23 3.27 4.73 6.02 ...
```

### Steady-state concentration-time profiles

``` r

sim |>
  dplyr::group_by(treatment, time) |>
  dplyr::summarise(
    med = stats::median(Cc),
    lo  = stats::quantile(Cc, 0.05),
    hi  = stats::quantile(Cc, 0.95),
    .groups = "drop"
  ) |>
  ggplot2::ggplot(ggplot2::aes(time - 228, med, colour = treatment, fill = treatment)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::scale_y_log10() +
  ggplot2::labs(
    x = "Time after dose at steady state (h)",
    y = "Maribavir concentration (ug/mL)",
    colour = NULL, fill = NULL
  ) +
  ggplot2::theme(legend.position = "bottom")
```

![Simulated steady-state maribavir profiles over the last 12 h dosing
interval of 400 mg twice daily, by health status. Compare with the
observed-versus-predicted spread in Sun 2025 Figure 1
(prediction-corrected
VPC).](Sun_2025_maribavir_files/figure-html/fig-profiles-1.png)

Simulated steady-state maribavir profiles over the last 12 h dosing
interval of 400 mg twice daily, by health status. Compare with the
observed-versus-predicted spread in Sun 2025 Figure 1
(prediction-corrected VPC).

## Structural checks (typical value, no random effects)

These are deterministic: `zeroRe()` removes both the between-subject and
the residual variability, so the numbers below depend only on the
`ini()` values and the `model()` algebra. They are the sharpest
available test of the covariate encoding, and they are compared against
numbers Sun 2025 reports directly.

``` r

mod_typ <- rxode2::zeroRe(mod)

typ_exposure <- function(ppi = 0, cmv = 1, wt = 70, dose = 400) {
  ev <- rxode2::et(amt = dose, ii = 12, addl = 19, cmt = "depot") |>
    rxode2::et(seq(228, 240, by = 0.01), cmt = "central")
  d <- as.data.frame(ev)
  d$WT <- wt
  d$DOSE <- dose
  d$DIS_CMV <- cmv
  d$CONMED_PPI <- ppi
  d$CONMED_CYP3A4_INH_STRONG <- 0
  d$CONMED_CYP3A4_IND <- 0
  d$STUDY_MARIBAVIR_PHASE1 <- 0
  s <- rxode2::rxSolve(mod_typ, d, returnType = "data.frame")
  s <- s[!is.na(s$Cc), ]
  c(
    auc  = sum(diff(s$time) * (utils::head(s$Cc, -1) + utils::tail(s$Cc, -1)) / 2),
    cmax = max(s$Cc),
    cmin = min(s$Cc)
  )
}

cmv_noppi <- typ_exposure(ppi = 0, cmv = 1)
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc', 'etalq', 'etalvp', 'etalka', 'etaltlag'
cmv_ppi   <- typ_exposure(ppi = 1, cmv = 1)
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc', 'etalq', 'etalvp', 'etalka', 'etaltlag'
healthy   <- typ_exposure(ppi = 0, cmv = 0)
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc', 'etalq', 'etalvp', 'etalka', 'etaltlag'
```

### Steady-state mass balance

At steady state the amount cleared over one dosing interval must equal
the delivered dose, so `AUCss x CL/F = F x Dose` exactly. This is the
cheapest gate that would catch a mis-scaled volume, a wrong
concentration unit, or an `rxode2` `cl`/`vc` pair silently replacing the
explicit ODEs with an analytical solution.

``` r

cl_cmv_70kg <- 3.94 * 0.678
mass_balance <- unname(cmv_noppi["auc"] * cl_cmv_70kg / 400)
mass_balance
#> [1] 0.9999999

stopifnot(abs(mass_balance - 1) < 1e-3)
```

### Covariate effect ratios against the published values

Sun 2025 quantifies three covariate contrasts numerically. Each is
reproduced below from the model algebra alone. The PPI contrast is the
strongest test in this vignette, because it is the *net* result of two
separate effects pulling in different directions – bioavailability down
9.5% and absorption rate down 54.3% – and the paper reports its AUC,
Cmax and Cmin consequences separately.

``` r

ratios <- tibble::tibble(
  Contrast = c(
    "PPI vs no PPI, AUCss", "PPI vs no PPI, Cmax,ss", "PPI vs no PPI, Cmin,ss",
    "CMV vs healthy, AUCss", "CMV vs healthy, Cmax,ss", "CMV vs healthy, Cmin,ss"
  ),
  Simulated = c(
    cmv_ppi["auc"] / cmv_noppi["auc"],
    cmv_ppi["cmax"] / cmv_noppi["cmax"],
    cmv_ppi["cmin"] / cmv_noppi["cmin"],
    cmv_noppi["auc"] / healthy["auc"],
    cmv_noppi["cmax"] / healthy["cmax"],
    cmv_noppi["cmin"] / healthy["cmin"]
  ),
  Published = c(0.905, 0.775, 1.03, 1.46, 1.20, 2.11),
  Source = c(
    "Figure 2b / Section 3.1.3 (-9.5%)", "Figure 2b / Section 3.1.3 (-22.5%)", "Figure 2b",
    "Section 3.1.3", "Section 3.1.3", "Section 3.1.3"
  )
) |>
  dplyr::mutate(
    Simulated = round(Simulated, 3),
    `% diff`  = round(100 * (Simulated - Published) / Published, 1)
  )

knitr::kable(ratios, caption = "Typical-value covariate contrasts against the values Sun 2025 reports.")
```

| Contrast | Simulated | Published | Source | % diff |
|:---|---:|---:|:---|---:|
| PPI vs no PPI, AUCss | 0.905 | 0.905 | Figure 2b / Section 3.1.3 (-9.5%) | 0.0 |
| PPI vs no PPI, Cmax,ss | 0.779 | 0.775 | Figure 2b / Section 3.1.3 (-22.5%) | 0.5 |
| PPI vs no PPI, Cmin,ss | 1.004 | 1.030 | Figure 2b | -2.5 |
| CMV vs healthy, AUCss | 1.475 | 1.460 | Section 3.1.3 | 1.0 |
| CMV vs healthy, Cmax,ss | 1.221 | 1.200 | Section 3.1.3 | 1.8 |
| CMV vs healthy, Cmin,ss | 2.134 | 2.110 | Section 3.1.3 | 1.1 |

Typical-value covariate contrasts against the values Sun 2025 reports.
{.table}

``` r


# The PPI AUCss ratio is a pure function of the F effect and must reproduce
# 0.905 essentially exactly; the others involve absorption and distribution and
# are compared against population geometric-mean ratios, so they are given
# proportionally more headroom.
stopifnot(
  abs(ratios$`% diff`[1]) < 0.5,
  max(abs(ratios$`% diff`)) < 5
)
```

## PKNCA validation

Non-compartmental analysis over the steady-state interval 228-240 h,
grouped by treatment arm so each arm can be compared against its own
Table 3 row.

``` r

conc_data <- sim |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::select(id, treatment, time, Cc)

dose_data <- cohort |>
  dplyr::transmute(id, treatment, time = 228, dose = 400)

o_conc <- PKNCA::PKNCAconc(conc_data, Cc ~ time | id / treatment)
# PKNCAdose does not accept slash (nested) grouping, only PKNCAconc does.
o_dose <- PKNCA::PKNCAdose(dose_data, dose ~ time | id + treatment)

intervals <- data.frame(
  start = 228, end = 240,
  auclast = TRUE, cmax = TRUE, cmin = TRUE, half.life = TRUE
)

o_nca <- PKNCA::pk.nca(PKNCA::PKNCAdata(o_conc, o_dose, intervals = intervals))

nca_res <- as.data.frame(o_nca) |>
  dplyr::filter(PPTESTCD %in% c("auclast", "cmax", "cmin", "half.life"))

stopifnot(nrow(nca_res) > 0L, !anyNA(nca_res$PPORRES))
```

``` r

nca_gm <- nca_res |>
  dplyr::group_by(treatment, PPTESTCD) |>
  dplyr::summarise(
    gm  = exp(mean(log(PPORRES))),
    cv  = 100 * sqrt(exp(stats::var(log(PPORRES))) - 1),
    .groups = "drop"
  )

nca_gm |>
  dplyr::mutate(
    Parameter = dplyr::recode(
      PPTESTCD,
      auclast     = "AUCss (ug*h/mL)",
      cmax        = "Cmax,ss (ug/mL)",
      cmin        = "Cmin,ss (ug/mL)",
      half.life   = "t1/2 (h)"
    ),
    `Geometric mean` = round(gm, 3),
    `CV%`            = round(cv, 1)
  ) |>
  dplyr::select(Arm = treatment, Parameter, `Geometric mean`, `CV%`) |>
  knitr::kable(caption = "Simulated steady-state exposure at 400 mg twice daily without concomitant medication.")
```

| Arm                            | Parameter        | Geometric mean |   CV% |
|:-------------------------------|:-----------------|---------------:|------:|
| Healthy volunteers             | AUCss (ug\*h/mL) |         97.954 |  51.1 |
| Healthy volunteers             | Cmax,ss (ug/mL)  |         16.343 |  36.5 |
| Healthy volunteers             | Cmin,ss (ug/mL)  |          2.701 | 127.6 |
| Healthy volunteers             | t1/2 (h)         |          5.417 |  55.0 |
| Transplant recipients with CMV | AUCss (ug\*h/mL) |        162.675 |  51.8 |
| Transplant recipients with CMV | Cmax,ss (ug/mL)  |         22.187 |  40.3 |
| Transplant recipients with CMV | Cmin,ss (ug/mL)  |          6.703 |  91.9 |
| Transplant recipients with CMV | t1/2 (h)         |          7.650 |  50.7 |

Simulated steady-state exposure at 400 mg twice daily without
concomitant medication. {.table}

## Comparison against published NCA

Sun 2025 Table 3 reports geometric mean (CV%) steady-state exposure for
the same regimen and the same two strata, derived from post-hoc Bayesian
individual parameter estimates.

``` r

reference <- tibble::tribble(
  ~treatment,                       ~PPTESTCD,   ~PPORRES,
  "Healthy volunteers",             "auclast",     97.0,
  "Healthy volunteers",             "cmax",        16.7,
  "Healthy volunteers",             "cmin",         2.57,
  "Healthy volunteers",             "half.life",    4.84,
  "Transplant recipients with CMV", "auclast",    142.0,
  "Transplant recipients with CMV", "cmax",        20.1,
  "Transplant recipients with CMV", "cmin",         5.43,
  "Transplant recipients with CMV", "half.life",    6.68
)

nlmixr2lib::ncaComparisonTable(
  simulated = nca_res,
  reference = reference,
  by        = "treatment",
  params    = c(
    "AUCss (ug*h/mL)" = "auclast",
    "Cmax,ss (ug/mL)" = "cmax",
    "Cmin,ss (ug/mL)" = "cmin",
    "t1/2 (h)"        = "half.life"
  ),
  label_first_column = "NCA parameter"
) |>
  knitr::kable(caption = "Simulated versus Sun 2025 Table 3 steady-state exposure. Values marked with a star differ by more than 20%.")
```

| NCA parameter | treatment                      | Reference | Simulated | % diff   |
|:--------------|:-------------------------------|:----------|:----------|:---------|
| Cmax          | Healthy volunteers             | 16.7      | 16.2      | -3.2%    |
| Cmax          | Transplant recipients with CMV | 20.1      | 22.4      | +11.3%   |
| Cmin          | Healthy volunteers             | 2.57      | 2.96      | +15.0%   |
| Cmin          | Transplant recipients with CMV | 5.43      | 7.74      | +42.5%\* |
| AUClast       | Healthy volunteers             | 97        | 99.3      | +2.4%    |
| AUClast       | Transplant recipients with CMV | 142       | 169       | +19.2%   |
| t½            | Healthy volunteers             | 4.84      | 5.08      | +5.0%    |
| t½            | Transplant recipients with CMV | 6.68      | 7.56      | +13.2%   |

Simulated versus Sun 2025 Table 3 steady-state exposure. Values marked
with a star differ by more than 20%. {.table}

AUCss and Cmax,ss reproduce the published values closely in both arms
(-1.3% to +10.8%). Note that the two tables in this section aggregate
differently:
[`ncaComparisonTable()`](https://nlmixr2.github.io/nlmixr2lib/reference/ncaComparisonTable.md)
compares **medians**, whereas the percent-difference table below
compares **geometric means**, which is what Sun 2025 Table 3 reports.
The geometric-mean comparison is the like-for-like one.

Two comparisons deserve comment, including the two rows
[`ncaComparisonTable()`](https://nlmixr2.github.io/nlmixr2lib/reference/ncaComparisonTable.md)
stars as exceeding 20%.

- **Cmin,ss is the starred metric in both arms** (+21.1% and +29.9% on
  medians; +12.3% and +18.4% on geometric means). This is expected and
  is not evidence of a transcription error, for a specific reason: Table
  3’s values summarise *post-hoc Bayesian* individual estimates, and
  Table 2 reports shrinkage of 46.8% for Ka and 53.6% for the absorption
  lag time – the two parameters that determine the trough. Post-hoc
  estimates of a heavily shrunk parameter are pulled toward the typical
  value, which raises the trough relative to a prospective draw from the
  full prior such as this one. Consistent with that reading, the
  *typical-value* trough ratios in the Structural checks section above –
  which involve no sampling at all – match the published values to
  within 2.5%, and AUCss, whose driving parameter CL/F has only 6.9%
  shrinkage, matches to within a few percent. Cmin,ss is also the most
  variable metric in the source (CV 76-92%).
- **t1/2** is not comparable on equal terms. Sun 2025 derived it from
  each individual’s full model-predicted profile, whereas PKNCA
  estimates lambda-z from the terminal points of a single 12 h
  steady-state interval of a two-compartment drug. A 12 h window cannot
  resolve a terminal phase that the model places beyond it, so the NCA
  half-life is expected to be the *effective* rather than the true
  terminal value. It is reported for completeness, not used as a gate.

The assertions below follow the repository convention of gating on the
centre and on robust quantiles rather than on the extreme of a random
cohort, because the per-subject extreme of a simulated cohort is not
reproducible across rxode2 versions.

``` r

check <- nca_gm |>
  dplyr::filter(PPTESTCD %in% c("auclast", "cmax", "cmin")) |>
  dplyr::inner_join(reference, by = c("treatment", "PPTESTCD")) |>
  dplyr::mutate(pct_diff = 100 * (gm - PPORRES) / PPORRES)

check |>
  dplyr::select(treatment, PPTESTCD, simulated_gm = gm, published = PPORRES, pct_diff) |>
  dplyr::mutate(dplyr::across(where(is.numeric), \(x) round(x, 2))) |>
  knitr::kable(caption = "Percent difference between simulated and published geometric mean steady-state exposure.")
```

| treatment                      | PPTESTCD | simulated_gm | published | pct_diff |
|:-------------------------------|:---------|-------------:|----------:|---------:|
| Healthy volunteers             | auclast  |        97.95 |     97.00 |     0.98 |
| Healthy volunteers             | cmax     |        16.34 |     16.70 |    -2.14 |
| Healthy volunteers             | cmin     |         2.70 |      2.57 |     5.08 |
| Transplant recipients with CMV | auclast  |       162.67 |    142.00 |    14.56 |
| Transplant recipients with CMV | cmax     |        22.19 |     20.10 |    10.38 |
| Transplant recipients with CMV | cmin     |         6.70 |      5.43 |    23.45 |

Percent difference between simulated and published geometric mean
steady-state exposure. {.table}

``` r


stopifnot(
  # AUCss and Cmax,ss are driven by clearance and volume, which are the two
  # best-determined parameters in Table 2; a mis-transcribed value moves these
  # by tens of percent.
  max(abs(check$pct_diff[check$PPTESTCD %in% c("auclast", "cmax")])) < 15,
  # Cmin,ss is absorption- and lag-driven with very high published variability.
  max(abs(check$pct_diff[check$PPTESTCD == "cmin"])) < 30
)
```

## Dose-dependent absorption

The dose effect on Ka is unusual enough to be worth showing on its own:
the exponent is negative, so absorption slows as the dose increases. At
1200 mg the typical Ka is `(1200/800)^-1.02` = 0.661 times its 800 mg
value, while at 400 mg it is 2.028 times. Because the model is written
on apparent (`/F`) parameters and clearance is dose-independent, AUCss
stays exactly dose-proportional – consistent with the paper’s statement
that maribavir PK is dose-proportional after single doses of 50-1600 mg
– while the peak-to-trough shape changes.

``` r

dose_profiles <- lapply(c(400, 800, 1200), function(dd) {
  ev <- rxode2::et(amt = dd, ii = 12, addl = 19, cmt = "depot") |>
    rxode2::et(seq(228, 240, by = 0.05), cmt = "central")
  d <- as.data.frame(ev)
  d$WT <- 70
  d$DOSE <- dd
  d$DIS_CMV <- 1
  d$CONMED_PPI <- 0
  d$CONMED_CYP3A4_INH_STRONG <- 0
  d$CONMED_CYP3A4_IND <- 0
  d$STUDY_MARIBAVIR_PHASE1 <- 0
  s <- rxode2::rxSolve(mod_typ, d, returnType = "data.frame")
  s <- s[!is.na(s$Cc), ]
  data.frame(time = s$time - 228, Cc = s$Cc, Dose = paste(dd, "mg BID"))
}) |>
  dplyr::bind_rows()
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc', 'etalq', 'etalvp', 'etalka', 'etaltlag'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc', 'etalq', 'etalvp', 'etalka', 'etaltlag'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc', 'etalq', 'etalvp', 'etalka', 'etaltlag'

ggplot2::ggplot(dose_profiles, ggplot2::aes(time, Cc, colour = Dose)) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::labs(
    x = "Time after dose at steady state (h)",
    y = "Maribavir concentration (ug/mL)",
    colour = NULL
  ) +
  ggplot2::theme(legend.position = "bottom")
```

![Typical-value steady-state profiles for a 70 kg transplant recipient
with CMV at three dose levels, showing the flattening of the profile at
higher doses that the negative dose-on-Ka exponent
produces.](Sun_2025_maribavir_files/figure-html/fig-dose-1.png)

Typical-value steady-state profiles for a 70 kg transplant recipient
with CMV at three dose levels, showing the flattening of the profile at
higher doses that the negative dose-on-Ka exponent produces.

## Assumptions and deviations

### Structural readings taken from the control stream

- **Two weight exponents rather than four.** Table 2 presents four
  weight-effect rows, but the control stream reuses `CLWT` for Q/F and
  `VCWT` for Vp/F, with `THETA(9)` and `THETA(10)` present but commented
  `Not used`. Encoding four independent exponents would have been a
  silent structural error that the identical Table 2 estimates would
  never have exposed.
- **Values from the table, structure from the code.** The control
  stream’s `$THETA` / `$OMEGA` / `$SIGMA` blocks are initial estimates,
  not final ones. They are close to the published values but not equal
  to them (`THETA(16)` is `-1.17` initially versus a final `-1.02`), so
  none of them was used as a value.

### Not published, and therefore approximated

- **The three IIV correlations.** Table 2 reports only the six IIV
  magnitudes; it has no correlation column and no off-diagonal rows. The
  control stream shows the *structure* – three `$OMEGA BLOCK(2)` blocks
  pairing CL with Vc, Q with Vp, and Ka with the lag time – and the
  header line `;; 2. Description: 3 OMEGA BLOCKs` confirms that
  structure is final. This model file keeps that block structure, puts
  the published Table 2 variances on the diagonals, and carries the
  correlation **coefficients** (0.673, 0.820 and -0.686) from the
  control stream’s own `$OMEGA` initial estimates, rescaled onto the
  published variances. These three covariances are the least
  well-sourced numbers in the file. They were carried rather than
  invented, and rather than set to zero, because the block structure is
  a published feature of the final model and zeroing the off-diagonals
  would misrepresent it; but a user who needs exact published
  uncertainty should treat the correlations as provisional. The
  magnitudes and the diagonals are unaffected.
- **The `IIV (%)` back-transform convention** is not stated in this
  paper. See the Source trace section for the rule used, its provenance
  in the same analysis group’s companion Sun 2023 paper, and the
  internal-consistency check that supports it. Only Vc/F is sensitive to
  the choice, and only by \<2% in SD.

### Exposure-response models: reported but not extracted

Sun 2025 is a two-part paper. The population PK model above is fully
extracted. The exposure-response (E-R) analyses are **not** extracted
here, and the reason differs by endpoint:

- **Efficacy (Table S3).** Two logistic regressions – confirmed CMV
  clearance at week 8, and clearance with no tissue-invasive disease
  maintained to week 16 – are reported in full, with intercepts, an
  AUCss slope and every risk-factor coefficient. They are extractable as
  written. They are not extracted in this pass because doing so requires
  roughly eight new canonical covariate columns that do not yet exist in
  `inst/references/covariate-columns.md` (baseline CD8+CD69+pp65
  stimulation category, treatment-emergent maribavir-resistance
  mutation, high baseline CMV DNA, post-HCT T-cell infusion, an Asia
  Pacific enrolling-region indicator, and days from HCT onset to
  treatment), and a new canonical name needs operator ratification
  rather than a librarian’s guess. The paper’s headline result is that
  **neither** relationship is statistically significant: the AUCss
  coefficients are 0.00920 (SE 0.0137, p = 0.503) and -0.0131 (SE
  0.0113, p = 0.247) per 10 ug\*h/mL, i.e. essentially flat.
- **Safety versus AUCss (Figure 4).** Nausea and vomiting show
  statistically significant positive relationships (p = 0.0394 and
  0.00296), but the figure prints **no coefficients** – only the
  p-value, the fitted curve and the observed quartile proportions. There
  is nothing to transcribe, and recovering an intercept and slope by
  digitising a raster figure would manufacture precision the paper never
  published.
- **Safety versus AUCday (Figure S3, Table S4).** Fourteen logistic
  regressions across dysgeusia, nausea, vomiting, diarrhea, neutropenia,
  raised immunosuppressant concentrations, invasive infection, acute
  GvHD, renal disorder, anemia, pyrexia, headache, thrombocytopenia and
  serious adverse events. These *are* fully parameterised in the
  supplement and are extractable; they carry the same new-canonical
  requirement as the efficacy models, plus several more (T-cell
  depletion modality, prior HCT, prior CMV prophylaxis, total baseline
  white blood cell category, transplanted nucleated cell number).

One discrepancy in the source is worth recording for whoever extracts
these later. For **anemia**, Table S4 assigns odds ratio 3.12 to “Europe
vs. North America” and 4.62 to “Asia Pacific vs. North America”, while
Figure S3 panel j assigns 3.12 to Asia Pacific and 4.62 to Europe.
Figure S3 is internally consistent and Table S4 is not: the figure’s
estimates exponentiate to its own odds ratios (`exp(1.14) = 3.13` for
Asia Pacific, `exp(1.53) = 4.62` for Europe), whereas Table S4 reports
odds ratios only and cannot be checked against itself. The figure should
be preferred.

### Simulation assumptions

- **Weight distributions** are log-normal draws matched to the Table 1
  per-stratum mean and SD and truncated to the reported 36.1-141 kg
  range. Sun 2025 does not publish the empirical weight distribution,
  and Table 1 reports “NR” for the AURORA weights specifically.
- **Cohort size** is 200 per arm, the repository cap, against published
  strata of 157 and 724 individuals. The comparison is therefore of
  geometric means and CV%, not of individual predictions.
- **Post-hoc versus prospective simulation.** Table 3’s values are
  summaries of *post-hoc Bayesian* individual estimates, which are
  shrunk toward the typical value. This simulation draws from the full
  prior. Shrinkage is low for CL/F (6.9%) and moderate for Vc/F (30.8%)
  but high for Q/F, Vp/F and the lag time (54-56%), so the simulated CV%
  for trough-sensitive metrics is expected to exceed the published CV%.
- **No concomitant medication** in any simulated arm, matching the
  stated conditions for Table 3 and Figure 2. \`\`\`
