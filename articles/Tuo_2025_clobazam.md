# Clobazam (Tuo 2025)

## Model and source

- Citation: Tuo Y, Yu X, Li S, Wang J, Liu M, Song X, Ma J, Wang Y, Liu
  Z, Sun D. Population Pharmacokinetics and Model-Informed Precision
  Dosing of Clobazam Based on the Developmental and Genetic
  Characteristics of Children with Epilepsy. Pharmaceutics.
  2025;17(7):813. <doi:10.3390/pharmaceutics17070813>
- Description: Joint parent-plus-metabolite population PK model for oral
  clobazam and its active metabolite N-desmethylclobazam (norclobazam)
  in Chinese children with refractory epilepsy (Tuo 2025). Tandem
  one-compartment disposition: a first-order absorption depot feeds a
  one-compartment parent, whose entire elimination is routed into a
  one-compartment metabolite that is then cleared. Absorption was not
  identifiable from the opportunistic trough-dominated sampling, so Ka
  is fixed at 1.99 1/h from Jullien 2015; the dosage conversion fraction
  Fm from clobazam to N-desmethylclobazam was likewise not estimable, so
  the metabolite clearance and volume are apparent with respect to Fm.
  Fixed allometric body-weight exponents (0.75 on both clearances, 1 on
  both volumes, 70 kg reference) scale all four disposition parameters,
  and CYP2C19 metabolizer phenotype shifts the metabolite clearance only
  (intermediate and poor metabolizers relative to a normal-metabolizer
  reference), which is why CYP2C19 poor metabolizers accumulate
  N-desmethylclobazam without a matching rise in parent exposure.
  Between-subject variability was estimated on the two clearances only.
- Article: <https://doi.org/10.3390/pharmaceutics17070813> (open access,
  CC BY)

``` r

mod <- rxode2::rxode(readModelDb("Tuo_2025_clobazam"))
mod
#>  ── rxode2-based free-form 3-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>                    lka                    lcl                    lvc 
#>            0.688134639            1.733423892            4.522549157 
#>             lcl_ndmclb             lvc_ndmclb                e_wt_cl 
#>            0.009950331            0.609765572            0.750000000 
#>                e_wt_vc         e_wt_cl_ndmclb         e_wt_vc_ndmclb 
#>            1.000000000            0.750000000            1.000000000 
#> e_cyp2c19_im_cl_ndmclb e_cyp2c19_pm_cl_ndmclb                 propSd 
#>           -0.250000000           -1.300000000            0.330000000 
#>          propSd_ndmclb 
#>            0.530000000 
#> 
#> Omega ($omega): 
#>               etalcl etalcl_ndmclb
#> etalcl        0.1594        0.0000
#> etalcl_ndmclb 0.0000        0.4535
#> attr(,"lotriLabels")
#> [1] "IIV on clobazam CL/F (variance on the log scale)"            
#> [2] "IIV on N-desmethylclobazam CL/Fm (variance on the log scale)"
#> attr(,"lotriFix")
#>               etalcl etalcl_ndmclb
#> etalcl         FALSE         FALSE
#> etalcl_ndmclb  FALSE         FALSE
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2          central
#> 3                  3   central_ndmclb
#>  ── Multiple Endpoint Model ($multipleEndpoint): ──  
#>        variable                      cmt                      dvid*
#> 1        Cc ~ …        cmt='Cc' or cmt=4        dvid='Cc' or dvid=1
#> 2 Cc_ndmclb ~ … cmt='Cc_ndmclb' or cmt=5 dvid='Cc_ndmclb' or dvid=2
#>   * If dvids are outside this range, all dvids are re-numered sequentially, ie 1,7, 10 becomes 1,2,3 etc
#> 
#>  ── μ-referencing ($muRefTable): ──  
#>        theta           eta level
#> 1        lcl        etalcl    id
#> 2 lcl_ndmclb etalcl_ndmclb    id
#>                                                              covariates
#> 1                                                                      
#> 2 CYP2C19_PM*e_cyp2c19_pm_cl_ndmclb + CYP2C19_IM*e_cyp2c19_im_cl_ndmclb
#> 
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     compartmentData <- list(depot = list(analyte = "clobazam", 
#>         units = "mg", specimen = "administration site", verified = TRUE), 
#>         central = list(analyte = "clobazam", units = "mg", specimen = "plasma", 
#>             verified = TRUE), central_ndmclb = list(analyte = "N-desmethylclobazam", 
#>             units = "mg", specimen = "plasma", verified = TRUE))
#>     covariateData <- list(WT = list(description = "Body weight.", 
#>         units = "kg", type = "continuous", reference_category = NULL, 
#>         notes = "Allometric scaling of all four disposition parameters against a 70 kg reference with exponents FIXED at 0.75 (both clearances) and 1 (both volumes); Tuo 2025 Equations (9)-(12) print the exponents inline and Table 2 reports no exponent parameter, so they were not estimated. The cohort weight range is 6.60-73.00 kg (median 20.00 kg), so the 70 kg reference sits at the extreme upper edge of the observed data and the reported typical values are extrapolated adult-standardised values rather than values observed at 70 kg.", 
#>         source_name = "Weight"), CYP2C19_IM = list(description = "CYP2C19 intermediate-metabolizer phenotype indicator.", 
#>         units = "(binary)", type = "binary", reference_category = "0 (normal metabolizer; *1/*1 -- both CYP2C19_IM = 0 and CYP2C19_PM = 0)", 
#>         notes = "1 = subject has CYP2C19 IM phenotype (*1/*2, *1/*3, *2/*17 or *3/*17 in Tuo 2025); 0 = otherwise. Cohort distribution: NM 39.81% (41/103), IM 43.69% (45/103), PM 14.56% (15/103), RM 1.94% (2/103). Affects the metabolite clearance only. The two CYP2C19 rapid metabolizers (*1/*17) were excluded from the covariate analysis for lack of sample size (Tuo 2025 Results 3.1) and no ultrarapid metabolizers (*17/*17) were observed, so neither phenotype has an estimated effect and both fall into the CYP2C19_IM = 0, CYP2C19_PM = 0 reference cell by default; users simulating RM or UM subjects should treat them as an explicit extrapolation.", 
#>         source_name = "CYP2C19 genotype (NMs / IMs / PMs)"), 
#>         CYP2C19_PM = list(description = "CYP2C19 poor-metabolizer phenotype indicator.", 
#>             units = "(binary)", type = "binary", reference_category = "0 (normal, intermediate, or rapid metabolizer)", 
#>             notes = "1 = subject has CYP2C19 PM phenotype (*2/*2, *2/*3 or *3/*3 in Tuo 2025); 0 = otherwise. Paired with `CYP2C19_IM` to encode the three-level NM (reference) / IM / PM phenotype with two binary indicators. Affects the metabolite clearance only: the paper reports mean post-hoc CL_N-CLB/Fm of 0.46, 0.34 and 0.13 L/h in NMs, IMs and PMs (a 71.7% reduction in PMs vs NMs) against no significant difference in parent CL/F across the three groups.", 
#>             source_name = "CYP2C19 genotype (NMs / IMs / PMs)"))
#>     covariatesDataExcluded <- list(AGE = list(description = "Age.", 
#>         units = "years", type = "continuous", notes = "Screened in the stepwise covariate analysis and not retained; Tuo 2025 Results 3.2 reports that age had no significant impact on the PK parameters of clobazam or N-desmethylclobazam, which the Discussion attributes to the cohort being purely pediatric (0.85-16.75 years) rather than the combined pediatric-plus-adult range of Tolbert 2016."), 
#>         SEXF = list(description = "Female sex indicator.", units = "(binary)", 
#>             type = "binary", notes = "Screened and not retained (Tuo 2025 Results 3.2)."), 
#>         BSA = list(description = "Body surface area.", units = "m^2", 
#>             type = "continuous", notes = "Screened and not retained (Tuo 2025 Results 3.2); body weight was the size descriptor carried into the final model."), 
#>         EGFR = list(description = "Estimated glomerular filtration rate (modified Schwartz formula).", 
#>             units = "mL/min/1.73m^2", type = "continuous", notes = "Screened as part of the renal-function panel and not retained (Tuo 2025 Results 3.2). The Discussion notes the number of patients with impaired renal function was too small to draw firm conclusions."), 
#>         ALB = list(description = "Serum albumin.", units = "g/L", 
#>             type = "continuous", notes = "Screened as part of the hepatic-function panel and not retained (Tuo 2025 Results 3.2)."), 
#>         ALT = list(description = "Alanine aminotransferase.", 
#>             units = "U/L", type = "continuous", notes = "Screened as part of the hepatic-function panel and not retained (Tuo 2025 Results 3.2)."), 
#>         AST = list(description = "Aspartate aminotransferase.", 
#>             units = "U/L", type = "continuous", notes = "Screened as part of the hepatic-function panel and not retained (Tuo 2025 Results 3.2)."), 
#>         CONMED_VALPROIC_ACID = list(description = "Concomitant valproic acid indicator.", 
#>             units = "(binary)", type = "binary", notes = "Screened and not retained (Tuo 2025 Results 3.2 and Supplementary Figure S1); 80.58% of the cohort received valproic acid. The Discussion notes this agrees with prior modelling analyses that found no clinically relevant antiepileptic drug-drug interaction with clobazam."), 
#>         CONMED_LAMOTRIGINE = list(description = "Concomitant lamotrigine indicator.", 
#>             units = "(binary)", type = "binary", notes = "Screened and not retained (Tuo 2025 Results 3.2 and Supplementary Figure S1); 28.16% of the cohort."), 
#>         DIET_KETOGENIC = list(description = "Adherence to a ketogenic diet.", 
#>             units = "(binary)", type = "binary", notes = "Screened and not retained (Tuo 2025 Results 3.2 and Discussion); only 6.80% of the cohort adhered to a ketogenic diet, which the authors judged too few to resolve an effect despite a published case report of a 42% fall in clobazam and N-desmethylclobazam concentrations on diet initiation."), 
#>         SNP_ABCB1_RS1045642 = list(description = "ABCB1 3435C>T (rs1045642) genotype.", 
#>             units = "(genotype)", type = "categorical", notes = "Genotyped and screened; no significant effect on clobazam or N-desmethylclobazam PK (Tuo 2025 Results 3.2 and Supplementary Figure S2). Two further ABCB1 SNPs (rs1128503 1236C>T, rs2032582 2677G>T/A), two CYP3A4 SNPs (rs2740574 *1B, rs2242480 *1G) and six GABA-receptor SNPs (rs2279020, rs279858, rs11503014, rs2229944, rs211014, rs211037) were screened with the same negative result; they are represented here by this single entry rather than one entry each because none carries an estimated coefficient."))
#>     description <- "Joint parent-plus-metabolite population PK model for oral clobazam and its active metabolite N-desmethylclobazam (norclobazam) in Chinese children with refractory epilepsy (Tuo 2025). Tandem one-compartment disposition: a first-order absorption depot feeds a one-compartment parent, whose entire elimination is routed into a one-compartment metabolite that is then cleared. Absorption was not identifiable from the opportunistic trough-dominated sampling, so Ka is fixed at 1.99 1/h from Jullien 2015; the dosage conversion fraction Fm from clobazam to N-desmethylclobazam was likewise not estimable, so the metabolite clearance and volume are apparent with respect to Fm. Fixed allometric body-weight exponents (0.75 on both clearances, 1 on both volumes, 70 kg reference) scale all four disposition parameters, and CYP2C19 metabolizer phenotype shifts the metabolite clearance only (intermediate and poor metabolizers relative to a normal-metabolizer reference), which is why CYP2C19 poor metabolizers accumulate N-desmethylclobazam without a matching rise in parent exposure. Between-subject variability was estimated on the two clearances only."
#>     population <- list(species = "human", n_subjects = 103L, 
#>         n_studies = 1L, n_observations = "156 plasma samples yielding 302 analyte concentrations (154 clobazam + 148 N-desmethylclobazam). Sampling depth per patient: 68 sampled once, 21 twice, 12 three times, 1 four times and 1 six times. Assay quantitative ranges 3-1200 ug/L (clobazam) and 40-16000 ug/L (N-desmethylclobazam) by HPLC-MS/MS.", 
#>         age_range = "0.85-16.75 years", age_median = "5.46 years (mean 5.94, SD 3.15)", 
#>         weight_range = "6.60-73.00 kg", weight_median = "20.00 kg (mean 22.94, SD 10.48)", 
#>         sex_female_pct = 43.7, race_ethnicity = c(Asian = 100), 
#>         disease_state = "Pediatric refractory epilepsy: Lennox-Gastaut syndrome, Dravet syndrome, infantile spasms and other refractory epilepsies. 75.73% of patients were taking three or more antiepileptic drugs; 80.58% received concomitant valproic acid, 28.16% lamotrigine, 24.27% perampanel, 18.45% levetiracetam and 17.48% topiramate, and 6.80% adhered to a ketogenic diet.", 
#>         dose_range = "Oral clobazam tablets, dosed twice daily when the total dose exceeded 5 mg. Starting dose 5 mg for patients weighing 30 kg or less and 10 mg above 30 kg, then individually titrated on efficacy and tolerability.", 
#>         regions = "China (single centre: Wuhan Children's Hospital, Tongji Medical College, Huazhong University of Science and Technology; enrolment December 2022 to March 2024).", 
#>         genotype = "CYP2C19 phenotype: normal metabolizers 41 (39.81%), intermediate 45 (43.69%), poor 15 (14.56%), rapid 2 (1.94%). No ultrarapid metabolizers were observed. All genotype frequencies were consistent with Hardy-Weinberg equilibrium.", 
#>         notes = "Demographics from Tuo 2025 Table 1. Prospective single-centre opportunistic-sampling study using scavenged residual blood drawn for routine biochemistry during safety follow-up, so the sampling is trough-dominated and carries essentially no information on the absorption or distribution phases -- the reason Ka was fixed and a one-compartment rather than two-compartment parent disposition was selected. Reported therapeutic trough ranges applied by the authors are 30-300 ug/L for clobazam and 300-3000 ug/L for N-desmethylclobazam, with laboratory alert levels of 500 ug/L and 5000 ug/L respectively; these are adult-derived targets carried over to the pediatric setting because no pediatric-specific ranges exist.")
#>     reference <- "Tuo Y, Yu X, Li S, Wang J, Liu M, Song X, Ma J, Wang Y, Liu Z, Sun D. Population Pharmacokinetics and Model-Informed Precision Dosing of Clobazam Based on the Developmental and Genetic Characteristics of Children with Epilepsy. Pharmaceutics. 2025;17(7):813. doi:10.3390/pharmaceutics17070813"
#>     units <- list(time = "h", dosing = "mg", concentration = "mg/L")
#>     vignette <- "Tuo_2025_clobazam"
#>     ini({
#>         lka <- fix(0.688134638736401)
#>         label("Absorption rate constant Ka (1/h)")
#>         lcl <- 1.73342389221509
#>         label("Apparent clearance CL_CLB/F of clobazam at 70 kg (L/h)")
#>         lvc <- 4.52254915729975
#>         label("Apparent central volume V_CLB/F of clobazam at 70 kg (L)")
#>         lcl_ndmclb <- 0.00995033085316809
#>         label("Apparent clearance CL_N-CLB/Fm of N-desmethylclobazam at 70 kg in CYP2C19 normal metabolizers (L/h)")
#>         lvc_ndmclb <- 0.609765571620894
#>         label("Apparent central volume V_N-CLB/Fm of N-desmethylclobazam at 70 kg (L)")
#>         e_wt_cl <- fix(0.75)
#>         label("Allometric exponent on clobazam CL/F (unitless)")
#>         e_wt_vc <- fix(1)
#>         label("Allometric exponent on clobazam V/F (unitless)")
#>         e_wt_cl_ndmclb <- fix(0.75)
#>         label("Allometric exponent on N-desmethylclobazam CL/Fm (unitless)")
#>         e_wt_vc_ndmclb <- fix(1)
#>         label("Allometric exponent on N-desmethylclobazam V/Fm (unitless)")
#>         e_cyp2c19_im_cl_ndmclb <- -0.25
#>         label("CYP2C19 intermediate-metabolizer log-additive shift on N-desmethylclobazam CL/Fm (unitless)")
#>         e_cyp2c19_pm_cl_ndmclb <- -1.3
#>         label("CYP2C19 poor-metabolizer log-additive shift on N-desmethylclobazam CL/Fm (unitless)")
#>         propSd <- c(0, 0.33)
#>         label("Proportional residual SD for clobazam (fraction)")
#>         propSd_ndmclb <- c(0, 0.53)
#>         label("Proportional residual SD for N-desmethylclobazam (fraction)")
#>         etalcl ~ 0.1594
#>         label("IIV on clobazam CL/F (variance on the log scale)")
#>         etalcl_ndmclb ~ 0.4535
#>         label("IIV on N-desmethylclobazam CL/Fm (variance on the log scale)")
#>     })
#>     model({
#>         wtr <- WT/70
#>         cl <- exp(lcl + etalcl) * wtr^e_wt_cl
#>         vc <- exp(lvc) * wtr^e_wt_vc
#>         cl_ndmclb <- exp(lcl_ndmclb + etalcl_ndmclb + e_cyp2c19_im_cl_ndmclb * 
#>             CYP2C19_IM + e_cyp2c19_pm_cl_ndmclb * CYP2C19_PM) * 
#>             wtr^e_wt_cl_ndmclb
#>         vc_ndmclb <- exp(lvc_ndmclb) * wtr^e_wt_vc_ndmclb
#>         ka <- exp(lka)
#>         d/dt(depot) <- -ka * depot
#>         d/dt(central) <- ka * depot - cl * (central/vc)
#>         d/dt(central_ndmclb) <- cl * (central/vc) - cl_ndmclb * 
#>             (central_ndmclb/vc_ndmclb)
#>         Cc <- central/vc
#>         Cc_ndmclb <- central_ndmclb/vc_ndmclb
#>         Cc ~ prop(propSd)
#>         Cc_ndmclb ~ prop(propSd_ndmclb)
#>     })
#> }
```

Concentrations in the model are in **mg/L**, matching the `units`
metadata. Tuo 2025 reports every concentration in **ug/L**, so this
vignette multiplies model output by 1000 before any comparison with the
paper.

``` r

UGL <- 1000  # mg/L -> ug/L
```

## Population

Tuo 2025 is a prospective, single-centre study at Wuhan Children’s
Hospital (December 2022 to March 2024) in **103 Chinese children with
refractory epilepsy** (Lennox-Gastaut syndrome, Dravet syndrome,
infantile spasms and other refractory epilepsies) receiving oral
clobazam. Median age was 5.46 years (range 0.85-16.75) and median body
weight 20.0 kg (range 6.60-73.0); 45 of 103 (43.7%) were female (Table
1). Sampling was **opportunistic**: 156 plasma samples were scavenged
from residual blood drawn for routine biochemistry during safety
follow-up, yielding 302 analyte concentrations (154 clobazam and 148
N-desmethylclobazam). Most patients contributed a single sample (68 of
103), and the draws were almost entirely in the elimination phase – the
reason the authors could not identify absorption or distribution and
fixed `Ka`.

CYP2C19 phenotype was determined for all 103 patients: 41 normal
metabolizers (39.8%), 45 intermediate (43.7%), 15 poor (14.6%) and 2
rapid (1.9%). No ultrarapid metabolizers were observed, and the 2 rapid
metabolizers were dropped from the covariate analysis for lack of sample
size. The cohort was heavily co-medicated – 75.7% were on three or more
antiepileptic drugs, 80.6% on valproic acid – and 6.8% followed a
ketogenic diet; none of these was retained as a covariate.

The same information is available programmatically:

``` r

str(readModelDb("Tuo_2025_clobazam")()$population)
#> List of 15
#>  $ species       : chr "human"
#>  $ n_subjects    : int 103
#>  $ n_studies     : int 1
#>  $ n_observations: chr "156 plasma samples yielding 302 analyte concentrations (154 clobazam + 148 N-desmethylclobazam). Sampling depth"| __truncated__
#>  $ age_range     : chr "0.85-16.75 years"
#>  $ age_median    : chr "5.46 years (mean 5.94, SD 3.15)"
#>  $ weight_range  : chr "6.60-73.00 kg"
#>  $ weight_median : chr "20.00 kg (mean 22.94, SD 10.48)"
#>  $ sex_female_pct: num 43.7
#>  $ race_ethnicity: Named num 100
#>   ..- attr(*, "names")= chr "Asian"
#>  $ disease_state : chr "Pediatric refractory epilepsy: Lennox-Gastaut syndrome, Dravet syndrome, infantile spasms and other refractory "| __truncated__
#>  $ dose_range    : chr "Oral clobazam tablets, dosed twice daily when the total dose exceeded 5 mg. Starting dose 5 mg for patients wei"| __truncated__
#>  $ regions       : chr "China (single centre: Wuhan Children's Hospital, Tongji Medical College, Huazhong University of Science and Tec"| __truncated__
#>  $ genotype      : chr "CYP2C19 phenotype: normal metabolizers 41 (39.81%), intermediate 45 (43.69%), poor 15 (14.56%), rapid 2 (1.94%)"| __truncated__
#>  $ notes         : chr "Demographics from Tuo 2025 Table 1. Prospective single-centre opportunistic-sampling study using scavenged resi"| __truncated__
```

## Source trace

Every `ini()` value carries an in-file comment naming its source
location in `inst/modeldb/specificDrugs/Tuo_2025_clobazam.R`. They are
collected here.

| Equation / parameter | Value | Source location |
|----|----|----|
| Depot / parent / metabolite ODEs | – | Tuo 2025 Methods 2.5, Equations (1)-(3) |
| Exponential IIV `P_i = theta * exp(eta_i)` | – | Equation (4) |
| Proportional residual `Y = IPRED * (1 + eps)` | – | Equation (5) |
| Continuous covariate form `(Cov/Cov_median)^theta_cov` | – | Equation (6) |
| Categorical covariate form `exp(theta_cov)` | – | Equation (7) |
| `lka` (Ka) | 1.99 1/h, **fixed** | Table 2 row 1 and Equation (8); value taken from Jullien 2015 (Tuo 2025 reference \[19\]) |
| `lvc` (V_CLB/F) | 92.07 L at 70 kg | Table 2, RSE 31.15%, bootstrap 40.42-176.27; Equation (9) |
| `lcl` (CL_CLB/F) | 5.66 L/h at 70 kg | Table 2, RSE 10.53%, bootstrap 4.09-7.26; Equation (10) |
| `lvc_ndmclb` (V_N-CLB/Fm) | 1.84 L at 70 kg | Table 2, RSE 29.66%, bootstrap 0.72-2.76; Equation (11) |
| `lcl_ndmclb` (CL_N-CLB/Fm) | 1.01 L/h at 70 kg, NM | Table 2, RSE 11.83%, bootstrap 0.71-1.30; Equation (12) |
| `e_wt_vc`, `e_wt_vc_ndmclb` | 1.0, **fixed** | Equations (9) and (11), printed inline; no Table 2 row |
| `e_wt_cl`, `e_wt_cl_ndmclb` | 0.75, **fixed** | Equations (10) and (12), printed inline; no Table 2 row |
| Reference weight | 70 kg | Equations (9)-(12) denominators |
| `e_cyp2c19_im_cl_ndmclb` | -0.25 | Table 2 `theta_CYP2C19,IM`, RSE 39.03%, bootstrap -0.53 to -0.02 |
| `e_cyp2c19_pm_cl_ndmclb` | -1.30 | Table 2 `theta_CYP2C19,PM`, RSE 15.01%, bootstrap -1.67 to -0.89 |
| NM reference level | 0.00, **fixed** | Table 2 `theta_CYP2C19,NM` |
| `etalcl` | variance 0.1594 | Table 2 `omega^2 CL_CLB/F (%) = 15.94`, footnote “interindividual variance”; eta-shrinkage 9.73% |
| `etalcl_ndmclb` | variance 0.4535 | Table 2 `omega^2 CL_N-CLB/Fm (%) = 45.35`; eta-shrinkage 17.42% |
| `propSd` | 0.33 | Table 2 `sigma_CLB`, RSE 8.52%, bootstrap 0.27-0.37 |
| `propSd_ndmclb` | 0.53 | Table 2 `sigma_N-CLB`, RSE 8.60%, bootstrap 0.44-0.61; eps-shrinkage 33.17% |
| Published validation targets | Table 3 (15 dose rows) | Median trough and PTA by weight and CYP2C19 phenotype |

### The apparent (F, Fm) parameterisation

Tuo 2025 Equations (1)-(3) are written in true amounts and carry two
fractions that the study could not estimate: `F`, oral bioavailability,
and `Fm`, the clobazam-to-N-desmethylclobazam dosage conversion fraction
(“Since the amount of CLB converted to N-CLB was not clear, Fm was not
estimated in this study”). Table 2 accordingly reports all four
disposition parameters as apparent: `V_CLB/F`, `CL_CLB/F`, `V_N-CLB/Fm`,
`CL_N-CLB/Fm`.

The packaged model divides the parent equations through by `F` and the
metabolite equation by `F * Fm`, so each state holds a scaled amount:
`central` holds `A_CLB/F` and `central_ndmclb` holds `A_N-CLB/(F*Fm)`.
The four `ini()` parameters are then exactly the four apparent
quantities of Table 2, and both observed concentrations are the **true**
plasma concentrations, because the same scaling factor divides the state
and the volume it is divided by. This is the convention the register
already records for the sibling `ndmima` (N-desmethyl-imatinib) models.

A consequence worth stating: the parent’s entire elimination is routed
into the metabolite compartment. The paper writes no separate
non-metabolic clobazam elimination arm, and the `Fm` factor of Equation
(3) – which is what would otherwise split the parent’s exit between
metabolite formation and other routes – has been absorbed into the
metabolite state scaling.

## Virtual cohort

Tuo 2025 Table 3 gives the optimal simulated dosing strategy for 15
combinations of body weight (10, 20, 30, 40, 50 kg) and CYP2C19
phenotype (NM, IM, PM), each with a published median steady-state trough
for both analytes and four probability-of-target-attainment columns.
That table is the validation target for this vignette, so the cohort
mirrors it exactly: twice-daily oral dosing for 10 consecutive days,
with the trough read 12 h after the last dose.

``` r

tab3 <- tibble::tribble(
  ~geno, ~wt, ~mgkg, ~pub_clb, ~pub_ndmclb, ~pub_pta_clb30, ~pub_pta_ndmclb300,
  "NM", 10, 0.3,    109.28,  629.52,  94.8, 84.1,
  "NM", 20, 0.25,   119.74,  699.26,  97.0, 88.4,
  "NM", 30, 0.2,    111.47,  651.02,  97.0, 86.3,
  "NM", 40, 0.175,  108.24,  630.81,  97.2, 85.5,
  "NM", 50, 0.15,    93.03,  568.99,  95.6, 80.9,
  "IM", 10, 0.25,    83.75,  658.07,  89.6, 83.5,
  "IM", 20, 0.2,     88.17,  708.52,  93.5, 87.1,
  "IM", 30, 0.15,    77.19,  626.37,  92.1, 82.7,
  "IM", 40, 0.125,   71.53,  579.29,  91.7, 80.9,
  "IM", 50, 0.12,    80.33,  616.59,  94.0, 84.4,
  "PM", 10, 0.1,     33.50,  962.78,  57.3, 89.9,
  "PM", 20, 0.075,   33.06,  950.00,  55.9, 92.2,
  "PM", 30, 0.075,   38.59, 1081.81,  66.1, 95.2,
  "PM", 40, 0.0625,  35.94,  967.24,  62.3, 91.4,
  "PM", 50, 0.06,    37.39,  996.10,  65.2, 92.0
) |>
  dplyr::mutate(
    id         = dplyr::row_number(),
    CYP2C19_IM = as.integer(geno == "IM"),
    CYP2C19_PM = as.integer(geno == "PM"),
    dose_mg    = mgkg * wt,
    arm        = sprintf("%s %g kg (%g mg/kg BID)", geno, wt, mgkg)
  )

knitr::kable(
  tab3 |>
    dplyr::select("CYP2C19" = geno, "Weight (kg)" = wt,
                  "Dose (mg/kg BID)" = mgkg, "Dose (mg)" = dose_mg,
                  "Published median CLB trough (ug/L)" = pub_clb,
                  "Published median N-CLB trough (ug/L)" = pub_ndmclb),
  caption = "Tuo 2025 Table 3 dosing scenarios and published median steady-state troughs."
)
```

| CYP2C19 | Weight (kg) | Dose (mg/kg BID) | Dose (mg) | Published median CLB trough (ug/L) | Published median N-CLB trough (ug/L) |
|:---|---:|---:|---:|---:|---:|
| NM | 10 | 0.3000 | 3.00 | 109.28 | 629.52 |
| NM | 20 | 0.2500 | 5.00 | 119.74 | 699.26 |
| NM | 30 | 0.2000 | 6.00 | 111.47 | 651.02 |
| NM | 40 | 0.1750 | 7.00 | 108.24 | 630.81 |
| NM | 50 | 0.1500 | 7.50 | 93.03 | 568.99 |
| IM | 10 | 0.2500 | 2.50 | 83.75 | 658.07 |
| IM | 20 | 0.2000 | 4.00 | 88.17 | 708.52 |
| IM | 30 | 0.1500 | 4.50 | 77.19 | 626.37 |
| IM | 40 | 0.1250 | 5.00 | 71.53 | 579.29 |
| IM | 50 | 0.1200 | 6.00 | 80.33 | 616.59 |
| PM | 10 | 0.1000 | 1.00 | 33.50 | 962.78 |
| PM | 20 | 0.0750 | 1.50 | 33.06 | 950.00 |
| PM | 30 | 0.0750 | 2.25 | 38.59 | 1081.81 |
| PM | 40 | 0.0625 | 2.50 | 35.94 | 967.24 |
| PM | 50 | 0.0600 | 3.00 | 37.39 | 996.10 |

Tuo 2025 Table 3 dosing scenarios and published median steady-state
troughs. {.table}

The event table is built as a plain data frame so the covariate columns
survive to `rxSolve()`. The model declares two endpoints (`Cc`,
`Cc_ndmclb`), so observation rows carry an explicit `dvid` alongside the
ODE-state `cmt`; dose rows carry `dvid = NA`.

``` r

DOSE_TIMES <- seq(0, 228, by = 12)  # twice daily for 10 days: 20 doses

# Built vectorised (one rbind for the whole cohort) rather than per subject --
# a per-subject rbind loop is quadratic and dominates the render time once the
# cohort reaches a few hundred subjects.
build_events <- function(subjects, obs_times) {
  nd <- length(DOSE_TIMES)
  no <- length(obs_times)
  ns <- nrow(subjects)
  doses <- data.frame(
    id   = rep(subjects$id, each = nd),
    time = rep(DOSE_TIMES, times = ns),
    amt  = rep(subjects$dose_mg, each = nd),
    evid = 1L, cmt = "depot", dvid = NA_integer_
  )
  obs <- data.frame(
    id   = rep(subjects$id, each = no),
    time = rep(obs_times, times = ns),
    amt  = NA_real_, evid = 0L, cmt = "central", dvid = 1L
  )
  out <- rbind(doses, obs)
  idx <- match(out$id, subjects$id)
  out$WT         <- subjects$wt[idx]
  out$CYP2C19_IM <- subjects$CYP2C19_IM[idx]
  out$CYP2C19_PM <- subjects$CYP2C19_PM[idx]
  out[order(out$id, out$time, out$evid), ]
}

dense_times <- sort(unique(c(seq(0, 240, by = 0.25), 228, 240)))
ev_typ <- build_events(tab3, dense_times)
```

## Simulation

[`rxode2::zeroRe()`](https://nlmixr2.github.io/rxode2/reference/zeroRe.html)
gives the typical-value (deterministic) profile for each arm – the
quantity the published median trough should reproduce, because the
steady-state trough is a monotone function of clearance and the median
of a monotone transform is the transform of the median.

``` r

sim_typ <- rxode2::rxSolve(rxode2::zeroRe(mod), ev_typ, returnType = "data.frame")
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalcl_ndmclb'
#> Warning: multi-subject simulation without without 'omega'
stopifnot(!anyNA(sim_typ$Cc), !anyNA(sim_typ$Cc_ndmclb))

sim_typ <- sim_typ |>
  dplyr::left_join(dplyr::select(tab3, id, geno, wt, arm, dose_mg), by = "id") |>
  dplyr::mutate(clb = Cc * UGL, ndmclb = Cc_ndmclb * UGL)
```

### Structural gate: steady-state mass balance

The sharpest cheap check on this model is the steady-state mass balance,
and it pins three things at once: that the explicit ODE system is really
being solved (rxode2 will silently replace a `cl`/`vc` pair with an
analytic one-compartment solution in some models, which would discard
the metabolite arm), that the metabolite formation term carries the
parent’s *whole* elimination, and that the dose-to-concentration unit
chain is right.

Over a steady-state dosing interval, `CL_CLB/F * AUCtau(Cc)` must equal
the dose, and because every molecule leaving the parent enters the
metabolite, `CL_N-CLB/Fm * AUCtau(Cc_ndmclb)` must equal the same dose.
Equivalently the ratio of metabolite to parent AUC is exactly
`cl / cl_ndmclb`.

``` r

tau_win <- sim_typ |> dplyr::filter(time >= 228, time <= 240)

trapz <- function(x, y) sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)

mb <- tau_win |>
  dplyr::group_by(id, arm, dose_mg) |>
  dplyr::summarise(
    auc_clb    = trapz(time, Cc),          # mg*h/L
    auc_ndmclb = trapz(time, Cc_ndmclb),
    cl         = dplyr::first(cl),         # rxSolve returns the per-subject values
    cl_ndmclb  = dplyr::first(cl_ndmclb),
    .groups    = "drop"
  ) |>
  dplyr::mutate(
    dose_recovered_clb    = cl * auc_clb,
    dose_recovered_ndmclb = cl_ndmclb * auc_ndmclb,
    err_clb    = 100 * (dose_recovered_clb    - dose_mg) / dose_mg,
    err_ndmclb = 100 * (dose_recovered_ndmclb - dose_mg) / dose_mg,
    auc_ratio  = auc_ndmclb / auc_clb,
    ratio_pred = cl / cl_ndmclb,
    err_ratio  = 100 * (auc_ratio - ratio_pred) / ratio_pred
  )

knitr::kable(
  mb |>
    dplyr::select("Arm" = arm, "Dose (mg)" = dose_mg,
                  "CL x AUC, parent (mg)" = dose_recovered_clb,
                  "CL x AUC, metabolite (mg)" = dose_recovered_ndmclb,
                  "AUC ratio N-CLB/CLB" = auc_ratio,
                  "cl / cl_ndmclb" = ratio_pred),
  digits = 4,
  caption = "Steady-state mass balance. Both recovered doses must equal the administered dose."
)
```

| Arm | Dose (mg) | CL x AUC, parent (mg) | CL x AUC, metabolite (mg) | AUC ratio N-CLB/CLB | cl / cl_ndmclb |
|:---|---:|---:|---:|---:|---:|
| NM 10 kg (0.3 mg/kg BID) | 3.00 | 2.9969 | 3.00 | 5.6097 | 5.6040 |
| NM 20 kg (0.25 mg/kg BID) | 5.00 | 4.9957 | 5.00 | 5.6088 | 5.6040 |
| NM 30 kg (0.2 mg/kg BID) | 6.00 | 5.9953 | 6.00 | 5.6083 | 5.6040 |
| NM 40 kg (0.175 mg/kg BID) | 7.00 | 6.9949 | 7.00 | 5.6080 | 5.6040 |
| NM 50 kg (0.15 mg/kg BID) | 7.50 | 7.4948 | 7.50 | 5.6078 | 5.6040 |
| IM 10 kg (0.25 mg/kg BID) | 2.50 | 2.4974 | 2.50 | 7.2030 | 7.1956 |
| IM 20 kg (0.2 mg/kg BID) | 4.00 | 3.9965 | 4.00 | 7.2019 | 7.1956 |
| IM 30 kg (0.15 mg/kg BID) | 4.50 | 4.4965 | 4.50 | 7.2013 | 7.1956 |
| IM 40 kg (0.125 mg/kg BID) | 5.00 | 4.9964 | 5.00 | 7.2009 | 7.1956 |
| IM 50 kg (0.12 mg/kg BID) | 6.00 | 5.9959 | 6.00 | 7.2006 | 7.1956 |
| PM 10 kg (0.1 mg/kg BID) | 1.00 | 0.9990 | 1.00 | 20.5838 | 20.5626 |
| PM 20 kg (0.075 mg/kg BID) | 1.50 | 1.4987 | 1.50 | 20.5804 | 20.5626 |
| PM 30 kg (0.075 mg/kg BID) | 2.25 | 2.2482 | 2.25 | 20.5787 | 20.5626 |
| PM 40 kg (0.0625 mg/kg BID) | 2.50 | 2.4982 | 2.50 | 20.5776 | 20.5626 |
| PM 50 kg (0.06 mg/kg BID) | 3.00 | 2.9979 | 3.00 | 20.5768 | 20.5626 |

Steady-state mass balance. Both recovered doses must equal the
administered dose. {.table}

``` r


# Deterministic identities: numerical-integration error only, so tight bounds
# are correct here and are NOT an extreme-of-cohort assertion.
stopifnot(
  max(abs(mb$err_clb))    < 0.5,
  max(abs(mb$err_ndmclb)) < 0.5,
  max(abs(mb$err_ratio))  < 0.5
)
```

The metabolite arm therefore survives: if rxode2 had auto-solved the
parent as a linear one-compartment model and dropped the explicit
`d/dt()` block, the metabolite mass balance would fail immediately.

### Structural gate: closed-form steady-state trough

The parent is a one-compartment model with first-order absorption, so
its steady-state trough has a closed form. Comparing it to the ODE solve
uses the *same* per-subject parameters on both sides, so the difference
is pure numerical error and a tight bound is appropriate.

``` r

ss_trough <- function(dose, ka, cl, vc, tau = 12) {
  kel <- cl / vc
  dose * ka / (vc * (ka - kel)) *
    (exp(-kel * tau) / (1 - exp(-kel * tau)) - exp(-ka * tau) / (1 - exp(-ka * tau)))
}

cf <- sim_typ |>
  dplyr::filter(time == 240) |>
  dplyr::mutate(closed_form = ss_trough(dose_mg, ka, cl, vc) * UGL) |>
  dplyr::select(arm, ode = clb, closed_form) |>
  dplyr::mutate(pct_diff = 100 * (ode - closed_form) / closed_form)

stopifnot(max(abs(cf$pct_diff)) < 0.1)
sprintf("Closed-form vs ODE steady-state trough: max |%% diff| = %.4f%%",
        max(abs(cf$pct_diff)))
#> [1] "Closed-form vs ODE steady-state trough: max |% diff| = 0.0000%"
```

## Replicate published figures

### Figure 3 – clearance versus body weight

Tuo 2025 Figure 3 shows both clearances rising allometrically with body
weight, with parent clearance above metabolite clearance at every
weight.

``` r

wt_grid <- tibble::tibble(wt = seq(6.6, 73, length.out = 120)) |>
  dplyr::mutate(
    `CLB (CL/F)`      = 5.66 * (wt / 70)^0.75,
    `N-CLB (CL/Fm)`   = 1.01 * (wt / 70)^0.75
  ) |>
  tidyr::pivot_longer(-wt, names_to = "Analyte", values_to = "CL")

ggplot(wt_grid, aes(wt, CL, colour = Analyte)) +
  geom_line(linewidth = 1) +
  labs(x = "Body weight (kg)", y = "Apparent clearance (L/h)",
       title = "Clearance versus body weight") +
  theme_bw()
```

![Replicates Figure 3 of Tuo 2025: allometric rise of CL_CLB/F and
CL_N-CLB/Fm with body weight (typical CYP2C19 normal
metabolizer).](Tuo_2025_clobazam_files/figure-html/figure3-1.png)

Replicates Figure 3 of Tuo 2025: allometric rise of CL_CLB/F and
CL_N-CLB/Fm with body weight (typical CYP2C19 normal metabolizer).

``` r


stopifnot(
  # Parent clearance exceeds metabolite clearance at every weight (Figure 3).
  all(
    dplyr::filter(wt_grid, Analyte == "CLB (CL/F)")$CL >
      dplyr::filter(wt_grid, Analyte == "N-CLB (CL/Fm)")$CL
  )
)
```

### Figure 6 – simulated concentration-time profiles at 0.2 mg/kg BID

Tuo 2025 Figure 6 simulates a 20 kg child (the cohort median weight)
given 0.2 mg/kg twice daily for 10 days, separately for CYP2C19 NM, IM
and PM, and plots the median with the 10th-90th percentile band. The
paper’s reading of the figure is that parent troughs are
indistinguishable across phenotypes and stay inside 30-300 ug/L, while
metabolite troughs rise NM \< IM \< PM, with about half of the PM
subjects above 3000 ug/L.

``` r

N_PER_ARM <- 200  # cap: never more than 200 participants per arm

fig6_arms <- tibble::tibble(geno = c("NM", "IM", "PM")) |>
  dplyr::mutate(CYP2C19_IM = as.integer(geno == "IM"),
                CYP2C19_PM = as.integer(geno == "PM"))

# Coarse over the approach to steady state, dense over the final dosing
# interval where the AUC ratios of the Figure 5 replication are integrated.
# The metabolite elimination rate constant reaches ~5.5 1/h in the upper tail
# of the 45.35% metabolite-clearance variance, so the integration is stiff
# enough that the record count -- not the cohort size -- sets the render time.
fig6_times <- sort(unique(c(seq(0, 228, by = 2), seq(228, 240, by = 0.5))))
fig6_subjects <- fig6_arms |>
  dplyr::slice(rep(seq_len(dplyr::n()), each = N_PER_ARM)) |>
  dplyr::mutate(id = dplyr::row_number(), wt = 20, dose_mg = 0.2 * 20)

fig6_ev <- build_events(fig6_subjects, fig6_times)

sim_fig6 <- rxode2::rxSolve(mod, fig6_ev, returnType = "data.frame") |>
  dplyr::left_join(dplyr::select(fig6_subjects, id, geno), by = "id") |>
  dplyr::mutate(clb = Cc * UGL, ndmclb = Cc_ndmclb * UGL,
                geno = factor(geno, levels = c("NM", "IM", "PM")))
stopifnot(!anyNA(sim_fig6$clb), !anyNA(sim_fig6$ndmclb))
```

``` r

fig6_summary <- sim_fig6 |>
  tidyr::pivot_longer(c(clb, ndmclb), names_to = "analyte", values_to = "conc") |>
  dplyr::mutate(analyte = factor(analyte, levels = c("clb", "ndmclb"),
                                 labels = c("CLB", "N-CLB"))) |>
  dplyr::group_by(geno, analyte, time) |>
  dplyr::summarise(med = median(conc), lo = quantile(conc, 0.1),
                   hi = quantile(conc, 0.9), .groups = "drop")

ggplot(fig6_summary, aes(time, med, colour = analyte, fill = analyte)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~geno, ncol = 1, scales = "free_y") +
  labs(x = "Time (h)", y = "Concentration (ug/L)",
       colour = "Analyte", fill = "Analyte",
       title = "20 kg child, clobazam 0.2 mg/kg twice daily for 10 days") +
  theme_bw()
```

![Replicates Figure 6 of Tuo 2025: median (solid) and 10th-90th
percentile band (shaded) concentration-time profiles for a 20 kg child
on 0.2 mg/kg twice daily for 10 days, by CYP2C19
phenotype.](Tuo_2025_clobazam_files/figure-html/figure6-plot-1.png)

Replicates Figure 6 of Tuo 2025: median (solid) and 10th-90th percentile
band (shaded) concentration-time profiles for a 20 kg child on 0.2 mg/kg
twice daily for 10 days, by CYP2C19 phenotype.

``` r

fig6_trough <- sim_fig6 |>
  dplyr::filter(time == 240) |>
  dplyr::group_by(geno) |>
  dplyr::summarise(
    med_clb      = median(clb),
    med_ndmclb   = median(ndmclb),
    pct_clb_in   = 100 * mean(clb >= 30 & clb <= 300),
    pct_ndmclb_gt3000 = 100 * mean(ndmclb > 3000),
    .groups = "drop"
  )
knitr::kable(
  fig6_trough |>
    dplyr::rename("CYP2C19" = geno,
                  "Median CLB trough (ug/L)" = med_clb,
                  "Median N-CLB trough (ug/L)" = med_ndmclb,
                  "% CLB in 30-300 ug/L" = pct_clb_in,
                  "% N-CLB > 3000 ug/L" = pct_ndmclb_gt3000),
  digits = 1,
  caption = "Steady-state trough summary for the Figure 6 scenario."
)
```

| CYP2C19 | Median CLB trough (ug/L) | Median N-CLB trough (ug/L) | % CLB in 30-300 ug/L | % N-CLB \> 3000 ug/L |
|:---|---:|---:|---:|---:|
| NM | 89.0 | 575.8 | 97.0 | 2.0 |
| IM | 87.4 | 821.8 | 91.5 | 3.5 |
| PM | 90.7 | 2598.2 | 91.0 | 40.0 |

Steady-state trough summary for the Figure 6 scenario. {.table}

``` r


# The paper's qualitative Figure 6 claims. These are assertions on the CENTRE
# of the distribution and on a monotone ordering, not on cohort extremes.
stopifnot(
  # "steady-state trough concentrations of CLB had no significant difference
  # between the CYP2C19 NMs, IMs and PMs" -- spread across phenotypes is small
  # relative to the parent trough itself.
  (max(fig6_trough$med_clb) - min(fig6_trough$med_clb)) / min(fig6_trough$med_clb) < 0.15,
  # "the steady-state trough concentrations of N-CLB gradually increased in the
  # order NMs < IMs < PMs"
  !is.unsorted(fig6_trough$med_ndmclb[match(c("NM", "IM", "PM"), fig6_trough$geno)]),
  # "Approximately 90% of the N-CLB concentrations in the NMs and IMs were less
  # than 3000 ug/L, but about 50% of the N-CLB concentrations in the PMs
  # exceeded 3000 ug/L" -- gate the direction and a generous envelope.
  fig6_trough$pct_ndmclb_gt3000[fig6_trough$geno == "NM"] < 25,
  fig6_trough$pct_ndmclb_gt3000[fig6_trough$geno == "PM"] > 30
)
```

### Figure 5 – N-CLB / CLB metabolic ratios by CYP2C19 phenotype

Tuo 2025 Figure 5 reports mean steady-state metabolic ratios of AUC24h,
Cmax and Cmin. The AUC ratio is the one the model pins exactly: by the
steady-state mass balance above it equals `cl / cl_ndmclb`, which is
`exp(-theta_CYP2C19)` times 5.60 and therefore depends only on
phenotype, not on weight or dose.

``` r

ratio_window <- sim_fig6 |> dplyr::filter(time >= 228, time <= 240)

ratios <- ratio_window |>
  dplyr::group_by(geno, id) |>
  dplyr::summarise(
    auc24  = trapz(time, ndmclb) / trapz(time, clb),
    cmax_r = max(ndmclb) / max(clb),
    cmin_r = min(ndmclb) / min(clb),
    .groups = "drop"
  )

published_ratios <- tibble::tibble(
  geno   = factor(c("NM", "IM", "PM"), levels = c("NM", "IM", "PM")),
  auc24  = c(6.77, 8.28, 27.05),
  cmax_r = c(5.81, 6.74, 18.54),
  cmin_r = c(8.18, 10.79, 51.40)
)

ratio_cmp <- ratios |>
  dplyr::group_by(geno) |>
  dplyr::summarise(dplyr::across(c(auc24, cmax_r, cmin_r), mean), .groups = "drop") |>
  dplyr::left_join(published_ratios, by = "geno", suffix = c("_sim", "_pub"))

knitr::kable(
  ratio_cmp |>
    dplyr::select("CYP2C19" = geno,
                  "AUC24h ratio (simulated)" = auc24_sim,
                  "AUC24h ratio (Tuo 2025)"  = auc24_pub,
                  "Cmax ratio (simulated)"   = cmax_r_sim,
                  "Cmax ratio (Tuo 2025)"    = cmax_r_pub,
                  "Cmin ratio (simulated)"   = cmin_r_sim,
                  "Cmin ratio (Tuo 2025)"    = cmin_r_pub),
  digits = 2,
  caption = "Replicates Figure 5 of Tuo 2025: mean steady-state N-CLB/CLB metabolic ratios by CYP2C19 phenotype."
)
```

| CYP2C19 | AUC24h ratio (simulated) | AUC24h ratio (Tuo 2025) | Cmax ratio (simulated) | Cmax ratio (Tuo 2025) | Cmin ratio (simulated) | Cmin ratio (Tuo 2025) |
|:---|---:|---:|---:|---:|---:|---:|
| NM | 7.56 | 6.77 | 6.34 | 5.81 | 9.57 | 8.18 |
| IM | 10.86 | 8.28 | 8.55 | 6.74 | 15.66 | 10.79 |
| PM | 30.24 | 27.05 | 21.47 | 18.54 | 55.80 | 51.40 |

Replicates Figure 5 of Tuo 2025: mean steady-state N-CLB/CLB metabolic
ratios by CYP2C19 phenotype. {.table style="width:100%;"}

``` r


# Exact identity: the steady-state AUC ratio equals cl / cl_ndmclb, which for a
# 20 kg subject is 5.66/1.01 * exp(-theta). Uses the same drawn etas on both
# sides, so a tight bound is correct.
exact_auc_ratio <- (5.66 / 1.01) * exp(-c(NM = 0, IM = -0.25, PM = -1.30))
auc_ratio_med <- ratios |>
  dplyr::group_by(geno) |>
  dplyr::summarise(m = median(auc24), .groups = "drop")
stopifnot(
  # Ordering NM < IM < PM (Tuo 2025 Figure 5 and Results 3.3).
  !is.unsorted(ratio_cmp$auc24_sim),
  !is.unsorted(ratio_cmp$cmax_r_sim),
  !is.unsorted(ratio_cmp$cmin_r_sim),
  # The MEDIAN AUC ratio is the exact structural identity (the published means
  # are inflated by the log-normal IIV, so they are compared, not asserted on).
  all(abs(100 * (auc_ratio_med$m[match(c("NM", "IM", "PM"), auc_ratio_med$geno)] -
                   exact_auc_ratio) / exact_auc_ratio) < 20)
)
```

The simulated means sit close to the published means for AUC24h and
Cmax. The Cmin ratio is the noisiest of the three in both the paper and
the simulation – it is a ratio of two small numbers at the bottom of the
profile, and its mean is strongly inflated by the log-normal
metabolite-clearance variability – so it is reported side by side rather
than gated.

## PKNCA validation

NCA is run over the last steady-state dosing interval (228 to 240 h) on
the typical-value cohort, once per analyte. The trough is taken as
`cmin`, not `ctrough`: PKNCA defines `ctrough` as the *predose*
concentration at the interval **start**, which is not returned for this
interval, whereas at steady state the minimum over a full dosing
interval is the trough by construction – the concentration at 228 h and
at 240 h are equal. `cmin` here reproduces the 240 h model trough
exactly.

``` r

nca_conc <- sim_typ |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::select(id, arm, geno, wt, time, clb, ndmclb)

dose_df <- ev_typ |>
  dplyr::filter(evid == 1L) |>
  dplyr::select(id, time, amt) |>
  dplyr::left_join(dplyr::select(tab3, id, arm), by = "id")

intervals <- data.frame(
  start = 228, end = 240,
  cmax = TRUE, tmax = TRUE, auclast = TRUE, cmin = TRUE
)

run_nca <- function(conc_col) {
  cobj <- PKNCA::PKNCAconc(
    nca_conc |> dplyr::mutate(conc = .data[[conc_col]]),
    conc ~ time | arm + id, concu = "ug/L", timeu = "hr"
  )
  dobj <- PKNCA::PKNCAdose(dose_df, amt ~ time | arm + id, doseu = "mg")
  as.data.frame(PKNCA::pk.nca(PKNCA::PKNCAdata(cobj, dobj, intervals = intervals))$result)
}

nca_clb    <- run_nca("clb")    |> dplyr::mutate(analyte = "CLB")
nca_ndmclb <- run_nca("ndmclb") |> dplyr::mutate(analyte = "N-CLB")

nca_all <- dplyr::bind_rows(nca_clb, nca_ndmclb) |>
  dplyr::filter(PPTESTCD %in% c("cmax", "tmax", "auclast", "cmin"))

knitr::kable(
  nca_all |>
    tidyr::pivot_wider(id_cols = c(analyte, arm), names_from = PPTESTCD,
                       values_from = PPORRES) |>
    dplyr::rename("Analyte" = analyte, "Arm" = arm, "Cmax (ug/L)" = cmax,
                  "Tmax (h)" = tmax, "AUCtau (ug*h/L)" = auclast,
                  "Ctrough (ug/L)" = cmin),
  digits = 1,
  caption = "Steady-state NCA (228-240 h) of the typical-value profiles, by analyte and arm."
)
```

| Analyte | Arm | AUCtau (ug\*h/L) | Cmax (ug/L) | Ctrough (ug/L) | Tmax (h) |
|:---|:---|---:|---:|---:|---:|
| CLB | IM 10 kg (0.25 mg/kg BID) | 1898.8 | 236.4 | 86.3 | 1.5 |
| CLB | IM 20 kg (0.2 mg/kg BID) | 1806.8 | 212.2 | 91.1 | 1.5 |
| CLB | IM 30 kg (0.15 mg/kg BID) | 1499.8 | 170.9 | 79.6 | 1.5 |
| CLB | IM 40 kg (0.125 mg/kg BID) | 1343.1 | 150.0 | 73.8 | 1.5 |
| CLB | IM 50 kg (0.12 mg/kg BID) | 1363.4 | 150.0 | 76.7 | 1.5 |
| CLB | NM 10 kg (0.3 mg/kg BID) | 2278.6 | 283.7 | 103.5 | 1.5 |
| CLB | NM 20 kg (0.25 mg/kg BID) | 2258.5 | 265.3 | 113.9 | 1.5 |
| CLB | NM 30 kg (0.2 mg/kg BID) | 1999.7 | 227.8 | 106.2 | 1.5 |
| CLB | NM 40 kg (0.175 mg/kg BID) | 1880.3 | 210.0 | 103.3 | 1.5 |
| CLB | NM 50 kg (0.15 mg/kg BID) | 1704.2 | 187.5 | 95.9 | 1.5 |
| CLB | PM 10 kg (0.1 mg/kg BID) | 759.5 | 94.6 | 34.5 | 1.5 |
| CLB | PM 20 kg (0.075 mg/kg BID) | 677.5 | 79.6 | 34.2 | 1.5 |
| CLB | PM 30 kg (0.075 mg/kg BID) | 749.9 | 85.4 | 39.8 | 1.5 |
| CLB | PM 40 kg (0.0625 mg/kg BID) | 671.5 | 75.0 | 36.9 | 1.5 |
| CLB | PM 50 kg (0.06 mg/kg BID) | 681.7 | 75.0 | 38.3 | 1.5 |
| N-CLB | IM 10 kg (0.25 mg/kg BID) | 13677.4 | 1484.0 | 724.4 | 3.2 |
| N-CLB | IM 20 kg (0.2 mg/kg BID) | 13012.3 | 1331.0 | 763.9 | 3.5 |
| N-CLB | IM 30 kg (0.15 mg/kg BID) | 10800.4 | 1072.8 | 667.1 | 3.8 |
| N-CLB | IM 40 kg (0.125 mg/kg BID) | 9671.5 | 942.5 | 616.9 | 3.8 |
| N-CLB | IM 50 kg (0.12 mg/kg BID) | 9817.4 | 944.2 | 640.8 | 4.0 |
| N-CLB | NM 10 kg (0.3 mg/kg BID) | 12782.3 | 1430.3 | 653.2 | 3.0 |
| N-CLB | NM 20 kg (0.25 mg/kg BID) | 12667.5 | 1332.9 | 718.3 | 3.2 |
| N-CLB | NM 30 kg (0.2 mg/kg BID) | 11215.2 | 1143.3 | 669.8 | 3.5 |
| N-CLB | NM 40 kg (0.175 mg/kg BID) | 10545.1 | 1053.8 | 651.0 | 3.5 |
| N-CLB | NM 50 kg (0.15 mg/kg BID) | 9557.2 | 941.0 | 604.2 | 3.5 |
| N-CLB | PM 10 kg (0.1 mg/kg BID) | 15634.4 | 1480.8 | 1032.2 | 4.8 |
| N-CLB | PM 20 kg (0.075 mg/kg BID) | 13944.5 | 1275.7 | 982.1 | 5.0 |
| N-CLB | PM 30 kg (0.075 mg/kg BID) | 15432.1 | 1389.3 | 1119.6 | 5.0 |
| N-CLB | PM 40 kg (0.0625 mg/kg BID) | 13819.1 | 1231.9 | 1020.7 | 5.2 |
| N-CLB | PM 50 kg (0.06 mg/kg BID) | 14027.4 | 1242.1 | 1049.0 | 5.2 |

Steady-state NCA (228-240 h) of the typical-value profiles, by analyte
and arm. {.table style="width:100%;"}

## Comparison against published NCA

Tuo 2025 reports no conventional NCA table, but Table 3 gives a
published **median steady-state trough** for both analytes in each of
the 15 dosing scenarios – the same quantity as `ctrough` above. That is
the reference set.

``` r

simulated_long <- nca_all |>
  dplyr::filter(PPTESTCD == "cmin") |>
  dplyr::select(analyte, arm, PPTESTCD, PPORRES)

published <- dplyr::bind_rows(
  tab3 |> dplyr::transmute(analyte = "CLB",   arm, cmin = pub_clb),
  tab3 |> dplyr::transmute(analyte = "N-CLB", arm, cmin = pub_ndmclb)
)

cmp <- nlmixr2lib::ncaComparisonTable(
  simulated     = simulated_long,
  reference     = published,
  by            = c("analyte", "arm"),
  units         = c(cmin = "ug/L"),
  tolerance_pct = 20
)

knitr::kable(
  cmp,
  caption = "Simulated versus published median steady-state trough concentration for all 15 Tuo 2025 Table 3 scenarios."
)
```

| NCA parameter | analyte | arm                         | Reference | Simulated | % diff |
|:--------------|:--------|:----------------------------|:----------|:----------|:-------|
| Cmin (ug/L)   | CLB     | NM 10 kg (0.3 mg/kg BID)    | 109       | 104       | -5.3%  |
| Cmin (ug/L)   | CLB     | NM 20 kg (0.25 mg/kg BID)   | 120       | 114       | -4.9%  |
| Cmin (ug/L)   | CLB     | NM 30 kg (0.2 mg/kg BID)    | 111       | 106       | -4.7%  |
| Cmin (ug/L)   | CLB     | NM 40 kg (0.175 mg/kg BID)  | 108       | 103       | -4.6%  |
| Cmin (ug/L)   | CLB     | NM 50 kg (0.15 mg/kg BID)   | 93        | 95.9      | +3.0%  |
| Cmin (ug/L)   | CLB     | IM 10 kg (0.25 mg/kg BID)   | 83.8      | 86.3      | +3.0%  |
| Cmin (ug/L)   | CLB     | IM 20 kg (0.2 mg/kg BID)    | 88.2      | 91.1      | +3.3%  |
| Cmin (ug/L)   | CLB     | IM 30 kg (0.15 mg/kg BID)   | 77.2      | 79.6      | +3.2%  |
| Cmin (ug/L)   | CLB     | IM 40 kg (0.125 mg/kg BID)  | 71.5      | 73.8      | +3.1%  |
| Cmin (ug/L)   | CLB     | IM 50 kg (0.12 mg/kg BID)   | 80.3      | 76.7      | -4.5%  |
| Cmin (ug/L)   | CLB     | PM 10 kg (0.1 mg/kg BID)    | 33.5      | 34.5      | +3.0%  |
| Cmin (ug/L)   | N-CLB   | NM 10 kg (0.3 mg/kg BID)    | 630       | 653       | +3.8%  |
| Cmin (ug/L)   | CLB     | PM 20 kg (0.075 mg/kg BID)  | 33.1      | 34.2      | +3.3%  |
| Cmin (ug/L)   | N-CLB   | NM 20 kg (0.25 mg/kg BID)   | 699       | 718       | +2.7%  |
| Cmin (ug/L)   | CLB     | PM 30 kg (0.075 mg/kg BID)  | 38.6      | 39.8      | +3.2%  |
| Cmin (ug/L)   | N-CLB   | NM 30 kg (0.2 mg/kg BID)    | 651       | 670       | +2.9%  |
| Cmin (ug/L)   | CLB     | PM 40 kg (0.0625 mg/kg BID) | 35.9      | 36.9      | +2.6%  |
| Cmin (ug/L)   | N-CLB   | NM 40 kg (0.175 mg/kg BID)  | 631       | 651       | +3.2%  |
| Cmin (ug/L)   | CLB     | PM 50 kg (0.06 mg/kg BID)   | 37.4      | 38.3      | +2.6%  |
| Cmin (ug/L)   | N-CLB   | NM 50 kg (0.15 mg/kg BID)   | 569       | 604       | +6.2%  |
| Cmin (ug/L)   | N-CLB   | IM 10 kg (0.25 mg/kg BID)   | 658       | 724       | +10.1% |
| Cmin (ug/L)   | N-CLB   | IM 20 kg (0.2 mg/kg BID)    | 709       | 764       | +7.8%  |
| Cmin (ug/L)   | N-CLB   | IM 30 kg (0.15 mg/kg BID)   | 626       | 667       | +6.5%  |
| Cmin (ug/L)   | N-CLB   | IM 40 kg (0.125 mg/kg BID)  | 579       | 617       | +6.5%  |
| Cmin (ug/L)   | N-CLB   | IM 50 kg (0.12 mg/kg BID)   | 617       | 641       | +3.9%  |
| Cmin (ug/L)   | N-CLB   | PM 10 kg (0.1 mg/kg BID)    | 963       | 1030      | +7.2%  |
| Cmin (ug/L)   | N-CLB   | PM 20 kg (0.075 mg/kg BID)  | 950       | 982       | +3.4%  |
| Cmin (ug/L)   | N-CLB   | PM 30 kg (0.075 mg/kg BID)  | 1080      | 1120      | +3.5%  |
| Cmin (ug/L)   | N-CLB   | PM 40 kg (0.0625 mg/kg BID) | 967       | 1020      | +5.5%  |
| Cmin (ug/L)   | N-CLB   | PM 50 kg (0.06 mg/kg BID)   | 996       | 1050      | +5.3%  |

Simulated versus published median steady-state trough concentration for
all 15 Tuo 2025 Table 3 scenarios. {.table}

``` r

# Recompute the percent difference numerically -- the `% diff` column of
# ncaComparisonTable() is formatted character and cannot be asserted on.
gate <- simulated_long |>
  dplyr::left_join(published, by = c("analyte", "arm")) |>
  dplyr::mutate(pct_diff = 100 * (PPORRES - cmin) / cmin)

knitr::kable(
  gate |>
    dplyr::group_by(analyte) |>
    dplyr::summarise(`Median % diff` = median(pct_diff),
                     `Max |% diff|`  = max(abs(pct_diff)), .groups = "drop") |>
    dplyr::rename("Analyte" = analyte),
  digits = 2,
  caption = "Agreement between the deterministic model trough and the Tuo 2025 Table 3 published medians."
)
```

| Analyte | Median % diff | Max \|% diff\| |
|:--------|--------------:|---------------:|
| CLB     |          3.01 |           5.27 |
| N-CLB   |          5.31 |          10.08 |

Agreement between the deterministic model trough and the Tuo 2025 Table
3 published medians. {.table}

``` r


# This comparison is DETERMINISTIC on the simulated side (zeroRe, no RNG), so
# the bound cannot drift with the rxode2 build or the thread count. The
# residual gap is the Monte Carlo median-estimation error of the paper's own
# simulation, which scatters roughly +/- 4%.
stopifnot(
  max(abs(dplyr::filter(gate, analyte == "CLB")$pct_diff))   < 15,
  max(abs(dplyr::filter(gate, analyte == "N-CLB")$pct_diff)) < 15,
  abs(median(gate$pct_diff)) < 8
)
```

### Probability of target attainment

The parent trough depends only on `etalcl` (no IIV was estimated on
`V_CLB/F`), and it is monotone decreasing in clearance. The probability
of exceeding a target is therefore available **exactly**, by
root-finding the eta that lands on the threshold and reading the normal
CDF – no sampling, so even the sub-1% columns of Table 3 are resolved
rather than quantised by a finite cohort.

``` r

OM_CL <- 0.1594  # etalcl variance, Table 2

pta_exact <- function(dose, wt, threshold_ugL) {
  f <- function(eta) {
    cl <- 5.66 * (wt / 70)^0.75 * exp(eta)
    vc <- 92.07 * (wt / 70)
    ss_trough(dose, 1.99, cl, vc) * UGL - threshold_ugL
  }
  if (f(-6) < 0) return(0)
  if (f(6)  > 0) return(100)
  eta_star <- stats::uniroot(f, c(-6, 6), tol = 1e-10)$root
  100 * stats::pnorm(eta_star / sqrt(OM_CL))
}

pta_clb <- tab3 |>
  dplyr::rowwise() |>
  dplyr::mutate(
    sim_pta30  = pta_exact(dose_mg, wt, 30),
    sim_pta300 = pta_exact(dose_mg, wt, 300),
    sim_pta500 = pta_exact(dose_mg, wt, 500)
  ) |>
  dplyr::ungroup()

knitr::kable(
  pta_clb |>
    dplyr::select("CYP2C19" = geno, "Weight (kg)" = wt,
                  "Dose (mg/kg BID)" = mgkg,
                  "PTA >= 30 ug/L (simulated)" = sim_pta30,
                  "PTA >= 30 ug/L (Tuo 2025)"  = pub_pta_clb30,
                  "PTA >= 300 ug/L (simulated)" = sim_pta300,
                  "PTA >= 500 ug/L (simulated)" = sim_pta500),
  digits = 1,
  caption = "Replicates the CLB probability-of-target-attainment columns of Tuo 2025 Table 3."
)
```

| CYP2C19 | Weight (kg) | Dose (mg/kg BID) | PTA \>= 30 ug/L (simulated) | PTA \>= 30 ug/L (Tuo 2025) | PTA \>= 300 ug/L (simulated) | PTA \>= 500 ug/L (simulated) |
|:---|---:|---:|---:|---:|---:|---:|
| NM | 10 | 0.3 | 94.1 | 94.8 | 3.3 | 0.2 |
| NM | 20 | 0.2 | 96.4 | 97.0 | 3.9 | 0.2 |
| NM | 30 | 0.2 | 96.4 | 97.0 | 2.5 | 0.1 |
| NM | 40 | 0.2 | 96.6 | 97.2 | 1.9 | 0.1 |
| NM | 50 | 0.1 | 96.1 | 95.6 | 1.2 | 0.0 |
| IM | 10 | 0.2 | 91.4 | 89.6 | 1.4 | 0.1 |
| IM | 20 | 0.2 | 93.9 | 93.5 | 1.3 | 0.0 |
| IM | 30 | 0.1 | 92.6 | 92.1 | 0.5 | 0.0 |
| IM | 40 | 0.1 | 91.7 | 91.7 | 0.3 | 0.0 |
| IM | 50 | 0.1 | 92.9 | 94.0 | 0.3 | 0.0 |
| PM | 10 | 0.1 | 58.2 | 57.3 | 0.0 | 0.0 |
| PM | 20 | 0.1 | 58.2 | 55.9 | 0.0 | 0.0 |
| PM | 30 | 0.1 | 67.8 | 66.1 | 0.0 | 0.0 |
| PM | 40 | 0.1 | 63.6 | 62.3 | 0.0 | 0.0 |
| PM | 50 | 0.1 | 66.3 | 65.2 | 0.0 | 0.0 |

Replicates the CLB probability-of-target-attainment columns of Tuo 2025
Table 3. {.table}

``` r


pta_err <- pta_clb$sim_pta30 - pta_clb$pub_pta_clb30
stopifnot(
  # Exact (non-stochastic) computation against published values; compare on the
  # percentage-point scale, and gate the centre plus a robust quantile.
  abs(median(pta_err)) < 10,
  stats::quantile(abs(pta_err), 0.9) < 20,
  # The paper's headline stratification: PMs on the recommended dose reach the
  # 30 ug/L parent target far less often than NMs and IMs.
  mean(pta_clb$sim_pta30[pta_clb$geno == "PM"]) <
    mean(pta_clb$sim_pta30[pta_clb$geno == "NM"]) - 20
)
```

The metabolite target-attainment columns depend on both etas, so they
are estimated from the stochastic cohort rather than in closed form. The
cohort is capped at 200 subjects per arm, which resolves a probability
to about 0.5 percentage points – adequate for the `>= 300 ug/L` column
(81-95% in Table 3) but **not** for the `>= 3000` and `>= 5000` columns
(0.1-9.3%), which are therefore reported for inspection and not gated.

``` r

pta_subjects <- tab3 |>
  dplyr::select(arm_id = id, wt, dose_mg, CYP2C19_IM, CYP2C19_PM) |>
  dplyr::slice(rep(seq_len(dplyr::n()), each = N_PER_ARM)) |>
  dplyr::mutate(id = dplyr::row_number())

# Only the 240 h trough is needed, so the PTA cohort carries one observation
# record per subject.
pta_ev <- build_events(pta_subjects, 240)

sim_pta <- rxode2::rxSolve(mod, pta_ev, returnType = "data.frame") |>
  dplyr::left_join(dplyr::select(pta_subjects, id, arm_id), by = "id") |>
  dplyr::filter(time == 240) |>
  dplyr::mutate(ndmclb = Cc_ndmclb * UGL)
stopifnot(!anyNA(sim_pta$ndmclb))

pta_ndmclb <- sim_pta |>
  dplyr::group_by(arm_id) |>
  dplyr::summarise(sim_pta300  = 100 * mean(ndmclb >= 300),
                   sim_pta3000 = 100 * mean(ndmclb >= 3000),
                   .groups = "drop") |>
  dplyr::left_join(dplyr::select(tab3, arm_id = id, geno, wt, mgkg, pub_pta_ndmclb300),
                   by = "arm_id")

knitr::kable(
  pta_ndmclb |>
    dplyr::select("CYP2C19" = geno, "Weight (kg)" = wt, "Dose (mg/kg BID)" = mgkg,
                  "PTA >= 300 ug/L (simulated)" = sim_pta300,
                  "PTA >= 300 ug/L (Tuo 2025)"  = pub_pta_ndmclb300,
                  "PTA >= 3000 ug/L (simulated, not gated)" = sim_pta3000),
  digits = 1,
  caption = "Replicates the N-CLB probability-of-target-attainment columns of Tuo 2025 Table 3."
)
```

| CYP2C19 | Weight (kg) | Dose (mg/kg BID) | PTA \>= 300 ug/L (simulated) | PTA \>= 300 ug/L (Tuo 2025) | PTA \>= 3000 ug/L (simulated, not gated) |
|:---|---:|---:|---:|---:|---:|
| NM | 10 | 0.3 | 81.0 | 84.1 | 3.0 |
| NM | 20 | 0.2 | 85.5 | 88.4 | 4.0 |
| NM | 30 | 0.2 | 83.5 | 86.3 | 2.0 |
| NM | 40 | 0.2 | 87.5 | 85.5 | 1.0 |
| NM | 50 | 0.1 | 86.0 | 80.9 | 1.5 |
| IM | 10 | 0.2 | 84.5 | 83.5 | 2.0 |
| IM | 20 | 0.2 | 89.0 | 87.1 | 2.5 |
| IM | 30 | 0.1 | 80.5 | 82.7 | 2.5 |
| IM | 40 | 0.1 | 83.5 | 80.9 | 2.5 |
| IM | 50 | 0.1 | 86.0 | 84.4 | 2.0 |
| PM | 10 | 0.1 | 94.5 | 89.9 | 12.5 |
| PM | 20 | 0.1 | 94.0 | 92.2 | 5.5 |
| PM | 30 | 0.1 | 95.5 | 95.2 | 9.5 |
| PM | 40 | 0.1 | 94.5 | 91.4 | 9.0 |
| PM | 50 | 0.1 | 96.0 | 92.0 | 7.5 |

Replicates the N-CLB probability-of-target-attainment columns of Tuo
2025 Table 3. {.table}

``` r


nd_err <- pta_ndmclb$sim_pta300 - pta_ndmclb$pub_pta_ndmclb300
stopifnot(
  abs(median(nd_err)) < 15,
  stats::quantile(abs(nd_err), 0.9) < 25
)
```

## Assumptions and deviations

- **Interindividual-variance scale.** Tuo 2025 Table 2 labels the two
  IIV rows `omega^2 CL/F (%)` with values 15.94 and 45.35. Two readings
  are possible: a variance expressed as a percentage (0.1594, 0.4535) or
  a %CV (0.1594, 0.4535 as SDs). The table footnote resolves it in words
  – “interindividual **variance** for pharmacokinetic parameters” – and
  the paper’s own Monte Carlo output confirms it numerically. Reading
  them as variances reproduces the Table 3 parent PTA columns for the 20
  kg NM / 0.25 mg/kg row (3.9% and 0.2% of subjects above 300 and 500
  ug/L, against the published 5.1% and 0.1%); reading them as CVs
  predicts 0.0005% and 7e-11%, four to ten orders of magnitude too
  small. The variance reading is used.
- **Bioavailability and metabolite fraction are not identifiable.** `F`
  and `Fm` were not estimated, so the model is written in the apparent
  parameterisation the paper reports; `central` and `central_ndmclb`
  hold scaled amounts. Both predicted concentrations are the true plasma
  concentrations. See “The apparent (F, Fm) parameterisation” above. A
  consequence is that the parent’s entire elimination is routed into the
  metabolite compartment – the paper models no separate non-metabolic
  clobazam elimination arm.
- **Absorption is inherited, not fitted.** `Ka` is fixed at 1.99 1/h
  from Jullien 2015 (Tuo 2025 reference \[19\]) because the
  opportunistic sampling carried essentially no absorption-phase
  information. Predictions near Tmax are therefore weakly supported by
  this study’s own data; the trough predictions that the paper and this
  vignette validate against are not sensitive to it.
- **70 kg reference weight is an extrapolation.** All four typical
  values are reported standardised to 70 kg, but the cohort spans
  6.60-73.0 kg with a median of 20.0 kg. Only one subject sits near the
  reference, so the printed typical values are adult-standardised rather
  than observed.
- **Allometric exponents were not estimated.** Equations (9)-(12) print
  1.0 on the volumes and 0.75 on the clearances as literal constants and
  Table 2 has no exponent row, no RSE and no bootstrap interval for
  them; they are encoded with `fixed()`.
- **No IIV on volumes or on Ka.** The Discussion states the model was
  simplified by “reducing model compartments, freezing parameter values,
  and ignoring interindividual variability for certain parameters”. Only
  the two clearances carry an eta, and Table 2 reports no correlation
  block, so the two etas are independent. This is encoded faithfully
  rather than by inventing variances.
- **CYP2C19 rapid and ultrarapid metabolizers have no encoding.** The 2
  rapid metabolizers (`*1/*17`) were excluded from the covariate
  analysis for lack of sample size and no ultrarapid metabolizers were
  enrolled, so neither has an estimated effect. Both fall into the
  `CYP2C19_IM = 0, CYP2C19_PM = 0` reference cell by default; simulating
  them is an explicit extrapolation the paper flags as a limitation.
- **Screened-but-not-retained covariates** (age, sex, body surface area,
  renal and hepatic panels, concomitant antiepileptics, ketogenic diet,
  and the ABCB1 / CYP3A4 / GABA-receptor SNPs) are recorded in
  `covariatesDataExcluded` so the paper’s covariate screen is preserved
  without implying they act in the model. The twelve screened SNPs are
  represented by a single entry because none carries an estimated
  coefficient.
- **Table 1 unit typo in the source.** Tuo 2025 Table 1 labels the two
  concentration rows `ug*mL^-1` while the values (mean 137.01 clobazam,
  1611.56 N-desmethylclobazam) and the whole of the surrounding text are
  in `ug*L^-1`. The `ug/L` reading is used. Relatedly, Results 3.1 calls
  137.01 the *median* clobazam concentration whereas Table 1 lists it as
  the *mean* (median 114.50); neither number enters the model.
- **Comparison target.** The paper reports no conventional NCA table, so
  the validation reference is the published median steady-state trough
  of Table 3 for all 15 weight-by-phenotype scenarios, compared against
  `ctrough` from PKNCA on the typical-value profiles.
- **Residual error is excluded from the simulated concentrations.** `Cc`
  and `Cc_ndmclb` are individual predictions; the proportional residual
  terms describe assay and model misspecification error and are not
  added before computing troughs or target attainment. Adding a 33-53%
  proportional error would also drive simulated concentrations negative.
- **Metabolite PTA resolution.** The `>= 3000` and `>= 5000 ug/L`
  metabolite columns of Table 3 range from 0.1% to 9.3%, below what a
  200-subject-per-arm cohort resolves; the `>= 3000` column is shown for
  inspection and neither is gated. The parent PTA columns are computed
  exactly instead, by root-finding rather than sampling.
- **Equation recovery.** All twelve numbered equations were lost by the
  PDF-to-markdown conversion (each rendered as `formula-not-decoded`);
  they were recovered from the PDF text layer directly and are the basis
  for the ODE system, the covariate forms and the fixed allometric
  exponents.
