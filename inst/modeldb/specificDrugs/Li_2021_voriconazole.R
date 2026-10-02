Li_2021_voriconazole <- function() {
  description <- "Joint parent-metabolite population pharmacokinetic model for oral voriconazole and voriconazole N-oxide in Chinese immunocompromised patients with invasive fungal infection (Li 2021). Voriconazole: one compartment with first-order absorption (ka and F fixed) and parallel linear plus Michaelis-Menten elimination, where the Michaelis-Menten pathway forms the N-oxide and is inhibited by the N-oxide concentration through a fixed Imax/IC50 term. Voriconazole N-oxide: one compartment with first-order elimination. CYP2C19 metabolizer phenotype (IM, PM vs NM) acts exponentially on Vmax."
  reference <- paste(
    "Li S, Wu S, Gong W, Cao P, Chen X, Liu W, Xiang L, Wang Y, Huang J.",
    "Application of Population Pharmacokinetic Analysis to Characterize",
    "CYP2C19 Mediated Metabolic Mechanism of Voriconazole and Support Dose",
    "Optimization. Front Pharmacol. 2021;12:730826 (published 3 January 2022).",
    "doi:10.3389/fphar.2021.730826. PMCID PMC8762230.",
    sep = " "
  )
  vignette <- "Li_2021_voriconazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "voriconazole", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "voriconazole", units = "mg", specimen = "plasma", verified = TRUE),
    central_noxvori = list(analyte = "voriconazole N-oxide", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CYP2C19_IM = list(
      description = "CYP2C19 intermediate-metabolizer phenotype indicator; exponential effect on the Michaelis-Menten Vmax of voriconazole",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal metabolizer, the NM reference group, when both CYP2C19_IM and CYP2C19_PM are 0)",
      notes = "1 = CYP2C19 intermediate metabolizer. Li 2021 Methods 'Genotype Analysis' assigns phenotype per the CPIC guideline: NM = *1/*1; IM = one loss-of-function allele with one normal allele (*1/*2, *1/*3, *2/*17); PM = *2/*2, *2/*3, *3/*3. No *17 allele was detected in the cohort, so no RM/UM subjects were analysed. 32 of the 75 genotyped patients were IM (Results 'Patient Characteristics'). Paired with CYP2C19_PM; both 0 is the NM reference (theta_NM = 0 FIX in Table 2).",
      source_name = "CYP2C19 phenotype (IM)"
    ),
    CYP2C19_PM = list(
      description = "CYP2C19 poor-metabolizer phenotype indicator; exponential effect on the Michaelis-Menten Vmax of voriconazole",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal metabolizer, the NM reference group, when both CYP2C19_IM and CYP2C19_PM are 0)",
      notes = "1 = CYP2C19 poor metabolizer (*2/*2, *2/*3, *3/*3 per Li 2021 Methods 'Genotype Analysis'). 16 of the 75 genotyped patients were PM. Paired with CYP2C19_IM; both 0 is the NM reference.",
      source_name = "CYP2C19 phenotype (PM)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 78L,
    n_studies = 1L,
    age_range = "14-70 years",
    age_median = "36.5 years",
    weight_range = "44-111 kg",
    weight_median = "64 kg",
    sex_female_pct = 26.9,
    race_ethnicity = c(Asian = 100),
    disease_state = "Immunocompromised patients receiving oral voriconazole for possible, probable or proven invasive fungal infection. Four patients were adolescents (14-17 years), the rest adults. CYP2C19 phenotype (75 genotyped): NM 27, IM 32, PM 16; three patients had no genotype. Co-medication with proton-pump inhibitors in 46.7% and glucocorticoids in 32.7% of samples; neither was a significant covariate.",
    dose_range = "200 mg oral voriconazole twice daily without a loading dose.",
    regions = "China (Union Hospital, Wuhan; single centre, February 2017 - July 2018)",
    notes = "Retrospective TDM study. 214 voriconazole and 213 voriconazole N-oxide plasma concentrations (LC-MS/MS, LLOQ 0.5 ng/mL for both), up to 8 per patient, almost all pre-dose troughs; sampling 23-4,223 h after the first dose. Demographics from Li 2021 Table 1 (57 male : 21 female). Phoenix NLME 8.2, FOCE-ELS."
  )

  # Implementation notes (the vignette 'Assumptions and deviations' section
  # carries the full justification):
  #
  # * Structure (Li 2021 Equations 8-13). The Michaelis-Menten clearance
  #     CLnonlin = Vmax / (C1 + km) * (1 - Imax * C2 / (IC50 + C2))
  #   is the voriconazole-to-N-oxide formation pathway and is inhibited by
  #   the N-oxide concentration C2; CL1 is voriconazole clearance by all
  #   other routes.
  #
  # * Equation 10 multiplies the formation flux by a factor k_n that is not
  #   defined or valued anywhere in the paper or its supplement, and is not
  #   one of the 11/13 estimated parameters (Results 'PopPK Model
  #   Development'). The prose defines the N-oxide input rate as "the same
  #   as the conversion rate from VCZ to VNO" and both analytes are in mg/L,
  #   so k_n is taken as the stoichiometric molecular-weight ratio
  #   M(N-oxide) / M(voriconazole) = 365.31 / 349.31 = 1.0458 (voriconazole
  #   C16H14F3N5O plus one oxygen). k_n = 1 would lower N-oxide
  #   concentrations by 4.4% and leave voriconazole essentially unchanged.
  #
  # * Table 2 reports the IIV as 'omega (%)' and its footnote defines omega
  #   as the 'square root of interindividual variance', so each value / 100
  #   is the SD of eta and is squared here. Table 3's Monte Carlo median
  #   day-20 troughs adjudicate this: they are reproduced to within about 5%
  #   under the SD reading and are 25-32% too high if 240.77% is read as a CV.

  ini({
    # Absorption, fixed from the literature (Li 2021 Methods 'Base Model')
    lka <- fixed(log(1.1)); label("Absorption rate constant ka (1/h)") # Li 2021 Table 2 'Ka (h-1) 1.1 Fix' (from Pascual 2012 / Wang 2014)
    lfdepot <- fixed(log(0.895)); label("Oral bioavailability F (fraction)") # Li 2021 Table 2 'F 0.895 Fix'

    # Voriconazole disposition
    lvc <- log(207.29); label("Voriconazole volume of distribution V1 (L)") # Li 2021 Table 2 'V1 (L) 207.29'
    lcl <- log(1.91); label("Voriconazole clearance other than the N-oxide pathway CL1 (L/h)") # Li 2021 Table 2 'CL1 (L/h) 1.91'
    lvmax <- log(18.80); label("Maximum N-oxidation rate Vmax for a CYP2C19 NM (mg/h)") # Li 2021 Table 2 'Vmax (mg/h) 18.80'
    lkm <- fixed(log(1.15)); label("Michaelis-Menten constant km (mg/L)") # Li 2021 Table 2 'Km 1.15 Fix'; units mg/L per Results text below Table 2
    imax <- fixed(0.75); label("Maximal fractional inhibition of the N-oxidation clearance by the N-oxide Imax (fraction)") # Li 2021 Table 2 'Imax 0.75 Fix'
    lic50 <- fixed(log(14.6)); label("N-oxide concentration giving half-maximal inhibition IC50 (mg/L)") # Li 2021 Table 2 'IC50 (mg/L) 14.6 Fix'

    # Voriconazole N-oxide disposition
    lvc_noxvori <- log(10.01); label("Voriconazole N-oxide volume of distribution V2 (L)") # Li 2021 Table 2 'V2 (L) 10.01'
    lcl_noxvori <- log(4.65); label("Voriconazole N-oxide clearance CL2 (L/h)") # Li 2021 Table 2 'CL2 (L/h) 4.65'

    # CYP2C19 phenotype on Vmax: Vmax = 18.80 * exp(theta_CYP2C19), theta_NM = 0 FIX
    e_cyp2c19_im_vmax <- -0.31; label("CYP2C19 IM effect on Vmax, exponential (unitless)") # Li 2021 Table 2 'theta IM -0.31' and footnote a
    e_cyp2c19_pm_vmax <- -0.61; label("CYP2C19 PM effect on Vmax, exponential (unitless)") # Li 2021 Table 2 'theta PM -0.61' and footnote a

    # IIV: Table 2 omega (%) = 100 * SD of eta; variances = (omega / 100)^2
    etalvc ~ 5.79702 # Li 2021 Table 2 'omega V1 (%) 240.77' -> 2.4077^2
    etalcl ~ 0.00362404 # Li 2021 Table 2 'omega CL1 (%) 6.02' -> 0.0602^2
    etalcl_noxvori ~ 0.06538249 # Li 2021 Table 2 'omega CL2 (%) 25.57' -> 0.2557^2
    etalvmax ~ 0.04464769 # Li 2021 Table 2 'omega Vmax (%) 21.13' -> 0.2113^2

    # Residual error: proportional for both analytes
    propSd <- 0.4697; label("Voriconazole proportional residual error (fraction)") # Li 2021 Table 2 'VCZ-sigma (%) 46.97'
    propSd_noxvori <- 0.2793; label("Voriconazole N-oxide proportional residual error (fraction)") # Li 2021 Table 2 'VNO-sigma (%) 27.93'
  })

  model({
    # Stoichiometric factor k_n of Li 2021 Equation 10 (value not printed;
    # see the implementation notes above): N-oxide / voriconazole molar mass
    mw_ratio <- 365.31 / 349.31

    # Individual parameters
    ka <- exp(lka)
    fdepot <- exp(lfdepot)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl)
    vmax <- exp(lvmax + e_cyp2c19_im_vmax * CYP2C19_IM + e_cyp2c19_pm_vmax * CYP2C19_PM + etalvmax)
    km <- exp(lkm)
    ic50 <- exp(lic50)
    vc_noxvori <- exp(lvc_noxvori)
    cl_noxvori <- exp(lcl_noxvori + etalcl_noxvori)

    # Concentrations (Equations 12-13)
    Cc <- central / vc
    Cc_noxvori <- central_noxvori / vc_noxvori

    # N-oxide-inhibited Michaelis-Menten formation clearance (Equation 11)
    cl_nonlin <- vmax / (Cc + km) * (1 - imax * Cc_noxvori / (ic50 + Cc_noxvori))

    # ODEs (Equations 8-10)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl / vc * central - cl_nonlin / vc * central
    d/dt(central_noxvori) <- mw_ratio * cl_nonlin / vc * central - cl_noxvori / vc_noxvori * central_noxvori

    f(depot) <- fdepot

    Cc ~ prop(propSd)
    Cc_noxvori ~ prop(propSd_noxvori)
  })
}
