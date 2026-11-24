Daryani_2016_fluorouracil_scaled_adult <- function() {
  description <- paste(
    "Translational simulation model (not fitted) used by Daryani 2016 to",
    "choose the pediatric phase I dose of intravenous bolus 5-fluorouracil",
    "(5-FU) for recurrent ependymoma. The adult two-compartment",
    "Michaelis-Menten plasma model of Terret et al. (2000), with Vmax and the",
    "central volume normalised to body surface area, drives the mouse tumor",
    "extracellular fluid (tECF) model and the mouse tumor-growth inhibition",
    "model of Daryani 2016 (fraction unbound changed to the human value 0.9).",
    "Typical-value, deterministic model: no IIV and no residual error. Tumor",
    "growth is the fold change in a bioluminescence-equivalent tumor burden",
    "from its value at the first dose. Dose in mg; concentrations in uM."
  )
  reference <- paste(
    "Daryani VM, Patel YT, Tagen M, Turner DC, Carcaboso AM, Atkinson JM,",
    "Gajjar A, Gilbertson RJ, Wright KD, Stewart CF. Translational",
    "Pharmacokinetic-Pharmacodynamic Modeling and Simulation: Optimizing",
    "5-Fluorouracil Dosing in Children With Pediatric Ependymoma.",
    "CPT Pharmacometrics Syst Pharmacol. 2016;5(4):211-221.",
    "doi:10.1002/psp4.12075.",
    "Plasma PK: Terret C, et al. Dose and time",
    "dependencies of 5-fluorouracil pharmacokinetics. Clin Pharmacol Ther.",
    "2000;68:270-279, doi:10.1067/mcp.2000.109352, as BSA-normalised in Daryani 2016 Supplementary",
    "Table 1."
  )
  vignette <- "Daryani_2016_fluorouracil"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "umol/L (plasma and tECF); tumor burden is a fold change relative to time zero (unitless)"
  )

  compartmentData <- list(
    central = list(analyte = "5-fluorouracil", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "5-fluorouracil", units = "mg", specimen = "tissue", verified = TRUE),
    ecf = list(
      analyte = "5-fluorouracil (unbound, tumor extracellular fluid; predicted from the mouse model)",
      units = "mg",
      specimen = "tumor",
      verified = TRUE
    ),
    int_tumor = list(
      analyte = "5-fluorouracil (intracellular tumor; predicted from the mouse model)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    cycling_cells = list(
      analyte = "proliferating tumor cells (bioluminescence-equivalent fold change)",
      units = "unitless",
      specimen = "tumor",
      verified = TRUE
    ),
    damaged_cells1 = list(
      analyte = "drug-damaged tumor cells, transit 1",
      units = "unitless",
      specimen = "tumor",
      verified = TRUE
    ),
    damaged_cells2 = list(
      analyte = "drug-damaged tumor cells, transit 2",
      units = "unitless",
      specimen = "tumor",
      verified = TRUE
    ),
    damaged_cells3 = list(
      analyte = "drug-damaged tumor cells, transit 3",
      units = "unitless",
      specimen = "tumor",
      verified = TRUE
    )
  )

  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NA,
      notes = paste(
        "Vmax, the central volume and the tECF volume and Vmax are linear in",
        "BSA (Supplementary Table 1 per-m^2 values times BSA). With dosing in",
        "mg/m^2 x BSA every concentration is independent of BSA."
      ),
      source_name = "BSA"
    ),
    STUDY_EXP23 = list(
      description = "1 = use the exponential tumor growth rate of mouse efficacy experiments two and three; 0 = that of experiment one",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Carried over from the mouse tumor-growth model (Daryani 2016",
        "Equations 7-8). The paper does not say which experiment's growth rate",
        "its pediatric simulations (Figure 4d-f) used; 0 is the default."
      ),
      source_name = "EXPER"
    )
  )

  population <- list(
    species = "human (simulated pediatric patients; plasma PK from adults)",
    n_subjects = "Not fitted: simulation model. Terret 2000 adult plasma PK, Daryani 2016 mouse tumor PK-PD.",
    n_studies = 0L,
    disease_state = "Recurrent pediatric ependymoma (simulation target population)",
    dose_range = "400, 500 and 650 mg/m^2 5-FU intravenous bolus weekly x 4 (Figure 4)",
    regions = "USA (St. Jude Children's Research Hospital)",
    notes = paste(
      "Qualified by simulating a 12-hour 80 mg/m^2/h 5-FU infusion in",
      "children (mean steady-state 4.1 uM versus 6.7 uM observed by",
      "Daryani 2016 ref 11). The tumor PK and PD parameters are the mouse",
      "estimates held at their values; only the plasma fraction unbound was",
      "changed from 0.366 (mouse) to 0.9 (human)."
    )
  )

  ini({
    # Terret 2000 adult plasma PK, BSA-normalised (Supplementary Table 1).
    lvmax <- fixed(log(5773)); label("Maximum plasma elimination rate per m^2 BSA (umol/h/m^2)") # Supplementary Table 1 Vmax,ap scaled 5773 (Terret 751 mg/h/m^2)
    lkm <- fixed(log(42.8)); label("Michaelis-Menten constant for plasma elimination (umol/L)") # Supplementary Table 1 Km,ap scaled 42.8 (Terret 5.57 mg/L)
    lvc <- fixed(log(7.4)); label("Central volume per m^2 BSA (L/m^2)") # Supplementary Table 1 V1,ap scaled 7.4 (Terret 12.7 L)
    lk12 <- fixed(log(5.35)); label("Central-to-peripheral rate constant (1/h)") # Supplementary Table 1 K12,ap = 5.35 (Terret)
    lk21 <- fixed(log(5.69)); label("Peripheral-to-central rate constant (1/h)") # Supplementary Table 1 K21,ap = 5.69 (Terret)

    # Tumor ECF link (mouse Table 1 values as held in Supplementary Material 9).
    fu <- fixed(0.9); label("Fraction unbound in human plasma (fraction)") # Results 'Application of scaled-adult model': 5-FU plasma fraction unbound 0.9 (ref 23)
    lvecf <- fixed(log(0.001)); label("Tumor ECF volume per m^2 BSA (L/m^2)") # Table 1 V3,m = 0.001, used on the same per-m^2 scale as V1
    lkin_ecf <- fixed(log(0.00426)); label("Rate constant from plasma (unbound) into tumor ECF (1/h)") # Supplementary Material 9 K13 = 0.00426; Table 1 K13,m = 0.0043
    lkout_ecf <- fixed(log(3.24)); label("Rate constant from tumor ECF back to plasma (1/h)") # Table 1 K31,m = 3.24
    lkin_int_tumor <- fixed(log(1.82)); label("Rate constant from tumor ECF into intracellular tumor (1/h)") # Table 1 K34,m = 1.82
    lkout_int_tumor <- fixed(log(0.334)); label("Rate constant from intracellular tumor back to tumor ECF (1/h)") # Table 1 K43,m = 0.334
    lvmax_ecf <- fixed(log(0.00633)); label("Maximum elimination rate from tumor ECF per m^2 BSA (umol/h/m^2)") # Supplementary Material 9 K30 = 0.00633; Table 1 Vmax,tumor,m = 0.0063
    lkm_ecf <- fixed(log(0.0107)); label("Michaelis-Menten constant for tumor ECF elimination (umol/L)") # Supplementary Material 9 KM2 = 0.0107; Table 1 prints 0.012

    # Mouse tumor growth and 5-FU effect (Table 1 'Efficacy study').
    lkgrow <- fixed(log(0.158)); label("Exponential tumor growth rate, mouse experiment one (1/h)") # Table 1 theta_Kexp,m,1 = 0.158
    e_study_exp23_kgrow <- fixed(0.102); label("Additive change in exponential growth rate for mouse experiments two and three (1/h)") # Table 1 theta_Kexp,m,2/3 = 0.102
    lkgrow_lin <- fixed(log(214000)); label("Linear tumor growth rate (1/h)") # Table 1 Klin,m = 214,000
    psi <- fixed(0.0936); label("Exponent switching tumor growth from exponential to linear (unitless)") # Table 1 Psi_m = 0.0936
    lkmax <- fixed(log(0.658)); label("Maximum 5-FU kill rate of proliferating tumor cells (1/h)") # Table 1 Kmax,m = 0.658
    lec50 <- fixed(log(2.12)); label("Tumor ECF 5-FU concentration giving half-maximal kill rate (umol/L)") # Table 1 IC50,m = 2.12
    hill <- fixed(1.22); label("Hill coefficient of the kill-rate Emax model (unitless)") # Table 1 Hm = 1.22
    lktr <- fixed(log(0.946)); label("Transit rate constant through the damaged-cell compartments (1/h)") # Table 1 Kdel,m = 0.946
  })

  model({
    # 5-FU molecular weight (g/mol); Supplementary Table 1 converts
    # 751 mg/h/m^2 to 5773 umol/h/m^2 with this factor.
    mw <- 130.08

    vmax <- exp(lvmax) * BSA
    km <- exp(lkm)
    vc <- exp(lvc) * BSA
    k12 <- exp(lk12)
    k21 <- exp(lk21)

    vecf <- exp(lvecf) * BSA
    kin_ecf <- exp(lkin_ecf)
    kout_ecf <- exp(lkout_ecf)
    kin_int_tumor <- exp(lkin_int_tumor)
    kout_int_tumor <- exp(lkout_int_tumor)
    vmax_ecf <- exp(lvmax_ecf) * BSA
    km_ecf <- exp(lkm_ecf)

    kgrow <- exp(lkgrow) + e_study_exp23_kgrow * STUDY_EXP23
    kgrow_lin <- exp(lkgrow_lin)
    kmax <- exp(lkmax)
    ec50 <- exp(lec50)
    ktr <- exp(lktr)

    Cc <- central / vc * 1000 / mw
    Cecf <- ecf / vecf * 1000 / mw

    Cecf_pd <- Cecf
    if (Cecf_pd < 0) {
      Cecf_pd <- 0
    }
    kkill <- kmax * Cecf_pd^hill / (ec50^hill + Cecf_pd^hill)
    growth_rate <- kgrow * cycling_cells / (1 + (kgrow / kgrow_lin * cycling_cells)^psi)^(1 / psi)

    d/dt(central) <- -vmax * mw / 1000 * Cc / (km + Cc) - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(ecf) <- kin_ecf * fu * central - kout_ecf * ecf - kin_int_tumor * ecf + kout_int_tumor * int_tumor -
      vmax_ecf * mw / 1000 * Cecf / (km_ecf + Cecf)
    d/dt(int_tumor) <- kin_int_tumor * ecf - kout_int_tumor * int_tumor

    d/dt(cycling_cells) <- growth_rate - kkill * cycling_cells
    d/dt(damaged_cells1) <- kkill * cycling_cells - ktr * damaged_cells1
    d/dt(damaged_cells2) <- ktr * damaged_cells1 - ktr * damaged_cells2
    d/dt(damaged_cells3) <- ktr * damaged_cells2 - ktr * damaged_cells3
    cycling_cells(0) <- 1

    tumor_size <- cycling_cells + damaged_cells1 + damaged_cells2 + damaged_cells3
  })
}
