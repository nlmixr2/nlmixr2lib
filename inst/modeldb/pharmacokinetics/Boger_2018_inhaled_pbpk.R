Boger_2018_inhaled_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body inhalation, partial-differential-equation lung, 169",
    "ODEs; bespoke MATLAB model deposited as Appendix S2). The first",
    "inhalation PBPK model to treat the lung as a continuous heterogeneous",
    "organ and the inhaled powder as a distribution that is continuous in",
    "both lung depth and particle size. Regional deposition (typical-path",
    "Findeisen-Landahl model over a Weibel airway tree, plus Cheng",
    "extrathoracic deposition) sets the initial particle field; particles",
    "then dissolve by the Nernst-Brunner equation into the epithelial",
    "lining fluid, are swept proximally by mucociliary clearance and",
    "swallowed, or permeate epithelium and sub-epithelium into the",
    "bronchial or pulmonary circulation of a seven-tissue whole-body PBPK",
    "model. The spatial coordinate is discretised by the method of lines",
    "into 12 tracheobronchial and 8 alveolar slabs, and the particle size",
    "coordinate into 8 size classes carrying a conserved number density and",
    "z-moment each. There is no drug: the compound is a hypothetical",
    "neutral small molecule (MW 250) relying on transcellular transport,",
    "used to demonstrate the framework in three case studies (mucociliary",
    "clearance, particle size distribution, overdosing). Deterministic",
    "typical-value simulation model: the paper reports no IIV and no",
    "residual-error model."
  )
  reference <- paste(
    "Boger E, Wigstrom O. A Partial Differential Equation Approach to",
    "Inhalation Physiologically Based Pharmacokinetic Modeling. CPT",
    "Pharmacometrics Syst Pharmacol. 2018;7(10):638-646.",
    "doi:10.1002/psp4.12344.",
    sep = " "
  )
  vignette <- "Boger_2018_inhaled_pbpk"
  units <- list(time = "h", dosing = "nmol", concentration = "nM")

  # No covariates: the paper simulates a single 70 kg reference adult and
  # varies only the solubility, particle size distribution and lung-deposited
  # dose between case studies (Table 1).
  covariateData <- list()

  compartmentData <- list(
    a_spleen = list(analyte = "generic neutral small molecule", units = "nmol", specimen = "tissue", verified = TRUE),
    a_rapidly_perfused = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    a_slowly_perfused = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    a_fat = list(analyte = "generic neutral small molecule", units = "nmol", specimen = "tissue", verified = TRUE),
    a_liver = list(analyte = "generic neutral small molecule", units = "nmol", specimen = "tissue", verified = TRUE),
    a_gut = list(analyte = "generic neutral small molecule", units = "nmol", specimen = "tissue", verified = TRUE),
    a_arterial = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "whole blood",
      verified = TRUE
    ),
    a_venous = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "whole blood",
      verified = TRUE
    ),
    depot = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    elf_tb_slab1 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    epithelium_tb_slab1 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    subepithelium_tb_slab1 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    elf_tb_slab2 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    epithelium_tb_slab2 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    subepithelium_tb_slab2 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    elf_tb_slab3 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    epithelium_tb_slab3 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    subepithelium_tb_slab3 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    elf_tb_slab4 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    epithelium_tb_slab4 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    subepithelium_tb_slab4 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    elf_tb_slab5 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    epithelium_tb_slab5 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    subepithelium_tb_slab5 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    elf_tb_slab6 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    epithelium_tb_slab6 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    subepithelium_tb_slab6 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    elf_tb_slab7 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    epithelium_tb_slab7 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    subepithelium_tb_slab7 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    elf_tb_slab8 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    epithelium_tb_slab8 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    subepithelium_tb_slab8 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    elf_tb_slab9 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    epithelium_tb_slab9 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    subepithelium_tb_slab9 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    elf_tb_slab10 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    epithelium_tb_slab10 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    subepithelium_tb_slab10 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    elf_tb_slab11 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    epithelium_tb_slab11 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    subepithelium_tb_slab11 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    elf_tb_slab12 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    epithelium_tb_slab12 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    subepithelium_tb_slab12 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    elf_alv_slab1 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    epithelium_alv_slab1 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    subepithelium_alv_slab1 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    elf_alv_slab2 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    epithelium_alv_slab2 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    subepithelium_alv_slab2 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    elf_alv_slab3 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    epithelium_alv_slab3 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    subepithelium_alv_slab3 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    elf_alv_slab4 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    epithelium_alv_slab4 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    subepithelium_alv_slab4 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    elf_alv_slab5 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    epithelium_alv_slab5 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    subepithelium_alv_slab5 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    elf_alv_slab6 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    epithelium_alv_slab6 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    subepithelium_alv_slab6 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    elf_alv_slab7 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    epithelium_alv_slab7 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    subepithelium_alv_slab7 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    elf_alv_slab8 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    epithelium_alv_slab8 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    subepithelium_alv_slab8 = list(
      analyte = "generic neutral small molecule",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    particles_tb_a_slab1 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_a_slab1 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_b_slab1 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_b_slab1 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_c_slab1 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_c_slab1 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_d_slab1 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_d_slab1 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_e_slab1 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_e_slab1 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_f_slab1 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_f_slab1 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_g_slab1 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_g_slab1 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_h_slab1 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_h_slab1 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_a_slab2 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_a_slab2 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_b_slab2 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_b_slab2 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_c_slab2 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_c_slab2 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_d_slab2 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_d_slab2 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_e_slab2 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_e_slab2 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_f_slab2 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_f_slab2 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_g_slab2 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_g_slab2 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_h_slab2 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_h_slab2 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_a_slab3 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_a_slab3 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_b_slab3 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_b_slab3 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_c_slab3 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_c_slab3 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_d_slab3 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_d_slab3 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_e_slab3 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_e_slab3 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_f_slab3 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_f_slab3 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_g_slab3 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_g_slab3 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_h_slab3 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_h_slab3 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_a_slab4 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_a_slab4 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_b_slab4 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_b_slab4 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_c_slab4 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_c_slab4 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_d_slab4 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_d_slab4 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_e_slab4 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_e_slab4 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_f_slab4 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_f_slab4 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_g_slab4 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_g_slab4 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_h_slab4 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_h_slab4 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_a_slab5 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_a_slab5 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_b_slab5 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_b_slab5 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_c_slab5 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_c_slab5 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_d_slab5 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_d_slab5 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_e_slab5 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_e_slab5 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_f_slab5 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_f_slab5 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_g_slab5 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_g_slab5 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_h_slab5 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_h_slab5 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_a_slab6 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_a_slab6 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_b_slab6 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_b_slab6 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_c_slab6 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_c_slab6 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_d_slab6 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_d_slab6 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_e_slab6 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_e_slab6 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_f_slab6 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_f_slab6 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_g_slab6 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_g_slab6 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_h_slab6 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_h_slab6 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_a_slab7 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_a_slab7 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_b_slab7 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_b_slab7 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_c_slab7 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_c_slab7 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_d_slab7 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_d_slab7 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_e_slab7 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_e_slab7 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_f_slab7 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_f_slab7 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_g_slab7 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_g_slab7 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_h_slab7 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_h_slab7 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_a_slab8 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_a_slab8 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_b_slab8 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_b_slab8 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_c_slab8 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_c_slab8 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_d_slab8 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_d_slab8 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_e_slab8 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_e_slab8 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_f_slab8 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_f_slab8 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_g_slab8 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_g_slab8 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_h_slab8 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_h_slab8 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_a_slab9 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_a_slab9 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_b_slab9 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_b_slab9 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_c_slab9 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_c_slab9 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_d_slab9 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_d_slab9 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_e_slab9 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_e_slab9 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_f_slab9 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_f_slab9 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_g_slab9 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_g_slab9 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_h_slab9 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_h_slab9 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_a_slab10 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_a_slab10 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_b_slab10 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_b_slab10 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_c_slab10 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_c_slab10 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_d_slab10 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_d_slab10 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_e_slab10 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_e_slab10 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_f_slab10 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_f_slab10 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_g_slab10 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_g_slab10 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_h_slab10 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_h_slab10 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_a_slab11 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_a_slab11 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_b_slab11 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_b_slab11 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_c_slab11 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_c_slab11 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_d_slab11 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_d_slab11 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_e_slab11 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_e_slab11 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_f_slab11 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_f_slab11 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_g_slab11 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_g_slab11 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_h_slab11 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_h_slab11 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_a_slab12 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_a_slab12 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_b_slab12 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_b_slab12 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_c_slab12 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_c_slab12 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_d_slab12 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_d_slab12 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_e_slab12 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_e_slab12 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_f_slab12 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_f_slab12 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_g_slab12 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_g_slab12 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_tb_h_slab12 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_tb_h_slab12 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_a_slab1 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_a_slab1 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_b_slab1 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_b_slab1 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_c_slab1 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_c_slab1 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_d_slab1 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_d_slab1 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_e_slab1 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_e_slab1 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_f_slab1 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_f_slab1 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_g_slab1 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_g_slab1 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_h_slab1 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_h_slab1 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_a_slab2 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_a_slab2 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_b_slab2 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_b_slab2 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_c_slab2 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_c_slab2 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_d_slab2 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_d_slab2 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_e_slab2 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_e_slab2 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_f_slab2 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_f_slab2 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_g_slab2 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_g_slab2 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_h_slab2 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_h_slab2 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_a_slab3 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_a_slab3 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_b_slab3 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_b_slab3 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_c_slab3 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_c_slab3 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_d_slab3 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_d_slab3 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_e_slab3 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_e_slab3 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_f_slab3 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_f_slab3 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_g_slab3 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_g_slab3 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_h_slab3 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_h_slab3 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_a_slab4 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_a_slab4 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_b_slab4 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_b_slab4 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_c_slab4 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_c_slab4 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_d_slab4 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_d_slab4 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_e_slab4 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_e_slab4 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_f_slab4 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_f_slab4 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_g_slab4 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_g_slab4 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_h_slab4 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_h_slab4 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_a_slab5 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_a_slab5 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_b_slab5 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_b_slab5 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_c_slab5 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_c_slab5 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_d_slab5 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_d_slab5 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_e_slab5 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_e_slab5 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_f_slab5 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_f_slab5 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_g_slab5 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_g_slab5 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_h_slab5 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_h_slab5 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_a_slab6 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_a_slab6 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_b_slab6 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_b_slab6 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_c_slab6 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_c_slab6 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_d_slab6 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_d_slab6 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_e_slab6 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_e_slab6 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_f_slab6 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_f_slab6 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_g_slab6 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_g_slab6 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_h_slab6 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_h_slab6 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_a_slab7 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_a_slab7 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_b_slab7 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_b_slab7 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_c_slab7 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_c_slab7 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_d_slab7 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_d_slab7 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_e_slab7 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_e_slab7 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_f_slab7 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_f_slab7 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_g_slab7 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_g_slab7 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_h_slab7 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_h_slab7 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_a_slab8 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_a_slab8 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_b_slab8 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_b_slab8 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_c_slab8 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_c_slab8 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_d_slab8 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_d_slab8 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_e_slab8 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_e_slab8 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_f_slab8 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_f_slab8 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_g_slab8 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_g_slab8 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    particles_alv_h_slab8 = list(
      analyte = "undissolved drug particles",
      units = "particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    pmass_alv_h_slab8 = list(
      analyte = "undissolved drug particles",
      units = "dm^2 particles per dm lung depth",
      specimen = "epithelial lining fluid",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 0,
    notes = paste(
      "No subjects: a theoretical modelling exercise, not a fit to data.",
      "System physiology is a single 70 kg reference adult with a 5.2",
      "L/minute cardiac output; tissue volumes and blood flows are from",
      "Brown et al. and Bernareggi and Rowland (Table S1). The airway tree",
      "is Weibel symmetric-branching geometry scaled to a functional",
      "residual capacity of 3,000 mL (Table S3), with epithelial lining",
      "fluid, epithelium and sub-epithelium layer heights from ICRP66,",
      "Patton and Byron, and Mariassy (Table S2). Deposition assumes tidal",
      "breathing with a 1,000 mL tidal volume, 15 breaths/minute, a 1:1",
      "inspiratory:expiratory ratio and no breath hold. The compound is a",
      "hypothetical neutral small molecule (MW 250, Table 2).",
      sep = " "
    )
  )

  ini({
    # ---------------- drug-specific parameters (Table 2) ----------------
    lcl <- fixed(log(         70)); label("Log blood clearance (L/h)")                        # Table 2, CL_b = 70 L/hour (1 L/h/kg x 70 kg)
    bp <- fixed(1); label("Blood/plasma ratio")                                      # Table 2, blood/plasma ratio = 1
    fup <- fixed(0.75); label("Unbound fraction in plasma")                          # Table 2, f_u,p = 0.75
    fuelf <- fixed(1); label("Unbound fraction in the epithelial lining fluid")      # Table 2, f_u,fluid = 1
    foral <- fixed(       0.05); label("Oral bioavailability of the swallowed fraction")      # Appendix S2 +drug/loadData.m (data.F); Table 2 prints F = 0.20
    lka <- fixed(log(0.6)); label("Log oral absorption rate constant (1/h)")          # Table 2, k_a = 0.6 hour-1
    cs <- fixed(100); label("Aqueous solubility (nM)")                                # Table 1, case study 1 and 3: C_s = 100 nM (case study 2 uses 250 nM)
    peff <- fixed(0.01510974984); label("Effective airway permeability (dm/h)")                  # Table 2 P_app = 1.5e-6 cm/s, converted per Sjogren 2013 (Appendix S2 +drug/CalcPeffFromPapp.m)
    vdiff <- fixed(0.0003054532316); label("Aqueous diffusion coefficient (dm^2/h)")               # Table 2, v_diff = 8.5e-6 cm^2/s (Stokes-Einstein, MW 250)
    dens <- fixed( 4000000000); label("Particle density (nmol/dm^3)")                          # Table 2, rho = 1 g/cm^3 and MW = 250 g/mol
    vrr <- fixed(1); label("Stagnant-layer thickness relative to particle radius")    # Appendix S2 +drug/loadData.m (data.dr = 1)
    rmu <- fixed(    1.5e-05); label("Mean particle radius (dm)")                              # Case study 1/3 wide PSD: d ~ N(3, 0.6) um, so r ~ N(1.5, 0.3) um
    rsd <- fixed(      3e-06); label("SD of particle radius (dm)")                             # Case study 1/3 wide PSD (Results, case study 2 definition)
    mccvel <- fixed(       2.16); label("Mucociliary clearance velocity at the trachea (dm/h)") # Eq. 2, alpha_0 = 3.6 mm/minute (Yeates 1975)

    # ------- tissue-to-plasma partition coefficients (Table 2) -------
    kpspleen <- fixed( 3.88275569); label("Kp, spleen")                                        # Table 2, K_p,spleen = 3.9
    kprapid <- fixed(3.075006329); label("Kp, richly perfused tissue")                          # Table 2, K_p,richly = 3.1
    kpslow <- fixed(2.316650947); label("Kp, poorly perfused tissue")                           # Table 2, K_p,poorly = 2.3
    kpfat <- fixed(0.5558844999); label("Kp, adipose")                                           # Table 2, K_p,adipose = 0.56
    kpliver <- fixed(5.926939171); label("Kp, hepatic")                                         # Table 2, K_p,hepatic = 5.9
    kpgut <- fixed( 3.69975855); label("Kp, gut")                                               # Table 2, K_p,gut = 3.7
    kpulung <- fixed(6.5); label("Unbound tissue-plasma partition coefficient, lung")  # Table 2, K_p,u,lung = 6.5

    # ---------------- physiology (Table S1, 70 kg adult) ----------------
    vspleen <- fixed(      0.182); label("Spleen volume (L)")                                   # Table S1, 0.0026 x 70 kg
    vrapid <- fixed(      2.037); label("Richly perfused tissue volume (L)")                    # Table S1, 0.0291 x 70 kg
    vslow <- fixed(     43.862); label("Poorly perfused tissue volume (L)")                     # Table S1, 0.6266 x 70 kg
    vfat <- fixed(     14.994); label("Adipose volume (L)")                                     # Table S1, 0.2142 x 70 kg
    vliver <- fixed(      1.799); label("Liver volume (L)")                                     # Table S1, 0.0257 x 70 kg
    vgut <- fixed(      1.197); label("Gut volume (L)")                                         # Table S1, 0.0171 x 70 kg
    varterial <- fixed(      1.799); label("Arterial blood volume (L)")                         # Table S1, 0.0257 x 70 kg
    vvenous <- fixed(      3.598); label("Venous blood volume (L)")                             # Table S1, 0.0514 x 70 kg
    vlung <- fixed(      0.532); label("Total lung volume (L)")                                 # Table S1, 0.0076 x 70 kg
    qco <- fixed(        312); label("Cardiac output (L/h)")                                     # Table S1, 5.2 L/minute
    qspleen <- fixed(       6.24); label("Spleen blood flow (L/h)")                              # Table S1, 0.020 x cardiac output
    qrapid <- fixed(      109.2); label("Richly perfused tissue blood flow (L/h)")               # Table S1, 0.35 x cardiac output
    qslow <- fixed(      101.4); label("Poorly perfused tissue blood flow (L/h)")                # Table S1, 0.325 x cardiac output
    qfat <- fixed(       15.6); label("Adipose blood flow (L/h)")                                # Table S1, 0.050 x cardiac output
    qliver <- fixed(      18.72); label("Hepatic arterial blood flow (L/h)")                     # Table S1, 0.060 x cardiac output
    qgut <- fixed(      53.04); label("Gut blood flow (L/h)")                                    # Table S1, 0.17 x cardiac output
    qbr <- fixed(        7.8); label("Total bronchial blood flow (L/h)")                         # Table S1, 0.025 x cardiac output
  })

  model({
    # ---- derived quantities ----
    cl <- exp(lcl)
    ka <- exp(lka)
    # Kp,lung from the unbound partition coefficient and the plasma
    # unbound fraction (Appendix S2, +drug/loadData.m).
    kplung <- kpulung * fup
    # Nernst-Brunner shrinkage rate in squared radius (Eq. 10):
    # dz/dt = -2 * vdiff / (dens * vrr) * (cs - fuelf * C_elf).
    kdiss <- 2 * vdiff / (dens * vrr)
    # Particle-mass coefficient of Eq. S1: mass per unit depth =
    # pmc * sum_j pmass_j.
    pmc <- dens * 2 / 3 * pi

    # ---- particle size classes: K equal-width radius bins ----
    hr <- 7 * rsd / 8
    rba <- rmu + ( -3.5) * hr; hza <- 2 * rba * hr
    rbb <- rmu + ( -2.5) * hr; hzb <- 2 * rbb * hr
    rbc <- rmu + ( -1.5) * hr; hzc <- 2 * rbc * hr
    rbd <- rmu + ( -0.5) * hr; hzd <- 2 * rbd * hr
    rbe <- rmu + (  0.5) * hr; hze <- 2 * rbe * hr
    rbf <- rmu + (  1.5) * hr; hzf <- 2 * rbf * hr
    rbg <- rmu + (  2.5) * hr; hzg <- 2 * rbg * hr
    rbh <- rmu + (  3.5) * hr; hzh <- 2 * rbh * hr

    # ---- systemic concentrations ----
    cart <- a_arterial / varterial
    cven <- a_venous / vvenous
    csp <- a_spleen / vspleen
    cri <- a_rapidly_perfused / vrapid
    cpo <- a_slowly_perfused / vslow
    cad <- a_fat / vfat
    chep <- a_liver / vliver
    cgut <- a_gut / vgut

    # ================= tracheobronchial region =================
    # -- tb slab 1 (x =       0 dm) --
    vuf1 <- 6.9742109e-06; vup1 <- 3.505846516e-05; vus1 <- 0.003864809353
    cuf1 <- elf_tb_slab1 / vuf1
    cup1 <- epithelium_tb_slab1 / vup1
    cus1 <- subepithelium_tb_slab1 / vus1
    que1 <- peff * 0.1269173741; qus1 <- peff * 0.1280536963
    qu1 <- qbr * 0.00262906147
    gamu1 <- -kdiss * (cs - fuelf * cuf1)
    wu1a <- min(particles_tb_a_slab1, max(0, pmass_tb_a_slab1 / hza))
    wu1b <- min(particles_tb_b_slab1, max(0, pmass_tb_b_slab1 / hzb))
    wu1c <- min(particles_tb_c_slab1, max(0, pmass_tb_c_slab1 / hzc))
    wu1d <- min(particles_tb_d_slab1, max(0, pmass_tb_d_slab1 / hzd))
    wu1e <- min(particles_tb_e_slab1, max(0, pmass_tb_e_slab1 / hze))
    wu1f <- min(particles_tb_f_slab1, max(0, pmass_tb_f_slab1 / hzf))
    wu1g <- min(particles_tb_g_slab1, max(0, pmass_tb_g_slab1 / hzg))
    wu1h <- min(particles_tb_h_slab1, max(0, pmass_tb_h_slab1 / hzh))
    yintu1 <- wu1a + wu1b + wu1c + wu1d + wu1e + wu1f + wu1g + wu1h
    pmu1 <- pmc * (pmass_tb_a_slab1 +
      pmass_tb_b_slab1 +
      pmass_tb_c_slab1 +
      pmass_tb_d_slab1 +
      pmass_tb_e_slab1 +
      pmass_tb_f_slab1 +
      pmass_tb_g_slab1 +
      pmass_tb_h_slab1)
    disu1 <- -pmc * gamu1 * yintu1 * 0.1025204545

    # -- tb slab 2 (x = 0.205041 dm) --
    vuf2 <- 1.19325047e-05; vup2 <- 5.99458921e-05; vus2 <- 0.006574800996
    cuf2 <- elf_tb_slab2 / vuf2
    cup2 <- epithelium_tb_slab2 / vup2
    cus2 <- subepithelium_tb_slab2 / vus2
    que2 <- peff * 0.2171263438; qus2 <- peff * 0.2188440296
    qu2 <- qbr * 0.004464114264
    gamu2 <- -kdiss * (cs - fuelf * cuf2)
    wu2a <- min(particles_tb_a_slab2, max(0, pmass_tb_a_slab2 / hza))
    wu2b <- min(particles_tb_b_slab2, max(0, pmass_tb_b_slab2 / hzb))
    wu2c <- min(particles_tb_c_slab2, max(0, pmass_tb_c_slab2 / hzc))
    wu2d <- min(particles_tb_d_slab2, max(0, pmass_tb_d_slab2 / hzd))
    wu2e <- min(particles_tb_e_slab2, max(0, pmass_tb_e_slab2 / hze))
    wu2f <- min(particles_tb_f_slab2, max(0, pmass_tb_f_slab2 / hzf))
    wu2g <- min(particles_tb_g_slab2, max(0, pmass_tb_g_slab2 / hzg))
    wu2h <- min(particles_tb_h_slab2, max(0, pmass_tb_h_slab2 / hzh))
    yintu2 <- wu2a + wu2b + wu2c + wu2d + wu2e + wu2f + wu2g + wu2h
    pmu2 <- pmc * (pmass_tb_a_slab2 +
      pmass_tb_b_slab2 +
      pmass_tb_c_slab2 +
      pmass_tb_d_slab2 +
      pmass_tb_e_slab2 +
      pmass_tb_f_slab2 +
      pmass_tb_g_slab2 +
      pmass_tb_h_slab2)
    disu2 <- -pmc * gamu2 * yintu2 * 0.2050409091

    # -- tb slab 3 (x = 0.410082 dm) --
    vuf3 <- 1.100557391e-05; vup3 <- 5.526708895e-05; vus3 <- 0.00604178535
    cuf3 <- elf_tb_slab3 / vuf3
    cup3 <- epithelium_tb_slab3 / vup3
    cus3 <- subepithelium_tb_slab3 / vus3
    que3 <- peff * 0.2002463229; qus3 <- peff * 0.2016961644
    qu3 <- qbr * 0.004103136302
    gamu3 <- -kdiss * (cs - fuelf * cuf3)
    wu3a <- min(particles_tb_a_slab3, max(0, pmass_tb_a_slab3 / hza))
    wu3b <- min(particles_tb_b_slab3, max(0, pmass_tb_b_slab3 / hzb))
    wu3c <- min(particles_tb_c_slab3, max(0, pmass_tb_c_slab3 / hzc))
    wu3d <- min(particles_tb_d_slab3, max(0, pmass_tb_d_slab3 / hzd))
    wu3e <- min(particles_tb_e_slab3, max(0, pmass_tb_e_slab3 / hze))
    wu3f <- min(particles_tb_f_slab3, max(0, pmass_tb_f_slab3 / hzf))
    wu3g <- min(particles_tb_g_slab3, max(0, pmass_tb_g_slab3 / hzg))
    wu3h <- min(particles_tb_h_slab3, max(0, pmass_tb_h_slab3 / hzh))
    yintu3 <- wu3a + wu3b + wu3c + wu3d + wu3e + wu3f + wu3g + wu3h
    pmu3 <- pmc * (pmass_tb_a_slab3 +
      pmass_tb_b_slab3 +
      pmass_tb_c_slab3 +
      pmass_tb_d_slab3 +
      pmass_tb_e_slab3 +
      pmass_tb_f_slab3 +
      pmass_tb_g_slab3 +
      pmass_tb_h_slab3)
    disu3 <- -pmc * gamu3 * yintu3 * 0.2050409091

    # -- tb slab 4 (x = 0.615123 dm) --
    vuf4 <- 1.099684268e-05; vup4 <- 5.522321019e-05; vus4 <- 0.006036950921
    cuf4 <- elf_tb_slab4 / vuf4
    cup4 <- epithelium_tb_slab4 / vup4
    cus4 <- subepithelium_tb_slab4 / vus4
    que4 <- peff * 0.200087439; qus4 <- peff * 0.2015359242
    qu4 <- qbr * 0.004108724065
    gamu4 <- -kdiss * (cs - fuelf * cuf4)
    wu4a <- min(particles_tb_a_slab4, max(0, pmass_tb_a_slab4 / hza))
    wu4b <- min(particles_tb_b_slab4, max(0, pmass_tb_b_slab4 / hzb))
    wu4c <- min(particles_tb_c_slab4, max(0, pmass_tb_c_slab4 / hzc))
    wu4d <- min(particles_tb_d_slab4, max(0, pmass_tb_d_slab4 / hzd))
    wu4e <- min(particles_tb_e_slab4, max(0, pmass_tb_e_slab4 / hze))
    wu4f <- min(particles_tb_f_slab4, max(0, pmass_tb_f_slab4 / hzf))
    wu4g <- min(particles_tb_g_slab4, max(0, pmass_tb_g_slab4 / hzg))
    wu4h <- min(particles_tb_h_slab4, max(0, pmass_tb_h_slab4 / hzh))
    yintu4 <- wu4a + wu4b + wu4c + wu4d + wu4e + wu4f + wu4g + wu4h
    pmu4 <- pmc * (pmass_tb_a_slab4 +
      pmass_tb_b_slab4 +
      pmass_tb_c_slab4 +
      pmass_tb_d_slab4 +
      pmass_tb_e_slab4 +
      pmass_tb_f_slab4 +
      pmass_tb_g_slab4 +
      pmass_tb_h_slab4)
    disu4 <- -pmc * gamu4 * yintu4 * 0.2050409091

    # -- tb slab 5 (x = 0.820164 dm) --
    vuf5 <- 1.173552429e-05; vup5 <- 5.895694517e-05; vus5 <- 0.006466676859
    cuf5 <- elf_tb_slab5 / vuf5
    cup5 <- epithelium_tb_slab5 / vup5
    cus5 <- subepithelium_tb_slab5 / vus5
    que5 <- peff * 0.2135424458; qus5 <- peff * 0.2152354362
    qu5 <- qbr * 0.004427560208
    gamu5 <- -kdiss * (cs - fuelf * cuf5)
    wu5a <- min(particles_tb_a_slab5, max(0, pmass_tb_a_slab5 / hza))
    wu5b <- min(particles_tb_b_slab5, max(0, pmass_tb_b_slab5 / hzb))
    wu5c <- min(particles_tb_c_slab5, max(0, pmass_tb_c_slab5 / hzc))
    wu5d <- min(particles_tb_d_slab5, max(0, pmass_tb_d_slab5 / hzd))
    wu5e <- min(particles_tb_e_slab5, max(0, pmass_tb_e_slab5 / hze))
    wu5f <- min(particles_tb_f_slab5, max(0, pmass_tb_f_slab5 / hzf))
    wu5g <- min(particles_tb_g_slab5, max(0, pmass_tb_g_slab5 / hzg))
    wu5h <- min(particles_tb_h_slab5, max(0, pmass_tb_h_slab5 / hzh))
    yintu5 <- wu5a + wu5b + wu5c + wu5d + wu5e + wu5f + wu5g + wu5h
    pmu5 <- pmc * (pmass_tb_a_slab5 +
      pmass_tb_b_slab5 +
      pmass_tb_c_slab5 +
      pmass_tb_d_slab5 +
      pmass_tb_e_slab5 +
      pmass_tb_f_slab5 +
      pmass_tb_g_slab5 +
      pmass_tb_h_slab5)
    disu5 <- -pmc * gamu5 * yintu5 * 0.2050409091

    # -- tb slab 6 (x =  1.0252 dm) --
    vuf6 <- 1.3050832e-05; vup6 <- 6.561098322e-05; vus6 <- 0.007237342316
    cuf6 <- elf_tb_slab6 / vuf6
    cup6 <- epithelium_tb_slab6 / vup6
    cus6 <- subepithelium_tb_slab6 / vus6
    que6 <- peff * 0.2375040969; qus6 <- peff * 0.2396668275
    qu6 <- qbr * 0.005033401726
    gamu6 <- -kdiss * (cs - fuelf * cuf6)
    wu6a <- min(particles_tb_a_slab6, max(0, pmass_tb_a_slab6 / hza))
    wu6b <- min(particles_tb_b_slab6, max(0, pmass_tb_b_slab6 / hzb))
    wu6c <- min(particles_tb_c_slab6, max(0, pmass_tb_c_slab6 / hzc))
    wu6d <- min(particles_tb_d_slab6, max(0, pmass_tb_d_slab6 / hzd))
    wu6e <- min(particles_tb_e_slab6, max(0, pmass_tb_e_slab6 / hze))
    wu6f <- min(particles_tb_f_slab6, max(0, pmass_tb_f_slab6 / hzf))
    wu6g <- min(particles_tb_g_slab6, max(0, pmass_tb_g_slab6 / hzg))
    wu6h <- min(particles_tb_h_slab6, max(0, pmass_tb_h_slab6 / hzh))
    yintu6 <- wu6a + wu6b + wu6c + wu6d + wu6e + wu6f + wu6g + wu6h
    pmu6 <- pmc * (pmass_tb_a_slab6 +
      pmass_tb_b_slab6 +
      pmass_tb_c_slab6 +
      pmass_tb_d_slab6 +
      pmass_tb_e_slab6 +
      pmass_tb_f_slab6 +
      pmass_tb_g_slab6 +
      pmass_tb_h_slab6)
    disu6 <- -pmc * gamu6 * yintu6 * 0.2050409091

    # -- tb slab 7 (x = 1.23025 dm) --
    vuf7 <- 1.477199092e-05; vup7 <- 7.43280738e-05; vus7 <- 0.008255333613
    cuf7 <- elf_tb_slab7 / vuf7
    cup7 <- epithelium_tb_slab7 / vup7
    cus7 <- subepithelium_tb_slab7 / vus7
    que7 <- peff * 0.2688653618; qus7 <- peff * 0.2717024467
    qu7 <- qbr * 0.005982635903
    gamu7 <- -kdiss * (cs - fuelf * cuf7)
    wu7a <- min(particles_tb_a_slab7, max(0, pmass_tb_a_slab7 / hza))
    wu7b <- min(particles_tb_b_slab7, max(0, pmass_tb_b_slab7 / hzb))
    wu7c <- min(particles_tb_c_slab7, max(0, pmass_tb_c_slab7 / hzc))
    wu7d <- min(particles_tb_d_slab7, max(0, pmass_tb_d_slab7 / hzd))
    wu7e <- min(particles_tb_e_slab7, max(0, pmass_tb_e_slab7 / hze))
    wu7f <- min(particles_tb_f_slab7, max(0, pmass_tb_f_slab7 / hzf))
    wu7g <- min(particles_tb_g_slab7, max(0, pmass_tb_g_slab7 / hzg))
    wu7h <- min(particles_tb_h_slab7, max(0, pmass_tb_h_slab7 / hzh))
    yintu7 <- wu7a + wu7b + wu7c + wu7d + wu7e + wu7f + wu7g + wu7h
    pmu7 <- pmc * (pmass_tb_a_slab7 +
      pmass_tb_b_slab7 +
      pmass_tb_c_slab7 +
      pmass_tb_d_slab7 +
      pmass_tb_e_slab7 +
      pmass_tb_f_slab7 +
      pmass_tb_g_slab7 +
      pmass_tb_h_slab7)
    disu7 <- -pmc * gamu7 * yintu7 * 0.2050409091

    # -- tb slab 8 (x = 1.43529 dm) --
    vuf8 <- 1.795418281e-05; vup8 <- 9.049884083e-05; vus8 <- 0.01019217703
    cuf8 <- elf_tb_slab8 / vuf8
    cup8 <- epithelium_tb_slab8 / vup8
    cus8 <- subepithelium_tb_slab8 / vus8
    que8 <- peff * 0.3268808029; qus8 <- peff * 0.3312930846
    qu8 <- qbr * 0.008421930238
    gamu8 <- -kdiss * (cs - fuelf * cuf8)
    wu8a <- min(particles_tb_a_slab8, max(0, pmass_tb_a_slab8 / hza))
    wu8b <- min(particles_tb_b_slab8, max(0, pmass_tb_b_slab8 / hzb))
    wu8c <- min(particles_tb_c_slab8, max(0, pmass_tb_c_slab8 / hzc))
    wu8d <- min(particles_tb_d_slab8, max(0, pmass_tb_d_slab8 / hzd))
    wu8e <- min(particles_tb_e_slab8, max(0, pmass_tb_e_slab8 / hze))
    wu8f <- min(particles_tb_f_slab8, max(0, pmass_tb_f_slab8 / hzf))
    wu8g <- min(particles_tb_g_slab8, max(0, pmass_tb_g_slab8 / hzg))
    wu8h <- min(particles_tb_h_slab8, max(0, pmass_tb_h_slab8 / hzh))
    yintu8 <- wu8a + wu8b + wu8c + wu8d + wu8e + wu8f + wu8g + wu8h
    pmu8 <- pmc * (pmass_tb_a_slab8 +
      pmass_tb_b_slab8 +
      pmass_tb_c_slab8 +
      pmass_tb_d_slab8 +
      pmass_tb_e_slab8 +
      pmass_tb_f_slab8 +
      pmass_tb_g_slab8 +
      pmass_tb_h_slab8)
    disu8 <- -pmc * gamu8 * yintu8 * 0.2050409091

    # -- tb slab 9 (x = 1.64033 dm) --
    vuf9 <- 2.854957941e-05; vup9 <- 0.0001447915723; vus9 <- 0.01708541125
    cuf9 <- elf_tb_slab9 / vuf9
    cup9 <- epithelium_tb_slab9 / vup9
    cus9 <- subepithelium_tb_slab9 / vus9
    que9 <- peff * 0.5203218461; qus9 <- peff * 0.5327078374
    qu9 <- qbr * 0.02762984831
    gamu9 <- -kdiss * (cs - fuelf * cuf9)
    wu9a <- min(particles_tb_a_slab9, max(0, pmass_tb_a_slab9 / hza))
    wu9b <- min(particles_tb_b_slab9, max(0, pmass_tb_b_slab9 / hzb))
    wu9c <- min(particles_tb_c_slab9, max(0, pmass_tb_c_slab9 / hzc))
    wu9d <- min(particles_tb_d_slab9, max(0, pmass_tb_d_slab9 / hzd))
    wu9e <- min(particles_tb_e_slab9, max(0, pmass_tb_e_slab9 / hze))
    wu9f <- min(particles_tb_f_slab9, max(0, pmass_tb_f_slab9 / hzf))
    wu9g <- min(particles_tb_g_slab9, max(0, pmass_tb_g_slab9 / hzg))
    wu9h <- min(particles_tb_h_slab9, max(0, pmass_tb_h_slab9 / hzh))
    yintu9 <- wu9a + wu9b + wu9c + wu9d + wu9e + wu9f + wu9g + wu9h
    pmu9 <- pmc * (pmass_tb_a_slab9 +
      pmass_tb_b_slab9 +
      pmass_tb_c_slab9 +
      pmass_tb_d_slab9 +
      pmass_tb_e_slab9 +
      pmass_tb_f_slab9 +
      pmass_tb_g_slab9 +
      pmass_tb_h_slab9)
    disu9 <- -pmc * gamu9 * yintu9 * 0.2050409091

    # -- tb slab 10 (x = 1.84537 dm) --
    vuf10 <- 7.906312043e-05; vup10 <- 0.0004048111709; vus10 <- 0.05019504822
    cuf10 <- elf_tb_slab10 / vuf10
    cup10 <- epithelium_tb_slab10 / vup10
    cus10 <- subepithelium_tb_slab10 / vus10
    que10 <- peff * 1.443266165; qus10 <- peff * 1.500815121
    qu10 <- qbr * 0.1563165917
    gamu10 <- -kdiss * (cs - fuelf * cuf10)
    wu10a <- min(particles_tb_a_slab10, max(0, pmass_tb_a_slab10 / hza))
    wu10b <- min(particles_tb_b_slab10, max(0, pmass_tb_b_slab10 / hzb))
    wu10c <- min(particles_tb_c_slab10, max(0, pmass_tb_c_slab10 / hzc))
    wu10d <- min(particles_tb_d_slab10, max(0, pmass_tb_d_slab10 / hzd))
    wu10e <- min(particles_tb_e_slab10, max(0, pmass_tb_e_slab10 / hze))
    wu10f <- min(particles_tb_f_slab10, max(0, pmass_tb_f_slab10 / hzf))
    wu10g <- min(particles_tb_g_slab10, max(0, pmass_tb_g_slab10 / hzg))
    wu10h <- min(particles_tb_h_slab10, max(0, pmass_tb_h_slab10 / hzh))
    yintu10 <- wu10a + wu10b + wu10c + wu10d + wu10e + wu10f + wu10g + wu10h
    pmu10 <- pmc * (pmass_tb_a_slab10 +
      pmass_tb_b_slab10 +
      pmass_tb_c_slab10 +
      pmass_tb_d_slab10 +
      pmass_tb_e_slab10 +
      pmass_tb_f_slab10 +
      pmass_tb_g_slab10 +
      pmass_tb_h_slab10)
    disu10 <- -pmc * gamu10 * yintu10 * 0.2050409091

    # -- tb slab 11 (x = 2.05041 dm) --
    vuf11 <- 0.0002767656165; vup11 <- 0.0004580807835; vus11 <- 0.0205146674
    cuf11 <- elf_tb_slab11 / vuf11
    cup11 <- epithelium_tb_slab11 / vup11
    cus11 <- subepithelium_tb_slab11 / vus11
    que11 <- peff * 6.373688624; qus11 <- peff * 6.488726383
    qu11 <- qbr * 0.1101472592
    gamu11 <- -kdiss * (cs - fuelf * cuf11)
    wu11a <- min(particles_tb_a_slab11, max(0, pmass_tb_a_slab11 / hza))
    wu11b <- min(particles_tb_b_slab11, max(0, pmass_tb_b_slab11 / hzb))
    wu11c <- min(particles_tb_c_slab11, max(0, pmass_tb_c_slab11 / hzc))
    wu11d <- min(particles_tb_d_slab11, max(0, pmass_tb_d_slab11 / hzd))
    wu11e <- min(particles_tb_e_slab11, max(0, pmass_tb_e_slab11 / hze))
    wu11f <- min(particles_tb_f_slab11, max(0, pmass_tb_f_slab11 / hzf))
    wu11g <- min(particles_tb_g_slab11, max(0, pmass_tb_g_slab11 / hzg))
    wu11h <- min(particles_tb_h_slab11, max(0, pmass_tb_h_slab11 / hzh))
    yintu11 <- wu11a + wu11b + wu11c + wu11d + wu11e + wu11f + wu11g + wu11h
    pmu11 <- pmc * (pmass_tb_a_slab11 +
      pmass_tb_b_slab11 +
      pmass_tb_c_slab11 +
      pmass_tb_d_slab11 +
      pmass_tb_e_slab11 +
      pmass_tb_f_slab11 +
      pmass_tb_g_slab11 +
      pmass_tb_h_slab11)
    disu11 <- -pmc * gamu11 * yintu11 * 0.2050409091

    # -- tb slab 12 (x = 2.25545 dm) --
    vuf12 <- 0.003508095001; vup12 <- 0.009102683329; vus12 <- 0.08170383473
    cuf12 <- elf_tb_slab12 / vuf12
    cup12 <- epithelium_tb_slab12 / vup12
    cus12 <- subepithelium_tb_slab12 / vus12
    que12 <- peff * 118.2029603; qus12 <- peff * 124.5352618
    qu12 <- qbr * 0.6667357366
    gamu12 <- -kdiss * (cs - fuelf * cuf12)
    wu12a <- min(particles_tb_a_slab12, max(0, pmass_tb_a_slab12 / hza))
    wu12b <- min(particles_tb_b_slab12, max(0, pmass_tb_b_slab12 / hzb))
    wu12c <- min(particles_tb_c_slab12, max(0, pmass_tb_c_slab12 / hzc))
    wu12d <- min(particles_tb_d_slab12, max(0, pmass_tb_d_slab12 / hzd))
    wu12e <- min(particles_tb_e_slab12, max(0, pmass_tb_e_slab12 / hze))
    wu12f <- min(particles_tb_f_slab12, max(0, pmass_tb_f_slab12 / hzf))
    wu12g <- min(particles_tb_g_slab12, max(0, pmass_tb_g_slab12 / hzg))
    wu12h <- min(particles_tb_h_slab12, max(0, pmass_tb_h_slab12 / hzh))
    yintu12 <- wu12a + wu12b + wu12c + wu12d + wu12e + wu12f + wu12g + wu12h
    pmu12 <- pmc * (pmass_tb_a_slab12 +
      pmass_tb_b_slab12 +
      pmass_tb_c_slab12 +
      pmass_tb_d_slab12 +
      pmass_tb_e_slab12 +
      pmass_tb_f_slab12 +
      pmass_tb_g_slab12 +
      pmass_tb_h_slab12)
    disu12 <- -pmc * gamu12 * yintu12 * 0.1025204545

    # ================= alveolar region =================
    # -- alv slab 1 (x = 2.25545 dm) --
    vlf1 <- 2.118316689e-06; vlp1 <- 1.093287077e-05; vls1 <- 0.0007765911139
    clf1 <- elf_alv_slab1 / vlf1
    clp1 <- epithelium_alv_slab1 / vlp1
    cls1 <- subepithelium_alv_slab1 / vls1
    qle1 <- peff * 6.053090035; qls1 <- peff * 6.060894197
    ql1 <- qco * 0.003303802586
    gaml1 <- -kdiss * (cs - fuelf * clf1)
    wl1a <- min(particles_alv_a_slab1, max(0, pmass_alv_a_slab1 / hza))
    wl1b <- min(particles_alv_b_slab1, max(0, pmass_alv_b_slab1 / hzb))
    wl1c <- min(particles_alv_c_slab1, max(0, pmass_alv_c_slab1 / hzc))
    wl1d <- min(particles_alv_d_slab1, max(0, pmass_alv_d_slab1 / hzd))
    wl1e <- min(particles_alv_e_slab1, max(0, pmass_alv_e_slab1 / hze))
    wl1f <- min(particles_alv_f_slab1, max(0, pmass_alv_f_slab1 / hzf))
    wl1g <- min(particles_alv_g_slab1, max(0, pmass_alv_g_slab1 / hzg))
    wl1h <- min(particles_alv_h_slab1, max(0, pmass_alv_h_slab1 / hzh))
    yintl1 <- wl1a + wl1b + wl1c + wl1d + wl1e + wl1f + wl1g + wl1h
    pml1 <- pmc * (pmass_alv_a_slab1 +
      pmass_alv_b_slab1 +
      pmass_alv_c_slab1 +
      pmass_alv_d_slab1 +
      pmass_alv_e_slab1 +
      pmass_alv_f_slab1 +
      pmass_alv_g_slab1 +
      pmass_alv_h_slab1)
    disl1 <- -pmc * gaml1 * yintl1 *     0.00525

    # -- alv slab 2 (x = 2.26595 dm) --
    vlf2 <- 6.063049221e-06; vlp2 <- 3.12939954e-05; vls2 <- 0.002307952975
    clf2 <- elf_alv_slab2 / vlf2
    clp2 <- epithelium_alv_slab2 / vlp2
    cls2 <- subepithelium_alv_slab2 / vls2
    qle2 <- peff * 111.4578363; qls2 <- peff * 111.6217573
    ql2 <- qco * 0.009818578749
    gaml2 <- -kdiss * (cs - fuelf * clf2)
    wl2a <- min(particles_alv_a_slab2, max(0, pmass_alv_a_slab2 / hza))
    wl2b <- min(particles_alv_b_slab2, max(0, pmass_alv_b_slab2 / hzb))
    wl2c <- min(particles_alv_c_slab2, max(0, pmass_alv_c_slab2 / hzc))
    wl2d <- min(particles_alv_d_slab2, max(0, pmass_alv_d_slab2 / hzd))
    wl2e <- min(particles_alv_e_slab2, max(0, pmass_alv_e_slab2 / hze))
    wl2f <- min(particles_alv_f_slab2, max(0, pmass_alv_f_slab2 / hzf))
    wl2g <- min(particles_alv_g_slab2, max(0, pmass_alv_g_slab2 / hzg))
    wl2h <- min(particles_alv_h_slab2, max(0, pmass_alv_h_slab2 / hzh))
    yintl2 <- wl2a + wl2b + wl2c + wl2d + wl2e + wl2f + wl2g + wl2h
    pml2 <- pmc * (pmass_alv_a_slab2 +
      pmass_alv_b_slab2 +
      pmass_alv_c_slab2 +
      pmass_alv_d_slab2 +
      pmass_alv_e_slab2 +
      pmass_alv_f_slab2 +
      pmass_alv_g_slab2 +
      pmass_alv_h_slab2)
    disl2 <- -pmc * gaml2 * yintl2 *      0.0105

    # -- alv slab 3 (x = 2.27645 dm) --
    vlf3 <- 1.348342041e-05; vlp3 <- 6.958045711e-05; vls3 <- 0.003561309935
    clf3 <- elf_alv_slab3 / vlf3
    clp3 <- epithelium_alv_slab3 / vlp3
    cls3 <- subepithelium_alv_slab3 / vls3
    qle3 <- peff * 631.4909041; qls3 <- peff * 631.9330213
    ql3 <- qco * 0.0151506562
    gaml3 <- -kdiss * (cs - fuelf * clf3)
    wl3a <- min(particles_alv_a_slab3, max(0, pmass_alv_a_slab3 / hza))
    wl3b <- min(particles_alv_b_slab3, max(0, pmass_alv_b_slab3 / hzb))
    wl3c <- min(particles_alv_c_slab3, max(0, pmass_alv_c_slab3 / hzc))
    wl3d <- min(particles_alv_d_slab3, max(0, pmass_alv_d_slab3 / hzd))
    wl3e <- min(particles_alv_e_slab3, max(0, pmass_alv_e_slab3 / hze))
    wl3f <- min(particles_alv_f_slab3, max(0, pmass_alv_f_slab3 / hzf))
    wl3g <- min(particles_alv_g_slab3, max(0, pmass_alv_g_slab3 / hzg))
    wl3h <- min(particles_alv_h_slab3, max(0, pmass_alv_h_slab3 / hzh))
    yintl3 <- wl3a + wl3b + wl3c + wl3d + wl3e + wl3f + wl3g + wl3h
    pml3 <- pmc * (pmass_alv_a_slab3 +
      pmass_alv_b_slab3 +
      pmass_alv_c_slab3 +
      pmass_alv_d_slab3 +
      pmass_alv_e_slab3 +
      pmass_alv_f_slab3 +
      pmass_alv_g_slab3 +
      pmass_alv_h_slab3)
    disl3 <- -pmc * gaml3 * yintl3 *      0.0105

    # -- alv slab 4 (x = 2.28695 dm) --
    vlf4 <- 4.148341334e-05; vlp4 <- 0.0002140150633; vls4 <- 0.00590746216
    clf4 <- elf_alv_slab4 / vlf4
    clp4 <- epithelium_alv_slab4 / vlp4
    cls4 <- subepithelium_alv_slab4 / vls4
    qle4 <- peff *  2738.70014; qls4 <- peff * 2739.435159
    ql4 <- qco * 0.02513174361
    gaml4 <- -kdiss * (cs - fuelf * clf4)
    wl4a <- min(particles_alv_a_slab4, max(0, pmass_alv_a_slab4 / hza))
    wl4b <- min(particles_alv_b_slab4, max(0, pmass_alv_b_slab4 / hzb))
    wl4c <- min(particles_alv_c_slab4, max(0, pmass_alv_c_slab4 / hzc))
    wl4d <- min(particles_alv_d_slab4, max(0, pmass_alv_d_slab4 / hzd))
    wl4e <- min(particles_alv_e_slab4, max(0, pmass_alv_e_slab4 / hze))
    wl4f <- min(particles_alv_f_slab4, max(0, pmass_alv_f_slab4 / hzf))
    wl4g <- min(particles_alv_g_slab4, max(0, pmass_alv_g_slab4 / hzg))
    wl4h <- min(particles_alv_h_slab4, max(0, pmass_alv_h_slab4 / hzh))
    yintl4 <- wl4a + wl4b + wl4c + wl4d + wl4e + wl4f + wl4g + wl4h
    pml4 <- pmc * (pmass_alv_a_slab4 +
      pmass_alv_b_slab4 +
      pmass_alv_c_slab4 +
      pmass_alv_d_slab4 +
      pmass_alv_e_slab4 +
      pmass_alv_f_slab4 +
      pmass_alv_g_slab4 +
      pmass_alv_h_slab4)
    disl4 <- -pmc * gaml4 * yintl4 *      0.0105

    # -- alv slab 5 (x = 2.29745 dm) --
    vlf5 <- 0.0001263306547; vlp5 <- 0.0006516619444; vls5 <- 0.01098217376
    clf5 <- elf_alv_slab5 / vlf5
    clp5 <- epithelium_alv_slab5 / vlp5
    cls5 <- subepithelium_alv_slab5 / vls5
    qle5 <- peff * 9513.371249; qls5 <- peff * 9514.825618
    ql5 <- qco * 0.04672076904
    gaml5 <- -kdiss * (cs - fuelf * clf5)
    wl5a <- min(particles_alv_a_slab5, max(0, pmass_alv_a_slab5 / hza))
    wl5b <- min(particles_alv_b_slab5, max(0, pmass_alv_b_slab5 / hzb))
    wl5c <- min(particles_alv_c_slab5, max(0, pmass_alv_c_slab5 / hzc))
    wl5d <- min(particles_alv_d_slab5, max(0, pmass_alv_d_slab5 / hzd))
    wl5e <- min(particles_alv_e_slab5, max(0, pmass_alv_e_slab5 / hze))
    wl5f <- min(particles_alv_f_slab5, max(0, pmass_alv_f_slab5 / hzf))
    wl5g <- min(particles_alv_g_slab5, max(0, pmass_alv_g_slab5 / hzg))
    wl5h <- min(particles_alv_h_slab5, max(0, pmass_alv_h_slab5 / hzh))
    yintl5 <- wl5a + wl5b + wl5c + wl5d + wl5e + wl5f + wl5g + wl5h
    pml5 <- pmc * (pmass_alv_a_slab5 +
      pmass_alv_b_slab5 +
      pmass_alv_c_slab5 +
      pmass_alv_d_slab5 +
      pmass_alv_e_slab5 +
      pmass_alv_f_slab5 +
      pmass_alv_g_slab5 +
      pmass_alv_h_slab5)
    disl5 <- -pmc * gaml5 * yintl5 *      0.0105

    # -- alv slab 6 (x = 2.30795 dm) --
    vlf6 <- 0.0005593265354; vlp6 <- 0.002884891048; vls6 <- 0.02337101051
    clf6 <- elf_alv_slab6 / vlf6
    clp6 <- epithelium_alv_slab6 / vlp6
    cls6 <- subepithelium_alv_slab6 / vls6
    qle6 <- peff * 46293.57826; qls6 <- peff * 46296.95989
    ql6 <- qco * 0.09942581572
    gaml6 <- -kdiss * (cs - fuelf * clf6)
    wl6a <- min(particles_alv_a_slab6, max(0, pmass_alv_a_slab6 / hza))
    wl6b <- min(particles_alv_b_slab6, max(0, pmass_alv_b_slab6 / hzb))
    wl6c <- min(particles_alv_c_slab6, max(0, pmass_alv_c_slab6 / hzc))
    wl6d <- min(particles_alv_d_slab6, max(0, pmass_alv_d_slab6 / hzd))
    wl6e <- min(particles_alv_e_slab6, max(0, pmass_alv_e_slab6 / hze))
    wl6f <- min(particles_alv_f_slab6, max(0, pmass_alv_f_slab6 / hzf))
    wl6g <- min(particles_alv_g_slab6, max(0, pmass_alv_g_slab6 / hzg))
    wl6h <- min(particles_alv_h_slab6, max(0, pmass_alv_h_slab6 / hzh))
    yintl6 <- wl6a + wl6b + wl6c + wl6d + wl6e + wl6f + wl6g + wl6h
    pml6 <- pmc * (pmass_alv_a_slab6 +
      pmass_alv_b_slab6 +
      pmass_alv_c_slab6 +
      pmass_alv_d_slab6 +
      pmass_alv_e_slab6 +
      pmass_alv_f_slab6 +
      pmass_alv_g_slab6 +
      pmass_alv_h_slab6)
    disl6 <- -pmc * gaml6 * yintl6 *      0.0105

    # -- alv slab 7 (x = 2.31845 dm) --
    vlf7 <- 0.002495939728; vlp7 <- 0.01287295348; vls7 <- 0.06351171644
    clf7 <- elf_alv_slab7 / vlf7
    clp7 <- epithelium_alv_slab7 / vlp7
    cls7 <- subepithelium_alv_slab7 / vls7
    qle7 <- peff * 213470.4355; qls7 <- peff * 213480.0477
    ql7 <- qco * 0.2701938888
    gaml7 <- -kdiss * (cs - fuelf * clf7)
    wl7a <- min(particles_alv_a_slab7, max(0, pmass_alv_a_slab7 / hza))
    wl7b <- min(particles_alv_b_slab7, max(0, pmass_alv_b_slab7 / hzb))
    wl7c <- min(particles_alv_c_slab7, max(0, pmass_alv_c_slab7 / hzc))
    wl7d <- min(particles_alv_d_slab7, max(0, pmass_alv_d_slab7 / hzd))
    wl7e <- min(particles_alv_e_slab7, max(0, pmass_alv_e_slab7 / hze))
    wl7f <- min(particles_alv_f_slab7, max(0, pmass_alv_f_slab7 / hzf))
    wl7g <- min(particles_alv_g_slab7, max(0, pmass_alv_g_slab7 / hzg))
    wl7h <- min(particles_alv_h_slab7, max(0, pmass_alv_h_slab7 / hzh))
    yintl7 <- wl7a + wl7b + wl7c + wl7d + wl7e + wl7f + wl7g + wl7h
    pml7 <- pmc * (pmass_alv_a_slab7 +
      pmass_alv_b_slab7 +
      pmass_alv_c_slab7 +
      pmass_alv_d_slab7 +
      pmass_alv_e_slab7 +
      pmass_alv_f_slab7 +
      pmass_alv_g_slab7 +
      pmass_alv_h_slab7)
    disl7 <- -pmc * gaml7 * yintl7 *      0.0105

    # -- alv slab 8 (x = 2.32895 dm) --
    vlf8 <- 0.006203039848; vlp8 <- 0.03199211534; vls8 <- 0.1246415646
    clf8 <- elf_alv_slab8 / vlf8
    clp8 <- epithelium_alv_slab8 / vlp8
    cls8 <- subepithelium_alv_slab8 / vls8
    qle8 <- peff * 536048.4349; qls8 <- peff * 536068.4136
    ql8 <- qco * 0.5302547453
    gaml8 <- -kdiss * (cs - fuelf * clf8)
    wl8a <- min(particles_alv_a_slab8, max(0, pmass_alv_a_slab8 / hza))
    wl8b <- min(particles_alv_b_slab8, max(0, pmass_alv_b_slab8 / hzb))
    wl8c <- min(particles_alv_c_slab8, max(0, pmass_alv_c_slab8 / hzc))
    wl8d <- min(particles_alv_d_slab8, max(0, pmass_alv_d_slab8 / hzd))
    wl8e <- min(particles_alv_e_slab8, max(0, pmass_alv_e_slab8 / hze))
    wl8f <- min(particles_alv_f_slab8, max(0, pmass_alv_f_slab8 / hzf))
    wl8g <- min(particles_alv_g_slab8, max(0, pmass_alv_g_slab8 / hzg))
    wl8h <- min(particles_alv_h_slab8, max(0, pmass_alv_h_slab8 / hzh))
    yintl8 <- wl8a + wl8b + wl8c + wl8d + wl8e + wl8f + wl8g + wl8h
    pml8 <- pmc * (pmass_alv_a_slab8 +
      pmass_alv_b_slab8 +
      pmass_alv_c_slab8 +
      pmass_alv_d_slab8 +
      pmass_alv_e_slab8 +
      pmass_alv_f_slab8 +
      pmass_alv_g_slab8 +
      pmass_alv_h_slab8)
    disl8 <- -pmc * gaml8 * yintl8 *     0.00525

    # ---- mucociliary clearance velocity (Eq. S7; alpha <= 0) ----
    alp1 <- mccvel * (       -1)
    alp2 <- mccvel * (-0.91478489)
    alp3 <- mccvel * (-0.81834906)
    alp4 <- mccvel * (-0.71221726)
    alp5 <- mccvel * (-0.59791426)
    alp6 <- mccvel * (-0.47696481)
    alp7 <- mccvel * (-0.35089368)
    alp8 <- mccvel * (-0.21789601)
    alp9 <- mccvel * (-0.068290643)
    alp10 <- mccvel * (-0.023651828)
    alp11 <- mccvel * (-0.0060322285)
    alp12 <- mccvel * (       -0)
    # Eq. 5: drug leaving the lung to the gut at the proximal boundary.
    mgut <- -alp1 * pmu1

    # ================= ODEs =================
    # Systemic PBPK, Eqs S8-S16 (R = bp).
    d/dt(a_spleen) <- qspleen * (cart - bp * csp / kpspleen)
    d/dt(a_rapidly_perfused) <- qrapid * (cart - bp * cri / kprapid)
    d/dt(a_slowly_perfused) <- qslow * (cart - bp * cpo / kpslow)
    d/dt(a_fat) <- qfat * (cart - bp * cad / kpfat)
    d/dt(a_liver) <- qliver * cart + qgut * bp * cgut / kpgut +
      qspleen * bp * csp / kpspleen -
      (qliver + qspleen + qgut) * bp * chep / kpliver
    d/dt(a_gut) <- qgut * (cart - bp * cgut / kpgut) + foral * ka * depot
    d/dt(depot) <- -ka * depot + mgut
    d/dt(a_arterial) <- ql1 * bp * cls1 / kplung +
      ql2 * bp * cls2 / kplung +
      ql3 * bp * cls3 / kplung +
      ql4 * bp * cls4 / kplung +
      ql5 * bp * cls5 / kplung +
      ql6 * bp * cls6 / kplung +
      ql7 * bp * cls7 / kplung +
      ql8 * bp * cls8 / kplung -
      (qspleen + qrapid + qslow + qfat + qliver + qgut + qbr) * cart
    d/dt(a_venous) <- qu1 * bp * cus1 / kplung +
      qu2 * bp * cus2 / kplung +
      qu3 * bp * cus3 / kplung +
      qu4 * bp * cus4 / kplung +
      qu5 * bp * cus5 / kplung +
      qu6 * bp * cus6 / kplung +
      qu7 * bp * cus7 / kplung +
      qu8 * bp * cus8 / kplung +
      qu9 * bp * cus9 / kplung +
      qu10 * bp * cus10 / kplung +
      qu11 * bp * cus11 / kplung +
      qu12 * bp * cus12 / kplung -
      qco * cven +
      qrapid * bp * cri / kprapid +
      qslow * bp * cpo / kpslow +
      qfat * bp * cad / kpfat +
      (qliver + qspleen + qgut) * bp * chep / kpliver -
      cl * cven

    # Tracheobronchial slabs (Eqs 7-9); perfused by bronchial blood flow.
    d/dt(elf_tb_slab1) <- -fuelf * que1 * cuf1 + que1 * cup1 / kpulung + disu1
    d/dt(epithelium_tb_slab1) <- fuelf * que1 * cuf1 - que1 * cup1 / kpulung -
      qus1 * cup1 / kpulung + qus1 * cus1 / kpulung
    d/dt(subepithelium_tb_slab1) <- qus1 * cup1 / kpulung - qus1 * cus1 / kpulung +
      qu1 * cart - qu1 * bp * cus1 / kplung
    d/dt(elf_tb_slab2) <- -fuelf * que2 * cuf2 + que2 * cup2 / kpulung + disu2
    d/dt(epithelium_tb_slab2) <- fuelf * que2 * cuf2 - que2 * cup2 / kpulung -
      qus2 * cup2 / kpulung + qus2 * cus2 / kpulung
    d/dt(subepithelium_tb_slab2) <- qus2 * cup2 / kpulung - qus2 * cus2 / kpulung +
      qu2 * cart - qu2 * bp * cus2 / kplung
    d/dt(elf_tb_slab3) <- -fuelf * que3 * cuf3 + que3 * cup3 / kpulung + disu3
    d/dt(epithelium_tb_slab3) <- fuelf * que3 * cuf3 - que3 * cup3 / kpulung -
      qus3 * cup3 / kpulung + qus3 * cus3 / kpulung
    d/dt(subepithelium_tb_slab3) <- qus3 * cup3 / kpulung - qus3 * cus3 / kpulung +
      qu3 * cart - qu3 * bp * cus3 / kplung
    d/dt(elf_tb_slab4) <- -fuelf * que4 * cuf4 + que4 * cup4 / kpulung + disu4
    d/dt(epithelium_tb_slab4) <- fuelf * que4 * cuf4 - que4 * cup4 / kpulung -
      qus4 * cup4 / kpulung + qus4 * cus4 / kpulung
    d/dt(subepithelium_tb_slab4) <- qus4 * cup4 / kpulung - qus4 * cus4 / kpulung +
      qu4 * cart - qu4 * bp * cus4 / kplung
    d/dt(elf_tb_slab5) <- -fuelf * que5 * cuf5 + que5 * cup5 / kpulung + disu5
    d/dt(epithelium_tb_slab5) <- fuelf * que5 * cuf5 - que5 * cup5 / kpulung -
      qus5 * cup5 / kpulung + qus5 * cus5 / kpulung
    d/dt(subepithelium_tb_slab5) <- qus5 * cup5 / kpulung - qus5 * cus5 / kpulung +
      qu5 * cart - qu5 * bp * cus5 / kplung
    d/dt(elf_tb_slab6) <- -fuelf * que6 * cuf6 + que6 * cup6 / kpulung + disu6
    d/dt(epithelium_tb_slab6) <- fuelf * que6 * cuf6 - que6 * cup6 / kpulung -
      qus6 * cup6 / kpulung + qus6 * cus6 / kpulung
    d/dt(subepithelium_tb_slab6) <- qus6 * cup6 / kpulung - qus6 * cus6 / kpulung +
      qu6 * cart - qu6 * bp * cus6 / kplung
    d/dt(elf_tb_slab7) <- -fuelf * que7 * cuf7 + que7 * cup7 / kpulung + disu7
    d/dt(epithelium_tb_slab7) <- fuelf * que7 * cuf7 - que7 * cup7 / kpulung -
      qus7 * cup7 / kpulung + qus7 * cus7 / kpulung
    d/dt(subepithelium_tb_slab7) <- qus7 * cup7 / kpulung - qus7 * cus7 / kpulung +
      qu7 * cart - qu7 * bp * cus7 / kplung
    d/dt(elf_tb_slab8) <- -fuelf * que8 * cuf8 + que8 * cup8 / kpulung + disu8
    d/dt(epithelium_tb_slab8) <- fuelf * que8 * cuf8 - que8 * cup8 / kpulung -
      qus8 * cup8 / kpulung + qus8 * cus8 / kpulung
    d/dt(subepithelium_tb_slab8) <- qus8 * cup8 / kpulung - qus8 * cus8 / kpulung +
      qu8 * cart - qu8 * bp * cus8 / kplung
    d/dt(elf_tb_slab9) <- -fuelf * que9 * cuf9 + que9 * cup9 / kpulung + disu9
    d/dt(epithelium_tb_slab9) <- fuelf * que9 * cuf9 - que9 * cup9 / kpulung -
      qus9 * cup9 / kpulung + qus9 * cus9 / kpulung
    d/dt(subepithelium_tb_slab9) <- qus9 * cup9 / kpulung - qus9 * cus9 / kpulung +
      qu9 * cart - qu9 * bp * cus9 / kplung
    d/dt(elf_tb_slab10) <- -fuelf * que10 * cuf10 + que10 * cup10 / kpulung + disu10
    d/dt(epithelium_tb_slab10) <- fuelf * que10 * cuf10 - que10 * cup10 / kpulung -
      qus10 * cup10 / kpulung + qus10 * cus10 / kpulung
    d/dt(subepithelium_tb_slab10) <- qus10 * cup10 / kpulung - qus10 * cus10 / kpulung +
      qu10 * cart - qu10 * bp * cus10 / kplung
    d/dt(elf_tb_slab11) <- -fuelf * que11 * cuf11 + que11 * cup11 / kpulung + disu11
    d/dt(epithelium_tb_slab11) <- fuelf * que11 * cuf11 - que11 * cup11 / kpulung -
      qus11 * cup11 / kpulung + qus11 * cus11 / kpulung
    d/dt(subepithelium_tb_slab11) <- qus11 * cup11 / kpulung - qus11 * cus11 / kpulung +
      qu11 * cart - qu11 * bp * cus11 / kplung
    d/dt(elf_tb_slab12) <- -fuelf * que12 * cuf12 + que12 * cup12 / kpulung + disu12
    d/dt(epithelium_tb_slab12) <- fuelf * que12 * cuf12 - que12 * cup12 / kpulung -
      qus12 * cup12 / kpulung + qus12 * cus12 / kpulung
    d/dt(subepithelium_tb_slab12) <- qus12 * cup12 / kpulung - qus12 * cus12 / kpulung +
      qu12 * cart - qu12 * bp * cus12 / kplung

    # Alveolar slabs (Eqs 7-9); perfused by the whole cardiac output.
    d/dt(elf_alv_slab1) <- -fuelf * qle1 * clf1 + qle1 * clp1 / kpulung + disl1
    d/dt(epithelium_alv_slab1) <- fuelf * qle1 * clf1 - qle1 * clp1 / kpulung -
      qls1 * clp1 / kpulung + qls1 * cls1 / kpulung
    d/dt(subepithelium_alv_slab1) <- qls1 * clp1 / kpulung - qls1 * cls1 / kpulung +
      ql1 * cven - ql1 * bp * cls1 / kplung
    d/dt(elf_alv_slab2) <- -fuelf * qle2 * clf2 + qle2 * clp2 / kpulung + disl2
    d/dt(epithelium_alv_slab2) <- fuelf * qle2 * clf2 - qle2 * clp2 / kpulung -
      qls2 * clp2 / kpulung + qls2 * cls2 / kpulung
    d/dt(subepithelium_alv_slab2) <- qls2 * clp2 / kpulung - qls2 * cls2 / kpulung +
      ql2 * cven - ql2 * bp * cls2 / kplung
    d/dt(elf_alv_slab3) <- -fuelf * qle3 * clf3 + qle3 * clp3 / kpulung + disl3
    d/dt(epithelium_alv_slab3) <- fuelf * qle3 * clf3 - qle3 * clp3 / kpulung -
      qls3 * clp3 / kpulung + qls3 * cls3 / kpulung
    d/dt(subepithelium_alv_slab3) <- qls3 * clp3 / kpulung - qls3 * cls3 / kpulung +
      ql3 * cven - ql3 * bp * cls3 / kplung
    d/dt(elf_alv_slab4) <- -fuelf * qle4 * clf4 + qle4 * clp4 / kpulung + disl4
    d/dt(epithelium_alv_slab4) <- fuelf * qle4 * clf4 - qle4 * clp4 / kpulung -
      qls4 * clp4 / kpulung + qls4 * cls4 / kpulung
    d/dt(subepithelium_alv_slab4) <- qls4 * clp4 / kpulung - qls4 * cls4 / kpulung +
      ql4 * cven - ql4 * bp * cls4 / kplung
    d/dt(elf_alv_slab5) <- -fuelf * qle5 * clf5 + qle5 * clp5 / kpulung + disl5
    d/dt(epithelium_alv_slab5) <- fuelf * qle5 * clf5 - qle5 * clp5 / kpulung -
      qls5 * clp5 / kpulung + qls5 * cls5 / kpulung
    d/dt(subepithelium_alv_slab5) <- qls5 * clp5 / kpulung - qls5 * cls5 / kpulung +
      ql5 * cven - ql5 * bp * cls5 / kplung
    d/dt(elf_alv_slab6) <- -fuelf * qle6 * clf6 + qle6 * clp6 / kpulung + disl6
    d/dt(epithelium_alv_slab6) <- fuelf * qle6 * clf6 - qle6 * clp6 / kpulung -
      qls6 * clp6 / kpulung + qls6 * cls6 / kpulung
    d/dt(subepithelium_alv_slab6) <- qls6 * clp6 / kpulung - qls6 * cls6 / kpulung +
      ql6 * cven - ql6 * bp * cls6 / kplung
    d/dt(elf_alv_slab7) <- -fuelf * qle7 * clf7 + qle7 * clp7 / kpulung + disl7
    d/dt(epithelium_alv_slab7) <- fuelf * qle7 * clf7 - qle7 * clp7 / kpulung -
      qls7 * clp7 / kpulung + qls7 * cls7 / kpulung
    d/dt(subepithelium_alv_slab7) <- qls7 * clp7 / kpulung - qls7 * cls7 / kpulung +
      ql7 * cven - ql7 * bp * cls7 / kplung
    d/dt(elf_alv_slab8) <- -fuelf * qle8 * clf8 + qle8 * clp8 / kpulung + disl8
    d/dt(epithelium_alv_slab8) <- fuelf * qle8 * clf8 - qle8 * clp8 / kpulung -
      qls8 * clp8 / kpulung + qls8 * cls8 / kpulung
    d/dt(subepithelium_alv_slab8) <- qls8 * clp8 / kpulung - qls8 * cls8 / kpulung +
      ql8 * cven - ql8 * bp * cls8 / kplung

    # Particle number density, advected by MCC in conservative upwind
    # form (Eq. 1). alpha <= 0, so the upwind face value is the distal
    # neighbour; dividing by this slab's own cell width makes the chain
    # telescope to the single boundary flux alp1 * pmu1.
    d/dt(particles_tb_a_slab1) <- -(alp2 * particles_tb_a_slab2 - alp1 * particles_tb_a_slab1) / 0.1025204545
    d/dt(particles_tb_b_slab1) <- -(alp2 * particles_tb_b_slab2 - alp1 * particles_tb_b_slab1) / 0.1025204545
    d/dt(particles_tb_c_slab1) <- -(alp2 * particles_tb_c_slab2 - alp1 * particles_tb_c_slab1) / 0.1025204545
    d/dt(particles_tb_d_slab1) <- -(alp2 * particles_tb_d_slab2 - alp1 * particles_tb_d_slab1) / 0.1025204545
    d/dt(particles_tb_e_slab1) <- -(alp2 * particles_tb_e_slab2 - alp1 * particles_tb_e_slab1) / 0.1025204545
    d/dt(particles_tb_f_slab1) <- -(alp2 * particles_tb_f_slab2 - alp1 * particles_tb_f_slab1) / 0.1025204545
    d/dt(particles_tb_g_slab1) <- -(alp2 * particles_tb_g_slab2 - alp1 * particles_tb_g_slab1) / 0.1025204545
    d/dt(particles_tb_h_slab1) <- -(alp2 * particles_tb_h_slab2 - alp1 * particles_tb_h_slab1) / 0.1025204545
    d/dt(particles_tb_a_slab2) <- -(alp3 * particles_tb_a_slab3 - alp2 * particles_tb_a_slab2) / 0.2050409091
    d/dt(particles_tb_b_slab2) <- -(alp3 * particles_tb_b_slab3 - alp2 * particles_tb_b_slab2) / 0.2050409091
    d/dt(particles_tb_c_slab2) <- -(alp3 * particles_tb_c_slab3 - alp2 * particles_tb_c_slab2) / 0.2050409091
    d/dt(particles_tb_d_slab2) <- -(alp3 * particles_tb_d_slab3 - alp2 * particles_tb_d_slab2) / 0.2050409091
    d/dt(particles_tb_e_slab2) <- -(alp3 * particles_tb_e_slab3 - alp2 * particles_tb_e_slab2) / 0.2050409091
    d/dt(particles_tb_f_slab2) <- -(alp3 * particles_tb_f_slab3 - alp2 * particles_tb_f_slab2) / 0.2050409091
    d/dt(particles_tb_g_slab2) <- -(alp3 * particles_tb_g_slab3 - alp2 * particles_tb_g_slab2) / 0.2050409091
    d/dt(particles_tb_h_slab2) <- -(alp3 * particles_tb_h_slab3 - alp2 * particles_tb_h_slab2) / 0.2050409091
    d/dt(particles_tb_a_slab3) <- -(alp4 * particles_tb_a_slab4 - alp3 * particles_tb_a_slab3) / 0.2050409091
    d/dt(particles_tb_b_slab3) <- -(alp4 * particles_tb_b_slab4 - alp3 * particles_tb_b_slab3) / 0.2050409091
    d/dt(particles_tb_c_slab3) <- -(alp4 * particles_tb_c_slab4 - alp3 * particles_tb_c_slab3) / 0.2050409091
    d/dt(particles_tb_d_slab3) <- -(alp4 * particles_tb_d_slab4 - alp3 * particles_tb_d_slab3) / 0.2050409091
    d/dt(particles_tb_e_slab3) <- -(alp4 * particles_tb_e_slab4 - alp3 * particles_tb_e_slab3) / 0.2050409091
    d/dt(particles_tb_f_slab3) <- -(alp4 * particles_tb_f_slab4 - alp3 * particles_tb_f_slab3) / 0.2050409091
    d/dt(particles_tb_g_slab3) <- -(alp4 * particles_tb_g_slab4 - alp3 * particles_tb_g_slab3) / 0.2050409091
    d/dt(particles_tb_h_slab3) <- -(alp4 * particles_tb_h_slab4 - alp3 * particles_tb_h_slab3) / 0.2050409091
    d/dt(particles_tb_a_slab4) <- -(alp5 * particles_tb_a_slab5 - alp4 * particles_tb_a_slab4) / 0.2050409091
    d/dt(particles_tb_b_slab4) <- -(alp5 * particles_tb_b_slab5 - alp4 * particles_tb_b_slab4) / 0.2050409091
    d/dt(particles_tb_c_slab4) <- -(alp5 * particles_tb_c_slab5 - alp4 * particles_tb_c_slab4) / 0.2050409091
    d/dt(particles_tb_d_slab4) <- -(alp5 * particles_tb_d_slab5 - alp4 * particles_tb_d_slab4) / 0.2050409091
    d/dt(particles_tb_e_slab4) <- -(alp5 * particles_tb_e_slab5 - alp4 * particles_tb_e_slab4) / 0.2050409091
    d/dt(particles_tb_f_slab4) <- -(alp5 * particles_tb_f_slab5 - alp4 * particles_tb_f_slab4) / 0.2050409091
    d/dt(particles_tb_g_slab4) <- -(alp5 * particles_tb_g_slab5 - alp4 * particles_tb_g_slab4) / 0.2050409091
    d/dt(particles_tb_h_slab4) <- -(alp5 * particles_tb_h_slab5 - alp4 * particles_tb_h_slab4) / 0.2050409091
    d/dt(particles_tb_a_slab5) <- -(alp6 * particles_tb_a_slab6 - alp5 * particles_tb_a_slab5) / 0.2050409091
    d/dt(particles_tb_b_slab5) <- -(alp6 * particles_tb_b_slab6 - alp5 * particles_tb_b_slab5) / 0.2050409091
    d/dt(particles_tb_c_slab5) <- -(alp6 * particles_tb_c_slab6 - alp5 * particles_tb_c_slab5) / 0.2050409091
    d/dt(particles_tb_d_slab5) <- -(alp6 * particles_tb_d_slab6 - alp5 * particles_tb_d_slab5) / 0.2050409091
    d/dt(particles_tb_e_slab5) <- -(alp6 * particles_tb_e_slab6 - alp5 * particles_tb_e_slab5) / 0.2050409091
    d/dt(particles_tb_f_slab5) <- -(alp6 * particles_tb_f_slab6 - alp5 * particles_tb_f_slab5) / 0.2050409091
    d/dt(particles_tb_g_slab5) <- -(alp6 * particles_tb_g_slab6 - alp5 * particles_tb_g_slab5) / 0.2050409091
    d/dt(particles_tb_h_slab5) <- -(alp6 * particles_tb_h_slab6 - alp5 * particles_tb_h_slab5) / 0.2050409091
    d/dt(particles_tb_a_slab6) <- -(alp7 * particles_tb_a_slab7 - alp6 * particles_tb_a_slab6) / 0.2050409091
    d/dt(particles_tb_b_slab6) <- -(alp7 * particles_tb_b_slab7 - alp6 * particles_tb_b_slab6) / 0.2050409091
    d/dt(particles_tb_c_slab6) <- -(alp7 * particles_tb_c_slab7 - alp6 * particles_tb_c_slab6) / 0.2050409091
    d/dt(particles_tb_d_slab6) <- -(alp7 * particles_tb_d_slab7 - alp6 * particles_tb_d_slab6) / 0.2050409091
    d/dt(particles_tb_e_slab6) <- -(alp7 * particles_tb_e_slab7 - alp6 * particles_tb_e_slab6) / 0.2050409091
    d/dt(particles_tb_f_slab6) <- -(alp7 * particles_tb_f_slab7 - alp6 * particles_tb_f_slab6) / 0.2050409091
    d/dt(particles_tb_g_slab6) <- -(alp7 * particles_tb_g_slab7 - alp6 * particles_tb_g_slab6) / 0.2050409091
    d/dt(particles_tb_h_slab6) <- -(alp7 * particles_tb_h_slab7 - alp6 * particles_tb_h_slab6) / 0.2050409091
    d/dt(particles_tb_a_slab7) <- -(alp8 * particles_tb_a_slab8 - alp7 * particles_tb_a_slab7) / 0.2050409091
    d/dt(particles_tb_b_slab7) <- -(alp8 * particles_tb_b_slab8 - alp7 * particles_tb_b_slab7) / 0.2050409091
    d/dt(particles_tb_c_slab7) <- -(alp8 * particles_tb_c_slab8 - alp7 * particles_tb_c_slab7) / 0.2050409091
    d/dt(particles_tb_d_slab7) <- -(alp8 * particles_tb_d_slab8 - alp7 * particles_tb_d_slab7) / 0.2050409091
    d/dt(particles_tb_e_slab7) <- -(alp8 * particles_tb_e_slab8 - alp7 * particles_tb_e_slab7) / 0.2050409091
    d/dt(particles_tb_f_slab7) <- -(alp8 * particles_tb_f_slab8 - alp7 * particles_tb_f_slab7) / 0.2050409091
    d/dt(particles_tb_g_slab7) <- -(alp8 * particles_tb_g_slab8 - alp7 * particles_tb_g_slab7) / 0.2050409091
    d/dt(particles_tb_h_slab7) <- -(alp8 * particles_tb_h_slab8 - alp7 * particles_tb_h_slab7) / 0.2050409091
    d/dt(particles_tb_a_slab8) <- -(alp9 * particles_tb_a_slab9 - alp8 * particles_tb_a_slab8) / 0.2050409091
    d/dt(particles_tb_b_slab8) <- -(alp9 * particles_tb_b_slab9 - alp8 * particles_tb_b_slab8) / 0.2050409091
    d/dt(particles_tb_c_slab8) <- -(alp9 * particles_tb_c_slab9 - alp8 * particles_tb_c_slab8) / 0.2050409091
    d/dt(particles_tb_d_slab8) <- -(alp9 * particles_tb_d_slab9 - alp8 * particles_tb_d_slab8) / 0.2050409091
    d/dt(particles_tb_e_slab8) <- -(alp9 * particles_tb_e_slab9 - alp8 * particles_tb_e_slab8) / 0.2050409091
    d/dt(particles_tb_f_slab8) <- -(alp9 * particles_tb_f_slab9 - alp8 * particles_tb_f_slab8) / 0.2050409091
    d/dt(particles_tb_g_slab8) <- -(alp9 * particles_tb_g_slab9 - alp8 * particles_tb_g_slab8) / 0.2050409091
    d/dt(particles_tb_h_slab8) <- -(alp9 * particles_tb_h_slab9 - alp8 * particles_tb_h_slab8) / 0.2050409091
    d/dt(particles_tb_a_slab9) <- -(alp10 * particles_tb_a_slab10 - alp9 * particles_tb_a_slab9) / 0.2050409091
    d/dt(particles_tb_b_slab9) <- -(alp10 * particles_tb_b_slab10 - alp9 * particles_tb_b_slab9) / 0.2050409091
    d/dt(particles_tb_c_slab9) <- -(alp10 * particles_tb_c_slab10 - alp9 * particles_tb_c_slab9) / 0.2050409091
    d/dt(particles_tb_d_slab9) <- -(alp10 * particles_tb_d_slab10 - alp9 * particles_tb_d_slab9) / 0.2050409091
    d/dt(particles_tb_e_slab9) <- -(alp10 * particles_tb_e_slab10 - alp9 * particles_tb_e_slab9) / 0.2050409091
    d/dt(particles_tb_f_slab9) <- -(alp10 * particles_tb_f_slab10 - alp9 * particles_tb_f_slab9) / 0.2050409091
    d/dt(particles_tb_g_slab9) <- -(alp10 * particles_tb_g_slab10 - alp9 * particles_tb_g_slab9) / 0.2050409091
    d/dt(particles_tb_h_slab9) <- -(alp10 * particles_tb_h_slab10 - alp9 * particles_tb_h_slab9) / 0.2050409091
    d/dt(particles_tb_a_slab10) <- -(alp11 * particles_tb_a_slab11 - alp10 * particles_tb_a_slab10) / 0.2050409091
    d/dt(particles_tb_b_slab10) <- -(alp11 * particles_tb_b_slab11 - alp10 * particles_tb_b_slab10) / 0.2050409091
    d/dt(particles_tb_c_slab10) <- -(alp11 * particles_tb_c_slab11 - alp10 * particles_tb_c_slab10) / 0.2050409091
    d/dt(particles_tb_d_slab10) <- -(alp11 * particles_tb_d_slab11 - alp10 * particles_tb_d_slab10) / 0.2050409091
    d/dt(particles_tb_e_slab10) <- -(alp11 * particles_tb_e_slab11 - alp10 * particles_tb_e_slab10) / 0.2050409091
    d/dt(particles_tb_f_slab10) <- -(alp11 * particles_tb_f_slab11 - alp10 * particles_tb_f_slab10) / 0.2050409091
    d/dt(particles_tb_g_slab10) <- -(alp11 * particles_tb_g_slab11 - alp10 * particles_tb_g_slab10) / 0.2050409091
    d/dt(particles_tb_h_slab10) <- -(alp11 * particles_tb_h_slab11 - alp10 * particles_tb_h_slab10) / 0.2050409091
    d/dt(particles_tb_a_slab11) <- -(alp12 * particles_tb_a_slab12 - alp11 * particles_tb_a_slab11) / 0.2050409091
    d/dt(particles_tb_b_slab11) <- -(alp12 * particles_tb_b_slab12 - alp11 * particles_tb_b_slab11) / 0.2050409091
    d/dt(particles_tb_c_slab11) <- -(alp12 * particles_tb_c_slab12 - alp11 * particles_tb_c_slab11) / 0.2050409091
    d/dt(particles_tb_d_slab11) <- -(alp12 * particles_tb_d_slab12 - alp11 * particles_tb_d_slab11) / 0.2050409091
    d/dt(particles_tb_e_slab11) <- -(alp12 * particles_tb_e_slab12 - alp11 * particles_tb_e_slab11) / 0.2050409091
    d/dt(particles_tb_f_slab11) <- -(alp12 * particles_tb_f_slab12 - alp11 * particles_tb_f_slab11) / 0.2050409091
    d/dt(particles_tb_g_slab11) <- -(alp12 * particles_tb_g_slab12 - alp11 * particles_tb_g_slab11) / 0.2050409091
    d/dt(particles_tb_h_slab11) <- -(alp12 * particles_tb_h_slab12 - alp11 * particles_tb_h_slab11) / 0.2050409091
    d/dt(particles_tb_a_slab12) <- 0
    d/dt(particles_tb_b_slab12) <- 0
    d/dt(particles_tb_c_slab12) <- 0
    d/dt(particles_tb_d_slab12) <- 0
    d/dt(particles_tb_e_slab12) <- 0
    d/dt(particles_tb_f_slab12) <- 0
    d/dt(particles_tb_g_slab12) <- 0
    d/dt(particles_tb_h_slab12) <- 0
    # No mucociliary clearance in the alveolar region (Eq. S7).
    d/dt(particles_alv_a_slab1) <- 0
    d/dt(particles_alv_b_slab1) <- 0
    d/dt(particles_alv_c_slab1) <- 0
    d/dt(particles_alv_d_slab1) <- 0
    d/dt(particles_alv_e_slab1) <- 0
    d/dt(particles_alv_f_slab1) <- 0
    d/dt(particles_alv_g_slab1) <- 0
    d/dt(particles_alv_h_slab1) <- 0
    d/dt(particles_alv_a_slab2) <- 0
    d/dt(particles_alv_b_slab2) <- 0
    d/dt(particles_alv_c_slab2) <- 0
    d/dt(particles_alv_d_slab2) <- 0
    d/dt(particles_alv_e_slab2) <- 0
    d/dt(particles_alv_f_slab2) <- 0
    d/dt(particles_alv_g_slab2) <- 0
    d/dt(particles_alv_h_slab2) <- 0
    d/dt(particles_alv_a_slab3) <- 0
    d/dt(particles_alv_b_slab3) <- 0
    d/dt(particles_alv_c_slab3) <- 0
    d/dt(particles_alv_d_slab3) <- 0
    d/dt(particles_alv_e_slab3) <- 0
    d/dt(particles_alv_f_slab3) <- 0
    d/dt(particles_alv_g_slab3) <- 0
    d/dt(particles_alv_h_slab3) <- 0
    d/dt(particles_alv_a_slab4) <- 0
    d/dt(particles_alv_b_slab4) <- 0
    d/dt(particles_alv_c_slab4) <- 0
    d/dt(particles_alv_d_slab4) <- 0
    d/dt(particles_alv_e_slab4) <- 0
    d/dt(particles_alv_f_slab4) <- 0
    d/dt(particles_alv_g_slab4) <- 0
    d/dt(particles_alv_h_slab4) <- 0
    d/dt(particles_alv_a_slab5) <- 0
    d/dt(particles_alv_b_slab5) <- 0
    d/dt(particles_alv_c_slab5) <- 0
    d/dt(particles_alv_d_slab5) <- 0
    d/dt(particles_alv_e_slab5) <- 0
    d/dt(particles_alv_f_slab5) <- 0
    d/dt(particles_alv_g_slab5) <- 0
    d/dt(particles_alv_h_slab5) <- 0
    d/dt(particles_alv_a_slab6) <- 0
    d/dt(particles_alv_b_slab6) <- 0
    d/dt(particles_alv_c_slab6) <- 0
    d/dt(particles_alv_d_slab6) <- 0
    d/dt(particles_alv_e_slab6) <- 0
    d/dt(particles_alv_f_slab6) <- 0
    d/dt(particles_alv_g_slab6) <- 0
    d/dt(particles_alv_h_slab6) <- 0
    d/dt(particles_alv_a_slab7) <- 0
    d/dt(particles_alv_b_slab7) <- 0
    d/dt(particles_alv_c_slab7) <- 0
    d/dt(particles_alv_d_slab7) <- 0
    d/dt(particles_alv_e_slab7) <- 0
    d/dt(particles_alv_f_slab7) <- 0
    d/dt(particles_alv_g_slab7) <- 0
    d/dt(particles_alv_h_slab7) <- 0
    d/dt(particles_alv_a_slab8) <- 0
    d/dt(particles_alv_b_slab8) <- 0
    d/dt(particles_alv_c_slab8) <- 0
    d/dt(particles_alv_d_slab8) <- 0
    d/dt(particles_alv_e_slab8) <- 0
    d/dt(particles_alv_f_slab8) <- 0
    d/dt(particles_alv_g_slab8) <- 0
    d/dt(particles_alv_h_slab8) <- 0

    # Particle z-moment: same advection, plus the Eq. 10 dissolution
    # sink, which is handed to the epithelial lining fluid unchanged.
    d/dt(pmass_tb_a_slab1) <- -(alp2 * pmass_tb_a_slab2 - alp1 * pmass_tb_a_slab1) / 0.1025204545 +
      gamu1 * wu1a
    d/dt(pmass_tb_b_slab1) <- -(alp2 * pmass_tb_b_slab2 - alp1 * pmass_tb_b_slab1) / 0.1025204545 +
      gamu1 * wu1b
    d/dt(pmass_tb_c_slab1) <- -(alp2 * pmass_tb_c_slab2 - alp1 * pmass_tb_c_slab1) / 0.1025204545 +
      gamu1 * wu1c
    d/dt(pmass_tb_d_slab1) <- -(alp2 * pmass_tb_d_slab2 - alp1 * pmass_tb_d_slab1) / 0.1025204545 +
      gamu1 * wu1d
    d/dt(pmass_tb_e_slab1) <- -(alp2 * pmass_tb_e_slab2 - alp1 * pmass_tb_e_slab1) / 0.1025204545 +
      gamu1 * wu1e
    d/dt(pmass_tb_f_slab1) <- -(alp2 * pmass_tb_f_slab2 - alp1 * pmass_tb_f_slab1) / 0.1025204545 +
      gamu1 * wu1f
    d/dt(pmass_tb_g_slab1) <- -(alp2 * pmass_tb_g_slab2 - alp1 * pmass_tb_g_slab1) / 0.1025204545 +
      gamu1 * wu1g
    d/dt(pmass_tb_h_slab1) <- -(alp2 * pmass_tb_h_slab2 - alp1 * pmass_tb_h_slab1) / 0.1025204545 +
      gamu1 * wu1h
    d/dt(pmass_tb_a_slab2) <- -(alp3 * pmass_tb_a_slab3 - alp2 * pmass_tb_a_slab2) / 0.2050409091 +
      gamu2 * wu2a
    d/dt(pmass_tb_b_slab2) <- -(alp3 * pmass_tb_b_slab3 - alp2 * pmass_tb_b_slab2) / 0.2050409091 +
      gamu2 * wu2b
    d/dt(pmass_tb_c_slab2) <- -(alp3 * pmass_tb_c_slab3 - alp2 * pmass_tb_c_slab2) / 0.2050409091 +
      gamu2 * wu2c
    d/dt(pmass_tb_d_slab2) <- -(alp3 * pmass_tb_d_slab3 - alp2 * pmass_tb_d_slab2) / 0.2050409091 +
      gamu2 * wu2d
    d/dt(pmass_tb_e_slab2) <- -(alp3 * pmass_tb_e_slab3 - alp2 * pmass_tb_e_slab2) / 0.2050409091 +
      gamu2 * wu2e
    d/dt(pmass_tb_f_slab2) <- -(alp3 * pmass_tb_f_slab3 - alp2 * pmass_tb_f_slab2) / 0.2050409091 +
      gamu2 * wu2f
    d/dt(pmass_tb_g_slab2) <- -(alp3 * pmass_tb_g_slab3 - alp2 * pmass_tb_g_slab2) / 0.2050409091 +
      gamu2 * wu2g
    d/dt(pmass_tb_h_slab2) <- -(alp3 * pmass_tb_h_slab3 - alp2 * pmass_tb_h_slab2) / 0.2050409091 +
      gamu2 * wu2h
    d/dt(pmass_tb_a_slab3) <- -(alp4 * pmass_tb_a_slab4 - alp3 * pmass_tb_a_slab3) / 0.2050409091 +
      gamu3 * wu3a
    d/dt(pmass_tb_b_slab3) <- -(alp4 * pmass_tb_b_slab4 - alp3 * pmass_tb_b_slab3) / 0.2050409091 +
      gamu3 * wu3b
    d/dt(pmass_tb_c_slab3) <- -(alp4 * pmass_tb_c_slab4 - alp3 * pmass_tb_c_slab3) / 0.2050409091 +
      gamu3 * wu3c
    d/dt(pmass_tb_d_slab3) <- -(alp4 * pmass_tb_d_slab4 - alp3 * pmass_tb_d_slab3) / 0.2050409091 +
      gamu3 * wu3d
    d/dt(pmass_tb_e_slab3) <- -(alp4 * pmass_tb_e_slab4 - alp3 * pmass_tb_e_slab3) / 0.2050409091 +
      gamu3 * wu3e
    d/dt(pmass_tb_f_slab3) <- -(alp4 * pmass_tb_f_slab4 - alp3 * pmass_tb_f_slab3) / 0.2050409091 +
      gamu3 * wu3f
    d/dt(pmass_tb_g_slab3) <- -(alp4 * pmass_tb_g_slab4 - alp3 * pmass_tb_g_slab3) / 0.2050409091 +
      gamu3 * wu3g
    d/dt(pmass_tb_h_slab3) <- -(alp4 * pmass_tb_h_slab4 - alp3 * pmass_tb_h_slab3) / 0.2050409091 +
      gamu3 * wu3h
    d/dt(pmass_tb_a_slab4) <- -(alp5 * pmass_tb_a_slab5 - alp4 * pmass_tb_a_slab4) / 0.2050409091 +
      gamu4 * wu4a
    d/dt(pmass_tb_b_slab4) <- -(alp5 * pmass_tb_b_slab5 - alp4 * pmass_tb_b_slab4) / 0.2050409091 +
      gamu4 * wu4b
    d/dt(pmass_tb_c_slab4) <- -(alp5 * pmass_tb_c_slab5 - alp4 * pmass_tb_c_slab4) / 0.2050409091 +
      gamu4 * wu4c
    d/dt(pmass_tb_d_slab4) <- -(alp5 * pmass_tb_d_slab5 - alp4 * pmass_tb_d_slab4) / 0.2050409091 +
      gamu4 * wu4d
    d/dt(pmass_tb_e_slab4) <- -(alp5 * pmass_tb_e_slab5 - alp4 * pmass_tb_e_slab4) / 0.2050409091 +
      gamu4 * wu4e
    d/dt(pmass_tb_f_slab4) <- -(alp5 * pmass_tb_f_slab5 - alp4 * pmass_tb_f_slab4) / 0.2050409091 +
      gamu4 * wu4f
    d/dt(pmass_tb_g_slab4) <- -(alp5 * pmass_tb_g_slab5 - alp4 * pmass_tb_g_slab4) / 0.2050409091 +
      gamu4 * wu4g
    d/dt(pmass_tb_h_slab4) <- -(alp5 * pmass_tb_h_slab5 - alp4 * pmass_tb_h_slab4) / 0.2050409091 +
      gamu4 * wu4h
    d/dt(pmass_tb_a_slab5) <- -(alp6 * pmass_tb_a_slab6 - alp5 * pmass_tb_a_slab5) / 0.2050409091 +
      gamu5 * wu5a
    d/dt(pmass_tb_b_slab5) <- -(alp6 * pmass_tb_b_slab6 - alp5 * pmass_tb_b_slab5) / 0.2050409091 +
      gamu5 * wu5b
    d/dt(pmass_tb_c_slab5) <- -(alp6 * pmass_tb_c_slab6 - alp5 * pmass_tb_c_slab5) / 0.2050409091 +
      gamu5 * wu5c
    d/dt(pmass_tb_d_slab5) <- -(alp6 * pmass_tb_d_slab6 - alp5 * pmass_tb_d_slab5) / 0.2050409091 +
      gamu5 * wu5d
    d/dt(pmass_tb_e_slab5) <- -(alp6 * pmass_tb_e_slab6 - alp5 * pmass_tb_e_slab5) / 0.2050409091 +
      gamu5 * wu5e
    d/dt(pmass_tb_f_slab5) <- -(alp6 * pmass_tb_f_slab6 - alp5 * pmass_tb_f_slab5) / 0.2050409091 +
      gamu5 * wu5f
    d/dt(pmass_tb_g_slab5) <- -(alp6 * pmass_tb_g_slab6 - alp5 * pmass_tb_g_slab5) / 0.2050409091 +
      gamu5 * wu5g
    d/dt(pmass_tb_h_slab5) <- -(alp6 * pmass_tb_h_slab6 - alp5 * pmass_tb_h_slab5) / 0.2050409091 +
      gamu5 * wu5h
    d/dt(pmass_tb_a_slab6) <- -(alp7 * pmass_tb_a_slab7 - alp6 * pmass_tb_a_slab6) / 0.2050409091 +
      gamu6 * wu6a
    d/dt(pmass_tb_b_slab6) <- -(alp7 * pmass_tb_b_slab7 - alp6 * pmass_tb_b_slab6) / 0.2050409091 +
      gamu6 * wu6b
    d/dt(pmass_tb_c_slab6) <- -(alp7 * pmass_tb_c_slab7 - alp6 * pmass_tb_c_slab6) / 0.2050409091 +
      gamu6 * wu6c
    d/dt(pmass_tb_d_slab6) <- -(alp7 * pmass_tb_d_slab7 - alp6 * pmass_tb_d_slab6) / 0.2050409091 +
      gamu6 * wu6d
    d/dt(pmass_tb_e_slab6) <- -(alp7 * pmass_tb_e_slab7 - alp6 * pmass_tb_e_slab6) / 0.2050409091 +
      gamu6 * wu6e
    d/dt(pmass_tb_f_slab6) <- -(alp7 * pmass_tb_f_slab7 - alp6 * pmass_tb_f_slab6) / 0.2050409091 +
      gamu6 * wu6f
    d/dt(pmass_tb_g_slab6) <- -(alp7 * pmass_tb_g_slab7 - alp6 * pmass_tb_g_slab6) / 0.2050409091 +
      gamu6 * wu6g
    d/dt(pmass_tb_h_slab6) <- -(alp7 * pmass_tb_h_slab7 - alp6 * pmass_tb_h_slab6) / 0.2050409091 +
      gamu6 * wu6h
    d/dt(pmass_tb_a_slab7) <- -(alp8 * pmass_tb_a_slab8 - alp7 * pmass_tb_a_slab7) / 0.2050409091 +
      gamu7 * wu7a
    d/dt(pmass_tb_b_slab7) <- -(alp8 * pmass_tb_b_slab8 - alp7 * pmass_tb_b_slab7) / 0.2050409091 +
      gamu7 * wu7b
    d/dt(pmass_tb_c_slab7) <- -(alp8 * pmass_tb_c_slab8 - alp7 * pmass_tb_c_slab7) / 0.2050409091 +
      gamu7 * wu7c
    d/dt(pmass_tb_d_slab7) <- -(alp8 * pmass_tb_d_slab8 - alp7 * pmass_tb_d_slab7) / 0.2050409091 +
      gamu7 * wu7d
    d/dt(pmass_tb_e_slab7) <- -(alp8 * pmass_tb_e_slab8 - alp7 * pmass_tb_e_slab7) / 0.2050409091 +
      gamu7 * wu7e
    d/dt(pmass_tb_f_slab7) <- -(alp8 * pmass_tb_f_slab8 - alp7 * pmass_tb_f_slab7) / 0.2050409091 +
      gamu7 * wu7f
    d/dt(pmass_tb_g_slab7) <- -(alp8 * pmass_tb_g_slab8 - alp7 * pmass_tb_g_slab7) / 0.2050409091 +
      gamu7 * wu7g
    d/dt(pmass_tb_h_slab7) <- -(alp8 * pmass_tb_h_slab8 - alp7 * pmass_tb_h_slab7) / 0.2050409091 +
      gamu7 * wu7h
    d/dt(pmass_tb_a_slab8) <- -(alp9 * pmass_tb_a_slab9 - alp8 * pmass_tb_a_slab8) / 0.2050409091 +
      gamu8 * wu8a
    d/dt(pmass_tb_b_slab8) <- -(alp9 * pmass_tb_b_slab9 - alp8 * pmass_tb_b_slab8) / 0.2050409091 +
      gamu8 * wu8b
    d/dt(pmass_tb_c_slab8) <- -(alp9 * pmass_tb_c_slab9 - alp8 * pmass_tb_c_slab8) / 0.2050409091 +
      gamu8 * wu8c
    d/dt(pmass_tb_d_slab8) <- -(alp9 * pmass_tb_d_slab9 - alp8 * pmass_tb_d_slab8) / 0.2050409091 +
      gamu8 * wu8d
    d/dt(pmass_tb_e_slab8) <- -(alp9 * pmass_tb_e_slab9 - alp8 * pmass_tb_e_slab8) / 0.2050409091 +
      gamu8 * wu8e
    d/dt(pmass_tb_f_slab8) <- -(alp9 * pmass_tb_f_slab9 - alp8 * pmass_tb_f_slab8) / 0.2050409091 +
      gamu8 * wu8f
    d/dt(pmass_tb_g_slab8) <- -(alp9 * pmass_tb_g_slab9 - alp8 * pmass_tb_g_slab8) / 0.2050409091 +
      gamu8 * wu8g
    d/dt(pmass_tb_h_slab8) <- -(alp9 * pmass_tb_h_slab9 - alp8 * pmass_tb_h_slab8) / 0.2050409091 +
      gamu8 * wu8h
    d/dt(pmass_tb_a_slab9) <- -(alp10 * pmass_tb_a_slab10 - alp9 * pmass_tb_a_slab9) / 0.2050409091 +
      gamu9 * wu9a
    d/dt(pmass_tb_b_slab9) <- -(alp10 * pmass_tb_b_slab10 - alp9 * pmass_tb_b_slab9) / 0.2050409091 +
      gamu9 * wu9b
    d/dt(pmass_tb_c_slab9) <- -(alp10 * pmass_tb_c_slab10 - alp9 * pmass_tb_c_slab9) / 0.2050409091 +
      gamu9 * wu9c
    d/dt(pmass_tb_d_slab9) <- -(alp10 * pmass_tb_d_slab10 - alp9 * pmass_tb_d_slab9) / 0.2050409091 +
      gamu9 * wu9d
    d/dt(pmass_tb_e_slab9) <- -(alp10 * pmass_tb_e_slab10 - alp9 * pmass_tb_e_slab9) / 0.2050409091 +
      gamu9 * wu9e
    d/dt(pmass_tb_f_slab9) <- -(alp10 * pmass_tb_f_slab10 - alp9 * pmass_tb_f_slab9) / 0.2050409091 +
      gamu9 * wu9f
    d/dt(pmass_tb_g_slab9) <- -(alp10 * pmass_tb_g_slab10 - alp9 * pmass_tb_g_slab9) / 0.2050409091 +
      gamu9 * wu9g
    d/dt(pmass_tb_h_slab9) <- -(alp10 * pmass_tb_h_slab10 - alp9 * pmass_tb_h_slab9) / 0.2050409091 +
      gamu9 * wu9h
    d/dt(pmass_tb_a_slab10) <- -(alp11 * pmass_tb_a_slab11 - alp10 * pmass_tb_a_slab10) / 0.2050409091 +
      gamu10 * wu10a
    d/dt(pmass_tb_b_slab10) <- -(alp11 * pmass_tb_b_slab11 - alp10 * pmass_tb_b_slab10) / 0.2050409091 +
      gamu10 * wu10b
    d/dt(pmass_tb_c_slab10) <- -(alp11 * pmass_tb_c_slab11 - alp10 * pmass_tb_c_slab10) / 0.2050409091 +
      gamu10 * wu10c
    d/dt(pmass_tb_d_slab10) <- -(alp11 * pmass_tb_d_slab11 - alp10 * pmass_tb_d_slab10) / 0.2050409091 +
      gamu10 * wu10d
    d/dt(pmass_tb_e_slab10) <- -(alp11 * pmass_tb_e_slab11 - alp10 * pmass_tb_e_slab10) / 0.2050409091 +
      gamu10 * wu10e
    d/dt(pmass_tb_f_slab10) <- -(alp11 * pmass_tb_f_slab11 - alp10 * pmass_tb_f_slab10) / 0.2050409091 +
      gamu10 * wu10f
    d/dt(pmass_tb_g_slab10) <- -(alp11 * pmass_tb_g_slab11 - alp10 * pmass_tb_g_slab10) / 0.2050409091 +
      gamu10 * wu10g
    d/dt(pmass_tb_h_slab10) <- -(alp11 * pmass_tb_h_slab11 - alp10 * pmass_tb_h_slab10) / 0.2050409091 +
      gamu10 * wu10h
    d/dt(pmass_tb_a_slab11) <- -(alp12 * pmass_tb_a_slab12 - alp11 * pmass_tb_a_slab11) / 0.2050409091 +
      gamu11 * wu11a
    d/dt(pmass_tb_b_slab11) <- -(alp12 * pmass_tb_b_slab12 - alp11 * pmass_tb_b_slab11) / 0.2050409091 +
      gamu11 * wu11b
    d/dt(pmass_tb_c_slab11) <- -(alp12 * pmass_tb_c_slab12 - alp11 * pmass_tb_c_slab11) / 0.2050409091 +
      gamu11 * wu11c
    d/dt(pmass_tb_d_slab11) <- -(alp12 * pmass_tb_d_slab12 - alp11 * pmass_tb_d_slab11) / 0.2050409091 +
      gamu11 * wu11d
    d/dt(pmass_tb_e_slab11) <- -(alp12 * pmass_tb_e_slab12 - alp11 * pmass_tb_e_slab11) / 0.2050409091 +
      gamu11 * wu11e
    d/dt(pmass_tb_f_slab11) <- -(alp12 * pmass_tb_f_slab12 - alp11 * pmass_tb_f_slab11) / 0.2050409091 +
      gamu11 * wu11f
    d/dt(pmass_tb_g_slab11) <- -(alp12 * pmass_tb_g_slab12 - alp11 * pmass_tb_g_slab11) / 0.2050409091 +
      gamu11 * wu11g
    d/dt(pmass_tb_h_slab11) <- -(alp12 * pmass_tb_h_slab12 - alp11 * pmass_tb_h_slab11) / 0.2050409091 +
      gamu11 * wu11h
    d/dt(pmass_tb_a_slab12) <- gamu12 * wu12a
    d/dt(pmass_tb_b_slab12) <- gamu12 * wu12b
    d/dt(pmass_tb_c_slab12) <- gamu12 * wu12c
    d/dt(pmass_tb_d_slab12) <- gamu12 * wu12d
    d/dt(pmass_tb_e_slab12) <- gamu12 * wu12e
    d/dt(pmass_tb_f_slab12) <- gamu12 * wu12f
    d/dt(pmass_tb_g_slab12) <- gamu12 * wu12g
    d/dt(pmass_tb_h_slab12) <- gamu12 * wu12h
    d/dt(pmass_alv_a_slab1) <- gaml1 * wl1a
    d/dt(pmass_alv_b_slab1) <- gaml1 * wl1b
    d/dt(pmass_alv_c_slab1) <- gaml1 * wl1c
    d/dt(pmass_alv_d_slab1) <- gaml1 * wl1d
    d/dt(pmass_alv_e_slab1) <- gaml1 * wl1e
    d/dt(pmass_alv_f_slab1) <- gaml1 * wl1f
    d/dt(pmass_alv_g_slab1) <- gaml1 * wl1g
    d/dt(pmass_alv_h_slab1) <- gaml1 * wl1h
    d/dt(pmass_alv_a_slab2) <- gaml2 * wl2a
    d/dt(pmass_alv_b_slab2) <- gaml2 * wl2b
    d/dt(pmass_alv_c_slab2) <- gaml2 * wl2c
    d/dt(pmass_alv_d_slab2) <- gaml2 * wl2d
    d/dt(pmass_alv_e_slab2) <- gaml2 * wl2e
    d/dt(pmass_alv_f_slab2) <- gaml2 * wl2f
    d/dt(pmass_alv_g_slab2) <- gaml2 * wl2g
    d/dt(pmass_alv_h_slab2) <- gaml2 * wl2h
    d/dt(pmass_alv_a_slab3) <- gaml3 * wl3a
    d/dt(pmass_alv_b_slab3) <- gaml3 * wl3b
    d/dt(pmass_alv_c_slab3) <- gaml3 * wl3c
    d/dt(pmass_alv_d_slab3) <- gaml3 * wl3d
    d/dt(pmass_alv_e_slab3) <- gaml3 * wl3e
    d/dt(pmass_alv_f_slab3) <- gaml3 * wl3f
    d/dt(pmass_alv_g_slab3) <- gaml3 * wl3g
    d/dt(pmass_alv_h_slab3) <- gaml3 * wl3h
    d/dt(pmass_alv_a_slab4) <- gaml4 * wl4a
    d/dt(pmass_alv_b_slab4) <- gaml4 * wl4b
    d/dt(pmass_alv_c_slab4) <- gaml4 * wl4c
    d/dt(pmass_alv_d_slab4) <- gaml4 * wl4d
    d/dt(pmass_alv_e_slab4) <- gaml4 * wl4e
    d/dt(pmass_alv_f_slab4) <- gaml4 * wl4f
    d/dt(pmass_alv_g_slab4) <- gaml4 * wl4g
    d/dt(pmass_alv_h_slab4) <- gaml4 * wl4h
    d/dt(pmass_alv_a_slab5) <- gaml5 * wl5a
    d/dt(pmass_alv_b_slab5) <- gaml5 * wl5b
    d/dt(pmass_alv_c_slab5) <- gaml5 * wl5c
    d/dt(pmass_alv_d_slab5) <- gaml5 * wl5d
    d/dt(pmass_alv_e_slab5) <- gaml5 * wl5e
    d/dt(pmass_alv_f_slab5) <- gaml5 * wl5f
    d/dt(pmass_alv_g_slab5) <- gaml5 * wl5g
    d/dt(pmass_alv_h_slab5) <- gaml5 * wl5h
    d/dt(pmass_alv_a_slab6) <- gaml6 * wl6a
    d/dt(pmass_alv_b_slab6) <- gaml6 * wl6b
    d/dt(pmass_alv_c_slab6) <- gaml6 * wl6c
    d/dt(pmass_alv_d_slab6) <- gaml6 * wl6d
    d/dt(pmass_alv_e_slab6) <- gaml6 * wl6e
    d/dt(pmass_alv_f_slab6) <- gaml6 * wl6f
    d/dt(pmass_alv_g_slab6) <- gaml6 * wl6g
    d/dt(pmass_alv_h_slab6) <- gaml6 * wl6h
    d/dt(pmass_alv_a_slab7) <- gaml7 * wl7a
    d/dt(pmass_alv_b_slab7) <- gaml7 * wl7b
    d/dt(pmass_alv_c_slab7) <- gaml7 * wl7c
    d/dt(pmass_alv_d_slab7) <- gaml7 * wl7d
    d/dt(pmass_alv_e_slab7) <- gaml7 * wl7e
    d/dt(pmass_alv_f_slab7) <- gaml7 * wl7f
    d/dt(pmass_alv_g_slab7) <- gaml7 * wl7g
    d/dt(pmass_alv_h_slab7) <- gaml7 * wl7h
    d/dt(pmass_alv_a_slab8) <- gaml8 * wl8a
    d/dt(pmass_alv_b_slab8) <- gaml8 * wl8b
    d/dt(pmass_alv_c_slab8) <- gaml8 * wl8c
    d/dt(pmass_alv_d_slab8) <- gaml8 * wl8d
    d/dt(pmass_alv_e_slab8) <- gaml8 * wl8e
    d/dt(pmass_alv_f_slab8) <- gaml8 * wl8f
    d/dt(pmass_alv_g_slab8) <- gaml8 * wl8g
    d/dt(pmass_alv_h_slab8) <- gaml8 * wl8h

    # ---- default initial state: case study 1, MCC + wide PSD ----
    # (d ~ N(3, 0.6) um; LDD 100 ug; Cs 100 nM, Table 1). Override with
    # rxSolve(inits = ...) for another particle size distribution or
    # lung-deposited dose - the vignette derives those from the
    # Appendix S1 deposition model.
    depot(0) <- 60.57435708
    particles_tb_a_slab1(0) <- 0.0460090653
    pmass_tb_a_slab1(0) <- 1.55442346e-12
    particles_tb_b_slab1(0) <- 0.3487095539
    pmass_tb_b_slab1(0) <- 2.482512351e-11
    particles_tb_c_slab1(0) <- 1.291032248
    pmass_tb_c_slab1(0) <- 1.579951144e-10
    particles_tb_d_slab1(0) <- 2.343963533
    pmass_tb_d_slab1(0) <- 4.391360742e-10
    particles_tb_e_slab1(0) <- 2.100071805
    pmass_tb_e_slab1(0) <- 5.588241854e-10
    particles_tb_f_slab1(0) <- 0.930102332
    pmass_tb_f_slab1(0) <- 3.33561582e-10
    particles_tb_g_slab1(0) <- 0.2029085333
    pmass_tb_g_slab1(0) <- 9.434057881e-11
    particles_tb_h_slab1(0) <- 0.02165372898
    pmass_tb_h_slab1(0) <- 1.266819272e-11
    particles_tb_a_slab2(0) <- 0.03489129098
    pmass_tb_a_slab2(0) <- 1.178807718e-12
    particles_tb_b_slab2(0) <- 0.2570499488
    pmass_tb_b_slab2(0) <- 1.829974733e-11
    particles_tb_c_slab2(0) <- 0.942543362
    pmass_tb_c_slab2(0) <- 1.153474257e-10
    particles_tb_d_slab2(0) <- 1.705156157
    pmass_tb_d_slab2(0) <- 3.194570095e-10
    particles_tb_e_slab2(0) <- 1.526161634
    pmass_tb_e_slab2(0) <- 4.06108034e-10
    particles_tb_f_slab2(0) <- 0.6760686733
    pmass_tb_f_slab2(0) <- 2.424577689e-10
    particles_tb_g_slab2(0) <- 0.1476219855
    pmass_tb_g_slab2(0) <- 6.863557353e-11
    particles_tb_h_slab2(0) <- 0.01577464563
    pmass_tb_h_slab2(0) <- 9.22872227e-12
    particles_tb_a_slab3(0) <- 0.02947322283
    pmass_tb_a_slab3(0) <- 9.957574386e-13
    particles_tb_b_slab3(0) <- 0.2126357501
    pmass_tb_b_slab3(0) <- 1.513783807e-11
    particles_tb_c_slab3(0) <- 0.774024566
    pmass_tb_c_slab3(0) <- 9.47242798e-11
    particles_tb_d_slab3(0) <- 1.396570326
    pmass_tb_d_slab3(0) <- 2.616441774e-10
    particles_tb_e_slab3(0) <- 1.249124528
    pmass_tb_e_slab3(0) <- 3.323891092e-10
    particles_tb_f_slab3(0) <- 0.5535172303
    pmass_tb_f_slab3(0) <- 1.985072789e-10
    particles_tb_g_slab3(0) <- 0.1209663738
    pmass_tb_g_slab3(0) <- 5.624227595e-11
    particles_tb_h_slab3(0) <- 0.01294185327
    pmass_tb_h_slab3(0) <- 7.571439148e-12
    particles_tb_a_slab4(0) <- 0.02945036611
    pmass_tb_a_slab4(0) <- 9.949852205e-13
    particles_tb_b_slab4(0) <- 0.2124262634
    pmass_tb_b_slab4(0) <- 1.512292442e-11
    particles_tb_c_slab4(0) <- 0.7731997053
    pmass_tb_c_slab4(0) <- 9.462333425e-11
    particles_tb_d_slab4(0) <- 1.395031653
    pmass_tb_d_slab4(0) <- 2.613559106e-10
    particles_tb_e_slab4(0) <- 1.247725712
    pmass_tb_e_slab4(0) <- 3.320168877e-10
    particles_tb_f_slab4(0) <- 0.5528918374
    pmass_tb_f_slab4(0) <- 1.982829949e-10
    particles_tb_g_slab4(0) <- 0.1208289575
    pmass_tb_g_slab4(0) <- 5.617838539e-11
    particles_tb_h_slab4(0) <- 0.01292709823
    pmass_tb_h_slab4(0) <- 7.562806931e-12
    particles_tb_a_slab5(0) <- 0.03451822607
    pmass_tb_a_slab5(0) <- 1.166203661e-12
    particles_tb_b_slab5(0) <- 0.2533807943
    pmass_tb_b_slab5(0) <- 1.803853506e-11
    particles_tb_c_slab5(0) <- 0.9277926257
    pmass_tb_c_slab5(0) <- 1.135422468e-10
    particles_tb_d_slab5(0) <- 1.677365749
    pmass_tb_d_slab5(0) <- 3.142505417e-10
    particles_tb_e_slab5(0) <- 1.500730415
    pmass_tb_e_slab5(0) <- 3.99340846e-10
    particles_tb_f_slab5(0) <- 0.6646363292
    pmass_tb_f_slab5(0) <- 2.383577998e-10
    particles_tb_g_slab5(0) <- 0.1450969956
    pmass_tb_g_slab5(0) <- 6.746160119e-11
    particles_tb_h_slab5(0) <- 0.01550212684
    pmass_tb_h_slab5(0) <- 9.069289199e-12
    particles_tb_a_slab6(0) <- 0.044372308
    pmass_tb_a_slab6(0) <- 1.499125359e-12
    particles_tb_b_slab6(0) <- 0.3324586484
    pmass_tb_b_slab6(0) <- 2.36681987e-11
    particles_tb_c_slab6(0) <- 1.225527173
    pmass_tb_c_slab6(0) <- 1.49978675e-10
    particles_tb_d_slab6(0) <- 2.220398225
    pmass_tb_d_slab6(0) <- 4.159864034e-10
    particles_tb_e_slab6(0) <- 1.986903861
    pmass_tb_e_slab6(0) <- 5.287104606e-10
    particles_tb_f_slab6(0) <- 0.8791945402
    pmass_tb_f_slab6(0) <- 3.153045763e-10
    particles_tb_g_slab6(0) <- 0.1916577476
    pmass_tb_g_slab6(0) <- 8.910962268e-11
    particles_tb_h_slab6(0) <- 0.02043868544
    pmass_tb_h_slab6(0) <- 1.195734953e-11
    particles_tb_a_slab7(0) <- 0.05870833075
    pmass_tb_a_slab7(0) <- 1.983470128e-12
    particles_tb_b_slab7(0) <- 0.4466206452
    pmass_tb_b_slab7(0) <- 3.179555179e-11
    particles_tb_c_slab7(0) <-  1.65413264
    pmass_tb_c_slab7(0) <- 2.024309432e-10
    particles_tb_d_slab7(0) <- 3.000964404
    pmass_tb_d_slab7(0) <- 5.622236477e-10
    particles_tb_e_slab7(0) <- 2.685019846
    pmass_tb_e_slab7(0) <- 7.144774879e-10
    particles_tb_f_slab7(0) <- 1.187014042
    pmass_tb_f_slab7(0) <- 4.256975474e-10
    particles_tb_g_slab7(0) <- 0.2583992747
    pmass_tb_g_slab7(0) <- 1.201405222e-10
    particles_tb_h_slab7(0) <- 0.0275086045
    pmass_tb_h_slab7(0) <- 1.609350073e-11
    particles_tb_a_slab8(0) <- 0.09987813767
    pmass_tb_a_slab8(0) <- 3.374398487e-12
    particles_tb_b_slab8(0) <- 0.7596789506
    pmass_tb_b_slab8(0) <- 5.408261279e-11
    particles_tb_c_slab8(0) <- 2.810395812
    pmass_tb_c_slab8(0) <- 3.439331656e-10
    particles_tb_d_slab8(0) <- 5.089303831
    pmass_tb_d_slab8(0) <- 9.534691446e-10
    particles_tb_e_slab8(0) <- 4.542300842
    pmass_tb_e_slab8(0) <- 1.208695608e-09
    particles_tb_f_slab8(0) <- 2.002083969
    pmass_tb_f_slab8(0) <- 7.180051841e-10
    particles_tb_g_slab8(0) <- 0.4343282557
    pmass_tb_g_slab8(0) <- 2.0193719e-10
    particles_tb_h_slab8(0) <- 0.04605963725
    pmass_tb_h_slab8(0) <- 2.694650707e-11
    particles_tb_a_slab9(0) <- 0.3291643946
    pmass_tb_a_slab9(0) <- 1.11208705e-11
    particles_tb_b_slab9(0) <- 2.500411804
    pmass_tb_b_slab9(0) <- 1.780078325e-10
    particles_tb_c_slab9(0) <- 9.230587807
    pmass_tb_c_slab9(0) <- 1.12962924e-09
    particles_tb_d_slab9(0) <- 16.66193535
    pmass_tb_d_slab9(0) <- 3.121574536e-09
    particles_tb_e_slab9(0) <- 14.80595705
    pmass_tb_e_slab9(0) <- 3.939830471e-09
    particles_tb_f_slab9(0) <- 6.490086253
    pmass_tb_f_slab9(0) <- 2.327532534e-09
    particles_tb_g_slab9(0) <- 1.398814324
    pmass_tb_g_slab9(0) <- 6.503666989e-10
    particles_tb_h_slab9(0) <- 0.1472467661
    pmass_tb_h_slab9(0) <- 8.614453481e-11
    particles_tb_a_slab10(0) <- 0.5716259938
    pmass_tb_a_slab10(0) <- 1.931247352e-11
    particles_tb_b_slab10(0) <- 4.327512782
    pmass_tb_b_slab10(0) <- 3.080817205e-10
    particles_tb_c_slab10(0) <-  15.9060954
    pmass_tb_c_slab10(0) <- 1.946570558e-09
    particles_tb_d_slab10(0) <-  28.5188362
    pmass_tb_d_slab10(0) <- 5.342937121e-09
    particles_tb_e_slab10(0) <- 25.10028706
    pmass_tb_e_slab10(0) <- 6.679127557e-09
    particles_tb_f_slab10(0) <- 10.86689542
    pmass_tb_f_slab10(0) <- 3.897182817e-09
    particles_tb_g_slab10(0) <- 2.307382151
    pmass_tb_g_slab10(0) <- 1.072797502e-09
    particles_tb_h_slab10(0) <- 0.2387362818
    pmass_tb_h_slab10(0) <- 1.396691179e-10
    particles_tb_a_slab11(0) <- 0.8825559546
    pmass_tb_a_slab11(0) <- 2.981729083e-11
    particles_tb_b_slab11(0) <- 6.574073586
    pmass_tb_b_slab11(0) <- 4.680175434e-10
    particles_tb_c_slab11(0) <- 23.87157865
    pmass_tb_c_slab11(0) <- 2.921377685e-09
    particles_tb_d_slab11(0) <- 42.17493586
    pmass_tb_d_slab11(0) <- 7.901375386e-09
    particles_tb_e_slab11(0) <- 36.41312697
    pmass_tb_e_slab11(0) <- 9.689447742e-09
    particles_tb_f_slab11(0) <- 15.38832714
    pmass_tb_f_slab11(0) <- 5.51869893e-09
    particles_tb_g_slab11(0) <- 3.174072073
    pmass_tb_g_slab11(0) <- 1.475757533e-09
    particles_tb_h_slab11(0) <- 0.3175794377
    pmass_tb_h_slab11(0) <- 1.857951359e-10
    particles_tb_a_slab12(0) <- 2.751424875
    pmass_tb_a_slab12(0) <- 9.295731931e-11
    particles_tb_b_slab12(0) <- 17.34474544
    pmass_tb_b_slab12(0) <- 1.234796819e-09
    particles_tb_c_slab12(0) <- 58.39072386
    pmass_tb_c_slab12(0) <- 7.145792921e-09
    particles_tb_d_slab12(0) <- 98.03141091
    pmass_tb_d_slab12(0) <- 1.836595507e-08
    particles_tb_e_slab12(0) <- 80.83894126
    pmass_tb_e_slab12(0) <- 2.15110528e-08
    particles_tb_f_slab12(0) <- 32.57760665
    pmass_tb_f_slab12(0) <- 1.168327144e-08
    particles_tb_g_slab12(0) <- 6.383082397
    pmass_tb_g_slab12(0) <- 2.967759306e-09
    particles_tb_h_slab12(0) <- 0.6037925269
    pmass_tb_h_slab12(0) <- 3.532398553e-10
    particles_alv_a_slab1(0) <- 2.751424875
    pmass_alv_a_slab1(0) <- 9.295731931e-11
    particles_alv_b_slab1(0) <- 17.34474544
    pmass_alv_b_slab1(0) <- 1.234796819e-09
    particles_alv_c_slab1(0) <- 58.39072386
    pmass_alv_c_slab1(0) <- 7.145792921e-09
    particles_alv_d_slab1(0) <- 98.03141091
    pmass_alv_d_slab1(0) <- 1.836595507e-08
    particles_alv_e_slab1(0) <- 80.83894126
    pmass_alv_e_slab1(0) <- 2.15110528e-08
    particles_alv_f_slab1(0) <- 32.57760665
    pmass_alv_f_slab1(0) <- 1.168327144e-08
    particles_alv_g_slab1(0) <- 6.383082397
    pmass_alv_g_slab1(0) <- 2.967759306e-09
    particles_alv_h_slab1(0) <- 0.6037925269
    pmass_alv_h_slab1(0) <- 3.532398553e-10
    particles_alv_a_slab2(0) <- 3.709705906
    pmass_alv_a_slab2(0) <- 1.253329937e-10
    particles_alv_b_slab2(0) <- 23.18710822
    pmass_alv_b_slab2(0) <- 1.650722841e-09
    particles_alv_c_slab2(0) <-  77.6503285
    pmass_alv_c_slab2(0) <- 9.502762272e-09
    particles_alv_d_slab2(0) <- 129.7732639
    pmass_alv_d_slab2(0) <- 2.431271684e-08
    particles_alv_e_slab2(0) <- 106.5090664
    pmass_alv_e_slab2(0) <- 2.834181293e-08
    particles_alv_f_slab2(0) <- 42.70173204
    pmass_alv_f_slab2(0) <- 1.531407546e-08
    particles_alv_g_slab2(0) <- 8.319891556
    pmass_alv_g_slab2(0) <- 3.86826208e-09
    particles_alv_h_slab2(0) <- 0.7822513577
    pmass_alv_h_slab2(0) <- 4.576445453e-10
    particles_alv_a_slab3(0) <- 5.518234975
    pmass_alv_a_slab3(0) <- 1.864344308e-10
    particles_alv_b_slab3(0) <- 34.61794963
    pmass_alv_b_slab3(0) <- 2.464500516e-09
    particles_alv_c_slab3(0) <- 115.8316356
    pmass_alv_c_slab3(0) <- 1.417534887e-08
    particles_alv_d_slab3(0) <- 192.9547246
    pmass_alv_d_slab3(0) <- 3.614961541e-08
    particles_alv_e_slab3(0) <- 157.5887923
    pmass_alv_e_slab3(0) <- 4.193400828e-08
    particles_alv_f_slab3(0) <- 62.79202835
    pmass_alv_f_slab3(0) <- 2.251903645e-08
    particles_alv_g_slab3(0) <- 12.14709125
    pmass_alv_g_slab3(0) <- 5.647685685e-09
    particles_alv_h_slab3(0) <- 1.133071402
    pmass_alv_h_slab3(0) <- 6.628866048e-10
    particles_alv_a_slab4(0) <- 9.334575145
    pmass_alv_a_slab4(0) <- 3.153700798e-10
    particles_alv_b_slab4(0) <- 59.60597014
    pmass_alv_b_slab4(0) <- 4.243432835e-09
    particles_alv_c_slab4(0) <- 200.1570648
    pmass_alv_c_slab4(0) <- 2.449500267e-08
    particles_alv_d_slab4(0) <-  332.476614
    pmass_alv_d_slab4(0) <- 6.22887144e-08
    particles_alv_e_slab4(0) <- 269.7416215
    pmass_alv_e_slab4(0) <- 7.177761327e-08
    particles_alv_f_slab4(0) <- 106.4900294
    pmass_alv_f_slab4(0) <- 3.819040276e-08
    particles_alv_g_slab4(0) <- 20.37189242
    pmass_alv_g_slab4(0) <- 9.47173631e-09
    particles_alv_h_slab4(0) <- 1.876476429
    pmass_alv_h_slab4(0) <- 1.097804681e-09
    particles_alv_a_slab5(0) <- 17.92348213
    pmass_alv_a_slab5(0) <- 6.055476444e-10
    particles_alv_b_slab5(0) <- 117.0761683
    pmass_alv_b_slab5(0) <- 8.334817062e-09
    particles_alv_c_slab5(0) <- 393.8753418
    pmass_alv_c_slab5(0) <- 4.820203353e-08
    particles_alv_d_slab5(0) <- 649.0300969
    pmass_alv_d_slab5(0) <- 1.215942675e-07
    particles_alv_e_slab5(0) <- 519.0092865
    pmass_alv_e_slab5(0) <- 1.381071547e-07
    particles_alv_f_slab5(0) <- 200.9769875
    pmass_alv_f_slab5(0) <- 7.20761572e-08
    particles_alv_g_slab5(0) <- 37.57000691
    pmass_alv_g_slab5(0) <- 1.746785184e-08
    particles_alv_h_slab5(0) <- 3.371717754
    pmass_alv_h_slab5(0) <- 1.972573423e-09
    particles_alv_a_slab6(0) <- 44.62437773
    pmass_alv_a_slab6(0) <- 1.507641574e-09
    particles_alv_b_slab6(0) <- 300.9619348
    pmass_alv_b_slab6(0) <- 2.142590337e-08
    particles_alv_c_slab6(0) <- 1004.874431
    pmass_alv_c_slab6(0) <- 1.229754338e-07
    particles_alv_d_slab6(0) <- 1608.472532
    pmass_alv_d_slab6(0) <- 3.013435591e-07
    particles_alv_e_slab6(0) <- 1227.995971
    pmass_alv_e_slab6(0) <- 3.267668496e-07
    particles_alv_f_slab6(0) <-  446.830291
    pmass_alv_f_slab6(0) <- 1.602462585e-07
    particles_alv_g_slab6(0) <- 77.37934674
    pmass_alv_g_slab6(0) <- 3.597686229e-08
    particles_alv_h_slab6(0) <-  6.35554941
    pmass_alv_h_slab6(0) <- 3.718219842e-09
    particles_alv_a_slab7(0) <- 37.59696249
    pmass_alv_a_slab7(0) <- 1.270219252e-09
    particles_alv_b_slab7(0) <- 251.5533352
    pmass_alv_b_slab7(0) <- 1.790843568e-08
    particles_alv_c_slab7(0) <- 815.0217617
    pmass_alv_c_slab7(0) <- 9.974147176e-08
    particles_alv_d_slab7(0) <- 1238.175468
    pmass_alv_d_slab7(0) <- 2.31969272e-07
    particles_alv_e_slab7(0) <- 871.6438057
    pmass_alv_e_slab7(0) <- 2.319423738e-07
    particles_alv_f_slab7(0) <- 281.1311402
    pmass_alv_f_slab7(0) <- 1.008217533e-07
    particles_alv_g_slab7(0) <- 40.87418519
    pmass_alv_g_slab7(0) <- 1.900410114e-08
    particles_alv_h_slab7(0) <- 2.614699479
    pmass_alv_h_slab7(0) <- 1.529691118e-09
    particles_alv_a_slab8(0) <-           0
    pmass_alv_a_slab8(0) <-           0
    particles_alv_b_slab8(0) <-           0
    pmass_alv_b_slab8(0) <-           0
    particles_alv_c_slab8(0) <-           0
    pmass_alv_c_slab8(0) <-           0
    particles_alv_d_slab8(0) <-           0
    pmass_alv_d_slab8(0) <-           0
    particles_alv_e_slab8(0) <-           0
    pmass_alv_e_slab8(0) <-           0
    particles_alv_f_slab8(0) <-           0
    pmass_alv_f_slab8(0) <-           0
    particles_alv_g_slab8(0) <-           0
    pmass_alv_g_slab8(0) <-           0
    particles_alv_h_slab8(0) <-           0
    pmass_alv_h_slab8(0) <-           0

    # ---- outputs ----
    alung <- elf_tb_slab1 + epithelium_tb_slab1 + subepithelium_tb_slab1 +
      elf_tb_slab2 + epithelium_tb_slab2 + subepithelium_tb_slab2 +
      elf_tb_slab3 + epithelium_tb_slab3 + subepithelium_tb_slab3 +
      elf_tb_slab4 + epithelium_tb_slab4 + subepithelium_tb_slab4 +
      elf_tb_slab5 + epithelium_tb_slab5 + subepithelium_tb_slab5 +
      elf_tb_slab6 + epithelium_tb_slab6 + subepithelium_tb_slab6 +
      elf_tb_slab7 + epithelium_tb_slab7 + subepithelium_tb_slab7 +
      elf_tb_slab8 + epithelium_tb_slab8 + subepithelium_tb_slab8 +
      elf_tb_slab9 + epithelium_tb_slab9 + subepithelium_tb_slab9 +
      elf_tb_slab10 + epithelium_tb_slab10 + subepithelium_tb_slab10 +
      elf_tb_slab11 + epithelium_tb_slab11 + subepithelium_tb_slab11 +
      elf_tb_slab12 + epithelium_tb_slab12 + subepithelium_tb_slab12 +
      elf_alv_slab1 + epithelium_alv_slab1 + subepithelium_alv_slab1 +
      elf_alv_slab2 + epithelium_alv_slab2 + subepithelium_alv_slab2 +
      elf_alv_slab3 + epithelium_alv_slab3 + subepithelium_alv_slab3 +
      elf_alv_slab4 + epithelium_alv_slab4 + subepithelium_alv_slab4 +
      elf_alv_slab5 + epithelium_alv_slab5 + subepithelium_alv_slab5 +
      elf_alv_slab6 + epithelium_alv_slab6 + subepithelium_alv_slab6 +
      elf_alv_slab7 + epithelium_alv_slab7 + subepithelium_alv_slab7 +
      elf_alv_slab8 + epithelium_alv_slab8 + subepithelium_alv_slab8 +
      pmu1 * 0.1025204545 +
      pmu2 * 0.2050409091 +
      pmu3 * 0.2050409091 +
      pmu4 * 0.2050409091 +
      pmu5 * 0.2050409091 +
      pmu6 * 0.2050409091 +
      pmu7 * 0.2050409091 +
      pmu8 * 0.2050409091 +
      pmu9 * 0.2050409091 +
      pmu10 * 0.2050409091 +
      pmu11 * 0.2050409091 +
      pmu12 * 0.1025204545 +
      pml1 *     0.00525 +
      pml2 *      0.0105 +
      pml3 *      0.0105 +
      pml4 *      0.0105 +
      pml5 *      0.0105 +
      pml6 *      0.0105 +
      pml7 *      0.0105 +
      pml8 *     0.00525
    # Fig. 3 definition: all lung drug (undissolved particles, fluid,
    # epithelium and sub-epithelium) over the total lung volume.
    Clung <- alung / vlung
    # Plasma concentration (the paper's C_p). The paper reports no
    # residual-error model: this is a deterministic simulation model.
    Cc <- cven / bp
  })
}
