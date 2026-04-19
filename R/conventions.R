#' nlmixr2lib convention standards
#'
#' Internal register of canonical parameter, compartment, covariate, and
#' residual-error names used by [checkModelConventions()]. Mirrors the
#' authoritative documentation in
#' `.claude/skills/extract-literature-model/references/` (especially
#' `naming-conventions.md` and `covariate-columns.md`) and
#' `vignettes/create-model-library.Rmd`. Keep this list in sync with those
#' sources; when a new canonical covariate is ratified in `covariate-columns.md`
#' add it here in the same commit.
#'
#' @keywords internal
#' @noRd
.nlmixr2libConventions <- list(
  pkParams = c(
    "lka", "lcl", "lvc", "lvp", "lvp2", "lq", "lq2", "lfdepot"
  ),
  pkBareParams = c(
    "ka", "cl", "vc", "vp", "vp2", "q", "q2", "kel",
    "k12", "k21", "k13", "k31", "fdepot"
  ),
  compartments = c(
    "depot", "central", "peripheral1", "peripheral2", "effect"
  ),
  compartmentRegex = "^(transit|effect)[0-9]+$",
  observationVar = "Cc",
  residualError = c("propSd", "addSd"),
  transformPrefixes = c("l", "logit", "probit"),
  covEffectPattern = "^e_[A-Za-z0-9]+_[A-Za-z0-9]+$",
  requiredUnits = c("time", "dosing", "concentration"),
  requiredMetadata = c("description", "reference", "units"),
  deprecatedResidualError = c(
    "prop.err", "add.err", "propErr", "addErr",
    "err.prop", "err.add"
  ),
  deprecatedIivPrefixes = c("iiv_", "IIV_", "bsv_", "BSV_"),
  canonicalCovariates = list(
    WT = list(units = "kg", type = "continuous", aliases = character()),
    AGE = list(units = "years", type = "continuous", aliases = character()),
    LBM = list(units = "kg", type = "continuous", aliases = character()),
    SEXF = list(units = "(binary)", type = "binary",
                aliases = c("SEXM", "SEX")),
    CHILD = list(units = "(binary)", type = "binary", aliases = character()),
    ADOLESCENT = list(units = "(binary)", type = "binary",
                      aliases = character()),
    PAGE = list(units = "months", type = "continuous", aliases = character()),
    PNA = list(units = "months", type = "continuous", aliases = character()),
    GA = list(units = "weeks", type = "continuous", aliases = character()),
    eGFR = list(units = "mL/min/1.73 m^2", type = "continuous",
                aliases = character()),
    CREAT = list(units = "umol/L or mg/dL", type = "continuous",
                 aliases = c("CRE", "SCR")),
    ALB = list(units = "g/dL or g/L", type = "continuous",
               aliases = character()),
    hsCRP = list(units = "mg/L", type = "continuous",
                 aliases = c("CRPHS", "HSCRP")),
    RACE_BLACK = list(units = "(binary)", type = "binary",
                      aliases = c("BLACK")),
    RACE_BLACK_OTH = list(units = "(binary)", type = "binary",
                          aliases = c("BLACK_OTH")),
    RACE_ASIAN = list(units = "(binary)", type = "binary",
                      aliases = c("ASIAN")),
    RACE_ASIAN_AMIND_MULTI = list(units = "(binary)", type = "binary",
                                  aliases = c("ASIAN_AMIND_MULTI")),
    RACE_MULTI = list(units = "(binary)", type = "binary",
                      aliases = c("MULTIRACIAL")),
    RACE_OTHER = list(units = "(binary)", type = "binary", aliases = character()),
    ADA_POS = list(units = "(binary)", type = "binary", aliases = c("ADA")),
    FED = list(units = "(binary)", type = "binary", aliases = character()),
    TABLET = list(units = "(binary)", type = "binary", aliases = character()),
    RIA_ASSAY = list(units = "(binary)", type = "binary",
                     aliases = character()),
    FORM_NS0 = list(units = "(binary)", type = "binary", aliases = character()),
    FORM_CHO_PHASE2 = list(units = "(binary)", type = "binary",
                           aliases = character()),
    dilution = list(units = "(binary)", type = "binary",
                    aliases = c("DILUTION")),
    nonECZTRA = list(units = "(binary)", type = "binary",
                     aliases = c("NON_ECZTRA", "STUDY_NON_ECZTRA")),
    SEASON2 = list(units = "(binary)", type = "binary", aliases = character()),
    STUDY1 = list(units = "(binary)", type = "binary", aliases = character()),
    STUDY5 = list(units = "(binary)", type = "binary", aliases = character()),
    ooc1 = list(units = "(binary)", type = "binary", aliases = character()),
    ooc2 = list(units = "(binary)", type = "binary", aliases = character()),
    ooc3 = list(units = "(binary)", type = "binary", aliases = character()),
    ooc4 = list(units = "(binary)", type = "binary", aliases = character()),
    JAPANESE_HV = list(units = "(binary)", type = "binary",
                       aliases = character())
  )
)

.nlmixr2libCovariateAliasMap <- function() {
  conv <- .nlmixr2libConventions$canonicalCovariates
  out <- character()
  for (canon in names(conv)) {
    for (a in conv[[canon]]$aliases) {
      out[a] <- canon
    }
  }
  out
}
