<!-- badges: start -->
[![R-CMD-check](https://github.com/nlmixr2/nlmixr2lib/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/nlmixr2/nlmixr2lib/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/nlmixr2/nlmixr2lib/branch/main/graph/badge.svg)](https://app.codecov.io/gh/nlmixr2/nlmixr2lib?branch=main)
<!-- badges: end -->

# nlmixr2lib

This is a model library for `nlmixr2`.  The package allows a few ways to interact with the model library:


```r
# See all available models
modellib()
# Load the "PK_1cmt" model
modellib(name="PK_1cmt")
# Switch residual error to additive
modellib(name="PK_1cmt", reserr = "add")
# Add inter-individual variability on ka and v and switch residual error to
# additive and proportional
modellib(name="PK_1cmt", eta = c("lka", "lv"), reserr = "add+prop")
```

# Modifying models by piping

You may also modify any model from the library (or your own models) with a
piping interface.  The code below adds inter-individual variability on ka and v
and then switches residual error to additive and proportional.

```r
modellib(name="PK_1cmt") %>%
  addEta(c("lka", "lv") %>%
  addResErr("add+prop")
```

# Possible extensions

The `modellib` function is set-up in way that it can be easily
extended and used in other applications.  A possible extension could
be implementation in a shiny app. An app can be created to easily add
new models to the model library database (curated?), and directly make
these models available for other users.  I believe there can be added
value in having a base model library that can be easily extended by
the community this way.
