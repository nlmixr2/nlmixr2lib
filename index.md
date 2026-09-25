# nlmixr2lib

> ### Which version are you on?
>
> **This site documents the development version.** The articles here
> cover every model in the development library; the CRAN release
> contains far fewer, so a
> [`readModelDb()`](https://nlmixr2.github.io/nlmixr2lib/reference/readModelDb.md)
> call copied from an article will fail with `'name' not in database` if
> you are on CRAN.
>
> |                                    | models |
> |------------------------------------|--------|
> | CRAN `0.3.2` (released 2026-01-18) | 59     |
> | development `0.3.2.9000`           | 2912   |
>
> To get the models these articles describe:
>
> ``` r
>
> install.packages(
>   "nlmixr2lib",
>   repos = c("https://nlmixr2.r-universe.dev", getOption("repos"))
> )
> ```
>
> **Then restart R.**
> [`readModelDb()`](https://nlmixr2.github.io/nlmixr2lib/reference/readModelDb.md)
> and `modeldb` both resolve against the package as loaded, so a
> reinstall in a running session leaves the old model list in memory and
> the new models will still appear missing.
>
> Check it worked:
>
> ``` r
>
> packageVersion("nlmixr2lib")   # 0.3.2.9000
> nrow(nlmixr2lib::modeldb)      # 2912
> ```
>
> The r-universe repository is needed rather than plain
> `install_github()`: the development version requires
> `rxode2 (>= 5.1.8)`, which is not on CRAN yet, and r-universe serves
> prebuilt binaries so no compiler is required.

This is a model library for `nlmixr2`. The package allows a few ways to
interact with the model library:

``` r

# See all available models
modellib()
# Load the "PK_1cmt" model
modellib(name="PK_1cmt")
# Switch residual error to additive
modellib(name="PK_1cmt", reserr = "addSd")
# Add inter-individual variability on ka and v and switch residual error to
# additive and proportional
modellib(name="PK_1cmt", eta = c("lka", "lv"), reserr = c("addSd", "propSd"))
```

# Modifying models by piping

You may also modify any model from the library (or your own models) with
a piping interface. The code below adds inter-individual variability on
ka and v and then switches residual error to additive and proportional.

``` r
modellib(name="PK_1cmt") |>
  addEta(c("lka", "lv") |>
  addResErr(c("addSd", "propSd"))
```

# Possible extensions

The `modellib` function is set-up in way that it can be easily extended
and used in other applications. A possible extension could be
implementation in a shiny app. An app can be created to easily add new
models to the model library database (curated?), and directly make these
models available for other users. I believe there can be added value in
having a base model library that can be easily extended by the community
this way.
