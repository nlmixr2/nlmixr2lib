test_that("convertQuad", {

   r1 <-readModelDb("PK_2cmt_no_depot") |>
     addIndirectLin(stim="out") |>
     convertQuad()

   expect_equal(rxode2::modelExtract(r1, "d/dt(R)"),
                "d/dt(R) <- kin - kout * R * (1 + (Ek * Cc + Ek2 * Cc^2))")

    r1 <-readModelDb("PK_2cmt_no_depot") |>
      addDirectLin() |>
      convertQuad()

    expect_equal(rxode2::modelExtract(r1, "effect"),
                "effect <- Ek * Cc + Ek2 * Cc^2")

    r1 <- readModelDb("PK_2cmt_no_depot") |>
     addEffectCmtLin() |>
     convertQuad()

    expect_equal(rxode2::modelExtract(r1, "effect"),
                 "effect <- Ek * Ce + Ek2 * Ce^2")

})
