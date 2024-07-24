for (m in c("PK_1cmt_des", "PK_2cmt_des", "PK_3cmt_des")) {
  for (v in c("in", "out")) {
    test_that(paste0("addIndirectLin: ", m, "; stim=", v), {
      expect_error(readModelDb(m) |> addIndirectLin(stim=v), NA)
    })
    test_that(paste0("addIndirectLin: ", m, "; inhib=", v), {
      expect_error(readModelDb(m) |> addIndirectLin(inhib=v), NA)
    })
  }
}
