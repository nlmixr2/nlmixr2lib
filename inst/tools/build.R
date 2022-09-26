modeldb <- data.frame()

addDirToModelDb(dir = "inst")

if (!dir.exists("data")) {
  dir.create("data")
}

.md5 <- digest::digest(modeldb)
print(modeldb)
.len <- length(modeldb$name)
if (file.exists("data/modeldb.rda")) unlink("data/modeldb.rda")
save(modeldb, file="data/modeldb.rda", compress="bzip2", version=2, ascii=FALSE)

modeldb.Rd <- suppressWarnings(readLines(con = "man/modeldb.Rd"))
modeldb.Rd <- gsub("XXX", .len, modeldb.Rd)

file.out <- file("man/modeldb.Rd", "wb")
writeLines(modeldb.Rd, file.out)
close(file.out)

file.out <- file("src/dbver.h", "wb")
writeLines(c(paste0("#define __MD5__ \"", .md5, "\""),
             ""), file.out)
close(file.out)


.in <- suppressWarnings(readLines("src/Makevars.in"))
if (.Platform$OS.type == "windows" && !file.exists("src/Makevars.win")) {
  file.out <- file("src/Makevars.win", "wb")
  writeLines(.in, file.out)
  close(file.out)
} else {
  file.out <- file("src/Makevars", "wb")
  writeLines(.in, file.out)
  close(file.out)
}
