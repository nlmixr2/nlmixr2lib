# Load the xie_2019_agomelatine model from the paper supplement

# The code here was used to initially prepare the model. Additional
# simplifications and updates were manually performed following the code here
# for loading.

library(readtext)
library(nonmem2rx)

lines_raw <- readtext("xie_2019_agomelatine_BCP-85-1003-s001.docx")
lines_prep <- strsplit(lines_raw$text, split = "\n")[[1]]
lines <- lines_prep[-length(lines_prep)]

# Modify the file to be a valid NONMEM control stream
lines[1] <- paste("$PROB", lines[1])
input_line <- which(startsWith(lines, "$INPUT"))
data_line <- which(startsWith(lines, "$DATA"))

lines[input_line] <- gsub(x = lines[input_line], pattern = "PERIOD.*$", replacement = "PERIOD")

file_prep <- writeLines(con = "xie_2019_agomelatine_BCP-85-1003-s001.ctl", text = lines)
nonmem_model <- nonmem2rx("xie_2019_agomelatine_BCP-85-1003-s001.ctl")
