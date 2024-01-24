# This script is for converting the elementary
# effects to .csv values

fname = "2024-01-22_MorrisAnalysis_r100_50days_r-100.RData"
rval = "100"
load(fname)

date2save <- Sys.Date()
notes <- readline(prompt = "notes for filename: ")

# Kplas
save_fname = paste0("./results/",
                    date2save,
                    "_MA_ee",
                    "_var-", "Kplas",
                    "_r-", rval,
                    "_notes-", notes,
                    ".csv")
write.csv(x_Kplas$ee, file = save_fname)

# Kmuscle
save_fname = paste0("./results/",
                    date2save,
                    "_MA_ee",
                    "_var-", "Kmusc",
                    "_r-", rval,
                    "_notes-", notes,
                    ".csv")
write.csv(x_Kmuscle$ee, file = save_fname)