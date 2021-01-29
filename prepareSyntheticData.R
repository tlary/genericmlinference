library(data.table)

# read the generated data
y <- fread("y.csv")[2:20001]
X <- fread("X.csv")[2:20001]
treatment <- fread("treatment.csv")[2:20001]
assignment <- fread("assignment.csv")[2:20001]

# rename columns
for (dt in c("y", "X", "treatment", "assignment")) {
  setnames(get(dt), old="V1", new="id")
}
setnames(y, old="V2", new="outcome")
setnames(assignment, old="V2", new="treatment")

# combine to single data.table
dt <- merge(y, X, by="id")
dt <- merge(dt, assignment)

# save as .csv file
fwrite(dt, "syntheticData.csv")
