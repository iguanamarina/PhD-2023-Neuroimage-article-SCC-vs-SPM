install.packages("unix") 
library(unix)

rlimit_as(1e12)  #increases to ~12GB

rlimit_all()
