setwd('~/Documents/Harvey/research/Tyranni/tyranni_revision/sister_clade/neotropicalbirds/')
getwd()

source("R/packages.R")  # Load all the packages you need.
source("R/functions.R") # Load all the functions into your environment.
source("R/plan.R")      # Build your workflow plan data frame.

make(my_plan)
