# ============================================================== #
#                 Main script to run packages                    #
# ============================================================== #

rm(list = ls())

# options(warn = 0)
options(error=recover)
options(show.error.locations=TRUE)

# Set different seeds
seed <- floor(runif(1)*1e6)
seed <- 12345
set.seed(seed)

