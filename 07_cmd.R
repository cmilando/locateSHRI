load("test_objects.Rdata")

# obviously takes slower
system("R CMD SHLIB simann_constrain.f90")
#dyn.unload("simann_constrain.so")
dyn.load("simann_constrain.so")

oo <- .Fortran("simann_constrain",
               S = as.integer(S),
               magic_n = as.integer(length(which(S == 1))),
               np = as.integer(nrow(site_pairs_sub)),
               nsites = as.integer(length(unique(site_pairs_sub$site_id))),
               site_pairs_sub = aa,
               row_lookup = bb,
               penalty = 1,
               SCORE = 0.,
               cooling_rate = 0.95,
               cc = cc,
               nzones = nzones,
               verbose = as.integer(1)
) 
which(S == 1)
which(oo$S == 1)
