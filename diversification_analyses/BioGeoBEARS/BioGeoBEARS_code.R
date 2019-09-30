setwd('/Users/michaelharvey/Documents/research/Tyranni/v2/div_analyses/biogeobears')
getwd()

library(optimx)
library(FD)
library(snow)
library(parallel)
library(BioGeoBEARS)
library(ape)
library(phytools)

library(BioGeoBEARS)
source("http://phylo.wdfiles.com/local--files/biogeobears/cladoRcpp.R") # (needed now that traits model added; source FIRST!)
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_add_fossils_randomly_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_basics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_calc_transition_matrices_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_classes_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_detection_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_DNA_cladogenesis_sim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_extract_Qmat_COOmat_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_generics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_models_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_on_multiple_trees_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_plots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_readwrite_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_simulate_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_makePlots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stochastic_mapping_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stratified_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_univ_model_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_uppass_probs_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_loglike_sp_v01.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/get_stratified_subbranch_top_downpass_likelihoods_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/runBSM_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/stochastic_map_given_inputs.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/summarize_BSM_tables_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_traits_v1.R") # added traits model
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)

# Process geography file and tree
#tree.data <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_HowardMoore.tre')
#geog.data <- read.table("~/Documents/research/Tyranni/v2/div_analyses/biogeobears/Range_data_AOSHM_wwf_DEC_noheader.txt", header=FALSE, sep="\t")
#tree.data
#nrow(geog.data)
#keep <- geog.data[geog.data$V1 %in% tree.data$tip.label,]
#nrow(keep)
#length(tree.data$tip.label)
#range.missing <- tree.data$tip.label[!(tree.data$tip.label %in% geog.data$V1)]
#tree.trimmed <- drop.tip(tree.data, range.missing)
#tree.trimmed
#write.tree(tree.trimmed, "DEC_input.tre")
#write.table(cbind(as.character(keep$V1), keep$V2), "DEC_input.txt")

tree <- np('~/Documents/research/Tyranni/v2/div_analyses/biogeobears/DEC_input.tre')
geog <- np("~/Documents/research/Tyranni/v2/div_analyses/biogeobears/DEC_input.txt")
max_range_size <- 3

#######################################################
# Run DEC
#######################################################

# Intitialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = tree

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geog

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1
# (use more cores to speed it up; this requires
# library(parallel) and/or library(snow). The package "parallel" 
# is now default on Macs in R 3.0+, but apparently still 
# has to be typed on some Windows machines. Note: apparently 
# parallel works on Mac command-line R, but not R.app.
# BioGeoBEARS checks for this and resets to 1
# core with R.app)

# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
# I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
# but the results are imprecise and so I haven't explored it further.
# In a Bayesian analysis, it might work OK, but the ML point estimates are
# not identical.
# Also, I have not implemented all functions to work with force_sparse=TRUE.
# Volunteers are welcome to work on it!!
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DEC model
# (nothing to do; defaults)

# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object

# This contains the model object
BioGeoBEARS_run_object$BioGeoBEARS_model_object

# This table contains the parameters of the model 
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

# Run this to check inputs. Read the error messages if you get them!
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = "DEC_Tyranni.Rdata"
if (runslow)
    {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res    

    save(res, file=resfn)
    resDEC = res
    } else {
    # Loads to "res"
    load(resfn)
    resDEC = res
    }

#######################################################
# Plot ancestral states - DEC
#######################################################

pdf("DEC_Tyranni.pdf", width=9, height=7)

tree.data <- read.tree(tree)
pie.info <- results_object$ML_marginal_prob_each_state_at_branch_top_AT_node[(Ntip(tree.data)+1):nrow(results_object$ML_marginal_prob_each_state_at_branch_top_AT_node),]
tipstates.data <- read.table("~/Documents/research/Tyranni/v2/div_analyses/biogeobears/DEC_input_noheader.txt", row.names=1)
tipstates <- tipstates.data[,1]
names(tipstates) <- rownames(tipstates.data)
tipstates <- tipstates[tree.data$tip.label]
colors <- c("gray", "blue", "green",  "red", "lightseagreen", "mediumpurple2", "orange", "gray")
tipstates[tipstates==100] <- colors[2]
tipstates[tipstates==010] <- colors[3]
tipstates[tipstates==001] <- colors[4]
tipstates[tipstates==110] <- colors[5]
tipstates[tipstates==101] <- colors[3]
tipstates[tipstates==011] <- colors[6]
plot(tree.data, type="fan", tip.color=tipstates, cex=0.1)
par(lwd=0.01)
nodelabels(pie=pie.info, piecol=colors, border=NULL, cex=0.15)
ranges_list = areas_list_to_states_list_new(areas=areas, maxareas=length(areas),
include_null_range=TRUE, split_ABC=FALSE)
ranges = unlist(ranges_list)
ranges
legend("topleft", legend = c("Neotropical", "Nearctic", "Old World", "Neotropical/Nearctic", "Neotropical/Old World", "Nearctic/Old World"), col=colors[2:8], pch=16, bty="n", cex=0.8)

dev.off()


#######################################################
# Stochastic mapping on DEC M3b stratified with islands coming up
#######################################################
clado_events_tables = NULL
ana_events_tables = NULL
lnum = 0

#######################################################
# Get the inputs for Biogeographical Stochastic Mapping
# Note: this can be slow for large state spaces and trees, since 
# the independent likelihoods for each branch are being pre-calculated
# E.g., for 10 areas, this requires calculation of a 1024x1024 matrix
# for each branch.  On a tree with ~800 tips and thus ~1600 branches, this was about 1.6 gigs
# for storage of "BSM_inputs_file.Rdata".
# Update: 2015-09-23 -- now, if you used multicore functionality for the ML analysis,
# the same settings will be used for get_inputs_for_stochastic_mapping().
#######################################################
BSM_inputs_fn = "BSM_inputs_file.Rdata"
runInputsSlow = TRUE
if (runInputsSlow)
    {
    stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=resDEC)
    save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)
    } else {
    # Loads to "stochastic_mapping_inputs_list"
    load(BSM_inputs_fn)
    } # END if (runInputsSlow)

# Check inputs (doesn't work the same on unconstr)
names(stochastic_mapping_inputs_list)
stochastic_mapping_inputs_list$phy2
stochastic_mapping_inputs_list$COO_weights_columnar
stochastic_mapping_inputs_list$unconstr
set.seed(seed=as.numeric(Sys.time()))

runBSMslow = TRUE
if (runBSMslow == TRUE)
    {
    # Saves to: RES_clado_events_tables.Rdata
    # Saves to: RES_ana_events_tables.Rdata
    BSM_output = runBSM(resDEC, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxnum_maps_to_try=200, nummaps_goal=100, maxtries_per_branch=40000, save_after_every_try=TRUE, savedir=getwd(), seedval=12345, wait_before_save=0.01)

    RES_clado_events_tables = BSM_output$RES_clado_events_tables
    RES_ana_events_tables = BSM_output$RES_ana_events_tables
    } else {
    # Load previously saved...

    # Loads to: RES_clado_events_tables
    load(file="RES_clado_events_tables.Rdata")
    # Loads to: RES_ana_events_tables
    load(file="RES_ana_events_tables.Rdata")
    BSM_output = NULL
    BSM_output$RES_clado_events_tables = RES_clado_events_tables
    BSM_output$RES_ana_events_tables = RES_ana_events_tables
    } # END if (runBSMslow == TRUE)

# Extract BSM output
clado_events_tables = BSM_output$RES_clado_events_tables
ana_events_tables = BSM_output$RES_ana_events_tables
head(clado_events_tables[[1]])
head(ana_events_tables[[1]])
length(clado_events_tables)
length(ana_events_tables)

getranges_from_LagrangePHYLIP(lgdata_fn=BioGeoBEARS_run_object$geogfn)
include_null_range = TRUE
areanames = names(tipranges@df)
areas = areanames
max_range_size = 3

# Note: If you did something to change the states_list from the default given the number of areas, you would
# have to manually make that change here as well! (e.g., areas_allowed matrix, or manual reduction of the states_list)
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)

colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=max_range_size, plot_null_range=TRUE)

############################################
# Setup for painting a single stochastic map
############################################
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
stratified = FALSE
clado_events_table = clado_events_tables[[1]]
ana_events_table = ana_events_tables[[1]]

#######################################################
# Summarize stochastic map tables
#######################################################
length(clado_events_tables)
length(ana_events_tables)

head(clado_events_tables[[1]][,-20])
tail(clado_events_tables[[1]][,-20])

head(ana_events_tables[[1]])
tail(ana_events_tables[[1]])

areanames = names(tipranges@df)
actual_names = areanames
actual_names

# Get the dmat and times (if any)
dmat_times = get_dmat_times_from_res(res=resDEC, numstates=NULL)
dmat_times

# Extract BSM output
clado_events_tables = BSM_output$RES_clado_events_tables
ana_events_tables = BSM_output$RES_ana_events_tables

# Simulate the source areas
BSMs_w_sourceAreas = simulate_source_areas_ana_clado(resDEC, clado_events_tables, ana_events_tables, areanames)
clado_events_tables = BSMs_w_sourceAreas$clado_events_tables
ana_events_tables = BSMs_w_sourceAreas$ana_events_tables
# Count all anagenetic and cladogenetic events
counts_list = count_ana_clado_events(clado_events_tables, ana_events_tables, areanames, actual_names)

Nt2Na <- counts_list$all_dispersals_counts_cube[1,2,1:50]
hist(Nt2Na)
min(Nt2Na)
max(Nt2Na)
mean(Nt2Na)
sd(Nt2Na)
Na2Nt <- counts_list$all_dispersals_counts_cube[2,1,1:50]
min(Na2Nt)
max(Na2Nt)
mean(Na2Nt)
sd(Na2Nt)
hist(Na2Nt)

# Age of oldest event
H <- nodeHeights(tree.data)
nn <- findMRCA(tree.data, c("Xenotr_mexic_LSUMZ24529", "Empinax_oberh_L48075"))
max(H) - H[tree.data$edge==nn][1]

# Plot with nodes colored according to ancestral states