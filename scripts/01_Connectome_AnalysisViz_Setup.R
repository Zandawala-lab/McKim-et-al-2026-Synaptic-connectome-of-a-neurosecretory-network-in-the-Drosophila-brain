# -- General setup steps for running connectome analyses in R ------------------
# ------------------------------------------------------------------------------

# Assumes starting from scratch

# Install the main package needed ----------------------------------------------
# Natverse: https://natverse.org/install/
# Follow instructions on that page

# Run full install version
install.packages("natmanager")
natmanager::install('natverse') 

# Other packages to install
install.packages('rglplus')
install.packages('ggforce')
install.packages('readr')

install.packages('reticulate') # to use python in R; required for caveclient connectivity data
# if issues, see: https://natverse.org/fafbseg/articles/installing-cloudvolume-meshparty.html
simple_python("full")
# if still issues with using python and reticulate, follow steps here:
# https://natverse.org/fafbseg/reference/simple_python.html#details
# specify which Python you want to use with the RETICULATE_PYTHON environment variable
# set RETICULATE_PYTHON with usethis::edit_r_environ().
# RETICULATE_PYTHON = add_path_info_here/r-miniconda-arm64/envs/r-reticulate/bin/python

# Need to run it again to get cloudvolume properly installed (and potentially restart R)
simple_python("full")

# set flywire token for access and version for caveclient: https://natverse.org/fafbseg/
# How to check access: 
# https://natverse.org/fafbseg/reference/flywire_set_token.html
flywire_set_token()

# other packages needed
install.packages('neuprintr')             # hemibrain/maleCNS specific
natmanager::install(pkgs = "fafbseg")     # check setup with dr_fafbseg()
natmanager::install(pkgs = 'coconatfly')  # used for clustering analyses

remotes::install_github('natverse/bancr') # BANC
# when calling library('bancr') was warned to install:
download_jefferislab_registrations()

# other things to do/set:
# ensure neuprint account created and grab token for setup:
# https://natverse.org/neuprintr/#authentication

# banc token: https://github.com/natverse/bancr#installation
banc_set_token()




