#parameter file 

#for file full_path='/orange/eiken/Tcell_motility/Tcells_TRITC_20X_0p86px_1p45µm_5min_2021FEB09.nd2'

#get peak z
flux_lim=1000.#1221. #originally from get_spot_flux_cutoffs.ipynb, lowered by factor of 3
#spots with fluxes below this limit will not be analyzed

sig_lim=3.
n_points=3#round(z_levels[-1]*0.05) #### spots must have at least 5 % of total # of points
r1=10. #source radius
r2=15. #inner annulus radius
r3=20. #outer annulus radius

#track_spots_new
# flux_lim=400.#1221. #spots with fluxes below this limit will not be analyzed
d_lim=25. #spots with distances greater than this limit will not be counted as a potential match/ 
#brightest spot within this limit will be considered the match
num_times=4 #if a spot appears less than this limit/# of time slices, it will not be written to a file
num_empty=4 #if a spot doesn't appear in the next 4 time slices, stop looking for it
flux_drop_lim1=20.
flux_drop_lim2=230.
factor=1.5