# try a grid of params

h2_array='0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8'

# mode 1: infinitesimal
mode=infinitesimal
for h2 in $h2_array
do
  snakemake -s simulator.snmk --configfile config.yaml -p all_simulate --config beta_dist=$mode pve=$h2
done

# mode 2: spike and slab
modepre=spike_n_slab_
nonzerofrac_array='0.01 0.05 0.1 0.2 0.3 0.4 0.5'

for h2 in $h2_array
do
  for nonzerofrac in $nonzerofrac_array
  do
    snakemake -s simulator.snmk --configfile config.yaml -p all_simulate --config beta_dist=$mode$nonzerofrac pve=$h2
  done
done