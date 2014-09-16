import sys

poly = open(sys.argv[1])
div = open(sys.argv[2])

theta_avg_file = open(sys.argv[3], "w")
theta_mu_file = open(sys.argv[4], "w")

divergence = {}
#key:values
#gene_name : (num_syn diverged, rep diverged)
for line in div:
    sline = line.split()
    name = sline[0]
    num_syn = float(sline[1])
    num_div = float(sline[2])
    syn_frac = float(sline[3])
    div_frac = float(sline[4])
    
    divergence[name] = (syn_frac, div_frac)
    
print("theta_syn\tdivergence_syn\tratio_syn\ttheta_tot\tdivergence_tot\tratio_tot")
    
poly.readline()#get rid of header
#name samp_size num_syn_sites theta_syn num_poly_syn .... [16] num_rep_sites [17] theta_rep [18] num_poly_rep
lengths = []
polys = 0
for line in poly:
    sline = line.split()
    name = sline[0]
    syn_theta = float(sline[3])*float(sline[2])
    syn_theta_site = float(sline[3])
    rep_theta = float(sline[17])*float(sline[18])
    rep_theta_site = float(sline[17])
    
    #throw out loci with fewer than 50 syn sites
    #TODO add filter for total sites instead?
    if float(sline[2]) < 50:
        continue
    lengths.append(float(sline[2]))
    polys+= float(sline[4])
    
    theta_mu_file.write("%s\t%s\n" % (syn_theta, float(sline[2])))
    
    tot_theta = syn_theta+rep_theta
    tot_theta_site = tot_theta/(float(sline[2])+float(sline[18]))
    
    if name in divergence.keys():
        if divergence[name][0] == 0 or divergence[name][1] == 0:
            #no divergence, problem!
            sys.stderr.write("Divergence == 0 for %s\n" % name)
            continue
        print("%s\t%s\t%s\t%s\t%s\t%s" % (syn_theta_site, divergence[name][0], float(syn_theta_site)/(divergence[name][0]), tot_theta_site, divergence[name][1], float(tot_theta_site)/(divergence[name][1])))
    else:
        sys.stderr.write("Missing divergence info for gene - %s\n" % name)

N = 26

harmonic = 0

for i in range(1,N):
    harmonic += 1.0/i

sys.stderr.write("Harmonic mean: %.4f\n" % harmonic)
        
avg_theta = polys/harmonic/sum(lengths)
sys.stderr.write("Average theta: %.4f\n" % avg_theta)

for locus in lengths:
    t = locus*avg_theta
    theta_avg_file.write("%s\t%s\n" % (t, locus))
    
theta_avg_file.close()
theta_mu_file.close()