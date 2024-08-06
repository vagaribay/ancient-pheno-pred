import fastDTWF, torch, numpy as np

# population size: diploid/haploid
N, N_haploid = 10000, 20000

# region size
L = 20*25000

# per-generation mutation rate
mu = 1e-8
mu_torch = torch.tensor(mu, dtype=torch.float64)

# no strength of selection against the "1" allele in heterozygotes/homozygotes
s_het = torch.tensor(0, dtype=torch.float64)

# SFS assuming neutrality: xi = thetha * 1/i, x = (x1,...,xn-1)
# theta = 2Nμ; assuming haploids and a locus size of 500,000 bp (20 regions of 25000 bp)
# theta = 4NμL
#theta = 4*N*mu*L # 200

# sfs
#sfs = np.zeros(N+1)
#for i in range(1, N):
#    sfs[i] = round((theta)*(1 / i), 0)

#print(sfs)
#exit()

freqs = np.round(np.arange(0, 0.021, 0.001),3)
probs = np.arange(0, N_haploid+1) / N_haploid

out = open("/Users/valeriagby/desktop/ancient-pheno-prediction/scripts/fastDTWF/DTWF_seq_0_02_OUT.txt", "w+")

for frq in freqs:
    print(frq)
    idx = np.where(probs == frq)[0][0]
    # compute transition mass functions (TMF) --- the probability of transitioning from some number of alleles at one time point to a different number at some point in the future
    # we assume a population of 20000 haploids, and an initial frequency of freq
    transition_distribution = torch.zeros(N_haploid + 1, dtype=torch.float64)
    transition_distribution[idx] = 1.
    # compute the expected allele frequencies for each given allele frequency
    expected_allele_freqs = fastDTWF.wright_fisher_ps_mutate_first(pop_size=N_haploid, mu_0_to_1=mu_torch, mu_1_to_0=mu_torch, s_het=s_het)
    # "coarse grain" expected allele frequencies
    index_sets = fastDTWF.coarse_grain(p_list=expected_allele_freqs, pop_size=N_haploid, tv_sd=0.1, no_fix=False, sfs=False)
    
    # evolve the probabilities over g = 400 generations
    for g in range(401):
        transition_distribution = fastDTWF.mat_multiply(vec=transition_distribution, index_sets=index_sets, p_list=expected_allele_freqs, pop_size=N_haploid, row_eps=1e-8, no_fix=False, sfs=False)
        if g == 100 or g == 200 or g == 300 or g == 400:
            probs_T = transition_distribution.numpy()
            line = '\t'.join(map(str, probs_T))
            out.write(str(frq) + "\t" + str(g) + "\t" + line + '\n')

out.close()
