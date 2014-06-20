import numpy as np
from scipy import stats
from numpy import median

def main():
    import sys
    sc = float(sys.argv[1])
    #alist=[1,1,0]
    #alist = [ float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]) ]
    alist = map(float, sys.argv[2:])
    #print alist
    m, s = calc_mean_and_sigma(alist)
    if s == 0:
        print "0 variation"
    else:
        print "(%.2f - %.2f) / %.2f = %.2f" % ( sc, m, s, (sc - m) / s)
    return

def convert_p_values_to_z_scores(pvalues):
    a = np.random.normal(size=1000000) #1000000 # 10000000
    return map(lambda x: stats.scoreatpercentile(a, 100-(100*x/2.0)), pvalues)

def correct_pvalues_for_multiple_testing(pvalues, correction_type = "Benjamini-Hochberg"):
    """
    consistent with R - print correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1]) 
    """
    from numpy import array, empty
    pvalues = array(pvalues)
    n = pvalues.shape[0]
    new_pvalues = empty(n)
    n = float(n)
    if correction_type == "Bonferroni":
	new_pvalues = n * pvalues
    elif correction_type == "Bonferroni-Holm":
	values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
	values.sort()
	for rank, vals in enumerate(values):
	    pvalue, i = vals
	    new_pvalues[i] = (n-rank) * pvalue
    elif correction_type == "Benjamini-Hochberg":
	values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
	values.sort()
	values.reverse()
	new_values = []
	for i, vals in enumerate(values):
	    rank = n - i
	    pvalue, index = vals
	    new_values.append((n/rank) * pvalue)
	for i in xrange(0, int(n)-1): 
	    if new_values[i] < new_values[i+1]:
		new_values[i+1] = new_values[i]
	for i, vals in enumerate(values):
	    pvalue, index = vals
	    new_pvalues[index] = new_values[i]
	#for rank, vals in enumerate(values):
	    #pvalue, i = vals
	    #new_pvalues[i] = (n/(rank+1)) * pvalue
    else:
	raise ValueError("Unknown correction type: " + correction_type)
    return new_pvalues

def calc_mean_and_sigma(alist):
    return mean(alist), sigma(alist)

def mean(x):
    return np.mean(x)

def sigma(x):
    return np.std(x)

def correlation(x, y, cor_type="pearson"):
    # coef, p-val
    if cor_type == "pearson":
	coef, pval = np.ravel(stats.pearsonr(x, y))
    elif cor_type == "spearman":
	coef, pval = np.ravel(stats.spearmanr(x, y))[0]
    else:
	raise ValueError("Invalid correlation type!")
    return coef, pval

def statistical_test(x, y, test_type="wilcoxon", alternative="two-sided"):
    # test stat, p-val
    if test_type == "t":
	stat, pval = np.ravel(stats.ttest_ind(x, y, equal_var=False))
    elif test_type == "wilcoxon": # Requires equal size
	stat, pval = np.ravel(stats.wilcoxon(x, y))
    elif test_type == "mannwhitney": # returns one-sided by default
	stat, pval = np.ravel(stats.mannwhitneyu(x, y))
    elif test_type == "ks":
	stat, pval = np.ravel(stats.ks_2samp(x,y))
    else:
	raise ValueError("Invalid correlation type!")
    #return stat, pval
    # To convert p-value to one-way, it is inconsistent with R though
    if test_type == "wilcoxon":
	stat2 = median(x) - median(y)
	if stat2 >= 0:
	    if alternative == "greater":
		pval = pval / 2
	    elif alternative == "less":
		pval = 1 - pval / 2
	else:
	    if alternative == "greater":
		pval = 1 - pval / 2
	    elif alternative == "less":
		pval = pval / 2
    elif test_type == "mannwhitney":
	stat2 = median(x) - median(y)
	if alternative == "two-sided":
	    pval = (2 * pval)
	elif alternative == "less":
	    if stat2 >= 0: 
		pval = 1 - pval
	elif alternative == "greater":
	    if stat2 < 0:
		pval = 1 - pval
    elif alternative != "two-sided":
	raise ValueError("Not implemented!")
    return stat, pval

def hypergeometric_test(picked_good, picked_all, all_all, all_good):
    k = len(picked_good) 
    M = len(all_all) 
    n = len(all_good) 
    N = len(picked_all) 
    val = sum(stats.hypergeom.pmf(range(k,min(N,n)+1), M, n, N))
    return val

def density_estimation(occurences, possible_values):
    kde = stats.gaussian_kde(map(float, occurences))
    p = kde(possible_values)
    return p / sum(p)

def fisher_exact(tp, fp, fn, tn, alternative="two-sided"):
    oddsratio, pvalue = stats.fisher_exact([[tp, fp], [fn, tn]], alternative)
    return oddsratio, pvalue

if __name__ == "__main__":
    main()

