import numpy as np
from scipy import stats

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

def statistical_test(x, y, test_type="wilcoxon"):
    # test stat, p-val
    if test_type == "t":
	stat, pval = np.ravel(stats.ttest_ind(x, y))
    elif test_type == "wilcoxon":
	stat, pval = np.ravel(stats.wilcoxon(x, y))
    elif test_type == "mannwhitney":
	stat, pval = np.ravel(stats.mannwhitneyu(x, y))
    elif test_type == "kolmogorovsmirnov":
	stat, pval = np.ravel(stats.ks_2samp(x,y))
    else:
	raise ValueError("Invalid correlation type!")
    # To convert p-value to one-way need to apply the following
    # taking into account the direction of the comparison
    # if stat >= 0:
    # 	pval = pval[1] / 2
    # else:
    #	pval = 1 - pval[1] / 2
    return stat, pval

def density_estimation(self, occurences, possible_values):
    kde = stats.gaussian_kde(map(float, occurences))
    p = kde(possible_values)
    return p / sum(p)

if __name__ == "__main__":
    main()

