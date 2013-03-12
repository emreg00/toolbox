
def main():
    return

def get_subsamples(scores, n_fold=10000, n_sample=1000):
    from random import shuffle
    for i in xrange(n_fold):
	#if with_replacement:
	#	size = len(scores)-1
	#	selected = empty(n_sample)
	#	for i in xrange(n_sample): 
	#	    selected[i] = scores[randint(0,size)]
	shuffle(scores)
	selected = scores[:n_sample]
	yield selected
    return

def k_fold_cross_validation(X, K, randomize = False, replicable = None):
    """
    By John Reid (code.activestate.com)
    Generates K (training, validation) pairs from the items in X.

    Each pair is a partition of X, where validation is an iterable
    of length len(X)/K. So each training iterable is of length (K-1)*len(X)/K.

    If randomise is true, a copy of X is shuffled before partitioning,
    otherwise its order is preserved in training and validation.

    If replicable is not None, this number is used to create the same random splits at each call
    """
    #if randomize: from random import shuffle; X=list(X); shuffle(X)
    if randomize: 
	from random import shuffle, seed
	X=list(X)
	if replicable is not None: 
	    seed(replicable)
	shuffle(X)
    for k in xrange(K):
	training = [x for i, x in enumerate(X) if i % K != k]
	validation = [x for i, x in enumerate(X) if i % K == k]
	yield k+1, training, validation
    return

def generate_samples_from_list_without_replacement(elements, sample_size, n_folds = None, replicable = None):
    """
	Iteratively returns (yields) n_folds sublists of elements with a size of sample_size 
	n_folds: If None calculated to cover as much elements as possible
	replicable: If not None uses this replicable as the seed for random
    """
    from random import shuffle, seed
    if replicable is not None:
	seed(replicable)
    shuffle(elements)
    if n_folds is None:
	from math import ceil
	#n_folds = len(elements) / sample_size
	n_folds = int(ceil(float(len(elements)) / sample_size))
    for i in range(n_folds):
	if (i+1)*sample_size < len(elements):
	    yield elements[i*sample_size:(i+1)*sample_size]
	else:
	    yield elements[i*sample_size:]
    return


if __name__ == "__main__":
    main()

