
import network_utilities, text_utilities, wrappers
from time import strftime, sleep
import subprocess
import os, numpy, cPickle
import hashlib


def run_proximity_on_cluster(parameters, source_to_geneids, target_to_geneids, run_mode='array job'):
    """
    run_mode: array job | single job | run local | run cluster
    """
    network_file = parameters.get("network_file")
    n_random = int(parameters.get("n_random"))
    min_bin_size = int(parameters.get("n_node_in_bins"))
    seed = int(parameters.get("random_seed"))
    executable_path = parameters.get("executable_path")
    try:
	cluster_dir = parameters.get("cluster_dir")
	network_file = network_file.replace(parameters.get("base_dir"), cluster_dir) 
    except:
	cluster_dir = None
    qname = "all.q"
    delay = 10
    values = []
    source_to_md5 = {}
    md5_to_sources = {}
    n_start, n_end = 0, 700000 #15000 638952 
    if run_mode == "run cluster":
	i = n_start
	while i < n_end+1:
	    experiment_count = get_number_of_jobs_in_queues()
	    while experiment_count > 60: 
		sleep(delay)
		experiment_count = get_number_of_jobs_in_queues()
	    input_file = parameters.get("data_dir") + "/input/"
	    if not os.path.exists(input_file):
		continue
	    score_command = "-p %s/input/ -i %d -j %d" % (parameters.get("data_dir"), i, i + 1000)
	    os.system("sbatch -x node30 run_proximity.sh %s" % score_command)  
	    i += 1000
	return
    for source, geneids_source in source_to_geneids.iteritems(): 
	#print source, len(geneids_source)
	#! source = text_utilities.convert_to_R_string(source)
	md5 = hashlib.md5("-".join(sorted(geneids_source))).hexdigest()
	source_to_md5[source] = (md5, geneids_source)
	if md5 in md5_to_sources:
	    md5_to_sources.setdefault(md5, set()).add(source)
	    continue
	else:
	    md5_to_sources.setdefault(md5, set()).add(source)
	source = md5
	for target, geneids_target in target_to_geneids.iteritems():
	    #print target, len(geneids_target)
	    #! target = text_utilities.convert_to_R_string(target)
	    out_file = parameters.get("output_dir") + "/%s_%s.out" % (source, target) 
	    if cluster_dir is not None:
		out_file = out_file.replace(parameters.get("base_dir"), cluster_dir) 

	    if os.path.exists(out_file):
	    	continue
	    score_command = ' -x %d -m %d -n %d -e "%s" -o "%s" -s "%s" -t "%s"' % (seed, min_bin_size, n_random, network_file, out_file, ",".join(geneids_source), ",".join(geneids_target))
	    if run_mode != "array job" and run_mode != "run cluster":
		score_command = executable_path + score_command
	    if run_mode == "array job":
		print "%s" % (score_command.replace('"', ''))
		values.append(score_command.replace('"', ''))
	    elif run_mode == "single job":
		print "qsub -cwd -S /bin/bash -o out -e err -v PATH=$PATH -v PYTHONPATH=$PYTHONPATH -q %s -N %s_%s -b y %s" % (qname, source[:3], target[:3], score_command)
	    elif run_mode == "run local":
		print "%s" % score_command
		os.system(score_command)
	    elif run_mode == "run cluster":
		experiment_count = get_number_of_jobs_in_queues()
		while experiment_count > 60: 
		    sleep(delay)
		    experiment_count = get_number_of_jobs_in_queues()
		#print score_command
		#os.system("qsub -cwd -S /bin/bash -o out -e err -v PATH=$PATH -v PYTHONPATH=$PYTHONPATH -q %s -N %s_%s -b y %s" % (qname, source[:3], target[:3], score_command))
		#os.system("sbatch -x node30 run_proximity.sh -f ../data/input/%i.txt" % i)
		os.system("sbatch -x node30 run_proximity.sh %s" % score_command)  
	    else:
		raise ValueError("Unknown run_mode: %s" % run_mode)
    n = 0
    for md5, sources in md5_to_sources.iteritems():
	if len(sources) > 1:
	    n += len(sources) - 1
	    print len(sources)
	    for source in sources:
		val, targets = source_to_md5[source]
		for source2 in sources:
		    val, targets2 = source_to_md5[source2]
		    if targets != targets2:
			print targets, targets2
    print len(source_to_geneids), n, len(md5_to_sources)
    return values


def run_guild_on_cluster(parameters, target_to_geneids, run_mode='array job', method = 's'):
    """
    run_mode: array job | run local 
    method: netshort 's' | netrank 'r'
    """
    network_lcc_file = parameters.get("network_file")
    executable_path = parameters.get("guild_path")
    output_dir = parameters.get("output_dir") + "/"
    qname = "all.q"
    delay = 10
    network = wrappers.get_network(parameters.get("network_file"), only_lcc = True) # already using LCC file
    nodes = network.nodes()
    for target, geneids in target_to_geneids.iteritems():
	#print target, len(geneids_target)
	target = text_utilities.convert_to_R_string(target)
	target_to_score = dict((gene, 1.0) for gene in geneids)
	out_file = parameters.get("output_dir") + "/%s.n%s" % (target, method) 
	if os.path.exists(out_file):
	    continue
	if run_mode == "run local":
	    qName = None
	elif run_mode != "array job":
	    raise ValueError("Unknown run_mode: %s" % run_mode)
	score_command = wrappers.run_guild(target, target_to_score, nodes, network_lcc_file, output_dir, executable_path, background_score = 0.01, qname = qname, method = method) 
    return


def get_number_of_jobs_in_queues():
    #p1 = subprocess.Popen(["qstat -u eguney"], stdout=subprocess.PIPE)
    #p2 = subprocess.Popen(["wc", "-l"], stdin=p1.stdout, stdout=subprocess.PIPE)
    #experiment_count = int(p2.communicate()[0])
    #text = subprocess.check_output(["qstat", "-u", "eguney"])
    text = subprocess.check_output(["squeue", "-u", "emre"])
    experiment_count = len(text.split("\n")) - 1
    return experiment_count


def output_proximity_results(parameters, sources, targets, out_file, source_to_targets=None):
    f = open(out_file, 'w')
    if source_to_targets is not None:
	f.write("source\ttarget\tflag\tz\n")
    else:
	f.write("source\ttarget\tz\n")
    source_to_target_to_proximity = get_proximity_values(parameters, sources, targets)
    for source in sources: 
	for target in targets:
	    z = source_to_target_to_proximity[source][target]
	    if source_to_targets is not None:
		f.write("%s\t%s\t%d\t%s\n" % (source, target, target in source_to_targets[source], z))
	    else:
		f.write("%s\t%s\t%s\n" % (source, target, z))
    f.close()
    return


def get_proximity_values(parameters, sources, targets, dump_file=None):
    if dump_file is None:
	dump_file = parameters.get("proximity_file")  
    if os.path.exists(dump_file):
	try:
	    source_to_target_to_proximity, source_to_target_to_d = cPickle.load(open(dump_file))
	except: # For old dumps storing only z
	    raise ValueError("Update proximity dump to store d in addition to z!") # print
	    source_to_target_to_proximity = cPickle.load(open(dump_file))
	    source_to_target_to_d = None
	return source_to_target_to_proximity, source_to_target_to_d 
    source_to_target_to_proximity = {} # before source was stored as R string
    source_to_target_to_d = {} 
    f = open(dump_file + ".txt", 'w')
    f.write("source\ttarget\tz\td\n")
    for source in sources: 
	source_mod = text_utilities.convert_to_R_string(source) 
	source_to_target_to_proximity[source] = {}
	source_to_target_to_d[source] = {}
	for target in targets:
	    target_mod = text_utilities.convert_to_R_string(target)
	    #target_mod = target_mod.lower() 
	    out_file = parameters.get("output_dir") + "/%s_%s.out" % (source_mod, target_mod) 
	    if not os.path.exists(out_file):
		print "File not found:", out_file
		#continue 
		raise ValueError("Proximity values missing!")
	    z, d, m, s = open(out_file).readline().strip("\n").split()
	    source_to_target_to_proximity[source][target] = float(z)
	    source_to_target_to_d[source][target] = float(d)
	    f.write("%s\t%s\t%s\t%s\n" % (source, target, z, d))
    f.close()
    if dump_file is not None:
	cPickle.dump((source_to_target_to_proximity, source_to_target_to_d), open(dump_file, 'w'))
    return source_to_target_to_proximity, source_to_target_to_d


def get_guild_values(parameters, targets, source_to_genes, method='s', dump_file=None):
    if dump_file is None:
	dump_file = parameters.get("guild_file")  
    if os.path.exists(dump_file):
	source_to_target_to_score = cPickle.load(open(dump_file))
	return source_to_target_to_score 
    source_to_target_to_score = {} 
    f = open(dump_file + ".txt", 'w')
    f.write("source\ttarget\tscore\n")
    for target in targets:
	target_mod = text_utilities.convert_to_R_string(target)
	out_file = parameters.get("output_dir") + "/%s.n%s" % (target_mod, method) 
	if not os.path.exists(out_file):
	    print "File not found:", out_file
	    raise ValueError("GUILD values missing!")
	node_to_score = dict(line.strip("\n").split() for line in open(out_file).readlines())
	values = map(float, numpy.array(node_to_score.values()))
	m = numpy.mean(values)
	s = numpy.std(values)
	for source, genes in source_to_genes.iteritems(): 
	    score = numpy.mean([(float(node_to_score[gene]) - m) / s for gene in genes])
	    d = source_to_target_to_score.setdefault(source, {})
	    d[target] = score
	    f.write("%s\t%s\t%f\n" % (source, target, score))
    f.close()
    if dump_file is not None:
	cPickle.dump(source_to_target_to_score, open(dump_file, 'w'))
    return source_to_target_to_score


if __name__ == "__main__":
    main()

