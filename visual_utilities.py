import matplotlib.pyplot as plt


def histogram(values, bin_size=None, color=None, xlab=None, ylab="# of occurrence", title=None, xrot=None, xfont=None, mark_loc=None, mark_text="", out_file=None):
    fig = plt.figure()
    if mark_loc is not None:
	ax = fig.add_subplot(111)
    n, bins, patches = plt.hist(values, bins = bin_size, color = color) # color='g', alpha=0.75, normed=1, density=True
    if xlab is not None:
	plt.xlabel(xlab)
    plt.ylabel(ylab)
    if title is not None:
	plt.title(title)
    if mark_loc is not None:
	ax.annotate(mark_text, xy=(mark_loc[0], 0), xytext=(mark_loc[0], mark_loc[1]), arrowprops=dict(facecolor='red', shrink=0.05)) 
	#ax.set_xlim(min(values), max(values))
    #plt.axis([xmin, xmax, ymin, ymax])
    #plt.grid(True)
    if xrot is not None:
	plt.xticks(rotation=xrot)
    if xfont is not None:
	plt.xticks(fontsize=xfont)
    if out_file is not None:
	fig.savefig(out_file)
	plt.close(fig)
    else:
	plt.show()
    return


def venn_diagram(values, labels = ('A', 'B'), title = "", colors=None, out_file = None):
    if len(values) == 2:
	a, b = values
	aa = a - b
	bb = b - a 
	ab = a & b
	values = (len(aa), len(bb), len(ab))
    from matplotlib_venn import venn2
    fig = plt.figure() #figsize=(4,4))
    v = venn2(subsets = values, set_labels = labels)
    if colors is not None:
	v.get_patch_by_id('10').set_color(colors[0])
	v.get_patch_by_id('01').set_color(colors[1])
	try:
	    v.get_patch_by_id('11').set_color(colors[2])
	    v.get_patch_by_id('11').set_edgecolor('none')
	    v.get_patch_by_id('11').set_alpha(0.4)
	except:
	    pass
    plt.title(title)
    if out_file is not None:
	fig.savefig(out_file)
	plt.close(fig)
    else:
	plt.show()
    return

