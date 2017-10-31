import matplotlib.pyplot as plt


def histogram(values, bin_size=None, color=None, xlab=None, ylab="# of occurrence", title=None, xrot=None, xfont=None, out_file=None):
    fig = plt.figure()
    n, bins, patches = plt.hist(values, bins = bin_size, color = color) # color='g', alpha=0.75, normed=1, density=True
    if xlab is not None:
	plt.xlabel(xlab)
    plt.ylabel(ylab)
    if title is not None:
	plt.title(title)
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


