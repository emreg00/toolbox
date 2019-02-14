###
# Copied pasted code from previous projects
# Need to refactor significantly to create independent reusable functions
###

library(ggplot2)
library(reshape2)
library(scales)
library(RColorBrewer)
#source("heatmap3.R")

color.palette = c("blue", "orange", "green", "red", "grey20", "yellow") 
cb.palette <- c("blue", "red", "green", "grey20", "orange", "blue") #c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cb.bluered <- rep(c("blue", "red"), 20) 
cb.set <- brewer.pal(8,"Set1")
cb.pastel <- brewer.pal(8,"Pastel1")
cb.accent <- brewer.pal(8,"Accent")
cb.dark = brewer.pal(8, "Dark2")
cb.blues = brewer.pal(3,"Blues")
cb.reds = brewer.pal(3,"Reds")

txt.size = 18

### Statistical utilities
get.z.score<-function(val, values) {
    z = (val - mean(values)) / sd(values)
}


convert.p.value<-function(pval, or=1) {
    if ( or > 1 ) {
	z <- qnorm( pval / 2 );
    } else { 
	z <- -(qnorm( pval / 2 ));
    }
    return(z);
}


get.confidence.interval<-function(asample) {
    a = t.test(asample)
    m = mean(asample)
    s = sd(asample)
    q = qt(.975, df = length(asample)-1)
    #q = quantile(rt(1000000, length(asample)), probs=seq(0,1,by=0.025))["97.5%"]
    ci = m + q * s / sqrt(length(asample))
}


convert.z.score<-function(z, one.sided=NULL) {
    if(is.null(one.sided)) {
	#pval = pnorm(-abs(z), mean=mean(z), sd=sd(z))
	pval = pnorm(-abs(z));
	pval = 2 * pval
	pval[pval > 1] = 1
    } else if(one.sided=="-") {
	pval = pnorm(z);
    } else {
	pval = pnorm(-z);
    }
    return(pval);
}


get.similarity<-function(d, method="binary") {
    #"binary" "euclidian"
    if(startsWith(method, "jaccard")) {
	e = as.matrix(d)
	# Jaccard sum(a*b)/sum(a|b)
	f = e %*% t(e) # sum(a*b) 
	if(method == "jaccard") {
	    g = apply(e, 1, function(k) { rowSums(t(apply(e, 1, function(x) { x | k }))) }) # sum(a|b)
	}
	# Jaccard min
	if(method == "jaccard.min") {
	    g = apply(e, 1, function(k) { t(apply(e, 1, function(x) { min(sum(x), sum(k)) })) }) # min(sum(a), sum(b))
	}
	# Cosine
	#sum(a*b)/(sqrt(sum(a*a))*sqrt(sum(b*b)))
	jacc = f / g
    } else {
	jacc = dist(d, method = method, diag = TRUE, upper = TRUE) # binary manhattan euclidian
    }
    jacc = as.data.frame(as.matrix(jacc))
    rownames(jacc) = rownames(d)
    colnames(jacc) = rownames(d)

    if(method == "binary") { 
	jacc = 1 - jacc
    } else if(method == "euclidian") {
	jacc = (sqrt(ncol(d)) - jacc) / sqrt(ncol(d))
    }
    #View(jacc)
    #print(jacc)
    return(jacc)
}


### General visualization utilities
add.theme<-function(p, no.background = F, x.orthogonal = F, vertical.lines=F, text.size=18) {
    p = p + theme(plot.background = element_blank(), panel.border = element_blank(), line = element_line(size=1.3))
    if(no.background == T) {
	#p = p + theme_bw() 
	p = p + theme(panel.background = element_blank(), axis.ticks=element_blank(), axis.line = element_line(color = 'black', size=1.3)) 
	if(vertical.lines == T) {
	    p = p + theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour="gray", size=1.3))
	} else {
	    p = p + theme(panel.grid.major = element_blank())
	}
    } else {
	p = p + theme(panel.background = element_rect(fill = "ghostwhite"), panel.grid.major = element_line(colour="white", size=1.3), axis.line = element_blank(), axis.ticks = element_line(colour="darkgrey")) # lavender
    }
    if(x.orthogonal == T) {
	p = p + theme(axis.text.x = element_text(size = 8, angle=90, vjust=0.5, hjust=0.5))
    }
    p = p + theme(text = element_text(size = text.size+2), axis.text = element_text(color='black', size=text.size), axis.title.x = element_text(vjust=-0.5), axis.title.y = element_text(vjust=1.5), panel.grid.minor = element_blank()) 
    return(p)
}


print.plot<-function(p, out.file, landscape=F) {
    if(!is.null(out.file)) {
	if(landscape == T) {
	    svg(out.file, width = 8, height = 6, onefile = TRUE)
	} else {
	    svg(out.file) 
	}
	print(p)
	dev.off()
    }
}


draw.histogram<-function(d, variable, x.lab, y.lab, binwidth=NULL, x.scale=NULL, y.scale=NULL, text.size=18, coord.flip=F, out.file=NULL) {
    if(is.null(binwidth)) {
	p = ggplot(d, aes_string(variable)) + geom_histogram(color="grey95", fill=cb.palette[1], alpha = 0.8) 
	#x = d[[variable]]
	#binwidth = round(abs(max(x) - min(x)) / 30)
	#print(binwidth)
    } else {
	p = ggplot(d, aes_string(variable)) + geom_histogram(color="grey95", fill=cb.palette[1], alpha = 0.8, binwidth=binwidth) # aes(y = ..count..)
    }
    p = p + labs(x = x.lab, y = y.lab) 
    if(!is.null(x.scale)) {
	if(x.scale == "sqrt") {
	    p = p + scale_x_sqrt() 
	}
	if(x.scale == "log") {
	    p = p + scale_x_log10() + annotation_logticks(sides="b")
	}
	if(x.scale == "discrete") {
	    n = max(d[[variable]])
	    x.seq = c(1, seq(5, n+1, by=5))
	    p = p + scale_x_discrete(breaks = x.seq) 
	}
    }
    if(!is.null(y.scale)) {
       if(y.scale == "log") {
	    p = p + scale_y_log10() + annotation_logticks(sides="l")
	    #p = p + coord_trans(y="log")
	}
    }
    if(coord.flip==T) {
	p = p + coord_flip()
    }
    p = add.theme(p, no.background=F, text.size=text.size)

    print.plot(p, out.file)
    print(variable)
    print(summary(d[[variable]]))
    return(p)
}


draw.scatterplot<-function(d, variable, selected, x.text, y.text, x.log=F, y.log=F, var.color=NULL, regression.line=F, out.file=NULL) {
    p = ggplot(data = d, aes_string(x = variable, y = selected))
    if(is.null(var.color)) {
	p = p + geom_point(color=cb.palette[1], alpha=0.5)
    } else {
	p = p + geom_point(aes_string(color=var.color), alpha=0.5) + guides(color=guide_legend(var.color)) + scale_color_manual(values=cb.palette)  
    }
    if(x.log == T) {
	p = p + scale_x_log10() + annotation_logticks(sides="b") 
    }
    if(y.log == T) {
	p = p + scale_y_log10() + annotation_logticks(sides="l") # + coord_trans(ytrans = "log10") 
    }
    if(regression.line==T) {
	p = p + geom_smooth(method = "lm", se=FALSE, color="gray", size=1.3, alpha=0.5, formula = y ~ x)
	f = paste(selected, "~", variable)
	m = lm(f, d)
	r2 = format(summary(m)$r.squared, digits = 3)
	print(sprintf("R^2: %s", r2))
	lm_eqn = function(df){
	    m = lm(y ~ x, df);
	    eq <- substitute("~~italic(R)^2~"="~r2", #substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
		  #list(a = format(coef(m)[1], digits = 2), 
		  #b = format(coef(m)[2], digits = 2), 
		  list(r2 = format(summary(m)$r.squared, digits = 3)))
	    as.character(as.expression(eq));                 
	}
	#p = p + geom_text(aes(x = -3, y = 250, label = lm_eqn(d)), parse = TRUE) 
	#p = p + geom_text(aes(x = -3, y = 250, label = paste("R^2", "=", r2, sep="")), parse = TRUE) 
	p = p + annotate("text", x = mean(d[[variable]],na.rm=T)+1*sd(d[[variable]],na.rm=T), y = mean(d[[selected]], na.rm=T)+1*sd(d[[selected]],na.rm=T), label = paste("R^2", "=", r2, sep="")) # -4, 60
    }
    p = p + labs(x=x.text, y=y.text) 
    p = add.theme(p, no.background = T)
    print.plot(p, out.file)
    a = cor.test(d[,variable], d[,selected], method="spearman")
    print(sprintf("Cor.test: %.3f %e", a$estimate, a$p.value))
    return(p)
}


draw.violinplot<-function(d, variable, value, x.lab, y.lab, y.log=F, out.file=NULL) {
    if(y.log == T) {
	d[,value] = log10(d[[value]])
    }
    p = ggplot(data=d, aes_string(x=variable, y=value)) + geom_violin(alpha=0.5, color="gray") + geom_jitter(alpha=0.5, aes_string(color=variable), position = position_jitter(width = 0.1)) + stat_summary(fun.y="median", geom='point', color='black', size=10, alpha=0.5) + coord_flip() + scale_color_manual(values=cbPalette.set) 
    p = p + labs(y = y.lab, x="") + guides(color=F) 
    p = add.theme(p)

    #library(beanplot)
    #beanplot(y$z, x$z, overalline="median",  col=cb.palette[1:2], ylab="Closeness (z)", names=c("Symptomatic", "Rest"), ll=0.05)
    #beanplot(y$ri, x$ri, overalline="median",  col=cb.palette[1:2], ylab="Relative inefficacy (%)", names=c("Symptomatic", "Rest"), ll=0.05) #, boder=NA, frame.plot=F)
    print.plot(p, out.file)
    return(p)
}


draw.distribution.with.fit<-function(d, column="degree", xlab="Degree (k)", ylab="Cumulative count", out.file=NULL) {
    x = table(d[[column]])
    e = data.frame(k=as.numeric(names(x)), count=as.vector(x))
    e$count.cum = rev(cumsum(rev(e$count)))
    p = ggplot(data=e, aes(x=log10(k), y=log10(count.cum))) + geom_point(alpha=0.5) + scale_y_continuous(name = ylab, labels = math_format(10^.x)) + scale_x_continuous(name = xlab, labels = math_format(10^.x)) 
    p = p + geom_smooth(method=glm) 
    p = p + scale_color_manual(values=cb.palette) #+ theme_bw()
    p = add.theme(p)
    p = p + annotation_logticks(sides="lb")

    if(!is.null(out.file)) {
	
	svg(out.file) 
	print(p)
	dev.off()

    	out.file = paste(out.file, ".more", sep="")
	e$p.cum = e$count.cum / sum(e$count.cum)
	e$k.log = log(e$k)
	e$count.cum.log = log(e$count.cum)
	e$p.cum.log = log(e$p.cum)
	model.log = lm(count.cum.log ~ k.log, data=e)
	#model.log = lm(p.cum.log ~ k.log, data=e)
	model = glm(count.cum ~ k, family=quasi(link=power()), data=e)
	#model = glm(p.cum ~ k, family=quasi(link=power()), data=e)
	k.seq = seq(1, max(e$k), by=0.05)
	fit.log = predict(model.log, data.frame(k.log=log(k.seq)), type="response")
	fit = predict.glm(model, data.frame(k=k.seq), type="response")
	svg(out.file) 
	par(mfrow=c(2,2))
	# Linear fit on log log
	x = as.vector(coef(model.log))
	b = x[2]
	a = x[1]
	plot(count.cum.log~k.log, data=e)
	abline(a, b, col=2)
	# Power law fit
	a = as.vector(coef(model))
	b = x[2]
	a = x[1]
	plot(count.cum~k, data=e, log="xy")
	abline(log(a), b, col=2)

	library("poweRlaw")
	e.pl = displ$new(d$degree)
	est = estimate_xmin(e.pl) 
	e.pl$setXmin(est) 
	plot(e.pl)
	lines(e.pl, col=2)
	est = estimate_xmin(e.pl, pars = seq(1, 2.5, 0.1)) 
	e.pl$setXmin(est) 
	lines(e.pl, col=3)

	dev.off()
    }
    return(p)
}


draw.roc.curve<-function(d, titles, value.column="score", label.column="flag", method="", plot.type="roc") {
    if(plot.type == "roc") { 
	# Get ROC using all instances in ROC (average over all repetitions and folds)
	pred = prediction(d[[value.column]], d[[label.column]])
	perf = performance(pred, "tpr", "fpr")
	if(i == 1) {
	    plot(perf, col=cb.dark[i], type="l", lwd=3, xlab="False Positive Rate", ylab="True Positive Rate") #, ylim=c(0,1) # cex.lab=1.5, cex.axis=1.5, bty="n", 
	} else {
	    p = plot(perf, lwd=3, add=T, col=cb.dark[i])
	}
	#out.file = paste0(data.dir, "/roc.dat")
	#write(mapply(paste, perf@x.values, perf@y.values), out.file)
	perf = performance(pred, "auc")
	auc = 100*perf@y.values[[1]]
	label = title
	#label = paste0(title, sprintf("(AUC: %.0f%%)", auc))
	# Get ROC per repetition
	#a = plyr::ddply(d, .(rep), function(x) { pred=prediction(split(x[[value.column]], x$fold), split(x[[label.column], x$fold)); perf=performance(pred, "fpr", "tpr"); mapply(paste, perf@x.values, perf@y.values) }) # folds contain different number of points
	a = plyr::ddply(d, .(rep), function(x) { pred=prediction(x[[value.column]], x[[score.column]]); perf=performance(pred, "fpr", "tpr"); mapply(paste, perf@x.values, perf@y.values) })
	b = do.call(rbind, by(a, a$rep, rbind, ""))
	#out.file = paste0(data.dir, "/roc_per_rep.dat")
	#write.table(b[,2], out.file, row.names=F, quote=F, col.names=F)
    } else {
	pred = prediction(d[[value.column]], d[[label.column]])
	perf = performance(pred, "sens")
	#output.dir = paste0(data.dir, "/") 
	#out.file = paste0(data.dir, "/sens.dat")
	#write(mapply(paste, perf@x.values, perf@y.values), out.file)
	p = plot(perf, col=cb.dark[1], type="l", lwd=3, xlab="Fraction of module genes", ylab="Sensitivity / Specificity") 
	perf = performance(pred, "spec")
	#out.file = paste0(data.dir, "/spec.dat")
	#write(mapply(paste, perf@x.values, perf@y.values), out.file)
	p = plot(perf, lwd=3, add=T, col=cb.dark[2])
	legend(quantile(perf@x.values[[1]])[4], quantile(perf@y.values[[1]])[4], legend=c("sensitivity", "specificity"), col=cb.dark[1:2], lty=rep(1,2), lwd=rep(3, 2), bty="n") #, cex=1.5) 
    }
    if(plot.type == "roc") { 
	legend("bottomright", legend=label, col=cb.dark[1:1], lty=rep(1,1), lwd=rep(3, 1), bty="n") #, cex=1.5) # -3, 60
    }
    return(p)
}



### Expression data visualization
volcano.plot<-function(de, cutoff.logFC=1, cutoff.adj.P.Val=0.05, title=NULL, out.file=NULL) {
    de$threshold = as.factor(abs(de$logFC) > cutoff.logFC & de$adj.P.Val < cutoff.adj.P.Val) 
    #de$threshold = as.factor(abs(de$logFC) > cutoff.logFC & de$P.Value < cutoff.adj.P.Val) # P-value based cutoffing
    p = ggplot(data=de, aes(x=logFC, y=-log2(P.Value), colour=threshold)) + geom_point(alpha=0.4, size=1.75) 
    p = p + xlim(c(-5, 5)) + ylim(c(0, 15)) + scale_colour_manual(values=cb.dark) # c(-10, 10)
    p = p + guides(colour=F) + labs(y = "-log2 P-value", x = "log2 Fold Change")
    #gene1 = de[1:3,]
    # g + geom_text(aes(x=gene1$logFC, y=-log10(gene1$P.Value),                      label=gene1$ID, size=1.2), colour="black")
    if(!is.null(title)) {
	p = p + ggtitle(title)
    }
    if(!is.null(out.file)) {
	if(endsWith(out.file, "jpg")) {
	    jpeg(out.file, width = 4, height = 4, units = 'in', res=300)
	} else {
	    pdf(out.file)
	}
	print(p)
	dev.off()
    } else {
	print(p)
    }
    return(p)
}


get.sample.ordering <- function(sample.mapping, states, title=NULL, out.file=NULL) {
    groups = split(sample.mapping$sample, sample.mapping$type)
    samples.ordered = c()
    for(state in states) {
	#arrays<-sample.mapping[sample.mapping$Type == state,]$Sample
	samples = as.vector(groups[[state]])
	#samples = samples(order(samples))
	samples.ordered = c(samples.ordered, samples)
	#print(c(state, samples))
    }
    return(samples.ordered)
}


visualize.array.correlation.grouped.by.sample.type <- function(expr, sample.mapping, states=NULL, title=NULL, out.file=NULL) {
    library(gplots)
    if(!is.null(states)) {
	samples.ordered = get.sample.ordering(sample.mapping, states) # when case / control scheme is used
	d = cor(expr)
	e = d[samples.ordered, samples.ordered]
    }
    e = expr[, colnames(expr) %in% sample.mapping$sample]
    e = cor(e)
    #print(dim(e))
    #val.cols <- heat.colors(5) 
    val.cols <- cb.blues # greenred
    library(gdata)
    map = lapply(mapLevels(sample.mapping$type), function(x) { x+1 })
    #val.cols.side = ifelse(sample.mapping$type == "N/A", "2", "3")
    val.cols.side = as.character(map[sample.mapping$type])
    if(is.null(out.file)) {
	p = heatmap(e, Rowv=NA, Colv=NA, col=val.cols, revC=T) 
	#heatmap(e, Rowv=NA, Colv=NA, revC=T, symm=T, margins=c(6,6), col=val.cols, RowSideColors=val.cols.side, ColSideColors=val.cols.side)
	# key=T, scale="none", cexRow=1, cexCol=1, #symkey=T, symbreaks=T, KeyValueName="Correlation" RowSideColors=val.cols.side
	#heatmap.3
    } else {
	pdf(out.file) 
	if(F) { # order based on type
	    p = heatmap.2(e, dendrogram="col", Rowv=F, Colv=F, symm=T, revC=F, na.color="grey", col=val.cols, ColSideColors=val.cols.side, margins=c(6,6), trace="none", main=title, density.info="none", keysize=1, key.title=F) # main="Sample correlation"   
	} else { # order based on dendongram
	    p = heatmap.2(e, dendrogram="col", Rowv=T, Colv=T, symm=T, revC=T, na.color="grey", col=val.cols, ColSideColors=val.cols.side, margins=c(6,6), trace="none", main=title, density.info="none", keysize=1, key.title=F) 
	}
	dev.off()
    }
    # For gene vs sample heatmap
    # heatmap.2(as.matrix(expr), dendrogram="col", Rowv=T, Colv=T, na.color="grey", col=val.cols, ColSideColors=val.cols.side, trace="none", density.info="none", keysize=1, key.title=F)
    return(p)
}


visualize.array.distribution.grouped.by.sample.type <- function(expr, sample.mapping, states, title=NULL, out.file=NULL) {
    if(!is.null(states)) {
	samples.ordered = get.sample.ordering(sample.mapping, states) # when case / control scheme is used
	e = expr[, colnames(expr) %in% samples.ordered]
    } else {
	e = expr[, colnames(expr) %in% sample.mapping$sample]
    }
    d = melt(e)
    #p = boxplot(e, col=factor(sample.mapping$type))
    p = ggplot(data=d, aes(variable, value)) + geom_boxplot(aes(fill=sample.mapping[variable,"type"]))
    p = p + scale_fill_manual(values=cb.dark) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top")
    p = p + guides(fill=guide_legend(title.position="left", title="Type")) + labs(y="log2 Expression", x="Samples")
    if(!is.null(title)) {
	p = p + ggtitle(title)
    }
    if(!is.null(out.file)) {
	if(endsWith(out.file, "jpg")) {
	    jpeg(out.file, width = 4, height = 4, units = 'in', res=300)
	} else {
	    pdf(out.file)
	}
	print(p)
	dev.off()
    } else {
	print(p)
    }
    return(p)
}


visualize.pca.grouped.by.sample.type <- function(expr, sample.mapping, states, title=NULL, out.file=NULL) {
    library(ggbiplot)
    if(!is.null(states)) {
	samples.ordered = get.sample.ordering(sample.mapping, states) # when case / control scheme is used
	e = expr[, colnames(expr) %in% samples.ordered]
    } else {
	e = expr[, colnames(expr) %in% sample.mapping$sample]
    }
    pca <- prcomp(t(e), center = TRUE, scale. = TRUE) 
    print(summary(pca))
    p = ggbiplot(pca, choices=c(1,2), obs.scale = 1, var.scale = 1, groups = sample.mapping$type, labels=colnames(e), ellipse=F, circle=F, var.axes=F)
    #p = p + scale_color_discrete(name = '') + theme(legend.direction = 'horizontal', legend.position = 'top')
    p = p + scale_color_manual(values=cb.dark) + guides(color=guide_legend(title.position="left", title="Type")) + theme(legend.position="top")
    if(!is.null(title)) {
	p = p + ggtitle(title)
    }
    if(!is.null(out.file)) {
	if(endsWith(out.file, "jpg")) {
	    jpeg(out.file, width = 4, height = 4, units = 'in', res=300)
	} else {
	    pdf(out.file, width=10, height=6)
	}
	print(p)
	dev.off()
    } else {
	print(p)
    }
}


draw.heatmap<-function(d, plot.type="tile", similarity.conversion=NULL, txt.size=txt.size, d.label=NULL, id.col=NULL, flip=F) {
    if(!is.null(similarity.conversion)) {
	d = get.similarity(d, similarity.conversion) #"binary" "euclidian"
    }
    if(plot.type == "tile") {
	# Code refactored from https://journocode.com/2016/03/13/similarity-and-distance-in-data-part-2/
	if(is.null(id.col)) {
	    d$name = rownames(d)
	    id.col = "name"
	}
	f = melt(d, id.vars = id.col)	
	if(!is.null(similarity.conversion)) {
	    # Order by name
	    labels = rownames(d)[order(rownames(d))]
	    f[[id.col]] = factor(f[[id.col]], labels)
	    f$variable = factor(f$variable, rev(labels))
	    f = plyr::arrange(f, variable, plyr::desc(!!sym(id.col)))
	    if(!is.null(d.label)) { 
		stop("External labels are not supported for similarity!")
	    }
	} else {
	    labels = levels(f$variable[order(f$variable)])
	    f$variable = factor(f$variable, rev(labels))
	    if(!is.null(d.label)) { 
		d.label[[id.col]] = rownames(d.label)
		d.label = melt(d.label, id.vars = id.col)
		f$label = d.label$value
	    }
	}
	#print(f)
	if(!is.null(d.label)) { 
	    p = ggplot(f, aes_string(id.col, "variable", label="label"))
	} else {
	    p = ggplot(f, aes_string(id.col, "variable"))
	}
	p = p + geom_tile(aes(fill=value), colour = "white")
	p = p + scale_fill_gradient(low="white", high="red") #low = "#b7f7ff", high = "#0092a3")
	p = p + theme_light(base_size = txt.size) + labs(x = "", y = "") + guides(fill=guide_legend(title=NULL)) 
	p = p + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
	if(flip) {
	    p = p + coord_flip() 
	}
	if(!is.null(d.label)) { 
	    p = p + geom_text(size = txt.size/2)
	}
	p = p + theme(axis.ticks = element_blank(), axis.text.x = element_text(size = txt.size * 0.8, angle = 310, hjust = 0), 
		  axis.text.y = element_text(size = txt.size * 0.8))
    } else if(plot.type == "simplot") { 
	p = simplot.mod(d)
    } else if(plot.type == "pheatmap") { 
	library(pheatmap)
	p = pheatmap(d, fontsize = txt.size)
    } else if(plot.type == "heatmap2") { 
	library(gplots)
	val.cols = c("white", cb.reds)
	a = colnames(d)
	i = regexpr("([0-9]+)", a)
	k = lapply(a, function(x) { k=which(a==x); substr(x, i[k][1], (i[k]+attr(i, "match.length")[k]-1)) } )
	k = as.numeric(k)
	val.cols.side = color.palette[k]
	title = ""
	margins = c(16,16) # c(6,6)
	p = heatmap.2(as.matrix(d), dendrogram="col", Rowv=T, Colv=T, symm=T, revC=T, na.color="grey", col=val.cols, ColSideColors=val.cols.side, margins=margins, trace="none", main=title, density.info="none", keysize=1, key.title=F) 
    }
    return(p)
}


simplot.mod<-function(sim, xlab = "", ylab = "", color.low = "white", color.high = "red", labs = TRUE, digits = 2, labs.size = 5, font.size = 10) 
{
    # Modified version of the simplot function in DOSE package
    sim.df <- as.data.frame(sim)
    rn <- row.names(sim.df)
    sim.df <- cbind(ID = rownames(sim.df), sim.df)
    #print(sim.df)
    sim.df <- melt(sim.df, na.rm=T)
    #print(sim.df)
    sim.df[, 1] <- factor(sim.df[, 1], levels = rev(rn))
    if (labs == TRUE) {
        sim.df$label <- as.character(round(sim.df$value, digits))
    }
    variable <- ID <- value <- label <- NULL
    if (labs == TRUE) {
        p <- ggplot(sim.df, aes(variable, ID, fill = value, label = label))
    }
    else {
	p <- ggplot(sim.df, aes(variable, ID, fill = value))
    }
    p <- p + geom_tile(color = "black") + scale_fill_gradient(low = color.low, high = color.high, name="Similarity") 
    p <- p + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))  # x = c(0.1, 0)
    if (labs == TRUE) 
        p <- p + geom_text(size = labs.size)
    p <- p + theme(axis.text.x = element_text(size = font.size), axis.text.y = element_text(size = font.size))
    #p <- p + theme_dose(font.size)
    #p <- p + theme(axis.text.x = element_text(hjust = 0, angle = -90)) + theme(axis.text.y = element_text(hjust = 0))
    #p <- p + theme(legend.title = element_blank())
    #p <- p + xlab(xlab) + ylab(ylab)
    #p <- p + theme(axis.text.x = element_text(vjust = 0.5))
    #p <- p + theme_bw()
    #p <- p + theme_classic()
    p <- p + xlab(NULL) + ylab(NULL) 
    p <- p + theme(axis.text.x=element_blank(), axis.ticks=element_blank(), panel.border=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # axis.text.y=element_blank(), 
    #p <- p + theme(legend.justification = c(1, 0), legend.position = c(0.6, 0.1), legend.direction = "horizontal")
    p <- p + guides(fill=F)
    #p <- p + annotate("text", x = rev((1:nrow(sim))-2), y = 1:nrow(sim), label = rownames(sim))
    return(p)
}


### BELOW NOT ADOPTED / NOT TESTED

draw.two.sided.barplots<-function() {
    file.name = paste(data.dir, prefix, ".atc", sep="")
    d = read.table(file.name, header=T, sep="\t")
    d$category = paste(d$category, " (", d$atc, ")", sep="")
    e = melt(d[,c(2,3,4)], "category")
    #e$value.close = ifelse(e$variable=="close", e$value, -e$value)
    for(i in e$category) {
	val = d[d$category==i,]$close-d[d$category==i,]$not.close
	e[e$category==i,"value.diff"] = val - 0.1*d[d$category==i,]$not.close
    }
    #e <- transform(e, category = reorder(category, value.close))
    e.close = transform(e[e$variable=="close",], category=reorder(category, value.diff))
    e.not.close = e[e$variable=="not.close",]
    #ggplot(data=e, aes(category, value, group=variable)) + geom_bar(aes(fill=variable), stat="identity", position="dodge") + coord_flip() + theme_bw()
    p = ggplot() + geom_bar(data=e.close, aes(category, value), fill="red", stat="identity", alpha=0.5) + geom_bar(data=e.not.close, aes(category, -value), fill="blue", stat="identity", alpha=0.5) + coord_flip() + theme_bw()
    p = p + theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank())
    p = p + theme(axis.ticks=element_blank(), axis.text.y = element_text(size = 8)) + scale_y_discrete(limits=seq(-15,15,by=5), labels=abs) # mesh69
    #p = p + theme(axis.ticks=element_blank(), axis.text.y = element_text(size = 7)) + scale_y_discrete(labels=abs, limits=seq(-30,30,by=10)) # ,breaks=seq(-30,30,by=10) # mesh299
    p = p + labs(x=NULL, y=NULL) #, axis.line = element_line(color = 'black')) #x="ATC code", y="Number of drugs") # Number of occurences of ATC terms 
    out.file = paste(img.dir, "atc_enrichment_mesh69_ndfrt_m2m_subset_seeds_100000.svg", sep="")
    #out.file = paste(img.dir, "atc_enrichment_mesh299_ndfrt_m2m_subset_seeds.svg", sep="")
    svg(out.file)
    print(p)
    dev.off()
}


