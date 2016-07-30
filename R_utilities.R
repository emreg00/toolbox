###
# Copied pasted code from previous projects
# Need to refactor significantly to create independent reusable functions
###

library(ggplot2)
library(reshape2)
library(scales)
library(RColorBrewer)
library(beanplot)
#source("heatmap3.R")

cbPalette <- c("blue", "red", "green", "grey20", "orange", "blue") #c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette2 <- rep(c("blue", "red"), 20) 
cbPalette3 <- brewer.pal(9,"Blues")
cbPalette.set <- brewer.pal(8,"Set1")
cbPalette.pastel <- brewer.pal(8,"Pastel1")
cbPalette.accent <- brewer.pal(8,"Accent")

### Analytical functions
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


convert.z.score<-function(z, one.sided=NULL) {
    if(is.null(one.sided)) {
	pval = pnorm(-abs(z));
	pval = 2 * pval
    } else if(one.sided=="-") {
	pval = pnorm(z);
    } else {
	pval = pnorm(-z);
    }
    return(pval);
}

### Visualization methods
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
	p = ggplot(d, aes_string(variable)) + geom_histogram(color="grey95", fill=cbPalette[1], alpha = 0.8) 
	#x = d[[variable]]
	#binwidth = round(abs(max(x) - min(x)) / 30)
	#print(binwidth)
    } else {
	p = ggplot(d, aes_string(variable)) + geom_histogram(color="grey95", fill=cbPalette[1], alpha = 0.8, binwidth=binwidth) # aes(y = ..count..)
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
	p = p + geom_point(color=cbPalette[1], alpha=0.5)
    } else {
	p = p + geom_point(aes_string(color=var.color), alpha=0.5) + guides(color=guide_legend(var.color)) + scale_color_manual(values=cbPalette)  
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

    #beanplot(y$z, x$z, overalline="median",  col=cbPalette[1:2], ylab="Closeness (z)", names=c("Symptomatic", "Rest"), ll=0.05)
    #beanplot(y$ri, x$ri, overalline="median",  col=cbPalette[1:2], ylab="Relative inefficacy (%)", names=c("Symptomatic", "Rest"), ll=0.05) #, boder=NA, frame.plot=F)
    print.plot(p, out.file)
    return(p)
}


draw.degree.distribution.with.fit<-function(d, column="degree", out.file=NULL) {
    x = table(d[[column]])
    e = data.frame(k=as.numeric(names(x)), count=as.vector(x))
    e$count.cum = rev(cumsum(rev(e$count)))
    p = ggplot(data=e, aes(x=log10(k), y=log10(count.cum))) + geom_point(alpha=0.5) + scale_y_continuous(name = "Cumulative count", labels = math_format(10^.x)) + scale_x_continuous(name = "Degree (k)", labels = math_format(10^.x)) 
    p = p + geom_smooth(method=glm) 
    p = p + scale_color_manual(values=cbPalette) #+ theme_bw()
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


### BELOW NOT ADOPTED / NOT TESTED
### Statistical methods
get.Fishers.enrichment<-function(file.name) {

    d=read.table(file.name, header=T)

    method = "z" #"z.tissue" #"z.pathway" #"pval" 
    cutoff = -1.28 #-1.61 #-0.84 #-2 #0.05 #-2 

    #d$m.e.z = scale(-d$m.e)
    # e = d
    e = d[!is.na(d[[method]]),]
    x=e[e$flag=="True", ]
    #y=e[e$flag=="False", ]
    y=e[e$flag!="True", ]

    #method = "z" #"z.tissue" 

    a = sum(x[[method]]<=cutoff, na.rm=T)
    b = sum(y[[method]]<=cutoff, na.rm=T) 
    #a = sum(x[[method]]>=cutoff, na.rm=T) # n.overlap
    #b = sum(y[[method]]>=cutoff, na.rm=T) 
    contigency = matrix(c(a, nrow(x)-a, b, nrow(y)-b), 2, 2)
    print(contigency)
    print(sum(contigency))
    s = fisher.test(contigency)
    print(s)

    print(sum(x[[method]]<=cutoff & x[["n.overlap"]]>0))

    return();
}


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


### Expression data analysis
visualize.array.correlation.grouped.by.sample.type <- function(expr, sample.mapping, states, out.file) {
    groups = split(sample.mapping$sample, sample.mapping$type)
    d<-cor(expr)
    samples.ordered = c()
    for(state in states) {
	#arrays<-sample.mapping[sample.mapping$Type == state,]$Sample
	samples = as.vector(groups[[state]])
	#samples = samples(order(samples))
	samples.ordered = c(samples.ordered, samples)
	print(c(state, samples))
	e<-d[samples.ordered, samples.ordered]
    }
    val.cols <- heat.colors(5) #brewer.pal(9,"Blues") # greenred
    svg(out.file) 
    heatmap(e, Rowv=NA, Colv=NA, col=val.cols, scale="none", revC=T) #labCol=F, labRow=F, , margins=c(1,1), keysize=0.9, colsep=seperators, rowsep=seperators, sepcolor="white")
    #heatmap.3(e, labCol=colnames(e), dendrogram="none", Rowv=F, Colv=F, revC=F, scale="none", na.color="grey", col=val.cols, RowSideColors=row.cols, ColSideColors=column.cols, margins=c(6,20), cexRow=0.8, cexCol=0.8, trace="none", xlab="Drug-disease closeness (z)", NumRowSideColors=0.5, NumColSideColors=0.5, symbreaks=symetry, key=T, symkey=symetry, density.info="none") # KeyValueName="Log expression")
    dev.off()
}


find.differentially.expressed.genes<-function(expr, sample.mapping, states, gene.mapping, out.file) {
    library(limma)
    design = model.matrix(~ 0 + sample.mapping$Type)
    colnames(design) = gsub(" ", "_", states) #colnames(design))
    fit = lmFit(expr, design)
    ref = states[1]
    contrast = unlist(sapply(states, function (x) { if(ref != x) { paste(x, ref, sep = "-") } }))
    contrast = gsub(" ", "_", contrast)
    #contrast = c()
    #for(i in 1:length(states)) { 
    #	for(j in 1:length(states)) { 
    #	    if(i<j) { 
    #		contrast = c(contrast, paste(gsub(" ", "_", states[j]), gsub(" ", "_", states[i]), sep = "-"))
    #	    }
    #	}
    #}
    cont.matrix = makeContrasts(contrasts=contrast, levels=design)
    fit2 = contrasts.fit(fit, cont.matrix)
    fit2 = eBayes(fit2)
    results = decideTests(fit2, p.value=0.05)
    vennDiagram(results)

    d = topTable(fit2, coef=1, adjust="BH", sort.by="B", number=50000)
    #d = d[d$P.Value <= 0.05,]
    d$geneid = as.vector(gene.mapping[as.integer(rownames(d)), "Gene.ID"])
    d$gene = as.vector(gene.mapping[as.integer(rownames(d)), "Gene.symbol"])
    write.table(d, file=out.file, row.names=F, sep="\t", quote=F) 
    #d$type = ifelse(d$logFC>0, "up", "down")
    #d = d[d$adj.P.Val <= 0.05,c("gene", "type")]
    #write.table(d, file=paste(out.file, ".fdr5", sep=""), row.names=F, sep="\t", quote=F) #, col.names=F)
    #return(fit2)
    out.file = sub("\\.dat", "_volcano.svg", out.file)
    draw.volcano.plot(d, out.file)
}


draw.volcano.plot <- function(d, out.file) {
    d$flag = ifelse(abs(d$logFC) >= log(1.2, 2), ifelse(d$logFC >= log(1.2, 2), "up", "down"), "no change")
    p = ggplot(data=d, aes(logFC, -log(adj.P.Val,2))) + geom_point(aes(color=factor(flag))) + scale_color_manual(values=c("green", "grey60", "red")) + geom_hline(yintercept=-log(0.05,2), linetype="dashed", color="black") + theme_bw()
    svg(out.file)
    print(p)
    dev.off()
}



find.differentially.expressed.genes.sam<-function(expr, sample.mapping, states, out.file=NULL, state.background=NULL, adjust.method='BH', cutoff=0.2) {
    library(samr)
    if(is.null(state.background)) {
	state.background = "case"
    }
    expr = expr[,as.vector(sample.mapping$sample)]
    x = as.matrix(expr)
    y = ifelse(sample.mapping$type == state.background, 1, 2)

    samfit<-invisible(SAM(x, y, resp.type="Two class unpaired", fdr.output=cutoff))
    d = rbind(samfit$siggenes.table$genes.up, samfit$siggenes.table$genes.lo)
    d = data.frame(GeneID=d[,2], logFC=as.double(d[,6]), adj.P.Val=as.double(d[,7])/100)
    d = d[order(d$adj.P.Val),]

    if(!is.null(out.file)) {
	write.table(d, file=out.file, row.names=F, quote=F, sep="\t")
    }
    return(d) 
}


find.differentially.expressed.genes.welch<-function(expr, sample.mapping, states, out.file=NULL, state.background=NULL, adjust.method='BH', cutoff=0.2) {
    if(is.null(state.background)) {
	state.background = "case"
    }
    expr = expr[,as.vector(sample.mapping$sample)]
    samples.background = which(sample.mapping$type == state.background)
    samples = setdiff(1:ncol(expr), samples.background)
    p.values = apply(expr, 1, function(x) { b = t.test(x[samples], x[samples.background]); return(b$p.value) })
    fc = apply(expr, 1, function(x) { log(mean(x[samples]) / mean(x[samples.background])) })

    d = data.frame(GeneID=names(p.values), logFC=as.vector(fc), P.Value=as.vector(p.values), adj.P.Val=p.adjust(p.values, method=adjust.method))
    d = d[order(d$adj.P.Val),]

    if(!is.null(out.file)) {
	write.table(d, file=out.file, row.names=F, quote=F, sep="\t")
    }
    return(d) 
}


find.differentially.expressed.genes.limma<-function(expr, sample.mapping, states, out.file=NULL, state.background=NULL, adjust.method='BH', cutoff=0.2) {
    library(limma)
    #if(ncol(expr) != nrow(sample.mapping)) {
    #	print(c("Warning: inconsistent dimenstions!", ncol(expr), nrow(sample.mapping)))
    #}
    expr = expr[,as.vector(sample.mapping$sample)]
    design = model.matrix(~ 0 + sample.mapping$type)
    colnames(design) = gsub(" ", "_", states) #colnames(design))
    fit = lmFit(expr, design)
    #contrast = unlist(sapply(states, function (x) { if(ref != x) { paste(ref, x, sep = "-") } }))
    contrast = c()
    if(!is.null(state.background)) {
	# W.r.t. background
	for(i in 1:length(states)) { 
	    if(states[i] != state.background)
		contrast = c(contrast, paste(gsub(" ", "_", states[i]), gsub(" ", "_", state.background), sep = "-"))
	}
    } else {
	# All vs all DE test
	for(i in 1:length(states)) { 
	    for(j in 1:length(states)) { 
		if(i<j) { 
		    contrast = c(contrast, paste(gsub(" ", "_", states[i]), gsub(" ", "_", states[j]), sep = "-"))
		}
	    }
	}
    }
    cont.matrix = makeContrasts(contrasts=contrast, levels=design)
    fit2 = contrasts.fit(fit, cont.matrix)
    fit2 = eBayes(fit2)
    #topTable(fit2, coef=1, adjust="BH") #adjust="fdr", sort.by="B",
    d = topTable(fit2, coef=1, adjust=adjust.method, sort.by="B", number=50000)
    d$GeneID = rownames(d)
    n = ncol(d)
    d = d[,c(n, 1:(n-1))]
    if(!is.null(out.file)) {
	write.table(d, file=out.file, row.names=F, quote=F, sep="\t")
	#results = decideTests(fit2, adjust.method=adjust.method, p.value=cutoff) # none
	#png(paste(out.file, ".png", sep=""))
	#vennDiagram(results)
	#dev.off()
    }
    return(d) #(d[d$adj.P.Val<=cutoff,])
}


convert.probe.to.gene.expression<-function(expr, gene.mapping, selection.function=NULL) {
    if(nrow(gene.mapping) != nrow(expr)) {
	print("Warning: dimension inconsisitency in gene annotation!")
	gene.mapping = as.data.frame(gene.mapping[gene.mapping$Probe %in% rownames(expr),]) #gene.mapping[rownames(gene.mapping) %in% rownames(expr),])
	expr = expr[rownames(expr) %in% gene.mapping$Probe,] #expr[rownames(gene.mapping),]
    }
    gene.mapping = factor(gene.mapping[,"Gene"])
    if(is.null(selection.function)) {
	selection.function<-function(asample){ # max
	   #return(tapply(abs(asample), gene.mapping, max)) #mean)) 
       return(tapply(col.values, mapping, function(x) { x[which.max(abs(x))] }))
	}
	#variances = apply(expr, 1, var)
	#get.max<-function(e) {
	#    idx = which.max(e[,2])
	#    return(e[idx,1])
	#}
	#selection.function<-function(asample) { # max.variance
	#    return(by(data.frame(asample, variances), gene.mapping, get.max))
	#}
    }
    expr.gene<-apply(expr, 2, selection.function)
    expr.gene = expr.gene[rownames(expr.gene) != "", ] # Filter the one from no gene probes
    return(expr.gene)
}


get.data.set<-function(gds.id, output.dir) {
    data.set = NULL
    file.name = paste(output.dir, gds.id, ".annot.gz", sep="")
    file.name2 = paste(output.dir, gds.id, ".soft.gz", sep="")
    if(file.exists(file.name)) {
	data.set = getGEO(filename=file.name)
    } else if(file.exists(file.name2)) {
	data.set = getGEO(filename=file.name2)
    } else {
	file.name = paste(output.dir, gds.id, ".annot", sep="")
	file.name2 = paste(output.dir, gds.id, ".soft", sep="")
	if(file.exists(file.name)) {
	    data.set = getGEO(filename=file.name)
	} else if(file.exists(file.name2)) {
	    data.set = getGEO(filename=file.name2)
	} else {
	    # Get GDS file
	    data.set = getGEO(gds.id, destdir=output.dir, AnnotGPL=T)
	}
    }
    return(data.set)
}


get.platform.annotation<-function(data.set, probe.conversion, output.dir) {
    gds.id = Meta(data.set)$platform
    data.set = get.data.set(gds.id, output.dir)
    d = Table(data.set)
    print(colnames(d)) 
    #print(probe.conversion)
    print(length(d[,probe.conversion]))
    gene.mapping = data.frame(Probe = d[,"ID"], Gene = d[,probe.conversion]) #row.names = d[,"ID"],  "Gene.Symbol",
    return(gene.mapping)
    # Alternative using bioconductor annotation packages
    #source("http://www.bioconductor.org/biocLite.R")
    #biocLite("hgu133a")
    #mget("121_at",hgu133aSYMBOL)
    #mget("121_at",hgu133aUNIGENE)
}


### NOT VERIFIED 
expressionPreprocess<-function() {
    expr.file = paste(data.dir, "lung_expression.csv", sep="") # Manually assigned a header line containing sample name and call
    expr.gene.file = paste(data.dir, "lung_normalized.dat", sep="")
    expr<-read.table(expr.file, header=T, sep = "\t", quote = "\"", row.names=1, check.names = F, comment.char = "!")
    call.indices = seq(2, 44, by=2)
    filter.by.call<-function(x) {
        #if(all(x[call.indices] == "P" | x[call.indices] == "M")) {
        #    return(as.double(x[-call.indices]))
        #} 
        #return(rep(NA, length(call.indices)))
        missing.indices = x[call.indices] == "A" | x[call.indices] == "NC"
        val = as.double(x[-call.indices])
        val[missing.indices] = NA
        return(val)
    }
    e = apply(expr, 1, filter.by.call)
    e = t(e)
    colnames(e) = gsub("#", "", gsub(" ", ".", colnames(expr[-call.indices]))) # R compatible colnames
    expr = na.omit(e)
    #library(impute)
    #expr = impute.knn(e)$data
    #expr = expr[,-c(18,21,22)]
    # Normalize probe level expression data
    expr.norm<-normalizeArrays(expr)
    # Get probe to gene mapping (choose gene ids that are in the interactome)
    gene.mapping<-getGeneMapping(expr.norm)
    # Map probe level data to gene level data
    expr.gene<-convertProbetoGeneExpression(expr.norm, gene.mapping)
    i = which(rownames(expr.gene)=="NA")
    expr.gene = expr.gene[-i,]
    write.table(expr.gene, expr.gene.file, quote=F)
}

convertProbetoGeneExpression.MaxVar<-function(d, gene.mapping) {
    variances<-apply(d, 1, var)
    #indices<-tapply(variances, gene.mapping, which.max)
    getMax<-function(e) {
        #print(c(length(e), dim(e)))
        idx = which.max(e[,2])
        return(e[idx,1])
    }
    maxVariance<-function(asample) {
        return(by(data.frame(asample, variances), gene.mapping, getMax))
    }
    e<-apply(d, 2, maxVariance)
    png("after_mapping_max_var.png")
    boxplot(e)
    dev.off()
    return(e);
}

convertProbetoGeneExpression<-function(d, gene.mapping) {
    #gene.mapping<-factor(gene.mapping)
    averageProbesOfGenes<-function(asample){
       return(tapply(asample, gene.mapping, max)) # before it was mean
    }
    e<-apply(d, 2, averageProbesOfGenes)
    png("after_mapping.png")
    boxplot(e)
    dev.off()
    #rownames(e)<-gene.mapping
    return(e);
}

getGeneMapping<-function(expr) {
    annotation.file = paste(data.dir, "GPL96_trimmed.annot", sep="")
    network.file<-paste(data.dir, "../9606/network_no_tap_geneid.sif", sep="")
    annot<-read.table(annotation.file, header=T, sep="\t", row.names=1, check.names = F)
    netw<-read.table(network.file, sep=" ")
    geneids.network<-levels(factor(c(netw$V1,netw$V3)))
    #gene.mapping<-factor(annot[rownames(expr),"Gene ID"])
    gene.mapping<-rep(NA,nrow(expr))
    for(i in 1:nrow(expr)) {
        #print(as.vector(annot[rownames(expr)[i], "Gene ID"]))
        words = unlist(strsplit(as.vector(annot[rownames(expr)[i], "Gene ID"]), "///"))
        val = "NA"
        for(word in words) {
            if(is.element(word, geneids.network)) {
                val = word
                break;
            }
        }
        gene.mapping[i] = val
    }
    return(gene.mapping)
}

normalizeArrays<-function(d) {
    library(affy)
    library(affyPLM)
    d.colnames<-colnames(d)
    d.rownames<-rownames(d)
    e<-as.matrix(d)
    png("before_normalization.png")
    boxplot(d)
    dev.off()
    e<-rma.background.correct(e)
    e<-normalize.quantiles(e)
    e<-log2(e)
    colnames(e)<-d.colnames
    rownames(e)<-d.rownames
    png("after_normalization.png")
    boxplot(e)
    dev.off()
    #library(beadarray)
    #plotMAXY(e, arrays=1:3)
    return(e);
}

# Illimuna
normalizeArrays<-function(expr.file, expr.val.file, tf.file, merged.file, merged.norm.file, merged.norm.gene.file, merged.norm.tf.file) {
    library(beadarray)
    data<-readBeadSummaryData(expr.file, skip=0, columns=list(exprs = "AVG_Signal"))
    data.val<-readBeadSummaryData(expr.val.file, skip=0, columns=list(exprs = "AVG_Signal"))
    data.merged<-combine(data, data.val)
    d<-exprs(data)
    e<-exprs(data.val)
    data.merged<-data.merged[rownames(d)[rownames(d) %in% rownames(e)]]
    #merged<-merge(d, e, by.x = "row.names", by.y = "row.names", all = FALSE)
    #write.table(merged, merged.file)
    data.merged.norm <- normaliseIllumina(data.merged, method="quantile", transform="log2")
    merged.norm <- exprs(data.merged.norm)
    write.table(merged.norm, merged.norm.file)
    #plotMAXY(merged.norm, arrays=1:3)
    gene.levels<-factor(fData(data.merged.norm)$TargetID)
    averageProbesOfGenes<-function(asample){ 
       return(tapply(asample,gene.levels,mean)) 
    }
    merged.norm.gene<-apply(merged.norm, 2, averageProbesOfGenes)
    #write.table(merged.norm.gene, merged.norm.gene.file)
    tfs<-read.table(tf.file)
    merged.norm.tf<-merged.norm.gene[rownames(merged.norm.gene) %in% levels(tfs[,1]),]
    write.table(merged.norm.tf, merged.norm.tf.file)
}

