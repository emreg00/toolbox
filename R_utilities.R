###
# Copied pasted code from previous projects
# Need to refactor significantly to create indipendent reusable functions
###

library(ggplot2)
library(reshape)
library(scales)
library(RColorBrewer)
source("heatmap3.R")

cbPalette <- c("blue", "red", "green", "grey20", "orange") #c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cPalette <- brewer.pal(9,"Blues")

main<-function() {
}

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

draw.scatterplots<-function(file.name) {
    d = read.table(file.name, header=T)
    #var_size = "n.target"
    e = d
    e$z = ifelse(e$z > 6, 6, e$z)
    e$Category = ifelse(e$flag=="True", "Known", "Unknown")
    selected = "z"

    variable = "n.target"
    p = ggplot(data = e, aes_string(x = variable, y = selected)) + geom_point(aes_string(color="Category"), alpha=0.5) + theme_bw() + guides(color=guide_legend("Category"), alpha=F) + labs(x="Number of drug targets", y="Normalized distance (z)") + scale_color_manual(values=cbPalette[1:2])  # size = var_size, size=guide_legend("N target"), 
    p = p + theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank())
    p = p + theme(axis.ticks=element_blank(), axis.line = element_line(color = 'black'))
    out.file = paste(img.dir, "scatter_", sub("\\.", "_", variable), ".svg", sep="")
    svg(out.file) #, width=10, height=6)
    print(p)
    dev.off()
    a = cor.test(e[,variable], e[,selected], method="spearman")
    print(c(a$estimate, a$p.value))

    variable = "n.disease"
    p = ggplot(data = e, aes_string(x = variable, y = selected)) + geom_point(aes_string(color="Category"), alpha=0.5) + theme_bw() + guides(color=guide_legend("Category"), alpha=F) + labs(x="Number of disease genes", y="Normalized distance (z)") + scale_color_manual(values=cbPalette[1:2])  
    p = p + theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank())
    p = p + theme(axis.ticks=element_blank(), axis.line = element_line(color = 'black'))
    out.file = paste(img.dir, "scatter_", sub("\\.", "_", variable), ".svg", sep="")
    svg(out.file) 
    print(p)
    dev.off()
    a = cor.test(e[,variable], e[,selected], method="spearman")
    print(c(a$estimate, a$p.value))
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

disease.drug.barplots<-function(file.name) {

    method = "z" #"z.tissue"
    cutoff = -0.84 #-2

    d=read.table(paste(file.name, sep=""), header=T)
    e = d[!is.na(d[[method]]),]
    e = e[e$flag=="True",]
    a = split(e$group, e$disease)
    y = sapply(a, length)
    e = e[e[[method]]<=cutoff,]
    b = split(e$group, e$disease)
    x = sapply(b, length)

    out.file = paste(img.dir, "disease_drug_barplot_ndfrt_m2m_100000_subset.svg", sep="") 
    #out.file = paste(img.dir, "disease_drug_barplot_mesh299_ndfrt_m2m_subset.svg", sep="") # mesh299
    svg(out.file, width=8, height=6)
    #par(family = "Arial")
    par(mar=c(12, 4, 4, 5) + 0.1)
    n.max = max(unlist(lapply(a, length)))
    n.max = ceiling(n.max / 10) * 10
    a<-barplot(rbind(x, y), beside=T, ylab="Number of drugs", xaxt="n", width=1, ylim=c(0,n.max))
    x.labels = lapply(names(x), function(x){ return(gsub("[.]", " ", x)) })
    print(x.labels)
    text(colMeans(a),  par("usr")[3], labels=x.labels, adj=c(1, 1), srt=45, cex=0.9, xpd=T) # mesh69
    #text(colMeans(a),  par("usr")[3], labels=x.labels, adj=c(1, 0.5), srt=90, cex=0.4, xpd=T) # mesh299
    par(new=T)
    plot(colMeans(a), 100*x/y, col=2, xaxt="n", yaxt="n", xlab="", ylab="", type='p', pch=4, bty="n", ylim=c(0,100))
    axis(4, xpd=T, col=2, col.axis=2)
    mtext("Ratio (%)", side=4, line=3, col=2)
    legend("topright", c("Close", "All", "Ratio"), pch=c(15, 15, 4), pt.cex=2, col=c("grey30", "grey70", 2), lty=c(0,0,0), bty="n")
    dev.off()
}

enrichment.incompleteness.plots<-function(file.name) {
    method = "z" 
    cutoff = -2 
    i = 1
    f = data.frame()
    for(suffix in c("HI4/mesh69/values/drug-individual/seeds/closest_None.dat", "STRING9/mesh69/values/drug-individual/seeds/closest_None.dat", "PPI2012/mesh69/values/drug-individual/seeds/closest_None_random50.dat", "PPI2012/mesh69/values/drug-individual/seeds/closest_None.dat", "PPI2012/mesh69/values/drug-individual/seeds/closest_None_random25.dat", "PPI2012/mesh69/values/drug-individual/seeds/closest_None_random10.dat")) {
	d=read.table(paste(file.name, suffix, sep=""), header=T)
	if(i < 4) {
	    e = d[!is.na(d[[method]]),]
	    x=e[e$flag=="True", ]
	    y=e[e$flag!="True", ]
	    a = sum(x[[method]]<=cutoff, na.rm=T)
	    b = sum(y[[method]]<=cutoff, na.rm=T) 
	    contigency = matrix(c(a, nrow(x)-a, b, nrow(y)-b), 2, 2)
	    s = fisher.test(contigency)
	    or = data.frame(value=s$estimate, variable=i)
	    if(i == 1) {
		p = ggplot(or, aes(x=variable, y=value)) + geom_point(aes(color=3))
	    } else {
		p = p + geom_point(data=or, aes(color=3))
	    }
	} else {
	    m = mean(d$or)
	    s = sd(d$or)
	    or = data.frame(value=m, variable=i)
	    p = p + geom_point(data=or, aes(color=3)) 
	    f[i-3,"value"] = m
	    f[i-3,"variable"] = i
	    f[i-3,"s"] = s
	}
	print(or)
	i = i + 1
    }
    p = p + geom_errorbar(data=f, aes(x=variable, ymax = value + s, ymin = value - s), width=0.4)  
    p = p + labs(x="Interaction network", y="Odds Ratio") + theme_bw() + guides(color=F) + xlim("Y2H", "Functional", "PPI", "PPI_50", "PPI_25", "PPI_10")
    out.file = paste(img.dir, "OR_comparison.svg", sep="")
    svg(out.file, width = 6, height = 6, onefile = TRUE)
    print(p)
    dev.off()
    return();
}


drug.disease.heatmap<-function(diseases, suffix) {
    out.file = paste(img.dir, "distance_heatmap.svg", suffix, sep="")
    file.name = paste(data.dir, prefix, ".values", suffix, sep="")
    e = read.table(file.name, header=T, sep="\t", check.names=F)
    e = e[,diseases]
    file.name = paste(data.dir, prefix, ".rowmapping", suffix, sep="")
    f = read.table(file.name, header=T, sep="\t", check.names=F)
    f = f[,diseases]
    indices = rowSums(f) > 0
    e = e[indices,]
    f = f[indices,]
    #e = d[d$flag == "True",]
    #drugs = split(e$group, e$disease)
    #n.drug=tapply(e$group, factor(e$disease), function(x){length(levels(factor(x)))})

    # To order f w.r.t. first disease then name
    indices = order(f[,2],f[,1],rownames(f))
    # Order labels alphabetically
    #indices = order(rownames(e))
    f = f[indices,]
    e = e[indices,]

    val.cols = brewer.pal(7,"RdBu") #bluered(75)
    cols = rep(c("darkolivegreen3", "darkorchid3", "wheat3", "aquamarine3", "yellow3", "orchid3", "cadetblue3", "chocolate3", "green", "yellow", "orange3", "slateblue3"), 3) #  # brewer.pal(8,"Dark2") 
    n = ncol(f)
    g = c()
    ##column.cols = c()
    for(i in 1:n) {
	g<-cbind(g, ifelse(f[,i]==1, cols[i], "lightgray"))
	##column.cols<-c(column.cols, rep(cols[i], d[colnames(f)[i], "n"]+1))
    }
    g<-g[,n:1]
    colnames(g) = colnames(f)
    row.cols = t(as.matrix(g))
    ##column.cols = cbind(cols, cols)
    ##colnames(column.cols) = rep(NA,2)
    column.cols = as.matrix(cols[1:n])
    colnames(column.cols) = c("Disease") #NA)
    
    # For d (distance) symbreaks=F, symkey=F
    symetry = T
    if(suffix == ".d") {
	symetry = F
    }
    svg(out.file, width = 6, height = 8, onefile = TRUE)
    heatmap.3(e, labCol=colnames(e), dendrogram="none", Rowv=F, Colv=F, revC=F, scale="none", na.color="grey", col=val.cols, RowSideColors=row.cols, ColSideColors=column.cols, margins=c(6,20), cexRow=0.8, cexCol=0.8, trace="none", xlab="Drug-disease closeness (z)", NumRowSideColors=0.5, NumColSideColors=0.5, symbreaks=symetry, key=T, symkey=symetry, density.info="none") # KeyValueName="Log expression")
    dev.off()
}

jaccard.heatmaps<-function() {
    out.file = paste(img.dir, "drug_jaccard_heatmap.svg", sep="")
    file.name = paste(output.dir, "drug_jaccard.dat", sep="")
    d=read.table(file.name, header=T, sep="\t")
    p = ggplot(d, aes(x=x, y=y)) + geom_tile(aes(fill=value)) + coord_equal() + scale_fill_gradient(low="white", high="steelblue", na.value = "lightgrey") + theme_bw()
    p = p + theme(axis.text.y = element_text(size = 6), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) + xlab(NULL) + ylab("Drugs")
    svg(out.file, width = 6, height = 6, onefile = TRUE)
    print(p)
    dev.off()

    for(prefix in c("disease_jaccard", "seed_jaccard", "disease_drug_jaccard")) { 
	out.file = paste(img.dir, prefix, "_heatmap.svg", sep="")
	file.name = paste(output.dir, prefix, ".dat", sep="")
	d=read.table(file.name, header=T, sep="\t")
	p = ggplot(d, aes(x=x, y=y)) + geom_tile(aes(fill=value)) + coord_equal() + scale_fill_gradient(low="white", high="steelblue", na.value = "lightgrey") + theme_bw()
	p = p + theme(axis.text.y = element_text(size = 8), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) + xlab(NULL) + ylab("Diseases")
	svg(out.file, width = 6, height = 6, onefile = TRUE)
	print(p)
	dev.off()
    }

    # n = length(levels(d$x))
    # a = matrix(nrow=n, ncol=n)
    # rownames(a)=levels(d$x)
    # colnames(a)=levels(d$y)
    # for(i in 1:nrow(d)) { a[d[i,]$x, d[i,]$y] = d[i,"value"] }
    # a[upper.tri(a)] = NA
    # diag(a) = NA
    # heatmap(a, Rowv=NA, Colv=NA, scale="none",revC=T, col=brewer.pal(7,"Blues"))

    out.file = paste(img.dir, "disease_drug_heatmap.svg", sep="")
    file.name = paste(output.dir, "disease_drug.dat", sep="")
    d=read.table(file.name, header=T, sep="\t")
    p = ggplot(d, aes(x=y, y=x)) + geom_tile(aes(fill=3)) + guides(fill=F) + theme_bw()   
    p = p + theme(axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 5, angle=90, vjust=0.5, hjust=0.5), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank()) + xlab("Drugs") + ylab("Diseases")
    svg(out.file, width = 6, height = 6, onefile = TRUE)
    print(p)
    dev.off()

    for(suffix in c("biogps", "ictnet", "lage")) {
	out.file = paste(img.dir, "tissue_info_", suffix, ".svg", sep="")
	file.name = paste(output.dir, "tissue_info_", suffix, ".dat", sep="")
	d=read.table(file.name, header=T, sep="\t")
	p = ggplot(d, aes(x=y, y=x)) + geom_tile(aes(fill=3)) + guides(fill=F) + theme_bw()
	p = p + theme(axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8, angle=90, vjust=0.5, hjust=0.5), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank()) + xlab("Tissues") + ylab("Diseases")
	svg(out.file, width = 6, height = 6, onefile = TRUE)
	print(p)
	dev.off()
    }
}

draw.separation.barplot<-function(file.name, out.file, selected=NULL) {
    d = read.table(file.name, sep=" ", header=T)
    if(is.null(selected)) {
	selected = colnames(d)[-(1:2)]
    }
    if(length(selected) == 3) {
	e = d
	e$variable = e[[selected[2]]]
	e$value = e[[selected[3]]]
	#x.seq = seq(0.5,3.5,0.75)
	#e = transform(e, variable = x.seq)
	e <- transform(e, variable = reorder(variable, 1:6)) # ggplot keeps the order in the data frame c(2,1,3,4,5,6))) 
	p = ggplot(data = e, aes(x = variable, y = value)) + geom_bar(aes(fill=method), width = 0.6, position="dodge", stat="identity", alpha=0.5) # + scale_x_discrete(breaks = NA)
	p = p + scale_fill_manual(values=cbPalette[1:2]) + coord_cartesian(ylim=c(0.4,0.8)) # ylim(c(0.5, 1)): drops values below the limit
	p = p + theme_bw() + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=0.5)) + labs(x=NULL, y="AUC (%)") + guides(fill=F) #guide_legend("Negative selection")) # size=7, 
	p = p + theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank())
	p = p + theme(axis.ticks=element_blank(), axis.line = element_line(color = 'black'))
    } 
}


compare.values<-function(file.name, out.file, alternative=NULL, xlabel="") {
    d=read.table(file.name, header=T)
    x=d[d$variable=="True","value"]
    y=d[d$variable=="False","value"]
    print(c(length(x), length(y)))
    a=wilcox.test(x, y, alternative=alternative)
    print(c(mean(x), mean(y)))
    print(c(median(x), median(y), a$p.value))
    p = ggplot(d, aes(value, fill = variable)) + geom_histogram(alpha = 0.4, aes(y = ..density..), position = 'identity') #, binwidth=1)
    p = p + labs(x=xlabel, y="Density")
    #p = p + geom_density(aes(color = variable), alpha=0.0)
    p = p + theme_bw()
    svg(out.file, width=10, height=6) 
    print(p)
    dev.off()
}

histogram.drug.info<-function() {
    file.name = paste(output.dir, "target.dat", sep="") 
    e = read.table(file.name, header=T)
    n = max(e$n.degree)
    p = ggplot(e, aes(n.degree)) + geom_histogram(fill="blue", color="black", alpha = 0.5, binwidth=1) + theme_bw() 
    p = p + ylab("Number of drugs") + xlab("Average degree of targets") + scale_x_discrete(breaks=seq(0,n,by=50)) 
    p = p + theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank())
    p = p + theme(axis.ticks=element_blank(), axis.line = element_line(color = 'black')) #, axis.text.x=element_text(seq(0,n,by=50)))
    print(summary(e$n.degree))

    out.file = paste(img.dir, "degree.svg", sep="")
    svg(out.file, width = 8, height = 6, onefile = TRUE)
    print(p)
    dev.off()
}

get.ml.model<-function(file.name) {

    d=read.table(file.name, header=T)

    d$class = ifelse(d$flag=="True", 1, 0)
    fit <- glm(class ~ z + z.disease + z.group + m.p.z, data = d, family = "binomial")

    library(rpart)
    x=d[d$flag=="True", ]
    y=d[d$flag=="False", ]
    z=y[sample(nrow(y))[1:nrow(x)],]
    e=as.data.frame(rbind(x,z))
    #fit <- rpart(class ~ z + z.disease + z.group, method="class", data=e)
    #fit <- rpart(class ~ z + z.disease + z.group + m.e + m.n + p.m.n + pval, method="class", data=e)
    fit <- rpart(class ~ z + z.disease + z.group + m.e + m.n.j + m.p.j + p.m.n + n + m + d.m + d.n + pval, method="class", data=d)

    #printcp(fit) 
    #plotcp(fit) 
    #summary(fit) 

    #plot(fit, uniform=TRUE, main="Classification Tree for Kyphosis")
    #text(fit, use.n=TRUE, all=TRUE, cex=.8)

    pfit<- prune(fit, cp=fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])

    printcp(pfit) 
    plotcp(pfit) 
    summary(pfit) 

    # plot the pruned tree
    plot(pfit, uniform=TRUE, main="Pruned Classification Tree for Kyphosis")
    text(pfit, use.n=TRUE, all=TRUE, cex=.8)

    library(randomForest)
    e=na.omit(e)
    #fit <- randomForest(class ~ z + z.disease + z.group + m.e + m.n.j + m.p.j + p.m.n + d.m + d.n + pval, method="class", data=e)
    fit <- randomForest(class ~ z + z.disease + z.group + m.e + m.n.j + m.p.j + p.m.n + n + m + d.m + d.n + pval, method="class", data=e)
    print(fit)
    importance(fit)
}

draw.inflammation.heatmaps<-function() {
    suffices = c("closest") #, "shortest", "kernel", "center", "jorg.individual")
    pdf(paste(img.dir, "inflammation_z.pdf", sep=""))
    for(suffix in suffices) {
	file.name = paste(output.dir, "PPI2011/inflammation/values/disease-disease/seeds/", suffix, "_None.dat", sep="")
	f = read.table(file.name, header=T)
	p = ggplot(data=f, aes(group, disease)) + geom_tile(aes(fill=z)) + labs(x=NULL, y=NULL, title=suffix) 
	p = p + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) 
	#p = p + scale_fill_gradient(low = "red", high = "grey")
	#p = p + scale_fill_gradient(low = "steelblue", high = "white") 
	p = p + scale_fill_gradient(low = "steelblue", high = "red")
	#print(p)
	p = ggplot(data=f, aes(group, disease)) + geom_point(aes(color=z), size=10) + labs(x=NULL, y=NULL, title=suffix) 
	p = p + theme_bw() + theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank())
	p = p + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) 
	p = p + scale_color_gradient(low = "steelblue", high = "red")
	print(p)
    }
    dev.off()
}




get.classifier.plots.helper<-function(file.name, suffix, suffix2, outfix) {
    file.name = paste(file.name, suffix, ".classifier", suffix2, sep="")
    f = read.table(file.name)
    # Note that z's are inverted (-z)
    #f$V1 = -f$V1 #, label.ordering=c(1,0)) this would not help
    pred = prediction(f$V1, f$V2) 
    perf = performance(pred, "auc")
    print(sprintf("%s %.3f", "auc", perf@y.values[[1]]))

    perf = performance(pred, "tpr", "fpr")
    out.file = paste(img.dir, outfix, "_roc", suffix, suffix2, ".svg", sep="")
    #svg(out.file) #!
    plot(perf)
    #dev.off()

    perf = performance(pred, "sens", "spec")
    #a = perf@x.values[[1]] + perf@y.values[[1]]
    #i = which.max(a)
    #print(perf@alpha.values[[1]][i])

    out.file = paste(img.dir, outfix, suffix, suffix2, ".svg", sep="")
    #svg(out.file) #!
    perf = performance(pred, "sens")
    a = perf@y.values[[1]] 
    plot(perf, col="blue", bty="n", xlab="Normalized distance (z)", ylab="Performance") # -perf@x.values[[1]], 100*perf@y.values[[1]]
    perf = performance(pred, "spec")
    plot(perf, col="red", bty="n", add=T)
    legend("right", legend=c("Sensitivity", "Specificity"), col=c("blue","red"), lty=c(1,1), bty="n")
    b = perf@y.values[[1]]
    #dev.off()

    i = which.max(a + b)
    #i = which(a == b)
    #i = which.max(abs(a - b) <= 0.005)
    print(sprintf("%s %.3f", "cutoff", perf@x.values[[1]][i]))
    print(sprintf("%s %.3f %s %.3f", "sens", a[i], "spec", b[i]))

    #i = which(a>=0.326)
    #print(max(perf@x.values[[1]][i]))

    #i = which(a>=0.626)
    #print(max(perf@x.values[[1]][i]))
}


visualize.array.correlation.grouped.by.sample.type <- function(expr, sample.mapping, states, out.file) {
    groups = split(sample.mapping$Sample, sample.mapping$Type)
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
    val.cols <- brewer.pal(9,"Blues") # greenred
    svg(out.file) 
    heatmap(e, Rowv=NA, Colv=NA, col=val.cols, scale="none", revC=T) #labCol=F, labRow=F, , margins=c(1,1), keysize=0.9, colsep=seperators, rowsep=seperators, sepcolor="white")
    dev.off()
}


find.differentially.expressed.genes<-function(expr, sample.mapping, states, out.file) {
    library(limma)
    design = model.matrix(~ 0 + sample.mapping$Type)
    colnames(design) = gsub(" ", "_", states) #colnames(design))
    fit = lmFit(expr, design)
    #contrast = unlist(sapply(states, function (x) { if(ref != x) { paste(ref, x, sep = "-") } }))
    contrast = c()
    for(i in 1:length(states)) { 
	for(j in 1:length(states)) { 
	    if(i<j) { 
		contrast = c(contrast, paste(gsub(" ", "_", states[i]), gsub(" ", "_", states[j]), sep = "-"))
	    }
	}
    }
    cont.matrix = makeContrasts(contrasts=contrast, levels=design)
    fit2 = contrasts.fit(fit, cont.matrix)
    fit2 = eBayes(fit2)
    topTable(fit2, coef=1, adjust="BH") #adjust="fdr", sort.by="B",
    write.table(topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=50000), file=out.file, row.names=F, sep="\t") 
    results = decideTests(fit2, p.value=0.05)
    vennDiagram(results)
    return(fit2)
}


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


main()

