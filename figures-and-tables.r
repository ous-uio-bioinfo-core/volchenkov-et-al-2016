


### Load libs, make outputfolders
library(limma)
library(RColorBrewer)
library(gplots)

figuredir = paste("figures-and-tables", sep="")
if(!file.exists(figuredir))dir.create(figuredir)

outputdir = paste("auxiliary-results", sep="")
if(!file.exists(outputdir))dir.create(outputdir)
genelistdir = paste(outputdir, "/genelists", sep="")
if(!file.exists(genelistdir))dir.create(genelistdir)
plotdir= paste(outputdir, "/plots", sep="")
if(!file.exists(plotdir))dir.create(plotdir)


### Read data and sampleannotation. Should maybe be done directly from GEO rep when availible?

# The datafile from GenomeStudio
datafil <- "not_in_github/SampleGeneProfileQuantileNormalizedImputed.txt"
dataset = read.ilmn(files=datafil , other.columns="Detection", probeid="TargetID$")

# Read sampleannot
safile="sampleannotation.xls"
sa <- read.table(file=safile, header=TRUE, sep="\t", stringsAsFactors=TRUE)
# table(sa$sampleid %in% colnames(dataset)) # all there
# table(colnames(dataset) %in% sa$sampleid )

sa = sa[match(colnames(dataset), sa$sampleid ), ]
# make some combo groups
sa$combocell = factor(paste(sa$celltype,sa$addition, sep=""))
sa$combogroup = factor(paste(sa$oxygen, sa$combocell,  sep="_"))
sa$comboignoreCSE = factor(paste(sa$oxygen, sa$celltype,  sep="_"))

# library(GEOquery)  
# matrixeset=getGEO("GSE90882") # did not work. Repo was not accessible as long as it is private



### Find DEGs between a lot of groups using LIMMA 
ds=log2(dataset$E)

# The "diff" contrasts are interaction (using first method from limma guide)
# (B_ThCSE - FS_ThCSE) - (B_beads - FS_beads) #  answers this question:
# Which genes respond to hypoxia(B) differently in ThCSE compared to bead cells

toptables = list()
fac=factor(sa$combogroup)
design = model.matrix(~0 + fac)
colnames(design)=levels(fac)
cont.matrix = makeContrasts ( 
	"B_beads - FS_beads",
	"B_ThCSE - FS_ThCSE", 
	"B_Th - FS_Th",
	diff_ThCSE_beads="(B_ThCSE - FS_ThCSE) - (B_beads - FS_beads)",
	diff_Th_beads="(B_Th - FS_Th) - (B_beads - FS_beads)",
	diff_ThCSE_Th="(B_ThCSE - FS_ThCSE) - (B_Th - FS_Th)",
	"FS_Th - FS_ThCSE",
	"B_Th - B_ThCSE",
	"FS_beads - FS_ThCSE", 
	"B_beads - B_ThCSE",
	"B_beads - B_Th", 
	"FS_Th - FS_beads",
	levels=design)
# not all contrasts are reported in the figures later!

fit <- lmFit(ds, design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2) 
fit2[["genes"]] <- data.frame( geneSymbol=rownames(fit))

for(t in 1:dim(cont.matrix)[2])
{
	this_coef=dimnames(cont.matrix)[2][[1]][t]
	coefname = gsub(" - ", "_vs_", this_coef)
	tab = topTable(fit2, adjust='fdr', coef=this_coef , number=dim(ds)[1])	
	write.table(tab, file=paste(genelistdir, "/", coefname, ".txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE)
	print(paste( "Printing ", this_coef, ", FDR<0.05 ", sum(tab$adj.P.Val<0.05)   ))
	toptables[[this_coef]]=tab	
}

# selected groupings to be counted in venn diagrams
vennnames = list(
	effectspergroup = c("B_beads - FS_beads", "B_Th - FS_Th", "B_ThCSE - FS_ThCSE"),
	diff_Th_beads = c("B_Th - FS_Th", "B_beads - FS_beads", "diff_Th_beads"),
	diff_ThCSE_beads = c("B_ThCSE - FS_ThCSE", "B_beads - FS_beads", "diff_ThCSE_beads"),
	diff_ThCSE_Th = c("B_ThCSE - FS_ThCSE", "B_Th - FS_Th", "diff_ThCSE_Th")
)

# Figure not included in manuscript
pdf(file = paste(plotdir, "/venndiffcount.pdf", sep=""), width = 12, height = 12)
par(mfrow=c(2,2))
venninclude = "both"
results <- decideTests(fit2)
vennDiagram(results[,vennnames[[1]]], include=venninclude, cex=1)
vennDiagram(results[,vennnames[[2]]], include=venninclude, cex=1)
vennDiagram(results[,vennnames[[3]]], include=venninclude, cex=1)
vennDiagram(results[,vennnames[[4]]], include=venninclude, cex=1)
dev.off()

# Figure 4b
pdf(file = paste(figuredir, "/figure4b.pdf", sep=""), width = 12, height = 12)
results <- decideTests(fit2)
vennDiagram(results[,c("B_beads - FS_beads", "B_Th - FS_Th", "B_ThCSE - FS_ThCSE")], include="both", cex=1)
dev.off()


# the effect of oxygen
contrasts = c("oxygenFS - oxygenB")
design = model.matrix(~0 + oxygen + combocell, data=sa)
cont.matrix = makeContrasts ( contrasts=contrasts, levels=design)
fit <- lmFit(ds, design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2) 
fit2[["genes"]] <- data.frame( geneSymbol=rownames(fit))
coefname = "FS_vs_B"
this_coef=dimnames(cont.matrix)[2][[1]][1]
tab = topTable(fit2, adjust='fdr', coef=this_coef , number=dim(ds)[1])	
write.table(tab, file=paste(genelistdir, "/", coefname, ".txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE)
print(paste( "Printing ", this_coef, ", FDR<0.05 ", sum(tab$adj.P.Val<0.05)   ))
toptables[[this_coef]]=tab	
results =data.frame(results, decideTests(fit2))


### combining Th and ThCSE and look for differences in response to hypoxia
# assmuning Th and ThCSE is very similar and combining them will increase sample size
design = model.matrix(~0 +comboignoreCSE , data=sa)
colnames(design)=levels(sa$comboignoreCSE)
cont.matrix = makeContrasts ( 
	hyp_ThignoreCSE="B_Th - FS_Th",
	diff_ThignoreCSE_beads="(B_Th - FS_Th) - (B_beads - FS_beads)",	
	levels=design)

fit <- lmFit(ds, design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2) 

fit2[["genes"]] <- data.frame( geneSymbol=rownames(fit))
for(t in 1:dim(cont.matrix)[2])
{
	this_coef=dimnames(cont.matrix)[2][[1]][t]
	coefname = gsub(" - ", "_vs_", this_coef)
	tab = topTable(fit2, adjust='fdr', coef=this_coef , number=dim(ds)[1])	
	write.table(tab, file=paste(genelistdir, "/", coefname, ".txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE)
	print(paste( "Printing ", this_coef, ", FDR<0.05 ", sum(tab$adj.P.Val<0.05)   ))
	toptables[[this_coef]]=tab	
}
results =data.frame(results, decideTests(fit2))
## no genes found so maybe Th and ThCSE are not som similar.


summarytab = data.frame( comparison=names(results), up=colSums(results==1), down=colSums(results==-1))
write.csv(summarytab, file=paste(genelistdir, "/", "summarycounts.csv", sep=""), row.names=FALSE)


####### some heatmaps for selected genes
heatmapwrapper = function(dsclust, salabels, main="", plotlegend=TRUE, plotlabRow=TRUE, dendrogram = "both", inset=0)
{
	labRow=rownames(dsclust)
	if(plotlabRow==FALSE)
	{
		labRow=NA
	}	
	colpalette = c("green", "orange", "blue", "brown", "yellow", "magenta")
	heatmap.2( dsclust, scale="row", ColSideColors=colpalette[salabels], trace="none",
						 labRow=labRow, col=bluered(25) ,
						 margins=c(4,8), key=TRUE, main=main, dendrogram=dendrogram)
	#col=greenred(25)[-c(10,12,14,16)]
	#margins=c(8,1),
	if(plotlegend)
	{
		legend("topright",  xpd=TRUE,    
					 #inset=c(-0.08,-0.15),  # location of the legend on the heatmap plot
					 inset=inset,  # location of the legend on the heatmap plot
					 legend = unique(salabels), # category labels
					 fill = colpalette[unique(salabels)],  # color key
					 bty = "n"
		)		
	}	
}



# Figure 4a
pdf(paste(figuredir, "/figure4a.pdf", sep=""))
comp="B_Th - FS_Th"
genes = as.character(toptables[[comp]][ toptables[[comp]]$"adj.P.Val"<0.05 ,"geneSymbol"])
salabels = sa$combogroup
heatmapwrapper(ds[genes,],
							 sa$combogroup,
							 main=paste("DEGs ", comp, sep=""), inset=c(-0.08,-0.15))

dev.off()



# Heatmaps not included in manuscript
heatmapcomps = c("diff_Th_beads", "FS_Th - FS_ThCSE", "B_Th - B_ThCSE", 
								 "FS_beads - FS_ThCSE", "B_beads - B_ThCSE", "FS_Th - FS_beads",
								 "B_Th - FS_Th")
for(comp in heatmapcomps)
{
	
	pdf(paste(plotdir, "/heatmap_", comp, ".pdf", sep=""))
	genes = as.character(toptables[[comp]][ toptables[[comp]]$"adj.P.Val"<0.05 ,"geneSymbol"])
	salabels = sa$combogroup
	heatmapwrapper(ds[genes,],
								 sa$combogroup,
								 main=paste("DEGs ", comp, sep=""), inset=c(-0.08,-0.15))
	dev.off()
}

# Heatmap not included in manuscript
pdf(paste(plotdir, "/heatmap_25percvariablegenes.pdf", sep=""))
salabels = sa$combogroup
genevar = apply(ds, MARGIN=1, var)
dsclust = ds[order(genevar, decreasing=TRUE),]
colnames(dsclust) = paste("s_", colnames(dsclust), sep="")
sel = 1:(round(nrow(dsclust)/4))
heatmapwrapper(dsclust[sel,],
							 sa$combogroup,
							 main=paste("25% most variable genes", sep=""), dendrogram="column", plotlabRow=FALSE, inset=c(-0.08,-0.15))

dev.off()

# Figure 5
metabolicgenes = unique(scan("metabolicgenes.txt", what="character"))
pdf(paste(figuredir, "/figure5.pdf", sep=""))
heatmapwrapper(ds[metabolicgenes,],
							 sa$combogroup,
							 main=paste("Genes involved in metabolism", sep=""), dendrogram="both", plotlabRow=TRUE, inset=c(-0.08,-0.15))
dev.off()

# Figure 2
# Hierarchical cluster with heatmap to visualize genes that differ between the three FS groups (pairwise)
FScomps= c("FS_beads - FS_ThCSE", "FS_Th - FS_ThCSE", "FS_Th - FS_beads")
FS_DEG=c()
for(comp in FScomps)
{
	a=toptables[[comp]]$adj.P.Val < 0.05085  #NB odd cutoff, why was this used?
	FS_DEG = c(FS_DEG, as.character(toptables[[comp]]$geneSymbol[a]))
}

a = sa$combogroup %in% c("FS_beads", "FS_Th", "FS_ThCSE")
pdf(paste(figuredir, "/figure2.pdf", sep=""), height = 12, width=7)
heatmapwrapper(ds[FS_DEG,a],
							 sa$combogroup[a],
							 main=paste("DEG between 21% O2", sep=""), dendrogram="both", plotlabRow=TRUE,
							 inset=c(-0.03, -0.04))
dev.off()


# data file for GEO
#x = data.frame(ID_REF=rownames(dataset$E), dataset$E, check.names = FALSE)
#write.csv(x, file="not_in_github/geosubmission/matrix_table.csv", row.names = FALSE, quote=FALSE)


sink("sessionInfo.txt")
print(sessionInfo())
sink()