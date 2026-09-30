# Code to generate Supplemental Figure 5 and Supplemental Figure 6
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(data.table)
library(stringr)
library(RColorBrewer)

setwd("/mnt/d/Research/WhitesonCollab/Data/")

samplelist = c('C1I','C1M','C2I','C2M','C3I','C3M','S1I','S1M','S2I','S2M','S3I','S3M')
timelist = c('t1','t2','t3','t4','t5','t10','t15','t20')

calc_ht_size<-function(ht, unit='inch'){
	pdf(NULL)
	ht = draw(ht)
	w = ComplexHeatmap:::width(ht)
	w = convertX(w, unit, valueOnly = TRUE)
	h = ComplexHeatmap:::height(ht)
	h = convertY(h, unit, valueOnly = TRUE)
	dev.off()
	c(w, h)
}

darkpal = brewer.pal(8, 'Dark2')

allelelist = list()
for(i in 1:length(isomuts)){
	allelelist[[names(isomuts)[i]]] = paste(isomuts[[i]]$V1, isomuts[[i]]$V2, isomuts[[i]]$V4, isomuts[[i]]$V5, sep="_")
}

unialleles = unique(unlist(allelelist))
allelemat = array(0, dim=c(length(unialleles),length(isomuts)))
rownames(allelemat) = unialleles
colnames(allelemat) = names(isomuts)
for(i in 1:length(isomuts)){
	allelemat[allelelist[[i]],i] = 1
}

#Now, iterate through each endpoint isolate and each allele, find the earliest time that the allele appeared.
#Matrix setup: rows = alleles. Columns =  population

#Updated the following to include results from the endpoints in the same matrix
isopopmat = array(0,dim=c(length(unialleles),96))
rownames(isopopmat) = unialleles
colnames(isopopmat) = names(vcflist_trimmed)

#rebuild mutlist
mutlist_trimmed = list()
for(i in 1:length(vcflist_trimmed)){
	mutlist_trimmed[[i]] = paste(vcflist_trimmed[[i]]$V1,vcflist_trimmed[[i]]$V2,vcflist_trimmed[[i]]$V4,vcflist_trimmed[[i]]$V5, sep='_')
}
names(mutlist_trimmed) = names(vcflist_trimmed)
allvcf_trimmed=mutlist_trimmed
for(i in 1:length(unialleles)){
	for(j in 1:length(allvcf_trimmed)){
		if(unialleles[i]%in%allvcf_trimmed[[j]]){isopopmat[i,j]=1}
	}
}

#Add the endpoints in
combomat = cbind(isopopmat, allelemat)
#manually add in combined columns for S1M, S1I, S2I, S3M, S3I
combomat = cbind(combomat, combomat[,'S1M_bob']+combomat[,'S1M_car'])
colnames(combomat)[ncol(combomat)] = 'S1M'
combomat = cbind(combomat, combomat[,'S1I_bob']+combomat[,'S1I_car'])
colnames(combomat)[ncol(combomat)] = 'S1I'
combomat = cbind(combomat, combomat[,'S2I_bob']+combomat[,'S2I_car'])
colnames(combomat)[ncol(combomat)] = 'S2I'
combomat = cbind(combomat, combomat[,'S3M_bob']+combomat[,'S3M_car'])
colnames(combomat)[ncol(combomat)] = 'S3M'
combomat = cbind(combomat, combomat[,'S3I_bob']+combomat[,'S3I_car'])
colnames(combomat)[ncol(combomat)] = 'S3I'
#Then remove the original separated columns
combomat = combomat[,-c(103:108,110:113)]
#And correct the order for the sisters columns at the end
combomat = combomat[,c(1:102,105,104,106,103,108,107)]


# Load host data
hostdata = read.csv("./HostRange/AppelmansTS_ScoresAncesctal_and_Endpoints_IsoandCocktail_vs_TS_20_replicates.csv", header=T)
hostdata = hostdata[1:4,]
rownames(hostdata) = hostdata[,1]
hostdata = hostdata[,-1]
colnames(hostdata)[1:12] = paste0(str_extract(colnames(hostdata)[1:12],'[CS][123][IM]'))
colnames(hostdata)[13:24] = paste0(str_extract(colnames(hostdata)[13:24],'[CS][123][IM]'),'_t20')

combomat_subset = combomat[,c(samplelist, paste0(samplelist,'_t20'))]
#remove mutations in t0muts
combomat_subset = combomat_subset[-which(rownames(combomat_subset)%in%t0muts_trimmed),]

#condense repetitive rows to make it easier to read
#convert each row to a character vector. then use table to find rows that repeat more than 3 times.

comborows = character()
for(i in 1:nrow(combomat_subset)){
	comborows[i] = paste(combomat_subset[i,], collapse="")
}

combotab = table(comborows)
keep = numeric()
namelist = names(combotab)
for(i in 1:length(namelist)){keep[i] = grep(namelist[i],comborows)[1]}
combomat_blocked = combomat_subset[keep,]

#Updates to the host range gene heatmap:
#transpose (samples on rows instead of columns)
#reorder rows to match figure 6
combomat_reorder = combomat_blocked[,c('C1M','C2M','C3M','C1I','C2I','C3I',
										'S1M','S2M','S3M','S1I','S2I','S3I',
										'C1M_t20','C2M_t20','C3M_t20','C1I_t20','C2I_t20','C3I_t20',
										'S1M_t20','S2M_t20','S3M_t20','S1I_t20','S2I_t20','S3I_t20')]

combomat_reorder_subset = combomat_reorder[-which(apply(combomat_reorder,1,sum)==1),]
# update row labels to use corrected positions, then reorder them as done previously
newrows = character()
rowsplit = strsplit(rownames(combomat_reorder_subset), '_')
bobgene = grep('terminase', tempgtf_bob$product)[1]
bobterm = tempgtf_bob$start[bobgene]
cargene = grep('terminase', newgtf_car$product)[1]
carterm = newgtf_car$start2[cargene]
for(i in 1:length(rowsplit)){
	if(rowsplit[[i]][1] == 'Bob2'){
		if(as.numeric(rowsplit[[i]][2]) >= bobterm){
			newpos = as.numeric(rowsplit[[i]][2]) - bobterm + 1
		}
		if(as.numeric(rowsplit[[i]][2]) < bobterm){
			newpos = nbob - (bobterm - as.numeric(rowsplit[[i]][2]))
		}
		newrows[i] = paste(c('Bob',newpos,rowsplit[[i]][3],rowsplit[[i]][4]), collapse='_')
	}
	if(rowsplit[[i]][1] == 'Car'){
		newpos = ncar - as.numeric(rowsplit[[i]][2]) + 1
		if(newpos >= carterm){
			newpos = newpos - carterm + 1
		}
		if(newpos < carterm){
			newpos = ncar - (carterm - newpos)
		}
		newrows[i] = paste(c('Car',newpos,rowsplit[[i]][3],rowsplit[[i]][4]), collapse='_')
	}
}
rownames(combomat_reorder_subset) = newrows
combomat_bob = combomat_reorder_subset[grep('Bob',rownames(combomat_reorder_subset)),]
combomat_car = combomat_reorder_subset[grep('Car',rownames(combomat_reorder_subset)),]

#reorder rows of each to match nucleotide position
bobposmat = as.data.frame(cbind(rownames(combomat_bob),str_extract(rownames(combomat_bob),'[0-9]+') %>% as.numeric()))
carposmat = as.data.frame(cbind(rownames(combomat_car),str_extract(rownames(combomat_car),'[0-9]+') %>% as.numeric()))
colnames(bobposmat) = colnames(carposmat) = c('Variant','Position')
bobposmat$Position = as.numeric(bobposmat$Position)
carposmat$Position = as.numeric(carposmat$Position)
bobposmat = arrange(bobposmat, Position)
carposmat = arrange(carposmat, Position)
combomat_bob = combomat_bob[bobposmat$Variant,]
combomat_car = combomat_car[carposmat$Variant,]

#Update heatmaps with rearranged positions
bobgtf = tempgtf_bob
cargtf = newgtf_car

#next, match up the gtf to the positions in bobposmat and carposmat and get functions. Might need to rerun blast
bobposmat$product = bobposmat$category = rep(NA, nrow(bobposmat))
carposmat$product = carposmat$category = rep(NA, nrow(carposmat))
for(i in 1:nrow(bobposmat)){
	temp = which(bobgtf$newstart < bobposmat$Position[i] & bobgtf$newend > bobposmat$Position[i])
	if(length(temp) > 0){
		bobposmat$product[i] = bobgtf$product[temp]
		bobposmat$category[i] = bobgtf$category[temp]
	}
	if(length(temp) == 0){
		bobposmat$product[i] = 'intergenic'
		bobposmat$category[i] = 'intergenic'
	}
}
for(i in 1:nrow(carposmat)){
	temp = which(cargtf$newstart < carposmat$Position[i] & cargtf$newend > carposmat$Position[i])
	if(length(temp) > 0){
		carposmat$product[i] = cargtf$product[temp]
		carposmat$category[i] = cargtf$category[temp]
	}
	if(length(temp) == 0){
		carposmat$product[i] = 'intergenic'
		carposmat$category[i] = 'intergenic'
	}
}

#specify colors to define different functional categories within the column names. Match the colors from the genome maps and add a new color for "intergenic"

bobposmat$color = rep(NA, nrow(bobposmat))
bobposmat$color[bobposmat$category == 'assembly'] = darkpal[1]
bobposmat$color[bobposmat$category == 'capsid'] = darkpal[2]
bobposmat$color[bobposmat$category == 'hypothetical'] = 'black'
bobposmat$color[bobposmat$category == 'other'] = darkpal[3]
bobposmat$color[bobposmat$category == 'recombination'] = darkpal[4]
bobposmat$color[bobposmat$category == 'replication'] = darkpal[5]
bobposmat$color[bobposmat$category == 'tail'] = darkpal[6]
bobposmat$color[bobposmat$category == 'intergenic'] = 'dodgerblue3'

bobposmat$fontface = rep('bold',nrow(bobposmat))
bobposmat$fontface[bobposmat$color == 'black'] = 'plain'

#Prior to plotting, update labels to match manuscript nomenclature
colnames(combomat_bob) = str_replace(colnames(combomat_bob), 'C', 'D')
colnames(combomat_bob) = str_replace(colnames(combomat_bob), 'I', 'P')

ht1 = Heatmap(t(combomat_bob),
	row_order = colnames(combomat_bob),
	cluster_columns = FALSE,
	show_column_dend = FALSE,
	#top_annotation = ta,
	#right_annotation = ha3,
	rect_gp = gpar(col = "black", lwd = 0.5),
	cluster_rows = FALSE,
	col = c('white','red'),
	row_names_gp = gpar(fontsize=12, fontface='bold'),
	column_names_gp = gpar(fontsize=12, fontface=bobposmat$fontface, col = bobposmat$color),
	column_names_side = "bottom",
	width=nrow(combomat_bob)*unit(5,'mm'),
	height=ncol(combomat_bob)*unit(4,'mm'),
	show_heatmap_legend = FALSE,
	#name = "Success Rate"
	)

#Since cousins don't have Car, let's set those values to NA
combomat_car[,grep('^C',colnames(combomat_car))] = NA
carposmat$color = rep(NA, nrow(carposmat))
carposmat$color[carposmat$category == 'assembly'] = darkpal[1]
carposmat$color[carposmat$category == 'capsid'] = darkpal[2]
carposmat$color[carposmat$category == 'hypothetical'] = 'black'
carposmat$color[carposmat$category == 'other'] = darkpal[3]
carposmat$color[carposmat$category == 'recombination'] = darkpal[4]
carposmat$color[carposmat$category == 'replication'] = darkpal[5]
carposmat$color[carposmat$category == 'tail'] = darkpal[6]
carposmat$color[carposmat$category == 'intergenic'] = 'dodgerblue3'
carposmat$fontface = rep('bold',nrow(carposmat))
carposmat$fontface[carposmat$color == 'black'] = 'plain'

#Prior to plotting, update labels to match manuscript nomenclature
colnames(combomat_car) = str_replace(colnames(combomat_car), 'C', 'D')
colnames(combomat_car) = str_replace(colnames(combomat_car), 'I', 'P')
ht2 = Heatmap(t(combomat_car),
	row_order = colnames(combomat_car),
	cluster_columns = FALSE,
	show_column_dend = FALSE,
	#top_annotation = ta,
	#right_annotation = ha3,
	rect_gp = gpar(col = "black", lwd = 0.5),
	cluster_rows = FALSE,
	col = c('white','red'),
	row_names_gp = gpar(fontsize=12, fontface='bold'),
	column_names_gp = gpar(fontsize=12, fontface=carposmat$fontface, col = carposmat$color),
	column_names_side = "bottom",
	width=nrow(combomat_car)*unit(5,'mm'),
	height=ncol(combomat_car)*unit(4,'mm'),
	show_heatmap_legend = FALSE,
	#name = "Success Rate"
	)

#Prior to plotting, update labels to match manuscript nomenclature
colnames(hostdata) = str_replace(colnames(hostdata), 'C', 'D')
colnames(hostdata) = str_replace(colnames(hostdata), 'I', 'P')
hostdata_subset = hostdata[,colnames(combomat_bob)]
colscheme = c(
  "Resistance" = "#eff3ff",
  "Partial Suppression" = "#6baed6",
  "Full Suppression" = "#0d2c39"
)

hostcols = structure(c("#eff3ff","#6baed6","#0d2c39"), names=c("1","2","3"))


ht3 = Heatmap(t(hostdata_subset),
	row_order = colnames(hostdata_subset),
	cluster_columns = FALSE,
	show_column_dend = FALSE,
	#top_annotation = ta,
	#right_annotation = ha3,
	rect_gp = gpar(col = "black", lwd = 0.5),
	cluster_rows = FALSE,
	col = hostcols,
	row_names_gp = gpar(fontsize=12, fontface='bold'),
	column_names_gp = gpar(fontsize=12, fontface='bold'),
	column_names_side = "bottom",
	width=nrow(hostdata_subset)*unit(5,'mm'),
	height=ncol(hostdata_subset)*unit(4,'mm'),
	show_heatmap_legend = FALSE,
	#name = "Success Rate"
	)

figsize = calc_ht_size(ht3+ht1+ht2)
pdf('./hostrange_mutations_heatmap_FigureSF4.pdf',height=figsize[2],width=figsize[1])
print(ht3+ht1+ht2)
dev.off()

plotdata = t(hostdata[,c('CCS4','SDS1','Bob','Bill','Car','Ancestral.Cousins','Ancestral.Sisters')])

#Update labels
rownames(plotdata)[6:7] = c('Ancestral.Dissimilar', 'Ancestral.Similar')
ht4 = Heatmap(plotdata,
	row_order = rownames(plotdata),
	cluster_columns = FALSE,
	show_column_dend = FALSE,
	#top_annotation = ha11,
	#right_annotation = ha3,
	rect_gp = gpar(col = "black", lwd = 0.5),
	cluster_rows = FALSE,
	col = hostcols,
	row_names_gp = gpar(fontsize=12, fontface='bold'),
	column_names_gp = gpar(fontsize=12, fontface='bold'),
	column_names_side = "bottom",
	width=ncol(plotdata)*unit(5,'mm'),
	height=nrow(plotdata)*unit(4,'mm'),
	show_heatmap_legend = FALSE,
)

figsize=calc_ht_size(ht4)
pdf('./ancestral_hostmat_FigureSF4_inset.pdf',height=figsize[2],width=figsize[1])
print(ht4)
dev.off()

#Analysis for Supplemental Figure 6: correlation between isolate and population host ranges (usiing any suppression)
#Want to do a comparison of isolates relative to their endpoint population for host range

allhost = read.csv("./HostData/HostRangeSimplified1.csv")

#Update the columns to be easier to read. Only want to compare the endpoints
keep = grep('Evolved',colnames(allhost))
keep = c(keep, grep('T20', colnames(allhost)))
allhost_subset = allhost[,c(1,keep)]
rownames(allhost_subset) = allhost_subset$X
allhost_subset = allhost_subset[,-1]

colnames(allhost_subset) = str_extract(colnames(allhost_subset), '[SC][123].+')
drop = grep('R2', colnames(allhost_subset))
allhost_subset = allhost_subset[,-drop]

allhost_subset_part = allhost_subset
allhost_subset_part[allhost_subset == 1] = 0
allhost_subset_part[allhost_subset %in% c(2,3)] = 1
corvec_part = c(
	cor(allhost_subset_part$S1M, allhost_subset_part$S1M_T20),
	cor(allhost_subset_part$S2M, allhost_subset_part$S2M_T20),
	cor(allhost_subset_part$S3M, allhost_subset_part$S3M_T20),
	cor(allhost_subset_part$S1I, allhost_subset_part$S1I_T20),
	cor(allhost_subset_part$S2I, allhost_subset_part$S2I_T20),
	cor(allhost_subset_part$S3I, allhost_subset_part$S3I_T20),
	cor(allhost_subset_part$C1M, allhost_subset_part$C1M_T20),
	cor(allhost_subset_part$C2M, allhost_subset_part$C2M_T20),
	cor(allhost_subset_part$C3M, allhost_subset_part$C3M_T20),
	cor(allhost_subset_part$C1I, allhost_subset_part$C1I_T20),
	cor(allhost_subset_part$C2I, allhost_subset_part$C2I_T20),
	cor(allhost_subset_part$C3I, allhost_subset_part$C3I_T20)
)

samplelist = c('S1M','S2M','S3M','S1I','S2I','S3I','C1M','C2M','C3M','C1I','C2I','C3I')
names(corvec_part) = samplelist

cor.df = data.frame(ID = samplelist, cor_part = 0, cor_part = 0, p_part = 0, p_part = 0, cocktail = str_extract(samplelist, '[SC]'), hostmix = str_extract(samplelist, '[IM]'))

for(i in 1:length(samplelist)){
	allhost_subset_part = allhost_subset
	allhost_subset_part[allhost_subset == 3] = 1
	allhost_subset_part[allhost_subset != 3] = 0
	temptest = cor.test(allhost_subset_part[[samplelist[i]]], allhost_subset_part[[paste0(samplelist[i],'_T20')]])
	cor.df$cor_part[i] = temptest$estimate
	cor.df$p_part[i] = temptest$p.value
	
	allhost_subset_part = allhost_subset
	allhost_subset_part[allhost_subset %in% c(2,3)] = 1
	allhost_subset_part[allhost_subset == 1] = 0
	temptest = cor.test(allhost_subset_part[[samplelist[i]]], allhost_subset_part[[paste0(samplelist[i],'_T20')]])
	cor.df$cor_part[i] = temptest$estimate
	cor.df$p_part[i] = temptest$p.value
}

cor.df$ID = str_replace(cor.df$ID, 'C','D')
cor.df$cocktail = str_replace(cor.df$cocktail, 'C','Dissimilar')
cor.df$cocktail[cor.df$cocktail == 'S'] = 'Similar'
cor.df$hostmix = str_replace(cor.df$hostmix, 'I', 'Parallel')
cor.df$hostmix[cor.df$hostmix == 'M'] = 'Mixed'
corplot = ggplot(cor.df, aes(x=cocktail, y=cor_part, fill=hostmix))+geom_boxplot()+theme_bw()+
	ylab('Correlation')+
	xlab('Phage Relatedness')+
	labs(fill = 'Host Mixing')+
	ylim(c(min(cor.df[,c('cor_part','cor_part')]),max(cor.df[,c('cor_part','cor_part')])))+
	theme(text = element_text(face='bold', size=12))

png('../hostrange_correlation_iso_pop.png', height=400, width=500)
print(corplot)
dev.off()

pdf('../hostrange_correlation_iso_pop.pdf', height=4, width=5)
print(corplot)
dev.off()