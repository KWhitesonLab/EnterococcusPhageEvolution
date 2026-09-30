# Code to generate mutation plots (Figure 4, Supplemental Figure 2, Supplemental Figure 3)
##Note: the following requires steps to generate tempgtf_bob, newgtf_car, and tempgtf_ccs4, as described in Figure_2c.txt
# Note: this code includes steps that refer to naming conventions that were changed for the final manuscript. Previously, we referred to the "Dissimilar" phage cocktails as "cousins" and the "Similar" cocktails as "sisters" or "siblings". For the cousins, the sample IDs started with C. In the final figures, these are updated to begin with "D", and other labels are also corrected to agree with the manuscript naming scheme.

library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(data.table)
library(stringr)
library(RColorBrewer)

samplelist = c('C1I','C1M','C2I','C2M','C3I','C3M','S1I','S1M','S2I','S2M','S3I','S3M')
timelist = c('t1','t2','t3','t4','t5','t10','t15','t20')

setwd('/mnt/d/Research/WhitesonCollab/Data/RawReads/trimmed')
dirlist = dir()
vcflist_trimmed = list()

for(i in 1:length(samplelist)){for(j in 1:length(timelist)){
	if(length(grep("C",samplelist[i])) > 0){
		vcflist_trimmed[[paste0(samplelist[i],'_',timelist[j])]] =
		read.table(paste0('cousins_',timelist[j],'/',samplelist[i],'_on_combo_trimmed_output.g.vcf.gz'), skip=29, header=F, sep='\t')
	}
	if(length(grep("S",samplelist[i])) > 0){
		vcflist_trimmed[[paste0(samplelist[i],'_',timelist[j])]] =
		read.table(paste0('sisters_',timelist[j],'/',samplelist[i],'_on_combo_trimmed_output.g.vcf.gz'), skip=29, header=F, sep='\t')
	}
}}

for(i in 1:length(vcflist_trimmed)){
	vcflist_trimmed[[i]]$label = paste0(vcflist_trimmed[[i]][,1],vcflist_trimmed[[i]][,2],vcflist_trimmed[[i]][,3],vcflist_trimmed[[i]][,4],vcflist_trimmed[[i]][,5])
}

#Add t0 which has a different naming convention
t0cousins_vcf_trimmed = list()
setwd("cousins_t0")
t0files = dir(pattern='trimmed_output.g.vcf.gz$')
for(i in 1:length(t0files)){
	t0cousins_vcf_trimmed[[t0files[i]]] = read.table(t0files[i], skip=29, header=F, sep='\t')
}

cousins_background = t0cousins_vcf_trimmed[[1]][,1:5]
for(i in 2:length(t0cousins_vcf_trimmed)){
	cousins_background = rbind(cousins_background, t0cousins_vcf_trimmed[[i]][,1:5])
}
cousins_background = distinct(cousins_background)
cousins_background$label = paste0(cousins_background[,1],cousins_background[,2],cousins_background[,3],cousins_background[,4],cousins_background[,5])

t0sisters_vcf_trimmed = list()
setwd("../sisters_t0")
t0files = dir(pattern='trimmed_output.g.vcf.gz$')
for(i in 1:length(t0files)){
	t0sisters_vcf_trimmed[[t0files[i]]] = read.table(t0files[i], skip=29, header=F, sep='\t')
}

sisters_background = t0sisters_vcf_trimmed[[1]][,1:5]
for(i in 2:length(t0sisters_vcf_trimmed)){
	sisters_background = rbind(sisters_background, t0sisters_vcf_trimmed[[i]][,1:5])
}
sisters_background = distinct(sisters_background)
sisters_background$label = paste0(sisters_background[,1],sisters_background[,2],sisters_background[,3],sisters_background[,4],sisters_background[,5])

t0muts_trimmed = list()
for(i in 1:length(t0cousins_vcf_trimmed)){
	t0muts_trimmed[[names(t0cousins_vcf_trimmed)[i]]] = paste(t0cousins_vcf_trimmed[[i]]$V1,t0cousins_vcf_trimmed[[i]]$V2,t0cousins_vcf_trimmed[[i]]$V4,t0cousins_vcf_trimmed[[i]]$V5, sep='_')
}
for(i in 1:length(t0sisters_vcf_trimmed)){
	t0muts_trimmed[[names(t0sisters_vcf_trimmed)[i]]] = paste(t0sisters_vcf_trimmed[[i]]$V1,t0sisters_vcf_trimmed[[i]]$V2,t0sisters_vcf_trimmed[[i]]$V4,t0sisters_vcf_trimmed[[i]]$V5, sep='_')
}

t0muts_trimmed = unique(unlist(t0muts_trimmed))

vcftrim = list()
for(i in 1:length(vcflist_trimmed)){
	temp=str_extract(vcflist_trimmed[[i]]$V8,'DP=[0-9]+\\;') %>% str_replace_all('DP=','') %>% str_replace_all('\\;','')
	temp=as.numeric(temp)
	vcftrim[[i]] = vcflist_trimmed[[i]]
	drop = which(temp<5)
	if(length(drop)>0){
		vcftrim[[i]] = vcftrim[[i]][-drop,]
	}
}
names(vcftrim) = names(vcflist_trimmed)

setwd("../../../variants")
nocar = c('S2M','C1I','C1M','C2I','C2M','C3I','C3M')
hasboth = c('S1I','S1M','S2I','S3I','S3M')
isomuts = list()
for(i in 1:length(samplelist)){
	if(samplelist[i]%in%hasboth){
		isomuts[[paste0(samplelist[i],'_bob')]] = read.table(paste0(samplelist[i],'_bob.vcf.gz'), skip=29, header=F, sep='\t')
		isomuts[[paste0(samplelist[i],'_bob')]] = isomuts[[paste0(samplelist[i],'_bob')]][grep(';AF=1',isomuts[[paste0(samplelist[i],'_bob')]]$V8),]
		isomuts[[paste0(samplelist[i],'_car')]] = read.table(paste0(samplelist[i],'_car.vcf.gz'), skip=29, header=F, sep='\t')
		isomuts[[paste0(samplelist[i],'_car')]] = isomuts[[paste0(samplelist[i],'_car')]][grep(';AF=1',isomuts[[paste0(samplelist[i],'_car')]]$V8),]
	}
	if(samplelist[i]%in%nocar){
		isomuts[[samplelist[i]]] = read.table(paste0(samplelist[i],'_output.g.vcf.gz'), skip=29, header=F, sep='\t')
		isomuts[[samplelist[i]]] = subset(isomuts[[samplelist[i]]], V1!='Bob1')
		isomuts[[samplelist[i]]] = isomuts[[samplelist[i]]][grep(';AF=1',isomuts[[samplelist[i]]]$V8),]
	}
}

bobmutpos = carmutpos = ccs4mutpos = list()

for(i in 1:length(vcflist_trimmed)){
	tempbob = grep('Bob', vcflist_trimmed[[i]]$V1)
	if(length(tempbob)>0){
		bobmutpos[[names(vcflist_trimmed)[i]]] = vcflist_trimmed[[i]]$V2[tempbob]
	}
	tempcar = grep('Car', vcflist_trimmed[[i]]$V1)
	if(length(tempcar)>0){
		carmutpos[[names(vcflist_trimmed)[i]]] = vcflist_trimmed[[i]]$V2[tempcar]
	}
}

#repeat for isolates
isobobmutpos = isocarmutpos = isoccs4mutpos = list()

for(i in 1:length(isomuts)){
	tempbob = grep('Bob', isomuts[[i]]$V1)
	if(length(tempbob)>0){
		isobobmutpos[[names(isomuts)[i]]] = isomuts[[i]]$V2[tempbob]
	}
	tempcar = grep('Car', isomuts[[i]]$V1)
	if(length(tempcar)>0){
		isocarmutpos[[names(isomuts)[i]]] = isomuts[[i]]$V2[tempcar]
	}
}

#clean up any background
t0bob = t0muts_trimmed[grep('Bob', t0muts_trimmed)]
t0bob = str_extract_all(t0bob, '[0-9][0-9]+')
t0bob = as.numeric(unlist(t0bob))
for(i in 1:length(bobmutpos)){
	drops = which(bobmutpos[[i]]%in%t0bob)
	if(length(drops)>0){
		bobmutpos[[i]] = bobmutpos[[i]][-drops]
	}
}

t0car = t0muts_trimmed[grep('Car', t0muts_trimmed)]
t0car = as.numeric(unlist(str_extract_all(t0car, '[0-9][0-9]+')))
for(i in 1:length(carmutpos)){
	drops = which(carmutpos[[i]]%in%t0car)
	if(length(drops)>0){
		carmutpos[[i]] = carmutpos[[i]][-drops]
	}
}

for(i in 1:length(isobobmutpos)){
	drops = which(isobobmutpos[[i]]%in%t0bob)
	if(length(drops)>0){
		isobobmutpos[[i]] = isobobmutpos[[i]][-drops]
	}
}

for(i in 1:length(isocarmutpos)){
	drops = which(isocarmutpos[[i]]%in%t0car)
	if(length(drops)>0){
		isocarmutpos[[i]] = isocarmutpos[[i]][-drops]
	}
}

#Just need a two column matrix listing all of the results
bobmutdf = c()
for(i in 1:length(bobmutpos)){
	nmut = length(bobmutpos[[i]])
	tempname = rep(names(bobmutpos)[i], nmut)
	tempmat = cbind(rep(tempname, nmut), bobmutpos[[i]])
	bobmutdf = rbind(bobmutdf, tempmat)
}
bobmutdf = as.data.frame(bobmutdf)
colnames(bobmutdf) = c('Sample','Position')
bobmutdf$Position = as.numeric(bobmutdf$Position)
bobmutdf$Sample = factor(bobmutdf$Sample, levels=names(bobmutpos))
bobmutdf$Time = str_extract(as.character(bobmutdf$Sample), 't[0-9]+')
bobmutdf$Time = as.numeric(str_replace(bobmutdf$Time,'t',''))
bobmutdf$Replicate = str_extract(as.character(bobmutdf$Sample), '[CS][123][IM]')
bobmutdf$Treatment = bobmutdf$Replicate
bobmutdf$Treatment[grep('C', bobmutdf$Treatment)] = 'Cousins'
bobmutdf$Treatment[grep('S', bobmutdf$Treatment)] = 'Sisters'

bobref = readLines("/mnt/d/Research/WhitesonCollab/Data/GoogleDrive/Bob/2023-08-11_Bob_flyeAssembly_selectedContig.fasta")
bobref=c(">Bob",paste(bobref[2:length(bobref)], collapse=""))

carref = readLines("/mnt/d/Research/WhitesonCollab/Data/GoogleDrive/Car/2023-08-11_Car_flyeAssembly_selectedContig.fasta")
carref = c(">Car", paste(carref[2:length(carref)], collapse=""))

ccs4ref = readLines("../GoogleDrive/CCS4/2023-08-25_CCS4_flyeAssembly_selectedContig.fasta")
ccs4ref = c(">CCS4", paste(ccs4ref[2:length(ccs4ref)], collapse=''))

nbob = nchar(bobref[2])
ncar = nchar(carref[2])
nccs4 = nchar(ccs4ref[2])

carmutdf = c()
for(i in 1:length(carmutpos)){
	nmut = length(carmutpos[[i]])
	tempname = rep(names(carmutpos)[i], nmut)
	tempmat = cbind(rep(tempname, nmut), carmutpos[[i]])
	carmutdf = rbind(carmutdf, tempmat)
}
carmutdf = as.data.frame(carmutdf)
colnames(carmutdf) = c('Sample','Position')
carmutdf$Position = as.numeric(carmutdf$Position)
carmutdf$Sample = factor(carmutdf$Sample, levels=names(carmutpos))
carmutdf$Time = str_extract(as.character(carmutdf$Sample), 't[0-9]+')
carmutdf$Time = as.numeric(str_replace(carmutdf$Time,'t',''))
carmutdf$Replicate = str_extract(as.character(carmutdf$Sample), '[CS][123][IM]')
carmutdf$Treatment = carmutdf$Replicate
carmutdf$Treatment[grep('C', carmutdf$Treatment)] = 'Cousins'
carmutdf$Treatment[grep('S', carmutdf$Treatment)] = 'Sisters'

#isolates
isobobmutdf = c()
for(i in 1:length(isobobmutpos)){
	nmut = length(isobobmutpos[[i]])
	tempname = rep(names(isobobmutpos)[i], nmut)
	tempmat = cbind(rep(tempname, nmut), isobobmutpos[[i]])
	isobobmutdf = rbind(isobobmutdf, tempmat)
}
isobobmutdf = as.data.frame(isobobmutdf)
colnames(isobobmutdf) = c('Sample','Position')
isobobmutdf$Position = as.numeric(isobobmutdf$Position)
isobobmutdf$Sample = factor(isobobmutdf$Sample, levels=names(isobobmutpos))
isobobmutdf$Time = "t20 Isolate"
isobobmutdf$Replicate = str_extract(as.character(isobobmutdf$Sample), '[CS][123][IM]')
isobobmutdf$Treatment = isobobmutdf$Replicate
isobobmutdf$Treatment[grep('C', isobobmutdf$Treatment)] = 'Cousins'
isobobmutdf$Treatment[grep('S', isobobmutdf$Treatment)] = 'Sisters'
#combine with bobmutdf
allbobmut = rbind(bobmutdf, isobobmutdf)
allbobmut$Time = factor(allbobmut$Time, levels=c(1,2,3,4,5,10,15,20,'t20 Isolate'))

#Rebuild, shifting positions to match the terminase-shifted version of the genome map
termgene = grep('terminase', tempgtf_bob$product)[1]
termstart = tempgtf_bob$start[termgene]
allbobmut$NewPosition = 0
allbobmut$NewPosition[allbobmut$Position >= termstart] = allbobmut$Position[allbobmut$Position >= termstart] - termstart + 1
allbobmut$NewPosition[allbobmut$Position < termstart] = nbob - (termstart - allbobmut$Position[allbobmut$Position < termstart])
#Replot with NewPosition for x-axis
samplelist = c('C1M','C2M','C3M','C1I','C2I','C3I','S1M','S2M','S3M','S1I','S2I','S3I')

#Update labels for manuscript. C -> D. I -> P
allbobmut$Replicate = as.character(allbobmut$Replicate)
allbobmut$Replicate = str_replace(allbobmut$Replicate, 'C', 'D')
allbobmut$Replicate = str_replace(allbobmut$Replicate, 'I', 'P')
allbobmut$Replicate = factor(allbobmut$Replicate, levels = c('D1M','D2M','D3M','D1P','D2P','D3P','S1M','S2M','S3M','S1P','S2P','S3P'))

tempgtf_bob$molecule = 'Bob'
genome_bob = ggplot(tempgtf_bob, aes(xmin=newstart, xmax=newend, y=molecule, forward=orientation, fill = category,size=sizeclass))+
	geom_gene_arrow(arrow_body_height = grid::unit(2, 'mm'), arrowhead_width = grid::unit(2,'mm'))+
	scale_fill_manual(values=c(darkpal[1:2],'grey',darkpal[3:6]))+
	scale_size_manual(values=c(0.4,0.1))+
	guides(size = 'none')+
	ylab('')+
	theme_genes()+
	theme(legend.position = 'none',
	text=element_text(size=12),
	axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank())

pdf('../bob_genome_map.pdf', height = 1, width=7.5)
print(genome_bob)
dev.off()

pdf('../bob_genome_map_narrow.pdf', height = 1, width=5)
print(genome_bob)
dev.off()

pdf('../bob_genome_map_narrower.pdf', height = 1, width=3.75)
print(genome_bob)
dev.off()
bobmut_sim = ggplot(subset(allbobmut, Treatment=='Sisters'), aes(x=NewPosition, y=Time))+facet_grid(Replicate~.)+
	geom_point(size=0.1)+
	xlim(c(0,nbob))+
	xlab('Position')+
	ylab('Experimental Cycle')+
	theme_bw()+
	theme(text=element_text(face='bold', size=12),
	strip.text.y = element_text(angle = 0, size = 15),
	plot.margin = unit(c(1, 1, 2, 1), "cm"))

pdf('../bobmutplot_wIsolates_similar_Figure4a.pdf', height=8, width=7.5)
print(bobmut_sim)
dev.off()

pdf('../bobmutplot_wIsolates_similar_Figure4a_small.pdf', height=5.5, width=3.75)
print(bobmut_sim+
	theme(text=element_text(face='bold', size=8),
	panel.spacing = unit(2, 'pt'),
	strip.text.y = element_text(angle = 0, size = 4)))
dev.off()

bobmut_dis = ggplot(subset(allbobmut, Treatment=='Cousins'), aes(x=NewPosition, y=Time))+facet_grid(Replicate~.)+
	geom_point(size=0.1)+
	xlim(c(0,nbob))+
	xlab('Position')+
	ylab('Experimental Cycle')+
	theme_bw()+
	theme(text=element_text(face='bold', size=12),
	strip.text.y = element_text(angle = 0, size=15),
	plot.margin = unit(c(1, 1, 2, 1), "cm"))

pdf('../bobmutplot_wIsolates_dissimilar_Figure4b.pdf', height=8, width=5)
print(bobmut_dis)
dev.off()

pdf('../bobmutplot_wIsolates_dissimilar_Figure4b_small2.pdf', height=5.5, width=3.75)
print(bobmut_dis+
	theme(text=element_text(face='bold', size=8),
	panel.spacing = unit(2, 'pt'),
	strip.text.y = element_text(angle = 0, size = 4)))
dev.off()

#car
isocarmutdf = c()
for(i in 1:length(isocarmutpos)){
	nmut = length(isocarmutpos[[i]])
	tempname = rep(names(isocarmutpos)[i], nmut)
	tempmat = cbind(rep(tempname, nmut), isocarmutpos[[i]])
	isocarmutdf = rbind(isocarmutdf, tempmat)
}
isocarmutdf = as.data.frame(isocarmutdf)
colnames(isocarmutdf) = c('Sample','Position')
isocarmutdf$Position = as.numeric(isocarmutdf$Position)
isocarmutdf$Sample = factor(isocarmutdf$Sample, levels=names(isocarmutpos))
isocarmutdf$Time = "t20 Isolate"
isocarmutdf$Replicate = str_extract(as.character(isocarmutdf$Sample), '[CS][123][IM]')
isocarmutdf$Treatment = isocarmutdf$Replicate
isocarmutdf$Treatment[grep('C', isocarmutdf$Treatment)] = 'Cousins'
isocarmutdf$Treatment[grep('S', isocarmutdf$Treatment)] = 'Sisters'

#combine with carmutdf
allcarmut = rbind(carmutdf, isocarmutdf)
allcarmut$Time = factor(allcarmut$Time, levels=c(1,2,3,4,5,10,15,20,'t20 Isolate'))
allcarmut$Replicate = factor(allcarmut$Replicate, levels=samplelist[7:12])
allcarmut$NewPosition = ncar - allcarmut$Position + 1 #This should reverse the values to match the first step in newgtf_car
termgene = grep('terminase', newgtf_car$product)[1]
termstart = newgtf_car$start2[termgene]  #use start2 to match the reversed start values
allcarmut$NewPosition2 = 0
allcarmut$NewPosition2[allcarmut$NewPosition >= termstart] = allcarmut$NewPosition[allcarmut$NewPosition >= termstart] - termstart + 1
allcarmut$NewPosition2[allcarmut$NewPosition < termstart] = ncar - (termstart - allcarmut$NewPosition[allcarmut$NewPosition < termstart])
allcarmut$Replicate = as.character(allcarmut$Replicate)
allcarmut$Replicate = str_replace(allcarmut$Replicate, 'I', 'P')
allcarmut$Replicate = factor(allcarmut$Replicate, levels = c('S1M','S2M','S3M','S1P','S2P','S3P'))

newgtf_car$molecule = 'Car'
genome_car = ggplot(newgtf_car, aes(xmin=newstart, xmax=newend, y=molecule, forward=orientation, fill = category,size=sizeclass))+
	geom_gene_arrow(arrow_body_height = grid::unit(2, 'mm'), arrowhead_width = grid::unit(2,'mm'))+
	scale_fill_manual(values=c(darkpal[1:2],'grey',darkpal[3:6]))+
	scale_size_manual(values=c(0.4,0.1))+
	guides(size = 'none')+
	ylab('')+
	theme_genes()+
	theme(legend.position = 'none',
	text=element_text(size=12),
	axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank())

pdf('../car_genome_map.pdf', height = 1, width = 7.5)
print(genome_car)
dev.off()

pdf('../car_genome_map_narrow.pdf', height = 1, width = 3.75)
print(genome_car)
dev.off()
carmutplot = ggplot(allcarmut, aes(x=NewPosition2, y=Time))+facet_grid(Replicate~.)+
	geom_point(size=0.1)+
	xlim(c(0,nbob))+
	xlab('Position')+
	ylab('Experimental Cycle')+theme_bw()+
	theme(text=element_text(face='bold'),
	strip.text.y = element_text(angle = 0),
	plot.margin = unit(c(1, 1, 2, 1), "cm"))
	
pdf('../carmutplot_wIsolates_similar_FigureSF2a.pdf', height=8, width=7.5)
print(carmutplot)
dev.off()

pdf('../carmutplot_wIsolates_similar_FigureSF2a_small.pdf', height=5.5, width=3.75)
print(carmutplot+
	theme(text=element_text(face='bold', size=8),
	panel.spacing = unit(2, 'pt'),
	strip.text.y = element_text(angle = 0, size = 4)))
dev.off()

test = (genome_car + theme(plot.margin = margin(0,1,-5,1)))/carmutplot
#CCS4 for Cousins
ccs4mutpos = list()  #none detected in the isolates?

for(i in 1:length(vcflist_trimmed)){
	tempccs4= grep('CCS4', vcflist_trimmed[[i]]$V1)
	if(length(tempccs4)>0){
		ccs4mutpos[[names(vcflist_trimmed)[i]]] = vcflist_trimmed[[i]]$V2[tempccs4]
	}
}

#clean up any background
t0ccs4 = t0muts_trimmed[grep('CCS4', t0muts_trimmed)]
t0ccs4 = str_extract_all(t0ccs4, '[0-9][0-9]+')
t0ccs4 = as.numeric(unlist(t0ccs4))
for(i in 1:length(ccs4mutpos)){
	drops = which(ccs4mutpos[[i]]%in%t0ccs4)
	if(length(drops)>0){
		ccs4mutpos[[i]] = ccs4mutpos[[i]][-drops]
	}
}

ccs4mutdf = c()
for(i in 1:length(ccs4mutpos)){
	nmut = length(ccs4mutpos[[i]])
	tempname = rep(names(ccs4mutpos)[i], nmut)
	tempmat = cbind(rep(tempname, nmut), ccs4mutpos[[i]])
	ccs4mutdf = rbind(ccs4mutdf, tempmat)
}
ccs4mutdf = as.data.frame(ccs4mutdf)
colnames(ccs4mutdf) = c('Sample','Position')
ccs4mutdf$Position = as.numeric(ccs4mutdf$Position)
ccs4mutdf$Sample = factor(ccs4mutdf$Sample, levels=names(ccs4mutpos))
ccs4mutdf$Time = str_extract(as.character(ccs4mutdf$Sample), 't[0-9]+')
ccs4mutdf$Time = as.numeric(str_replace(ccs4mutdf$Time,'t',''))
ccs4mutdf$Replicate = str_extract(as.character(ccs4mutdf$Sample), '[CS][123][IM]')
ccs4mutdf$Treatment = ccs4mutdf$Replicate

#Shift values as in the genome map
termgene = grep('terminase', tempgtf_ccs4$product)[1]
termstart = tempgtf_ccs4 $start[termgene]
ccs4mutdf$NewPosition = 0
ccs4mutdf$NewPosition[ccs4mutdf$Position >= termstart] = ccs4mutdf$Position[ccs4mutdf$Position >= termstart] - termstart + 1
ccs4mutdf$NewPosition[ccs4mutdf$Position < termstart] = nccs4 - (termstart - ccs4mutdf$Position[ccs4mutdf$Position < termstart])
#add in empty rows for the isolate and times 1, 2, 4
ccs4mutdf$Time = factor(ccs4mutdf$Time, levels = c(1, 2, 3, 4, 5, 10, 15, 20, 't20 Isolate'))
ccs4mutdf$Replicate = factor(ccs4mutdf$Replicate, levels=samplelist[1:6])

#As above, update the labels for the manuscript
#Update labels for manuscript
ccs4mutdf$Replicate = as.character(ccs4mutdf$Replicate)
ccs4mutdf$Replicate = str_replace(ccs4mutdf$Replicate, 'C', 'D')
ccs4mutdf$Replicate = str_replace(ccs4mutdf$Replicate, 'I', 'P')
ccs4mutdf$Replicate = factor(ccs4mutdf$Replicate, levels = c('D1M','D2M','D3M','D1P','D2P','D3P'))

tempgtf_ccs4$molecule = 'CCS4'
genome_ccs4 = ggplot(tempgtf_ccs4, aes(xmin=newstart, xmax=newend, y=molecule, forward=orientation, fill = category,size=sizeclass))+
	geom_gene_arrow(arrow_body_height = grid::unit(2, 'mm'), arrowhead_width = grid::unit(2,'mm'))+
	scale_fill_manual(values=c(darkpal[1:2],'grey',darkpal[3:6]))+
	scale_size_manual(values=c(0.4,0.1))+
	guides(size = 'none')+
	ylab('')+
	theme_genes()+
	theme(legend.position = 'none',
	text=element_text(size=12),
	axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank())

pdf('../ccs4_genome_map.pdf', height=1, width=7.5)
print(genome_ccs4)
dev.off()

pdf('../ccs4_genome_map_narrow.pdf', height=1, width=5)
print(genome_ccs4)
dev.off()

pdf('../ccs4_genome_map_narrower_SF6.pdf', height=1, width=3.75)
print(genome_ccs4)
dev.off()

ccs4mutplot = ggplot(ccs4mutdf, aes(x=NewPosition, y=Time))+facet_grid(Replicate~., drop=FALSE)+
	geom_point(size=0.1)+
	xlim(c(0,nccs4))+
	xlab('Position')+
	ylab('Experimental Cycle')+theme_bw()+
	scale_y_discrete(drop=FALSE)+
	theme(text=element_text(face='bold'),
	strip.text.y = element_text(angle = 0),
	plot.margin = unit(c(1, 1, 2, 1), "cm"))
	
pdf('../ccs4mutplot_wIsolates_dissimilar_FigureSF4b.pdf', height=8, width=5)
print(ccs4mutplot)
dev.off()

pdf('../ccs4mutplot_wIsolates_dissimilar_FigureSF4b_small.pdf', height=5.5, width=3.75)
print(ccs4mutplot+
	theme(text=element_text(face='bold', size=8),
	panel.spacing = unit(2, 'pt'),
	strip.text.y = element_text(angle = 0, size = 4)))
dev.off()


#Add mutation summary boxplots
mut_bob = data.table(subset(distinct(allbobmut), Time==20))[,.N, by=.(Sample,Replicate)]
mut_car = data.table(subset(distinct(allcarmut), Time==20))[,.N, by=.(Sample,Replicate)]
mut_ccs4 = data.table(subset(distinct(ccs4mutdf), Time==20))[,.N, by=.(Sample,Replicate)]
temp = rbind(mut_bob, mut_car)
temp = rbind(temp, mut_ccs4)
temp$Treatment = str_extract(temp$Replicate, '[SD]')
temp$Treatment = ifelse(temp$Treatment=='D','Dissimilar','Similar')
temp$HostCondition = as.character(temp$HostCondition)
temp$HostCondition = str_extract(temp$Replicate, '[PM]')
temp$HostCondition = ifelse(temp$HostCondition == 'P','Parallel','Mixed')

## Add individual points ##
pdf('../mutation_summary_boxplot_wpoints_FigureSF3.pdf', height=4, width=6)
ggplot(temp, aes(x=Treatment, y=N, fill=HostCondition))+geom_boxplot() + 
geom_jitter(position=position_jitterdodge(jitter.width=0.2))+
ylab('# Mutations')+
theme_bw()+
theme(text=element_text(face='bold', size=15))
dev.off()