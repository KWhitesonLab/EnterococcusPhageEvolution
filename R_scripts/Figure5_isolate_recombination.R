#Code to generate Figure 5: Recombination between Bob and Car in endpoint isolates.
#Note: some steps are redundant with the beginning of the code to generate Figure 1c.

#Load libraries
library(gggenes)
library(stringr)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(patchwork)

#gbk2gtf parses a genbank file to gtf format.
gbk2gtf <- function(gbkfile, name=NULL){
  gbk = readLines(gbkfile)
  if(is.null(name)){
    name = str_replace(gbk[grep('DEFINITION', gbk)], 'DEFINITION  ', '')
  }
  cdspos = grep("CDS ", gbk)
  cds = gbk[cdspos]
  cds = str_replace_all(cds," ","")
  cds = str_replace_all(cds,"CDS","")
  comps = grep("complement", cds)
  cds = str_replace_all(cds,"complement\\(","")
  cds = str_replace_all(cds,"\\)","")
  cdssplit = strsplit(cds,"\\.\\.")
  prodpos = grep("product=",gbk)
  products = gbk[prodpos]
  products = str_extract(products,'\".+') %>% str_replace_all('\\\"','')
  genelist = paste0('gp_',c(1:length(cds)))
  starts = stops = numeric()
  for(i in 1:length(cdssplit)){
    starts[i] = as.numeric(cdssplit[[i]][1])
    stops[i] = as.numeric(cdssplit[[i]][2])
  }

  #gtf columns: start
  gtf = data.frame(molecule = name, gene=genelist, start = starts,
                   end = stops, strand = 'forward', orientation = 0, product = NA)
  for(i in 1:length(starts)){
      if(i %in% comps){
        gtf$strand[i] = 'reverse'
        gtf$orientation[i] = 1
      }
  }
  #Add products
  for(i in 1:nrow(gtf)){
	gtf$product[i] = products[which(prodpos > cdspos[i])[1]]
  }

  return(gtf)
}

darkpal = brewer.pal(8, 'Dark2')

#add maps for isolate genomes so that we can overlay these on recombination plots
#do not worry about rearranging genomes in any way for these, just do the quick annotations
#First, write out the gtfs with products for annotating manually
setwd("/mnt/d/Research/WhitesonCollab/Data/GoogleDrive")
samplelist = c('S1M','S1I','S2M','S2I','S3M','S3I')  #just need sisters for now since no recombination in the cousins
for(i in 1:length(samplelist)){
	tempgb = dir(samplelist[i], pattern='.gb')
	temp = gbk2gtf(paste0(samplelist[i],'/',tempgb))
	write.table(temp, paste0(samplelist[i],'/gtf_wproducts.txt'), row.names=F, quote=F, sep='\t')
}

#Update recombination figures with terminase reordering and inverted gene orientation where needed
genome_map <- function(gtffile, invert=FALSE, molecule){
	darkpal = brewer.pal(8, 'Dark2')
	tempgtf = read.csv(gtffile, header=T, sep='\t')
	tempgtf$length = tempgtf$end - tempgtf$start
	tempgtf$sizeclass = factor(ifelse(tempgtf$length > 500, 'large', 'small'))
	tempgtf$molecule = molecule
	tempgtf$molecule = str_replace(tempgtf$molecule, "I", "P") #correction for manuscript
	if(invert == FALSE){
		termgene = grep('terminase',tempgtf$product)[1]
		tempgtf$newstart = tempgtf$newend = 0
		tempgtf$newstart[termgene:nrow(tempgtf)] = tempgtf$start[termgene:nrow(tempgtf)] - tempgtf$start[termgene] + 1
		endshift = tempgtf$newstart[nrow(tempgtf)] + tempgtf$length[nrow(tempgtf)]
		tempgtf$newstart[1:(termgene - 1)] = tempgtf$start[1:(termgene - 1)] + endshift
		tempgtf$newend = tempgtf$newstart + tempgtf$length
	}
	if(invert == TRUE){
		tempgtf = tempgtf[nrow(tempgtf):1,]
		termgene = grep('terminase',tempgtf$product)[1]
		tempgtf$start2 = tempgtf$end2 = 0
		tempgtf$start2[1] = 1
		tempgtf$end2[1] = tempgtf$start2[1]+tempgtf$length[1]
		for(i in 2:nrow(tempgtf)){
			tempgtf$start2[i] = tempgtf$end2[i-1] + (tempgtf$start[i-1] - tempgtf$end[i])
			tempgtf$end2[i] = tempgtf$start2[i] + tempgtf$length[i]
		}
		tempgtf$newstart = tempgtf$newend = 0
		tempgtf$newstart[termgene:nrow(tempgtf)] = tempgtf$start2[termgene:nrow(tempgtf)] - tempgtf$start2[termgene] + 1
		endshift = tempgtf$newstart[nrow(tempgtf)] + tempgtf$length[nrow(tempgtf)]
		tempgtf$newstart[1:(termgene - 1)] = tempgtf$start2[1:(termgene - 1)] + endshift
		tempgtf$newend = tempgtf$newstart + tempgtf$length
		tempgtf$orientation = 1 - tempgtf$orientation
	}
	tempplot = ggplot(tempgtf, aes(xmin=newstart, xmax=newend, y=molecule, forward=orientation, fill = category,size=sizeclass))+
	geom_gene_arrow(arrow_body_height = grid::unit(2, 'mm'), arrowhead_width = grid::unit(2,'mm'))+
	scale_fill_manual(values=c(darkpal[1:2],'grey',darkpal[3:6]))+
	scale_size_manual(values=c(0.4,0.1))+
	guides(size = 'none')+
	theme_genes()+
	theme(legend.position = 'none')

	results = list()
	results$tempgtf = tempgtf
	results$plot = tempplot
	return(results)
}

setwd("/mnt/d/Research/WhitesonCollab/Data/alignments")

#Need to update recomplot to reorder as above

recomplot<-function(sampleID, reorder=TRUE){
	#easiest to build the new gtf first with code above
	tempgtf = genome_map(paste0('/mnt/d/Research/WhitesonCollab/Data/GoogleDrive/',sampleID,'/gtf_wproducts.txt'), molecule=sampleID, invert = reorder)$tempgtf
	bobfiles = dir(pattern = 'bob_regions.sort.depth')
	bobfile = bobfiles[grep(sampleID, bobfiles)]
	carfiles = dir(pattern = 'car_regions.sort.depth')
	carfile = carfiles[grep(sampleID, carfiles)]

	bobdepth = read.table(bobfile)
	cardepth = read.table(carfile)
	bobdepth$V1 = 'Bob'
	cardepth$V1 = 'Car'
	depth = rbind(bobdepth, cardepth)
	colnames(depth) = c('Reference', 'Position', 'Depth')
	depth$NewPosition = 0
	maxpos = max(depth$Position)
	termgene = grep('terminase', tempgtf$product)[1]

	if(reorder == FALSE){
		termpos = tempgtf$start[termgene]
	}
	if(reorder == TRUE){
		depth$Position = maxpos - depth$Position + 1
		termpos = tempgtf$start2[termgene]
	}
	depth$NewPosition[depth$Position >= termpos] = depth$Position[depth$Position >= termpos] - termpos + 1
	depth$NewPosition[depth$Position < termpos] = maxpos + depth$Position[depth$Position < termpos] - termpos + 1

	tempplot = ggplot(depth, aes(x=NewPosition, y=Depth, color=Reference))+
		geom_line(size=0.25)+
		xlab('Position')+
		scale_color_manual(values = c('dodgerblue','red'))+
		theme_bw()+
		theme(text = element_text(face='bold', size=15))

	results = list()
	results$depth = depth
	results$plot = tempplot
	return(results)
}


#Separate files for each panel in Figure 5, since each is paired with its own reference genome map.
#Note: file names reflect original nomenclature that referred to the "Parallel" condition (P) in the manuscript as "Individual" (I).

genome_S1P = genome_map('../GoogleDrive/S1I/gtf_wproducts.txt', invert=TRUE, molecule='S1I')$plot+
	ylab('')+
	theme(text=element_text(size=12),
	axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank())
rplot_S1P = recomplot('S1I', reorder=TRUE)$plot+theme(legend.position='none')+xlab('')

pdf('../S1P_recombination_plot_Figure5.pdf', height=2.5, width=7)
genome_S1P/rplot_S1P
dev.off()

genome_S1M = genome_map('../GoogleDrive/S1M/gtf_wproducts.txt', invert=TRUE, molecule='S1M')$plot+
	ylab('')+
	theme(text=element_text(size=12),
	axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank())
rplot_S1M = recomplot('S1M', reorder=TRUE)$plot+theme(legend.position='none')+xlab('')

pdf('../S1M_recombination_plot_Figure5.pdf', height=2.5, width=7)
genome_S1M/rplot_S1M
dev.off()

genome_S3M = genome_map('../GoogleDrive/S3M/gtf_wproducts.txt', invert=TRUE, molecule='S3M')$plot+
	ylab('')+
	theme(text=element_text(size=12),
	axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank())
rplot_S3M = recomplot('S3M', reorder=TRUE)$plot+theme(legend.position='none')+xlab('')
pdf('../S3M_recombination_plot_Figure5.pdf', height=2.5, width=7)
genome_S3M/rplot_S3M
dev.off()

genome_S2P = genome_map('../GoogleDrive/S2I/gtf_wproducts.txt', invert=FALSE, molecule='S2I')$plot+
	ylab('')+
	theme(text=element_text(size=12),
	axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank())
rplot_S2P = recomplot('S2I', reorder=FALSE)$plot+theme(legend.position='none')+xlab('')

pdf('../S2P_recombination_plot_Figure5.pdf', height=2.5, width=7)
genome_S2P/rplot_S2P
dev.off()

genome_S3P = genome_map('../GoogleDrive/S3I/gtf_wproducts.txt', invert=TRUE, molecule='S3I')$plot+
	ylab('')+
	theme(text=element_text(size=12),
	axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank())
rplot_S3P = recomplot('S3I', reorder=TRUE)$plot+theme(legend.position='none')+xlab('')

pdf('../S3P_recombination_plot_Figure5.pdf', height=2.5, width=7)
genome_S3P/rplot_S3P
dev.off()


#Also want to do S2M but it has no Car

genome_S2M = genome_map('../GoogleDrive/S2M/gtf_wproducts.txt', invert=FALSE, molecule='S2M')$plot+
	ylab('')+
	theme(text=element_text(size=12),
	axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank())

sampleID='S2M'
tempgtf = genome_map(paste0('/mnt/d/Research/WhitesonCollab/Data/GoogleDrive/',sampleID,'/gtf_wproducts.txt'), molecule=sampleID, invert = FALSE)$tempgtf
depth = read.table('S2M_plasmidsaurus_unmasked.sort.depth')
colnames(depth) = c('Reference', 'Position', 'Depth')
depth=subset(depth, Reference=='Bob2')
depth$Reference='Bob'
depth$NewPosition = 0
maxpos = max(depth$Position)
termgene = grep('terminase', tempgtf$product)[1]
termpos = tempgtf$start[termgene]
depth$NewPosition[depth$Position >= termpos] = depth$Position[depth$Position >= termpos] - termpos + 1
depth$NewPosition[depth$Position < termpos] = maxpos + depth$Position[depth$Position < termpos] - termpos + 1

rplot_S2M = ggplot(depth, aes(x=NewPosition, y=Depth, color=Reference))+
	geom_line(size=0.25)+
	xlab('')+
	scale_color_manual(values = c('dodgerblue'))+
	theme_bw()+
	theme(text = element_text(face='bold', size=15),
	legend.position='none')
	
pdf('../S2M_recombination_plot_Figure5.pdf', height=2.5, width=7)
genome_S2M/rplot_S2M
dev.off()


