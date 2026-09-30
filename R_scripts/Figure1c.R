# Code to generate Figure 1c

library(gggenes)
library(stringr)
library(RColorBrewer)
library(dplyr)
library(ggplot2)

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

#Set color palette
darkpal = brewer.pal(8, 'Dark2')

#Use gbk2gtf to parse the original .gb files. Then set up the categories in Excel.
bobgtf = gbk2gtf('./Bob/Enterococcus phage vB_OCPT_Bob 2023-08-11_Bob_flyeAssembly_selectedContig_annotated.gb')
cargtf = gbk2gtf('./Car/Enterococcus phage vB_OCPT_Car 2023-08-11_Car_flyeAssembly_selectedContig_annotated.gb')
billgtf = gbk2gtf('./Bill/Enterococcus phage vB_OCPT_Bill 2023-08-25_Bill_flyeAssembly_selectedContig_annotated.gb')
sds1gtf = gbk2gtf('./SDS1/Enterococcus phage vB_OCPT_SDS1 2023-08-11_SDS1_flyeAssembly_selectedContig_annotated.gb')
ccs4gtf=gbk2gtf('./CCS4/Enterococcus phage vB_OCPT_CCS4 2023-08-25_CCS4_flyeAssembly_selectedContig_annotated.gb')

write.table(bobgtf, './Bob/gtf_wproducts.txt', row.names=F, quote=F, sep='\t')
write.table(cargtf, './Car/gtf_wproducts.txt', row.names=F, quote=F, sep='\t')
write.table(billgtf, './Bill/gtf_wproducts.txt', row.names=F, quote=F, sep='\t')
write.table(sds1gtf, './SDS1/gtf_wproducts.txt', row.names=F, quote=F, sep='\t')
write.table(ccs4gtf, './CCS4/gtf_wproducts.txt', row.names=F, quote=F, sep='\t')

#Update category designations to simplify to fewer general groups
#Note: Each genome is reoriented to begin from the first terminase. "termgene" below is the gene index of that terminase
#Bob
tempgtf_bob = read.table("./Bob/gtf_wproducts.txt", header=T, sep='\t')
tempgtf_bob$category[tempgtf_bob$category %in% c('baseplate')] = 'tail'
tempgtf_bob$category[tempgtf_bob$category %in% c('enzyme')] = 'other'
tempgtf_bob$category[tempgtf_bob$category %in% c('integration')] = 'recombination'
tempgtf_bob$category[is.na(tempgtf_bob$category)] = 'hypothetical'
tempgtf_bob$length = tempgtf_bob$end - tempgtf_bob$start
tempgtf_bob$sizeclass = factor(ifelse(tempgtf_bob$length > 500, 'large', 'small'))
termgene = 8
tempgtf_bob$newstart = tempgtf_bob$newend = 0
tempgtf_bob$newstart[termgene:nrow(tempgtf_bob)] = tempgtf_bob$start[termgene:nrow(tempgtf_bob)] - tempgtf_bob$start[termgene] + 1
endshift = tempgtf_bob$newstart[nrow(tempgtf_bob)] + tempgtf_bob$length[nrow(tempgtf_bob)]
tempgtf_bob$newstart[1:(termgene - 1)] = tempgtf_bob$start[1:(termgene - 1)] + endshift
tempgtf_bob$newend = tempgtf_bob$newstart + tempgtf_bob$length

#Car
tempgtf_car = read.csv("./Car/gtf_wproducts.txt", header=T, sep='\t')
tempgtf_car$category[tempgtf_car$category %in% c('baseplate')] = 'tail'
tempgtf_car$category[tempgtf_car$category %in% c('enzyme')] = 'other'
tempgtf_car$category[tempgtf_car$category %in% c('integration')] = 'recombination'
tempgtf_car$category[is.na(tempgtf_car$category)] = 'hypothetical'
tempgtf_car$length = tempgtf_car$end - tempgtf_car$start
tempgtf_car$sizeclass = factor(ifelse(tempgtf_car$length > 500, 'large', 'small'))

#Car's nucleotide sequence was reversed relative to the orientation for Bob. Correct this to align them better
#Step 1. Reverse the row/gene order
#Step 2. Reset row 1 to start at 1. Adjust end to be Start + length
#Step 3. For successive genes, the new start should be the difference between the original end and the start of next position

newgtf_car = tempgtf_car[nrow(tempgtf_car):1,]
newgtf_car$start2 = newgtf_car$end2 = 0
newgtf_car$start2[1] = 1
newgtf_car$end2[1] = newgtf_car$start2[1]+newgtf_car$length[1]
for(i in 2:nrow(newgtf_car)){
	newgtf_car$start2[i] = newgtf_car$end2[i-1] + (newgtf_car$start[i-1] - newgtf_car$end[i])
	newgtf_car$end2[i] = newgtf_car$start2[i] + newgtf_car$length[i]
}
termgene = grep('terminase',newgtf_car$product)[1]
newgtf_car$newstart = newgtf_car$newend = 0
newgtf_car$newstart[termgene:nrow(newgtf_car)] = newgtf_car$start2[termgene:nrow(newgtf_car)] - newgtf_car$start2[termgene] + 1
endshift = newgtf_car$newstart[nrow(newgtf_car)] + newgtf_car$length[nrow(newgtf_car)]
newgtf_car$newstart[1:(termgene - 1)] = newgtf_car$start2[1:(termgene - 1)] + endshift
newgtf_car$newend = newgtf_car$newstart + newgtf_car$length
newgtf_car$orientation = 1 - newgtf_car$orientation  #key step to ensure the orientation is forward.

#Bill
tempgtf_bill = read.csv('./Bill/gtf_wproducts.txt', sep='\t')
tempgtf_bill$length = tempgtf_bill$end - tempgtf_bill$start
tempgtf_bill$sizeclass = factor(ifelse(tempgtf_bill$length > 500, 'large', 'small'))
termgene = 93
tempgtf_bill$newstart = tempgtf_bill$newend = 0
tempgtf_bill$newstart[termgene:nrow(tempgtf_bill)] = tempgtf_bill$start[termgene:nrow(tempgtf_bill)] - tempgtf_bill$start[termgene] + 1
endshift = tempgtf_bill$newstart[nrow(tempgtf_bill)] + tempgtf_bill$length[nrow(tempgtf_bill)]
tempgtf_bill$newstart[1:(termgene - 1)] = tempgtf_bill$start[1:(termgene - 1)] + endshift
tempgtf_bill$newend = tempgtf_bill$newstart + tempgtf_bill$length

#combined
sisters_gtf = rbind(tempgtf_bob, newgtf_car[,colnames(tempgtf_bob)])
sisters_gtf = rbind(sisters_gtf, tempgtf_bill)

#SDS1
tempgtf_sds1 = read.csv('./SDS1/gtf_wproducts.txt', sep='\t')
tempgtf_sds1$length = tempgtf_sds1$end - tempgtf_sds1$start
tempgtf_sds1$sizeclass = factor(ifelse(tempgtf_sds1$length > 500, 'large', 'small'))
termgene = 14
tempgtf_sds1$newstart = tempgtf_sds1$newend = 0
tempgtf_sds1$newstart[termgene:nrow(tempgtf_sds1)] = tempgtf_sds1$start[termgene:nrow(tempgtf_sds1)] - tempgtf_sds1$start[termgene] + 1
endshift = tempgtf_sds1$newstart[nrow(tempgtf_sds1)] + tempgtf_sds1$length[nrow(tempgtf_sds1)]
tempgtf_sds1$newstart[1:(termgene - 1)] = tempgtf_sds1$start[1:(termgene - 1)] + endshift
tempgtf_sds1$newend = tempgtf_sds1$newstart + tempgtf_sds1$length
newgtf_sds1 = tempgtf_sds1[nrow(tempgtf_sds1):1,]
newgtf_sds1$start2 = newgtf_sds1$end2 = 0
newgtf_sds1$start2[1] = 1
newgtf_sds1$end2[1] = newgtf_sds1$start2[1]+newgtf_sds1$length[1]
for(i in 2:nrow(newgtf_sds1)){
	newgtf_sds1$start2[i] = newgtf_sds1$end2[i-1] + (newgtf_sds1$start[i-1] - newgtf_sds1$end[i])
	newgtf_sds1$end2[i] = newgtf_sds1$start2[i] + newgtf_sds1$length[i]
}
termgene = grep('terminase',newgtf_sds1$product)[1]
newgtf_sds1$newstart = newgtf_sds1$newend = 0
newgtf_sds1$newstart[termgene:nrow(newgtf_sds1)] = newgtf_sds1$start2[termgene:nrow(newgtf_sds1)] - newgtf_sds1$start2[termgene] + 1
endshift = newgtf_sds1$newstart[nrow(newgtf_sds1)] + newgtf_sds1$length[nrow(newgtf_sds1)]
newgtf_sds1$newstart[1:(termgene - 1)] = newgtf_sds1$start2[1:(termgene - 1)] + endshift
newgtf_sds1$newend = newgtf_sds1$newstart + newgtf_sds1$length
newgtf_sds1$orientation = 1 - newgtf_sds1$orientation

#CCS4
tempgtf_ccs4 = read.csv('./CCS4/gtf_wproducts.txt', header=T, sep='\t')
tempgtf_ccs4$length = tempgtf_ccs4$end - tempgtf_ccs4$start
tempgtf_ccs4$sizeclass = factor(ifelse(tempgtf_ccs4$length > 500, 'large', 'small'))
termgene = 14
tempgtf_ccs4$newstart = tempgtf_ccs4$newend = 0
tempgtf_ccs4$newstart[termgene:nrow(tempgtf_ccs4)] = tempgtf_ccs4$start[termgene:nrow(tempgtf_ccs4)] - tempgtf_ccs4$start[termgene] + 1
endshift = tempgtf_ccs4$newstart[nrow(tempgtf_ccs4)] + tempgtf_ccs4$length[nrow(tempgtf_ccs4)]
tempgtf_ccs4$newstart[1:(termgene - 1)] = tempgtf_ccs4$start[1:(termgene - 1)] + endshift
tempgtf_ccs4$newend = tempgtf_ccs4$newstart + tempgtf_ccs4$length

cousins_gtf = rbind(tempgtf_bob, tempgtf_ccs4)
cousins_gtf = rbind(cousins_gtf, newgtf_sds1[,colnames(cousins_gtf)])

all_gtf = distinct(rbind(sisters_gtf, cousins_gtf))
all_gtf$molecule = str_extract(all_gtf$molecule, '_[BCS].+') %>% str_replace('_','')
all_gtf$molecule = factor(all_gtf$molecule, levels = c('Bill','Car','Bob','CCS4','SDS1'))
all_genomes_plot = ggplot(all_gtf, aes(xmin=newstart, xmax=newend, y=molecule, forward=orientation, fill = category,size=sizeclass))+
  facet_wrap(vars(molecule), ncol=1, scale='free_y')+
  geom_gene_arrow(arrow_body_height = grid::unit(2, 'mm'), arrowhead_width = grid::unit(2,'mm'))+
  scale_fill_manual(values=c(darkpal[1:2],'grey',darkpal[3:6]))+
  scale_size_manual(values=c(0.4,0.1))+
  guides(size = 'none')+
  ylab('')+
  labs(fill = 'Product')+
  theme_genes()+
  theme(legend.position = 'bottom',
  text = element_text(face='bold', size=15))
pdf('../all_genome_maps_annotated_Figure1c.pdf', height=3, width=8)
print(all_genomes_plot)
dev.off()

