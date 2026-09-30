# Code to generate Figure 2C-D
library(ggplot2)
library(stringr)
library(khroma)

#Set color palette
colrange = as.character(color('bright')(6))

#Phage coverage over time
setwd("/mnt/d/Research/WhitesonCollab/Data/RawReads")

#the directory consists of folders in the format sisters_tX and cousins_tX, where X can be 0, 1, 2, 3, 4, 5, 10, 15, 20, representing the time points.
poplist = dir()
poplist = poplist[-grep('trimmed', poplist)] #remove a directory that is not related to time
timelist = str_extract(poplist, 't[0-9]+')
covtimemat = character()

for(i in 1:length(poplist)){
	setwd(poplist[i])
	filelist = dir(pattern='_on_combo.*.sort.depth')
	#tempsamples = str_extract(filelist, '[CS][123][IM]')
	temptime = timelist[i]
	for(j in 1:length(filelist)){
		tempdepth = read.table(filelist[j])
		temprefs = unique(tempdepth$V1)
		tempres = character()
		for(k in 1:length(temprefs)){
			tempres[1] = filelist[j]
			tempres[2] = temptime
			tempres[3] = temprefs[k]
			tempres[4] = mean(tempdepth$V3[tempdepth$V1==temprefs[k]])
			covtimemat = rbind(covtimemat, tempres)
		}
	}
	setwd('..')
}
rownames(covtimemat) = NULL
colnames(covtimemat) = c('Sample','Time','Reference','MeanCov')

#remove the Bob1 results (The original assembly had an extraneous contig)
covtimemat = as.data.frame(covtimemat)
covtimemat = subset(covtimemat, Reference!='Bob1')

#rename Bob2 to Bob
covtimemat$Reference[covtimemat$Reference=='Bob2'] = 'Bob'
covtimemat[,1] = str_replace_all(covtimemat[,1], 'Phage_Cousins_R','C')
covtimemat[,1] = str_replace_all(covtimemat[,1], 'Phage_Sisters_R','S')
covtimemat[,1] = str_extract(covtimemat[,1],'[CS][123].')
covtimemat[,1] = str_replace_all(covtimemat[,1],'_','')
covtimemat[,1] = factor(covtimemat[,1], levels = c('C1','C2','C3','C1M','C2M','C3M','C1I','C2I','C3I','S1','S2','S3','S1M','S2M','S3M','S1I','S2I','S3I'))
covtimemat[,2] = factor(covtimemat[,2], levels = paste0('t',c(0,1,2,3,4,5,10,15,20)))
covtimemat[,3] = factor(covtimemat[,3], levels = c('Bob','Car','Bill','CCS4','SDS1'))
covtimemat[,4] = as.numeric(covtimemat[,4])

#To get the t0 into the other plots, let's make copies of each, one relabeled with I and one with M, since it is the same
covtimemat$Sample[covtimemat$Sample%in%c('C1','C2','C3','S1','S2','S3')] = paste0(covtimemat$Sample[covtimemat$Sample%in%c('C1','C2','C3','S1','S2','S3')],'I')
temp = subset(covtimemat, Time=='t0')
temp$Sample = str_replace_all(temp$Sample, 'I','M')
covtimemat = rbind(covtimemat, temp)

#Last, correct the label scheme to match the final version used for the manuscript. C -> D. I -> P
covtimemat$Sample = as.character(covtimemat$Sample)
covtimemat$Sample = str_replace(covtimemat$Sample, 'C', 'D')
covtimemat$Sample = str_replace(covtimemat$Sample, 'I', 'P')
covtimemat$Sample = factor(covtimemat$Sample, levels = c('S1M','S2M','S3M','S1P','S2P','S3P','D1M','D2M','D3M','D1P','D2P','D3P'))

#Last, each coverage estimate is adjusted based on the size of each phage reference genome.
#Pull in references and get the genome lengths
bobref = readLines("/mnt/d/Research/WhitesonCollab/Data/GoogleDrive/Bob/2023-08-11_Bob_flyeAssembly_selectedContig.fasta")
bobref=c(">Bob",paste(bobref[2:length(bobref)], collapse=""))

carref = readLines("/mnt/d/Research/WhitesonCollab/Data/GoogleDrive/Car/2023-08-11_Car_flyeAssembly_selectedContig.fasta")
carref = c(">Car", paste(carref[2:length(carref)], collapse=""))

billref = readLines("../GoogleDrive/Bill/2023-08-25_Bill_flyeAssembly_selectedContig.fasta")
billref = c(">Bill", paste(billref[2:length(billref)], collapse=""))

ccs4ref = readLines("../GoogleDrive/CCS4/2023-08-25_CCS4_flyeAssembly_selectedContig.fasta")
ccs4ref = c(">CCS4", paste(ccs4ref[2:length(ccs4ref)], collapse=''))

sds1ref = readLines("../GoogleDrive/SDS1/2023-08-11_SDS1_flyeAssembly_selectedContig.fasta")
sds1ref = c(">SDS1", paste(sds1ref[2:length(sds1ref)], collapse=''))

sizevec = c(nchar(bobref[2]), nchar(carref[2]), nchar(billref[2]), nchar(ccs4ref[2]), nchar(sds1ref[2]))
names(sizevec) = c('Bob','Car','Bill','CCS4','SDS1')

covtimemat$NormCov = 150*covtimemat$MeanCov/sizevec[covtimemat$Reference] #multiplied by 150 to rescale by the read length
p1=ggplot(covtimemat, aes(x=Time, y=NormCov, fill=Reference))+
 facet_wrap(vars(Sample), ncol=3)+
 ylab('Relative Coverage (Normalized by genome size)')+
 xlab('Experimental Cycle')+
 geom_bar(stat='identity', position='fill')+
 scale_fill_manual(values=colrange)+
 theme_bw()+
 theme(text=element_text(face='bold', size=12), axis.text = element_text(size=8))

pdf('../coverage_summary_by_time_wt0_Figure2cd_normalized.pdf',height=6, width=7.25)
print(p1)
dev.off()
