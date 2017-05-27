##########################################################################
#FSig-SV.R
#Khurana Lab, Weill Cornell  Medical Collage

#randomly shifts structural variant set 1000 times
#shifts within arm and avoids gaps
#pericentric SVs remain pericentric 

#for each SV, first find arm
#choose a random point on the arm (in "allowed region" i.e. not in a gap) 
#to be the new start of the SV
#add width to find new end 

#then overlap the observed and permuted SV against CDS or promoter or enhancer to compute p-values

#NOTE- This is a an extremely time-consuming process
#########################################################################

##!/usr/bin/env Rscript

#load libraries
args<-commandArgs(TRUE)
library(data.table)
library(doParallel)
library(ggplot2)

options(scipen=999)

#load SV input data set

input_SV = fread(args[1])
names(input_SV) = c("chr","start","end","sv_type","sample")

#Chromosome arm lengths
chr_arm_lengths = readRDS("./input/chr_arm_lengths.rds")
input_SV = merge(input_SV,chr_arm_lengths,by="chr")

#determine arm of each SV, add column to SV data table
input_SV$arm = ifelse(input_SV$end<=input_SV$p_end,'p',ifelse(input_SV$start>=input_SV$p_end,'q',ifelse(input_SV$start<=input_SV$p_end & input_SV$end>=input_SV$p_end,'b','z')))

#allowed regions/arm lengths/widths, split by chr and arm  
genome_chr_regions = readRDS("./input/chr_regions_arm_width.rds")
chr_regions = split(genome_chr_regions, list(genome_chr_regions$chr,genome_chr_regions$arm2))

#shuffles the SVs; restricts to arm/avoids gaps
random_shuffle = function(chr,arm,start,end,p_end,q_end){
        new_start = c()
        new_end = c()
        
        #index for allowed regions
        reg = paste(chr,arm,sep=".") #e.g. chr5.q
        
        #indices for p-allowed regions, q-allowed regions
        #these are only used when SV arm is 'b' (pericentric)
        reg_p = paste(chr,"p",sep=".")
        reg_q = paste(chr,"q",sep=".")

        #allowed regions
        allowed = (chr_regions[[reg]])[,c("start","end")]
        allowed_p = (chr_regions[[reg_p]])[,c("start","end")]
        allowed_q = (chr_regions[[reg_q]])[,c("start","end")]
        
        #arbitrary value starts loops
        y = 1000000000
 
        #used to calculate width in loops
        s_orig = start
        e_orig = end
        
        if(arm!='b'){
                #continue until SV end is within arm/in an allowed region
                arm_end = ifelse(arm=='p',p_end,q_end)
                while((start+y >= arm_end) & length(which(allowed$start <= (start+y) & allowed$end >= (start+y))) == 0){
                        p_allowed = allowed[sample(nrow(allowed), 1),]
                        rand_range = c(p_allowed[[1]],p_allowed[[2]])
                        rand = round(runif(1,rand_range[1],rand_range[2]))
                        y = e_orig-s_orig #set y to width of SV 
                        start = rand
                    }
                new_start = c(new_start,start)
                new_end = c(new_end,start+y)
                }
        #if arm is b, SV end should be in an allowed region on q
        #SV start should be in an allowed region on p
        if(arm=='b'){
                while(start+y >= q_end & length(which(allowed_q$start <= (start+y) & allowed_q$end >= (start+y))) == 0){
                b_allowed = allowed_p[sample(nrow(allowed_p), 1),]
                rand_range = c(b_allowed[[1]],b_allowed[[2]])
                rand_2 = round(runif(1,rand_range[1],rand_range[2]))   
                y = e_orig - s_orig
                start = rand_2
                }
                new_start = c(new_start,start)
                new_end = c(new_end,start+y)
        }
        #this data frame has the same number of rows as the original SV set
        #the columns will be added to the original set 
        return(data.frame(new_start,new_end))
}

#CDS
cds = fread("./input/gencode.v16.cds.bed")
names(cds) = c("chr","start","end","gene")
cds = setkey(cds,chr,start,end)

#get the samples recurrence for each genes affected in the observed SV data set
input_SV = setkey(input_SV, chr, start, end) #overlap by chr,start,end
cds = setkey(cds,chr,start,end)
print("observed")
print("starting overlap")
x = foverlaps(cds,input_SV,nomatch=0) #exclude elements w/ no match#print(x)
print("finished overlap")
x = x[,.(gene,sample)] #select columns by name 

sample_recurrence=x[, .(num_samples = length(unique(sample))), by = gene]
sample_recurrence$ran_count = c(rep(0,length(sample_recurrence$gene)))

#promoter
pro = fread("./input/gencode.v16.promoter.bed")
names(pro) = c("chr","start","end","gene")
pro = setkey(pro,chr,start,end)

#get the samples recurrence for each genes affected in the observed SV data set
input_SV = setkey(input_SV, chr, start, end) #overlap by chr,start,end
pro = setkey(pro,chr,start,end)
print("observed")
print("starting overlap")
x = foverlaps(pro,input_SV,nomatch=0) #exclude elements w/ no match#print(x)
print("finished overlap")
x = x[,.(gene,sample)] #select columns by name 

sample_pro_recurrence=x[, .(num_samples = length(unique(sample))), by = gene]
sample_pro_recurrence$ran_count = c(rep(0,length(sample_pro_recurrence$gene)))

#enhancer
drm = fread("./input/encode.v16.drm.bed")
names(drm) = c("chr","start","end","gene")
drm = setkey(drm,chr,start,end)

#get the samples recurrence for each genes affected in the observed SV data set
input_SV = setkey(input_SV, chr, start, end) #overlap by chr,start,end
drm = setkey(drm,chr,start,end)
print("observed")
print("starting overlap")
x = foverlaps(drm,input_SV,nomatch=0) #exclude elements w/ no match#print(x)
print("finished overlap")
x = x[,.(gene,sample)] #select columns by name 

sample_enh_recurrence=x[, .(num_samples = length(unique(sample))), by = gene]
sample_enh_recurrence$ran_count = c(rep(0,length(sample_enh_recurrence$gene)))


print("permutation")
#generate random replicates of SV data set (x1000)
iter<-args[2]
for(i in 1:iter){
print(i)

#applies to each row (avoid for loop!)
v=input_SV[,random_shuffle(chr,arm,start,end,p_end,q_end), by = 1:nrow(input_SV)]

#combine new columns with old SV set
input_SVi = cbind(input_SV,v)

#get the samples recurrence for each genes affected in the random SV data set and compares samples recurrence to observed counts
#if samples affected >= original samples affected, add 1 to count for that gene
    input_SVi = setkey(input_SVi,chr,new_start,new_end)
    print("starting overlap")
    x = foverlaps(cds,input_SVi,nomatch=0)
    print("finished overlap")
    x = x[,.(gene,sample)]
    y2=x[, .(num_samples = length(unique(sample))), by = gene]
    z = merge(x=sample_recurrence,y=y2,by.x = "gene", by.y = "gene", all.x = TRUE)
    z[num_samples.y >= num_samples.x]$ran_count = z[num_samples.y >= num_samples.x]$ran_count+1
    sample_recurrence = z[ , .(gene, num_samples.x,ran_count)]
    names(sample_recurrence)[names(sample_recurrence)=="num_samples.x"]="num_samples"

    x = foverlaps(pro,input_SVi,nomatch=0)
    print("finished overlap")
    x = x[,.(gene,sample)]
    y2=x[, .(num_samples = length(unique(sample))), by = gene]
    z = merge(x=sample_pro_recurrence,y=y2,by.x = "gene", by.y = "gene", all.x = TRUE)
    z[num_samples.y >= num_samples.x]$ran_count = z[num_samples.y >= num_samples.x]$ran_count+1
    sample_pro_recurrence = z[ , .(gene, num_samples.x,ran_count)]
    names(sample_pro_recurrence)[names(sample_pro_recurrence)=="num_samples.x"]="num_samples"

    x = foverlaps(drm,input_SVi,nomatch=0)
    print("finished overlap")
    x = x[,.(gene,sample)]
    y2=x[, .(num_samples = length(unique(sample))), by = gene]
    z = merge(x=sample_enh_recurrence,y=y2,by.x = "gene", by.y = "gene", all.x = TRUE)
    z[num_samples.y >= num_samples.x]$ran_count = z[num_samples.y >= num_samples.x]$ran_count+1
    sample_enh_recurrence = z[ , .(gene, num_samples.x,ran_count)]
    names(sample_enh_recurrence)[names(sample_enh_recurrence)=="num_samples.x"]="num_samples"

}

sample_recurrence$p.value <- (sample_recurrence$ran_count)/iter
sample_recurrence$p.corretion <- p.adjust(sample_recurrence$p.value, "BH", length(sample_recurrence$p.value))
sample_recurrence_significant <- subset(sample_recurrence, sample_recurrence$p.corretion <= 0.01)

sample_pro_recurrence$p.value <- (sample_pro_recurrence$ran_count)/iter
sample_pro_recurrence$p.corretion <- p.adjust(sample_pro_recurrence$p.value, "BH", length(sample_pro_recurrence$p.value))
sample_pro_recurrence_significant <- subset(sample_pro_recurrence, sample_pro_recurrence$p.corretion <= 0.01)

sample_enh_recurrence$p.value <- (sample_enh_recurrence$ran_count)/iter
sample_enh_recurrence$p.corretion <- p.adjust(sample_enh_recurrence$p.value, "BH", length(sample_enh_recurrence$p.value))
sample_enh_recurrence_significant <- subset(sample_enh_recurrence, sample_enh_recurrence$p.corretion <= 0.01)

#save
print("saving")
write.csv(sample_recurrence, file="./output/Gene_affected_CDS_pvalue.csv", row.name=TRUE)  
write.csv(sample_recurrence_significant, file="./output/Gene_significant_affected_CDS.csv", row.name=TRUE)

write.csv(sample_pro_recurrence, file="./output/Gene_affected_PRO_pvalue.csv", row.name=TRUE)
write.csv(sample_pro_recurrence_significant, file="./output/Gene_significant_affected_PRO.csv", row.name=TRUE)

write.csv(sample_enh_recurrence, file="./output/Gene_affected_ENH_pvalue.csv", row.name=TRUE)
write.csv(sample_enh_recurrence_significant, file="./output/Gene_significant_affected_ENH.csv", row.name=TRUE)

  
