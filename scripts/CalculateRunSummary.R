args <- commandArgs(trailingOnly = TRUE)

setwd(args[1])

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(rjson)

FixConcatNames <- function(a){
    a$Group.1 = as.character(a$Group.1)
    concatAlleleIndex = grep(",",a$Group.1)
    if ( length(concatAlleleIndex) > 0 ){
        for (i in 1:length(concatAlleleIndex)){
            alleles = strsplit(as.character(a$Group.1[concatAlleleIndex[i]]),",")[[1]]
            for (j in 1:length(alleles)){
                if (j == 1){
                    a$Group.1[concatAlleleIndex[i]] = alleles[j] #update the name in the original position
                }else if (is.na(match(alleles[j],a$Group.1))){
                    a=rbind(a,c(as.character(alleles[j]),as.character(a[concatAlleleIndex[i],-1])))
                }
            }
        }
    }
    return(a)
}

AccountMissedMarkers <- function(a, hotspotTable){
    a$Group.1 = as.character(a$Group.1)
    missedMarkers = which(!hotspotTable$V4 %in% a$Group.1)
    if ( length(missedMarkers) > 0 ){
        for (i in 1:length(missedMarkers)){
            a=rbind(a,c(as.character(hotspotTable$V4[missedMarkers[i]]),0,0,0,length(unique(a$Barcode)),0) )
        }
    }
    return (a)    
}



# Getting the variantCaller and coverageAnalysis ID, if it's 0, pick the latest run. IDs are incrementally created so the highest ID belongs to the latest run.
if (args[2] == 0 || is.na(args[2])){
    vids = list.files(paste0("../"),pattern="variantCaller_out.")
    vid = sort(as.numeric(gsub(".*\\.(\\d+)$",'\\1',vids)))[length(vids)]
}else{
    vid=args[2]
}

if (args[3] == 0 || is.na(args[3])){
    cids = list.files(paste0("../"),pattern="coverageAnalysis_out.")
    cid = sort(as.numeric(gsub(".*\\.(\\d+)$",'\\1',cids)))[length(cids)]
}else{
    cid=args[3]
}

if (args[4] == 0 || is.na(args[4])){
    fids = list.files(paste0("../"),pattern="HIDfDiagnosticsTool-RC2")
    fid = sort(as.numeric(gsub(".*\\.(\\d+)$",'\\1',fids)))[length(fids)]
}else{
    fid=args[4]
}

# Parse startplugin.json to get some constants about the run e.g. Run name etc.
startConfig = fromJSON(file = paste0("../variantCaller_out.",vid,"/startplugin.json"))
runName = startConfig$expmeta$run_name

hotspotFiles = regionFiles = c()

if ( !is.null(startConfig$pluginconfig$meta$targetloci) && startConfig$pluginconfig$meta$targetloci != "" ){
    hotspotFiles[1] = startConfig$pluginconfig$meta$targetloci
    regionFiles[1] = startConfig$pluginconfig$meta$targetregions
}else{
    if (!is.null(startConfig$plan$barcodedSamples)){
        for (i in 1:length(startConfig$plan$barcodedSamples)){
            barcodes = startConfig$plan$barcodedSamples[[i]]$barcodes
            for (j in 1:length(barcodes)){
                hotspotFiles[i] = startConfig$plan$barcodedSamples[[i]]$barcodeSampleInfo[barcodes[j]][[1]]$hotSpotRegionBedFile
                regionFiles[i] = startConfig$plan$barcodedSamples[[i]]$barcodeSampleInfo[barcodes[j]][[1]]$targetRegionBedFile
            }
        }
    }
}
##Get HotSpot File (Legacy)
#if (!is.null(startConfig$plan$regionfile)){
#    hotspotFile = startConfig$plan$regionfile
#}else if (!is.null(startConfig$planconfig$meta$targetloci)){
#    hotspotFile = startConfig$planconfig$meta$targetloci
#}
#hotspotTable=read.table(hotspotFile,sep="\t",skip=1)

hotspotFiles = unique(hotspotFiles)
for (i in 1:length(hotspotFiles)){
    if ( i == 1){
        hotspotTable=read.table(hotspotFiles[i],sep="\t",skip=1)
    }else{
        tmp = read.table(hotspotFiles[i],sep="\t",skip=1)
        hotspotTable = rbind(hotspotTable, tmp)
    }
}

hotspotTable$ref = gsub("REF=([ATGC]*);.*","\\1",hotspotTable$V7)
hotspotTable$obs = gsub(".*OBS=([ATGC]*);.*","\\1",hotspotTable$V7)

marker_number_snp=marker_number_del=marker_number_ins=marker_number_mnp=marker_number_complex=0;
for (i in 1:nrow(hotspotTable)){
    hotspotTable$obs[i] = sub(paste0("^",hotspotTable$ref[i]),"",hotspotTable$obs[i]) #Fix potential INS to be mis-categorized as COMPLEX
    hotspotTable$ref[i] = sub(paste0("^",hotspotTable$obs[i]),"",hotspotTable$ref[i]) #Fix potential DEL to be mis-categorized as COMPLEX
    
    lref = nchar(hotspotTable$ref[i])
    lobs = nchar(hotspotTable$obs[i])
    if ( lref == 1 &&  lobs == 1){
        marker_number_snp=marker_number_snp+1
    }else if (lref == 0){
        marker_number_ins=marker_number_ins+1
    }else if (lobs == 0){
        marker_number_del=marker_number_del+1
    }else if ( lref == lobs && lref > 1 ){
        marker_number_mnp=marker_number_mnp+1
    }else{
        marker_number_complex=marker_number_complex+1
    }
}

##Get Region File
#if (!is.null(startConfig$plan$bedfile)){
#    regionFile = startConfig$plan$bedfile
#}
#regionTable=read.table(regionFile,sep="\t",skip=1)

regionFiles = unique(regionFiles)
for (i in 1:length(regionFiles)){
    if ( i == 1){
        regionTable=read.table(regionFiles[i],sep="\t",skip=1)
    }else{
        tmp = read.table(regionFiles[i],sep="\t",skip=1)
        regionTable = rbind(regionTable, tmp)
    }
}


# We only include Hotspots in the analysis, so filter out the rest of the data
a=read.table(paste0("../variantCaller_out.",vid,"/",runName,".xls"),sep="\t",header=T)
a1=a[which(a$Allele.Source == "Hotspot"),]

#Marker Call Rate
b=aggregate(a1$Allele.Call,by=list(a1$Allele.Name),unlist(table))
b=as.data.frame(cbind(as.character(b[,1]),b[,2][,1:ncol(b[,2])]))
names(b)[1] = "Group.1"
btmp = setNames(data.frame(matrix(ncol = 5, nrow = nrow(b))), c("Group.1", "Heterozygous", "Homozygous", "Absent", "No Call"))
btmp[,match(names(b),names(btmp))] = b
btmp_no_data = which(!names(btmp) %in% names(b))
if (length(btmp_no_data) > 0 ){
    btmp[, btmp_no_data] = 0   
}
b=btmp

b$CallRate = apply(b,1,function(x){sum(as.numeric(x[-c(1,5)]))/sum(as.numeric(x[-1]))})

b[] <- lapply(b, as.character)

b = FixConcatNames(b)
b = AccountMissedMarkers(b, hotspotTable)
b$CallRate = as.numeric(b$CallRate)

#Sample Call Rate
c=aggregate(a1$Allele.Call,by=list(a1$Barcode),table)
c=as.data.frame(cbind(as.character(c[,1]),c[,2][,1:ncol(c[,2])]))
names(c)[1] = "Group.1"
ctmp = setNames(data.frame(matrix(ncol = 5, nrow = nrow(c))), c("Group.1", "Heterozygous", "Homozygous", "Absent", "No Call"))
ctmp[,match(names(c),names(ctmp))] = c
ctmp_no_data = which(!names(ctmp) %in% names(c))
if (length(ctmp_no_data) > 0 ){
    ctmp[, ctmp_no_data] = 0   
}
c=ctmp
c$CallRate = apply(c,1,function(x){sum(as.numeric(x[-c(1,5)]))/sum(as.numeric(x[-1]))})

#Coverage By Sample
d=aggregate(a1$Original.Coverage,by=list(a1$Barcode),FUN = function(x) tryCatch({t.test(x)$conf.int[1:2]},error = function(e) {c(NaN,NaN)}))
d1=aggregate(a1$Original.Coverage,by=list(a1$Barcode),mean)

#Coverage By Marker
e=aggregate(a1$Original.Coverage,by=list(a1$Allele.Name),FUN = function(x) tryCatch({t.test(x)$conf.int[1:2]},error = function(e) {c(NaN,NaN)}))
e1=aggregate(a1$Original.Coverage,by=list(a1$Allele.Name),mean)

e1 = FixConcatNames(e1)
e1 = AccountMissedMarkers(e1, hotspotTable)

out_prefix_name = paste0("Run_Summary_",vid,"_",cid,"_",fid)

cutoff = 0.95
number_markers=nrow(hotspotTable)
number_amplicons=nrow(regionTable)
samples_run=dim(c)[1]
samples_passed_threshold=sum(c$CallRate>=cutoff)
samples_call_rate=paste0(round(mean(c$CallRate)*100,0),"%")
markers_passed_threshold=sum(b$CallRate>=cutoff)
markers_call_rate=paste0(round(mean(b$CallRate)*100,0),"%")

index_markers_less_10 = b$CallRate>=0 & b$CallRate<0.1
index_samples_less_10 = c$CallRate>=0 & c$CallRate<0.1
markers_less_10=as.character(b$Group.1[index_markers_less_10])
samples_less_10=as.character(c$Group.1[index_samples_less_10])

poor_markers_tbl = cbind(e1[match(markers_less_10,e1$Group.1),],b[match(markers_less_10,b$Group.1),"CallRate"])
#poor_markers_tbl[,-1] = round(poor_markers_tbl[,-1],2)
write.table(poor_markers_tbl,paste0(out_prefix_name,"_poor_markers.txt"),sep="\t",quote=F,col.names=F,row.names=F)

poor_samples_tbl = cbind(d1[match(samples_less_10,d1$Group.1),],c[match(samples_less_10,c$Group.1),"CallRate"])
#poor_samples_tbl[,-1] = round(poor_samples_tbl[,-1],2)
write.table(poor_samples_tbl,paste0(out_prefix_name,"_poor_samples.txt"),sep="\t",quote=F,col.names=F,row.names=F)

marker_call_rate_coverage = cbind(e1[match(b$Group.1,e1$Group.1),],b[,"CallRate"])
names(marker_call_rate_coverage)=c("Allele Name","Mean Depth","Call Rate")

marker_call_rate_coverage2 = marker_call_rate_coverage
marker_call_rate_coverage2$"Call Rate" = 100*marker_call_rate_coverage2$"Call Rate"
marker_call_rate_coverage_melt <- melt(marker_call_rate_coverage2, id.var="Allele Name")
marker_call_rate_coverage_melt$value = as.numeric(marker_call_rate_coverage_melt$value)
write.table(marker_call_rate_coverage,paste0(out_prefix_name,"_marker_call_data.txt"),sep="\t",quote=F,col.names=T,row.names=F)

x=data.frame()
x[1,1] = number_markers
x$V1 = as.character(x$V1)
x[nrow(x)+1,1] = marker_number_complex
x[nrow(x)+1,1] = marker_number_del
x[nrow(x)+1,1] = marker_number_ins
x[nrow(x)+1,1] = marker_number_mnp
x[nrow(x)+1,1] = marker_number_snp
x[nrow(x)+1,1] = number_amplicons
x[nrow(x)+1,1] = samples_run
x[nrow(x)+1,1] = samples_passed_threshold
x[nrow(x)+1,1] = samples_call_rate
x[nrow(x)+1,1] = markers_passed_threshold
x[nrow(x)+1,1] = markers_call_rate
x[nrow(x)+1,1] = sum(b$CallRate==1)
x[nrow(x)+1,1] = sum(b$CallRate<1 & b$CallRate>=0.98)
x[nrow(x)+1,1] = sum(b$CallRate>=0.95 & b$CallRate<0.98)
x[nrow(x)+1,1] = sum(b$CallRate>=0.90 & b$CallRate<0.95)
x[nrow(x)+1,1] = sum(b$CallRate>=0.80 & b$CallRate<0.90)
x[nrow(x)+1,1] = sum(b$CallRate>=0.50 & b$CallRate<0.80)
x[nrow(x)+1,1] = sum(b$CallRate>=0.10 & b$CallRate<0.50)
x[nrow(x)+1,1] = sum(index_markers_less_10) - sum(b$CallRate==0)
x[nrow(x)+1,1] = sum(b$CallRate==0)
x[nrow(x)+1,1] = sum(c$CallRate==1)
x[nrow(x)+1,1] = sum(c$CallRate<1 & c$CallRate>=0.98)
x[nrow(x)+1,1] = sum(c$CallRate>=0.95 & c$CallRate<0.98)
x[nrow(x)+1,1] = sum(c$CallRate>=0.90 & c$CallRate<0.95)
x[nrow(x)+1,1] = sum(c$CallRate>=0.80 & c$CallRate<0.90)
x[nrow(x)+1,1] = sum(c$CallRate>=0.50 & c$CallRate<0.80)
x[nrow(x)+1,1] = sum(c$CallRate>=0.10 & c$CallRate<0.50)
x[nrow(x)+1,1] = sum(index_samples_less_10) - sum(c$CallRate==0)
x[nrow(x)+1,1] = sum(c$CallRate==0)
x[nrow(x)+1,1] = paste0(round(mean(a1$Original.Coverage),0),"X")

row.names(x)=c(
               "Number Markers",
               "Number COMPLEX Markers",
               "Number DEL Markers",
               "Number INS Markers",
               "Number MNP Markers",
               "Number SNP Markers",
               "Number Amplicons",
               "Samples Run",
               "Samples Passed >=95% CR",
               "Samples CR",
               "Markers Passed >=95% CR",
               "Markers CR",
               "# Markers = 100% CR",
               "# Markers 98% - 100% CR",
               "# Markers 95% - 98% CR",
               "# Markers 90% - 95% CR",
               "# Markers 80% - 90% CR",
               "# Markers 50% - 80% CR",
               "# Markers 10% - 50% CR",
               "# Markers 0% - 10% CR",
               "# Markers = 0% CR",
               "# Samples = 100% CR",
               "# Samples 98% - 100% CR",
               "# Samples 95% - 98% CR",
               "# Samples 90% - 95% CR",
               "# Samples 80% - 90% CR",
               "# Samples 50% - 80% CR",
               "# Samples 10% - 50% CR",
               "# Samples 0% - 10% CR",
               "# Samples = 0% CR",
               "Mean Coverage"
               )
write.table(x,paste0(out_prefix_name,".txt"),sep="\t",quote=F,col.names=F)

xls_files_cov = list.files(paste0("../coverageAnalysis_out.",cid),pattern="*bc_summary.xls$",full.names=T)
xls_file_cov=xls_files_cov[1]
acov=read.table(xls_file_cov,sep="\t",header=T)
acov$On.Target = gsub("%$","",acov$On.Target)
acov$Uniformity = gsub("%$","",acov$Uniformity)
acov_subset = subset(acov, select=c(Barcode.ID,Mean.Depth,On.Target, Uniformity))
acov_subset = cbind(acov_subset,"CallRate"=100*c[match(acov_subset$Barcode.ID,c$Group.1), "CallRate"])
names(acov_subset) = c('Barcode ID','Mean Depth', 'On Target', 'Uniformity', 'Call Rate')
acov_melt <- melt(acov_subset, id.var="Barcode ID")
acov_melt$value = as.numeric(acov_melt$value)
acov_melt$variable = factor(acov_melt$variable, levels=c("Call Rate","Uniformity","On Target","Mean Depth"))

max_labels = 96

countface <- 0
breaks_fun <- function(x) {
    countface <<- countface + 1L
    
    breaksy1 = seq(0, 100, by = 10)
    #breaksy2 = seq(0, max(x), by = if (max(x)/1000 > 1) ceiling((max(x)/10)/100)*100 else if (max(x)/600 > 1)  200 else if (max(x)/100 > 1) 100 else 10)
    breaksy2 = seq(0, max(x,na.rm=T), by = if (max(x,na.rm=T)/100 > 1) ceiling((max(x,na.rm=T)/10)/100)*100 else 10)
    switch(
        countface,
        breaksy1,
        breaksy1,
        breaksy1,
        breaksy2
    )
}

labels_fun <- function(x) {
    x[seq(0,length(x),2)]=""
    return(x)
}

p1=ggplot(acov_melt, aes(x = `Barcode ID`, y = value)) + geom_dotplot(binaxis='y', stackdir='center',aes(color = variable)) +
facet_grid(variable ~ ., scales = "free_y") + theme(legend.position = "none") +
    #theme_minimal() +
    ggtitle("Overall Sample Performance") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6))
    p11 = ggplot_build(p1)
    labels = p11$panel$ranges[[1]]$x.labels
    labels2=labels
    if (length(labels) > max_labels){
        labels2 = labels[seq(1,length(labels),by=ceiling(length(labels)/max_labels))]
    }
    
p1=p1+scale_x_discrete(breaks=labels2, labels=as.character(labels2)) +
scale_y_continuous(breaks = breaks_fun, labels=labels_fun) +
expand_limits(y = 0)
ggsave(paste0(out_prefix_name,"_1.png"),width=12, height=5)

labels = unique(b$Group.1)
labels2=labels
if (length(labels) > max_labels){
    labels2 = labels[seq(1,length(labels),by=ceiling(length(labels)/max_labels))]
}
p2=ggplot(a1,aes(Allele.Name,Frequency,colour = Allele.Call)) +
    geom_point(size = 1) +
    ggtitle("Allele Frequency Across Samples") +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6)) +
    scale_x_discrete(breaks=labels2, labels=as.character(labels2)) +
    scale_color_discrete(labels = c("Homozygous Ref", "Heterozygous", "Homozygous Variant", "No Call")) +
    labs(x = "Allele Name", y="Variant Allele Frequency (%)", colour = "Allele Call")
ggsave(paste0(out_prefix_name,"_2.png"),width=12, height=5)

breaks = data.frame()
for (i in 1:10){
    cut=10*i
    breaks[i,1] = (2*cut - 10)/2
    if (i == 1)
        breaks[i,2] = sum(100*b$CallRate>=cut-10 & 100*b$CallRate<=cut)
    else
        breaks[i,2] = sum(100*b$CallRate>cut-10 & 100*b$CallRate<=cut)
    breaks[i,3] = paste0((cut-10),"-",cut,"%")
}

p4=ggplot(breaks, aes(x=breaks$V1, y=breaks$V2, fill=factor(breaks$V3))) +
    ggtitle("Number of Markers wrt Call Rate") +
    guides(fill = guide_legend(reverse = TRUE)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=breaks$V2), vjust=-0.3, size=3.5) +
    scale_x_continuous(breaks = seq(0,100,10)) +
    scale_fill_brewer(palette="RdYlGn") +
    labs(x = "Call Rate Range", y="Number of Markers", fill="CR Range")
ggsave(paste0(out_prefix_name,"_3.png"),width=12, height=5)


breaks_sam = data.frame()
for (i in 1:10){
    cut=10*i
    breaks_sam[i,1] = (2*cut - 10)/2
    if (i == 1)
        breaks_sam[i,2] = sum(100*c$CallRate>=cut-10 & 100*c$CallRate<=cut)
    else
        breaks_sam[i,2] = sum(100*c$CallRate>cut-10 & 100*c$CallRate<=cut)
    breaks_sam[i,3] = paste0((cut-10),"-",cut,"%")
}

p5=ggplot(breaks_sam, aes(x=breaks_sam$V1, y=breaks_sam$V2, fill=factor(breaks_sam$V3))) +
    ggtitle("Number of Samples wrt Call Rate") +
    guides(fill = guide_legend(reverse = TRUE)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=breaks_sam$V2), vjust=-0.3, size=3.5) +
    scale_x_continuous(breaks = seq(0,100,10)) +
    scale_fill_brewer(palette="RdYlGn") +
    labs(x = "Call Rate Range", y="Number of Samples", fill="CR Range")
ggsave(paste0(out_prefix_name,"_4.png"),width=12, height=5)


p6=ggplot(marker_call_rate_coverage_melt, aes(x = `Allele Name`, y = value)) + geom_dotplot(binaxis='y', stackdir='center',aes(color = variable)) +
    ggtitle("Overall Marker Performance") +
    facet_grid(variable ~ ., scales = "free_y") + theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6))
    p11 = ggplot_build(p6)
    labels = p11$panel$ranges[[1]]$x.labels
    labels2=labels
    if (length(labels) > max_labels){
        labels2 = labels[seq(1,length(labels),by=ceiling(length(labels)/max_labels))]
    }
p6=p6+scale_x_discrete(breaks=labels2, labels=as.character(labels2))
ggsave(paste0(out_prefix_name,"_5.png"),width=12, height=5)