# Goes through .paf of reads to an assembly and outputs assignment of reads to regions

library(dplyr)
library(data.table)

# Read in defined regions of your reference genome
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
regions <- read.table("chrom_regions_3D7REF.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE) 

# paf input
hits <- fread("ANC.30kb-3D7REF.paf", select = c(1:14),fill=T) ## fill=T b/c some lines have end cols missing 
names(hits) = c("query.name","query.length" ,"query.start","query.end" ,"strand", "target.name","target.length" ,"target.start" ,"target.end","matching.bases", "alignment.length","mapping.quality", "type.alignment" , "N.minimizers")


### FILTER OUT LOW-QUALITY ALIGNMENTS 
# remove reads with low mapping quality 
hist(hits$mapping.quality)
hits<- hits%>%filter(mapping.quality > 10)
# remove short alignments
quantile(hits$alignment.length)
hits <- hits %>% filter(alignment.length > 10000)
# remove short reads 
hits <- hits %>% filter(query.length > 30000)
# add a %id column and filter by identity 
hits$identity <- hits$matching.bases/hits$alignment.length
hist(hits$identity)
hits <- hits %>% filter(identity > 0.65)

### sort by read, then alignment length 
hits <- hits[with(hits, order(query.name,alignment.length)), ] # most reads have one hit. Those with 2 tend to be triangles (inverted duplications)

### save names of reads, and create empty output table
reads <- unique(hits$query.name) 




### FOR EACH READ, FIND ITS REGION.
out <- data.frame(matrix(ncol = 2, nrow = length(reads)))
colnames(out) <- c("read","assigned_region")
out$read <- reads

for (i in 1:length(reads)){ # takes <10 sECONDS for 1000 reads
  # get the read's hits 
  rhits <-  hits %>% filter (query.name == reads[[i]])

  # if there are at least two hits, keep two with largest bitscore (should already be in order)
  # if there are fewer than two hits, will address with else statement later
  if (dim(rhits)[[1]] >=2){
        'proceed'
        rhits2 <- rhits[1:2,] # get first two hits only 
        
        # only proceed if both hits are from same chrom OR one is much (4x) better [longer] than the other
        if (rhits2[1,6] == rhits2[2,6] || rhits2[1,11] > 3.9*rhits2[2,11]){
          'proceed'
          
          # and that chrom is in regions
          if(rhits2[1,6] %in% regions$chr){
                      regionschr <- subset(regions,regions$chr==as.character(rhits[1,6]))
                      
                      # does the top hit overlap with anything in regionschr? #
                      # get the coordinates of the best hit 
                      besthit_start = min(rhits2[1,8],rhits2[1,9]) 
                      besthit_end = max(rhits2[1,8],rhits2[1,9])
                      
                      #add a column to regionschr for 'frac', the fraction of the besthit read contained in each possible region (row)
                      regionschr$frac = NA
                      # calculate the frac from besthit start and end 
                      for (row in 1:(dim(regionschr)[[1]])){
                        region_start = regionschr[row,2]
                        region_end=regionschr[row,3]
                        if (region_start>besthit_end || besthit_start>region_end) {
                          frac = 0} else {
                            os = max(besthit_start,region_start)
                            oe = min (besthit_end,region_end)
                            frac = (oe-os)/(besthit_end-besthit_start)
                          }
                        regionschr[row,5] = frac
                      }
                      
                      # only proceed if frac!=0 
                      match <- regionschr %>% filter (frac >0)
                      
                      # allow >1 row as long as region is the same 
                      if ( length(unique(as.list(match$region))) == 1  ){
                        out[i,2] = match[1,4]
                      }
                      else {out[i,2] = NA}
                      
          } # end is chrom in regions test 
        } # end from same chrom test 
        
   ## else if there aren't at least two hits 
  } else if (dim(rhits)[[1]] == 1){
      'proceed'
      
      # ..if it's big and on of the chroms we're considering
      if (rhits$alignment.length > 10000 && rhits$target.name %in% regions$chr){ 
        'proceed'
        regionschr <- subset(regions,regions$chr==as.character(rhits[1,6]))
        
        # does the top hit overlap with anything in regionschr? #
        # get the coordinates of the best hit 
        besthit_start = min(rhits[1,8],rhits[1,9]) 
        besthit_end = max(rhits[1,8],rhits[1,9])
        
        #add a column to regionschr for 'frac', the fraction of the besthit read contained in each possible region (row)
        regionschr$frac = NA
        # calculate the frac from besthit start and end 
        for (row in 1:(dim(regionschr)[[1]])){
          region_start = regionschr[row,2]
          region_end=regionschr[row,3]
          if (region_start>besthit_end || besthit_start>region_end) {
            frac = 0} else {
              os = max(besthit_start,region_start)
              oe = min (besthit_end,region_end)
              frac = (oe-os)/(besthit_end-besthit_start)
            }
          regionschr[row,5] = frac
        }
        
        # only save if there is one row where frac!=0 
        match <- regionschr %>% filter (frac >0)
        
        # allow >1 row as long as region is the same 
        if ( length(unique(as.list(match$region))) == 1  ){
          out[i,2] = match[1,4]
        }
      }
    }
   
}


###
goodout<-subset(out,is.na(out$assigned_region)==F)

summary <- goodout %>%
  group_by(assigned_region) %>%
  summarise(n = n())

write.table(goodout,file="ANC-REF.Chrom-RegionAssignments.csv",sep=",",col.names=TRUE,row.names=FALSE)
