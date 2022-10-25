### This script uses a .paf alignment between reads and a reference
### to identify reads that are large inverted duplicaations 

library(data.table)
library(dplyr)

setwd("/Users/emily/Desktop/forgithub/inverted-duplications") 
paf <- fread("ANC-psDNA-3D7ref.paf", select = c(1:14),fill=T) ## fill=T b/c some lines have end cols missing 

# how many reads to start?
length(unique(paf$query.name)) 

# remove alignments with low mapping quality 
hist(paf$mapping.quality)
paf <- subset(paf,paf$mapping.quality >=10)

# add a %id column
paf$identity <- paf$matching.bases/paf$alignment.length
hist(paf$identity) # bimodal 


### pull out reads with 2 alignments -- one high quality, one low quality  
triangles = c() # save the names of the reads
putative <- data.frame(matrix(ncol = 0, nrow = 0))
reads = unique(paf$query.name)

## for every read,
for (i in reads){
  ## get its alignments
  hits <- paf %>% filter(query.name == i)
  ## pick the 2 longest alignments 
  if (dim(hits)[[1]] > 2) {
    hits <- hits[with(hits, order(-alignment.length)), ]
    hits <- hits[1:2,] }
  ## has to have at least 2 to be considered
  if (dim(hits)[[1]] == 2) {
    ## they should be on same chrom; one + strand, one - strand; and at least one with reasonable identity 
    if (hits[1,6] == hits[2,6] && hits[1,5] != hits[2,5] && max(hits[1,15],hits[2,15])>=0.5 ) { 
      
      ## the query bits should not overlap much? ->  most don't, and those that do still have good overlap of target bits; so I'm not going to use this filter 
      ## the target bits SHOULD overlap 
      start1 = as.numeric(hits[1,8])
      end1 = as.numeric(hits[1,9])
      start2 = as.numeric(hits[2,8])
      end2 = as.numeric(hits[2,9])
      
      if (start2 > end1 || start1 > end2) {pass = FALSE # if they don't overlap
      } else { # how much does second overlap with the first?
        
        overlapstart = max(start1,start2)
        overlapend = min(end1,end2)
        overlaplen = as.numeric(overlapend - overlapstart)
        
        interval1len = as.numeric(hits[1,9] - hits[1,8])
        interval2len = as.numeric(hits[2,9] - hits[2,8])
        
        ratio1 = overlaplen/interval1len
        ratio2 = overlaplen/interval2len
        
        # overlaps = c(overlaps, min(ratio1,ratio2))
        
        if (min(ratio1,ratio2) > .20) {  ## if the overlap is sufficient, do the length test 
          #overlapreads = c(overlapreads,i)
          len1 = as.numeric(hits[1,11]); len2=as.numeric(hits[2,11])
          if ( min(len1,len2) >= 0.15*max(len1,len2)  ) {
            triangles=c(triangles,i)
            putative = rbind(putative,hits) } ## save and be done!
          
        } # end min ratio loop
      } # end else/making sure they overlap  
      
      
      ## when the overlap is low, one alignment is typically a lot shorter than the other 
      ## so maybe do a simpler ask: that one alignment is at least 15% of the other
      # -> looking through the first 100 reads this produces, the arms partially overlap in 100% of cases. so this seems fine 
      # -> but actually, for gDNA, there is 1 read in first 16 where there is 0 overlap
      
      
    } # end if on same chrom and strand
  } # end if dim(hits)==2
  
} # end for every read


length(triangles) 


# identity on each side 
putative$mean = (putative$query.end+putative$query.start)/2
putative$x = putative$mean/putative$query.length
putative$identity = putative$identity*100
left <- subset(putative,putative$x <0.5)
right <- subset(putative,putative$x >0.5)
breaks=seq(0,100,by=5)
hist(left$identity,main="",ylim=c(0,350),xlab="",ylab="",breaks=breaks,col=rgb(0, 0, 255, max = 255, alpha = 125) )
hist(right$identity,add=T,breaks=breaks,col=rgb(0, 255, 0, max = 255, alpha = 125))
title(xlab="Identity to 3D7 reference (%)",line=2.25)
title(ylab="ecDNA Reads",line=2.25)
legend(3, 350, legend=c("First arm", "Second arm"),
       col=c(rgb(0, 0, 255, max = 255, alpha = 125), rgb(0, 255, 0, max = 255, alpha = 125)), pch=15,  cex=0.99, box.lty=0)
