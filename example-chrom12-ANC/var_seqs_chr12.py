import csv,sys
csv.field_size_limit(sys.maxsize) #needed to read in long strings
 
    
# Part 1: get the coordinates of genes from a region and save them in 'segments'

interstart = 0 # placeholder for start position of intergenic region 
segments = {} 


with open("12var_edited.csv","r") as csvfile:
    reader = csv.reader(csvfile) 
    header = next(reader)
    for row in reader:
                
        if (interstart != 0): # save the intergenic region following the previous gene 
            name = 'post-'+intername
            segments[name] = (interstart, int(row[2])-1)  
            
        segments[row[0]] = (int(row[2]),int(row[3])) # save the current gene 
                
        if int(row[2]) >= int(row[3]):  print('ERROR') # make sure 'start' < 'stop
           
        interstart = int(row[3])+1 
        intername = row[0]

        
        

# Part 2: get sequences from reference, split into chunks, and write to file 

chrom = 0 # placeholder for reference chromosome 
seq = '' #chrom sequence
outfile = open("12var_seqs_chunks_edited.fasta","w")
targets = ['Pf3D7_12_v3'] # reference chromosome name 


with open("3D7_PlasmoDB.fasta","r") as csvfile: 
    reader = csv.reader(csvfile)
    for row in reader:
        if row[0][0] == ">":

            if (chrom in targets): 
                
                for seg in segments.keys():
                    start = segments[seg][0]-1
                    end = segments[seg][1]
                    seqtowrite = seq[start:end]
                    

                    if len(seqtowrite) < 500: ## DESIRED LENGTH OF BLOCKS. Note that the last segment will have length up to 2X-1 
                        outfile.write(">"+seg+'\n')
                        outfile.write(seq[start:end]+'\n') 

                    else:
                        n = 500 ##  LENGTH OF BLOCKS
                        chunks = [seqtowrite[i:i+n] for i in range(0, len(seqtowrite), n)]
                        if len(chunks[-1]) < 500:
                            chunks[-2] = chunks[-2] + chunks[-1]
                            chunks.pop()
                        # write each block with an informative header 
                        for c in range(len(chunks)):
                            outfile.write('>'+seg+'.'+str(c+1)+'of'+str(len(chunks))+'\n')
                            outfile.write(chunks[c]+'\n')

            # restart seq and chrom
            chrom = row[0][1:].split(' ')[0]           
            seq = ''
        else:
            seq = seq + row[0]
            
            
# In case the target chrom is the last one in the reference, write its segments:
if (chrom in targets): 

    for seg in segments.keys():
        start = segments[seg][0]-1
        end = segments[seg][1]
        seqtowrite = seq[start:end]

        if len(seqtowrite) < 500: ## DESIRED LENGTH OF BLOCKS
            outfile.write(">"+seg+'\n')
            print(seg)
            outfile.write(seq[start:end]+'\n') 

        else:
            n = 500 ##  LENGTH OF BLOCKS
            chunks = [seqtowrite[i:i+n] for i in range(0, len(seqtowrite), n)]
            if len(chunks[-1]) < 500:
                chunks[-2] = chunks[-2] + chunks[-1]
                chunks.pop()
            # write each block with an informative header 
            for c in range(len(chunks)):
                outfile.write('>'+seg+'.'+str(c+1)+'of'+str(len(chunks))+'\n')
                outfile.write(chunks[c]+'\n')

                
outfile.close()