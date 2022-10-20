# Supporting code and data for: "Antigenic diversity maintained on extrachromosomal DNA in haploid malaria parasites"

## Visualizing structural variation on large pieces of DNA

![alt text](https://github.com/emily-ebel/varSV/blob/main/ANC_PB_chr12.png)

<b>1. Define a region or regions to visualize.</b> <br><br>
Region files are created by you, based on annotated genes in a reference genome. Each chromosome with a region of interest needs its own region file. 
<br><br>For example, [12var.csv] contains all the <i>var</i> genes from chr12 in the 3D7 reference assembly of <i>Plasmodium falciparum</i>. We also included some adjacent, non-<i>var</i> genes in this region. 
<br><br>Note that the 'start' coordinate must be less than the 'stop' coordinate, but 'orientation' can be + or -. 


<br>
<br>
<br>


<b>2. BLAST reference sequences from these regions to your DNA of interest.</b><br><br>
[var_seqs_chr12.py] uses coordinates from [12var.csv] to pull sequences from the 3D7 reference genome and divide them into 500-bp blocks. Although this step is done separately for each chromosome, blocks from different chromosomes can be concatenated before the BLAST. For example, [allvar_seqs_chunks_040721.fasta] contains 500-bp blocks from all 37 <i>var</i> regions in the 3D7 genome.

Your DNA of interest could be some contigs you assembled; a set of long reads; or any other collection of sequences that are longer than the structural variation (SV) you expect. In this example, the DNA of interest is [ANC.30kb.fasta], which is a set of Nanopore reads longer than 30 kb from the Ancestor sample in our MA experiment. 

```
ncbi-blast+/2.7.1/makeblastdb -in ANC.30kb.fasta -parse_seqids -dbtype nucl -out ANC.30kb 
ncbi-blast+/2.7.1/blastn -query [allvar_seqs_chunks_040721.fasta] -task blastn -db ANC.30kb -out ANC.30kb.allvar.hit -evalue 0.00001 -word_size 8 -num_threads 6 -outfmt 6 
```

<br>
<br>
<br>

<b>3. Identify sequences in your DNA of interest that align to your defined regions.</b> 

To do this, first align your DNA of interest to your reference genome. 
 ```
cat ANC.30kb.fasta.gz | ./minimap2/minimap2 -x map-ont 3D7_PlasmoDB.fasta - > ANC.30kb-3D7REF.paf
```

Then, create a masterlist of region coordinates and names. Our example, [chrom_regions_3D7REF], lists all 37 <i>var</i> regions.

Finally, use [assign_reads_regions_paf_edited.R] to assign reads to your defined regions. The output file is [ANC-REF.Chrom-RegionAssignments.csv].
