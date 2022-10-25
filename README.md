## Supporting code and data for: "Antigenic diversity maintained on extrachromosomal DNA in haploid malaria parasites"

<br>

![alt text](https://github.com/emily-ebel/varSV/blob/main/chr10_telo2.png)
Four duplicated rifin genes visible on a single Nanopore read.
<br><br>
## Visualizing structural variation on large pieces of DNA

<b>1. Define a region or regions to visualize.</b> <br><br>
Region files are based on annotated genes in a reference genome. Each chromosome with at least one region of interest needs its own region file. The file can contain more than one region from the same chromosome, such as telomeric and internal <i>var</i> regions.
<br><br>For example, [12var_edited.csv] contains all the <i>var</i> genes from chr12 in the 3D7 reference assembly of <i>Plasmodium falciparum</i>. It also includes some adjacent, non-<i>var</i> genes from each region. Note that the 'start' coordinate must be less than the 'stop' coordinate, but you can still keep track of 'orientation'.


<br>
<br>


<b>2. BLAST reference sequences from these regions to your DNA of interest.</b><br><br>
[var_seqs_chr12_edited.ipynb] uses coordinates from [12var.csv] to pull sequences from the 3D7 reference genome and divide them into 500-bp blocks. Although this step is done separately for each chromosome, the resulting blocks from different chromosomes can be concatenated before the BLAST. For example, [allvar_seqs_chunks_040721.fasta] contains 500-bp blocks from all 35 <i>var</i> regions in the 3D7 genome.

Your DNA of interest could be some contigs you assembled; a set of long reads; or any other collection of sequences that are longer than the structural variation (SV) you expect. In this example, the DNA of interest is [ANC.30kb.fasta], which is a set of Nanopore reads longer than 30 kb from the Ancestor sample in our Mutation Accumulation experiment. 

```
ncbi-blast+/2.7.1/makeblastdb -in ANC.30kb.fasta -parse_seqids -dbtype nucl -out ANC.30kb 
ncbi-blast+/2.7.1/blastn -query [allvar_seqs_chunks_040721.fasta] -task blastn -db ANC.30kb -out ANC.30kb.allvar.hit -evalue 0.00001 -word_size 8 -num_threads 6 -outfmt 6 
```

Add a header to the output.
```
cat header.hit ANC.30kb.allvar.hit > ANC.30kb.allvar2.hit
```

<br>
<br>

<b>3. Assign each gene and intergenic sequence a color and level to display.</b><br><br>
See the example file, [segcol_all_070722]. 'Level' refers to the vertical height in the dot plot, with 1 at the bottom. The first gene in each region should have level 1; the next level 2; and so on. Intergenic segments (starting with 'post-') have the same level as the preceding gene. We manually selected pairs of darker and lighter colors for genes and intergenic sequence, respectively.

<br>
<br>

<b>4. Identify sequences in your DNA of interest that correspond to your defined regions.</b> 

In this example, most reads do not cover the <i>var</i> regions we are interested in. To identify the ones that do, first align the DNA of interest to the reference genome. 
```
cat ANC.30kb.fasta.gz | ./minimap2/minimap2 -x map-ont 3D7_PlasmoDB.fasta - > ANC.30kb-3D7REF.paf
```

Then, create a masterlist of region coordinates and names. Our example, [chrom_regions_3D7REF], lists all 35 <i>var</i> regions with up to 100 kb of padding in each direction. The padding helps pull out reads that partially align to repetitive <i>var</i> genes. 

Finally, use [assign_reads_regions_paf_edited.R] to assign reads to the defined regions. The example output file is [ANC-REF.Chrom-RegionAssignments.csv].

<br>
<br>

<b>5. Run the Shiny app.</b> 

The script with the Shiny code, [varSV.R], needs three files hardcoded in:<br>
<br>-The colors and levels for each gene and intergenic sequence [segcol_all_070722.csv]
<br>-The assignments of reads to regions [segcol_all_070722.csv]
<br>-The lengths of the reads [ANC.30kb.readlengths.txt], which can be easily generated with samtools: 
 ```
samtools faidx  ANC.30kb.fasta
cut -f1-2 ANC.30kb.fasta.fai > ANC.30kb.readlengths.txt
 ```
 
Once these files are correctly specified, run the app. Browse to select the .hit file from step 2, e.g. [ANC.30kb.allvar2.hit]. Once the data load, you can click from read to read; adjust the % identity slider; or select a new region from the dropdown menu. You can also save a clean .pdf version of the current image by clicking 'save.'
