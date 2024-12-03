# Final-Repository

NOX4 Gene Family Phylogenetic Analysis Pipeline
This repository contains the workflow and scripts used to analyze the phylogenetic relationships within the NOX4 gene family, starting from a human NOX4 protein sequence (NP_058627).

Prerequisites

NCBI E-utilities (ncbi-acc-download)
BLAST+ suite
Standard Unix tools (awk, grep, sed)
MUSCLE alignment tool
R (with required packages for MSA plotting)
T-Coffee
AlignBuddy
LaTeX (for PDF generation)
IQ-TREE
Newick Utilities (nw_display, nw_order, nw_reroot, nw_topology)
seqkit
gotree
ImageMagick (convert)
Notung 3.0-beta
Python 2.7
RecPhyloXML tools
thirdkind
ImageMagick
RPS-BLAST
Miller (mlr)
R with required libraries for domain plotting

Workflow
PART 1: Sequence Collection and BLAST Analysis

Sequence Retrieval
Download human NOX4 protein sequence in FASTA format:
ncbi-acc-download -F fasta -m protein "NP_058627"
BLAST Search
Perform BLAST search against custom protein database:
blastp -db ../allprotein.fas -query NP_058627.fa -outfmt 0 -max_hsps 1 -out NOX.blastp.typical.out
Filtering and Analysis
Filter BLAST results for significant hits (E-value < 1e-30):
awk '{if ($6< 1e-30)print $1 }' NOX.blastp.detail.out > NOX.blastp.detail.filtered.out
Count total number of filtered hits:
wc -l NOX.blastp.detail.filtered.out
Count species distribution in filtered results:
grep -o -E "^[A-Z].[a-z]+" NOX.blastp.detail.filtered.out | sort | uniq -c

PART 2: Multiple Sequence Alignment and Analysis

Generate Multiple Sequence Alignment
muscle -align ~/lab04-s1-mvito/toy.fas -output ~/lab04-s1-mvito/toy.al.fas
Plot MSA Visualization
Rscript --vanilla ~/lab04-$MYGIT/plotMSA.R ~/lab04-$MYGIT/NOX/NOX.homologs.al.fas
Generate PDF Output
pdflatex /home/bio312-user/lab04-s1-mvito/NOX/NOX.homologs.al.fas.tex
Alignment Analysis and Processing

View alignment:
alignbuddy -al ~/lab04-$MYGIT/NOX/NOX.homologs.al.fas
Trim all gaps:
alignbuddy -trm all ~/lab04-$MYGIT/NOX/NOX.homologs.al.fas | alignbuddy -al
Remove ambiguous positions:
alignbuddy -dinv 'ambig' ~/lab04-$MYGIT/NOX/NOX.homologs.al.fas | alignbuddy -al


Alignment Quality Assessment

Generate similarity matrix:
t_coffee -other_pg seq_reformat -in ~/lab04-$MYGIT/NOX/NOX.homologs.al.fas -output sim
Calculate average percent identity:
alignbuddy -pi ~/lab04-$MYGIT/NOX/NOX.homologs.al.fas | awk '(NR>2) { for (i=2;i<=NF;i++){ sum+=$i;num++} } END{ print(100*sum/num) }'

PART 3: Phylogenetic Tree Construction and Visualization

Sequence Preparation
Remove spaces from labels and filter duplicates:
sed 's/ /_/g' NOX.homologs.al.fas | seqkit grep -v -r -p "dupelabel" > NOX.homologsf.al.fas
Phylogenetic Tree Construction
Build tree with bootstrapping (1000 replicates):
iqtree -s NOX.homologsf.al.fas -bb 1000 -nt 2
Tree Visualization and Manipulation

Display basic tree:
nw_display NOX.homologsf.al.fas.treefile
Generate unrooted tree plot:
Rscript --vanilla plotUnrooted.R NOX.homologsf.al.fas.treefile NOX.homologsf.al.fas.treefile.pdf 0.4 15
Midpoint root the tree:
gotree reroot midpoint -i NOX.homologsf.al.fas.treefile -o NOX.homologsf.al.mid.treefile
Display ordered tree:
nw_order -c n NOX.homologsf.al.mid.treefile | nw_display -


Tree Format Conversion and Visualization

Generate SVG with branch lengths:
nw_order -c n NOX.homologsf.al.mid.treefile | nw_display -w 1000 -b 'opacity:0' -s > NOX.homologsf.al.mid.treefile.svg
Convert SVG to PDF:
convert NOX.homologsf.al.mid.treefile.svg NOX.homologsf.al.mid.treefile.pdf
Generate cladogram (topology only):
nw_order -c n NOX.homologsf.al.mid.treefile | nw_topology - | nw_display -s -w 1000 > NOX.homologsf.al.midCl.treefile.svg


Outgroup Rooting
Reroot using specific outgroups:
nw_reroot NOX.homologsf.al.fas.treefile H.sapiens_HBG1_hemoglobin_subunit_gamma1 H.sapiens_HBG2_hemoglobin_subunit_gamma2 H.sapiens_HBB_hemoglobin_subunit_beta H.sapiens_HBD_hemoglobin_subunit_delta > NOX.homologsf.outgroupbeta.treefile

PART 4: Phylogenetic Reconciliation Analysis

Prepare Gene Tree
Copy midpoint-rooted tree for reconciliation:
cp NOX.homologsf.al.mid.treefile NOX.homologs.al.mid.treefile
Reconciliation with Species Tree
Perform reconciliation using Notung:
java -jar Notung-3.0_24-beta.jar -s species.tre -g NOX.homologsf.pruned.treefile --reconcile --speciestag prefix --savepng --events --outputdir ./NOX/
Convert to RecPhyloXML
Convert Notung output to RecPhyloXML format:
python2.7 NOTUNGtoRecPhyloXML.py -g NOX.homologsf.pruned.treefile.rec.ntg --include.species
Visualization
Generate reconciliation visualization:
thirdkind -Iie -D 40 -f NOX.homologsf.al.mid.treefile.rec.ntg.xml -o NOX.homologsf.al.mid.treefile.rec.svg
Convert to PDF:
convert -density 150 NOX.homologsf.al.mid.treefile.rec.svg NOX.homologsf.al.mid.treefile.rec.pdf

PART 5: Protein Domain Analysis

Prepare Sequence File
Remove asterisks from sequence file:
sed 's/*//' NOX.homologs.fas > NOX.homologs.fas
Domain Search
Perform RPS-BLAST against Pfam database:
rpsblast -query NOX.homologs.fas -db Pfam/Pfam -out NOX.rps-blast.out -outfmt "6 qseqid qlen qstart qend evalue stitle" -evalue .0000000001
Tree and Domain Visualization
Copy midpoint-rooted tree:
cp NOX.homologsf.al.mid.treefile NOX/
Generate combined tree and domain visualization:
Rscript --vanilla plotTreeAndDomains.r NOX.homologsf.al.mid.treefile NOX.rps-blast.out NOX.tree.rps.pdf
Domain Analysis

View formatted BLAST results:
mlr --inidx --ifs "\t" --opprint cat NOX.rps-blast.out | tail -n +2
Count domains per sequence:
cut -f 1 NOX.rps-blast.out | sort | uniq -c
Count domain types:
cut -f 6 NOX.rps-blast.out | sort | uniq -c
Calculate and sort domain lengths:
awk "a=$4-$3;print $1,"\t",a;" NOX.rps-blast.out | sort -k2nr
Extract E-values:
cut -f 1,5 -d $'\t' NOX.rps-blast.out

Output Files

NP_058627.fa: Query sequence in FASTA format
NOX.blastp.typical.out: Standard BLAST output
NOX.blastp.detail.out: Detailed BLAST results
NOX.blastp.detail.filtered.out: Filtered list of significant hits
toy.al.fas: MUSCLE alignment output
NOX.homologs.al.fas.tex: LaTeX file for alignment visualization
NOX.homologs.al.fas.pdf: PDF visualization of the alignment
NOX.homologsf.al.fas: Filtered alignment file
NOX.homologsf.al.fas.treefile: IQ-TREE output tree
NOX.homologsf.al.fas.treefile.pdf: Unrooted tree visualization
NOX.homologsf.al.mid.treefile: Midpoint-rooted tree
NOX.homologsf.al.mid.treefile.svg/pdf: Tree visualizations
NOX.homologsf.al.midCl.treefile.svg/pdf: Cladogram visualizations
NOX.homologsf.outgroupbeta.treefile: Outgroup-rooted tree
NOX.homologsf.pruned.treefile.rec.ntg: Notung reconciliation output
NOX.homologsf.al.mid.treefile.rec.ntg.xml: RecPhyloXML format
NOX.homologsf.al.mid.treefile.rec.svg/pdf: Reconciliation visualizations
NOX.rps-blast.out: Domain search results
NOX.tree.rps.pdf: Combined tree and domain visualization

Notes

The pipeline uses an E-value threshold of 1e-30 to filter for significant homologs
Species distribution is analyzed using the first letter of genus and species (e.g., "H.sapiens")
The analysis uses a custom protein database (allprotein.fas)
Multiple sequence alignment is performed using MUSCLE
Alignment visualization uses custom R script (plotMSA.R)
Various alignment quality metrics are calculated using AlignBuddy and T-Coffee
Phylogenetic analysis uses IQ-TREE with 1000 bootstrap replicates
Trees are visualized in multiple formats: basic display, unrooted, midpoint-rooted, and outgroup-rooted
Two parallel threads (-nt 2) are used for tree construction
Various tree manipulations are performed using Newick Utilities
Both branch length trees and cladograms are generated
Reconciliation analysis requires a species tree (species.tre)
Notung is used for gene tree/species tree reconciliation
Visualization includes gene duplication and loss events
RecPhyloXML format allows for detailed reconciliation visualization
High-resolution PDFs are generated with 150 DPI
Domain analysis uses RPS-BLAST against Pfam database
Custom output format includes sequence IDs, lengths, positions, E-values, and titles
Domain visualization combines phylogenetic tree with domain architecture
Various analyses performed on domain distribution and characteristics
