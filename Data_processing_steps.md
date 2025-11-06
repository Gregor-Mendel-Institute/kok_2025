## Pipelines run
1. nf-core/MNase_seq
2. RNA_seq_processing
3. PRO_seq
4. nf-core/ChIP_seq
5. Long-read

## Files
| File | How it was generated | Explaination |
|------|---------|--------------|
|"annotation.gtf"| Long-read pipeline, followed by DEseq_antisense.R , followed by extract_de_anti.sh | DE antisense transcripts |
|"as_genes.bed"  | Long-read pipeline, followed by DEseq_antisense.R , followed by extract_de_anti.sh | DE antisense transcripts |


## Figures

### Main Figures:

| Figure | <div style="width:100px">Data</div>  | Processing       |	Intermediate | Analysis/Plot | Extra files for analysis | Status/Note | 
|--------|------------|------------------|---------------|---------------|-------------|-------
| 1B | MNase-seq | MNase_seq_nfcore | MNase_computematrix.sh | MNase_plotting.Rmd | plusonenuc_filtered_final_WT.bed (*) | Missing file
| 1C | MNase-seq | MNase_seq_nfcore	| MNase_computematrix.sh | MNase_plotting.Rmd | plusonenuc_filtered_final_WT.bed (*) | Missing file
| 1D | mRNA-seq | RNA_seq_processing.sh | NA | RNA_seq_remodelers.Rmd | Figure_1.xlsx |  Update to boxplot?, Missing file
| 2A | mRNA-seq | RNA_seq_processing.sh | NA | RNA_seq_CHD1.Rmd | Figure_2.xlsx | Missing file
| 2B | mRNA-seq | RNA_seq_processing.sh | NA | RNA_seq_CHD1.Rmd | Figure_2.xlsx | Missing file
| 2C | PRO-seq | PRO_seq_pipeline | NA | PRO_seq_CHD1.Rmd	| annotations.gtf | Missing file 
| 2D | PRO-seq | PRO_seq_pipeline | NA | PRO_seq_CHD1.Rmd	| annotations.gtf | Missing file 
| 3A | PRO-seq | PRO_seq_pipeline | NA | PRO_seq_CHD1.Rmd	| annotations.gtf | Missing file 
| 3B | PRO-seq | PRO_seq_pipeline | NA | PRO_seq_CHD1.Rmd	| annotations.gtf | Missing file 
| 3D | MNase-seq | MNase_seq_nfcore	|  MNase_computematrix.sh | MNase_plotting.Rmd | as_genes.bed (*) |  Missing file
| 3E | MNase-seq | MNase_seq_nfcore | MNase_computematrix.sh | PRO_seq_CHD1.Rmd | as_genes.bed (*)  |  Missing file
| 3F | MNase-seq | MNase_seq_nfcore, PRO_seq_pipeline | MNase_computematrix.sh | PRO_seq_CHD1.Rmd | as_genes.bed (*), annotations.gtf | Missing file
| 3G | NEW
| 4A | mRNA-seq |  RNA_seq_processing.sh | NA | RNA_seq_CHD1.Rmd |  Figure_2.xlsx, sexualreproductionlist.tsv, earlymeiotic.txt, middlemeiotic.txt, latemeiotic.txt | missing files
| 4B | mRNA-seq |  RNA_seq_processing.sh | NA | RNA_seq_CHD1.Rmd |  res_tib.csv, sexualreproductionlist.tsv, earlymeiotic.txt, middlemeiotic.txt, latemeiotic.txt | Missing files
| 4E | mRNA-seq | RNA_seq_processing.sh  | NA | RNA_seq_CHD1.Rmd |  earlymeiotic.txt, middlemeiotic.txt, latemeiotic.txt, nested_genes.csv, pcg_sig_up.csv | Missing files
| 4F | mRNA-seq | RNA_seq_processing.sh | NA | RNA_seq_CHD1.Rmd | res_tib.csv, earlymeiotic.txt, middlemeiotic.txt, latemeiotic.txt, nested_genes.csv | Missing files
| 5C | mRNA-seq | RNA_seq_processing.sh | NA | RNA_seq_tailswaps.Rmd | metadata.xlsx | Missing files
| 5D | mRNA-seq | RNA_seq_processing.sh | NA | RNA_seq_tailswaps.Rmd | metadata.xlsx | Missing files
| 7A | ChIP-seq | ChIP_seq_nfcore | ChIP_seq_computematrix.sh	| NA | annotation.bed (*) | Missing file
| 7B | ChIP-seq | ChIP_seq_nfcore | ChIP_seq_computematrix.sh | ChIP_seq_CHD1.Rmd | annotation.bed (*), Quartile_hrp1_scaled.tab, Quartile_hrp3_scaled.tab,Quartile_H2Bub_scaled.tab | Missing files
| 7C | ChIP-seq | ? New figure
| 7D | ChIP-seq | ? New figure
| 7E | ChIP-seq |  ChIP_seq_nfcore | ChIP_seq_computematrix.sh	| NA | annotation.bed | Missing file
| 7F | MNase-seq | MNase_seq_nfcore | MNase_computematrix.sh | MNase_plotting.Rmd | plusonenuc_filtered_final_WT.bed(*) | Updated figure?, Missing file
| 7G | mRNA-seq | RNA_seq_processing.sh | NA | RNA_seq_prf1.Rmd | metadata.xlsx | Missing file

(*): File is needed in Intermediate step 

### Supplementary Figures: 

| Figure | <div style="width:100px">Data</div>  | Processing       |	Intermediate | Analysis/Plot | Extra files for analysis | Status/Note | 
|--------|------------|------------------|---------------|---------------|-------------|-------
| S1 | ChIP-seq | ChIP_seq_nfcore | ChIP_seq_computematrix.sh | NA | annotation.bed (*) | Missing file
| S2B-S2F | RNA-seq | RNA_seq_processing.sh | NA | RNA_seq_CHD1.Rmd | Figure_2.xlsx | Missing file
| S3A | Iso-seq | Long_read | NA | supp_figure_pairs.R | NA | Done | 
| S3B | PRO-seq | PRO_seq_pipeline | NA | PRO_seq_CHD1.Rmd	| annotations.gtf | Missing file 
| S3E | MNase-seq | MNase_seq_nfcore | MNase_computematrix.sh ||| Please double check this!
| S3F | NEW
| S3G | NEW
| S4A | MNase-seq |MNase_seq_nfcore | MNase_computematrix.sh | NA|| Please double check this!
| S4B | MNase-seq | MNase_seq_nfcore | NA | GC_analysis.Rmd || Update figure?
| S4C | MNase-seq| MNase_seq_nfcore | MNase_computematrix.sh ||| Please double check this!
| S4D | Missing?
| S4E | Iso-seq | Long_read | NA | supp_figure_pairs.R | NA | Done
| S5 | NEW
| S6A | mRNA-seq | RNA_seq_processing.sh || RNA_seq_remodelers.Rmd | Figure_1.xlsx | Missing file
| S6B | mRNA-seq | RNA_seq_processing.sh || RNA_seq_remodelers.Rmd | Figure_1.xlsx, sexualreproductionlist.txt | Missing files
| S10A | ChIP-seq | ChIP_seq_nfcore | ChIP_seq_computematrix.sh |  ChIP_seq_CHD1.Rmd | annotation.bed (*) | Missing file
| S10B | ChIP-seq | ChIP_seq_nfcore | ChIP_seq_computematrix.sh |  ChIP_seq_CHD1.Rmd | annotation.bed (*) | Missing file
| S10D | NEW
| S10E | ChIP-seq | ChIP_seq_nfcore | ChIP_seq_computematrix.sh |  ChIP_seq_CHD1.Rmd | annotation.bed (*), mean_tpm.txt | Missing file
| S10F | mRNA-seq | RNA_seq_processing.sh | NA | RNA_seq_prf1.Rmd | metadata.xlsx | Missing file
| S10G | mRNA-seq | RNA_seq_processing.sh | NA | RNA_seq_prf1.Rmd | metadata.xlsx | Missing file
| S11 | Iso-seq, mRNA-seq| Long_read | sense_process.R  |  supp_figure_sense.R | NA | Done


