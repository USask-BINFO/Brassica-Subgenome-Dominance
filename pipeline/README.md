# Analysis pipeline

Scripts in run order. Paths at the top of each script point at the analysis working directory
(`expression_rebuild/`); adjust those if the layout differs.

| script | what it does |
|:-|:-|
| `01_build_saf.py`, `01b_build_saf_exon.py` | SAF annotations for featureCounts from each species' own annotation, keyed on transcript identifiers with a 95% match-rate guard |
| `02_build_gref.py` | the conserved homoeologous gene set and the per-subgenome count matrix |
| `02b_build_gref_Q10.py`, `02c_build_gref_span.py` | the same at MAPQ >= 10 and with gene-span counting, for the sensitivity analysis |
| `03_heb.R` | homoeologous expression bias: paired design, size factors held at 1, two one-sided tests so no pair is balanced by default |
| `03b_heb_Q10.R`, `03c_heb_span.R` | HEB under the two counting variants |
| `04_eld.R` | expression-level dominance, 12-pattern framework |
| `05_joint.R` | the joint model across tissues |
| `06_go.py` | GO annotation propagated over `is_a` and `part_of`, then enrichment |
| `07_summed.py` | summed subgenome expression |
| `08_ltr_vs_bias.py` | flanking LTR density against bias, *napus* only (superseded by 12) |
| `09_kaks.py` | K~a~/K~s~ per subgenome from SynMap K~n~ and K~s~ |
| `10_diploid_kaks.py` | the diploid-diploid control, protein-guided codon alignment and Nei-Gojobori |
| `11_peak_ks.py` | K~s~ per subgenome under both estimators; shows the published values are medians |
| `12_ltr_all_species.py` | LTR density against bias in all three allotetraploids, plus the diploids |
| `12b_ltr_coge_sensitivity.py` | the same from the complete CoGe annotations, raising the join rate |
| `13_nigra_library_workaround.py` | cross-host replication of each subgenome, and why background normalisation fails |
| `14_run_repeatmasker.sh` | re-annotates the three diploids on one common TE library; validates against the existing calls first |
| `15_nigra_compare.py` | the validation check and the resulting diploid comparison |
| `16_exact_analyzecoge.py` | faithful port of the published synteny pipeline; reproduces mean gene pairs per block exactly and identifies the undocumented settings |
| `17_4dtv.py`, `17b_4dtv_clustalw.py` | 4DTv under our aligner and under clustalw + pal2nal, on identical pairs |
| `18_inherited_vs_novel.R` | inheritance of bias: regression of allotetraploid bias on progenitor bias |
| `19_cytonuclear.py` | the maternal test, using the fact that the maternal subgenome is subgenome 1 in both testable species |
| `20_kaks_per_gene.py` | K~a~/K~s~ against bias within pairs |
| `21_hierarchy.py` | transitivity of the three contests, and the structural orders |
| `22_heb_inference.py` | exact binomial tests and Wilson intervals for every contest at every threshold |
| `23_centromere.py` | centromere distance against bias |
| `24_er_control.py` | the expression-rate control that removes the raw K~a~/K~s~ association |

## Added after the first release

| script | what it produces |
|:-|:-|
| `25_retention.py` | gene retention per Br-subgenome layer, over the 27,203 Arabidopsis anchors (Figure S12a) |
| `26_divergence_collinearity.py` | per-pair K~s~, K~n~ and percent identity from the SynMap blocks, and chromosomal collinearity (Figures S11, S13) |
| `27_pergene_bias_and_ltr_layers.py` | per-gene bias correlation between samples, and LTR density per Br-subgenome layer (Figures S14, S15) |
| `28_retention_syntenic.py` | retention measured as the share of each progenitor's own genes still syntenic, in 50-gene bins (Figure S12b) |
| `29_table_s1_alignment.py` | Table S1, streamed from the 39 BAMs with pysam, flags tallied as `samtools flagstat` does |
| `30_tables_s2_s9.py` | Tables S2, S9 and S10 |
| `32_diploid_density_one_pipeline.py` | gene-proximal LTR density of all three diploids from the single RepeatMasker pass (Table S7) |

`14_run_repeatmasker.sh` and `15_nigra_compare.py` annotate the three diploid assemblies on the
common TE library and check the run against the original calls; `32` then reads that output.

## Order

`01` -> `02` -> `03`/`04`/`05` give the expression results. `09`-`17` give the evolutionary
metrics. `18`-`24` are the analyses added for this version, and `25`-`32` those added after the first release: inheritance, the gene-level tests and
their controls. `22` should be run before quoting any dominance ratio, because it supplies the
intervals.
