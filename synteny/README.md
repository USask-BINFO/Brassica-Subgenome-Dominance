# Collinear blocks between each Allo-subgenome and its diploid progenitor

WGDI collinearity output, one file per comparison, gene order rather than base-pair
coordinates. These are the blocks behind the counts in the Methods and behind the
block-number comparison in the Results.

| file | Allo-subgenome | progenitor | blocks | median anchors | IQR | smallest |
|:-|:-|:-|-|-|-|-|
| `An_vs_progenitor.collinearity.gz` | A~n~ | *rapa* v3.0 | 1,467 | 8 | 6-20 | 5 |
| `Aj_vs_progenitor.collinearity.gz` | A~j~ | *rapa* v3.0 | 1,513 | 8 | 6-21 | 5 |
| `Bj_vs_progenitor.collinearity.gz` | B~j~ | *nigra* NI100 v2 | 1,944 | 8 | 6-22 | 5 |
| `Bc_vs_progenitor.collinearity.gz` | B~c~ | *nigra* NI100 v2 | 1,899 | 8 | 6-20 | 5 |
| `Cn_vs_progenitor.collinearity.gz` | C~n~ | *oleracea* v2.1 | 1,704 | 9 | 6-22 | 5 |
| `Cc_vs_progenitor.collinearity.gz` | C~c~ | *oleracea* v2.1 | 1,784 | 9 | 6-22 | 5 |

Each block begins with a header line carrying its score, p value, anchor count `N=`
and the two chromosomes, followed by one line per anchor pair. Gene identifiers are
the renamed forms used in that analysis, not the published accessions.

Block counts:

```bash
for f in *.collinearity.gz; do echo -n "$f "; zgrep -c '^# Alignment' "$f"; done
```

## Parameters

`collinearity.conf` is the WGDI configuration used, with the per-comparison file names
replaced by placeholders. WGDI v0.74, run on gene order. Anchors came from a BLASTN
search of coding sequences at an E value of 1e-5, with at most 20 repeats kept per gene.
The remaining settings are WGDI defaults.
