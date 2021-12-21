This workflow is related to the issue https://github.com/PedroMTQ/mantis/issues/37 in the Mantis repo

### Running workflow

Just run this outside the UniFunc folder, e.g.:

```python UniFunc/Workflows/Representative_function/Cluster_Representative_Function.py -h```

Use `-h` to check all available variables.

### How does it work?
1. yield orthogroups, genes and respective annotations
2. pre-process pre_process_annotations (here you can set whether hypothetical functions are removed or not - `keep_hypothetical`)
3. compare intra-orthogroup annotations in a pairwise manner
4. build cluster of annotations where similarity score is above `similarity_threshold`
5. Get similarity score of annotation cluster
6. get counts for all annotations
7. scale annotation clusters score and annotations counts with min-max scaling to obtain a 0-1 score for both variables. By default a weight of 0.5 is given to the missing annotations, but this can be changed in `unannotated_weight`
8. Average scaled annotation clusters score and scaled annotation counts (here you can choose whether to or not to take into account non-annotated genes with `remove_unannotated`)
9. output all annotations cluster with a score above `representative_threshold`

### OUTPUTS
- `OUTPUT_FUNC_SIMILARITY` contains all similarity comparisons between annotations:
the format is the following:
cluster_id,gene1,gene2,annotation1,annotation2,score

- `OUTPUT_FUNC_CLUSTER` contains all similarity comparisons between annotations:
the format is the following:
cluster_id,score,annotations

Keep in mind that the same cluster may have multiple representative annotations