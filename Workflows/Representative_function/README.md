This workflow is related to the issue https://github.com/PedroMTQ/mantis/issues/37 in the Mantis repo

### Running workflow

To run this do the following:

```unifunc cluster_function -h```

Use `-h` to check all available variables.

### How does it work?
1. parse orthogroups (i.e., cluster IDs), gene IDs and respective gene annotations
2. pre-process pre_process_annotations (here you can set whether hypothetical functions are removed or not - `keep_hypothetical`)
3. compare intra-orthogroup functional annotations in a pairwise manner
4. build clusters of functional annotations per orthogroup (similarity score is above `similarity_threshold`) and calculate the intra-cluster functional similarity, which should indicates how similar (in terms of function) the cluster is
5. obtain counts for all functional annotations per orthogroup, which indicates how many times a certain function appears in the orthogroup 
6. scale the functional clusters similarity (min-max scaling), so that is within the range 0 to 1. By default a weight of 0.5 is given to the missing annotations, but this can be changed in `unannotated_weight`
7. sum the functional annotations counts per cluster to obtain the total counts per cluster. Scale the cluster counts (min-max scaling) so that it is also within the range 0 to 1.
8. average the functional cluster similarity and counts to obtain a functional cluster score. Hhere you can choose whether to or not to take into account non-annotated genes with `remove_unannotated`.
9. output all annotations cluster with a score above `representative_threshold`


### OUTPUTS
- `OUTPUT_FUNC_SIMILARITY` contains all similarity comparisons between annotations:
the format is the following:
cluster_id,gene1,gene2,annotation1,annotation2,score

- `OUTPUT_FUNC_CLUSTER` contains all similarity comparisons between annotations:
the format is the following:
cluster_id,score,annotations

Keep in mind that the same cluster may have multiple representative annotations