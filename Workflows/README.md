This workflow is related to the issue https://github.com/PedroMTQ/mantis/issues/37 in the Mantis repo

How does it work?
1. yield orthogroups, genes and respective annotations
2. pre-process pre_process_annotations (here you can set whether hypothetical functions are removed or not, by setting REMOVE_HYPOTHETICAL to True or False, respectively)
3. compare intra-orthogroup annotations in a pairwise manner
4. build cluster of annotations where similarity score is above `SIMILARITY_THRESHOLD`
5. Get similarity score of annotation cluster
6. get counts for all annotations
7. scale annotation clusters score and annotations counts with min-max scaling to obtain a 0-1 score for both variables. By default a weight of 0.5 is given to the missing annotations, but this can be changed in `NO_ANNOT_WEIGHT`
8. Average scaled annotation clusters score and scaled annotation counts (here you can choose whether to or not to take into account non-annotated genes with `USE_UNANNOTATED`, by setting to True or False, respectively)
9. output all annotations cluster with a score above `FUNCTION_THRESHOLD`

### OUTPUTS
- `OUTPUT_FUNC_SIMILARITY` contains all similarity comparisons between annotations:
the format is the following:
cluster_id,gene1,gene2,annotation1,annotation2,score

- `OUTPUT_FUNC_CLUSTER` contains all similarity comparisons between annotations:
the format is the following:
cluster_id,score,annotations

Keep in mind that the same cluster may have multiple representative annotations