
import argparse
import re
import os
import sys
from pathlib import Path


from unifunc import source





class Cluster_Representative_Function():
    def __init__(self,
                 input_path,
                 output_folder,
                 verbose,
                 similarity_threshold,
                 keep_hypothetical,
                 remove_unannotated,
                 unannotated_weight,
                 representative_threshold,
                 output_without_representative,
                 ):
        self.unifunc = source.UniFunc()
        self.input_path = input_path
        self.output_folder = output_folder
        self.verbose = verbose
        self.similarity_threshold = similarity_threshold
        self.keep_hypothetical = keep_hypothetical
        self.remove_unannotated = remove_unannotated
        self.unannotated_weight = unannotated_weight
        self.representative_threshold = representative_threshold
        self.output_without_representative = output_without_representative

        # do not change this, it's just for processing steps
        self.ANNOTATION_SPLITTER = '#######'
        self.accession_pattern = re.compile('[\dA-Za-z_]+\.\d+')
        self.OUTPUT_FUNC_SIMILARITY = f'{output_folder}func_sim.tsv'
        self.OUTPUT_FUNC_CLUSTER = f'{output_folder}rep_func_cluster.tsv'
        Path(self.output_folder).mkdir(parents=True, exist_ok=True)

        self.get_representative_function()

    def pre_process_annotations(self,annotation_str):
        annotation_str=annotation_str.strip()
        accession_id=self.accession_pattern.match(annotation_str)
        if accession_id:
            accession_id=accession_id.group()
            annotation_str=annotation_str.replace(accession_id,'')
        if '[' in annotation_str:
            species_info=annotation_str.split('[')[1]
            species_info=f'[{species_info}'
            annotation_str=annotation_str.replace(species_info,'')
        annotation_str=annotation_str.strip()
        for bad_str in [
            'RecName:',
            ', partial',
            'Flags:',
            'Full=',
            'Precursor',
            'MULTISPECIES:',
            '; Short=FCP;',
        ]:
            annotation_str=annotation_str.replace(bad_str,'')
        if not self.keep_hypothetical:
            if 'hypothetical' in annotation_str: annotation_str=''
            if 'predicted' in annotation_str: annotation_str=''
        #if annotation_str in ['predicted protein','predicted protein, partial','hypothetical protein','hypothetical protein, partial',]:        annotation_str=''
        return annotation_str.strip()

    def read_clustered_annotations(self):
        #this yields all clusters, one by one - less memory footprint
        #FILE NEEDS TO BE SORTED BY CLUSTER_ID
        with open(self.input_path) as file:
            file.readline()
            temp=[]
            for line in file:
                line=line.strip('\n')
                line=line.split('\t')
                gene_id,cluster_id,annotation=line
                annotation=self.pre_process_annotations(annotation)
                if temp and cluster_id!=previous_cluster_id:
                    yield previous_cluster_id,temp
                    temp=[]
                temp.append([gene_id,annotation])
                previous_cluster_id=cluster_id
        if temp:
            yield previous_cluster_id,temp

    def compare_annotations(self):
        #this will yield the cluster_id, the genes (gene_id+annotations) a non-redundant score of all-vs-all annotations
        annotations_generator=self.read_clustered_annotations()
        with open(self.OUTPUT_FUNC_SIMILARITY,'w+') as outfile:
            for cluster_id,gene_annotation in annotations_generator:
                cluster_scores={}
                genes=[]
                for gene1,annotation1 in gene_annotation:
                    genes.append([gene1,annotation1])
                    for gene2,annotation2 in gene_annotation:
                        if gene1!=gene2 and annotation1 and annotation2:
                            annotation_str=self.ANNOTATION_SPLITTER.join(sorted([annotation1,annotation2]))
                            if annotation_str not in cluster_scores:
                                if annotation1!=annotation2:
                                    score = self.unifunc.get_similarity_score(annotation1, annotation2, only_return=True,verbose=self.verbose)
                                    line=[cluster_id,gene1,gene2,annotation1,annotation2,str(round(score,3)),'\n']
                                    line='\t'.join(line)
                                    outfile.write(line)
                                    cluster_scores[annotation_str]=score
                yield cluster_id,genes,cluster_scores

    def build_clusters(self,linked_annotations):
        #this will build all clusters of annotations with score above the threshold
        pool = set(map(frozenset, linked_annotations))
        groups = []
        while pool:
            groups.append(set(pool.pop()))
            while True:
                for candidate in pool:
                    if groups[-1] & candidate:
                        groups[-1] |= candidate
                        pool.remove(candidate)
                        break
                else:
                    break
        return groups


    def get_functional_clusters(self,cluster_vector):
        '''
        this will create clusters for annotations above the similarity threshold
        basically interconnecting clusters
        e.g.
        if these pairs are above the threshold:
            (a,e)
            (b,c)
            (c,d)
        then the resulting clusters will be:
            (a,e),(b,c,d)


        test={
            f'a{ANNOTATION_SPLITTER}b':0.5,
            f'b{ANNOTATION_SPLITTER}c':0.9,
            f'c{ANNOTATION_SPLITTER}d':0.5,
            f'b{ANNOTATION_SPLITTER}d':0.9,
            f'e{ANNOTATION_SPLITTER}a':0.95,
            f'f{ANNOTATION_SPLITTER}g':0.8,
        }
        get_functional_clusters(test)

        '''
        list_scores=list(cluster_vector.values())
        list_annotation_str=list(cluster_vector.keys())
        similar_annotations=[]
        for i in range(len(list_scores)):
            if list_scores[i]>=self.similarity_threshold:
                similar_annotations.append(list_annotation_str[i].split(self.ANNOTATION_SPLITTER))
        all_annotations=set()
        for annotation_str in list_annotation_str:
            annotation1,annotation2=annotation_str.split(self.ANNOTATION_SPLITTER)
            all_annotations.add(annotation1)
            all_annotations.add(annotation2)
        linked_annotations={i:set() for i in all_annotations}
        for annotation1 in all_annotations:
            linked_annotations[annotation1].add(annotation1)
        for annotation1,annotation2 in similar_annotations:
            linked_annotations[annotation1].add(annotation2)
            linked_annotations[annotation2].add(annotation1)
        linked_annotations=list(linked_annotations.values())
        clusters=self.build_clusters(linked_annotations)
        return clusters

    def get_functional_clusters_score(self,functional_clusters):
        res={}
        for cluster in functional_clusters:
            cluster_score=0
            cluster_str=self.ANNOTATION_SPLITTER.join(cluster)
            c=0
            for annotation1 in cluster:
                for annotation2 in cluster:
                    if annotation1 is not annotation2:
                        c+=1
                        score = self.unifunc.get_similarity_score(annotation1, annotation2, only_return=True,verbose=self.verbose)
                        cluster_score+=score
            if cluster_score:
                res[cluster_str]=cluster_score/c
            else:
                res[cluster_str]=1

        return res

    def clean_up_clusters(self,cluster_dict):
        single_annots=set()
        for cluster_str in cluster_dict:
            annotations=cluster_str.split(self.ANNOTATION_SPLITTER)
            if len(annotations)==1:
                single_annots.add(annotations[0])
        added=set()
        for cluster_str in cluster_dict:
            annotations=cluster_str.split(self.ANNOTATION_SPLITTER)
            if len(annotations)>1:
                added.update(annotations)
        for annot in single_annots:
            if annot in added:
                cluster_dict.pop(annot)
        for annot in cluster_dict:
            if not annot:
                cluster_dict[annot]*=self.unannotated_weight



    def get_clusters_scores(self,cluster_scores,genes):
        annotations_counts={}
        res = {}
        for gene,annotation in genes:
            if annotation not in annotations_counts:            annotations_counts[annotation]=0
            annotations_counts[annotation]+=1
        cluster_counts={}
        for annotation in annotations_counts:
            if not annotation:
                if not self.remove_unannotated:
                    cluster_counts[annotation] = annotations_counts[annotation]
                    cluster_scores[annotation] = 1
            else:
                cluster_counts[annotation] = annotations_counts[annotation]
                cluster_scores[annotation] = 1
        for cluster_str in cluster_scores:
            annotations=cluster_str.split(self.ANNOTATION_SPLITTER)
            cluster_counts[cluster_str]=0
            for annotation in annotations:
                cluster_counts[cluster_str]+=annotations_counts[annotation]
        self.clean_up_clusters(cluster_scores)
        self.clean_up_clusters(cluster_counts)
        if cluster_scores and cluster_counts:
            scaled_scores=self.scale_dict(cluster_scores)
            scaled_counts=self.scale_dict(cluster_counts)
            for cluster_str in cluster_scores:
                res[cluster_str]=(scaled_scores[cluster_str]+scaled_counts[cluster_str])/2
        if '' in res:
            res.pop('')
        return res


    def get_representative_function_cluster(self,cluster_vector,genes):
        functional_clusters=self.get_functional_clusters(cluster_vector)
        functional_clusters_score=self.get_functional_clusters_score(functional_clusters)
        clusters_score=self.get_clusters_scores(functional_clusters_score,genes)
        return clusters_score

    def scale_dict(self,input_dict):
        res={}
        min_dict=min(input_dict.values())
        max_dict=max(input_dict.values())
        for i in input_dict:
            res[i]=self.min_max_scale(input_dict[i],min_dict,max_dict)
        return res

    def min_max_scale(self,X, minX, maxX):
        if minX == maxX: return 1
        return (X - minX) / (maxX - minX)



    def get_representative_function(self):
        annotations_scores_generator=self.compare_annotations()
        with open(self.OUTPUT_FUNC_CLUSTER,'w+') as outfile:
            for cluster_id,genes,cluster_scores in annotations_scores_generator:
                temp=[]
                representative_function=self.get_representative_function_cluster(cluster_scores,genes)
                for f in representative_function:
                    if representative_function[f]>self.representative_threshold:
                        temp.append([str(round(representative_function[f],3)),f.split(self.ANNOTATION_SPLITTER)])
                for l in temp:
                    score=l[0]
                    annotations=l[1]
                    line=[cluster_id,score]
                    line.extend(annotations)
                    line.append('\n')
                    line='\t'.join(line)
                    outfile.write(line)
                if not temp:
                    if self.output_without_representative:
                        line=f'{cluster_id}\n'
                        outfile.write(line)
