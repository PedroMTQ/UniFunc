# UniFunc

UniFunc is a text mining tool that processes and analysis text similarity between a pair of protein function annotations.
It is mainly used as a cross-linking mechanism or redundancy elimination tool when processing annotations without any sort of database identifiers.


##### To update corpus:
1. Delete all files in `UniFunc/Resources/`
2. Go to https://www.uniprot.org/uniprot/?query=reviewed 
3. Search for all protein entries
4. Choose the columns `Entry`,`Protein names`,and `Function [CC]`
5. Apply columns
6. Download results in tab separated format
7. Check if download file has these 3 headers: `Entry	Protein names	Function [CC]`
8. Rename the downloaded file to `uniprot.tab` and move it `UniFunc/Resources/`
9. Go to http://geneontology.org/docs/download-ontology/
10. Download `go.obo`
11. Move the file`go.obo` to `UniFunc/Resources/`


Here's an overview of the UniFunc workflow:


![overview](Images/overview.png)
