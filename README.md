# UniFunc

UniFunc is a text mining tool that processes and analysis text similarity between a pair of protein function annotations.
It is mainly used as a cross-linking mechanism or redundancy elimination tool when processing annotations without any sort of database identifiers.

## Using UniFunc
UniFunc can be run in two modes:


The default mode returns the similarity score (float) between the provided strings, to run it use:
`python UniFunc "this is string1" "this is string2"`

The secondary mode requires the user to set a threshold (e.g. 0.95) with the argument `-t`, and True will be returned if the string similarity is above the threshold, and False otherwise. To run it use:
`python UniFunc string1 string2 -t 0.95`

To use verbose mode add the argument `-v`, to redirect output to a file, add the argument `-t file_path`

To run a sample execution use: `python UniFunc --example`


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
