import argparse
import os
import sys
from unifunc import source



def argv_cluster_representative_function():
    from Workflows.Representative_function.Cluster_Representative_Function import Cluster_Representative_Function
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    #just a placeholder
    parser.add_argument('workflow')

    parser.add_argument('-i', '--input_path',
                        help='[required]\tInput with a tsv of clusters formatted as: gene_id|cluster|annotation. The tsv should be sorted by cluster! Otherwise this code will not work properly.')
    parser.add_argument('-o', '--output_folder', help='[required]\tOutput folder path')

    parser.add_argument('-v', '--verbose', help='Verbose mode for UniFunc. Default is False', action='store_true')
    parser.add_argument('-st', '--similarity_threshold', help='Threshold for similarity analysis.Default is 0.8')
    parser.add_argument('-kh', '--keep_hypothetical',
                        help='Consider hypothetical annotations as not being annotations? Default is False',
                        action='store_true')
    parser.add_argument('-ra', '--remove_unannotated',
                        help='Consider unannotations (i.e., blank functions) in the scording system? Default is False',
                        action='store_true')
    parser.add_argument('-uw', '--unannotated_weight',
                        help='Weight for unannotations (i.e., blank functions) in the scording system? Default is 0.5')
    parser.add_argument('-rt', '--representative_threshold',
                        help='Score to consider a function representative? Default is 0.8')
    parser.add_argument('-owr', '--output_without_representative',
                        help='Output clusters without a representative function? Default is False',
                        action='store_true')
    args = parser.parse_args()
    input_path = args.input_path
    output_folder = args.output_folder
    verbose = args.verbose
    similarity_threshold = args.similarity_threshold
    keep_hypothetical = args.keep_hypothetical
    remove_unannotated = args.remove_unannotated
    unannotated_weight = args.unannotated_weight
    representative_threshold = args.representative_threshold
    output_without_representative = args.output_without_representative

    if not output_folder:
        output_folder = os.getcwd() + '/rep_cluster_out/'
    if not output_folder.endswith(source.SPLITTER):
        output_folder+=source.SPLITTER

    if similarity_threshold:
        similarity_threshold = float(similarity_threshold)
    else:
        similarity_threshold = 0.8

    if unannotated_weight:
        unannotated_weight = float(unannotated_weight)
    else:
        unannotated_weight = 0.5

    if representative_threshold:
        representative_threshold = float(representative_threshold)
    else:
        representative_threshold = 0.8

    if input_path:
        Cluster_Representative_Function(input_path=input_path,
                                        output_folder=output_folder,
                                        verbose=verbose,
                                        similarity_threshold=similarity_threshold,
                                        keep_hypothetical=keep_hypothetical,
                                        remove_unannotated=remove_unannotated,
                                        unannotated_weight=unannotated_weight,
                                        representative_threshold=representative_threshold,
                                        output_without_representative=output_without_representative,
                                        )

    else:
        print('Missing input path')


def main():
    if '--example' in sys.argv:
        source.run_example()
    elif 'cluster_function' in sys.argv:
        argv_cluster_representative_function()

    else:
        parser = argparse.ArgumentParser(description=' _   _         _  _____                     \n'
                                                     '| | | |       (_)|  ___|                    \n'
                                                     '| | | | _ __   _ | |_    _   _  _ __    ___ \n'
                                                     '| | | || \'_ \\ | ||  _|  | | | || \'_ \\  / __|\n'
                                                     '| |_| || | | || || |    | |_| || | | || (__ \n'
                                                     ' \\___/ |_| |_||_|\\_|     \\__,_||_| |_| \\___|, a functional annotation text similarity analysis tool.\n\n'+
                                                     'UniFunc can be run in two modes:\n' +
                                                     'The default mode returns the similarity score (float) between the provided strings, to run it use: '
                                                     '<python UniFunc "this is string1" "this is string2">\n' +
                                                     'The secondary mode requires the user to set a threshold (e.g. 0.95), and True will be returned if the string similarity is above the threshold, and False otherwise. To run it use: ' +
                                                     '<python UniFunc string1 string2 -t 0.95>\n' +
                                                     'To use verbose mode add <-v>, to redirect output to a file, add <-t file_path>'
                                         ,formatter_class=argparse.RawTextHelpFormatter)

        parser.add_argument('str1')
        parser.add_argument('str2')
        parser.add_argument('-v','--verbose',action='store_true',help='Verbose mode for UniFunc')
        parser.add_argument('-o','--output',help='File path for console output')
        parser.add_argument('-t','--threshold',help='Threshold for similarity analysis')
        args=parser.parse_args()
        str1=args.str1
        str2=args.str2
        verbose=args.verbose
        console_output = args.output
        threshold = args.threshold
        source.run_unifunc(str1=str1,str2=str2,verbose=verbose,console_output=console_output,threshold=threshold)


if __name__ == '__main__':
    main()

