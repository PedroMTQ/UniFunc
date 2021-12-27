import argparse
import os
from datetime import datetime
import sys

from UniFunc import run_unifunc,run_example


if __name__ == '__main__':
    if '--example' not in sys.argv:
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
        run_unifunc(str1=str1,str2=str2,verbose=verbose,console_output=console_output,threshold=threshold)
    else:
        run_example()
