#!/usr/bin/env python
"""
Process fastq files throughout a given directory with the specified enzymatic cut sites present.
Runs process_radtags from the Stacks software. Merging fastq files when appropriate.
"""

from __future__ import print_function

import argparse
import fnmatch
import os
import re
import subprocess
import sys

def main(rest_enzyme_list = None, input_dir = os.getcwd(), file_pattern='*_1.fastq*', output_dir='process_radtags_output', paired_end=False, overwrite_output_dir=False):

    enzyme_choices = ('apeKI', 'bamHI', 'dpnII', 'eaeI', 'ecoRI', 'ecoT22I', 'hindIII', 'mluCI', 'mseI', 'mspI', 'ndeI', 'nlaIII',
                      'notI', 'nsiI', 'pstI', 'sau3AI', 'sbfI', 'sexAI', 'sgrAI', 'sphI', 'taqI', 'xbaI')
    enzyme_choices_lower = map(str.lower, enzyme_choices)

    if not rest_enzyme_list:
        raise ValueError("No restriction enzyme specified; this is a required value.")
    if len(rest_enzyme_list) > 2:
        print("Too many enzymes specified ({0}). Truncating to first two enzymes (\"{1}\" and \"{2}\").".format(
                                                            len(rest_enzyme_list), rest_enzyme_list[0], rest_enzyme_list[1]), file=sys.stderr)
        rest_enzmye_list = rest_enzyme_list[:2]
    for rest_enzyme in rest_enzyme_list:
        if rest_enzyme.lower() not in enzyme_choices_lower:
            raise ValueError("Restriction enzyme \"{0}\" invalid. Valid values include:\n{1}".format(rest_enzyme, ", ".join(enzyme_choices)))

    files_to_process = []
    input_dir = os.path.abspath(input_dir)
    ## TODO this is a hack but fix this later so you don't have to write the same function twice in a row
    if paired_end:
        # Walk the directory, finding pairs of fastq files
        ## TODO add non-recursive method that doesn't walk the directory
        for root, dirs, files in os.walk(input_dir):
            read_1_files = fnmatch.filter(files, file_pattern)
            for filename in read_1_files:
                ## TODO Misses ".gzip" gzipped files
                m = re.match('(?P<file_name>.*)_1.fastq(?P<gzip_ext>.gz)', filename)
                if m:
                    read_1 = filename
                    read_2_string = "{0}_2.fastq{1}".format(m.group('file_name'), m.group('gzip_ext'))
                    read_2 = files.pop(files.index(read_2_string))
                    files_to_process.append({ "files"       : [ os.path.abspath(os.path.join(root,read_1)),
                                                                os.path.abspath(os.path.join(root,read_2))],
                                              "gzip"        :   m.group('gzip_ext'),
                                              "merge_file"  : "{0}.fastq".format(m.group('file_name')),})
    else:
        for root, dirs, files in os.walk(input_dir):
            matching_files = fnmatch.filter(files, file_pattern)
            for filename in matching_files:
                m = re.match('(?P<file_name>.*).fastq(?P<gzip_ext>.gz)', filename)
                files_to_process.append({   "files"         : [ os.path.abspath(os.path.join(root,filename))],
                                            "gzip"          :   m.group('gzip_ext'),
                                            "merge_file"    :   None})
    if not files_to_process:
        raise ValueError("No files matching pattern \"{0}\" found under directory \"{1}\". Exiting.".format(file_pattern, input_dir))

    # Create the output directory if it does not exist
    if os.path.exists(output_dir):
        if not overwrite_output_dir:
            raise OSError("Output directory \"{0}\" exists; will not overwrite.".format(output_dir))
    else:
        os.makedirs(output_dir)

    # Build the command
    for pair in files_to_process:
        cl =  ["process_radtags"]
        cl += ["-c", "-q", "--filter_illumina"]
        cl += ["-e", rest_enzyme_list[0]]
        if len(rest_enzyme_list) > 1:
            cl += ["--renz_2", rest_enzyme_list[1]]
        if pair["gzip"]:
            cl += ["-i", "gzfastq"]
        cl += ["-o", output_dir]
        if len(pair["files"]) == 2:
            cl += ["-1", pair["files"][0], "-2", pair["files"][1]]
        else:
            cl += ["-f", pair["files"][0]]
        subprocess.check_call(cl)

    # Merge individual reads into one file per sample
        if paired_end and len(pair["files"] > 1):
            files_to_merge = []
            merge_output_dir = os.path.join(output_dir, "merged")
            if not os.path.exists(merge_output_dir):
                os.makedirs(merge_output_dir)
            for file_to_merge in pair["files"]:
                file_to_merge = os.path.basename(file_to_merge)
                if pair["gzip"]:
                    file_to_merge = os.path.splitext(file_to_merge)[0]
                file_to_merge = os.path.join(output_dir, file_to_merge)
                files_to_merge.append(file_to_merge)
            with open(os.path.join(merge_output_dir, pair["merge_file"]), 'a+') as outf:
                # TODO: print this to the log file
                print("Merging files\n\t\"{0}\" and\n\t\"{1}\" into file\n\t\"{2}\"".format(files_to_merge[0], files_to_merge[1], outf.name), file=sys.stderr)
                for file in files_to_merge:
                    #print "input file is {0}".format(file)
                    with open(file) as inf:
                        for line in inf:
                            outf.write(line)
        _, columns = os.popen('stty size', 'r').read().split()
        columns = int(columns)
        if not columns:
            columns = 80
        ## TODO: print this to the log file too
        print("-"*columns + "\n\n", file=sys.stderr)

        ## TODO, maybe: implement denovo_map.pl execution on all these files

if __name__=="__main__":
    enzyme_choices = ('apeKI', 'bamHI', 'dpnII', 'eaeI', 'ecoRI', 'ecoT22I', 'hindIII', 'mluCI', 'mseI', 'mspI', 'ndeI', 'nlaIII',
                      'notI', 'nsiI', 'pstI', 'sau3AI', 'sbfI', 'sexAI', 'sgrAI', 'sphI', 'taqI', 'xbaI')
    parser = argparse.ArgumentParser(description="Run process_radtags on paired-end fastq files using the SciLifeLab naming conventions.")
    parser.add_argument("-o", "--output-dir", default="process_radtags_output",
                        help="The output directory in which to store processed files and logs. Default process_radtags_output.")
    parser.add_argument("-f", "--overwrite-output-dir", action="store_true",
                        help="Overwrite output directory if it already exists. Default is false.")
    parser.add_argument("-i", "--input-dir", default=os.getcwd(), help="The input directory to search for data. Default cwd.")
    parser.add_argument("-m", "--file-pattern", default="*_1.fastq*",
                        help="The pattern to use to locate files to be processed. Should correspond to read 1. Default \"*_1.fastq*\".")
    parser.add_argument("-p", "--paired-end", action="store_true", dest="paired_end",
                        help="Files are paired-end and should be merged after processing. Default false.")
    parser.add_argument("-e", "--rest-enzyme", dest="rest_enzyme_list", action="append", required=True,
                        help="""The restriction enzyme used to digest the genomic DNA; can be used up to twice.
                                Currently supported enzymes include:

                                    'apeKI', 'bamHI', 'dpnII', 'eaeI', 'ecoRI', 'ecoT22I',
                                    'hindIII', 'mluCI', 'mseI', 'mspI', 'ndeI', 'nlaIII',
                                    'notI', 'nsiI', 'pstI', 'sau3AI', 'sbfI', 'sexAI',
                                    'sgrAI', 'sphI', 'taqI', or 'xbaI'.

                                Default is ecoRI.
                                """)
    ## TODO add --quiet option that just prints a progress bar
    args = vars(parser.parse_args())
    main(**args)
