#!/usr/bin/env python
"""
Process fastq files throughout a given directory with the specified enzymatic cut sites present.
Runs process_radtags from the Stacks software. Merging fastq files when appropriate.


This entire script is more or less replaceable by

    parallel 'r1={}; process_radtags -c -q --filter_illumina -r -e ecoRI -i gzfastq -o process_radtags_output_parallel -1 $r1 -2 ${r1/_1./_2.}' ::: *_1.fastq.gz

and would probably never had been born if I knew about gnu parallel at its inception.
"""

from __future__ import print_function

import argparse
import fnmatch
import os
import re
import subprocess
import sys

def main(rest_enzyme_list = None, input_dir = os.getcwd(), file_pattern='*_1.fastq*',
         output_dir='process_radtags_output', overwrite_output_dir=False):
    enzyme_choices = ('apeKI', 'bamHI', 'claI', 'dpnII', 'eaeI', 'ecoRI',
                      'ecoT22I', 'hindIII', 'mluCI', 'mseI', 'mspI', 'ndeI',
                      'nheI', 'nlaIII', 'notI', 'nsiI', 'pstI', 'sau3AI',
                      'sbfI', 'sexAI', 'sgrAI', 'sphI', 'taqI', 'xbaI')
    enzyme_choices_lower = map(str.lower, enzyme_choices)
    if not rest_enzyme_list:
        raise ValueError("No restriction enzyme specified; this is a required value.")
    if len(rest_enzyme_list) > 2:
        print("Too many enzymes specified ({0}). Truncating to first two enzymes " \
              "(\"{1}\" and \"{2}\").".format(len(rest_enzyme_list), rest_enzyme_list[0],
                                              rest_enzyme_list[1]), file=sys.stderr)
        rest_enzmye_list = rest_enzyme_list[:2]
    for rest_enzyme in rest_enzyme_list:
        if rest_enzyme.lower() not in enzyme_choices_lower:
            raise ValueError("Restriction enzyme \"{0}\" invalid. Valid values " \
                             "include:\n{1}".format(rest_enzyme, ", ".join(enzyme_choices)))
            ## TODO call fails if i write e.g. ecori instead of ecoRI. Fix
    files_to_process = []
    input_dir = os.path.abspath(input_dir)
    naming_pattern = re.compile('(?P<file_name>.+?)(?P<read_number>_1){0,1}.fastq(?P<gzip_ext>.gz\w{0,2})')

    for root, dirs, files in os.walk(input_dir):
        matching_files = fnmatch.filter(files, file_pattern)
        for filename in matching_files:
            m = naming_pattern.match(filename)
            if m:
                sample_info = { "files": [ os.path.abspath(os.path.join(root, filename)) ],
                                "gzip" : m.group('gzip_ext')}
                print("Found sequencing data file:\t{}".format(sample_info['files'][0]), file=sys.stderr)
                if m.groupdict()["read_number"]:
                    # Read is one of a pair (paired-end data)
                    ## TODO this doesn't work if there is a non-standard file naming system
                    read_2 = os.path.abspath("{0}_2.fastq{1}".format(
                                     m.group('file_name'), m.group('gzip_ext')))
                    #read_2 = files.pop(files.index(read_2_string))
                    if os.path.exists(read_2):
                        sample_info['files'].append(read_2)
                        print("     ...with matching file:\t{}".format(read_2), file=sys.stderr)
                    merge_file = "{0}.fastq".format(m.group('file_name'))
                    sample_info["merge_file"] = "{0}.fastq".format(m.group('file_name'))
                    #print("Processing matching pair: {} and {}".format(read_1, read_2))
                files_to_process.append(sample_info)
    if not files_to_process:
        raise ValueError("No files matching pattern \"{0}\" found under directory" \
                         " \"{1}\". Exiting.".format(file_pattern, input_dir))
    # Create the output directory if it does not exist
    if os.path.exists(output_dir):
        if not overwrite_output_dir:
            raise OSError("Output directory \"{0}\" exists; " \
            "will not overwrite.".format(output_dir))
    else:
        os.makedirs(output_dir)
    # Build the command
    for sample_set in files_to_process:
        cl =  ["process_radtags"]
        # -c : clean reads (remove reads with N)
        # -q : discard reads with low quality scores
        # --filter_illumina : discard reads marked as failed by CASAVA (~redundant)
        cl += ["-c", "-q", "--filter_illumina"]
        # Rescue radtags
        cl += ["-r"]
        cl += ["-e", rest_enzyme_list[0]]
        if len(rest_enzyme_list) > 1:
            cl += ["--renz_2", rest_enzyme_list[1]]
        if sample_set["gzip"]:
            cl += ["-i", "gzfastq"]
        cl += ["-o", output_dir]
        if len(sample_set["files"]) > 1:
            cl += ["-1", sample_set["files"][0], "-2", sample_set["files"][1]]
        else:
            cl += ["-f", sample_set["files"][0]]
        print("Executing ctommand line: \"{}\"".format(" ".join(cl)), file=sys.stderr)
        subprocess.check_call(cl)

        # Merge individual reads into one file per sample
        if len(sample_set["files"]) > 1:
            files_to_merge = []
            merge_output_dir = os.path.join(output_dir, "merged")
            if not os.path.exists(merge_output_dir):
                os.makedirs(merge_output_dir)
            for file_to_merge in sample_set["files"]:
                file_to_merge = os.path.basename(file_to_merge)
                if sample_set["gzip"]:
                    file_to_merge = os.path.splitext(file_to_merge)[0]
                file_to_merge = os.path.join(output_dir, file_to_merge)
                files_to_merge.append(file_to_merge)
            with open(os.path.join(merge_output_dir, sample_set["merge_file"]), 'a+') as outf:
                print("Merging files\n\t\"{0}\" and\n\t\"{1}\" into file\n\t\"{2}\"".format(
                      files_to_merge[0], files_to_merge[1], outf.name), file=sys.stderr)
                for file in files_to_merge:
                    with open(file) as inf:
                        for line in inf:
                            outf.write(line)
        _, columns = os.popen('stty size', 'r').read().split()
        columns = int(columns)
        if not columns:
            columns = 80
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
    #parser.add_argument("-m", "--file-pattern", default="*_1.fastq*",
    #                    help="The pattern to use to locate files to be processed. Should correspond to read 1. Default \"*_1.fastq*\". " +
    #                         "Note that autodetection of paired-end files is not supported for custom patterns.")
    #parser.add_argument("-p", "--paired-end", action="store_true", dest="paired_end",
    #                    help="Files are paired-end and should be merged after processing. Default false.")
    parser.add_argument("-e", "--rest-enzyme", dest="rest_enzyme_list", action="append", required=True,
                        help="""The restriction enzyme used to digest the genomic DNA; can be used up to twice.
                                Currently supported enzymes include:

                                    'apeKI', 'bamHI', 'claI', 'dpnII', 'eaeI', 'ecoRI',
                                    'ecoT22I', 'hindIII', 'mluCI', 'mseI', 'mspI', 'ndeI',
                                    'nheI', 'nlaIII', 'notI', 'nsiI', 'pstI', 'sau3AI',
                                    'sbfI','sexAI, 'sgrAI', 'sphI', 'taqI', or 'xbaI'.
                                """)
    ## TODO add --quiet option that just prints a progress bar
    args = vars(parser.parse_args())
    main(**args)
