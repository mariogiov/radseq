#!/usr/bin/python

import argparse
import fnmatch
import os
import re
import subprocess

def main(rest_enzyme_list = None, input_dir = os.getcwd(), file_pattern='*_1.fastq*', output_dir='process_radtags_output', merge=True, overwrite_output_dir=False):
    enzyme_choices = ('apeKI', 'bamHI', 'dpnII', 'eaeI', 'ecoRI', 'ecoT22I', 'hindIII', 'mluCI', 'mseI', 'mspI', 'ndeI', 'nlaIII',
                      'notI', 'nsiI', 'pstI', 'sau3AI', 'sbfI', 'sexAI', 'sgrAI', 'sphI', 'taqI', 'xbaI')

    if not rest_enzyme_list:
        raise ValueError("No restriction enzyme specified; this is a required value.")
    if len(rest_enzyme_list) > 2:
        print >>sys.stderr, "Too many enzymes specified ({0}). Truncating to first two enzymes (\"{0}\" and \"{1}\").".format(
                                                            len(rest_enzyme_list), rest_enzyme_list[0], rest_enzyme_list[1])
        rest_enzmye_list = rest_enzyme_list[:2]
    # TODO: automatically format the input parameter to match the reqired input parameter for stacks
    #       e.g. ecori becomes ecoRI
    for rest_enzyme in rest_enzyme_list:
        if rest_enzyme not in enzyme_choices:
            raise ValueError("Restriction enzyme \"{0}\" invalid. Valid values include:\n{1}".format(rest_enzyme, ", ".join(enzyme_choices)))

    files_to_process = []
    # Walk the directory, finding pairs of fastq files
    for root, dirs, files in os.walk(input_dir):
        read_1_files = fnmatch.filter(files, file_pattern)
        for filename in read_1_files:
            m = re.match('(?P<file_name>.*)_1.fastq(?P<gzip_ext>.gz)', filename)
            if m:
                read_1 = filename
                read_2_string = "{0}_2.fastq{1}".format(m.group('file_name'), m.group('gzip_ext'))
                read_2 = files.pop(files.index(read_2_string))
                files_to_process.append({ "files"       : [ os.path.abspath(os.path.join(root,read_1)),
                                                            os.path.abspath(os.path.join(root,read_2))],
                                          "gzip"        :   m.group('gzip_ext'),
                                          "merge_file"  : "{0}.fastq".format(m.group('file_name')),})

    if not files_to_process:
        raise ValueError("No files matching pattern \"{0}\" found under directory {1}. Exiting.".format(file_pattern))

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
        #if merge:
        #    cl += ["--merge"]
        # Execute
        subprocess.check_call(cl)

    # Merge individual reads into one file per sample
        if merge:
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
                print("Merging files\n\t\"{0}\" and\n\t\"{1}\" into file\n\t\"{2}\"".format(files_to_merge[0], files_to_merge[1], outf.name))
                for file in files_to_merge:
                    #print "input file is {0}".format(file)
                    with open(file) as inf:
                        for line in inf:
                            outf.write(line)
            _, columns = os.popen('stty size', 'r').read().split()
            columns = int(columns)
            if not columns:
                columns = 80
            # TODO: print this to the log file too
            print("-"*columns + "\n\n")

        # TODO, maybe: implement denovo_map.pl execution on all these files

if __name__=="__main__":
    enzyme_choices = ('apeKI', 'bamHI', 'dpnII', 'eaeI', 'ecoRI', 'ecoT22I', 'hindIII', 'mluCI', 'mseI', 'mspI', 'ndeI', 'nlaIII',
                      'notI', 'nsiI', 'pstI', 'sau3AI', 'sbfI', 'sexAI', 'sgrAI', 'sphI', 'taqI', 'xbaI')
    parser = argparse.ArgumentParser(description="Run process_radtags on paired-end fastq files using the SciLifeLab naming conventions.")
    parser.add_argument("-o", "--output-dir", dest="output_dir", default="process_radtags_output",
                        help="The output directory in which to store processed files and logs. Default process_radtags_output.")
    parser.add_argument("-f", "--overwrite-output-dir", action="store_true", dest="overwrite_output_dir",
                        help="Overwrite output directory if it already exists. Default is false.")
    parser.add_argument("-i", "--input-dir", dest="input_dir", default=os.getcwd(), help="The input directory to search for data. Default cwd.")
    parser.add_argument("-p", "--file-pattern", dest="file_pattern", default="*_1.fastq*",
                        help="The pattern to use to locate files to be processed. Should correspond to read 1. Default \"*_1.fastq*\".")
    parser.add_argument("-n", "--no-merge", action="store_false", dest="merge",
                        help="Whether or not to merge processed reads from paired-end fastq files into one final output file. Default true.")
    parser.add_argument("-e", "--rest-enzyme", dest="rest_enzyme_list", action="append", required=True,
                        help="""The restriction enzyme used to digest the genomic DNA; can be used up to twice.
                                Currently supported enzymes include:

                                    'apeKI', 'bamHI', 'dpnII', 'eaeI', 'ecoRI', 'ecoT22I', 
                                    'hindIII', 'mluCI', 'mseI', 'mspI', 'ndeI', 'nlaIII', 
                                    'notI', 'nsiI', 'pstI', 'sau3AI', 'sbfI', 'sexAI', 
                                    'sgrAI', 'sphI', 'taqI', or 'xbaI'.

                                Default is ecoRI.
                                """)
    args = vars(parser.parse_args())
    main(**args)
