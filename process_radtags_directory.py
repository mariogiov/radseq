#!/usr/bin/python

#import argparse
import fnmatch
import os
#import pprint
import re
import subprocess

def main(rest_enzyme, root_path = os.getcwd(), pattern=r'*1.fastq*', output_dir='process_radtags_output', merge=False):
    files_to_process = []
    assert pattern is not None

    # Walk the directory, finding pairs of fastq files
    for root, dirs, files in os.walk(root_path):
        read_1_files = fnmatch.filter(files, pattern)
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

    # Create the output directory if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # Build the command
    for pair in files_to_process:
        cl =  ["process_radtags"]
        cl += ["-c", "-q", "--filter_illumina"]
        cl += ["-e", rest_enzyme]
        if pair["gzip"]:
            cl += ["-i", "gzfastq"]
        cl += ["-o", output_dir]
        if len(pair["files"]) == 2:
            cl += ["-1", pair["files"][0], "-2", pair["files"][1]]
        else:
            cl += ["-f", pair["files"][0]]
        if merge:
            cl += ["--merge"]
    # Execute
        subprocess.check_call(cl)

    # Merge individual reads into one file per sample
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

        with open(os.path.join(output_dir, pair["merge_file"]), 'a+') as outf:
            print("Merging files \"{0}\" and \"{1}\" into file \"{2}\"\n".format(files_to_merge[0], files_to_merge[1], outf.name))
            for file in files_to_merge:
                #print "input file is {0}".format(file)
                with open(file) as inf:
                    for line in inf:
                        outf.write(line)
        _, columns = os.popen('stty size', 'r').read().split()
        if not columns:
            columns = 80
        print("-"*columns + "\n\n")

        # TODO: implement denovo_map.pl execution on all these files

if __name__=="__main__":
    #parser = argparse.ArgumentParser(description="Run process_radtags on paired-end fastq files using the SciLifeLab naming conventions.")
    main(rest_enzyme="ecoRI")
