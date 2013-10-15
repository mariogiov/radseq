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
    for root, dirs, files in os.walk(root_path):
        read_1_files = fnmatch.filter(files, pattern)
        for filename in read_1_files:
            m = re.match('(?P<file_name>.*)_1.fastq(?P<gzip_ext>.gz)', filename)
            if m:
                read_1 = filename
                read_2_string = "{0}_2.fastq{1}".format(m.group('file_name'), m.group('gzip_ext'))
                read_2 = files.pop(files.index(read_2_string))
                files_to_process.append({ "files"   : [ os.path.abspath(os.path.join(root,read_1)),
                                                            os.path.abspath(os.path.join(root,read_2))],
                                          "gzip"        :   m.group('gzip_ext')})


    #pprint.pprint(files_to_process)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
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
        subprocess.check_call(cl)

if __name__=="__main__":
    #parser = argparse.ArgumentParser(description="Run process_radtags on paired-end fastq files using the SciLifeLab naming conventions.")
    main(rest_enzyme="ecoRI")
