#!/usr/bin/env python3

import argparse
import csv
import glob
import hashlib
import json
import logging
import os
import re
import subprocess

import requests


def parse_assembly_summary_lines(assembly_summary_lines, assembly_summary_header):
    """
    """
    assembly_summary = []

    for line in assembly_summary_lines[1:]:
        assembly_summary_record = {}
        line_split = line.strip().split('\t')
        if len(line_split) != len(assembly_summary_header):
            continue
        for i, header in enumerate(assembly_summary_header):
            # print(i, header)
            try:
                assembly_summary_record[header] = line_split[i]
            except IndexError as e:
                print(json.dumps(assembly_summary_header))
                print(json.dumps(line_split))
                exit(-1)
        assembly_summary.append(assembly_summary_record)
        
    return assembly_summary


def parse_md5checksums_lines(md5checksums_lines):
    md5checksums_by_filename = {}
    
    for line in md5checksums_lines:
        line_split = line.strip().split()
        if len(line_split) == 2:
            [md5, filename_with_leading_dotslash] = line_split
            filename = filename_with_leading_dotslash[2:]
            md5checksums_by_filename[filename] = md5

    return md5checksums_by_filename


def download_and_check(assembly, outdir):
    assembly_outdir = os.path.join(outdir, assembly['assembly_accession'])
    os.makedirs(assembly_outdir)
    ftp_path = assembly['ftp_path']
    file_suffixes = [
        'assembly_report.txt',
        'assembly_stats.txt',
        'genomic.fna.gz',
        'genomic.gff.gz',
    ]
    md5checksums_url = os.path.join(ftp_path, 'md5checksums.txt')
    r = requests.get(md5checksums_url, allow_redirects=True)
    with open(os.path.join(assembly_outdir, 'md5checksums.txt'), 'w') as f:
        f.write(r.text)

    md5checksums_lines = r.text.split('\n')
    md5checksums = parse_md5checksums_lines(md5checksums_lines)

    md5checksums_failed = []
    assembly_accession = assembly['assembly_accession']
    assembly_name = assembly['asm_name']
    chars_to_replace_with_underscores = [' ', '(', ')', '/']
    for c in chars_to_replace_with_underscores:
        assembly_name = assembly_name.replace(c, '_')
    if '__' in assembly_name:
        assembly_name = re.sub('_+', '_', assembly_name)
    for suffix in file_suffixes:
        filename = '_'.join([assembly_accession, assembly_name, suffix])
        url = os.path.join(ftp_path, filename)
        r = requests.get(url, allow_redirects=True)
        if filename.endswith('.gz'):
            with open(os.path.join(assembly_outdir, filename), 'wb') as f:
                f.write(r.content)
        else:
            with open(os.path.join(assembly_outdir, filename), 'w') as f:
                f.write(r.text)
        md5 = hashlib.md5(open(os.path.join(assembly_outdir, filename), 'rb').read()).hexdigest()
        if md5 != md5checksums[filename]:
            os.remove(os.path.join(assembly_outdir, filename))
            md5checksums_failed.append(filename)

    with open(os.path.join(assembly_outdir, 'md5checksums_failed.txt'), 'w') as f:
        for failed in md5checksums_failed:
            f.write(failed + '\n')

    gzip_files = glob.glob(os.path.abspath(os.path.join(assembly_outdir, '*.gz')))
    for gzip_file in gzip_files:
        gunzip_cmd = [
            'gunzip',
            gzip_file,
        ]
        subprocess.run(gunzip_cmd, capture_output=True, check=True)
    logging.info(json.dumps({"event_type": "download_complete", "assembly_accession": assembly_accession, "num_md5_checksums_failed": len(md5checksums_failed)}))


def main(args):
    logging.basicConfig(
        format='{"timestamp": "%(asctime)s.%(msecs)03d", "level": "%(levelname)s", "module", "%(module)s", "function_name": "%(funcName)s", "line_num", %(lineno)d, "message": %(message)s}',
        datefmt='%Y-%m-%dT%H:%M:%S',
        encoding='utf-8',
        level=logging.INFO,
    )

    assembly_summary_url = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt"

    r = requests.get(assembly_summary_url, allow_redirects=True)
    assembly_summary_str = r.text
    assembly_summary_lines = assembly_summary_str.split('\n')[1:]
    assembly_summary_header = assembly_summary_lines[0].strip().split('\t')
    assembly_summary_header[0] = assembly_summary_header[0].replace('#', '').strip()
    assembly_summary = parse_assembly_summary_lines(assembly_summary_lines[1:], assembly_summary_header)

    with open(os.path.join(args.outdir, 'assembly_summary_full.tsv'), 'w') as f:
        writer = csv.DictWriter(f, fieldnames=assembly_summary_header,
                                dialect='excel-tab', quoting=csv.QUOTE_MINIMAL)

    filter_criteria = {
        'assembly_level': lambda x: x == "Complete Genome",
        'version_status': lambda x: x == "latest",
        'refseq_category': lambda x: x == "representative genome" or x == "reference genome",
    }

    assembly_summary_filtered = []
    for assembly_summary_record in assembly_summary:
        criteria_met = []
        for k, v in filter_criteria.items():
            criteria_met.append(v(assembly_summary_record[k]))
        if all(criteria_met):
            assembly_summary_filtered.append(assembly_summary_record)


    with open(os.path.join(args.outdir, 'assembly_summary_filtered.tsv'), 'w') as f:
        writer = csv.DictWriter(f, fieldnames=assembly_summary_header, dialect='excel-tab', quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        for row in assembly_summary_filtered:
            writer.writerow(row)


    for assembly in assembly_summary_filtered:
        assembly_accession = assembly['assembly_accession']
        if not(os.path.exists(os.path.join(args.outdir, assembly_accession))):
            download_and_check(assembly, args.outdir)
    
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outdir')
    args = parser.parse_args()
    main(args)
