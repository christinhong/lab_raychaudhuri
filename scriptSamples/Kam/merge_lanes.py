#!/usr/bin/env python
# encoding: utf8
"""
merge_lanes.py
Kamil Slowikowski
July 17, 2015

Broad Technology labs provides sequencing data in a directory structure like
this:

    C1-123_2012-06-01_2012-06-03
    ├── 1
    │   ├── HCHT1BGXX.1.AGGCAGAA_AGAGGATA.unmapped.1.fastq.gz
    │   ├── HCHT1BGXX.1.AGGCAGAA_AGAGGATA.unmapped.2.fastq.gz
    │   ...
    ├── 2
    │   ├── HCHT1BGXX.2.AGGCAGAA_AGAGGATA.unmapped.1.fastq.gz
    │   ├── HCHT1BGXX.2.AGGCAGAA_AGAGGATA.unmapped.2.fastq.gz
    │   ...
    ...
    │
    └── info
        ├── logs
        │   ├── generateGetManifest.HCHT1BGXX.0.0.log
        │   ├── generateGetSite.HCHT1BGXX.0.0.log
        │   ├── ID_riker-12345-1234567890123.json
        │   ├── logs.tar.gz
        │   ├── sendFileEmail.sendGetSiteEmail.HCHT1BGXX.0.0.log
        │   └── workflow.log
        └── summary

This script reads the JSON file:

    info/logs/ID_riker-12345-1234567890123.json
    
and prints the Unix commands to merge FASTQ files from multiple lanes:

    gzip -cd 1/HCHT1BGXX.1.AGGCAGAA_AGAGGATA.unmapped.1.fastq.gz \
        HCHT1BGXX.2.AGGCAGAA_AGAGGATA.unmapped.1.fastq.gz \
        | pigz > sampleAlias.R1.fastq.gz

If we execute the commands, we will create new files like this:

    C1-166_2015-07-02_2015-07-06
    ├── sampleAlias.R1.fastq.gz
    └── sampleAlias.R2.fastq.gz

"""

import json
import os
import argparse

# We want to ignore this error:
#   IOError: [Errno 32] Broken pipe
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)

def main():
    args = get_args()

    x = json.loads(args.json.read())

    sampleAliases = get_sampleAliases(x)

    for sampleAlias in sampleAliases:

        for end in (1, 2):

            fastqs = filenames(x, sampleAlias, end)

            outform = '{sampleAlias}.R{end}.fastq.gz'
            
            d = {
                'sampleAlias': sampleAlias,
                'end': end
            }

            out = outform.format(**d)

            cmdform = 'gzip -cd {} | pigz > {}'

            print cmdform.format(' '.join(fastqs), out)


def get_args():
    parser = argparse.ArgumentParser(
        description="""Print commands needed to merge FASTQ files from
        multiple lanes."""
    )
    parser.add_argument(
        'json', metavar='FILE', type=argparse.FileType('r'),
        help='JSON file found in info/logs/*.json'
    )
    return parser.parse_args()


def filenames(js, sampleAlias, end = 1):
    """Returns a list of filenames corresponding to a sampleAlias.
    """
    fname = '{lane}/{flowcellBarcode}.{lane}.{barcode1}_{barcode2}.unmapped.{end}.fastq.gz'

    flowcellBarcode = js['workflow']['flowcellBarcode']

    lanes = js['workflow']['lanes']

    retval = []

    for lane_i, lane in enumerate(lanes, 1):
        libraries = lane['libraries']

        library_i = 0
        for i, x in enumerate(libraries):
            if x['sampleAlias'] == sampleAlias:
                library_i = i

        library = libraries[library_i]

        barcode1, barcode2 = library['molecularBarcodeSequences']

        d = {
            'flowcellBarcode': flowcellBarcode,
            'lane': lane_i,
            'barcode1': barcode1,
            'barcode2': barcode2,
            'end': end
        }

        filename = fname.format(**d)

        assert os.path.isfile(filename)

        retval.append(fname.format(**d))

    return retval


def get_sampleAliases(js):
    """Returns a list of sampleAliases.
    """
    retval = []
    for lane in js['workflow']['lanes']:
        for library in lane['libraries']:
            sampleAlias = library['sampleAlias']
            if not sampleAlias in retval:
                retval.append(sampleAlias)
    return retval


if __name__ == '__main__':
    main()

