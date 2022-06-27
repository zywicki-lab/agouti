#!/usr/bin/env python
# coding=utf-8

import os
import sys
import codecs
from agouti_pkg.eprint import eprint
import agouti_pkg.gffutils
import agouti_pkg.gffutils.inspect as inspect
from operator import attrgetter
from agouti_pkg.database import Database
from agouti_pkg.argument_parser import (parse_arguments,
                             create_attributes_and_features_dict)
from agouti_pkg.sequence_ontology import *
from agouti_pkg.processing_product import ProcessingProduct
from agouti_pkg.miscallaneous import *
from agouti_pkg.read_input import *
from agouti_pkg.header import *
from agouti_pkg.output_processing import prepare_output


def main(args):
    """Main function

    Arguments:
        args {argparse.Namespace} -- command line arguments parsed with\
            argparse
    """

    chromosomes_not_found = set()  ## chromosomes in the input file but not found in the reference annotations
    chromosomes_found = set()  ## chromosomes in the input file and found in the reference annotations

    num_of_bed_fields = -1
    lengths_dict = {}  # stores length of UTRs and CDS of a given transcript

    try:
        db = Database(agouti_pkg.gffutils.FeatureDB(args.database, keep_order=False),
                    args.database)  # reading the database from file
    except ValueError:
        eprint("ERROR: the database provided does not exists")
        sys.exit()

    # list all featuretypes existing in db
    featuretypes_from_db = list(db.database.featuretypes())

    try:

        attr_and_feat = create_attributes_and_features_dict  # abbr
        attributes_and_features = attr_and_feat(db, args, featuretypes_from_db)


    except IndexError:

        eprint("ERROR: arguments provided with --combine, --select_attributes or --select_features are incorrect")
        sys.exit()

    try:
        if (args.custom == "BED"):
            try:
                products, num_of_bed_fields = read_BED_file(args.bed,
                                                            num_of_bed_fields,
                                                            args.first_base_num, args.header_lines)
            except FileNotFoundError:
                eprint("ERROR: the input file does not exists")
                sys.exit()
        else:
            try:
                products, num_of_bed_fields = read_custom_format(args.bed, args.custom,
                                                                args.sep,
                                                                args.first_base_num,
                                                                num_of_bed_fields, args.header_lines)
            except FileNotFoundError:
                eprint("ERROR: the input file does not exists")
                sys.exit()
    except IndexError:
        eprint("ERROR: the input file has the wrong format")
        sys.exit()

    header = prepare_header(db, attributes_and_features, args,
                            num_of_bed_fields, args.header_lines)
    
    whole_output = f"{header[0].strip()}"

    for key, value in products.items():
        if (args.transcriptomic):
            try:
                id = value.coordinates[0]

                overlapping_features = [db.database[id]]
                g2t, cds_start = value.genomic_to_transcriptomic(
                    db, featuretypes_from_db)
                lengths_dict[value.coordinates[0]] = g2t

                out = prepare_output(lengths_dict, header, args,
                                     attributes_and_features, db, value,
                                     overlapping_features, g2t, cds_start)
                
                whole_output = f"{whole_output}\n{out}"

            except (agouti_pkg.gffutils.exceptions.FeatureNotFoundError):
                eprint("WARNING: Couldn't find transcript {}. Please make sure that it exists in your GTF/GFF3 file.".format(id))

        else:

            region = value.coordinates

            ## test whether chromosome is present in annotations
            if value.coordinates[0] not in chromosomes_found:
                if len(list(db.database.region(seqid=value.coordinates[0]))) == 0:
                    chromosomes_not_found.add(value.coordinates[0])
                else:
                    chromosomes_found.add(value.coordinates[0])
            ##
            

            if (not args.strand_specific):
                overlapping_features = list(
                    db.database.region(region=region,
                                       featuretype=attributes_and_features.keys(),
                                       completely_within=False))

            elif(args.strand_specific):
                overlapping_features = list(db.database.region(region=region,
                                                               strand=value.strand,
                                                               featuretype=attributes_and_features.keys(),
                                                               completely_within=False))

            out = prepare_output(lengths_dict, header, args,
                                 attributes_and_features, db, value,
                                 overlapping_features, region)
            
            whole_output = f"{whole_output}\n{out}"

    if args.statistics or args.stats_only:
        statistics(whole_output)
    
    if not args.stats_only:
        print(whole_output)

    if len(chromosomes_not_found):
        eprint(f"WARNING: The following chromosomes were not found in the annotations: {', '.join(list(chromosomes_not_found))}")
    return


if __name__ == "__main__":
    sys.exit(main(sys.argv))
