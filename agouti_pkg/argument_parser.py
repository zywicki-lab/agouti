import argparse
from agouti_pkg.sequence_ontology import *
import pickle
from agouti_pkg.eprint import eprint
import sys


def validate_arguments(args, parser):
    """ Additional argument validation. Applicable for
        mutually exclusive arguments, etc.

    Arguments:
        args {argparse.Namespace} -- argparse command-line arguments
        parser {ArgumentParser} - ArgumentParser object
    """

    if (args.transcriptomic and args.level == 1):
        parser.error("--level 1 cannot be combined with --transcriptomic")
    elif (not args.transcriptomic and args.annotate_relative_location):
        parser.error("--annotate_relative_location should be used with --transcriptomic")
    elif ((args.combine and (args.attributes or args.features))):
        parser.error("Use --select_features, --select_attributes OR --combine.")
    elif (args.offset != 0 and not args.transcriptomic):
        parser.error("--offset should be used with --transcriptomic")
    elif (args.sep != "\t" and args.custom == "BED"):
        parser.error("--separator must be used with --custom option")

    return


def parse_arguments(argv):
    """Parse command-line arguments for agouti

    Arguments:
        argv {list} -- list of command-line arguments (sys.argv)

    Returns:
        argparse.Namespace -- parsed arguments
    """

    parser = argparse.ArgumentParser('agouti', description='AGouTI - software for annotation of genomic and transcriptomic intervals',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subprasers = parser.add_subparsers(dest='command')
    create_db = subprasers.add_parser('create_db', help='create the database for agouti')
    create_db.add_argument('-a', '--annotation', type=str, help='input file with the genomic annotation in either GTF or GFF3',
                        required=True, dest='annotation')
    create_db.add_argument('-f', '--format', type=str, help='format of the input file with the genomic annotation',
                        required=True, choices=['GTF', 'GFF3'], dest='format')
    create_db.add_argument('-d', '--db', type=str, help='name for the output database',
                        required=True, dest='database')
    create_db.add_argument('-l', '--low-ram', help='enable low-memory mode of the database creation process (warning: slow!)',
                        action="store_true", dest='save_on_disk')
    create_db.add_argument('-i', '--infer_genes', help='infer gene features. Use only with GTF files that do not have lines describing genes (warning: slow!)',
                        action="store_true", dest='infer_genes')
    create_db.add_argument('-j', '--infer_transcripts', help='infer transcript features. Use only with GTF files that do not have lines describing transcripts (warning: slow!)',
                        action="store_true", dest='infer_transcripts')

    annotate = subprasers.add_parser('annotate', help='run annotation with agouti')
    annotate.add_argument('-i', '--input', type=str,
                        help='input file in BED or another column-based format (see --custom).', required=True,
                        dest='bed')
    annotate.add_argument('-d', '--database', type=str, help='database file created with the agouti create_db run mode',
                        required=True, dest='database')
    annotate.add_argument('-m', '--custom', type=str, help='the input text file is in custom format, other than BED. It should contain columns with information about feature id (id), chromosome (chr), start (s) and end (e) coordinates, and optionally about strand. User should provide the proper column indexes (starting from 1) in order: "id,chr,s,e[,strand]". The index of the strand column is optional. Example use: --custom 1,2,4,5,6 or --custom 1,2,4,5. The field separator used in a text file can be specified using the --separator option',
                        required=False, default="BED", dest='custom')
    annotate.add_argument('-p', '--separator', type=str, help='field separator to be used with the --custom option. Default is "\\t"', required=False, default="\t",
                        dest='sep')
    annotate.add_argument('-b', '--coordinates', type=int, help='indicate the coordinate system used in the input file (BED/CUSTOM). Either 0 (0-based coordinates) or 1 (1-based coordinates). Default is 0',
                        choices=[0, 1], required=False, default=0,
                        dest='first_base_num')
    annotate.add_argument('-n', '--header_lines', type=int, help='the number of header lines in the input file. Default is 0',
                        required=False, default=0, dest='header_lines')
    annotate.add_argument('-t', '--transcriptomic', action="store_true", help='transcriptomic annotation mode. In this mode, transcript IDs from the GTF/GFF3 are expected to be placed in the first column of provided BED file instead of chromosome names. Coordinates in this mode are assumed to reflect positions within the transcript')
    annotate.add_argument('-f', '--select_features', type=str, help='comma-separated list of feature names to be reported, e.g., "mRNA,CDS". Refer to [db_name].database.structure.txt file for a list of valid features. By default, all features are reported',
                        required=False, dest='features')
    annotate.add_argument('-a', '--select_attributes', type=str, help='comma-separated list of attribute names to be reported, e.g., "ID,description". Refer to [db_name].database.structure.txt file for a list of valid attributes. By default, all attributes are reported',
                        required=False, dest='attributes')
    annotate.add_argument('-c', '--combine', type=str, help='list of specific feature-attribute combinations to be reported. The combinations should be specified in the format: feature1-attribute1:attribute2,feature2-attribute1, e.g. "mRNA-ID:description,CDS-ID',
                        required=False, dest='combine')
    annotate.add_argument('-s', '--strand_specific', action="store_true",
                        help='strand-specific search')
    annotate.add_argument('-w', '--completly_within', action="store_true",
                        help='the annotated BED interval must be located entirely within the GTF/GFF3 feature. By default, any overlap is sufficient to trigger annotation')
    annotate.add_argument('-l', '--level', type=int, help='annotate results on a specific level (1 for gene level, 2 for mRNA, tRNA level, etc.). For available levels, refer to the tree-like representation of features in [db_name].database.structure.txt file. Please note that --level 1 cannot be combined with –transcriptomic mode. Default is 2',
                        choices=[1, 2], required=False, default=2,
                        dest='level')
    annotate.add_argument('-o', '--offset', type=int, help=argparse.SUPPRESS,
                        required=False, default=0, dest='offset')
    annotate.add_argument('-r', '--annotate_relative_location', action="store_true",
                        help='annotate the relative location of the interval within the feature. Designed to work with –transcriptomic mode')
    annotate.add_argument('--statistics', action="store_true", help='calculate additional feature statistics. Those will be displayed on the stderr',
                        dest='statistics')
    annotate.add_argument('--stats_only', action="store_true", help='calculate and display only feature statistics. No annotation will be performed.', dest='stats_only')
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()

    if args.command == "create_db":
        if ((args.infer_genes or args.infer_transcripts) and args.format == "GFF3"):
            parser.error("Use option --infer_genes and/or --infer_transcripts only with GTF file format")
    elif args.command == "annotate":
        args = parser.parse_args(argv[1:])
        validate_arguments(args, parser)

    return args


def create_attributes_and_features_dict(database, args, featuretypes_from_db):
    """Creates dictionary of features (keys) and attributes (values) provided \
        by the user via the command-line. These combinations of features and \
        attributes will be then reported and annotated.

    Arguments:
        database {modules.database.Database} -- object of Database class
        args {argparse.Namespace} -- command-line arguments
        featuretypes_from_db {list} -- list of feature types from the database

    Returns:
        dict -- dictionary of attributes and features to be annotated
    """

    attributes_and_features = {}

    if args.attributes and args.features:
        for feature in list(filter(None, args.features.strip().split(","))):
            attributes_and_features[feature.lower()] = tuple(
                map(lambda x: x.lower(), list(filter(None, args.attributes.strip().split(",")))))
    elif args.combine:
        combined_str = args.combine.strip().split(",")
        for combination in combined_str:
            c = combination.strip().split("-")
            attributes_and_features[c[0].lower()] = tuple(
                map(lambda x: x.lower(), list(filter(None, c[1].strip().split(":")))))
    else:
        try:
            with open('{}.attributes_and_features.pickle'.format(args.database), 'rb') as handle:
                temp_attributes_and_features = pickle.load(handle)
        except FileNotFoundError:
            eprint("Cannot locate file 'attributes_and_features.pickle', which should be created during database creation")
            sys.exit()
        if args.features:
            features_args_list = list(filter(None, args.features.strip().split(",")))
            for feature, attributes in temp_attributes_and_features.items():
                if feature in features_args_list:
                    attributes_and_features[feature.lower()] = attributes
        elif args.attributes:
            attributes_args_list = list(filter(None, args.attributes.strip().split(",")))
            for feature, attributes in temp_attributes_and_features.items():
                attr_list = []
                for attr in attributes:
                    if attr in attributes_args_list:
                        attr_list.append(attr)
                attr_list = list(filter(None, attr_list))
                if len(attr_list):
                    attributes_and_features[feature.lower()] = attr_list
        else:
            attributes_and_features = temp_attributes_and_features

    to_remove = []

    if (args.transcriptomic):
        for f in attributes_and_features.keys():
            if f in database.features_at_3_level and not (
                    f in cds_synonyms or f in UTR_synonyms
                    or f in three_prime_UTR_synonyms or f
                    in five_prime_UTR_synonyms):
                if (args.combine or args.features):
                    eprint("In your case (--transcriptomic), features at the 3rd level must be of CDS or UTR type")
                    sys.exit()
                else:
                    to_remove.append(f)

    for f in set(to_remove):
        del attributes_and_features[f]
    # whether attributes_and_features.keys() are at annotated level
    validity = False
    for a in attributes_and_features.keys():
        if args.level == 1:
            if a in database.features_at_1_level:
                validity = True
        elif args.level == 2:
            if a in database.features_at_2_level:
                validity = True
    if not validity:
        eprint("In your case, you need to specify at least one feature at annotated level")
        sys.exit()
    return attributes_and_features
