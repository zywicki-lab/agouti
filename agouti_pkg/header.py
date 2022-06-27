from agouti_pkg.read_input import read_header_line


def prepare_bed_header(bed_file, custom, transcriptomic, num_of_bed_fields, sep, header_line_num):
    """Prepares part of the header for the features from BED line

    Arguments:
        bed_file {str} -- name or path of the input file
        custom {str} -- custom file format description - as in --custom arg
        transcriptomic {bool} -- args.transcriptomic
        num_of_bed_fields {int} -- number of columns in the input file
                sep {str} -- separator

    Returns:
        [str] -- bed header of the output file
    """

    rhl = read_header_line(bed_file, custom, sep, header_line_num)
    bed_header = ""

    if custom != "BED":
        idx = 1
        custom_format = custom.strip().split(",")
        custom_format = list(map(int, custom_format))

        if not transcriptomic:
            bed_fields = ["feature_id", "chr", "feature_start", "feature_end",
                          "strand"]
        else:
            bed_fields = ["feature_id", "transcript_id", "feature_start",
                          "feature_end", "strand"]

        indices = [index for index, value in sorted(enumerate(custom_format),
                   key=lambda x: x[1])]
        temp = 0

        for i in range(1, num_of_bed_fields + 1):

            if i not in custom_format:

                if not rhl[0]:  # if there is no header
                    bed_header = "{}\tCUSTOM_field_{}".format(bed_header, idx)
                else:
                    bed_header = "{}\t{}".format(bed_header, rhl[1][idx - 1])
                idx = idx + 1

            else:
                bed_header = "{}\t{}".format(bed_header,
                                             bed_fields[(indices[temp])])
                temp = temp + 1
    else:
        bed_fields = []

        if not transcriptomic:
            bed_fields = ["chr", "feature_start", "feature_end",
                          "feature_name", "feature_score", "strand"]
        else:
            bed_fields = ["transcript_id", "feature_start", "feature_end",
                          "feature_name", "feature_score", "strand"]

        for field in range(0, num_of_bed_fields):

            if field < 6:
                bed_header = "{}\t{}".format(bed_header, bed_fields[field])
            else:
                idx = field - 6

                if not rhl[0]:  # if there is no header
                    bed_header = "{}\tBED_field_{}".format(bed_header, idx + 1)
                else:
                    bed_header = "{}\t{}".format(bed_header, rhl[1][idx])

    return bed_header


def prepare_header(db, attributes_and_features, args, num_of_bed_fields, header_line_num):
    """Prepare header for the output file

    Arguments:
        db {Database} -- object of the Database class
        attributes_and_features {dict} -- dictionary of attributes and\
            features that need to be annotated
        args {argparse} -- parsed command line argument using argparse
        num_of_bed_fields {int} -- number of columns in the input file

    Returns:
        tuple -- header of the output file, attributes choosen by the user to\
            be annotated and present in the database,\
            attributes choosen by the user
    """

    features_to_annotate = set(db.features_at_3_level)
    features_choosen_by_user = set(attributes_and_features.keys())
    intersection = features_to_annotate & features_choosen_by_user

    if args.level == 1:
        features_at_given_level = db.features_at_1_level
    else:
        features_at_given_level = db.features_at_2_level

    attributes_choosen_by_user = set()

    for f in attributes_and_features.keys():
        if (f in features_at_given_level):
            for a in attributes_and_features[f]:
                attributes_choosen_by_user.add(a)

    attributes_choosen_by_user = sorted(attributes_choosen_by_user)
    intersection = list(intersection)
    intersection.sort()
    bed_header = prepare_bed_header(args.bed, args.custom, args.transcriptomic,
                                    num_of_bed_fields, args.sep, header_line_num)

    if (args.level == 1):
        header = "{}\tannotated_gene_id\tannotated_featuretype\
                  \tannotated_gene_start\tannotated_gene_end\
                  \t{}\t{}".format(bed_header, "\t".join(intersection),
                                   "\t".join(attributes_choosen_by_user))
    elif (args.level == 2):

        if args.transcriptomic:
            header = "{}\tannotated_gene_id\tannotated_featuretype\
                      \tannotated_chromosome\tannotated_transcript_start\
                      \tannotated_transcript_end\t{}\t{}\
                      ".format(bed_header, "\t".join(intersection),
                               "\t".join(attributes_choosen_by_user))
        else:
            header = "{}\tannotated_transcript_id\tannotated_gene_id\
                      \tannotated_featuretype\tannotated_transcript_start\
                      \tannotated_transcript_end\t{}\t{}\
                      ".format(bed_header, "\t".join(intersection),
                               "\t".join(attributes_choosen_by_user))

    if (args.annotate_relative_location):
        header = header.strip() + "\tfeature_region"

    return header, list(intersection), list(attributes_choosen_by_user)
