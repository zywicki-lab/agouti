from agouti_pkg.miscallaneous import *


def prepare_output(lengths_dict, header, args, attributes_and_features,
                   database, processing_product, overlapping_features,
                   region, cds_start=None):
    """Prepares a single line of the output file. Parse ProcessingProduct\
        objects.

    Arguments:
        lengths_dict {dict} -- 5'UTR length, CDS length, 3'UTR length
        header {str} -- header
        args {argparse.Namespace} -- argparse command-line arguments
        attributes_and_features {dict} -- dictionary of attributes and\
            features that need to be annotated
        database {Database} -- Database
        processing_product {ProcessingProduct} -- object of ProcessingProduct
            class
        overlapping_features {list} -- list of overlapping features
        region {tuple} -- coordinates

    Keyword Arguments:
        cds_start {int} -- start of cds coordinate (default: {None})

    Returns:
        str -- single output line
    """

    out = ""

    if len(overlapping_features):  # if there are overlapping features

        for feature in overlapping_features:
            temp_out = "."
            temp_result = ""
            temp_featuretype="."
            temp_seqid="."
            temp_overlapping_feature_start="."
            temp_overlapping_feature_end="."

            if ((feature.strand != processing_product.strand) and
               args.strand_specific and not args.transcriptomic) or (args.transcriptomic and processing_product.strand == "-"):
                continue

            if (feature.featuretype in attributes_and_features.keys()):

                if not args.transcriptomic:

                    if ((args.completly_within and not completly_within(
                            region[1], region[2], feature)) or not
                            infer_feature_level(feature, database) ==
                            args.level):
                        continue

                    if args.level == 2:
                        gene_id = "{}\t".format((get_level1_parent(database,
                                                                   feature,
                                                                   id=True)))
                        out += "{bed}\t{featureid}\t{gene_id}\t{featuretype}\t\
                            {overlapping_feature_start}\t\
                            {overlapping_feature_end}\
                            ".format(bed=processing_product.bed_line.strip(),
                                     featureid=feature.id.strip(),
                                     gene_id=gene_id.strip(),
                                     featuretype=feature.featuretype.strip(),
                                     overlapping_feature_start=feature.start,
                                     overlapping_feature_end=feature.end)
                    else:
                        gene_id = ""
                        out += "{bed}\t{featureid}\t{featuretype}\t\
                                {overlapping_feature_start}\t\
                                {overlapping_feature_end}\
                                ".format(bed=processing_product.bed_line.strip(),
                                        featureid=feature.id.strip(),
                                        featuretype=feature.featuretype.strip(),
                                        overlapping_feature_start=feature.start,
                                        overlapping_feature_end=feature.end)

                    
                    out += "\t{}\n".format(
                        filter_overlapping_features(attributes_and_features,
                                                    args, database, feature,
                                                    overlapping_features,
                                                    header[1], header[2]))


                else:
                    comp_pos = completly_within_positive_coordinates  # abbr
                    temp = comp_pos(lengths_dict, feature, processing_product,
                                    cds_start)

                    if (not infer_feature_level(feature, database) ==
                        args.level or (args.completly_within and
                        processing_product.coords_outside_transcript ==
                        "just_one") or processing_product.coords_outside_transcript ==
                        "both") or (args.completly_within and temp[0] ==
                        "just_one") or temp == "both":
                        if args.level == 2:
                            temp_out = get_level1_parent(database, feature, id=True)
                            temp_featuretype=feature.featuretype
                            temp_seqid=feature.seqid
                            temp_overlapping_feature_start=feature.start
                            temp_overlapping_feature_end=feature.end
                            analyzed_set = set(attributes_and_features[feature.featuretype]).intersection(feature.attributes)
                            d = {}
                            for a in header[2]:
                                for attr in analyzed_set:
                                    if attr == a:
                                        d[attr] = ("").join(feature[attr])
                                # concatenation
                                temp_result = ""

                            for h in (header[1] + header[2]):
                                try:
                                    temp_result += "{}\t".format(d[h].strip())
                                except KeyError:
                                    temp_result += "{}\t".format(".")
                        continue

                    if args.level == 2:
                        gene_id = get_level1_parent(database, feature, id=True)
                    else:
                        gene_id = ""

                    out += "{bed}\t{gene_id}\t{featuretype}\t{seqid}\t\
                            {overlapping_feature_start}\t\
                            {overlapping_feature_end}\
                            ".format(bed=processing_product.bed_line.strip(),
                                     featuretype=feature.featuretype.strip(),
                                     gene_id=gene_id, seqid=feature.seqid.strip(),
                                     overlapping_feature_start=feature.start,
                                     overlapping_feature_end=feature.end)

                    filt_ov = filter_overlapping_features  # abbr
                    out += "\t{}\n".format(filt_ov(attributes_and_features,
                                                   args, database, feature, [],
                                                   header[1], header[2],
                                                   processing_product, region))

                    insert = processing_product.check_overlapping_feature__position(lengths_dict,
                                                                                    feature,
                                                                                    args.transcriptomic,
                                                                                    args.offset)

                    if (args.annotate_relative_location):
                        out = out.rstrip()
                        out += "\t{}\n".format(insert)

    if (out == ""):

        additional_fields = ""

        header_lengths = (len(header[1]) + len(header[2])) if not (
            args.annotate_relative_location) else (len(header[1]) +
                                                len(header[2]) + 1)

        if (args.transcriptomic):

            feature = database.database[processing_product.coordinates[0]]
            out += "{bed}\t{gene_id}\t{featuretype}\t{seqid}\t\
                    {overlapping_feature_start}\t{overlapping_feature_end}\
                    ".format(bed=processing_product.bed_line.strip(),
                            featuretype=temp_featuretype, gene_id=temp_out, seqid=temp_seqid,
                            overlapping_feature_start=temp_overlapping_feature_start,
                            overlapping_feature_end=temp_overlapping_feature_end)

            
            if args.level != 2 or ((feature.strand != processing_product.strand) and
               args.strand_specific and not args.transcriptomic) or (args.transcriptomic and processing_product.strand == "-"):
                out += "\t{}".format(".\t.\t.")
                for _ in range(0, len(header[2])):
                    out += "\t."

                out += "\n"
            else:
                out += "\t{}\n".format(temp_result)

            if (args.annotate_relative_location):
                out = out.rstrip()
                if not temp_featuretype == ".":
                    out += "\t{}\n".format(
                    processing_product.check_overlapping_feature__position(lengths_dict, feature,
                                                                           args.transcriptomic, args.offset))

        else:

            for _ in range(1, header_lengths):
                additional_fields += "\t."

            closest_gene = find_closest_gene(database, processing_product)
            find_closest_gene(database, processing_product)

            out += "{bed}\t{featureid}\t{gene_id}\t{featuretype}\t\
                    {overlapping_feature_start}\t{overlapping_feature_end}\
                    {add}".format(bed=processing_product.bed_line.strip(),
                                  featureid=closest_gene, featuretype=".",
                                  gene_id="intergenic",
                                  overlapping_feature_start=".",
                                  overlapping_feature_end=".",
                                  add=additional_fields)

    return out.strip()


