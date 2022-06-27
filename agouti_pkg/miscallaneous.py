from operator import attrgetter
from agouti_pkg.sequence_ontology import *
from io import StringIO
import pandas as pd
from agouti_pkg.eprint import eprint


def statistics(agouti_output : str):
    """Create statistics based on the agouti output

    Args:
        agouti_output (str): output of the agouti software
    """
    
    output = StringIO(agouti_output)
    df = pd.read_csv(output, sep="\t", index_col=False) 
    
    columns = df.columns
    
    #check the index of the first column for which stats will be calculated
    for i in range(0, len(columns)):
        if columns[i].startswith("annotated_") and not columns[i+1].startswith("annotated_"):
            n = i+1
            break
    cols = list(columns[n:])
    if "annotated_featuretype" in columns:
        cols.append("annotated_featuretype")
    eprint("##"*10)
    eprint("##STATISTICS")
    eprint("##"*10)
    for col in cols:
        try:
            df[col] = df[col].str.strip()
        except AttributeError:
            pass
        output = []
        col_stats = df[col].value_counts(dropna=False).to_dict()
        for key, value in col_stats.items():
            output.append(f"'{key}' - {value}")
        eprint(f"# statistics for the '{col}' column: {'; '.join(output)}")


def sum_features_length(features):
    """Returns the sum of the features lengths. Features should be of Feature\
        type

    Arguments:
        features {list} -- list or other iterable object of Features

    Returns:
        int -- sum of features lengths
    """

    length = 0
    for f in features:
        length += (f.end - f.start)
    return length


def get_level1_parent(db, child, id=False):
    """Returns 1st level parent or parent id of a given Feature.
        If child has no parrent child will be returned

    Arguments:
        db {modules.database.Database} -- Database object
        child {Feature} -- Feature object

    Keyword Arguments:
        id {bool} -- if True returns id of the object instead of Feature\
                     object (default: {False})

    Returns:
        Feature -- child Feature or id if id parameter is set to True
    """

    for p in db.database.parents(child):
        p_parents = db.database.parents(p)
        if sum(1 for a in db.database.parents(p)) == 0:
            return p.id if id else p
        elif sum(1 for c in db.database.parents(p)) == 1:
            for a in db.database.parents(p):
                if p.featuretype == a.featuretype:
                    return p.id if id else p
    return child.id if id else child


def completly_within(start, end, feature):
    """Checks if product described by start and end coordinates is completly\
       within feature.

    Arguments:
        start {int} -- start coordinate of annotated product
        end {int} -- end coordinate of annotated product
        feature {Feature} -- Feature object

    Returns:
        bool -- True if product is completly within Feature object,\
                False otherwise
    """

    return True if start >= feature.start and end <= feature.end else False


def get_min_start(featurelist):
    """Returns the smallest feature start coordinate from the list of Feature\
        objects

    Arguments:
        featurelist {list} -- list of Feature objects

    Returns:
        int -- min start coordinate from a given set
    """

    starts = []
    for feature in featurelist:
        starts.append(feature.start)
    return min(starts)


def get_max_end(featurelist):
    """Returns the biggest feature end coordinate from the list of Feature\
        objects

    Arguments:
        featurelist {list} -- list of Feature objects

    Returns:
        int -- max end coordinate from a given set
    """

    ends = []
    for feature in featurelist:
        ends.append(feature.start)
    return max(ends)


def distinguish_UTRs(utr, cds):
    """Distinguish 3' from 5' UTRs if not annotated

    Arguments:
        utr {list} -- list of UTRs
        cds {list} -- list of CDS

    Returns:
        tuple -- tuple of sets containing 3' and 5' UTRs
    """

    cds_start, cds_end = get_min_start(cds), get_max_end(cds)
    three_primes, five_primes = set(), set()
    for u in utr:
        if (((u.end <= cds_start) and u.strand == "+")
           or ((u.start >= cds_end) and u.strand == "-")):
            five_primes.add(u)
        elif (((u.start >= cds_end) and u.strand == "+")
              or ((u.end <= cds_start) and u.strand == "-")):
            three_primes.add(u)
    return (five_primes, three_primes)


def get_exons(feature, db):
    """Returns exons of a given transcript.

    Arguments:
        feature {Feature} -- transcript
        db {modules.database.Database} -- Database object

    Returns:
        list -- list of exons
    """

    utr5_children = db.get_children(feature, featuretypes=["exon"], level=1)
    exons_children = db.get_children(feature, featuretypes=["exon"], level=1)
    return list(exons_children)


def get_cds(feature, featuretypes_from_db, db):
    """Returns CDS of a given transcript

    Arguments:
        feature {Feature} -- transcript
        featuretypes_from_db {list} -- list of feature types from the database
        db {modules.database.Database} -- Database object

    Returns:
        tuple -- tuple of CDS
    """

    synonyms = cds_synonyms
    intersection = set(synonyms) & set(featuretypes_from_db)
    children = db.get_children(feature, featuretypes=intersection, level=1)
    return children


def get_three_prime_UTRs(feature, featuretypes_from_db, db):
    """Returns 3' UTR regions of a given transcript

    Arguments:
        feature {Feature} -- transcript
        featuretypes_from_db {list} -- list of feature types from the database
        db {modules.database.Database} -- Database object

    Returns:
        list -- list of 3 prime UTRs. Usually of length equal to 1
    """

    synonyms = [three_prime_UTR_synonyms, UTR_synonyms]
    intersection = set(synonyms[0]) & set(featuretypes_from_db)
    utr3_children = db.get_children(feature, featuretypes=intersection,
                                    level=1)
    intersection = set(synonyms[1]) & set(featuretypes_from_db)
    utr_children = db.get_children(feature, featuretypes=intersection, level=1)
    cds = get_cds(feature, featuretypes_from_db, db)
    five_primes, three_primes = [], []
    if len(utr_children) and len(cds) and len(intersection):
        five_primes, three_primes = distinguish_UTRs(utr_children, cds)

    result_list = []
    try:
        for a in (list(utr3_children) + list(three_primes)):
            if "+" in feature.strand:
                if not a.start < get_max_end(cds):
                    result_list.append(a)
            elif "-" in feature.strand:
                if not a.start > get_min_start(cds):
                    result_list.append(a)
    except ValueError:
        return []
    return result_list


def get_five_prime_UTRs(feature, featuretypes_from_db, db):
    """Returns 5' UTR regions of a given transcript

    Arguments:
        feature {Feature} -- transcript
        featuretypes_from_db {list} -- list of feature types from the database
        db {modules.database.Database} -- Database object

    Returns:
        list -- list of 3 prime UTRs. Usually of length equal to 1
    """

    synonyms = [five_prime_UTR_synonyms, UTR_synonyms]
    intersection = set(synonyms[0]) & set(featuretypes_from_db)
    utr5_children = db.get_children(
        feature, featuretypes=intersection, level=1)
    intersection = set(synonyms[1]) & set(featuretypes_from_db)
    utr_children = db.get_children(feature, featuretypes=intersection, level=1)
    cds = get_cds(feature, featuretypes_from_db, db)

    five_primes, three_primes = [], []
    if len(utr_children) and len(cds) and len(intersection):
        five_primes, three_primes = distinguish_UTRs(utr_children, cds)

    result_list = []
    try:
        for a in (list(utr5_children) + list(five_primes)):
            if "+" in feature.strand:
                if not a.start > get_min_start(cds):
                    result_list.append(a)
            elif "-" in feature.strand:
                if not a.start < get_max_end(cds):
                    result_list.append(a)
    except ValueError:
        return []
    return result_list


def handle_negative_coordinates(start, end):
    """Checks if coordinates are outside transcript boundaries

    Arguments:
        start {int} -- start coordinate
        end {int} -- end coordinate

    Returns:
        tuple -- updated start coordinate, updated end coordinate,
                 string describing which coordinates are outside
                 transcript boundaries
    """

    coord_outside_transcript = "No"

    if (start < 0 and end > 0):
        start_coord, end_coord = 0, end
        coord_outside_transcript = "just_one"
    elif (start < 0 and end < 0):
        start_coord, end_coord = 0, 0
        coord_outside_transcript = "both"
    else:
        start_coord, end_coord = start, end

    return start_coord, end_coord, coord_outside_transcript


def infer_feature_level(feature, db):
    """Infers feature level

    Arguments:
        feature {Feature} -- feature which level is to be inferred
        db {Database} -- object of the Database class

    Returns:
        int -- level
    """

    level = 1
    for parent in db.database.parents(feature):
        if (parent.featuretype != feature.featuretype):
            level = level + 1
    return level


def completly_within_positive_coordinates(lengths_dict, feature,
                                          processing_product, cds_start):
    """Checks if coordinates of ProcessingProduct object don't exceed\
        transcript length.

    Arguments:
        lengths_dict {dict} -- dictionary of utr and cds lengths
        feature {Feature} -- feature, transcript
        processing_product {ProcessingProduct} -- object of the\
            ProcessingProduct class
        cds_start {int} -- start of the CDS

    Returns:
        tuple -- ...
    """

    length = sum(lengths_dict[feature.id])
    five_prime_UTR_len = lengths_dict[feature.id][0]
    if (processing_product.coordinates[2] >
            length and processing_product.coordinates[1] > length):
        return ("both")
    elif (processing_product.coordinates[2] > length):
        return ("just_one", length, five_prime_UTR_len)
    return ("No", length, five_prime_UTR_len)


def add_quotation_marks(l):
    """Add quotation marks for each string in list

    Arguments:
        l {list} -- list of strings

    Returns:
        [list] -- returns list of strings with added quotation marks
    """

    return ["'" + x + "'" for x in l]


def check_position_within_transcript(args, lengths, num, processing_product):
    """Checks whether processing_product lies within 3' UTR, CDS or 5' UTR

    Arguments:
        args {argparse.Namespace} -- argparse command-line arguments
        lengths {tuple} -- 5'UTR length, CDS length, 3'UTR length
        num {int} -- index of the feature. Might be one 0, 1 or 2
        processing_product {ProcessingProduct} -- object of ProcessingProduct\
            class

    Returns:
        bool -- True if ProcessingProduct lies within 3'UTR, 5'UTR or CDS
    """
    
    if (args.completly_within):

        if (num == 0 and processing_product.coordinates[1] >= args.offset and
                processing_product.coordinates[2] < args.offset + lengths[0]):
            # if inside 5'UTR
            return True
        # if inside CDS
        elif (num == 1 and processing_product.coordinates[1] >= args.offset +
              lengths[0] and processing_product.coordinates[2] <=
              (args.offset + lengths[0] + lengths[1])):
            return True
        # if inside 3'UTR
        elif (num == 2 and processing_product.coordinates[1] > (args.offset +
              lengths[0] + lengths[1]) and processing_product.coordinates[2] <
              args.offset + lengths[0] + lengths[1] + lengths[2]):
            return True
        else:
            return False

    else:

        if (num == 0 and processing_product.coordinates[2] > args.offset
            and processing_product.coordinates[1] < args.offset + lengths[0]):
            return True
        elif (num == 1 and (processing_product.coordinates[1] in
              range(args.offset + lengths[0] + 1, args.offset + lengths[0] +
              lengths[1])   
              or processing_product.coordinates[2] in
              range(args.offset + lengths[0] + 1, args.offset +
              lengths[0] + lengths[1])
              or processing_product.coordinates[1] in
              range(args.offset + lengths[0] + 1, args.offset +
              lengths[0] + lengths[1])
              or (processing_product.coordinates[1] <=
              args.offset + lengths[0] and processing_product.coordinates[2] >=
              args.offset + lengths[0] + lengths[1]))):
            return True
        elif (num == 2 and processing_product.coordinates[2] > args.offset +
              lengths[0] + lengths[1] and processing_product.coordinates[1] <
              args.offset + lengths[0] + lengths[1] + lengths[2]):
            return True
        else:
            return False


def find_closest_gene(database, processing_product):
    """Find the closest gene to the processing_product

    Arguments:
        database {Database} -- Database object
        processing_product {ProcessingProduct} -- ProcessingProduct object

    Returns:
        [str] -- string to display in the output file
    """

    f1 = ", ".join(add_quotation_marks(list(database.features_at_1_level)))
    if "-" not in processing_product.strand:

        upstream = database.database.execute(
            'SELECT max(end), id from features where seqid == "{}" AND \
            featuretype IN ({}) AND end <= {}'.format(
                processing_product.coordinates[0], f1,
                processing_product.coordinates[1]))

        downstream = database.database.execute(
            'SELECT min(start), id from features where seqid == "{}" AND \
            featuretype IN ({}) AND start >= {}'.format(
                processing_product.coordinates[0], f1,
                processing_product.coordinates[2]))

    else:
        downstream = database.database.execute(
            'SELECT max(end), id from features where seqid == "{}" AND \
            featuretype IN ({}) AND end <= {}'.format(
                processing_product.coordinates[0], f1,
                processing_product.coordinates[1]))

        upstream = database.database.execute(
            'SELECT min(start), id from features where seqid == "{}" AND \
            featuretype IN ({}) AND start >= {}'.format(
                processing_product.coordinates[0], f1,
                processing_product.coordinates[2]))

    values = []

    for i in upstream.fetchone():
        values.append(i)

    for i in downstream.fetchone():
        values.append(i)

    if "-" not in processing_product.strand:

        return "closest gene upstream: {} - distance {} bp; closest gene downstream: {} - distance {} bp".format(values[1], processing_product.coordinates[1] -
                         values[0] if values[0] is not None else None,
                         values[3], values[2] -
                         processing_product.coordinates[2]
                         if values[2] is not None else None)
    else:

        return "closest gene upstream: {} - distance {} bp; closest gene downstream: {} - distance {} bp".format(values[1], values[0]-processing_product.coordinates[2]
                         if values[0] is not None else None, values[3],
                         processing_product.coordinates[1] - values[2]
                         if values[2] is not None else None)


def filter_overlapping_features(attributes_and_features, args, db, feature,
                                overlapping_features, header_features,
                                header_attr, processing_product=None,
                                lengths=None):
    """Filter overlapping features so that they meet criteria specified by\
        user.
       Parse the results to be displayed in the output file.

    Arguments:
        attributes_and_features {dict} -- dictionary of attributes and\
            features that need to be annotated
        args {argparse.Namespace} -- argparse command-line arguments
        db {Database} -- Database object
        feature {Feature} -- single overlapping Feature object
        overlapping_features {list} -- list of overlapping features
        header_features {list} -- features present in the header
        header_attr {list} -- attributes present in the header

    Keyword Arguments:
        processing_product {ProcessingProduct} -- object of the
            ProcessingProduct class (default: {None})
        lengths {dict} -- dictionary of utr and cds lengths (default: {None})

    Returns:
        str -- results to be displayed
    """

    children_ids = db.get_children(feature, ids=True)
    d = {}
    position = check_position_within_transcript  # abbreviation

    for h in (header_features + header_attr):
        d[h] = "."
    

    if (lengths and feature.featuretype in mRNA_synonyms):  # if transcriptomic

        if lengths[1] != 0:

            for f in d.keys():

                if f in three_prime_UTR_synonyms and lengths[2] != 0:

                    d[f] = "y" if position(args, lengths, 2,
                                            processing_product) else "."

                elif f in five_prime_UTR_synonyms and lengths[1] != 0:

                    d[f] = "y" if position(args, lengths, 0,
                                            processing_product) else "."

                elif f in UTR_synonyms:

                    d[f] = "y" if (
                        position(args, lengths, 2, processing_product) and
                        lengths[2] != 0) or (position(args, lengths, 0,
                                                processing_product) and
                                                lengths[0] != 1) else "."

                elif f in cds_synonyms:
                    d[f] = "y" if position(args, lengths, 1,
                                            processing_product) else "."

        else:

            for f in d.keys():

                if (f in three_prime_UTR_synonyms or
                   f in five_prime_UTR_synonyms or f in UTR_synonyms or
                   f in cds_synonyms):
                    d[f] = "NA"

    else:

        for o in overlapping_features:

            if (o.id in children_ids and
               o.featuretype in attributes_and_features.keys()):

                d[o.featuretype] = "y"

    # attributes
    analyzed_set = set(attributes_and_features[feature.featuretype]
                       ).intersection(list(feature.attributes))
    for a in header_attr:

        for attr in analyzed_set:
            if attr == a:

                d[attr] = ("").join(feature[attr])
    # concatenation
    result = ""

    for h in (header_features + header_attr):
        result += "{}\t".format(d[h])

    return result.strip()
