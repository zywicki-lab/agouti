from agouti_pkg.miscallaneous import *


class ProcessingProduct(object):
    """Feature from the BED file. One line in BED corresponds\
        to one ProcessingProduct"""

    def __init__(self, coordinates, processing_product, score, strand,
                 bed_line, first_base_num, coords_outside_transcript="No"):
        self.coordinates = coordinates
        if (first_base_num == 0):
            one_based_coordinates = (coordinates[0], coordinates[1]+1,
                                     coordinates[2])
            self.coordinates = one_based_coordinates

        self.processing_product = processing_product
        self.score = score
        self.strand = strand
        self.bed_line = bed_line
        self.coords_outside_transcript = coords_outside_transcript

    def convert_coordinates(self):
        """Converts coordinates into UCSC-like format

        Returns:
            str -- coordinates in UCSC-like format
        """

        coord_formated = "{chr}:{start}-{end}".format(
            chr=self.coordinates[0],
            start=self.coordinates[1],
            end=self.coordinates[2])

        return coord_formated

    def genomic_to_transcriptomic(self, database, featuretypes_from_db):
        """Converts genomic coordinates to transcriptomic coordinates of given\
            mRNA. Require two args. One of them is a sqlite3 database created \
            by gffutils, second is a generator object that yields featuretypes\
            from database (as strings)

        Arguments:
            database {modules.database.Database} -- Database object
            featuretypes_from_db {list} -- feature types in the database

        Returns:
            tuple -- first element is tuple of (5' utr len, cds len ,\
                     3' utr len), second is start coordinate of CDS within\
                     transcript
        """

        transcript_id = self.coordinates[0]

        mRNA = database.database[transcript_id]

        three_prime_UTRs = get_three_prime_UTRs(mRNA,
                                                list(featuretypes_from_db),
                                                database)

        five_prime_UTRs = get_five_prime_UTRs(mRNA, list(featuretypes_from_db),
                                              database)

        exons = get_exons(mRNA, database)
        cds = get_cds(mRNA, list(featuretypes_from_db), database)
        five_prime_utr_len = sum_features_length(five_prime_UTRs)
        three_prime_utr_len = sum_features_length(three_prime_UTRs)
        transcript_length = sum_features_length(exons)
        cds_len = transcript_length - five_prime_utr_len - three_prime_utr_len

        if len(cds) == 0:
            cds_start = None
        else:
            cds_start = min(cds, key=attrgetter('start')).start

        if (sum([cds_len, three_prime_utr_len, five_prime_utr_len]) !=
           transcript_length and len(cds)):

            left_side_of_cds_len = 0
            for exon in exons:
                if (cds_start >= exon.end):
                    left_side_of_cds_len += exon.end - exon.start
                elif (cds_start < exon.end and cds_start > exon.start):
                    left_side_of_cds_len += cds_start - exon.start

            five_prime_utr_len = left_side_of_cds_len
            three_prime_utr_len = (transcript_length - five_prime_utr_len
                                   - cds_len)

        lengths_tuple = (five_prime_utr_len, cds_len,
                         three_prime_utr_len)

        if lengths_tuple == (0, 0, 0):
            lengths_tuple = (0, transcript_length, 0)

        return lengths_tuple, cds_start

    def check_overlapping_feature__position(self, lengths_dict, feature,
                                            transcriptomic, offset=0):
        """Checks in which part of overlapping feature lies ProcessingProduct\
            object

        Arguments:
            feature {Feature} -- overlapping feature
            transcriptomic {bool} -- transcriptomic option (see help)

        Keyword Arguments:
            offset {int} -- offset option (see help) (default: {0})

        Returns:
            str -- localization within overlapping feature
        """
        
        start, end = min(self.coordinates[1], self.coordinates[2]), max(self.coordinates[1], self.coordinates[2]) 

        if (transcriptomic):
            feature_length = sum(lengths_dict[feature.id])
        else:
            feature_length = max(feature.stop, feature.start) - min(feature.stop, feature.start)

        pp_length = self.coordinates[2] - self.coordinates[1]

        if (offset):
            feature_start = offset
        else:
            feature_start = 0


        if start > feature_length and end > feature_length:
            loc = "downstream"
        elif self.coords_outside_transcript == "both":
            loc = "upstream"
        elif (start < feature_start + (0.25 * feature_length)
                and end < feature_start + (0.5 * feature_length)):
            loc = "5 prime"
        elif (start > feature_start + (0.25 * feature_length) and
              end < feature_start + (0.75 * feature_length)):
            loc = "middle"
        elif (start > feature_start + (0.5 * feature_length) and
              end > feature_start + (0.75 * feature_length)):
            loc = "3 prime"
        elif (start < feature_start + (0.25 * feature_length) and
              end > feature_start + (0.75 * feature_length) and
              pp_length < 0.9 * feature_length):
            loc = "whole"
        elif (start < feature_start + (0.25 * feature_length) and
              end > feature_start + (0.75 * feature_length) and
              pp_length >= 0.9 * feature_length):
            loc = "full"
        else:
            loc = "other"
        return loc
