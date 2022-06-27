# Sequence Ontology Synonyms
mRNA_synonyms = tuple(map(lambda x: x.lower(),
                          ["mRNA", "SO:0000234", "INSDC_feature:mRNA",
                           "messenger RNA", "protein_coding_transcript",
                           "transcript"]))
cds_synonyms = tuple(map(lambda x: x.lower(),
                         ["CDS", "SO:0000316", "coding sequence",
                          "coding_sequence", "INSDC_feature:CDS"]))
UTR_synonyms = tuple(
    map(lambda x: x.lower(), ["UTR", "SO:0000203", "untranslated region"]))
three_prime_UTR_synonyms = tuple(map(lambda x: x.lower(),
                                     ["three_prime_UTR", "three prime UTR",
                                      "three prime untranslated region",
                                      "3'UTR", "INSDC_feature:3'UTR",
                                      "SO:0000205"]))
five_prime_UTR_synonyms = tuple(map(lambda x: x.lower(),
                                    ["five_prime_UTR", "5'UTR",
                                     "five_prime_untranslated_region",
                                     "INSDC_feature:5'UTR", "SO:0000204"]))
