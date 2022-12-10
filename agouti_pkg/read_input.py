from agouti_pkg.miscallaneous import handle_negative_coordinates
from agouti_pkg.processing_product import ProcessingProduct
from agouti_pkg.eprint import eprint
import sys
import codecs
from agouti_pkg.argument_parser import parse_arguments


def read_custom_format(input_file, custom, sep, first_base_num,
                       num_of_bed_fields, header_line_num):
    """Reads and parse input file in CUSTOM format. Returns dictionary of\
        id (key) and ProcessingProduct's (value)

    Arguments:
        input_file {str} -- input file name or path
        custom {str} -- custom file format description - as in --custom arg
        sep {str} -- column separator
        first_base_num {int} -- 0 for 0-based coordinates or 1 for 1-based
        num_of_bed_fields {int} -- number of columns in the input file
        header_line_num {int} -- number of header lines from the top of the file

    Returns:
        tuple -- first element stands for dictionary based on input file.\
                Key stands for product ID, value is an object of\
                ProcessingProduct class. Second element is of int type and\
                describes number of columns in the input file
    """

    products = {}
    custom_format = custom.split(',')
    custom_format = list(filter(None, custom_format))
    try:
        custom_format = list(map(int, custom_format))
    except ValueError:
        eprint("ERROR: incorrect format provided with --custom flag. Invalid literal for int() with base 10")
        sys.exit()
    if (len(custom_format) < 4): 
        eprint("ERROR: incorrect format provided with --custom flag or wrong separator used")
        sys.exit()


    try:
        with open(input_file) as f:
            for index, line in enumerate(f): ###
                if header_line_num - index > 0:
                    if header_line_num - index > 1:
                        print(line.strip())
                    continue
                if not line.startswith("#"):
                    splitted_line = line.split(codecs.decode(
                                                       sep,
                                                       'unicode_escape'))

                    if (num_of_bed_fields == -1):
                        num_of_bed_fields = len(splitted_line)
                    line_updated = ""
                    if line.strip().endswith(sep):
                        line_updated = "{}.\n".format(line.strip())
                    else:
                        line_updated = line
                    

                    start_coord, end_coord, coord_outside_transcript =\
                        handle_negative_coordinates(int(splitted_line[
                                                     custom_format[2] - 1]),
                                                    int(splitted_line[
                                                     custom_format[3] - 1]))

                    coordinates = (splitted_line[custom_format[1] - 1],
                                   start_coord, end_coord)

                    processing_product = splitted_line[custom_format[0] - 1]
                    score = 0

                    if len(custom_format) == 4:
                        strand = "."
                    else:
                        strand = splitted_line[custom_format[4] - 1]

                    if coord_outside_transcript != "No":
                        product = ProcessingProduct(coordinates,
                                                    processing_product,
                                                    score, strand,
                                                    line_updated.strip().replace(sep,
                                                                         "\t"),
                                                    first_base_num,
                                                    coord_outside_transcript)
                    else:
                        product = ProcessingProduct(coordinates,
                                                    processing_product,
                                                    score, strand,
                                                    line_updated.strip().replace(sep,
                                                                         "\t"),
                                                    first_base_num)

                    products[processing_product] = product
        f.close()
    except IndexError:
        eprint("ERROR: incorrect format provided with --custom flag or wrong separator used")
        sys.exit()

    return products, num_of_bed_fields


def read_BED_file(bed_file, num_of_bed_fields, first_base_num, header_line_num):
    """Reads and parse input file in BED format. Returns dictionary of\
        id (key) and ProcessingProduct's (value)

    Arguments:
        bed_file {str} -- input file name or path
        num_of_bed_fields {int} -- number of columns in the input file
        first_base_num {int} -- 0 for 0-based coordinates or 1 for 1-based
        header_line_num {int} -- number of header lines from the top of the file

    Returns:
        [tuple] -- list of ProcessingProducts objects, number of columns in\
                   the input file
    """
    products = {}
    counter = 0
    with open(bed_file) as f: 
        for index, line in enumerate(f): ###
            if header_line_num - index > 0:
                continue
            if not line.startswith("#"):
                tab = line.strip().split()
                num_of_bed_fields = len(list(filter(None, tab))) if (
                    num_of_bed_fields == -1) else num_of_bed_fields
                start_coord, end_coord, coord_outside_transcript = (
                    handle_negative_coordinates(int(tab[1]), int(tab[2])))
                coordinates = (tab[0], start_coord, end_coord)
                while (len(tab) < 6):
                    # in case of less than 6 columns in BED,
                    # fill the rest (up to 6) with blank fields
                    tab.append(".")
                processing_product, score, strand = tab[3], tab[4], tab[5]
                if processing_product == ".":
                    processing_product = f"unnamed_agouti_feature_{counter}"
                    counter += 1
                if coord_outside_transcript != "No":
                    product = ProcessingProduct(coordinates,
                                                processing_product,
                                                score, strand, line.strip(),
                                                first_base_num,
                                                coord_outside_transcript)
                else:
                    product = ProcessingProduct(coordinates,
                                                processing_product,
                                                score, strand, line.strip(),
                                                0)
                products[processing_product] = product
    f.close()
    return products, num_of_bed_fields


def read_header_line(bed_file, custom, sep, header_line_num):
    """Read and parse header line from the input file

    Arguments:
        bed_file {str} -- input file name or path
        custom {str} -- custom file format description - as in --custom arg
        sep {str} -- separator
        header_line_num {int} -- number of header lines from the top of the file

    Returns:
        [tuple] -- header line, list of header fields
    """
    f = open(bed_file)
    lines = f.readlines()
    header = False
    field_list = []
    potential_header = ""
    if header_line_num != 0:
        potential_header = lines[header_line_num - 1]
        header = True
    f.close()
    if custom != "BED" and header:
        custom_format = list(
            map(int, list(filter(None, custom.split(',')))))
        sh = list(filter(None, potential_header.strip().split(codecs.decode(
                                                       sep,
                                                       'unicode_escape'))))
        for i in range(0, len(sh)):
            if (i + 1 not in custom_format):
                field_list.append(sh[i])
    elif custom == "BED":
        sh = list(filter(None, potential_header.strip().split()))
        if len(sh) > 6:
            field_list = sh[6:]
    return header, field_list