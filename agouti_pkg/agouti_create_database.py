#!/usr/bin/env python
# coding=utf-8

from agouti_pkg.eprint import eprint
import sys
import os
import pickle
import sqlite3
import agouti_pkg.gffutils
import agouti_pkg.gffutils.inspect as inspect
import argparse
from agouti_pkg.anytree import Node, RenderTree, importer
from agouti_pkg.anytree.exporter import DotExporter
from agouti_pkg.argument_parser import parse_arguments
import gzip

output_lines = []


def transform_func(x):
    x.featuretype = x.featuretype.lower()
    return x


def inspect_db(args, db):
    try:
        features = {}
        if args.annotation.endswith(".gz"):
            f = gzip.open(args.annotation, "rb")
        else:
            f = open(args.annotation)
        lines = f.readlines()
        for line in lines:
            if isinstance(line, bytes):
                line = line.decode()
            splitted = line.strip().split("\t")
            if (not line.startswith("#")):
                feature = splitted[2].lower()
                attributes = splitted[8].strip().replace('"', "").replace("'", "")
                if feature not in features.keys():
                    features[feature] = set()
                tab = attributes.split(db.dialect["field separator"])
                tab = list(filter(None, tab))
                for attribute in tab:
                    a = attribute.strip().split(db.dialect["keyval separator"])
                    a = list(filter(None, a))
                    features[feature].add(a[0])  # was a[0].lower()
        f.close()
    except IndexError:
        eprint("{}ERROR: The file {} has the wrong format".format(os.linesep,
                                                          args.annotation))
        os.system("rm -f {}".format(args.database))
        sys.exit()
    return features


def find_children(_parent, parent_child, prefix=""):
    pre = prefix+"    "
    output_lines.append("{}{}".format(prefix, _parent))
    for (parent, child) in parent_child:
        if _parent == parent:
            find_children(child, parent_child, pre)


def find_roots(x):
    """takes parent_child list of tuples as an argument. The list should be in\
        the format [(parent, child), (parent, child)], etc."""
    roots = set()
    parent_child_relation = list(zip(*x))
    for parent in parent_child_relation[0]:
        if (parent not in parent_child_relation[1]):
            roots.add(parent)
    return roots


def add_super_feature(parent_child):
    """Function adds super-feature if there are multiple first level features"""
    unzipped = list(zip(*parent_child))
    modified_parent_child = set()
    first_level_features = []
    for parent in unzipped[0]:
        if parent not in unzipped[1]:
            first_level_features.append(parent)
    if len(first_level_features) > 1:  # if more than one root
        for i in first_level_features:
            modified_parent_child.add((".", i))
    return parent_child.union(modified_parent_child)


def main(args):
    """Main function

    Arguments:
        argv {argparse.Namespace} -- command-line arguments parsed with\
            argparse
    """

    global output_lines

    try:
        if args.format == 'GFF3':
            if (not args.save_on_disk):
                db = agouti_pkg.gffutils.create_db(args.annotation, force=True,
                                        keep_order=False,
                                        sort_attribute_values=False,
                                        merge_strategy="create_unique",
                                        transform=transform_func,
                                        dbfn=":memory:",
                                        checklines=50)
            else:
                db = agouti_pkg.gffutils.create_db(args.annotation, force=True,
                                        keep_order=False,
                                        sort_attribute_values=False,
                                        merge_strategy="create_unique",
                                        transform=transform_func,
                                        dbfn=args.database,
                                        checklines=50)
        else:
            if (not args.save_on_disk):
                db = agouti_pkg.gffutils.create_db(args.annotation, force=True,
                                    keep_order=False,
                                    sort_attribute_values=False,
                                    merge_strategy="create_unique",
                                    disable_infer_genes= not args.infer_genes,
                                    disable_infer_transcripts=not args.infer_transcripts,
                                    transform=transform_func,
                                    dbfn=":memory:", checklines=50)
            else:
                db = agouti_pkg.gffutils.create_db(args.annotation, force=True,
                                    keep_order=False,
                                    sort_attribute_values=False,
                                    merge_strategy="create_unique",
                                    disable_infer_genes= not args.infer_genes,
                                    disable_infer_transcripts=not args.infer_transcripts,
                                    transform=transform_func,
                                    dbfn=args.database, checklines=50)

    except ValueError:
        eprint("{}ERROR: The file {} has the wrong format or does not exist".format(os.linesep, args.annotation))
        sys.exit()

    database_struct = open("{}.database.structure.txt".format(args.database), "w")

    database_struct.write("{}Attributes available for each feature type:{}\n".format(os.linesep, os.linesep))
    attributes_and_features = inspect_db(args, db)
    with open('{}.attributes_and_features.pickle'.format(args.database), 'wb') as handle:
        pickle.dump(attributes_and_features, handle,
                    protocol=pickle.HIGHEST_PROTOCOL)  # saving the dictionary
    for feature, attributes in attributes_and_features.items():
        database_struct.write("{}: {}\n".format(feature, ", ".join(set(attributes))))
    
    if args.format == "GTF":
        if (("gene" not in attributes_and_features.keys() and not args.infer_genes) or ("transcript" not in attributes_and_features.keys() and not args.infer_transcripts)):
            eprint("{}ERROR: The file {} seems to be missing gene or transcript features. Please use --infer_genes and/or --infer_transcripts".format(os.linesep, args.annotation))
            sys.exit()
    
    parent_child = set()
    # Find featuretypes without parents
    if args.format == "GFF3":
        for featuretype in db.featuretypes():
            children = set()
            for f in db.iter_by_parent_childs(featuretype, level=1):
                children.update(set([x.featuretype for x in f if x.featuretype != featuretype]))
            if len(children):
                for c in children:
                    parent_child.add((featuretype,c))
    
    if args.format == "GTF":
        parent_child.add(("gene", "transcript"))
        for i in db.featuretypes():
            if i not in ["gene", "transcript"]:
                parent_child.add(("transcript", i))
        
    
    file = open("{}.relations".format(args.database), "w")
    for relation in parent_child:
        file.write("{}\t{}\n".format(relation[0], relation[1]))
    file.close()
    output_lines = []

    database_struct.write(os.linesep)
    database_struct.write(os.linesep)

    parent_child = add_super_feature(parent_child)
   
    imp = importer.IndentedStringImporter()
    try:
        for r in find_roots(parent_child):
            output_lines = []
            find_children(r, parent_child)
            root = imp.import_(output_lines)
            tree = RenderTree(root.children[0])
            for pre, fill, node in tree:
                database_struct.write("%s%s\n" % (pre, node.name))
    except IndexError:
        eprint("{}ERROR: File {} has the wrong format. Couldn't create parent-child relations. If you are using GTF file format, please check if gene and transcript feature lines are present - if not, use --infer_genes option or consider using another annotation file. Check if transcript and gene names differ and are correctly formatted.".format(os.linesep, args.annotation))
        os.system("rm -f {}".format(args.database))
        sys.exit()

    if (not args.save_on_disk):
        bck = sqlite3.connect(args.database)
        with bck:
            db.conn.backup(bck)
        bck.close()
    print("-"*10)
    print("The pipeline finished successfully!")
    print("Available attributes and relations are available in the file {}.database.structure.txt".format(args.database))
    database_struct.close()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
