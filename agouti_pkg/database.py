from itertools import tee
import sys
from agouti_pkg.eprint import eprint


class Database(object):
    def __init__(self, database, db_name):
        self.database = database
        self.parent_child_relation = []
        self.features_at_1_level = []
        self.features_at_2_level = []
        self.features_at_3_level = []
        self.name = db_name
        self.parent_child_relation = set()
        self.create_parent_child_relation()
        self.find_featuretypes_at_given_level(3)
        self.find_featuretypes_at_given_level(2)
        self.find_featuretypes_at_given_level(1)

    def create_parent_child_relation(self):
        """Creates Parent -- Child like relations of Feature types
        """

        parent_child = set()
        try:
            file = open("{}.relations".format(self.name))
            lines = file.readlines()
            for line in lines:
                tab = line.strip().split("\t")
                parent_child.add((tab[0], tab[1]))
            file.close()
        except (IOError, IndexError):
            eprint("ERROR: Couldn't find the file {}.relations -> create the database again!".format(self.name))
            sys.exit()
        self.parent_child_relation = parent_child

    def find_featuretypes_at_given_level(self, level):
        """Returns Feature types at specified level

        Arguments:
            level {int} -- level of the feature (1, 2 or 3)
        """

        pc = self.parent_child_relation
        unzipped_pc = list(zip(*pc))
        featuretypes = set()
        l1, l2, l3 = set(), set(), set()
        for p in unzipped_pc[0]:  # finding features at level 1 -> no parents
            if (p not in unzipped_pc[1]):
                l1.add(p)
        for p, c in pc:
            if p in l1:
                l2.add(c)
        for p, c in pc:
            if p in l2:
                l3.add(c)
        if (level == 1):
            self.features_at_1_level = l1
        elif (level == 2):
            self.features_at_2_level = l2
        else:
            self.features_at_3_level = l3
            

    def get_children(self, parent, ids=False, featuretypes=(), level=None):
        """Get children or children ids of the given 'parent'

        Arguments:
            parent {Feature} -- parent Feature

        Keyword Arguments:
            ids {bool} -- True if ids instead of Features are to be returned\
                          (default: {False})
            featuretypes {tuple} -- children of only these featuretypes\
                                    will be returned (default: {()})
            level {int} -- return children at this level (default: {None})

        Returns:
            [list] -- list of chilldren Features or ids
        """

        children = []
        if len(featuretypes):
            child_iter = self.database.children(parent.id,
                                                featuretype=featuretypes,
                                                level=level)
        else:
            child_iter = self.database.children(parent.id, level=level)
        for child in child_iter:
            if(child.featuretype != parent.featuretype):
                if (not ids):
                    children.append(child)
                else:
                    children.append(child.id)
        return children
