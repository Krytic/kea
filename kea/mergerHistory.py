#
# Functions and classes to generate the merger tree histogry of galaxies
#
# Author: Max Briel
import pandas as pd

class node():
    """
    The node class will can hold data and will be in a tree structure.
    """
    def __init__(self):
        self.galaxyID = None
        self.data = {}
        self.children =  []
        self.parent = None

    def addChild(self, child):
        self.children.append(child)

    def addParent(self, parent):
        self.parent = parent

    def addData(self, data):
        """
        data should be a dictionary with the data you wish to add to the node
        """
        if isinstance(data, pd.Series):
            data = data.to_dict()

        for i in data:
            if i == "galaxyID":
                self.galaxyID = data[i]
            else:
                self.data[i] = data[i]

    @property
    def type(self):
        return getType(self)



def pprint_tree(node, file=None, _prefix="", _last=True):
    """
    print the tree beginning at the given node
    """
    print(_prefix, "\\- " if _last else "|- ", node.galaxyID, " (", node.data["snapnum"], ")", sep="", file=file)
    _prefix += " " if _last else "|  "
    child_count = len(node.children)
    for i, child in enumerate(node.children):
        _last = i == (child_count - 1)
        pprint_tree(child, file, _prefix, _last)


def findNode(nodeID, node):
    """
    find a node with a specific galaxyID. Returns that node
    """
    if node.galaxyID == nodeID:
        return node
    else:
        for i in node.children:
            x = findNode(nodeID, i)
            if x != None:
                return x


def getType(node):
    """
    Calculate the type of the given node based on the bulge-to-total mass ratio.

    M < 0.4: elliptical galaxy          : 0
    0.4 < M < 1.56: other galaxies      : 1
    M > 1.56: irregular & spiral        : 2
    """
    totalMass = node.coldGas + node.stellarMass + node.bulgeMass + node.hotGas
    if totalMass-node.bulgeMass < 0.4:
        return 0
    elif totalMass-node.bulgeMass > 1.56:
        return 1
    else:
        return 2

def buildHistory(cosmological_model):
    """
    Build the merger tree histories of galaxies existing today
    """
    root = None
    last = None
    parent = None
    root_list = []

    keys = cosmological_model.keys()

    for index, row in cosmological_model.iterrows():
        if (row["descendantId"] == -1):
            root = node()
            root.addData(row)
            root_list.append(root)
            last = root
        else:
            if (last.galaxyID == row["descendantId"]):
                item = node()
                item.addData(row)
                item.addParent(last)
                last.addChild(item)
                last = item
            else:
                parent = findNode(row["descendantId"], root)
                item = node()
                item.addData(row)
                item.addParent(parent)
                parent.addChild(item)
                last = item
