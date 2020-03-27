#
# Functions and classes to generate the merger tree history of galaxies
# Also contains functions to extract data from the tree structure.
#
# Author: Max Briel
import pandas as pd
import numpy as np

class node():
    """A node class to build a tree structure.

    Attributes
    ----------
    galaxyID : int
        The identifier of the node.
    data : dictionary
        A dictionary containing data
    children : array
        An array of child nodes
    parent : node
        The parent node

    """
    def __init__(self):
        self.galaxyID = None
        self.data = {}
        self.children =  []
        self.parent = None

    def addChild(self, child):
        """Add a child to the node

        Parameters
        ----------
        child : node
            Child to add to the node
        """

        self.children.append(child)
        return

    def addParent(self, parent):
        """ Add a parent to the node.

        Parameters
        ----------
        parent : node
            The parent of the current node.
        """
        self.parent = parent

    def addData(self, data):
        """ Function to add data to the node

        Parameters
        ----------
        data : dict or pandas Series
            A dictionary or a pandas Series of data to add to the node

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


def getSFR(galaxy):
    """ Extract the stellar formation rate given a root node of a merger tree.

    Parameters
    ----------
    galaxy : node
        A node of a galaxy

    Returns
    -------
    numpy array
        A numpy array of the stellar formation rate over time.

    """
    sfr = np.zeros(64)

    sfr[int(galaxy.data["snapnum"])] += galaxy.data["sfr"]

    for i in galaxy.children:
        sfr += getSFR(i)
    return sfr


def pprint_tree(node, file=None, _prefix="", _last=True):
    """ Print the node and its children; basically the tree.

    Parameters
    ----------
    node : node
        The starting node from where to print
    file : string
        output to a file
    _prefix : string
        A prefix for printing
    _last : Boolean
        Check if this is the last node

    """
    print(_prefix, "\\- " if _last else "|- ", node.galaxyID, " (", node.data["snapnum"], ")", sep="", file=file)
    _prefix += " " if _last else "|  "
    child_count = len(node.children)
    for i, child in enumerate(node.children):
        _last = i == (child_count - 1)
        pprint_tree(child, file, _prefix, _last)


def findNode(nodeID, node):
    """Find a node in a tree.

    Parameters
    ----------
    nodeID : int
        The node to find in the tree
    node : node
        A node with possible children attached.

    Returns
    -------
    node
        A node associated with the given **nodeID**

    """
    if node.galaxyID == nodeID:
        return node

    else:
        for i in node.children:
            x = findNode(nodeID, i)
            if x != None:
                return x


def getType(node):
    """ calculate the galaxy type according to ... paper
    with the bulge-to-total mass ratio.

    Parameters
    ----------
    node : node
        a node containing the required data:
            - coldGas
            - stellarMass
            - bulgeMass
            - hotGas

    Returns
    -------
    int
        An integer identifingthe type of galaxy.

    """
    totalMass = node.coldGas + node.stellarMass + node.bulgeMass + node.hotGas
    if totalMass-node.bulgeMass < 0.4:
        return 0
    elif totalMass-node.bulgeMass > 1.56:
        return 1
    else:
        return 2


def buildHistory(cosmological_model):
    """Build the merger tree histories of galaxies existing today.

    Parameters
    ----------
    cosmological_model : pandas DataFrame
        A pandas DataFrame containing the following information:
        - galaxyID
        - stellarMass
        - bulgeMass
        - sfr
        - hotGas
        - coldGas
        - descendantId

    Returns
    -------
    array
        An array of nodes with galaxies present in the now. These nodes contain
        the merger history of the galaxy in its children.

    """
    root = None
    last = None
    parent = None
    root_list = []

    keys = cosmological_model.keys()

    for index, row in cosmological_model.iterrows():
        if index % 1000 == 0:
            print(index)
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

    return root_list
