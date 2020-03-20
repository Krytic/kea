#
# Test for the merger history tests
#
#
from kea.mergerHistory import node, pprint_tree

def test_print():

    x = node()
    assert x.data == {}

def test_node():

    x = node()
    data = {"snapnum":10, "galaxyID":11}

    x.addData(data)

    assert x.data["snapnum"] == 10
    assert x.galaxyID == 11


def test_children():

    root = node()
    child1 = node()
    child2 = node()

    root.addChild(child1)
    root.addChild(child2)

    data = {"snapnum":10}
    root.addData(data)
    child1.addData(data)
    child2.addData(data)

    assert pprint_tree(root)  == None
