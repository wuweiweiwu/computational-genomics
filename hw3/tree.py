# class to represent a node in a Tree
# a node has its own id and a map of child nodes and respective distances
class TreeNode:
    def __init__(self, id):
        self.id = id
        self.children = None

    def add_child(self, node, distance):
        if not self.children:
            self.children = {}
        self.children[node] = distance
