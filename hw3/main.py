from tree import TreeNode
import sys
import random

# function to calculate the difference between two input id_sequences
# simply counts mismatches and return the percentage as a float
def calculate_difference(s1,s2):
    totalCount = len(s1)
    mismatchCount = 0
    for i in xrange(totalCount):
        if s1[i] != s2[i]:
            mismatchCount += 1
    if mismatchCount == 0:
        return 0
    return mismatchCount / float(totalCount)


# function to calculate the difference matrix from ids and the id sequence mapping
# uses the calculate_difference function
# returns the matrix
def get_difference_matrix(ids, id_sequences):
    count = len(ids)
    matrix = [[0] * count for _ in xrange(count)]
    for i, id1 in enumerate(ids):
        for j, id2 in enumerate(ids):
            matrix[i][j] = calculate_difference(id_sequences[id1], id_sequences[id2])
    return matrix


# writes the difference matrix to file
def write_distance_file(ids, matrix):
    count = len(ids)
    with open('genetic_distances.txt', 'w') as f:
        f.write('\t' + '\t'.join(ids) + '\n')
        for i in xrange(count):
            f.write(ids[i] + '\t' + '\t'.join(map(str, matrix[i])) + '\n')


# uses preorder traversal to traverse the tree and generate the edges file
def write_edge_file(root):
    traversal = []
    # recursive preorder traversal
    def preorder_traversal(node):
        if not node or not node.children:
            return
        for child, distance in node.children.items():
            # traverse root first
            traversal.append((node.id, child.id, distance))
            # then traverse children recursively
            preorder_traversal(child)
    preorder_traversal(root)
    # writes the traversal list to file
    with open('edges.txt', 'w') as f:
        for (parent, child, distance) in traversal:
            f.write(parent + '\t' + child + '\t' + str(distance) + '\n')


# uses post order traversal to write the newick file
def write_newick_file(original_ids, root):
    # recursive postorder traversal
    def postorder_traversal(node):
        if not node.children:
            # if it a leaf node then return the original id
            if 1 <= int(node.id) <= 61:
                return original_ids[int(node.id) - 1]
        vals = []
        # recursively traverse the children first
        for child, distance in node.children.items():
            vals.append(postorder_traversal(child) + ':' + str(distance))
        result = '(' + ','.join(vals) + ')'
        return result
    # get the root nodes children.
    # split it up into a 2 pair and 1 extra node for the root
    (first, fd), (second, sd), (third, td) = root.children.items()
    # recursively call the postorder traversal function to generate the newick format
    first_part = '(' + postorder_traversal(second) + ':' + str(sd) + ','+ postorder_traversal(third) + ':' + str(td) + ')'
    newick =  '(' + first_part + ':' + str(fd) + ');'
    with open('tree.txt', 'w') as f:
        f.write(newick)


# function to read the fna file. Returns the ids and the dictionary
# of relationships of id -> sequence for fast lookups
def read_fasta_file(filename):
    dictionary = {}
    ids = []
    with open(filename) as f:
        content = f.read().splitlines()
    currId = None
    for index, line in enumerate(content):
        if index % 2 == 0:
            currId = line[1:]
            ids.append(currId)
        else:
            dictionary[currId] = line
    return ids, dictionary


# implemetnation of the nei saitou algorithm.
# Takes in the original ids and the distance matrix
# need original id ordering
def nei_saitou(original_ids, original_distance_matrix):
    # global root value and nodeId
    global nodeId, root
    # nodeId starts at 120 and decrements
    nodeId = 120
    root = None
    # keeps track of child parent relationships for constructing the tree later
    child_parent = {}
    # keeps track of id -> index mappings for O(1) lookup
    id_index_map = {}
    for i, id in enumerate(original_ids):
        id_index_map[id] = str(i + 1)

    # recursive helper function that computes Q matrix, q values,
    # and the new distance matrix
    def helper(ids, distance_matrix):
        global nodeId, root
        N  = len(distance_matrix)

        # ending recursive case when the matrix is size 2x2
        if N <= 2:
            # add the relationship into our dictionary and also set the root
            node_1 = id_index_map[ids[0]] if ids[0] in id_index_map else ids[0]
            node_2 = id_index_map[ids[1]] if ids[1] in id_index_map else ids[1]
            child_parent[node_1] = (node_2, distance_matrix[0][1])
            root = node_2
            return

        # default q matrix initially all 0 (NxN)
        Q_matrix = [[0] * N for _ in xrange(N)]

        # variables to keep track of the minimum indices in the q matrix
        min_i = min_j = 0
        minVal = float('inf')

        for i in xrange(N):
            for j in xrange(N):
                if i == j:
                    continue
                # based on the equation in the book
                Q_matrix[i][j] = (N - 2) * distance_matrix[i][j] - sum(distance_matrix[i]) - sum(distance_matrix[j])
                # keep track of the minimum value and indices in the Q matrix
                if Q_matrix[i][j] < minVal:
                    min_i = i
                    min_j = j
                    minVal = Q_matrix[i][j]

        # u is the new node
        # distance from min_i node to u
        # Based on equation from the book
        q_i_u = 1 / 2.0 * distance_matrix[min_i][min_j] + 1 / (2.0 * (N - 2)) * (sum(distance_matrix[min_i]) - sum(distance_matrix[min_j]))
        # distance from min_j node to u
        q_j_u = distance_matrix[min_i][min_j] - q_i_u

        # store child parent relations
        # we want to get the node id. ranging from 1-120
        child_1 = id_index_map[ids[min_i]] if ids[min_i] in id_index_map else ids[min_i]
        child_2 = id_index_map[ids[min_j]] if ids[min_j] in id_index_map else ids[min_j]
        child_parent[child_1] = (str(nodeId), q_i_u)
        child_parent[child_2] = (str(nodeId), q_j_u)
        # remove the used ids from the id list and also add the new nodeid
        ids = [e for e in ids if e not in (ids[min_i], ids[min_j])] + [str(nodeId)]
        # atomically decrement the global nodeId
        nodeId -= 1

        # create a tmp distance matrix with new node at the end (index N)
        # used to calculate distances to the new node and then used to create the new smaller
        # distance matrix
        tmp_matrix = [[0] * (N + 1) for _ in xrange(N + 1)]
        # copy everything over
        for i in xrange(N):
            for j in xrange(N):
                tmp_matrix[i][j] = distance_matrix[i][j]
        # update the distances to the new node
        for k in xrange(N):
            # based on equations in the book
            tmp_matrix[N][k] = 1 / 2.0 * (distance_matrix[min_i][k] + distance_matrix[min_j][k] - distance_matrix[min_i][min_j])
            tmp_matrix[k][N] = tmp_matrix[N][k]

        # shrink the tmp matrix to create the new distance matrix
        new_distance_matrix = [[0] * (N - 1) for _ in xrange(N - 1)]

        # where to move values from tmp values to the new_distance_matrix
        i_pointer = j_pointer = 0
        for i in xrange(N + 1):
            # if it is the minimum indices, skip to the next iteration
            if i == min_i or i == min_j:
                continue
            j_pointer = 0
            for j in xrange(N + 1):
                if j == min_i or j == min_j:
                    continue
                # set the corresponding values in the new distance matrix
                new_distance_matrix[i_pointer][j_pointer] = tmp_matrix[i][j]
                j_pointer += 1
            i_pointer += 1

        # recursively call with the new distance matrix
        helper(ids, new_distance_matrix)

    # initial call to the recursive helper function with the original params
    helper(original_ids, original_distance_matrix)

    # construct tree from the child_parent relations
    rootNode = TreeNode(root)
    # queue for breadth first search
    queue = [rootNode]

    # bfs to construct the tree
    # will end when all nodes have been traversed
    while len(queue) > 0:
        tmp = []
        for node in queue:
            # for every parent node taht matches the current node, create
            # a new child node and add it to the queue and add it as a child
            for child, (parent, distance) in child_parent.items():
                if parent == node.id:
                    childNode = TreeNode(child)
                    node.add_child(childNode, distance)
                    tmp.append(childNode)
        # delete the node in the dictionary so we dont reuse it
        for node in tmp:
            del child_parent[node.id]
        # the queue is set to the next level
        queue = tmp

    # return the root of the tree
    return rootNode


# get the leaves under this node using dfs
# returns a set of node.ids
def get_leaves(node):
    leaves = set()
    def dfs(node):
        global count
        if not node or not node.children:
            leaves.add(node.id)
        else:
            for child in node.children.keys():
                dfs(child)
    dfs(node)
    return leaves


# return the dictionary of id -> set of leaves under that node
# using bfs with a queue
def get_partitions(root):
    queue = [root]
    partition = {}
    while len(queue):
        tmp = []
        for node in queue:
            # this is a leaf
            if not node.children:
                continue
            leaves = get_leaves(node)
            partition[node.id] = leaves
            for child in node.children.keys():
                tmp.append(child)
        queue = tmp
    return partition


# get the dfs order of the tree
def get_id_order(root):
    order = []
    def dfs(node):
        if not node or not node.children:
            return
        else:
            order.append(node.id)
            for child in node.children.keys():
                dfs(child)
    dfs(root)
    return order


# do the 100 bootstrap sampling
def bootstrap(original_root, ids, id_sequences):
    # original order of teh ids in the original tree (dfs)
    original_order = get_id_order(original_root)
    original_partitions = get_partitions(original_root)
    partition_count = [0] * 59
    for _ in xrange(100):
        # new bootstrap sequence
        bootstrap_sequences = {}
        # the columns to swap
        all_indices = [i for i in xrange(len(id_sequences['152801']))]
        indices = [random.choice(all_indices) for _ in xrange(len(all_indices))]
        for id in ids:
            # choose randomly the to reconstruct the sequence
            new_sequence = ''
            for index in indices:
                new_sequence += id_sequences[id][index]
            bootstrap_sequences[id] = new_sequence
        # do the nei_saitou and get a new tree
        distance_matrix = get_difference_matrix(ids, bootstrap_sequences)
        root = nei_saitou(ids, distance_matrix)
        # get the dictionary of partitions and ids under those partitions
        partitions = get_partitions(root)

        # compare partitions to see which trees make the same partition as the original
        for index, id in enumerate(original_order):
            if original_partitions[id] == partitions[id]:
                partition_count[index] += 1
    percentages = [count / 100.0 for count in partition_count]
    return percentages

# write the percentages in the bootstrap file
def write_bootstrap(percentages):
    with open('bootstrap.txt', 'w') as f:
        for percent in percentages:
            f.write(str(percent) + '\n')


def main(filename):
    # read the fna file and return the ids and id sequence mappings
    ids, id_sequences = read_fasta_file(filename)

    # generate the difference matrix
    distance_matrix = get_difference_matrix(ids, id_sequences)

    # write the distances.txt file (distance matrix)
    write_distance_file(ids, distance_matrix)

    # run the nei saitu algorithm. Returns the root of the tree
    root = nei_saitou(ids, distance_matrix)

    # write the edge file using preorder traversal
    write_edge_file(root)

    # write the newick file using postorder traversal
    write_newick_file(ids, root)

    # bootstrap bonus points
    percentages = bootstrap(root, ids, id_sequences)
    write_bootstrap(percentages)


if __name__ == '__main__':
    main(sys.argv[1])
