import numpy
import numba

# gah! nubma doesn't support recursion!
@numba.jit
def traverse_paths(
    node, dest_nodes, visited, npaths, path_nodes, path_length, n_edges, outgoing_edges
):
    visited[node] = True
    path_nodes[path_length] = node

    if dest_nodes[node]:
        # we have reached a destination node
        # increment the npaths counts for all the nodes in this path
        for i in range(path_length + 1):
            npaths[path_nodes[i]] += 1
    else:
        for i in range(n_edges[node]):
            neighb = outgoing_edges[node,i]
            if visited[neighb]:
                # this really shouldn't happen, right???
                continue
            traverse_paths(
                neighb,
                dest_nodes,
                visited,
                npaths,
                path_nodes,
                path_length + 1,
                n_edges,
                outgoing_edges,
            )
    visited[node] = False
    path_nodes[path_length] = -1


@numba.jit(nopython=True)
def traverse_paths_iterative(
    start_node,
    dest_nodes,
    visited,
    npaths,
    path_nodes,
    neighbor_count,
    n_edges,
    outgoing_edges,
):
    stack_depth = 0
    path_nodes[stack_depth] = start_node
    neighbor_count[stack_depth] = 0
    while stack_depth >= 0:
        node = path_nodes[stack_depth]
        visited[node] = True
        if dest_nodes[node]:
            for i in range(1, stack_depth + 1):
                npaths[path_nodes[i]] += 1

        if neighbor_count[stack_depth] < n_edges[node]:
            # descend
            last_stack_depth = stack_depth
            stack_depth += 1
            neighbor_count[stack_depth] = 0
            path_nodes[stack_depth] = outgoing_edges[
                node,
                neighbor_count[last_stack_depth]
            ]
            neighbor_count[last_stack_depth] += 1
        else:
            # ascend
            last_stack_depth = stack_depth
            stack_depth -= 1
            visited[node] = False
            path_nodes[last_stack_depth] = -1


@numba.jit(nopython=True)
def paths_to_destinations(
    start_node: int,
    edge_matrix: numpy.array,  # N x N 2D array of 0s and 1s
    dest_nodes: numpy.array,  # N 1D array of 0s and 1s
):
    N = edge_matrix.shape[0]
    visited = numpy.zeros((N,), dtype=numpy.int32)
    npaths = numpy.zeros((N,), dtype=numpy.int32)
    path_nodes = numpy.zeros((N,), dtype=numpy.int32)

    n_edges = numpy.zeros((N,), dtype=numpy.int32)
    outgoing_edges = numpy.zeros((N, N), dtype=numpy.int32)

    for i in range(N):
        for j in range(N):
            if edge_matrix[i, j]:
                outgoing_edges[i, n_edges[i]] = j
                n_edges[i] += 1

    # for i in range(n_edges[start_node]):
    #     neighb = outgoing_edges[start_node, i]
    #     traverse_paths(
    #         neighb, dest_nodes, visited, npaths, path_nodes, 0, n_edges, outgoing_edges
    #     )
    neighbor_count = numpy.zeros((N,), dtype=numpy.int32)
    print("traverse paths iterative")
    traverse_paths_iterative(
        start_node, dest_nodes, visited, npaths, path_nodes, neighbor_count, n_edges, outgoing_edges)
    

    return npaths


if __name__ == "__main__":
    edge_matrix = numpy.zeros((10, 10), dtype=numpy.int32)
    edge_matrix[0, 1] = 1
    edge_matrix[0, 2] = 1
    edge_matrix[0, 6] = 1
    edge_matrix[1, 4] = 1
    edge_matrix[1, 5] = 1
    edge_matrix[3, 9] = 1
    edge_matrix[4, 7] = 1
    edge_matrix[4, 8] = 1
    edge_matrix[5, 8] = 1
    edge_matrix[6, 8] = 1

    dest_nodes = numpy.zeros((10,), dtype=numpy.int32)
    dest_nodes[7] = 1
    dest_nodes[8] = 1
    dest_nodes[9] = 1

    npaths = paths_to_destinations(0, edge_matrix, dest_nodes)
    #print("npaths!", npaths)
    gold_npaths = numpy.array([0, 3, 0, 0, 2, 1, 1, 1, 3, 0], dtype=numpy.int32)
    numpy.testing.assert_equal(npaths, gold_npaths)
