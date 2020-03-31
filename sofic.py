import networkx as nx
import networkx.algorithms.isomorphism as iso

def is_labeled(G):
    """Returns True if each edge in G has a label attribute."""
    return all(label is not None for _, _, label in G.edges(data="label"))


def is_deterministic(G):
    """Returns True if G is deterministic """
    for vertex in G:
        labels = [label for (_, _, label) in G.out_edges(vertex, data="label")]
        if len(set(labels)) != len(labels):
            return False
    return True


def is_stranded(G, I):
    """Returns true if I is stranded (i.e. there are no edges starting at I or no edges ending at I)."""
    return (not G.out_edges(I)) or (not G.in_edges(I))


def is_essential(G):
    """Returns true if G is essential (i.e. no vertex in G is stranded)."""
    return all(not is_stranded(G, I) for I in G)


def trim(G):
    """Removes stranded vertices from G."""
    G = nx.MultiDiGraph(G)
    G.remove_nodes_from([I for I in G if is_stranded(G, I)])
    return G


def transition(G, I, w):
    """
    Given a deterministic graph G, run the transition associated 
    with the word w starting at I.

    If the transition is undefined, this returns None.
    
    If I is None, this returns None.
    """

    if I is None:
        return None
    
    for c in w:
        defined = False
        for (_, J, l) in G.out_edges(I, data="label"):
            if l == c:
                I = J
                defined = True
                break
        if not defined:
            return None
        
    return I


def table_generator(l):
    for i in range(len(l)):
        for j in range(i+1, len(l)):
            yield l[i], l[j]


def table_ordering(l, x, y):
    x_idx = l.index(x)
    y_idx = l.index(y)
    i, j = sorted([x_idx, y_idx])
    return l[i], l[j]


def minimize(g):
    g = trim(g)
    vertex_list = list(g)
    o = lambda x, y: table_ordering(vertex_list, x, y)
    marked = {}
    for i, j in table_generator(vertex_list):
        i_out_edges = g.out_edges(i, data="label")
        j_out_edges = g.out_edges(j, data="label")
        if set(l for _, _, l in i_out_edges) != set(l for _, _, l in j_out_edges):
            marked[i, j] = True
        else:
            marked[i, j] = False

    while True:
        new_marks = False
        
        for i, j in marked:
            if not marked[i, j]:        
                labels = [l for _, _, l in g.out_edges(i, data="label")]
                if any(marked.get(o(transition(g, i, l), transition(g, j, l)), True) for l in labels):
                    marked[i, j] = True
                    new_marks = True               
        
        if not new_marks:
            break
           
    if all(marked.values()):
        return g
        
    gp = nx.MultiDiGraph()
    gp.add_nodes_from([ij for ij, m in marked.items() if not m])
    for (i, j) in gp:
        for _, _, l in g.out_edges(i, data="label"):
            gp.add_edge((i, j), o(transition(g, i, l), transition(g, j, l)), label=l)
            
    return gp


def is_follower_separated(g):
    return is_label_isomorphic(g, minimize(g))


def is_label_isomorphic(g1, g2):
    return iso.is_isomorphic(g1, g2, edge_match=iso.categorical_multiedge_match("label", None))


def is_isomorphic(g1, g2):
    return iso.is_isomorphic(g1, g2)