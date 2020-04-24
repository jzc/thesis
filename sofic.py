import networkx as nx
import networkx.algorithms.isomorphism as iso
import random
from collections import deque
import heapq 
from datetime import datetime
from itertools import product


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
        i_outgoing_labels = set(l for _, _, l in g.out_edges(i, data="label"))
        j_outgoing_labels = set(l for _, _, l in g.out_edges(j, data="label"))
        
        # print(i_outgoing_labels, j_outgoing_labels, i_outgoing_labels==j_outgoing_labels)
        if i_outgoing_labels != j_outgoing_labels:
            marked[i, j] = True
        else:
            marked[i, j] = False

    # print(marked)
    while True:
        new_marks = False
        
        for i, j in marked:
            if not marked[i, j]:        
                labels = [l for _, _, l in g.out_edges(i, data="label")]
                # print([o(transition(g, i, l), transition(g, j, l)) for l in labels])
                if any(marked.get(o(transition(g, i, l), transition(g, j, l)), False) for l in labels):
                    marked[i, j] = True
                    new_marks = True               
        
        if not new_marks:
            break
           
    # print("eq")
    vertices = set(vertex_list)
    eq_classes = []
    eq_map = {}
    while vertices:
        it = iter(vertices)
        a = next(it)
        eq_class = set([a])
        for b in it:
            if not marked[o(a, b)]:
                eq_class.add(b)
        vertices -= eq_class  
        eq_class = iset(eq_class) 
        eq_classes.append(eq_class)
        for e in eq_class:
            eq_map[e] = eq_class
    # print(eq_classes)

    gp = nx.MultiDiGraph()
    gp.add_nodes_from(eq_classes)
    for I in gp:
        for i, j, l in g.out_edges(I, data="label"):
            if l not in {l for _, _, l in gp.out_edges(I, data="label")}:
                gp.add_edge(I, eq_map[j], label=l)
            
    return gp


def is_follower_separated(g):
    return is_label_isomorphic(g, minimize(g))


def is_label_isomorphic(g1, g2):
    return iso.is_isomorphic(g1, g2, edge_match=iso.categorical_multiedge_match("label", None))


def is_isomorphic(g1, g2):
    return iso.is_isomorphic(g1, g2)


def is_irreducible(G):
    return nx.is_strongly_connected(G)


def is_reducible(G):
    return is_essential(G) and not is_irreducible(G)


def graph_from_partial_fns(pfns):
    G = nx.MultiDiGraph()
    first = next(iter(pfns.values()))
    G.add_nodes_from(range(len(first)))

    for l, pfn in pfns.items():
        for I, J in enumerate(pfn):
            if J is not None:
                G.add_edge(I, J, label=l)

    return G


def random_graph_with_props(n, m, props):
    while True:
        G = random_graph(n, m)
        if all(prop(G) for prop in props):
            return G


def random_graph(n, m):
    pfns = {}

    for l in range(m):
        l = str(l)
        pfn = (random.randrange(-1, n) for _ in range(n))
        pfn = (x if x != -1 else None for x in pfn)
        pfn = list(pfn)
        pfns[l] = pfn

    return graph_from_partial_fns(pfns)


class iset(frozenset):
    def __repr__(self):
        return "{" f"{', '.join(repr(i) for i in self)}" "}"


def transition_subset(G, ss, w):
    ssp = (transition(G, I, w) for I in ss)
    ssp = (I for I in ssp if I is not None)
    return iset(ssp)


def subset_construction(G):
    Gp = nx.MultiDiGraph()
    init = iset(list(G))
    Gp.add_node(init)
    visited = set([init])
    q = deque([init])
    while q:
        current = q.popleft()
        # print(current)
        outgoing_labels = iset(l for _, _, l in G.out_edges(current, data="label"))
        # print(outgoing_labels)
        for l in outgoing_labels:
            # new = iset(j for _, j, lp, in G.out_edges(current, data="label") if l == lp)
            # print(l, new)
            new = transition_subset(G, current, l)
            if new not in visited:
                q.append(new)
                visited.add(new)
            Gp.add_edge(current, new, label=l)

    return Gp


def astar(G, source, target_set, h):
    visited = set([source])
    q = []
    heapq.heappush(q, (h(source), source))
    g = {source: 0}
    f = {source: h(source)}
    prev = {}
    n = 0

    while q:
        n += 1
        _, current = heapq.heappop(q)
        if current in target_set:
            return prev, current, n

        visited.add(current)

        for _, succ, l in G.out_edges(current, data="label"):
            g_succ = g[source] + 1
            if g_succ < g.get(succ, float("inf")):
                prev[succ] = current
                g[succ] = g_succ
                f[succ] = g_succ + h(succ)
                if succ not in visited:
                    heapq.heappush(q, (f[succ], succ))


def reconstruct_path(prev, target):
    p = []
    while target in prev:
        p.append(target)
        target = prev[target]
    
    p = list(reversed(p))
    x = (IJ for IJ in zip(p, p[1:]))
    x = ([l for _, K, l in G_ssc.out_edges(I, data="label") if K == J] for I, J in x)
    x = (ls[0] for ls in x)
    return ''.join(x)


def find_synchronizing_word(G, init=None):
    if init is None:
        init = iset(list(G))
    G_ssc = subset_construction(G)
    prev, target, _ = astar(G_ssc, init, {iset([i]) for i in G}, len)

    return reconstruct_path(prev, target)


def is_2_simple(G):
    sccs = list(nx.strongly_connected_components(G))
    if len(sccs) != 2:
        return False

    H1, H2 = G.subgraph(sccs[0]), G.subgraph(sccs[1])
    edges_H1_H2 = len([v for u, v in G.out_edges(H1) if v in H2])
    edges_H2_H1 = len([v for u, v in G.out_edges(H2) if v in H1])
    if ((edges_H1_H2 == 1 and edges_H2_H1 == 0) or 
        (edges_H2_H1 == 1 and edges_H1_H2 == 0)):
       return True

    return False

    
def write_graph(G):
    name = datetime.now().strftime("%Y-%m-%d--%H-%M-%S.xml")
    nx.write_graphml(G, f"graphs/{name}")


def find_shortest_synchronizing_word(G):
    G_ssc = subset_construction(G)
    source = iset(list(G))

    for u, v in nx.bfs_edges(G_ssc, source):
        G_ssc.nodes[v]["pred"] = u

    ws = {}
    for singleton in (s for s in G_ssc if len(s) == 1):
        c = singleton
        path = []
        no_path = False
        while c != source:
            path.append(c)
            if "pred" not in G_ssc.nodes[c]:
                no_path = True
                break
            else: 
                c = G_ssc.nodes[c]["pred"]

        if no_path:
            continue

        path.append(source)
        path = list(reversed(path))
        w = []
        for u, v in zip(path, path[1:]):
            w.append(next(iter(G_ssc[u][v].values()))["label"])
        ws[singleton] = ''.join(w)

    return ws


def make_scc_dag(G):
    Gp = nx.DiGraph()
    sccs = [iset(c) for c in nx.strongly_connected_components(G)]
    Gp.add_nodes_from(sccs)
    for ci in sccs:
        for cj in sccs:
            if ci != cj:
                edges_between_ci_cj = [(u, v) for u, v in G.out_edges(ci) if v in cj]
                if edges_between_ci_cj:
                    Gp.add_edge(ci, cj, edges_between=edges_between_ci_cj)

                edges_between_cj_ci = [(u, v) for u, v in G.out_edges(cj) if v in ci]
                if edges_between_cj_ci:
                    Gp.add_edge(cj, ci, edges_between=edges_between_cj_ci)

    return Gp

def add_kill_state(G, alphabet=None):
    Gp = nx.MultiDiGraph(G)
    Gp.add_node("K")
    
    if alphabet is None:
        alphabet = {l for _, _, l in G.edges(data="label")}

    for I in Gp:
        out_labels = {l for _, _, l in G.edges(I, data="label")}
        for a in (alphabet - out_labels):
            Gp.add_edge(I, "K", label=a)
        
    return Gp


def label_product(G, H):

    GH = nx.MultiDiGraph()
    GH.add_nodes_from(product(G, H))

    for (I, J) in GH:
        I_labels = {l for _, _, l in G.out_edges(I, data="label")}
        J_labels = {l for _, _, l in H.out_edges(J, data="label")}
        for l in (I_labels & J_labels):
            Ip = transition(G, I, l)
            Jp = transition(H, J, l)
            assert Ip is not None
            assert Jp is not None
            GH.add_edge((I, J), (Ip, Jp), label=l)

    return GH  


def first(it):
    return next(iter(it))


def is_subshift(G, H):
    alphabet = {l for _, _, l in G.edges(data="label")} | {l for _, _, l in H.edges(data="label")}
    Hk = add_kill_state(H, alphabet=alphabet)
    GH = label_product(G, Hk)
    paths = nx.shortest_path(GH)

    for I0 in G:
        X = list(H)
        w = ""
        I = I0
        while True:
            J = first(X)
            path_exists = False
            for A in G:
                if (A, "K") in paths[I, J]:
                    # there is a path starting at I labeled wp that 
                    # takes J to the kill state 
                    path_exists = True
                    path = paths[I, J][A, "K"]
                    wp = "".join([first(GH[a][b].values())["label"] for a, b in zip(path, path[1:])])
                    I = transition(G, I, wp)
                    X = transition_subset(H, X, wp)
                    w += wp
                    if len(X) == 0:
                        # w in language of G but not in language of H
                        print(I)
                        return w  
                    break

            if not path_exists:
                break

    return True

def sigma_star(sigma=["0", "1"]):
    c = 1
    while True:
        for s in product(sigma, repeat=c):
            yield "".join(s)
        c += 1

def enumerate_language(G):
    for w in sigma_star():
        if transition_subset(G, G, w):
            yield w

def enumerate_complement(G):
    for w in sigma_star():
        if not transition_subset(G, G, w):
            yield w

def take(x, n):
    return [i for i, _ in zip(x, range(n))]

# def minimize_automata(G, q0, F):  
