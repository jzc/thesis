import networkx as nx
import networkx.algorithms.isomorphism as iso
import random
from collections import deque
import heapq 
from datetime import datetime
from itertools import product, combinations


class iset(frozenset):
    def __repr__(self):
        return "{" f"{', '.join(repr(i) for i in self)}" "}"


def is_labeled(G):
    """Returns True if each edge in G has a label attribute."""
    return all(label is not None for _, _, label in G.edges(data="label"))


def _outgoing_labels_list(G, I):
    return [label for (_, _, label) in G.out_edges(I, data="label")]


def outgoing_labels(G, I):
    return iset(_outgoing_labels_list(G, I))


def is_deterministic(G):
    """Returns True if G is deterministic """
    for I in G:
        if len(outgoing_labels(G, I)) != len(_outgoing_labels_list(G, I)):
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
    
    for a in w:
        d = {b: J for _, J, b in G.out_edges(I, data="label")}
        if a not in d:
            return None
        else:
            I = d[a]
        
    return I


def make_inequivalence_table(G):
    marked = {}
    for (I, J) in combinations(G, 2):
        marked[iset([I, J])] = outgoing_labels(G, I) != outgoing_labels(G, J)

    while True:
        new_marks = False
        for s in marked:
            if not marked[s]:
                I, J = s
                x = outgoing_labels(G, I)
                x = (iset([transition(G, I, a), transition(G, J, a)]) for a in x)
                x = (marked.get(sp, False) for sp in x)
                if any(x):
                    marked[s] = True
                    new_marks = True
            
        if not new_marks:
            break
    
    return marked
            

def follower_separate(G):
    marked = make_inequivalence_table(G)
           
    vertices = set(G)
    eq_classes = []
    eq_map = {}
    while vertices:
        it = iter(vertices)
        a = next(it)
        eq_class = set([a])
        for b in it:
            if not marked[iset([a, b])]:
                eq_class.add(b)
        vertices -= eq_class  
        eq_class = iset(eq_class) 
        eq_classes.append(eq_class)
        for e in eq_class:
            eq_map[e] = eq_class
            
    Gp = nx.MultiDiGraph()
    Gp.add_nodes_from(eq_classes)
    for I in Gp:
        for _, J, a in G.out_edges(I, data="label"):
            if a not in outgoing_labels(Gp, I):
                Gp.add_edge(I, eq_map[J], label=a)
            
    return Gp


def is_follower_separated(G):
    return all(make_inequivalence_table(G).values())


def is_label_isomorphic(G1, G2):
    return iso.is_isomorphic(G1, G2, edge_match=iso.categorical_multiedge_match("label", None))


def is_isomorphic(G1, G2):
    return iso.is_isomorphic(G1, G2)


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


def transition_subset(G, ss, w):
    ssp = (transition(G, I, w) for I in ss)
    ssp = (I for I in ssp if I is not None)
    return iset(ssp)


def subset_construction(G, alphabet=None):
    if alphabet == iset():
        alphabet = iset(l for _, _, l in G.edges(data="label"))

    Gp = nx.MultiDiGraph()
    init = iset(list(G))
    Gp.add_node(init)
    visited = set([init])
    q = deque([init])
    while q:
        current = q.popleft()
        it = (outgoing_labels(G, current) if alphabet is None else alphabet)
        for a in it:
            new = transition_subset(G, current, a)
            if new not in visited:
                q.append(new)
                visited.add(new)
            Gp.add_edge(current, new, label=a)

    return Gp


def is_GH(G):
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


def get_GH(GH):
    sccs = list(nx.strongly_connected_components(GH))
    H1, H2 = GH.subgraph(sccs[0]), GH.subgraph(sccs[1])
    edges_H1_H2 = len([v for u, v in GH.out_edges(H1) if v in H2])
    edges_H2_H1 = len([v for u, v in GH.out_edges(H2) if v in H1])
    if (edges_H1_H2 == 1 and edges_H2_H1 == 0):
        return H1, H2
    elif (edges_H2_H1 == 1 and edges_H1_H2 == 0):
        return H2, H1
    else:
        return None

    
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


def find_shortest_GH_word(GH):
    res = get_GH(GH)
    assert res is not None
    G, H = res
    alphabet = iset(l for _, _, l in GH.edges(data="label"))
    P = label_product(
        subset_construction(GH),
        subset_construction(H, alphabet=alphabet)
    )
    # print(list(P))
    source = (iset(GH), iset(H))
    for u, v in nx.bfs_edges(P, source):
        P.nodes[v]["pred"] = u
    
    ws = {}
    for s in (s for s in P if s[1] == iset()):
        c = s
        path = []
        no_path = False
        while c != source:
            path.append(c)
            if "pred" not in P.nodes[c]:
                no_path = True
                break
            else:
                c = P.nodes[c]["pred"]

        if no_path:
            continue

        path.append(source)
        path = list(reversed(path))
        w = []
        for u, v in zip(path, path[1:]):
            w.append(first(P[u][v].values())["label"])
        ws[s] = ''.join(w)
    
    return ws


def make_scc_dag(G):
    Gp = nx.DiGraph()
    sccs = [iset(c) for c in nx.strongly_connected_components(G)]
    Gp.add_nodes_from(sccs)
    for ci in sccs:
        for cj in sccs:
            if ci != cj:
                edges_between_ci_cj = any(True for u, v in G.out_edges(ci) if v in cj)
                if edges_between_ci_cj:
                    Gp.add_edge(ci, cj, edges_between=edges_between_ci_cj)

                edges_between_cj_ci = any(True for u, v in G.out_edges(cj) if v in ci)
                if edges_between_cj_ci:
                    Gp.add_edge(cj, ci, edges_between=edges_between_cj_ci)

    return Gp


def add_kill_state(G, alphabet=None):
    Gp = nx.MultiDiGraph(G)
    Gp.add_node("K")
    
    if alphabet is None:
        alphabet = {l for _, _, l in G.edges(data="label")}

    for I in Gp:
        out_labels = outgoing_labels(G, I)
        for a in (alphabet - out_labels):
            Gp.add_edge(I, "K", label=a)
        
    return Gp


def label_product(G, H):
    GH = nx.MultiDiGraph()
    GH.add_nodes_from(product(G, H))

    for (I, J) in GH:
        for a in (outgoing_labels(G, I) & outgoing_labels(H, J)):
            Ip = transition(G, I, a)
            Jp = transition(H, J, a)
            assert Ip is not None
            assert Jp is not None
            GH.add_edge((I, J), (Ip, Jp), label=a)

    return GH  


def first(it):
    return next(iter(it))


def is_subshift(G, H):
    alphabet = {l for _, _, l in G.edges(data="label")} | {l for _, _, l in H.edges(data="label")}
    Hk = add_kill_state(H, alphabet=alphabet)
    GH = label_product(G, Hk)
    paths = nx.shortest_path(GH)
    I = first(G)
    X = list(H)
    w = "" 
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
                    return False  
                break

        if not path_exists:
            break

    return True


def find_shortest_separating_words(G, I, J):
    P = label_product(G, add_kill_state(G))
    paths = nx.shortest_path(P)
    paths = [paths[I, J][A, "K"] for A in G if (A, "K") in paths[I, J]]
    if paths:
        min_len = min([len(p) for p in paths])
        paths = [p for p in paths if len(p) == min_len]
        all_words = []
        for path in paths:
            words = [edge["label"] for edge in P[path[0]][path[1]].values()]
            for i in range(1, len(path)-1):
                labels = [edge["label"] for edge in P[path[i]][path[i+1]].values()]
                words = [word+label for word in words for label in labels]
            all_words.extend(words)        

        return all_words


def choice_graph(GH):
    res = get_GH(GH)
    assert res is not None
    G, H = res

    K = nx.MultiDiGraph()

    init = [(I, iset(H)) for I in G]
    K.add_nodes_from(init)
    q = deque(init)
    visited = set(init)

    while q:
        current = q.popleft()
        I, X = current
        it = (find_shortest_separating_words(GH, I, J) for J in X)
        it = (ws for ws in it if ws is not None)
        it = (w for ws in it for w in ws)
        for w in it:
            Ip = transition(GH, I, w)
            Xp = transition_subset(GH, X, w)
            new = (Ip, Xp)
            if new not in visited:
                q.append(new)
                visited.add(new)
            K.add_edge(current, new, label=w)

    return K