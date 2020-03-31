from sofic import *

G1 = nx.MultiDiGraph()
G1.add_nodes_from([1, 2])
G1.add_edge(1, 2)

# minimal even shift
G2 = nx.MultiDiGraph()
G2.add_nodes_from([1, 2])
G2.add_edge(1, 1, label="1")
G2.add_edge(1, 2, label="0")
G2.add_edge(2, 1, label="0")

# non-minimal even shift
G3 = nx.MultiDiGraph()
G3.add_nodes_from([1, 2, 3, 4])
G3.add_edge(1, 1, label="1")
G3.add_edge(1, 2, label="0")
G3.add_edge(2, 3, label="0")
G3.add_edge(3, 3, label="1")
G3.add_edge(3, 4, label="0")
G3.add_edge(4, 1, label="0")

# reducible, follower-separated even shift
G4 = nx.MultiDiGraph(G2)
G4.add_node(3)
G4.add_edge(3, 3, label="0")
G4.add_edge(3, 1, label="1")

# non-deterministic even shift
G5 = nx.MultiDiGraph()
G5.add_nodes_from([1, 2, 3])
G5.add_edge(1, 1, label="1")
G5.add_edge(1, 2, label="1")
G5.add_edge(2, 3, label="0")
G5.add_edge(3, 2, label="0")
G5.add_edge(3, 1, label="0")


def test_is_labeled():
    assert not is_labeled(G1)
    assert is_labeled(G2)
    assert is_labeled(G3)
    assert is_labeled(G4)
    assert is_labeled(G5)


def test_is_deterministic():
    assert is_deterministic(G2)
    assert is_deterministic(G3)
    assert is_deterministic(G4)
    assert not is_deterministic(G5)


def test_is_essential():
    assert not is_essential(G1)
    assert is_essential(G2)
    assert is_essential(G3)
    assert is_essential(G4)
    assert is_essential(G5)


def test_trim():
    assert len(trim(G1)) == 0
    assert len(trim(G2)) == len(G2)


def test_transition():
    assert transition(G2, 1, "1") == 1
    assert transition(G2, 1, "10") == 2
    assert transition(G2, 1, "100") == 1
    assert transition(G2, 1, "10000") == 1
    assert transition(G2, 2, "0") == 1
    assert transition(G2, 2, "00") == 2
    assert transition(G2, 2, "1") is None
    assert transition(G2, 1, "2") is None
    

def test_minimize():
    assert is_label_isomorphic(G2, minimize(G3))


def test_is_follower_separated():
    assert is_follower_separated(G2)
    assert not is_follower_separated(G3)
    assert is_follower_separated(G4)