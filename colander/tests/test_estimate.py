import pytest
import networkx as nx
from colander.estimate import (
    GreedyError,
    CycleSet,
    Cycle,
    get_cov,
    get_max_weight_edge_from_node,
)


def get_test_graph():
    G = nx.DiGraph()
    G.add_edges_from(
        [
            (0, 1, {"cov": 40}),
            (1, 2, {"cov": 20}),
            (1, 3, {"cov": 10}),
            (1, 4, {"cov": 10}),
            (2, 5, {"cov": 20}),
            (3, 5, {"cov": 10}),
            (4, 5, {"cov": 10}),
            (5, 6, {"cov": 40}),
            (6, 0, {"cov": 40}),
        ]
    )
    return G


def test_get_cov():
    g = get_test_graph()
    assert get_cov(g, (1, 3)) == 10
    assert get_cov(g, (2, 5)) == 20


def test_cycle_set_len():
    c = Cycle(((0, 1), (1, 2), (2, 5), (5, 6), (6, 0)), 33)
    cs = CycleSet()
    assert len(cs) == 0
    cs.add(c)
    assert len(cs) == 1


def test_cycle_set_cycle_coverage():
    c = Cycle(((0, 1), (1, 2), (2, 5), (5, 6), (6, 0)), 33)
    cs = CycleSet()
    assert cs.cycle_coverage((0, 1)) == 0
    cs.add(c)
    print(cs.cycles)
    assert cs.cycle_coverage((0, 1)) == 33


def test_cycle_set_conformity_score():
    c = Cycle(((0, 1), (1, 2), (2, 5), (5, 6), (6, 0)), 33)
    cs = CycleSet()
    g = get_test_graph()
    assert cs.conformity_score(g) == 6000
    cs.add(c)
    assert cs.conformity_score(g) == 885


def test_get_max_weight_edge_from_node():
    g = get_test_graph()
    assert get_max_weight_edge_from_node(g, 1, []) == (1, 2)
    assert get_max_weight_edge_from_node(g, 1, [(1, 2)]) in ((1, 3), (1, 4))


def test_trigger_greedyerror():
    g = get_test_graph()
    covdict = {}
    for e in g.edges:
        covdict[e] = {"cov": 0}
    nx.set_edge_attributes(g, covdict)
    with pytest.raises(GreedyError):
        get_max_weight_edge_from_node(g, 1, [])
