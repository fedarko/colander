import pytest
from colander.mock_data_generation.utils import (
    rand_nt,
    generate_random_sequence,
    Mutation,
)


def test_rand_nt():
    for i in range(100):
        n = rand_nt()
        assert n in ["A", "C", "G", "T"]


def test_rand_nt_exclude():
    for i in range(100):
        n = rand_nt(exclude="A")
        assert n in ["C", "G", "T"]


def test_generate_random_sequence():
    for i in range(100):
        seq = generate_random_sequence(i)
        assert len(seq) == i
        for nt in seq:
            assert nt in ["A", "C", "G", "T"]


def test_mutation():
    nts = ["A", "C", "G", "T"]
    for i in range(100):
        c = nts[i // 25]
        m = Mutation(i, c)
        if m.mtype == "insertion":
            assert m.__repr__() == (
                "insertion at pos {}: {} -> {}{}".format(i, c, m.new_nt, c)
            )
        elif m.mtype == "deletion":
            assert m.__repr__() == (
                "deletion at pos {}: {} -> ''".format(i, c)
            )
        else:
            assert m.__repr__() == (
                "mutation at pos {}: {} -> {}".format(i, c, m.new_nt)
            )
            assert m.new_nt != c


def test_invalid_mutation():
    m = Mutation(3, "A")
    m.mtype = "lol"
    with pytest.raises(ValueError):
        m.__repr__()
