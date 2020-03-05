import pytest
from colander.mock_data_generation.utils import (
    rand_nt,
    generate_random_sequence,
    Mutation,
    generate_mutations,
    add_mutations,
    generate_strain,
    generate_strains_from_genome,
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
        m.randomly_initialize()
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


def test_indeterminate_mutation():
    m = Mutation(3, "A")
    assert m.__repr__() == "indeterminate mutation"


def test_add_mutations():
    s = "ACGTAAAC"
    # Deterministically create mutations
    m1 = Mutation(0, "A")
    m1.make_insertion("T")

    m2 = Mutation(3, "T")
    m2.make_mutation("C")

    m3 = Mutation(7, "C")
    m3.make_deletion()
    assert add_mutations(s, [m1, m2, m3]) == "TACGCAAA"


def test_add_indeterminate_mutation():
    m = Mutation(0, "A")
    with pytest.raises(ValueError):
        add_mutations("ACTG", [m])


def test_gen_mutations():
    g = "ATCGAACGATAAACTAGACCCAA"
    hv_regions = [(3, 9)]
    hvmp = 1
    nmp = 0
    muts = generate_mutations(g, hv_regions, hvmp, nmp)
    assert [m.coordinate for m in muts] == list(range(3, 10))


def test_generate_strain():
    g = "AAAA"
    hv_regions = [(0, 1)]
    hvmp = 1
    nmp = 0
    for i in range(100):
        strain = generate_strain(g, 5, hv_regions, hvmp, nmp)
        # Since there's a normal mutation probability of 0, the areas outside
        # the hypervariable region (the first "AA") should be unmodified
        assert strain.seq.endswith("AA")
        assert strain.coverage == 5
        # Since there's a hypervariable region mutation probability of 1, there
        # should have been 2 mutations for each of the 2 bases in the HV region
        assert len(strain.originating_mutations) == 2


def test_generate_strains_from_genome():
    #     ***    -- the *s indicate the hypervariable region for this "genome"
    #    012345
    g = "GGGGGG"
    covs = [10, 1, 5]
    hv_regions = [(1, 3)]
    hvmp = 1
    nmp = 0
    for i in range(100):
        strains = generate_strains_from_genome(
            g,
            covs,
            hv_regions,
            hypervariable_mutation_probability=hvmp,
            normal_mutation_probability=nmp,
        )
        for s, cov in zip(strains, covs):
            assert s.seq.endswith("GG")
            assert len(s.originating_mutations) == 3
            assert s.coverage == cov
