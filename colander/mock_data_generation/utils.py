import random


def rand_nt():
    return random.choice(["A", "C", "T", "G"])


class Strain:
    def __init__(self, seq, coverage):
        self.seq = seq
        self.coverage = coverage


class Mutation:
    def __init__(self, coordinate, curr_nt):
        self.coordinate = coordinate
        self.curr_nt = curr_nt
        self.mtype = random.choice(["insertion", "deletion", "mutation"])
        if self.mtype == "insertion" or self.mtype == "mutation":
            self.new_nt = rand_nt()
        else:
            self.new_nt = None


def generate_random_sequence(length):
    """Returns a random DNA sequence of a specified length."""
    s = ""
    for i in range(length):
        s += rand_nt()
    return s


def generate_strain(genome, hv_regions, hvmp, nmp):
    """Creates and then applies mutations to a DNA sequence."""
    mutations_to_add = []
    in_hv = False
    curr_hv = -1
    for c in range(len(genome)):
        if in_hv:
            if c > hv_regions[curr_hv][1]:
                in_hv = False
                curr_hv += 1
        else:
            if c > hv_regions[curr_hv + 1][0]:
                in_hv = True

        mutation_threshold = hvmp if in_hv else nmp
        mp = random.random()
        if mp < mutation_threshold:
            # Add a random mutation
            mutations_to_add.append(Mutation(c))
    sg = genome
    # Indels mess up coordinates, so what we do is keep track of our "shift" --
    # a negative shift when we're at a given coordinate implies more deletions
    # than insertions to the left of that coordinate, and a positive shift
    # similarly implies more insertions than deletions to the left.
    # The reason we can get away with this is that, when creating mutations and
    # when applying mutations, we traverse the genome from left to right (...
    # i.e. 0 to len(genome)) consistently. So when applying a mutation at
    # coordinate C, we only need to worry about the mutations to the left of
    # that mutation.
    shift = 0
    for m in mutations_to_add:
        if m.mtype == "deletion":
            sg = sg[: m.coordinate + shift] + sg[m.coordinate + shift + 1 :]
            shift -= 1
        elif m.mtype == "insertion":
            sg = (
                sg[: m.coordinate + shift]
                + m.new_nt
                + sg[m.coordinate + shift :]
            )
            shift += 1
        elif m.mtype == "mutation":
            sg = (
                sg[: m.coordinate + shift]
                + m.new_nt
                + sg[m.coordinate + shift + 1 :]
            )
        else:
            raise ValueError("Invalid mutation mtype: {}".format(m.mtype))
    return sg


def generate_strains_from_genome(
    genome,
    strain_coverages,
    hypervariable_regions,
    hypervariable_mutation_probability=0.01,
    normal_mutation_probability=0.001,
):
    """
    Parameters
    ----------
        genome: str
            A sequence of DNA.

        strain_coverages: list of int
            A list where each entry is the coverage of a given strain in
            the mock assembly data. The number of entries in this list will be
            equal to the "true" number of strains in the mock data.

        hypervariable_regions: list of (int, int) tuples
            A list of regions in the genome in which mutations occur in strains
            with probability hypervariable_mutation_probability. In other
            regions in the genome, mutations occur with probability
            normal_mutation_probability.

        hypervariable_mutation_probability: float (optional)

        normal_mutation_probability: float (optional)

    Returns
    -------
        strains: list of Strain
            A list of Strain objects "derived" from the genome.
    """
    prev_region_end = 0
    for region in hypervariable_regions:
        if len(region) != 2:
            raise ValueError("HV region with =/= 2 coordinates.")
        for coord in region:
            if coord not in range(0, len(genome)):
                raise ValueError("Invalid coord in HV region.")
        if region[1] <= region[0]:
            raise ValueError("Backwards HV region.")
        if region[0] <= prev_region_end:
            raise ValueError("HV regions out of order and/or overlapping.")
        prev_region_end = region[1]

    out_strains = []
    for cov in strain_coverages:
        strain_seq = generate_strain(
            genome,
            hypervariable_regions,
            hypervariable_mutation_probability,
            normal_mutation_probability,
        )
        out_strains.append(Strain(strain_seq, cov))
    return out_strains
