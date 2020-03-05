import random

nts = set(["A", "C", "T", "G"])


def rand_nt(exclude=None):
    if exclude is not None:
        nts_to_consider = list(nts.difference(exclude))
    else:
        nts_to_consider = list(nts)
    return random.choice(nts_to_consider)


class Strain:
    def __init__(self, seq, coverage):
        self.seq = seq
        self.coverage = coverage


class Mutation:
    def __init__(self, coordinate, curr_nt):
        self.coordinate = coordinate
        self.curr_nt = curr_nt
        self.mtype = "indeterminate"
        self.new_nt = None

    def make_insertion(self, new_nt=None):
        self.mtype = "insertion"
        if new_nt is None:
            self.new_nt = rand_nt()
        else:
            self.new_nt = new_nt

    def make_deletion(self):
        self.mtype = "deletion"
        self.new_nt = None

    def make_mutation(self, new_nt=None):
        self.mtype = "mutation"
        if new_nt is None:
            self.new_nt = rand_nt(exclude=self.curr_nt)
        else:
            self.new_nt = new_nt

    def randomly_initialize(self):
        # Randomly determine mutation "type" then adjust
        r = random.random()
        if r < 1 / 3:
            self.make_insertion()
        elif r < 2 / 3:
            self.make_deletion()
        else:
            self.make_mutation()

    def __repr__(self):
        if self.mtype == "insertion":
            return "insertion at pos {}: {} -> {}{}".format(
                self.coordinate, self.curr_nt, self.new_nt, self.curr_nt
            )
        elif self.mtype == "deletion":
            return "deletion at pos {}: {} -> ''".format(
                self.coordinate, self.curr_nt
            )
        elif self.mtype == "mutation":
            return "mutation at pos {}: {} -> {}".format(
                self.coordinate, self.curr_nt, self.new_nt
            )
        else:
            return "indeterminate mutation"


def generate_random_sequence(length):
    """Returns a random DNA sequence of a specified length."""
    s = ""
    for i in range(length):
        s += rand_nt()
    return s


def add_mutations(seq, mutations):
    """Returns a modified DNA sequence based on a list of Mutation objects."""
    # Indels mess up coordinates, so what we do is keep track of our "shift" --
    # a negative shift when we're at a given coordinate implies more deletions
    # than insertions to the left of that coordinate, and a positive shift
    # similarly implies more insertions than deletions to the left.
    # The reason we can get away with this is that, when creating mutations and
    # when applying mutations, we traverse the genome from left to right (...
    # i.e. 0 to len(genome)) consistently. So when applying a mutation at
    # coordinate C, we only need to worry about the mutations to the left of
    # that mutation.
    seq2 = seq
    shift = 0
    print(seq2)
    for m in mutations:
        print(m)
        if m.mtype == "deletion":
            seq2 = (
                seq2[: m.coordinate + shift] + seq2[m.coordinate + shift + 1 :]
            )
            shift -= 1
        elif m.mtype == "insertion":
            seq2 = (
                seq2[: m.coordinate + shift]
                + m.new_nt
                + seq2[m.coordinate + shift :]
            )
            shift += 1
        elif m.mtype == "mutation":
            seq2 = (
                seq2[: m.coordinate + shift]
                + m.new_nt
                + seq2[m.coordinate + shift + 1 :]
            )
        else:
            raise ValueError("Invalid mutation mtype: {}".format(m.mtype))
        print(seq2)
    return seq2


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
            if curr_hv + 1 < len(hv_regions):
                if c > hv_regions[curr_hv + 1][0]:
                    in_hv = True

        mutation_threshold = hvmp if in_hv else nmp
        mp = random.random()
        if mp < mutation_threshold:
            # Add a random mutation
            m = Mutation(c, genome[c])
            m.randomly_initialize()
            mutations_to_add.append(m)
    return add_mutations(genome, mutations_to_add)


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
