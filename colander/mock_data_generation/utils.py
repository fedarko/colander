import random

nts = set(["A", "C", "T", "G"])


def rand_nt(exclude=None):
    if exclude is not None:
        nts_to_consider = list(nts.difference(exclude))
    else:
        nts_to_consider = list(nts)
    return random.choice(nts_to_consider)


class Strain:
    def __init__(self, seq, coverage, originating_mutations):
        self.seq = seq
        self.coverage = coverage
        self.originating_mutations = originating_mutations


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

    def make_substitution(self, new_nt=None):
        self.mtype = "substitution"
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
            self.make_substitution()

    def __repr__(self):
        if self.mtype == "insertion":
            return "insertion at pos {}: {} -> {}{}".format(
                self.coordinate, self.curr_nt, self.new_nt, self.curr_nt
            )
        elif self.mtype == "deletion":
            return "deletion at pos {}: {} -> ''".format(
                self.coordinate, self.curr_nt
            )
        elif self.mtype == "substitution":
            return "substitution at pos {}: {} -> {}".format(
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
    for m in mutations:
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
        elif m.mtype == "substitution":
            seq2 = (
                seq2[: m.coordinate + shift]
                + m.new_nt
                + seq2[m.coordinate + shift + 1 :]
            )
        else:
            raise ValueError("Invalid mutation mtype: {}".format(m.mtype))
    return seq2


def validate_genomic_regions(genome, regions):
    """Validates a list of "regions" (coordinates along a genome string)."""

    prev_region_end = -1
    for region in regions:
        if len(region) != 2:
            raise ValueError("Region with =/= 2 coordinates.")
        for coord in region:
            if coord not in range(0, len(genome)):
                raise ValueError("Out-of-range genomic coordinate in region.")
        if region[1] <= region[0]:
            raise ValueError("Backwards region.")
        if region[0] <= prev_region_end:
            raise ValueError("Regions out of order and/or overlapping.")
        prev_region_end = region[1]


def generate_mutations(genome, hv_regions, hvmp, nmp, only_subs=False):
    """Determines where to place mutations on a DNA sequence.

       If only_subs is True, this will only create substitution mutations.
       Otherwise, the type of mutation created (insertion / deletion /
       substitution) is random.
    """
    validate_genomic_regions(genome, hv_regions)
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
                if c >= hv_regions[curr_hv + 1][0]:
                    in_hv = True

        mutation_threshold = hvmp if in_hv else nmp
        mp = random.random()
        if mp < mutation_threshold:
            # Add a random mutation
            m = Mutation(c, genome[c])
            if only_subs:
                m.make_substitution()
            else:
                m.randomly_initialize()
            mutations_to_add.append(m)
    return mutations_to_add


def generate_strain(genome, coverage, hv_regions, hvmp, nmp):
    """Creates and then applies mutations to a DNA sequence.

       Produces a Strain object as output.
    """
    validate_genomic_regions(genome, hv_regions)
    mutations = generate_mutations(genome, hv_regions, hvmp, nmp)
    strain_seq = add_mutations(genome, mutations)
    return Strain(strain_seq, coverage, mutations)


def generate_strains_from_genome(
    genome,
    strain_coverages,
    hypervariable_regions,
    hypervariable_mutation_probability=0.01,
    normal_mutation_probability=0.001,
):
    """From a sequence of DNA ("genome"), generates a list of Strain objects.

    These strains' sequences are generated by applying random mutations
    (insertions, deletions, and substitutions) at each position along the
    initial genome. The probability of applying a mutation at any given
    position is dependent on the parameters to this function -- see the
    descriptions below.

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
            normal_mutation_probability. (The coordinates defining the
            boundaries of these regions are 0-indexed.)

            Note that you can totally pass in [] if you just want all parts of
            the genome to have the same mutation probability! That's fine.

        hypervariable_mutation_probability: float (optional)
            Probability of mutation at a position in a hypervariable region.

        normal_mutation_probability: float (optional)
            Probability of mutation at a position outside of a hypervariable
            region.

    Returns
    -------
        strains: list of Strain
            A list of Strain objects "derived" from the genome.
    """
    validate_genomic_regions(genome, hypervariable_regions)
    out_strains = []
    for cov in strain_coverages:
        out_strains.append(
            generate_strain(
                genome,
                cov,
                hypervariable_regions,
                hypervariable_mutation_probability,
                normal_mutation_probability,
            )
        )
    return out_strains


def shear_into_reads(seq, coverage, read_length, error_probability=0):
    """Breaks a sequence into reads.

       This function is fairly unrealistic, compared to actual shotgun
       sequencing. (e.g. reads are not always created from the same places in
       the genome...) However, for the purposes of simulation for a class
       project, this should be ok :)

    Parameters
    ----------
        seq: str

        coverage: int

        read_length: int

        error_probability: float
            Probability of an error (defined as just a substitution mutation).
            If you make this 0, the reads will be "error-free."
    Returns
    -------
        reads: list of str
    """
    reads = []
    for c in range(coverage):
        # Shear off read_length characters from the start of the sequence
        # repeatedly. There's probably a more efficient way to do this.
        # e.g. AAACCCGGGTTT, read_length = 5:
        #      AAACC
        #           CGGGT
        #                TT
        # When seq_copy gets to be less than read_length characters, this
        # will still work as intended -- it'll just add a read of length
        # len(seq) % read_length to reads.
        seq_copy = seq
        while len(seq_copy) > 0:
            reads.append(seq_copy[0:read_length])
            seq_copy = seq_copy[read_length:]

    # apply errors to reads
    for i in range(len(reads)):
        errored_read = add_mutations(
            reads[i],
            generate_mutations(
                reads[i], [], 0, error_probability, only_subs=True
            ),
        )
        reads[i] = errored_read
    return reads
