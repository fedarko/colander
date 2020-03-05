class Cycle:
    def __init__(self, edges):
        self.edges = edges


class CycleSet:
    def __init__(self):
        self.cycle2weight = {}

    def cycle_coverage(self, e):
        """Computes the CycleCoverage of an edge e relative to this CycleSet.

           This is defined as the sum of the weights of all cycles in the
           CycleSet that include e.
        """
        cc = 0
        for cycle, weight in self.cycle2weight.items():
            if e in cycle.edges:
                cc += weight
        return cc

    def conformity_score(self, G):
        """Computes the conformity score of this CycleSet for a graph G."""
        score = 0
        for e in G.edges:
            score += (e["cov"] - self.cycle_coverage(e)) ** 2
        return score


def estimate_strains(G, N):
    """Produces an optimal weighted cycle-set defined for G and N.

    Parameters
    ----------
        G: networkx.DiGraph
            A directed, strongly-connected edge-weighted graph.

        N: int
            A positive integer.

    Returns
    -------
        S: CycleSet object
            An (ideally optimal) set of N or fewer cycles in the graph G.
    """
