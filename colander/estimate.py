import networkx as nx
import copy


class GreedyError(Exception):
    """Used for errors in the greedy estimation algorithm."""

    pass


class Cycle:
    """Class representing a weighted Cycle in an arbitrary graph."""

    def __init__(self, edges, weight):
        self.edges = edges
        self.weight = weight

    def __str__(self):
        return "Cycle with wgt {} and edges {}".format(self.weight, self.edges)

    def __repr__(self):
        return "Cycle({}, {})".format(self.edges, self.weight)


class CycleSet:
    """Class representing a set of weighted cycles.

       ...Actually uses a list instead of a set internally. Sorry for ruining
       biology.
    """

    def __init__(self):
        self.cycles = []

    def __len__(self):
        return len(self.cycles)

    def __repr__(self):
        output = "CycleSet of {"
        for c in self.cycles:
            output += "\n\t" + str(c)
        output += "\n}"
        return output

    def add(self, cycle):
        self.cycles.append(cycle)

    def cycle_coverage(self, e):
        """Computes the CycleCoverage of an edge e relative to this CycleSet.

           This is defined as the sum of the weights of all cycles in the
           CycleSet that include e.
        """
        cc = 0
        for cycle in self.cycles:
            if e in cycle.edges:
                cc += cycle.weight
        return cc

    def conformity_score(self, G):
        """Computes the conformity score of this CycleSet for a graph G.

           Defined as the sum for each edge of the squared difference between
           the edge's coverage and the cycle coverage of this edge.

           Note that this makes no assumptions that all edges in G are in this
           CycleSet -- it's perfectly possible for an edge in G to be absent
           from this set, in which case a penalty on on the conformity score of
           (that edge's coverage) squared is incurred.
        """
        score = 0
        print("doin it")
        for e in G.edges:
            penalty = (get_cov(G, e) - self.cycle_coverage(e)) ** 2
            print("edge {} causes penalty {}".format(e, penalty))
            score += penalty
        return score


def get_cov(G, e):
    """Returns the cov attribute of an edge e in a graph G."""
    return G.get_edge_data(*e)["cov"]


def get_max_weight_edge_from_node(G, n, seen_edges):
    maxweightoutedge = None
    maxweight = 0
    for e in G.out_edges(n):
        if e not in seen_edges:
            ecov = get_cov(G, e)
            if ecov > maxweight:
                maxweight = ecov
                maxweightoutedge = e

    if maxweightoutedge is None:
        raise GreedyError("All out edges have cov <= 0.")

    return maxweightoutedge


def peel_max_weight_cycle(G):
    # start at a 1-in 1-out node to make sure the algorithm doesn't get "lost"
    starting_node = None
    for n in G.nodes:
        if G.in_degree(n) == G.out_degree(n) == 1:
            starting_node = n
            break
    if starting_node is None:
        raise ValueError(
            "Graph has no one-in-one-out nodes. This algorithm isn't that "
            "robust yet..."
        )

    prev_edge = get_max_weight_edge_from_node(G, starting_node, [])

    cycle_min_weight = get_cov(G, prev_edge)
    curr_node = prev_edge[1]
    cycle_edges = [prev_edge]
    while starting_node != curr_node:
        curr_node = prev_edge[1]
        # pick max weight outgoing edge that we haven't seen before in
        # cycle_edges and follow.
        next_edge = get_max_weight_edge_from_node(G, curr_node, cycle_edges)
        # Update the cycle
        cycle_edges.append(next_edge)
        # Decrease the "capacity" of this cycle if necessary
        cycle_min_weight = min(cycle_min_weight, get_cov(G, next_edge))

    return Cycle(cycle_edges, cycle_min_weight)


def remove_cycle_impact(G, cycle):
    covdict = {}
    for e in cycle.edges:
        covdict[e] = {"cov": get_cov(G, e) - cycle.weight}
    nx.set_edge_attributes(G, covdict)


def greedy_strain_estimation(G, N):
    """Tries to produce a good weighted cycle-set defined for G and N.

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
    Gc = copy.deepcopy(G)

    cycles = CycleSet()
    while len(cycles) < N:
        try:
            cycle = peel_max_weight_cycle(Gc)
        except GreedyError:
            # We seem to have exhausted the "capacity" of the graph
            break
        CycleSet.add(cycle)
        # Here's the critical thing: modify the graph to remove the "impact" of
        # the selected cycle.
        remove_cycle_impact(Gc, cycle)
    return cycles
