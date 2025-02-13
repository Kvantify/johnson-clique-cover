"""Tests containing the functionality necessary to reproduce the computational claims
in the paper.
"""

from collections import defaultdict
from itertools import combinations

import gurobipy as gp
import pytest
from gurobipy import GRB

from kvantify.johnsonclique.cover import (
    CliqueType,
    GurobiMaxDisjointCollection,
    GurobiSetCover,
    clique_generator_to_clique,
    clique_to_clique_generator,
)


# Table 1
# Here, we just show how to get a few values in the table.
@pytest.mark.parametrize("N,k,expected", [(7, 3, 9), (8, 4, 14)])
def test_table_1(N: int, k: int, expected: int) -> None:
    cliques = GurobiSetCover(N, k).cover()
    assert len(cliques) == expected


# Example 4.7, part 1
# Let us first test the claim that there must be as many
# cliques of type A as there are of type B.
def test_example_4_7() -> None:
    N = 8
    k = 4
    cover = GurobiSetCover(N, k)
    universe = dict(enumerate(map(frozenset, combinations(range(N), k))))
    possible_cliques = list(cover.maximal_cliques())

    model = gp.Model()
    model.Params.Threads = 16

    x = model.addVars(len(possible_cliques), vtype=GRB.BINARY, name="x")
    for e in universe.values():
        model.addConstr(
            gp.quicksum(
                x[i] for i in range(len(possible_cliques)) if e in possible_cliques[i]
            )
            >= 1,
            name=f"cover_{e}",
        )

    # Optimize on two levels: Slightly favor cliques of type B
    model.setObjective(
        gp.quicksum(
            x[i] * (1.01 if i < len(possible_cliques) // 2 else 1)
            for i in range(len(possible_cliques))
        ),
        GRB.MINIMIZE,
    )

    model.optimize()
    assert model.status == GRB.OPTIMAL
    assert model.objVal == pytest.approx(14 + 0.01 * 7)


# Example 4.7, part 2
# Now, check that if all of [8] is represented by cliques of
# type A, then the minimal cover has more than 14 elements.
def test_example_4_7_part_2() -> None:
    N = 8
    k = 4
    cover = GurobiSetCover(N, k)
    universe = dict(enumerate(map(frozenset, combinations(range(N), k))))
    possible_cliques = list(cover.maximal_cliques())

    model = gp.Model()
    model.Params.Threads = 16

    x = model.addVars(len(possible_cliques), vtype=GRB.BINARY, name="x")
    for e in universe.values():
        model.addConstr(
            gp.quicksum(
                x[i] for i in range(len(possible_cliques)) if e in possible_cliques[i]
            )
            >= 1,
            name=f"cover_{e}",
        )

    cliques_by_elements = defaultdict(list)
    for i, clique in enumerate(possible_cliques):
        # Only take cliques of type A
        gen = clique_to_clique_generator(clique, N, k)
        if gen[0] != CliqueType.A:
            continue
        for e in gen[1]:
            cliques_by_elements[e].append(i)

    # Require that each element is in at least one selected clique
    for e in range(N):
        print(cliques_by_elements[e])
        model.addConstr(gp.quicksum(x[i] for i in cliques_by_elements[e]) >= 1)

    model.setObjective(
        gp.quicksum(x[i] for i in range(len(possible_cliques))), GRB.MINIMIZE
    )

    model.optimize()
    assert model.status == GRB.OPTIMAL
    assert model.objVal == pytest.approx(15)


# Example 4.8, part 1
# Check that we get 66 cliques of each type
@pytest.mark.skip("this takes a while to run")
def test_example_4_8() -> None:
    N = 12
    k = 6
    cover = GurobiSetCover(N, k)
    universe = dict(enumerate(map(frozenset, combinations(range(N), k))))
    possible_cliques = list(cover.maximal_cliques())

    model = gp.Model()
    model.Params.Threads = 16

    x = model.addVars(len(possible_cliques), vtype=GRB.BINARY, name="x")
    for e in universe.values():
        model.addConstr(
            gp.quicksum(
                x[i] for i in range(len(possible_cliques)) if e in possible_cliques[i]
            )
            >= 1,
            name=f"cover_{e}",
        )

    # Optimize on two levels: Slightly favor cliques of type B
    model.setObjective(
        gp.quicksum(
            x[i] * (1.01 if i < len(possible_cliques) // 2 else 1)
            for i in range(len(possible_cliques))
        ),
        GRB.MINIMIZE,
    )

    model.optimize()
    assert model.status == GRB.OPTIMAL
    assert model.objVal == pytest.approx(66 + 0.01 * 66)


# Example 4.8, part 2
# Check that we can not have more than two elements belonging
# to strictly less than 30 cliques.
@pytest.mark.skip("this takes a while to run")
def test_example_4_8_part_2():
    N = 12
    k = 6
    cover = GurobiSetCover(N, k)
    universe = dict(enumerate(map(frozenset, combinations(range(N), k))))
    possible_cliques = list(cover.maximal_cliques())

    model = gp.Model()
    model.Params.Threads = 16

    x = model.addVars(len(possible_cliques), vtype=GRB.BINARY, name="x")
    for e in universe.values():
        model.addConstr(
            gp.quicksum(
                x[i] for i in range(len(possible_cliques)) if e in possible_cliques[i]
            )
            >= 1,
            name=f"cover_{e}",
        )
    zeros = []
    ones = []
    twos = []
    for i, clique in enumerate(possible_cliques):
        gen = clique_to_clique_generator(clique, N, k)
        if gen[0] == CliqueType.A:
            if 0 in gen[1]:
                zeros.append(i)
            if 1 in gen[1]:
                ones.append(i)
            if 2 in gen[1]:
                twos.append(i)
    model.addConstr(gp.quicksum(x[i] for i in zeros) <= 29)
    model.addConstr(gp.quicksum(x[i] for i in ones) <= 29)
    model.addConstr(gp.quicksum(x[i] for i in twos) <= 29)
    model.setObjective(
        gp.quicksum(x[i] for i in range(len(possible_cliques))), GRB.MINIMIZE
    )
    model.optimize()
    assert model.status == GRB.OPTIMAL
    assert model.objVal >= 133


# Example 4.8, part 3
# We check that there is no 132 clique solution in which one
# element belongs to between 1 and 14 cliques.
@pytest.mark.skip("this takes a while to run")
def test_example_4_8_part_3():
    N = 12
    k = 6
    cover = GurobiSetCover(N, k)
    universe = dict(enumerate(map(frozenset, combinations(range(N), k))))
    possible_cliques = list(cover.maximal_cliques())

    model = gp.Model()
    model.Params.Threads = 16

    x = model.addVars(len(possible_cliques), vtype=GRB.BINARY, name="x")
    for e in universe.values():
        model.addConstr(
            gp.quicksum(
                x[i] for i in range(len(possible_cliques)) if e in possible_cliques[i]
            )
            >= 1,
            name=f"cover_{e}",
        )
    zeros = []
    for i, clique in enumerate(possible_cliques):
        gen = clique_to_clique_generator(clique, N, k)
        if gen[0] == CliqueType.A:
            if 0 in gen[1]:
                zeros.append(i)
    model.addConstr(gp.quicksum(x[i] for i in zeros) >= 1)
    model.addConstr(gp.quicksum(x[i] for i in zeros) <= 14)
    model.setObjective(
        gp.quicksum(x[i] for i in range(len(possible_cliques))), GRB.MINIMIZE
    )
    model.optimize()
    assert model.status == GRB.OPTIMAL
    assert model.objVal >= 133


# Proposition 4.10
# We test that we can find at least 371 disjoint cliques
@pytest.mark.skip("The case k = 5 takes a while")
@pytest.mark.parametrize("k,expected", [(1, 1), (2, 2), (3, 4), (4, 14), (5, 36)])
def test_prop_4_10_equality(k, expected) -> None:
    cover = GurobiMaxDisjointCollection(2 * k, k, disp=True)
    res = cover.cover()
    assert len(res) == pytest.approx(expected)


# Proposition 4.10
# We test that we can find at least 371 disjoint cliques for k = 7
@pytest.mark.skip("The takes a while")
def test_prop_4_10() -> None:
    cover = GurobiMaxDisjointCollection(14, 6, 371)
    res = cover.cover()
    assert len(res) >= 371


# Theorem 4.12
# See the "lexicodes" subfolder


# Remark 4.17
# We test that for N = 10, k = 5, there smallest A-B-symmetric
# cover has 48 elements
def test_remark_4_17() -> None:
    N = 10
    k = 5
    cover = GurobiSetCover(N, k)
    universe = dict(enumerate(map(frozenset, combinations(range(N), k))))
    possible_cliques = list(cover.maximal_cliques())

    model = gp.Model()
    model.Params.Threads = 16

    x = model.addVars(len(possible_cliques), vtype=GRB.BINARY, name="x")
    for e in universe.values():
        model.addConstr(
            gp.quicksum(
                x[i] for i in range(len(possible_cliques)) if e in possible_cliques[i]
            )
            >= 1,
            name=f"cover_{e}",
        )

    pairs = []
    for i, clique in enumerate(possible_cliques):
        # Only take cliques of type A
        gen = clique_to_clique_generator(clique, N, k)
        if gen[0] != CliqueType.A:
            continue
        complement = set(range(N)) - gen[1]
        complement_gen = (CliqueType.B, complement)
        complement_clique = clique_generator_to_clique(complement_gen, N)
        complement_clique_index = possible_cliques.index(complement_clique)
        pairs.append((i, complement_clique_index))

    # If a clique is in the cover, its complement must be too
    for i, j in pairs:
        model.addConstr(x[i] == x[j])

    model.setObjective(
        gp.quicksum(x[i] for i in range(len(possible_cliques))), GRB.MINIMIZE
    )

    model.optimize()
    assert model.status == GRB.OPTIMAL
    assert model.objVal == pytest.approx(48)
