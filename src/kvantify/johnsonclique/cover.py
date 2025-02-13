from abc import ABC, abstractmethod
from enum import Enum
from itertools import combinations
from math import comb
from typing import Iterator

import gurobipy as gp
import numpy as np
from gurobipy import GRB
from scipy.optimize import linprog

Vertex = frozenset[int]
Clique = frozenset[Vertex]
CliqueCover = set[Clique]


class CliqueCoverSolver(ABC):
    def __init__(self, N: int, k: int):
        if N < 2:
            raise ValueError("N must be at least 2")
        if k < 1 or k > N - 1:
            raise ValueError("k must satisfy 0 < k < N")
        self.N = N
        if N < 2 * k:
            self.k = N - k
            self.flip = True
        else:
            self.k = k
            self.flip = False

    def maximal_cliques(
        self, A_only: bool = False, B_only: bool = False
    ) -> Iterator[Clique]:
        N, k = self.N, self.k
        if N == 2:
            yield frozenset({frozenset({0}), frozenset({1})})
            return
        if not B_only and k < N - 1:
            for overlap in map(set, combinations(range(N), k - 1)):
                yield frozenset(
                    frozenset(overlap | {x}) for x in set(range(N)) - overlap
                )
        if not A_only and k > 1:
            for overlap in map(set, combinations(range(N), k + 1)):
                yield frozenset(frozenset(overlap - {x}) for x in overlap)

    def flip_if_necessary(self, cover: CliqueCover) -> CliqueCover:
        if not self.flip:
            return cover
        flipped_cover = set()
        for clique in cover:
            flipped_clique: Clique = frozenset(
                frozenset(set(range(self.N)) - vertex) for vertex in clique
            )
            flipped_cover.add(flipped_clique)
        return flipped_cover

    @abstractmethod
    def cover(self) -> CliqueCover: ...


class TrivialCover(CliqueCoverSolver):
    def __init__(self, N: int, k: int):
        super().__init__(N, k)

    def cover(self) -> CliqueCover:
        N, k = self.N, self.k
        return self.flip_if_necessary(
            {frozenset([frozenset(s)]) for s in combinations(range(N), k)}
        )


class AllMaximalCliquesCover(CliqueCoverSolver):
    def __init__(self, N: int, k: int):
        super().__init__(N, k)

    def cover(self) -> CliqueCover:
        return self.flip_if_necessary(frozenset(self.maximal_cliques()))


class GreedySetCover(CliqueCoverSolver):
    def __init__(self, N: int, k: int):
        super().__init__(N, k)

    def cover(self) -> CliqueCover:
        N, k = self.N, self.k
        covered = set()
        cliques = set()
        while len(covered) != comb(N, k):
            # Pick the maximal clique leading to the greatest increase in
            # covered vertices
            clique = max(self.maximal_cliques(), key=lambda c: len(covered | c))
            covered |= clique
            cliques.add(clique)
        return self.flip_if_necessary(cliques)


class CliqueType(Enum):
    A = 0
    B = 1


CliqueGenerator = tuple[CliqueType, frozenset[int]]

CACHE: dict[tuple[int, int], list[CliqueGenerator]] = {
    (2, 1): [(CliqueType.A, set())],
    (3, 1): [(CliqueType.A, set())],
    (3, 2): [(CliqueType.B, {0, 1, 2})],
}


def clique_generator_to_clique(generator: CliqueGenerator, N: int) -> Clique:
    clique_type, base_set = generator
    return (
        frozenset(frozenset(base_set | {i}) for i in set(range(N)) - base_set)
        if clique_type == CliqueType.A
        else frozenset(frozenset(base_set - {i}) for i in base_set)
    )


def clique_generators(N: int, k: int) -> list[CliqueGenerator]:
    if generators := CACHE.get((N, k), []):
        return generators
    if k > 1:
        for generator in clique_generators(N - 1, k - 1):
            clique_type, base_set = generator
            generators.append((clique_type, base_set | {N - 1}))
    if k < N:
        generators += clique_generators(N - 1, k)
    CACHE[N, k] = generators
    return generators


def clique_to_clique_generator(clique: Clique, N: int, k: int) -> CliqueGenerator:
    intersection = frozenset.intersection(*clique)
    if intersection:
        return (CliqueType.A, intersection)
    else:
        union = frozenset.union(*clique)
        return (CliqueType.B, union)


class RecursiveCover(CliqueCoverSolver):
    def __init__(self, N: int, k: int):
        super().__init__(N, k)

    def cover(self) -> CliqueCover:
        N, k = self.N, self.k
        generators = clique_generators(N, k)
        cliques = frozenset(clique_generator_to_clique(g, N) for g in generators)
        return self.flip_if_necessary(cliques)


class GurobiSetCover(CliqueCoverSolver):
    def __init__(self, N: int, k: int, A_only: bool = False, disp: bool = False):
        super().__init__(N, k)
        self.A_only = A_only
        self.disp = disp

    def cover(self) -> CliqueCover:
        N, k = self.N, self.k
        universe = dict(enumerate(map(frozenset, combinations(range(N), k))))
        possible_cliques = list(self.maximal_cliques(self.A_only))

        model = gp.Model()
        if not self.disp:
            model.setParam("OutputFlag", 0)
        model.Params.Threads = 16

        x = model.addVars(len(possible_cliques), vtype=GRB.BINARY, name="x")
        for e in universe.values():
            model.addConstr(
                gp.quicksum(
                    x[i]
                    for i in range(len(possible_cliques))
                    if e in possible_cliques[i]
                )
                >= 1,
                name=f"cover_{e}",
            )

        model.setObjective(
            gp.quicksum(x[i] for i in range(len(possible_cliques))), GRB.MINIMIZE
        )

        # Optimize the model
        model.optimize()
        if model.status == GRB.OPTIMAL:
            selected_sets = [i for i in range(len(possible_cliques)) if x[i].X > 0.5]
        cliques = [possible_cliques[i] for i in selected_sets]
        gens = [clique_to_clique_generator(c, N, k) for c in cliques]
        cliques = frozenset(clique_generator_to_clique(g, N) for g in gens)
        return cliques


class GurobiMaxDisjointCollection(CliqueCoverSolver):
    def __init__(
        self,
        N: int,
        k: int,
        A_only: bool = False,
        disp: bool = False,
        best_obj_stop: int | None = None,
    ):
        super().__init__(N, k)
        self.A_only = A_only
        self.disp = disp
        self.best_obj_stop = best_obj_stop

    def cover(self) -> CliqueCover:
        N, k = self.N, self.k
        possible_cliques = list(self.maximal_cliques(self.A_only))
        possible_cliques = [
            c for c in possible_cliques if len(c) == max(k + 1, N - k + 1)
        ]

        model = gp.Model()
        if not self.disp:
            model.setParam("OutputFlag", 0)
        if self.best_obj_stop is not None:
            model.setParam("BestObjStop", self.best_obj_stop)
        model.Params.Threads = 16

        x = model.addVars(len(possible_cliques), vtype=GRB.BINARY, name="x")

        # Only disjoint sets
        for (i1, c1), (i2, c2) in combinations(enumerate(possible_cliques), 2):
            if not c1.isdisjoint(c2):
                model.addConstr(x[i1] + x[i2] <= 1)

        model.setObjective(
            gp.quicksum(x[i] for i in range(len(possible_cliques))), GRB.MAXIMIZE
        )
        model.optimize()
        if model.status == GRB.OPTIMAL:
            selected_sets = [i for i in range(len(possible_cliques)) if x[i].X > 0.5]
        cliques = [possible_cliques[i] for i in selected_sets]
        gens = [clique_to_clique_generator(c, N, k) for c in cliques]
        cliques = frozenset(clique_generator_to_clique(g, N) for g in gens)
        return cliques


def clique_generator_str(gen, N):
    s = "$A" if gen[0] == CliqueType.A else "B"
    s += "^{" + str(N) + "}" + r"_{\{" + str(set(gen[1])) + r"\}}"
    s += "$"
    return s


class ScipyMaxDisjointCollection(CliqueCoverSolver):
    def __init__(self, N: int, k: int, A_only: bool = False, disp: bool = False):
        super().__init__(N, k)
        self.A_only = A_only
        self.disp = disp

    def cover(self) -> CliqueCover:
        N, k = self.N, self.k
        possible_cliques = list(self.maximal_cliques(self.A_only))
        possible_cliques = [
            c for c in possible_cliques if len(c) == max(k + 1, N - k + 1)
        ]
        c = np.repeat(1.0, len(possible_cliques))
        overlap = []
        for (i1, c1), (i2, c2) in combinations(enumerate(possible_cliques), 2):
            if not c1.isdisjoint(c2):
                overlap.append((i1, i2))

        A = np.zeros((len(overlap), len(possible_cliques)))
        for j, (i1, i2) in enumerate(overlap):
            A[j, i1] = 1
            A[j, i2] = 1
        res = linprog(
            c=-c,
            A_ub=A,
            b_ub=np.ones(A.shape[0]),
            bounds=(0, 1),
            integrality=1,
            options={"disp": self.disp},
        )
        chosen_clique_indices = np.where(res.x > 0.5)[0]
        return self.flip_if_necessary(
            frozenset(possible_cliques[ci] for ci in chosen_clique_indices)
        )


class OptimalSetCover(CliqueCoverSolver):
    def __init__(self, N: int, k: int, A_only: bool = False, disp: bool = False):
        super().__init__(N, k)
        self.A_only = A_only
        self.disp = disp

    def cover(self) -> CliqueCover:
        N, k = self.N, self.k
        universe = dict(enumerate(map(frozenset, combinations(range(N), k))))
        universe_rev = {v: k for k, v in universe.items()}
        possible_cliques = list(self.maximal_cliques(self.A_only))
        c = np.repeat(1.0, len(possible_cliques))
        c[: len(c) // 2] += 0.00001
        # print(c)
        A = np.zeros((len(universe), len(possible_cliques)))
        for si, s in enumerate(possible_cliques):
            for x in s:
                A[universe_rev[x], si] = 1
        res = linprog(
            c,
            A_ub=-A,
            b_ub=-np.ones(A.shape[0]),
            bounds=(0, 1),
            integrality=1,
            options={"disp": self.disp},
        )
        chosen_clique_indices = np.where(res.x > 0.5)[0]
        return self.flip_if_necessary(
            frozenset(possible_cliques[ci] for ci in chosen_clique_indices)
        )

    def write_lp(self):
        N, k = self.N, self.k
        lp_str = ["MINIMIZE"]
        universe = dict(enumerate(map(frozenset, combinations(range(N), k))))
        universe_rev = {v: k for k, v in universe.items()}
        possible_cliques = list(self.maximal_cliques(self.A_only))
        A = np.zeros((len(universe), len(possible_cliques)))
        for si, s in enumerate(possible_cliques):
            for x in s:
                A[universe_rev[x], si] = 1
        lp_str += [" + ".join(f"s{i}" for i in range(len(possible_cliques)))]
        lp_str += ["SUBJECT TO"]

        for row in range(A.shape[0]):
            lp_str += [
                " + ".join(f"s{col}" for col in range(A.shape[1]) if A[row, col] != 0)
                + " >= 1"
            ]
        lp_str += ["BINARY"]
        lp_str += [" ".join(f"s{i}" for i in range(len(possible_cliques)))]
        lp_str += ["END"]
        with open("this_lp.lp", "w") as f:
            f.write("\n".join(lp_str))
