from abc import ABC, abstractmethod
from enum import Enum
from itertools import combinations
from math import comb
from typing import Iterator

import numpy as np
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

    def maximal_cliques(self) -> Iterator[Clique]:
        N, k = self.N, self.k
        if N == 2:
            yield frozenset({frozenset({0}), frozenset({1})})
            return
        if k < N - 1:
            for overlap in map(set, combinations(range(N), k - 1)):
                yield frozenset(
                    frozenset(overlap | {x}) for x in set(range(N)) - overlap
                )
        if k > 1:
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


class RecursiveCover(CliqueCoverSolver):
    def __init__(self, N: int, k: int):
        super().__init__(N, k)

    def cover(self) -> CliqueCover:
        N, k = self.N, self.k
        generators = clique_generators(N, k)
        cliques = frozenset(clique_generator_to_clique(g, N) for g in generators)
        return self.flip_if_necessary(cliques)


class OptimalSetCover(CliqueCoverSolver):
    def __init__(self, N: int, k: int):
        super().__init__(N, k)

    def cover(self) -> CliqueCover:
        N, k = self.N, self.k
        universe = dict(enumerate(map(frozenset, combinations(range(N), k))))
        universe_rev = {v: k for k, v in universe.items()}
        possible_cliques = list(self.maximal_cliques())
        c = np.repeat(1, len(possible_cliques))
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
        )
        chosen_clique_indices = np.where(res.x > 0.5)[0]
        return self.flip_if_necessary(
            frozenset(possible_cliques[ci] for ci in chosen_clique_indices)
        )
