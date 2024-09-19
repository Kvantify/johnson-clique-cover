from itertools import combinations, product

import pytest

from kvantify.johnsonclique import (
    AllMaximalCliquesCover,
    CliqueCover,
    GreedySetCover,
    OptimalSetCover,
    RecursiveCover,
    TrivialCover,
)


def is_cover(cover: CliqueCover, N: int, k: int) -> bool:
    covered = frozenset.union(*cover)
    return covered == set(map(frozenset, combinations(range(N), k)))


def all_cliques(cover: CliqueCover, N: int, k: int) -> bool:
    for subset in cover:
        for vertex in subset:
            if len(vertex) != k:
                return False
            if any(x < 0 or x > N - 1 for x in vertex):
                return False
        for v1, v2 in product(subset, subset):
            if len(v1 & v2) < k - 1:
                return False
    return True


all_solvers = [
    AllMaximalCliquesCover,
    GreedySetCover,
    OptimalSetCover,
    RecursiveCover,
    TrivialCover,
]


@pytest.mark.parametrize("solver", all_solvers)
def test_bad_inputs(solver):
    with pytest.raises(ValueError):
        solver(1, 1)
    with pytest.raises(ValueError):
        solver(0, 1)
    with pytest.raises(ValueError):
        solver(2, 0)
    with pytest.raises(ValueError):
        solver(2, 2)
    with pytest.raises(ValueError):
        solver(3, 5)


@pytest.mark.parametrize("solver", all_solvers)
def test_edge_case_k(solver):
    N, k = 4, 2
    cover = solver(N, k).cover()
    assert is_cover(cover, N, k)
    assert all_cliques(cover, N, k)
    N, k = 4, 3
    cover = solver(N, k).cover()
    assert is_cover(cover, N, k)
    assert all_cliques(cover, N, k)


@pytest.mark.parametrize(
    "solver", [GreedySetCover, AllMaximalCliquesCover, OptimalSetCover, RecursiveCover]
)
def test_many_2_1(solver):
    N, k = 2, 1
    cover = solver(N, k).cover()
    assert len(cover) == 1
    assert cover == {frozenset({frozenset({0}), frozenset({1})})}
    assert is_cover(cover, N, k)
    assert all_cliques(cover, N, k)


@pytest.mark.parametrize(
    "solver", [GreedySetCover, AllMaximalCliquesCover, OptimalSetCover, RecursiveCover]
)
def test_many_3_1(solver):
    N, k = 3, 1
    cover = solver(N, k).cover()
    assert len(cover) == 1
    assert cover == {frozenset({frozenset({0}), frozenset({1}), frozenset({2})})}
    assert is_cover(cover, N, k)
    assert all_cliques(cover, N, k)


@pytest.mark.parametrize(
    "solver", [GreedySetCover, AllMaximalCliquesCover, OptimalSetCover, RecursiveCover]
)
def test_many_3_2(solver):
    N, k = 3, 2
    cover = solver(N, k).cover()
    assert len(cover) == 1
    assert cover == {
        frozenset({frozenset({0, 1}), frozenset({0, 2}), frozenset({1, 2})})
    }
    assert is_cover(cover, N, k)
    assert all_cliques(cover, N, k)


@pytest.mark.parametrize("solver", [GreedySetCover, OptimalSetCover, RecursiveCover])
def test_many_4_2(solver):
    N, k = 4, 2
    cover = solver(N, k).cover()
    assert len(cover) == 2
    assert is_cover(cover, N, k)
    assert all_cliques(cover, N, k)


@pytest.mark.parametrize("solver", [GreedySetCover, OptimalSetCover, RecursiveCover])
def test_many_5_3(solver):
    N, k = 5, 3
    cover = solver(N, k).cover()
    assert len(cover) == 3
    assert is_cover(cover, N, k)
    assert all_cliques(cover, N, k)


def test_discrepancy_7_3():
    N, k = 7, 3

    cover = TrivialCover(N, k).cover()
    assert len(cover) == 35
    assert is_cover(cover, N, k)
    assert all_cliques(cover, N, k)

    cover = AllMaximalCliquesCover(N, k).cover()
    assert len(cover) == 56
    assert is_cover(cover, N, k)
    assert all_cliques(cover, N, k)

    cover = GreedySetCover(N, k).cover()
    assert len(cover) == 9
    assert is_cover(cover, N, k)
    assert all_cliques(cover, N, k)

    cover = RecursiveCover(N, k).cover()
    assert len(cover) == 10
    assert is_cover(cover, N, k)
    assert all_cliques(cover, N, k)

    cover = OptimalSetCover(N, k).cover()
    assert len(cover) == 9
    assert is_cover(cover, N, k)
    assert all_cliques(cover, N, k)


def test_recursive_outperforms_greedy_6_3():
    N, k = 6, 3
    cover = GreedySetCover(N, k).cover()
    assert len(cover) == 7
    assert is_cover(cover, N, k)
    assert all_cliques(cover, N, k)

    cover = RecursiveCover(N, k).cover()
    assert len(cover) == 6
    assert is_cover(cover, N, k)
    assert all_cliques(cover, N, k)

    cover = OptimalSetCover(N, k).cover()
    assert len(cover) == 6
    assert is_cover(cover, N, k)
    assert all_cliques(cover, N, k)
