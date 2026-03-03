# Johnson graph clique covers

This repository contains various methods for generating small clique covers of Johnson graphs.

## Background

This repository contains the code for the computer-assisted proofs in

* Jørgensen, S. F.. "On the clique covering numbers of Johnson graphs". *Des. Codes Cryptogr.* **93**, 3689–3705 (2025). [doi:10.1007/s10623-025-01663-3](https://doi.org/10.1007/s10623-025-01663-3). [arXiv:2502.15019](https://arxiv.org/abs/2502.15019).

Concretely, [`tests/test_reproduce_paper.py`](https://github.com/Kvantify/johnson-clique-cover/blob/master/tests/test_reproduce_paper.py) contains a collection of tests for each of the claims made in the above paper.

## Usage

The package contains a handful of different solvers with different performance/solution quality trade-offs.

For example, to find clique covers of the Johnson graph $J(9, 4)$, examples of different results are obtained as follows:

```python
from kvantify.johnsonclique import (
    TrivialCover,
    AllMaximalCliquesCover,
    GreedySetCover,
    RecursiveCover,
    OptimalSetCover,
)

N, k = 9, 4
cover = TrivialCover(N, k).cover()
print(len(cover))  # 126
cover = AllMaximalCliquesCover(N, k).cover()
print(len(cover))  # 210
cover = GreedySetCover(N, k).cover()
print(len(cover))  # 28
cover = RecursiveCover(N, k).cover()
print(len(cover))  # 35
cover = OptimalSetCover(N, k).cover()
print(len(cover))  # 25
```
