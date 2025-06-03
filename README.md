# Johnson graph clique covers

This repository contains various methods for generating small clique covers of Johnson graphs.

## Background

This repository contains the code for the computer-assisted proofs in

* S. F. Jørgensen. "On the clique covering numbers of Johnson graphs". *Designs, Codes and Cryptography* (to appear). [arXiv:2502.15019](https://arxiv.org/abs/2502.15019).

Concretely, [`tests/test_reproduce_paper.py`](https://github.com/Kvantify/johnson-clique-cover/blob/master/tests/test_reproduce_paper.py) contains a collection of tests for each of the claims made in the above paper.

### Follow-up work

The list below contains a collection of sources that build upon the above work.

* C. Rosin. "Using Reasoning Models to Generate Search Heuristics that Solve Open Instances of Combinatorial Design Problems". [arxiv:2505.23881](https://arxiv.org/abs/2505.23881).

This work describes a framework for automatically generating efficient implementations of solvers that in some cases improve on the bounds in the original paper. The implementations can be found [here](https://github.com/Constructive-Codes/CPro1/tree/main/designs/johnson-clique-cover).

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
