# Johnson graph clique covers

This repository contains various methods for generating small clique covers of Johnson graphs.

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