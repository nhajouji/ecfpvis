# Elliptic Curve Visualization

A Python library for computing with and visualizing elliptic curves, with a focus on their arithmetic over finite fields and connections to the theory of complex multiplication.

## Features

- **Elliptic curves over finite fields** — trace of Frobenius, point counting, isogeny classes
- **Mordell–Weil group structure** — generators and lattice representations
- **Complex multiplication** — quadratic forms, Hilbert and Atkin modular polynomials
- **Isogeny graphs** — visualizations of X₀(ℓ) modular curves and isogeny walks
- **Number theory utilities** — GCD, quadratic reciprocity, prime generation, ring-of-integer computations
- **Algebraic structures** — polynomials over ℤ and 𝔽ₚ, integer square matrices

## Project Structure

```
.
├── ecv.py                  # Main visualization module (matplotlib-based)
├── ecfp.py                 # Elliptic curves over finite fields (Fp)
├── qfs.py                  # Quadratic forms and class groups
├── nt.py                   # Number theory utilities
├── ringclasses.py          # Algebraic structures (polynomials, matrices)
├── modularpolynomials.py   # Hilbert and Atkin modular polynomial utilities
├── data/
│   ├── hilbpolys.json      # Precomputed Hilbert class polynomials
│   └── atkinpolys.json     # Precomputed Atkin modular polynomials
└── requirements.txt
```

## Requirements

- Python 3.10+
- [matplotlib](https://matplotlib.org/)
- [numpy](https://numpy.org/)

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/nhajouji/elliptic-curve-visualization.git
   cd elliptic-curve-visualization
   ```

2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. Make sure the `data/` folder is present in the project root (it contains precomputed polynomial data required by `modularpolynomials.py`).

## Usage

Import the visualization module and start exploring:

```python
from ecv import make_both_pics, IsogenyClassFp
import matplotlib.pyplot as plt

# Visualize an elliptic curve and its associated lattice
# fg = (f, g) coefficients of y² = x³ + fx + g
# ap = (a, p) trace of Frobenius and prime
# abc = quadratic form (a, b, c) for the CM lattice
make_both_pics(fg=(1, 0), ap=(1, 5), abc=(1, 0, 1))
plt.show()
```

You can also work directly with elliptic curves over finite fields:

```python
from ecfp import mw_gens, get_j_to_qfs_dict

# Get Mordell-Weil generators for a curve over Fp
gens = mw_gens(fg=(1, 0), p=17)

# Get a mapping from j-invariants to quadratic forms
j_to_qfs = get_j_to_qfs_dict(p=17)
```

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.
