# Ball-Pivoting Algorithm in one C++ file

![Reconstruction of the bunny model](bunny.png)

A small, complete implementation of the Ball-Pivoting Algorithm for turning an oriented
point cloud into a triangle mesh:

> Fausto Bernardini, Joshua Mittleman, Holly Rushmeier, Cláudio Silva and Gabriel Taubin.
> *The Ball-Pivoting Algorithm for Surface Reconstruction.* IEEE TVCG 5(4), 1999.

The library is one file of about 550 lines, `src/lib/bpa.cpp`, and follows the paper: a uniform grid of cells
of side 2ρ, a seed search that visits every cell once, the pivot that returns the first point
the rolling ball touches, and the advancing front with the `join` and `glue` operators of
section 4.4. It depends only on [glm](https://github.com/g-truc/glm) for vectors.

This is a fork of [bernhardmgruber/bpa](https://github.com/bernhardmgruber/bpa), whose
clarity it keeps. The changes are listed at the end; in short, the fork seeds more than once,
pivots as the paper describes, works in double precision with relative tolerances, and
returns indices. It was checked against [BPA.jl](https://github.com/ctsilva/BPA.jl),
Open3D, MeshLab and the IPOL reference implementation on the same inputs; the numbers are
below.

## Building

Needs CMake 3.16, a C++20 compiler and glm; Catch2 v3 for the tests, which are skipped
without it.

```sh
# macOS: brew install glm catch2      Debian/Ubuntu: apt install libglm-dev catch2
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/tests
```

## Using it

```
./build/bpa <points> <radius>[,<radius>...] [<mesh>] [--min-component N] [--seed-neighbors N]
```

`<points>` is either an `.xyz` file with `x y z nx ny nz` per line or a `.noff`/`.off` file
with a `NOFF` header (vertices with normals). The normals decide which side of the surface
the ball rolls on; they need not be unit length. `<radius>` is the ball radius in the units
of the points; around 1.5 times the mean point spacing is a reasonable start. Several
radii, comma-separated, give one pass each from the smallest (section 4.6 of the paper): the
boundary edges left by one ball resume pivoting with the next, larger one, so that gaps the
small ball could not cross are closed while the fine detail it captured is kept. The mesh is
written as OFF, with the input vertices in their order, or as binary STL if the name ends in
`.stl`.

From C++:

```cpp
#include <bpa.h>

std::vector<bpa::Point> points = ...;              // {glm::dvec3 pos, glm::dvec3 normal}
std::vector<bpa::Face> faces = bpa::reconstruct(points, radius);   // or a std::vector of radii
// faces[i] = {a, b, c}: indices into `points`, counter-clockwise seen from outside
```

`reconstruct` takes an optional `bpa::Options`: `minComponent` drops connected components
with fewer triangles than that, useful for scans whose overlapping layers leave stray
fragments; `seedNeighbors` (100) bounds the pair search around a seed candidate to its
nearest neighbours, and 0 restores the paper's unbounded search.

Every triangle admits an empty ball of the given radius on its outward side, so the output
is a subset of the alpha shape of the points; the mesh is orientable and edge-manifold by
construction. Points the ball cannot reach with the given radius, because they lie in a
gap wider than the ball or under another layer of samples, are left out, and the mesh has a
boundary there. The paper's remedy is a second pass with a larger radius, which is what the
radius list does.

## How it compares

The implementations below were run on the same inputs by the comparison harness of
[BPA.jl](https://github.com/ctsilva/BPA.jl) (`compare/`), which checks every output on its
own terms: each triangle must admit an empty ball of the radius on its outward side, and the
mesh must be orientable and edge-manifold. Synthetic surfaces, the Stanford bunny scans
(single, four and ten merged) and the 62 dragon scans (1.83 million points); a MacBook Air
with an Apple M5, single-threaded, times for the reconstruction alone.

| case | points | ρ | this code | BPA.jl | Open3D | IPOL (Digne) | MeshLab |
|---|---|---|---|---|---|---|---|
| bunny, 1 scan | 40 256 | 1.25 mm | 78 162 triangles, 0.09 s | 78 152, 0.12 s | 77 994, 0.35 s | 77 941, 0.51 s | 78 203, 0.28 s |
| bunny, 10 scans | 362 272 | 1.25 mm | 323 920, 1.1 s | 323 934, 1.6 s | 317 975, 37 s | 323 808, 82 s | 477 737, 295 s |
| dragon, 62 scans | 1 830 000 | 0.7 mm | 649 459, 4.5 s | 649 518, 6.6 s | 624 851, 831 s | 631 174, 1987 s | 2 545 632, 5.4 h |

On the eight synthetic inputs (sphere, plane, four tori including an exact lattice, two knot
radii) this code and BPA.jl produce identical triangle sets. On the scans no triangle of
either has a non-empty ball; the ten-scan bunny has 17 components with both and the dragon
103 against 101, where Open3D and the IPOL code fragment into hundreds (97 and 536; 1899 and
549). MeshLab's pivot has no empty-ball test, which is why it produces more triangles. The
per-case reports, with renderings, are in the BPA.jl repository under `compare/results/`.

## What changed from bernhardmgruber/bpa

Each change is one commit, with its measurements in the message.

- **Seeds again after each front is exhausted** (section 4.2 of the paper, fig. 5). The original
  seeded once, so a scan came out as whichever component the first seed grew: 1631 of the
  78 000 triangles of a bunny scan. The seed search keeps a cursor over the cells, skips
  cells that already hold a used point (fig. 4c), tries one candidate per cell and pairs it
  with its nearest 100 neighbours; without that bound the search under an already
  reconstructed sheet cost 34 s instead of 1.4 s on the merged bunny scans.
- **Pivots as the paper describes.** The ball now rolls over every point within reach and
  stops at the first one it touches; the point is tested afterwards (normal, unused or on
  the front, manifold), and the edge becomes a boundary edge if it fails. The original
  filtered the candidates first, so the ball could pass through points to reach an
  acceptable one; it excluded the opposite vertex of the edge, so a triangle whose ball
  contained that vertex could be accepted; and it folded far-side angles as α + π instead
  of 2π − α, which ordered them backwards. Cospherical hits, the rule on a lattice, are
  resolved so that the pivots into a quad agree on its diagonal.
- **Seeds face along their vertex normals**, as the paper requires, instead of along the
  average normal of the cell, which on a thin tube points anywhere.
- **Double precision and a relative emptiness tolerance.** The original tested emptiness
  with an absolute margin of 1e-4 on the squared distance, which is 4 % of a ball of radius
  0.05 and larger than a ball of radius 0.01.
- **Pivots without trigonometry.** The order in which the rolling ball reaches the
  candidates is read off the 2-D positions of the ball centre on its circle, without
  computing an angle: no `atan2` or `acos` per candidate. Same output, 17 to 30 % faster
  on the scans (BPA.jl's `pivot_contact`).
- **Several radii** (section 4.6): after a pass, the boundary edges whose triangles admit an
  empty ball of the next radius are put back on the front with that ball, and pivoting
  resumes. The grid is rebuilt for each radius; the front and the mesh persist.
- **Indices, files, options.** `reconstruct` returns index triples; the driver reads XYZ and
  NOFF and writes OFF or STL; `Options` adds `minComponent` and `seedNeighbors`; the tests
  check the meshes (closed sphere with Euler characteristic 2, outward-facing triangles,
  two components from two spheres, component dropping) rather than only that they are not
  empty; the debug dumps are gone.

The bunny model is provided by the
[Stanford University Computer Graphics Laboratory](http://graphics.stanford.edu/data/3Dscanrep/).

## License

Boost Software License 1.0, as the original. See `LICENSE.txt`.
