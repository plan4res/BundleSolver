# BundleSolver

BundleSolver class, which implements the CDASolver interface within the SMS++
framework. The class implements a NonDifferentiable Optimization Solver using
a "Generalized Bundle" algorithm; cf. e.g.

- [A. Frangioni "Generalized Bundle Methods" SIAM Journal on Optimization
  13(1), p. 117 - 156,
  2002](http://www.di.unipi.it/~frangio/abstracts.html#SIOPT02)

- [A. Frangioni, E. Gorgone "Generalized Bundle Methods for Sum-Functions
  with "Easy" Components: Applications to Multicommodity Network Design"
  Mathematical Programming 145(1), 133–161,
  2014](http://pages.di.unipi.it/frangio/abstracts.html#MP11c)

- [W. van Ackooij, A. Frangioni "Incremental Bundle Methods Using Upper
  Models" SIAM Journal on Optimization 28(1), 379–410,
  2018](http://pages.di.unipi.it/frangio/abstracts.html#SIOPT16)

- [A. Frangioni "Standard Bundle Methods: Untrusted Models and Duality" in
  Numerical Nonsmooth Optimization: State of the Art Algorithms, A.M. Bagirov,
  M. Gaudioso, N. Karmitsa, M. Mäkelä, S. Taheri (Eds.), 61—116, Springer,
  2020](http://pages.di.unipi.it/frangio/abstracts.html#NDOB18)


## Getting started

These instructions will let you build MCFBlock and MCFSolver on your system.

### Requirements

- The [SMS++ core library](https://gitlab.com/smspp/smspp) and its requirements.

- The [MILPSolver](https://gitlab.com/smspp/milpsolver) SMS++ module.

- The [NDOSolver/FiOracle project](https://gitlab.com/frangio68/ndosolver_fioracle_project)
  and its requirements (depending on the actual MPSolver built); note that this
  dependency is supposed to be removed down the line.

### Build and install with CMake

Configure and build the library with:

```sh
mkdir build
cd build
cmake ..
make
```
The library has the same configuration options of
[SMS++](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration).

You can also choose the following configuration options:

| Variable       | Description         | Default value |
|----------------|---------------------|---------------|
| `WHICH_OSI_QP` | Use CPLEX or GUROBI | GUROBI        |

You can set them with:

```sh
cmake <source-path> -D<var>=<value>
```

Optionally, install the library in the system with:

```sh
sudo make install
```

### Usage with CMake

After the library is built, you can use it in your CMake project with:

```cmake
find_package(BundleSolver)
target_link_libraries(<my_target> SMS++::BundleSolver)
```

### Build and install with makefiles

Carefully hand-crafted makefiles have also been developed for those unwilling
to use CMake. Makefiles build the executable in-source (in the same directory
tree where the code is) as opposed to out-of-source (in the copy of the
directory tree constructed in the build/ folder) and therefore it is more
convenient when having to recompile often, such as when developing/debugging
a new module, as opposed to the compile-and-forget usage envisioned by CMake.

Each executable using `BundleSolver` or `ParallelBundleSolver`,
such as the [tester for `PolyhedralFunction` comparing `BundleSolver` with
`MILPSolver`](https://gitlab.com/smspp/tests/-/tree/develop/PolyhedralFunction?ref_type=heads),
has to include a "main makefile" of the module, which typically is either
[makefile-c](makefile-c) including all necessary libraries comprised the
"core SMS++" one, or [makefile-s](makefile-s) including all necessary
libraries but not the "core SMS++" one (for the common case in which this is
used together with other modules that already include them). These in turn
recursively include all the required other makefiles, hence one should only
need to edit the "main makefile" for compilation type (C++ compiler and its
options) and it all should be good to go. In case some of the external
libraries are not at their default location, it should only be necessary to
create the `../extlib/makefile-paths` out of the
`extlib/makefile-default-paths-*` for your OS `*` and edit the relevant bits
(commenting out all the rest).

Check the [SMS++ installation wiki](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration#location-of-required-libraries)
for further details.

Note that the [NDOSolver/FiOracle
project](https://gitlab.com/frangio68/ndosolver_fioracle_project) has a similar
arrangement with its own extlib/ folder, but the `*_ROOT` values are set in the
SMS++ files and therefore are immediately available there, so there is no need
to separately edit the NDOSolver/FiOracle project ones (but there would be if
it were downloaded and compiled independenty).


## Getting help

If you need support, you want to submit bugs or propose a new feature, you can
[open a new issue](https://gitlab.com/smspp/bundlesolver/-/issues/new).


## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of
conduct, and the process for submitting merge requests to us.


## Authors

### Current Lead Authors

- **Antonio Frangioni**  
  Dipartimento di Informatica  
  Università di Pisa

- **Enrico Gorgone**  
  Dipartimento di Matematica ed Informatica  
  Università di Cagliari


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html) -
see the [LICENSE](LICENSE) file for details.


## Disclaimer

The code is currently provided free of charge under an open-source license.
As such, it is provided "*as is*", without any explicit or implicit warranty
that it will properly behave or it will suit your needs. The Authors of
the code cannot be considered liable, either directly or indirectly, for
any damage or loss that anybody could suffer for having used it. More
details about the non-warranty attached to this code are available in the
license description file.
