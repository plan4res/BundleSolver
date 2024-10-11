# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.4.4] - 2024-02-27

### Changed

- makefiles and Cmake files updated to new global SMS++ scheme

- documentation updated accordingly

- default OSISolver is Gurobi

- improved BundleSolverState by keeping all upper and lower values for each component

### Fixed

- fixed sloppy setting of OSI parameters

- fixed `has_*_solution()` and `has_*_feasible()`

## [0.4.3] - 2023-05-23

### Added

- handling the constant term Modification for hard
  components

- handling constant term for "easy" components

- support for OsiGrbSolverInterface

### Changed

- better printing in case of errors in a Fi

### Fixed

- missing check of max number of iterations in one case

- issues in put_State()

- issue in one case of process_outstanding_Modification (adding
  Linearization was not properly managed)

- management of numerical errors in the MPSolver

- typos and spaces

## [0.4.2] - 2022-06-28

### Added

- added new parameter intTrgtMng controlling different aspects of how models
  are used to pass targets and accuracy to the C05Function and to forecast
  the value (using the Lipschitz_constant)

### Changed

- generate the abstract representation in set_Block()

- minor improvements in logging and comments

- quite significant rehaul of t heuristics: reshuffled the bits in tSPar1
  (baiscally keeping compatibility with develop branch) and added the "reversal
  form of the poorman's quasi-Newton update". This requires the computation of
  the norm of the "aggregated representative subgradient", doing which called
  for some nontrivial changes in the computation of Alfa1 and ScPr1

### Fixed

- ParallelBundleSolver now properly counts the computed components; previously
  those of the ramp-down phase were not, which may lead to getting an error if
  the last step was not a SS and intFrcLstSS was true

- fixed get\_dual\_solution() when there is no easy component

## [0.4.1] - 2021-12-07

- improved Modification handling (no over-reacting to easy ones, important bugfix)

- better log, printing times for each component in verbosity 4 and higher

- fixed several flaws

- improved namespace handling and similar stuff

### Changed

## [0.4.0] - 2021-05-02

### Changed

Several major improvements:

- restructured stopping conditions, with a new one on the norm (INF, 1 or 2)
  of zStar, possibly relative to subgradient ; tStar anyway remains the driver
  of t-strategies.

- BundleSolver now detects when "zStar == 0" (see previous point) and uses
  this to algorithmically produce globally valid lower bounds.

- BundleSolver now checks if the objectives are (all) convex or concave and
  adapts accordingly (by sneakily converting the concave problem internally
  into a convex one and adapting the results on the fly when they go out).

- Box constraints are now handled, either by means of NNConstraint/BoxConstraint
  or as NN constraints "inherent" in the ColVariable.

- Changes in the costs or RHS/bounds in easy components now are allowed and managed.

- Changes in the 0-th component now allowed.

- Added hooks for perspective handling of changes in easy components.

- Complete rehaul of Modification in "easy" components.

- Better bad termination reporting.

- Allowing to switch away Cplex.

- Added ComputeConfig for hard and easy and NoEasy.

- Added individual Configurations.

### Fixed

- Many fixes throughout the code, too many to list.

## [0.3.0] - 2020-09-16

### Changed

- Support for better modification framework.

### Fixed

- Minor fixes.

## [0.2.0] - 2020-03-06

### Fixed

- Minor fixes.

## [0.1.2] - 2020-03-02

### Fixed

- Minor fixes.

## [0.1.1] - 2020-02-17

### Fixed

- Minor bug in CMake configuration file.

## [0.1.0] - 2020-02-10

### Added

- First test release.

[Unreleased]: https://gitlab.com/smspp/bundlesolver/-/compare/0.4.3...develop
[0.4.3]: https://gitlab.com/smspp/bundlesolver/-/compare/0.4.2...0.4.3
[0.4.2]: https://gitlab.com/smspp/bundlesolver/-/compare/0.4.1...0.4.2
[0.4.1]: https://gitlab.com/smspp/bundlesolver/-/compare/0.4.0...0.4.1
[0.4.0]: https://gitlab.com/smspp/bundlesolver/-/compare/0.3.0...0.4.0
[0.3.0]: https://gitlab.com/smspp/bundlesolver/-/compare/0.2.0...0.3.0
[0.2.0]: https://gitlab.com/smspp/bundlesolver/-/compare/0.1.2...0.2.0
[0.1.2]: https://gitlab.com/smspp/bundlesolver/-/compare/0.1.1...0.1.2
[0.1.1]: https://gitlab.com/smspp/bundlesolver/-/compare/0.1.0...0.1.1
[0.1.0]: https://gitlab.com/smspp/bundlesolver/-/tags/0.1.0
