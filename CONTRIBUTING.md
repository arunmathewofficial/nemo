# Contributing to NEMO

Thank you for your interest in NEMO.

## How the repositories work

- **This GitHub repository is public** and holds the released versions:
  `main` and the release tags (`vX.Y.Z`). The Python package
  [`nemo-astro`](https://pypi.org/project/nemo-astro/) is built from it.
- **Development happens in a private repository.** New features and
  changes that are not released yet are there. `main` here is updated from
  it at each release; it never receives changes from anywhere else.
- If you would like access to the development version (for example, to
  work on a feature together with us), contact the maintainers:
  [Arun Mathew](https://github.com/arunmathewofficial) and
  [Jonathan Mackey](https://github.com/jmackey-astro).

## Issues

Please report bugs, and suggest new features, as
[GitHub issues](https://github.com/arunmathewofficial/nemo/issues). For a
bug, the NEMO version (or commit), how you built it, and a small input
that shows the problem help a lot.

## Pull requests

Bug fixes and new features are welcome as pull requests against `main`.

- **Pull requests are not merged on GitHub.** We bring your commits into
  the development repository, test them there, and they appear on GitHub
  with the next release. We then close the pull request with a link to the
  release.
- **You stay the author of your commits.** Git records the author of each
  commit separately from the person who applies it, and we keep the
  author when we bring your commits in. GitHub credits you as the author
  once they appear in a release here. If we have to combine or rewrite
  commits, we add a `Co-authored-by:` line for you.
  - For GitHub to link the commits to your account, commit with an email
    address that is on your GitHub account (or your GitHub `noreply`
    address).
- `main` here can be behind the development version, so we may have to
  adapt your change before merging it. We will ask you if anything is
  unclear.
- Please:
  - keep a pull request to one change, and explain what it does and why;
  - add or update a test in `tests/` (C++) or `python/tests/` (Python) where
    possible;
  - format C++ code with `tools/format.sh` (clang-format 18).

By contributing you agree that your contribution is licensed under the
[BSD 3-Clause License](LICENSE) of this project.
