# Building and releasing the Python package

## What builds what

- **GitLab CI, job `Python`:** builds the sdist, then manylinux_2_28 x86_64
  wheels from it for CPython 3.10–3.15, and tests them. It runs the quick
  tests for every version and all tests for 3.13. Nothing is uploaded, and the
  files are kept as job artifacts for a week.
- **GitHub Actions, `.github/workflows/wheels.yml`:** builds the sdist and the
  wheels for Linux (x86_64, aarch64) and macOS (Intel, Apple silicon) with
  cibuildwheel, and tests them. The workflow runs on pull requests and on
  pushes to `main` of the GitHub repository. It uploads only when started by
  hand ("Run workflow") with `publish` set to `testpypi` or `pypi`.

The cibuildwheel settings (Python versions, platforms, test command) are in
`[tool.cibuildwheel]` in `pyproject.toml`. To build Linux wheels locally,
install docker or podman and run `pipx run cibuildwheel --platform linux`
(with podman, set `CIBW_CONTAINER_ENGINE=podman`).

## The package version

The version comes from `project(NEMO VERSION ...)` in the top-level
`CMakeLists.txt`. PyPI never accepts the same version twice, even after it
has been deleted, so increase the version for every upload. For a
Python-only fix, use a post-release, e.g. `1.4.0.post1`. That needs a change
to the version regex in `pyproject.toml`.

## First release

1. Create an account on https://pypi.org and on https://test.pypi.org (they
   are separate), and turn on two-factor authentication. Uploading is
   allowed only with 2FA.
2. On each site, under "Your projects", then "Publishing", add a *pending
   trusted publisher*: the project name, the GitHub owner and repository,
   workflow `wheels.yml`, and environment `pypi` (on test.pypi.org:
   `testpypi`). A pending publisher does not reserve the name; the first
   upload creates the project.
3. In the GitHub repository, under Settings, then Environments, create the
   environments `pypi` and `testpypi`. Requiring a reviewer for `pypi` adds a
   manual approval step.
4. Run the workflow with `publish = testpypi`, then check that it installs:
   `pip install -i https://test.pypi.org/simple/ --extra-index-url
   https://pypi.org/simple/ <name>`.
5. Run the workflow with `publish = pypi`.

The first upload to pypi.org fixes the project name for good. Choose the
name first: see below.

## Changing the package name

The distribution name (what `pip install` uses) and the import name can
differ, but keeping them alike is less confusing. The name currently
appears in these places:

- `pyproject.toml`: `name`, and `wheel.packages = ["python/<import name>"]`
- `python/<import name>/`: the package directory
- `python/CMakeLists.txt`: `install(... DESTINATION <import name>)`
- `python/src/bindings.cpp`: the module docstring
- `python/README.md`, `python/<import name>/__init__.py` (docstring),
  `python/tests/*.py` (`import ... as nemo`)
- the comment on `NEMO_PYTHON` in the top-level `CMakeLists.txt`

`grep -rn "nemo[-_]plasma" --exclude-dir=build-python .` lists them.
