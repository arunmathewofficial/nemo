# NEMO example: logging CVODE messages (by Claude AI)

A standalone project, like [examples/minimal](../minimal), that shows how to
send NEMO's CVODE (SUNDIALS) error messages to their own log file.  A hot,
dense H+He gas cell (73,000 K) cools over 5 steps of 1e14 s each.  The steps
are much longer than the cooling time, so CVODE often reaches its limit on
internal steps ("mxstep") and NEMO splits the step.

How it works:
- NEMO sends CVODE error messages to the spdlog logger named
  `NEMO::sundials_logger_name` (`"sundials"`).  The example creates that
  logger with a file sink, `sundials.log`.  Without it, the messages go to
  spdlog's default logger.
- mxstep errors are not logged one by one, because NEMO recovers from them.
  NEMO counts them (summed over all threads), and
  `NEMO::get_and_reset_cvode_mxstep_count()` returns and resets the count.
  The example logs the count once per step.
- Don't call `spdlog::shutdown()` while a `NEMO::nemo` object exists: it
  also removes the default logger, which NEMO's destructor logs through.
  The example calls `flush()` on the logger instead.

First install NEMO (see the top-level [README](../../README.md)), e.g. into
`~/.local`.  Then, from this directory:
```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH="$HOME/.local"
cmake --build build
cd build && ./nemo-sundials-logger
```

The output should be similar to this (the numbers depend on NEMO's cooling
tables):
```
Initial temperature: 73406 K
step 1: T = 5562.43 K
step 2: T = 5434.23 K
step 3: T = 5390.52 K
step 4: T = 5375.52 K
step 5: T = 5370.34 K
CVODE messages are in sundials.log
```
and `sundials.log` contains a line like
```
[2026-11-11 22:45:10.538] [sundials] [warning] step 1: 29 CVODE mxstep errors
```
