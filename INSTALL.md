# Installation

## Prerequisites

- Git
- Python 3.11+
- Optional: LaTeX distribution for manuscript builds

## Clone

```bash
git clone git@github.com:alexanderholman/mp1-1-sigma3-112-repro.git
cd mp1-1-sigma3-112-repro
```

## Update

```bash
git pull --ff-only
```

## MP1 tools

Install the pinned `mp1-tools` version recorded in `mp1_tools.lock`:

```bash
python -m pip install "git+ssh://git@github.com/alexanderholman/mp1-tools.git@e8b48cb25a7b810cf15c2599db27781dcba0cf6d"
```

## Remove

Delete the cloned directory when no longer needed.
