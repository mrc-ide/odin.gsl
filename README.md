## odin.gsl

<!-- badges: start -->
[![Project Status: Concept – Minimal or no implementation has been done yet, or the repository is only intended to be a limited example, demo, or proof-of-concept.](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
[![Build Status](https://travis-ci.com/mrc-ide/odin.gsl.svg?branch=master)](https://travis-ci.com/mrc-ide/odin.gsl)
[![codecov.io](https://codecov.io/github/mrc-ide/odin.gsl/coverage.svg?branch=master)](https://codecov.io/github/mrc-ide/odin.gsl?branch=master)
<!-- badges: end -->

## Process

The model was generated using

```
odin::odin("inst/examples/sir.R", workdir = "src")
```

which generates `src/sir.c` - from that file we keep:

just the rhs

## License

MIT © Imperial College of Science, Technology and Medicine
