# Contributing to Cysinger

Thank you for your interest in contributing to Cysinger! This document provides guidelines for contributing.

## How to Contribute

### Reporting Bugs

1. Check existing [Issues](https://github.com/go2bean/CysingerR/issues) to avoid duplicates
2. Open a new issue with:
   - A clear, descriptive title
   - Steps to reproduce the problem
   - Expected vs actual behavior
   - R version and OS information (`sessionInfo()`)

### Suggesting Features

Open an issue with the **"Feature Request"** label describing:
- The problem your feature would solve
- Your proposed solution
- Any alternative approaches you've considered

### Pull Requests

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/your-feature`)
3. Make your changes following the code style below
4. Add or update tests in `tests/testthat/`
5. Update documentation (roxygen2 comments)
6. Run `R CMD check` and ensure no errors
7. Submit a pull request

## Code Style

- Follow the [tidyverse style guide](https://style.tidyverse.org/)
- Use roxygen2 for documentation
- Use `rlang::.data` for tidy evaluation in ggplot2/dplyr
- Keep functions focused and well-documented

## Development Setup

```r
# Install development dependencies
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown"))

# Build and check
devtools::document()
devtools::test()
devtools::check()
```

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
