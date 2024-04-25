# geposan

Geposan is an R package for analyzing genes based on their position across
species. The package includes gene data from Ensembl. It provides multiple
methods to use that data to find genes that score well in comparison with a
set of reference genes.

## Installation

You can install the development version of geposan using:

```r
# install.packages("remotes")
remotes::install_github("johrpan/geposan")
```

See [this page](https://remotes.r-lib.org/reference/install_git.html) for more
information on this command.

## Publication

This method and its implementation have been peer-reviewed and published in
NAR Genomics and Bioinformatics. If you use the package in your research or
would like to refer to our methodology, please cite the following paper:

Elias F Projahn, Georg Fuellen, Michael Walter, Steffen Möller, Proposing
candidate genes under telomeric control based on cross-species position data,
NAR Genomics and Bioinformatics, Volume 6, Issue 2, June 2024, lqae037,
https://doi.org/10.1093/nargab/lqae037

## Graphical interface

Please also take a look at the interactive graphical web interface for geposan
which is available as a separate R package called
[geposanui](https://github.com/johrpan/geposanui).

## License

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU Affero General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the
[GNU Affero General Public License](https://www.gnu.org/licenses/agpl-3.0.html)
for more details.
