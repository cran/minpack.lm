".onLoad" <- function (lib, pack)
{
	library.dynam(pack, pack, lib)

}

".onAttach" <- function (lib, pack)
{
  cat("
This product includes software developed by the University of Chicago, as
Operator of Argonne National Laboratory.

See the LICENSE file distributed with the minpack.lm source code or
http://www.netlib.org/minpack/disclaimer for the full license.\n\n")

}
