# Testing Results for using par_bed in WDL
By using the following files, we can run a test to compare the output of our current WDL to the prior implementation, which didn't differentiate between PAR region and other regions.

Notice the difference when running the current WDL compared to the prior implementation. The mismatching GT fields in the X chromosome outside the PAR region, are now considered to be matches.
