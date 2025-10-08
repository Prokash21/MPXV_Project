# Load the package
library(GOplot)

# Get the version of the package
package_version <- packageVersion("GOplot")
print(paste("Version:", package_version))

# Get the citation information
package_citation <- citation("GOplot")
print(package_citation)

citation("lfcShrink")
citation("stats")
