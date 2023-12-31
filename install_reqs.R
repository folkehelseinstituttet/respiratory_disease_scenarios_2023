

drat:::add("ncov-ic")
install.packages(c("abind", "tidyverse", "data.table", "odin", "odin.dust", "dust", "mcstate", "doBy", "lhs", "truncnorm"))

options(
  repos = structure(c(
    SPLVERSE  = "https://docs.sykdomspulsen.no/drat/",
    CRAN      = "https://cran.rstudio.com"
  ))
)

install.packages(c("spldata", "splstyle"))


devtools::install_github(repo="https://github.com/Gulfa/metapop")
devtools::install_github(repo="https://github.com/Gulfa/metapopnorge")
