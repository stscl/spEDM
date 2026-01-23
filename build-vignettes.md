# NA

Due to the time-consuming computations involved in the vignettes of the
*spEDM* package, it is necessary to pre-build the vignettes prior to
package submission.

``` r
.prebuild_vignettes = \(name){
  out = paste0("vignettes/",name,".Rmd")
  inp = paste0(out,".orig")
  knitr::knit(inp,out)
}
```

- Compile all vignettes at once

``` r
# list vignette names
vignettes = paste0(c(paste0("main", 1:6), 
                     paste0("si", 1:2)), "_",
                   c("pkgintro", "ssr", "gccm", "gpc",
                     "gcmc", "scpcm", "slm", "sct"))

for (v in vignettes) {
  .prebuild_vignettes(v)
}
```

- Build vignettes separately

``` r
.prebuild_vignettes("main1_pkgintro")
.prebuild_vignettes("main2_ssr")
.prebuild_vignettes("main3_gccm")
.prebuild_vignettes("main4_gpc")
.prebuild_vignettes("main5_gcmc")
.prebuild_vignettes("main6_scpcm")
.prebuild_vignettes("si1_slm")
.prebuild_vignettes("si2_sct")
```
