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
for (v in c("PkgIntro", "SSR", "SLM", "SCT", "GCCM", "GCMC", "SCPCM")) {
  .prebuild_vignettes(v)
}
```

- Build vignettes separately

``` r
.prebuild_vignettes("PkgIntro")
.prebuild_vignettes("SSR")
.prebuild_vignettes("SLM")
.prebuild_vignettes("SCT")
.prebuild_vignettes("GCCM")
.prebuild_vignettes("GCMC")
.prebuild_vignettes("SCPCM")
```
