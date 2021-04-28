.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Written by Jing Zhang, Ph.D. Please direct questions to jzhangcad@gmail.com.

  To cite in publication:

    Zhang J, Hadj-Moussa H, Storey KB. 2020. Marine periwinkle stress-responsive microRNAs: a potential factor to reflect anoxia and freezing survival adaptations. GENOMICS. 2020 Jul 27: S0888-7543(20)30169-5. doi: 10.1016/j.ygeno.2020.07.036.
    Zhang J, Wallace SJ, Shiu MY, Smith I, Rhind SG, Langlois VS. 2017. Human hair follicle transcriptome profiling: a minimally invasive tool to assess molecular adaptations upon low-volume, high-intensity interval training. Physiological Reports. 5(23) pii: e13534. doi: 10.14814/phy2.13534.")
  # suppressPackageStartupMessages(require(pathview))
  return(TRUE)
}
