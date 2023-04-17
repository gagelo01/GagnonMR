gwasvcf::set_bcftools()
gwasvcf::set_plink()
test_that("intern function returns data.table", {
  uto <- GagnonMR:::test_object()
  expect_equal(intern_vcfpath_to_TwoSampelMR_region(vcf = uto$fn_exp1,
                                                    chrompos = uto$region["total"],
                                                    parameters = uto$parameters) %>% class,
               c("data.table", "data.frame"))
  expect_equal(intern_vcfpath_to_TwoSampelMR_region(vcf = uto$fn_exp1,
                                                    chrompos = uto$region["total"],
                                                    parameters = uto$parameters) %>% nrow %>%
           '>'(.,1), TRUE)
  expect_error(intern_vcfpath_to_TwoSampelMR_region(vcf = uto$fn_exp1,
                                                    chrompos = uto$region["total"],
                                                    parameters = uto$parameters), NA)
})

test_that("intern function gives right output", {
  uto <- GagnonMR:::test_object()
  expect_equal(intern_dt_null(vcffile_exp = uto$fn_exp1, vcffile_out = uto$fn_out, parameters = uto$parameters) %>% class,
               c("data.table", "data.frame"))
  expect_equal(intern_dt_null(vcffile_exp = uto$fn_exp1, vcffile_out = uto$fn_out, parameters = uto$parameters)$id.exposure,
               "bmi")
  expect_equal(intern_dt_null(vcffile_exp = uto$fn_exp1, vcffile_out = uto$fn_out, parameters = uto$parameters)$exposure,
               "bmi_ukbgiant")

})

test_that("different exposures returns", {
  uto <- GagnonMR:::test_object()
  expect_error(get_uni_cis(vcffile_exp = uto$fn_exp1,
              vcffile_out = uto$fn_out,
              chrompos = uto$region["subset"],
              parameters = uto$parameters), NA)

  expect_equal(get_uni_cis(vcffile_exp = uto$fn_exp1,
                           vcffile_out = uto$fn_out,
                           chrompos = uto$region["subset"],
                           parameters = uto$parameters) %>% nrow, 1)

  expect_equal(get_uni_cis(vcffile_exp = c(uto$fn_exp1, uto$fn_exp2),
                           vcffile_out = uto$fn_out,
                           chrompos = uto$region["subset"],
                           parameters = uto$parameters) %>% nrow, 2)

})

test_that("different exposure in get_coloc", {
  uto <- GagnonMR:::test_object()

  expect_equal(suppressWarnings(get_coloc(vcffile_exp = uto$fn_exp1,
                           vcffile_out = uto$fn_out,
                           chrompos = uto$region["subset"],
                           parameters = uto$parameters)) %>% nrow, 1)

  expect_error(get_coloc(vcffile_exp = c(uto$fn_exp1, uto$fn_exp2),
                           vcffile_out = uto$fn_out,
                           chrompos = uto$region["subset"],
                           parameters = uto$parameters))

})

test_that("different exposure in run_all_pqtl_analyses", {
  uto <- GagnonMR:::test_object()

  expect_equal(suppressWarnings(run_all_pqtl_analyses(vcffile_exp = uto$fn_exp1,
                                          vcffile_out = uto$fn_out,
                                          chrompos = uto$region["subset"],
                                     method_list = c("get_uni_cis", "get_coloc"),
                                          parameters = uto$parameters)) %>% nrow, 1)

  expect_equal(suppressWarnings(run_all_pqtl_analyses(vcffile_exp = c(uto$fn_exp1, uto$fn_exp2),
                         vcffile_out = uto$fn_out,
                         chrompos = uto$region["subset"],
                         method_list = c("get_uni_cis", "get_coloc"),
                         parameters = uto$parameters)) %>% nrow, 2)

})

test_that("coloc same", {
  uto <- GagnonMR:::test_object()
  dat_exposure <- GagnonMR:::intern_vcfpath_to_TwoSampelMR_region(vcf = uto$fn_exp1,
                                                                  chrompos = uto$region["subset"],
                                                                  parameters = uto$parameters)

  dat_outcome <- GagnonMR:::intern_vcfpath_to_TwoSampelMR_region(vcf = uto$fn_exp2,
                                                                 chrompos = uto$region["subset"],
                                                                 parameters = uto$parameters,
                                                                 type = "outcome")
  vout <- from_tsmr_to_coloc(dat_exposure = dat_exposure, dat_outcome = dat_outcome)
  expect_equal(suppressWarnings(get_coloc_intern(vout = vout,
                                dt_null = cbind(dat_exposure[,.(id.exposure, exposure)],
                                                dat_outcome[,.(id.outcome, outcome)])[1])),
               suppressWarnings(get_coloc(vcffile_exp = uto$fn_exp1,
                    vcffile_out = uto$fn_exp2,
                    chrompos = uto$region["subset"],
                    parameters = uto$parameters)))
})

