source('01_code/packages/packages.R')

#https://privefl.github.io/blog/a-guide-to-parallelism-in-r/
#https://nceas.github.io/oss-lessons/parallel-computing-in-r/parallel-computing-in-r.html
#https://blog.dominodatalab.com/multicore-data-science-r-python/
#https://github.com/ljdursi/beyond-single-core-R
#https://ljdursi.github.io/beyond-single-core-R/#/2

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)
foreach(i = 1:3, .combine = 'c') %dopar% {
    sqrt(i)
}
## [1] 1.000000 1.414214 1.732051
parallel::stopCluster(cl)
