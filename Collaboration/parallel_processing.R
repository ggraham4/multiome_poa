library(parallel)

test_function <- function(x){
  set.seed(x)
  random_data <- data.frame(col_1 = rnorm(100000),
                            col_2 = rnorm(100000))

  
  model <- lm(col_1~col_2, data = random_data)
  return()
  
}

#running serially
{
  start = (Sys.time())
  
  lapply(0:10000, test_function)
  
  end = Sys.time()
  print(end-start)
}
#2.36 mis

#running in parallel -- on mac
{
  start = (Sys.time())
  
  parallel::mclapply(X= 0:10000, FUN=test_function, mc.cores = detectCores()-1)
  
  end = Sys.time()
  print(end-start)
}
#35.5 s

#running in parallel -- on windows using parallelsugar wrapper on windows

#installing parallelsugar
library(devtools)
install_git('https://github.com/nathanvan/parallelsugar.git')
library(parallelsugar)
{
    start = (Sys.time())
  
  parallelsugar::mclapply(X= 0:10000, FUN=test_function, mc.cores = detectCores()-1)
  
  end = Sys.time()
  print(end-start)
  
}
#1.00 mins - mac chips are so fast