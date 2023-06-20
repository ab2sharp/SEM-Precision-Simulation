#!/bin/bash

nohup Rscript simpara__semRB_(nobs-p-iters).R 500 1 20 &
wait
nohup Rscript simpara__semRB_(nobs-p-iters).R 500 1 50 &
wait
nohup Rscript simpara__semRB_(nobs-p-iters).R 500 1 100 &
wait
nohup Rscript simpara__semRB_(nobs-p-iters).R 500 9 20 &
wait
nohup Rscript simpara__semRB_(nobs-p-iters).R 500 9 50 &
wait
nohup Rscript simpara__semRB_(nobs-p-iters).R 500 9 100 &
wait
nohup Rscript simpara__semRB_(nobs-p-iters).R 500 49 20 &
wait
nohup Rscript simpara__semRB_(nobs-p-iters).R 500 49 50 &
wait
nohup Rscript simpara__semRB_(nobs-p-iters).R 500 49 100 &
wait
nohup Rscript simpara__semRB_(nobs-p-iters).R 500 99 20 &
wait
nohup Rscript simpara__semRB_(nobs-p-iters).R 500 99 50 &
wait
nohup Rscript simpara__semRB_(nobs-p-iters).R 500 99 100 &
