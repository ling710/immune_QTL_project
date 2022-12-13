#!/bin/bash
#PBS -r n
#
#PBS -l walltime=12:00:00,nodes=1:ppn=1
#PBS -j oe
#PBS -o job-report
#
# --------------------------------------------
#
cd pwcoco/pwcoco-master/build/
module load tools/cmake/3.16.3-gcc.9.3

./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/1_exp.txt  --sum_stats2  oneK1K/1_out.txt  --out oneK1Kres/1 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/2_exp.txt  --sum_stats2  oneK1K/2_out.txt  --out oneK1Kres/2 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/3_exp.txt  --sum_stats2  oneK1K/3_out.txt  --out oneK1Kres/3 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/4_exp.txt  --sum_stats2  oneK1K/4_out.txt  --out oneK1Kres/4 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/5_exp.txt  --sum_stats2  oneK1K/5_out.txt  --out oneK1Kres/5 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/6_exp.txt  --sum_stats2  oneK1K/6_out.txt  --out oneK1Kres/6 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/7_exp.txt  --sum_stats2  oneK1K/7_out.txt  --out oneK1Kres/7 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/8_exp.txt  --sum_stats2  oneK1K/8_out.txt  --out oneK1Kres/8 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/9_exp.txt  --sum_stats2  oneK1K/9_out.txt  --out oneK1Kres/9 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/10_exp.txt  --sum_stats2  oneK1K/10_out.txt  --out oneK1Kres/10 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/11_exp.txt  --sum_stats2  oneK1K/11_out.txt  --out oneK1Kres/11 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/12_exp.txt  --sum_stats2  oneK1K/12_out.txt  --out oneK1Kres/12 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/13_exp.txt  --sum_stats2  oneK1K/13_out.txt  --out oneK1Kres/13 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/14_exp.txt  --sum_stats2  oneK1K/14_out.txt  --out oneK1Kres/14 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/15_exp.txt  --sum_stats2  oneK1K/15_out.txt  --out oneK1Kres/15 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/16_exp.txt  --sum_stats2  oneK1K/16_out.txt  --out oneK1Kres/16 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/17_exp.txt  --sum_stats2  oneK1K/17_out.txt  --out oneK1Kres/17 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/18_exp.txt  --sum_stats2  oneK1K/18_out.txt  --out oneK1Kres/18 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/19_exp.txt  --sum_stats2  oneK1K/19_out.txt  --out oneK1Kres/19 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/20_exp.txt  --sum_stats2  oneK1K/20_out.txt  --out oneK1Kres/20 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/21_exp.txt  --sum_stats2  oneK1K/21_out.txt  --out oneK1Kres/21 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/22_exp.txt  --sum_stats2  oneK1K/22_out.txt  --out oneK1Kres/22 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/23_exp.txt  --sum_stats2  oneK1K/23_out.txt  --out oneK1Kres/23 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/24_exp.txt  --sum_stats2  oneK1K/24_out.txt  --out oneK1Kres/24 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/25_exp.txt  --sum_stats2  oneK1K/25_out.txt  --out oneK1Kres/25 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/26_exp.txt  --sum_stats2  oneK1K/26_out.txt  --out oneK1Kres/26 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/27_exp.txt  --sum_stats2  oneK1K/27_out.txt  --out oneK1Kres/27 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/28_exp.txt  --sum_stats2  oneK1K/28_out.txt  --out oneK1Kres/28 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/29_exp.txt  --sum_stats2  oneK1K/29_out.txt  --out oneK1Kres/29 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/30_exp.txt  --sum_stats2  oneK1K/30_out.txt  --out oneK1Kres/30 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/31_exp.txt  --sum_stats2  oneK1K/31_out.txt  --out oneK1Kres/31 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/32_exp.txt  --sum_stats2  oneK1K/32_out.txt  --out oneK1Kres/32 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/33_exp.txt  --sum_stats2  oneK1K/33_out.txt  --out oneK1Kres/33 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/34_exp.txt  --sum_stats2  oneK1K/34_out.txt  --out oneK1Kres/34 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/35_exp.txt  --sum_stats2  oneK1K/35_out.txt  --out oneK1Kres/35 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/36_exp.txt  --sum_stats2  oneK1K/36_out.txt  --out oneK1Kres/36 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/37_exp.txt  --sum_stats2  oneK1K/37_out.txt  --out oneK1Kres/37 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/38_exp.txt  --sum_stats2  oneK1K/38_out.txt  --out oneK1Kres/38 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/39_exp.txt  --sum_stats2  oneK1K/39_out.txt  --out oneK1Kres/39 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/40_exp.txt  --sum_stats2  oneK1K/40_out.txt  --out oneK1Kres/40 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/41_exp.txt  --sum_stats2  oneK1K/41_out.txt  --out oneK1Kres/41 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/42_exp.txt  --sum_stats2  oneK1K/42_out.txt  --out oneK1Kres/42 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/43_exp.txt  --sum_stats2  oneK1K/43_out.txt  --out oneK1Kres/43 --out_cond
./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/44_exp.txt  --sum_stats2  oneK1K/44_out.txt  --out oneK1Kres/44 --out_cond