#PBS -N MPASMap 
#PBS -l walltime=2:00:00
#PBS -l nodes=1:ppn=1
#PBS -j oe

make

for ((i=0; i<12; i++)); do
	./MPASMap 0070_4.88364_578.19012_0.51473_227.95909 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0071_4.76373_659.84406_0.50628_279.50372 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0072_0.71478_881.02736_0.93875_178.47475 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0073_2.65927_529.08195_0.98952_218.40374 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0074_4.74674_692.31792_0.55049_159.21251 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0075_1.92645_426.88584_0.89721_247.26429 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0076_0.24047_1361.06602_0.91238_163.83314 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0077_2.51313_739.69077_0.49487_114.41521 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0078_4.34409_817.58109_0.90277_121.95627 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0079_3.43027_896.88548_0.27044_216.36624 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0080_3.49079_1106.38786_0.47725_173.44631 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0081_1.70972_522.88349_0.28171_207.02280 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0082_3.32575_306.11592_0.44649_263.19028 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0083_0.66710_865.73627_0.33224_293.92435 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0084_2.78137_303.29571_0.81286_247.29227 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0085_3.89023_979.79331_0.78842_156.21143 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0086_4.21449_1376.79154_0.34101_291.77479 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0087_1.24575_947.79456_0.31085_238.54343 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0088_4.61754_1395.80875_0.37381_236.79354 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0089_1.15061_507.89965_0.86718_277.42284 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0090_1.36644_1057.97160_0.73387_253.56073 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0091_2.77446_957.32246_0.82495_148.28125 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0092_4.39005_490.21656_0.50971_198.76377 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0093_0.25530_1083.39580_0.50553_162.80399 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0094_2.70176_1417.91279_0.52001_168.91967 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0095_4.52162_761.82800_0.98431_162.28460 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0096_0.00420_366.06495_0.50637_145.73566 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0097_0.03778_1319.40561_0.26651_237.70300 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0098_0.53889_1474.56055_0.81173_294.49235 $i
done
date
for ((i=0; i<12; i++)); do
    ./MPASMap 0099_4.52450_359.64880_0.86608_119.24767 $i
done
date