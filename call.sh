g++ main.cpp TN_KDE.cpp kernel.cpp init.cpp network.cpp edge.cpp library.cpp -std=c++11 -O3 -o main

dataset="JohnsCreek"
statistic="statistic.txt"
method="RFS"
kernel="Triangular"
prune="ON"
lixel_interval=10
spatial_bandwidth=1000
H=0
window_num=10

rm $statistic

./main "network."$dataset".txt" "Output/out."$dataset"_"$kernel"_l"$lixel_interval"_b"$spatial_bandwidth"_"$method"_H="$H".txt" $statistic $method $kernel $prune $lixel_interval $spatial_bandwidth $H $window_num \
0.10 0.80 \
0.11 0.81 \
0.12 0.82 \
0.13 0.83 \
0.14 0.84 \
0.15 0.85 \
0.16 0.86 \
0.17 0.87 \
0.18 0.88 \
0.19 0.89
