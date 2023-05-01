g++ main.cpp TN_KDE.cpp kernel.cpp init.cpp network.cpp edge.cpp library.cpp -std=c++11 -O3 -o main


dataset="JohnsCreek"
statistic="statistic.txt"
method="RFS"
kernel="Triangular"
lixel_interval=50
spatial_bandwidth=1000
depth=0
start_time=1400000000
end_time=1500000000

rm $statistic

./main "network."$dataset".txt" "out."$dataset"_"$kernel"_l"$lixel_interval"_b"$spatial_bandwidth"_"$method".txt" $statistic $method $kernel $lixel_interval $spatial_bandwidth $depth $start_time $end_time
