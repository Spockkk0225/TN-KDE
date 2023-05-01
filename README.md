# STNKDE

This is the code of the Temporal Network Kernel Density Estimation (TN-KDE).

### 1. Build & Run

To build and run this project, you just need to execute `call.sh` in the terminal.

This file includes:

```
first line: compiling command with c++11 and O3 optimization

dataset: network and event datasets(file name is network.XXX.txt)

statistic: statistics output file(parameters, time, etc.)

method: algorithms methods(RQS, SPS, RAS, RTS, RFS, DRFS)

kernel: kernel functions(Triangular, Cosine, Exponential)

lixel_interval: length of each lixel(10, 20, 30, 40, 50 in the paper)

spatial_bandwidth: bandwidth of distance(1000, 2000, 4000, 7000, 10000 in the paper)

h_value: tree depth used in DRFS

start_time/end_time (optional): specified querying time range, [20%,80%] as default if not given
```

### 2. Dataset

Since the dataset file is too large, we do not upload it. Here is the structure of dataset files:

```
<Number of Nodes> <Number of Edges>
(Edges 1) <Endpoint 1> <Endpoint 2> <Length>
(Edges 2) <Endpoint 1> <Endpoint 2> <Length>
...
(Edges m) <Endpoint 1> <Endpoint 2> <Length>

(Event 1) <Endpoint 1> <Endpoint 2> <Distance from Endpoint 1> <Timestamp>
...
(Event ?) <Endpoint 1> <Endpoint 2> <Distance from Endpoint 1> <Timestamp>
```

The first line is the number of nodes and the number of edges.

Then each line is an edge, including two endpoints ID and its length.

Finally each line is an event, including two endpoints of the edge, the distance away from the frist endpoint, and its time.

Part of `network.JohnsCreek.txt`:

```
3074 3471
1 2 178.64384457327444
2 3 168.02959873043565
2 552 315.866481429132
4 5 420.8691863406573
...
...
1645 2570 11.511241466084144 0.716155452365731 1251732711.0
1005 2555 120.40317486357394 5.985442657243191 1251732775.0
731 1637 252.06724715224766 203.03331062005964 1251733324.0
120 2279 71.13009119351476 6.07249890140156 1251733532.0
...
...
```

