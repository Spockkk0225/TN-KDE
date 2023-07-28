#include "TN_KDE.hpp"

int main(int argc, char** agrv) {
    clock_t start, end;

    Model model;
    InitParameters(model, argc, agrv);
    LoadNetwork(model);

    start = clock();
    AddLixels(model);
    end = clock();
    cout << "Preprocess Time : " << (double)(end - start) / CLOCKS_PER_SEC << " s " << endl;
    model.prepare_time = (double)(end - start) / CLOCKS_PER_SEC;

//    double init_time = 0;
//    double kde_time = 0;
    start = clock();
    for (int i = 0;i < model.windows_num; i++) {
        double start = clock();
        InitEdgeStructure(model, i);
        double end = clock();
        cout << "Init Time : " << (end - start) / CLOCKS_PER_SEC << " s " << endl;
        start = clock();
        TN_KDE(model);
        end = clock();
        cout << "KDE Time : " << (end - start) / CLOCKS_PER_SEC << " s " << endl;
//        Visualize(model, i);
    }
    end = clock();

//    cout << "Init Time : " << init_time / CLOCKS_PER_SEC << " s " << endl;
//    cout << "KDE Time : " << kde_time / CLOCKS_PER_SEC << " s " << endl;

    cout << "Running Time : " << (double)(end - start) / CLOCKS_PER_SEC << " s " << endl;
    model.running_time = (double)(end - start) / CLOCKS_PER_SEC;
//    cout << "Dijkstra Time : " << model.dijkstra_time << " s " << endl;

    Visualize(model, 0);
    Statistic(model);
//    cout << model.a << " " << model.b << " " << model.c << " " << model.d << endl;
    return 0;
}