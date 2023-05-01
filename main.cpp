#include "TN_KDE.hpp"

int main(int argc, char** agrv) {
    clock_t start, end;
    
    Model model;
    InitParameters(model, argc, agrv);
    LoadNetwork(model);

    start = clock();
    InitEdgeStructure(model);
    AddLixels(model);
    end = clock();
    cout << "Preprocess Time : " << (double)(end - start) / CLOCKS_PER_SEC << " s " << endl;
    model.prepare_time = (double)(end - start) / CLOCKS_PER_SEC;
    
    start = clock();
    TN_KDE(model);
    end = clock();
    cout << "Running Time : " << (double)(end - start) / CLOCKS_PER_SEC << " s " << endl;
    model.running_time = (double)(end - start) / CLOCKS_PER_SEC;
    cout << "Dijkstra Time : " << model.dijkstra_time / CLOCKS_PER_SEC << " s " << endl;
    
    Visualize(model);
}
