#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>  


using namespace std;

// 線性反插值 (利用 std::lower_bound 進行內插)
double linear_inverse_interpolation(const vector<double>& xs, const vector<double>& ys, double y_target) {
    for (size_t i = 1; i < ys.size(); ++i) {
        if ((ys[i - 1] - y_target) * (ys[i] - y_target) <= 0) { // 包含 y_target
            double t = (y_target - ys[i - 1]) / (ys[i] - ys[i - 1]);
            return xs[i - 1] + t * (xs[i] - xs[i - 1]);
        }
    }
    // fallback: 若找不到區段，使用最簡單線性內插 (外插)
    return xs[0];
}

double inverse_interpolation(vector<double> xs, vector<double> ys, double y_target = 0.0, double tol = 1e-6, int max_iter = 100) {
    double x_new = xs[0];
    for (int iter = 0; iter < max_iter; ++iter) {
        // 反插值估算新的 x
        x_new = linear_inverse_interpolation(xs, ys, y_target);
        double y_new = x_new - exp(-x_new);

        if (abs(y_new - y_target) < tol)
            return x_new;

        // 將新估值加入資料點中
        xs.push_back(x_new);
        ys.push_back(y_new);

        // 根據 y 排序（從小到大）
        vector<size_t> indices(xs.size());
        iota(indices.begin(), indices.end(), 0);

        sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
            return ys[a] < ys[b];
        });

        vector<double> xs_sorted, ys_sorted;
        for (size_t i : indices) {
            xs_sorted.push_back(xs[i]);
            ys_sorted.push_back(ys[i]);
        }
        xs = xs_sorted;
        ys = ys_sorted;
    }

    return x_new;
}

int main() {
    vector<double> x_values = {0.3, 0.4, 0.5, 0.6};
    vector<double> y_values;

    for (double x : x_values)
        y_values.push_back(x - exp(-x));

    double target_y = 0.0;
    double solution = inverse_interpolation(x_values, y_values, target_y);

    cout.precision(8);
    cout << "Approximate solution: x ≈ " << solution << endl;

    return 0;
}
