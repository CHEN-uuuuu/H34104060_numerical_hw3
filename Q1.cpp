#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

// 計算階乘
int factorial(int n) {
    int result = 1;
    for (int i = 2; i <= n; ++i)
        result *= i;
    return result;
}

// 計算誤差界限
double error_bound(const vector<double>& x_vals, double x_target, int degree) {
    double max_derivative = 1.0; // cos(x) 的高階導數最大值為 ±1
    double product_term = 1.0;

    for (int i = 0; i <= degree; ++i)
        product_term *= abs(x_target - x_vals[i]);

    return (max_derivative / factorial(degree + 1)) * product_term;
}

// Lagrange 插值函數
pair<double, double> lagrange_interpolation(const vector<double>& x_vals, const vector<double>& y_vals, double x_target, int degree) {
    int n = x_vals.size();

    // 找到最接近 x_target 的 (degree+1) 個點
    vector<int> indices(n);
    for (int i = 0; i < n; ++i)
        indices[i] = i;

    sort(indices.begin(), indices.end(), [&](int a, int b) {
        return abs(x_vals[a] - x_target) < abs(x_vals[b] - x_target);
    });

    vector<double> x_subset(degree + 1);
    vector<double> y_subset(degree + 1);
    for (int i = 0; i <= degree; ++i) {
        x_subset[i] = x_vals[indices[i]];
        y_subset[i] = y_vals[indices[i]];
    }

    // 建立 Lagrange 多項式
    double approx_value = 0.0;
    for (int i = 0; i <= degree; ++i) {
        double term = y_subset[i];
        for (int j = 0; j <= degree; ++j) {
            if (j != i)
                term *= (x_target - x_subset[j]) / (x_subset[i] - x_subset[j]);
        }
        approx_value += term;
    }

    // 計算誤差界限
    double error = error_bound(x_subset, x_target, degree);

    return {approx_value, error};
}

int main() {
    vector<double> x_values = {0.698, 0.733, 0.768, 0.803};
    vector<double> y_values;
    for (double x : x_values)
        y_values.push_back(cos(x));

    double x_target = 0.750;
    double true_value = cos(x_target);

    int max_degree = x_values.size() - 1;

    for (int deg = 1; deg <= 4; ++deg) {
        if (deg > max_degree) {
            cout << deg << " 次插值無法執行：資料點不足（只有 " << x_values.size() << " 個點）" << endl;
            continue;
        }

        auto [approx_value, error] = lagrange_interpolation(x_values, y_values, x_target, deg);
        double actual_error = abs(true_value - approx_value);

        cout << deg << " 次插值結果: P_" << deg << "(" << x_target << ") ≈ "
             << approx_value << ", 誤差界限 ≤ " << error
             << ", 真實誤差 = " << actual_error << endl;
    }

    return 0;
}
