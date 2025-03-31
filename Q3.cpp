#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

struct HermiteData {
    vector<double> z;
    vector<vector<double>> Q;
};

// 建立 Hermite 分差表
HermiteData computeHermiteTable(const vector<double>& T, const vector<double>& D, const vector<double>& V) {
    int n = T.size();
    int m = 2 * n;

    vector<double> z(m);
    vector<vector<double>> Q(m, vector<double>(m, 0.0));

    for (int i = 0; i < n; ++i) {
        z[2 * i] = z[2 * i + 1] = T[i];
        Q[2 * i][0] = Q[2 * i + 1][0] = D[i];
        Q[2 * i + 1][1] = V[i];

        if (i > 0) {
            Q[2 * i][1] = (D[i] - D[i - 1]) / (T[i] - T[i - 1]);
        }
    }

    for (int j = 2; j < m; ++j) {
        for (int i = 0; i < m - j; ++i) {
            Q[i][j] = (Q[i + 1][j - 1] - Q[i][j - 1]) / (z[i + j] - z[i]);
        }
    }

    return {z, Q};
}

// 計算 Hermite 多項式值（牛頓形式）
double evaluateHermite(const vector<double>& z, const vector<vector<double>>& Q, double x) {
    double result = Q[0][0];
    double term = 1.0;

    for (int j = 1; j < z.size(); ++j) {
        term *= (x - z[j - 1]);
        result += Q[0][j] * term;
    }

    return result;
}

// 計算 Hermite 導數值（使用數值微分）
double derivativeAt(const vector<double>& z, const vector<vector<double>>& Q, double x) {
    double h = 1e-5;
    double f_plus = evaluateHermite(z, Q, x + h);
    double f_minus = evaluateHermite(z, Q, x - h);
    return (f_plus - f_minus) / (2 * h);
}

int main() {
    vector<double> T = {0, 3, 5, 8, 13};
    vector<double> D = {0, 200, 375, 620, 990};
    vector<double> V = {75, 77, 80, 74, 72};

    HermiteData data = computeHermiteTable(T, D, V);

    double D_10 = evaluateHermite(data.z, data.Q, 10.0);
    double V_10 = derivativeAt(data.z, data.Q, 10.0);

    cout << fixed << setprecision(3);
    cout << "預測在 t = 10 秒時的位置 D(10): " << D_10 << " 英尺\n";
    cout << "預測在 t = 10 秒時的速度 V(10): " << V_10 << " 英尺/秒\n";

    // 超速檢查
    double speed_limit = 80.67;
    double exceed_time = -1;
    for (int i = 0; i <= 1000; ++i) {
        double t = 0 + (13.0 - 0) * i / 1000;
        double v = derivativeAt(data.z, data.Q, t);
        if (v > speed_limit) {
            exceed_time = t;
            break;
        }
    }

    if (exceed_time > 0)
        cout << "\n汽車在 t ≈ " << exceed_time << " 秒時，第一次超過速度限制 " << speed_limit << " ft/s\n";
    else
        cout << "\n汽車從未超過速度限制 " << speed_limit << " ft/s\n";

    // 最大速度與時間
    double max_speed = -1e9;
    double max_time = -1;

    for (int i = 0; i <= 1000; ++i) {
        double t = 0 + (13.0 - 0) * i / 1000;
        double v = derivativeAt(data.z, data.Q, t);
        if (v > max_speed) {
            max_speed = v;
            max_time = t;
        }
    }

    cout << "\n預測的最大速度 V_max: " << max_speed << " 英尺/秒，發生在 t ≈ " << max_time << " 秒\n";
    cout << "預測的最大速度 V_max: " << max_speed * 0.681818 << " 英里/小時\n";

    return 0;
}
