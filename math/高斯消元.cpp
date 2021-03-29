#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

const int maxn = 10;
const int maxq = 100;

int n;
double a[maxn][maxn + 1]; //ϵ������
double res[maxn] = {}; //�������
double query[maxq];

void solve() {
    int now; // ��ǰ������
    double temp; //���ڼ�¼��Ԫʱ������
    for (int j = 0; j < n; j++) {
        now = j;
        // �ҵ���ǰ�����ϵ��
        for (int i = j; i < n; i++)
            if (fabs(a[i][j]) > fabs(a[now][j]))
                now = i;
        if (now != j)
            //�н���
            for (int i = 0; i < n + 1; i++) {
                double t = a[j][i];
                a[j][i] = a[now][i];
                a[now][i] = t;
            }
        // ��Ԫ
        for (int i = j + 1; i < n; i++) {
            temp = a[i][j] / a[j][j];
            for (int k = j; k <= n; k++)
                a[i][k] -= a[j][k] * temp;
        }
    }
    //��ý��
    for (int i = n - 1; i >= 0; i--) {
        if (a[i][i] == 0) {
            res[i] = 0;
            continue;
        }
        res[i] = a[i][n];
        for (int j = i + 1; j < n; j++) {
            res[i] -= a[i][j] * res[j];
        }
        res[i] /= a[i][i];
    }
}

int main() {
    cin >> n;
    if (n == 1) cin >> a[0][0] >> a[0][1];
    else {
        for (int i = 0; i < n; i++) {
            double x, y;
            cin >> x >> y;
            a[i][n - 1] = 1; //������ϵ��
            a[i][n] = y;
            for (int j = 0; j < n - 1; j++) {
                a[i][j] = pow(x, n - j - 1);
            }
        }
    }
    int q;
    cin >> q;
    for (int i = 0; i < q; i++) cin >> query[i];
    if (n == 1) {
        double k = a[0][1] / a[0][0];
        for (int i = 0; i < q; i++) {
            double ans = query[i] * k;
            if (ans >= -0.005 && ans < 0) ans = 0;
            cout << fixed << setprecision(2) << ans << '\n';
        }
    } else {
        solve();
        for (int i = 0; i < q; i++) {
            double ans = 0;
            for (int j = 0; j < n; j++) {
                ans += pow(query[i], n - j - 1) * res[j];
            }
            if (ans >= -0.005 && ans < 0) ans = 0;
            cout << fixed << setprecision(2) << ans << '\n';
        }
    }
    return 0;
}
