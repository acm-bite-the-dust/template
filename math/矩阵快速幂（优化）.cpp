//矩阵快速幂
#include <vector>
#include <iostream>

using namespace std;
typedef long long ll;
const int maxn = 250;

struct mat {
    ll data[maxn][maxn] = {};
    int size;

    ll *operator[](int index) {
        return data[index];
    }
};

const ll mod = 1e9 + 7;

mat mul(mat &A, mat &B) {     //矩阵乘法
    mat C;
    C.size = A.size;
    for (int i = 0; i < A.size; i++) {
        for (int k = 0; k < A.size; k++) {
            for (int j = 0; j < A.size; j++) {
                C[i][j] = (C[i][j] + A[i][k] * B[k][j]) % mod;
            }
        }
    }
    return C;
}

mat matpow(mat A, ll n) {       //矩阵快速幂
    mat B;
    B.size = A.size;
    for (int i = 0; i < A.size; i++) {
        B[i][i] = 1;
    }
    while (n) {
        if (n & 1) B = mul(B, A);
        A = mul(A, A);
        n >>= 1;
    }
    return B;
}

int main() {
    int n;
    ll k;
    cin >> n >> k;
    mat A;
    A.size = n;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> A[i][j];
        }
    }
    A = matpow(A, k);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
    return 0;
}
