struct mat {
    int data[N][N] = {};
    int size{};

    int *operator[](int index) {
        return data[index];
    }
} G, e;

mat operator+(const mat &a, const mat &b) {
    mat ret;
    ret.size = a.size;
    for (int i = 1; i <= a.size; ++i) {
        for (int j = 1; j <= a.size; ++j) {
            ret.data[i][j] = a.data[i][j] + b.data[i][j];
        }
    }
    return ret;
}

mat mul(mat &A, mat &B) {
    mat C;
    C.size = A.size;
    for (int i = 1; i <= A.size; i++) {
        for (int k = 1; k <= A.size; k++) {
            for (int j = 1; j <= A.size; j++) {
                C[i][j] = (C[i][j] + A[i][k] * B[k][j]) % mod;
            }
        }
    }
    return C;
}

mat matpow(mat A, int n) {
    mat B;
    B.size = A.size;
    for (int i = 1; i <= A.size; i++) {
        B[i][i] = 1;
    }
    while (n) {
        if (n & 1) B = mul(B, A);
        A = mul(A, A);
        n >>= 1;
    }
    return B;
}

/*倍增法求解A^1 + A^2 + ... + A^n*/
mat pow_sum(const mat &a, int n) {
    if (n == 1) return a;
    mat tmp = pow_sum(a, n / 2);
    mat tt = matpow(a, n / 2);
    mat sum = tmp + mul(tmp, tt);
    ///若n为奇数，n/2 + n/2 = n-1, 因此sum需要加上A^(n)这一项
    if (n & 1) sum = sum + matpow(a, n);
    return sum;
}