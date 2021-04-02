//sqrt求单个欧拉函数
ll ph(ll x) {
    ll ret = x, tmp = x;
    for (ll i = 2; i * i <= x; i++) {
        if (tmp % i == 0) {
            ret = ret / i * (i - 1);
            while (tmp % i == 0) tmp /= i;
        }
    }
    if (tmp > 1) ret = ret / tmp * (tmp - 1);
    return ret;
}


//线性筛欧拉函数
int phi[N];
vector<int> prime;
bool isprime[N];
void init() {
    memset(isprime, 1, sizeof(isprime));
    phi[1] = 1, isprime[1] = false;
    for (int i = 2; i < N; ++i) {
        if (isprime[i]) {
            prime.push_back(i);
            phi[i] = i - 1;
        }
        for (int p : prime) {
            if (1ll * i * p > N) break;
            int now = i * p;
            isprime[now] = false;
            if (i % p == 0) {
                phi[now] = phi[i] * p;
                break;
            } else {
                phi[now] = phi[i] * (p - 1);
            }
        }
    }
}