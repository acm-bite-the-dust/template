void top() {
    for (int i = 1; i <= n; i++) {
        for (int v : G[i]) {
            in[v]++;
        }
    }
    queue<int> q;
    for (int i = 1; i <= n; i++) {
        if (!in[i]) {
            q.push(i);
        }
    }
    int tot = 0;
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        ans[++tot] = u;
        for (int v : G[u]) {
            in[v]--;
            if (!in[v]) {
                q.push(v);
            }
        }
    }
    if (tot != n) {
        printf("NO\n");
    } else {
        printf("YES\n");
        for (int i = 1; i <= n; i++) {
            printf("%d%c", ans[i], i == n ? '\n' : ' ');
        }
    }
}