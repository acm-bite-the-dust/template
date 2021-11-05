# 重载哈希用于unordered_set

```cpp
class my_hash {
public:
    ull operator()(const pair<ll, ll> &p) const {
        return (ull) p.first * P + (ull) p.second;
    }
};

//unorder_set<pair<ll, ll>, my_hash> s;
```

