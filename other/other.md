# 杂项

### mt19937

```cpp
mt19937 mt(chrono::steady_clock::now().time_since_epoch().count());
ll rng(ll l, ll r) {
    uniform_int_distribution<ll> uni(l, r);
    return uni(mt);
}
```



### dbg宏

```cpp
#define dbg(x...) \
    do { \
        cout << #x << " -> "; \
        err(x); \
    } while (0)

void err() {
    cout << endl;
}

template<class T, class... Ts>
void err(T arg, Ts &... args) {
    cout << arg << ' ';
    err(args...);
}
```



### getchars()快读

```cpp
//quick read
template<typename T>
inline void read(T &x) {
    int s = 1;
    x = 0;
    char ch = getchar();
    while (ch < '0' || ch > '9') {
        if (ch == '-') s = -1;
        ch = getchar();
    }
    while (ch >= '0' && ch <= '9') {
        x = (x << 3) + (x << 1) + (ch ^ 48);
        ch = getchar();
    }
    x *= s;
}

template<typename T, typename... Args>
inline void read(T &x, Args &... args) {
    read(x);
    read(args...);
}
```



### Linux对拍

```
while true; do
./make>tmp.in #出数据
./tmp<tmp.in>tmp.out #被测程序
./tmp2<tmp.in>tmp2.out #正确（暴力）程序
if diff tmp.out tmp2.out; then #比较两个输出文件
printf AC #结果相同显示AC
else
echo WA #结果不同显示WA，并退出
#cat tmp.out tmp2.out
exit 0
fi #if的结束标志,与C语言相反，0为真
done # while的结束标志

#BY NICK WONG 2014-08-29
#在终端下，进入当前目录，输入"sh ./nick.sh",（其中nick.sh为当前shell脚本名） '#'表示单行注释
#diff在两文件相同时返回空串
```



### windows对拍脚本

```
:loop
makedata.exe
K.exe
Kture.exe
fc a.out b.out
if %errorlevel%==0 goto loop
pause 
```

