## **问题**

有 ![[公式]](https://www.zhihu.com/equation?tex=n) 个有序对 ![[公式]](https://www.zhihu.com/equation?tex=%5Cleft%3Ca_i%2Cb_i%5Cright%3E) ，对于每个有序对，可以选择将 ![[公式]](https://www.zhihu.com/equation?tex=a_i) 加入集合 ![[公式]](https://www.zhihu.com/equation?tex=A) 或将 ![[公式]](https://www.zhihu.com/equation?tex=b_i) 加入集合 ![[公式]](https://www.zhihu.com/equation?tex=B) ，试最小化 ![[公式]](https://www.zhihu.com/equation?tex=%5Cmax%28A%29%2B%5Cmax%28B%29) 。（规定 ![[公式]](https://www.zhihu.com/equation?tex=%5Cmax%28%5Cvarnothing%29%3D0) ）

## **算法**

首先注意到对于 ![[公式]](https://www.zhihu.com/equation?tex=%5Cleft%3Ca_i%2Cb_i%5Cright%3E) 和 ![[公式]](https://www.zhihu.com/equation?tex=%5Cleft%3Ca_j%2Cb_j%5Cright%3E) ，如果在 ![[公式]](https://www.zhihu.com/equation?tex=a_i%5Cle+a_j) 的同时有 ![[公式]](https://www.zhihu.com/equation?tex=b_i%5Cle+b_j) ，那么就可以舍弃 ![[公式]](https://www.zhihu.com/equation?tex=%5Cleft%3Ca_i%2Cb_i%5Cright%3E) 。于是可以将有序对以 ![[公式]](https://www.zhihu.com/equation?tex=a_i) 为第一关键词、 ![[公式]](https://www.zhihu.com/equation?tex=b_i) 为第二关键词排序，这样可以保证 ![[公式]](https://www.zhihu.com/equation?tex=a_i) 是随 ![[公式]](https://www.zhihu.com/equation?tex=i) 单调不减的：

![[公式]](https://www.zhihu.com/equation?tex=%5Cbegin%7Bmatrix%7D+%5Cleft%3C1%2C4%5Cright%3E%5C%5C+%5Cleft%3C1%2C8%5Cright%3E%5C%5C+%5Cleft%3C2%2C3%5Cright%3E%5C%5C+%5Cleft%3C2%2C7%5Cright%3E%5C%5C+%5Cleft%3C3%2C1%5Cright%3E%5C%5C+%5Cleft%3C4%2C2%5Cright%3E+%5Cend%7Bmatrix%7D) 

显然， ![[公式]](https://www.zhihu.com/equation?tex=b_i) 应该随 ![[公式]](https://www.zhihu.com/equation?tex=i) 单调减。如何保证呢？注意到最后一个有序对肯定要保留，所以我们只需要从后往前遍历，划掉破坏 ![[公式]](https://www.zhihu.com/equation?tex=b_i) 单调性的元素即可：

![[公式]](https://www.zhihu.com/equation?tex=%5Cbegin%7Bmatrix%7D+%5Ccancel%7B%5Cleft%3C1%2C4%5Cright%3E%7D%5C%5C+%5Cleft%3C1%2C8%5Cright%3E%5C%5C+%5Ccancel%7B%5Cleft%3C2%2C3%5Cright%3E%7D%5C%5C+%5Cleft%3C2%2C7%5Cright%3E%5C%5C+%5Ccancel%7B%5Cleft%3C3%2C1%5Cright%3E%7D%5C%5C+%5Cleft%3C4%2C2%5Cright%3E+%5Cend%7Bmatrix%7D) 

当然，在编程中，删除数组中的元素还是太浪费时间了，我们转而把无需删除的有序对加入一个新的数组。只不过，因为我们是从后往前遍历，这样会导致顺序颠倒过来。

![[公式]](https://www.zhihu.com/equation?tex=%5Cbegin%7Bmatrix%7D%5Cleft%3C4%2C2%5Cright%3E%5C%5C%5Cleft%3C2%2C7%5Cright%3E%5C%5C%5Cleft%3C1%2C8%5Cright%3E%5Cend%7Bmatrix%7D) 

这时我们发现，如果我们选了某个 ![[公式]](https://www.zhihu.com/equation?tex=a_i) ，那么可以无代价地选择它**下方**的所有 ![[公式]](https://www.zhihu.com/equation?tex=a_i) ；如果我们选了某个 ![[公式]](https://www.zhihu.com/equation?tex=b_i) ，那么可以无代价地选择它**上方**的所有 ![[公式]](https://www.zhihu.com/equation?tex=b_i) 。相当于，我们只需要选择左边的**开头**和右边的**结尾**，而且它们应该是**相邻**的：

![img](https://pic4.zhimg.com/80/v2-8eefe70d19257a2a16588820e318fcb3_720w.jpg)

这样一来，我们一共就只有 ![[公式]](https://www.zhihu.com/equation?tex=m%2B1) 种选法了（ ![[公式]](https://www.zhihu.com/equation?tex=m) 为新数组的长度），把原来的指数级问题降成了线性。

## 代码

```cpp
int minsum(vector<pair<int, int>> &V)
{
    if (V.empty())
        return 0;
    sort(V.begin(), V.end());
    vector<pair<int, int>> U{V[V.size() - 1]};
    for (int i = V.size() - 2; i >= 0; --i)
        if (V[i].second > U[U.size() - 1].second)
            U.push_back(V[i]);
    int mi = min(U[0].first, U[U.size() - 1].second);
    for (int i = 0; i < U.size() - 1; ++i)
        mi = min(mi, U[i].second + U[i + 1].first);
    return mi;
}
```

## 例题

[CF1408D](https://link.zhihu.com/?target=https%3A//codeforces.com/contest/1408/problem/D)：通过计算，可求得每个强盗向上躲过探照灯需走 ![[公式]](https://www.zhihu.com/equation?tex=u_i) 格，向右躲过探照灯需走 ![[公式]](https://www.zhihu.com/equation?tex=r_i) 格，记为 ![[公式]](https://www.zhihu.com/equation?tex=%5Cleft%3Cu_i%2Cr_i%5Cright%3E) 。现在要让所有强盗向上、向右各走若干格躲过所有探照灯，也即选若干个 ![[公式]](https://www.zhihu.com/equation?tex=u_i) 求最大值，选剩余的 ![[公式]](https://www.zhihu.com/equation?tex=r_i) 求最大值，将这两个最大值的和最小化。这完美符合刚刚那个模型*（当然，因为这模型就是从这个题来的-w-）*。

[CF1435C](https://link.zhihu.com/?target=https%3A//codeforces.com/contest/1435/problem/C)：我们可以得到 ![[公式]](https://www.zhihu.com/equation?tex=n) 组 ![[公式]](https://www.zhihu.com/equation?tex=6) 个数，每组各选择一个，让选择的数的极差最小。可以固定其中一组，对于这组数中的每一个 ![[公式]](https://www.zhihu.com/equation?tex=v) ，都把其余各组中最大的 ![[公式]](https://www.zhihu.com/equation?tex=%5Cle+v) 的数和最小的 ![[公式]](https://www.zhihu.com/equation?tex=%5Cge+v) 的数找出来，分别求它们与 ![[公式]](https://www.zhihu.com/equation?tex=v) 的差（如果没有满足某个条件的数，记对应的差为`INF`）。这相当于是**向左扩展**和**向右扩展**。然后套模型。这样我们算出了选固定组每一个数时的最小极差，求一个最小值即可。