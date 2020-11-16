#include<iostream>
#include<cstring>
#include<queue>
#include<algorithm>

using namespace std;

const int maxn = 300;
const int maxs = 2e5;

int ans[maxn];

struct ACMachine {
    struct Tree//�ֵ���
    {
        int fail;//ʧ��ָ��
        int vis[26];//�ӽڵ��λ��
        int end;//���������ڵ��β�ĵ��ʱ��
    } AC[maxs * 10];//Trie��
    int cnt = 0;//Trie��ָ��

    inline void init(int x) {
        memset(AC[x].vis, 0, sizeof(AC[x].vis));
        AC[x].fail = 0;
        AC[x].end = 0;
    }

    inline void insert(string s, int id) {
        int l = s.length();
        int now = 0;//�ֵ����ĵ�ǰָ��
        for (int i = 0; i < l; ++i)//����Trie��
        {
            if (AC[now].vis[s[i] - 'a'] == 0)//Trie��û������ӽڵ�
            {
                AC[now].vis[s[i] - 'a'] = ++cnt;//�������
                init(cnt);
            }
            now = AC[now].vis[s[i] - 'a'];//���¹���
        }
        AC[now].end = id;//��ǵ��ʽ�β
    }

    void solve()//����failָ��
    {
        AC[0].fail = 0;//������־
        queue<int> Q;//����
        for (int i = 0; i < 26; ++i)//�ڶ����failָ����ǰ����һ��
        {
            if (AC[0].vis[i] != 0) {
                AC[AC[0].vis[i]].fail = 0;//ָ����ڵ�
                Q.push(AC[0].vis[i]);//ѹ�����
            }
        }
        while (!Q.empty())//BFS��failָ��
        {
            int u = Q.front();
            Q.pop();
            for (int i = 0; i < 26; ++i)//ö�������ӽڵ�
            {
                if (AC[u].vis[i] != 0)//��������ӽڵ�
                {
                    AC[AC[u].vis[i]].fail = AC[AC[u].fail].vis[i];
                    //�ӽڵ��failָ��ָ��ǰ�ڵ��
                    //failָ����ָ��Ľڵ����ͬ�ӽڵ�
                    Q.push(AC[u].vis[i]);//ѹ�����
                } else//����������ӽڵ�
                    AC[u].vis[i] = AC[AC[u].fail].vis[i];
                //��ǰ�ڵ������ӽڵ�ָ��
                //ǰ�ڵ�failָ�������ӽڵ�
            }
        }
    }

    int query(string s)//AC�Զ���ƥ��
    {
        int l = s.length();
        int now = 0, res = 0;
        for (int i = 0; i < l; ++i) {
            now = AC[now].vis[s[i] - 'a'];//����һ��
            for (int t = now; t; t = AC[t].fail)//ѭ�����
                ans[AC[t].end]++;
        }
        return res;
    }
} ACM;

string s[maxn];

int main() {
    ios_base::sync_with_stdio(false);
    int n;
    while (true) {
        scanf("%d", &n);
        if (n == 0) break;
        ACM.cnt = 0;
        ACM.init(0);
        for (int i = 1; i <= n; ++i) {
            s[i].resize(maxs);
            scanf("%s", &s[i][0]);
            s[i].resize(strlen(&s[i][0]));
            ans[i] = 0;
            ACM.insert(s[i], i);
        }
        ACM.solve();//���ʧ��ָ��
        string str;
        str.resize(1e6);
        scanf("%s", &str[0]);//�ı���
        str.resize(strlen(&str[0]));
        ACM.query(str);
        int resm = 0;
        for (int i = 1; i <= n; i++) resm = max(resm, ans[i]);
        printf("%d\n", resm);
        for (int i = 1; i <= n; i++) if (ans[i] == resm) printf("%s\n", s[i].c_str());
    }
    return 0;
}
