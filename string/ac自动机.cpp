#include<iostream>
#include<cstring>
#include<queue>
#include<algorithm>

using namespace std;

const int maxn = 1e5;
const int maxs = 1e5;

struct Result {
    int num;
    int pos;
} ans[maxn];//���е��ʵĳ��ִ���

bool operator<(Result a, Result b) {
    if (a.num != b.num)
        return a.num > b.num;
    else
        return a.pos < b.pos;
}

struct ACMachine {
    struct Tree//�ֵ���
    {
        int fail;//ʧ��ָ��
        int vis[26];//�ӽڵ��λ��
        int end;//���������ڵ��β�ĵ��ʱ��
    } AC[maxs];//Trie��
    int cnt = 0;//Trie��ָ��

    inline void Clean(int x) {
        memset(AC[x].vis, 0, sizeof(AC[x].vis));
        AC[x].fail = 0;
        AC[x].end = 0;
    }

    inline void Build(string s, int Num) {
        int l = s.length();
        int now = 0;//�ֵ����ĵ�ǰָ��
        for (int i = 0; i < l; ++i)//����Trie��
        {
            if (AC[now].vis[s[i] - 'a'] == 0)//Trie��û������ӽڵ�
            {
                AC[now].vis[s[i] - 'a'] = ++cnt;//�������
                Clean(cnt);
            }
            now = AC[now].vis[s[i] - 'a'];//���¹���
        }
        AC[now].end = Num;//��ǵ��ʽ�β
    }

    void Get_fail()//����failָ��
    {
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
                //��ǰ�ڵ������ӽڵ�ָ��ǰ�ڵ�failָ�������ӽڵ�
            }
        }
    }

    int AC_Query(string s)//AC�Զ���ƥ��
    {
        int l = s.length();
        int now = 0, ans = 0;
        for (int i = 0; i < l; ++i) {
            now = AC[now].vis[s[i] - 'a'];//����һ��
            for (int t = now; t; t = AC[t].fail)//ѭ�����
                ans[AC[t].end].num++;
        }
        return ans;
    }
} ACM;

int main() {
    string s[300];
    int n;
    while (true) {
        cin >> n;
        if (n == 0)break;
        ACM.cnt = 0;
        ACM.Clean(0);
        for (int i = 1; i <= n; ++i) {
            cin >> s[i];
            ans[i].num = 0;
            ans[i].pos = i;
            ACM.Build(s[i], i);
        }
        ACM.AC[0].fail = 0;//������־
        ACM.Get_fail();//���ʧ��ָ��
        cin >> s[0];//�ı���
        ACM.AC_Query(s[0]);
        sort(&ans[1], &ans[n + 1]);
        cout << ans[1].num << endl;
        cout << s[ans[1].pos] << endl;
        for (int i = 2; i <= n; ++i) {
            if (ans[i].num == ans[i - 1].num)
                cout << s[ans[i].pos] << endl;
            else
                break;
        }
    }
    return 0;
}
