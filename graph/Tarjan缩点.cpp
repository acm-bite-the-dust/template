const int N = 1e5 + 10;
vector<int> G[N];
int dfn[N], low[N];//dfn[u] -> u�������Ĵ��� low[u] -> u������dfn��Сֵ(������������)
int index;//�������
stack<int> S;
bool ins[N];//�Ƿ��ջ
int col[N], num_color;//Ⱦɫ

void Tarjan(int u) {
    dfn[u] = low[u] = ++index;
    S.push(u);//��ջ
    ins[u] = true;
    for (int i :G[u]) {
        int v = i;
        if (!dfn[v]) { //δ�����ʹ�
            Tarjan(v);
            low[u] = min(low[u], low[v]); //�Ұְ֣�����ͷ����С��
        } else if (ins[v]) { //�ѱ����ʹ�����ջ��,����Ҫ����; ������ջ��˵����Ӧǿ��ͨ�����������
            low[u] = min(low[u], dfn[v]); //�ж�˭�ǰְ�
        }
    }
    if (dfn[u] == low[u]) { //���ָ�����һ���Լ��ǰְ�(ĳǿ��ͨ�����н���һ�������ʵĽ������dfn[u] == low[u])
        num_color++;
        int tmp;
        do {
            tmp = S.top();
            col[tmp] = num_color; //��ջ��Ⱦɫ
            ins[tmp] = false;
            S.pop();
        } while (tmp != u);
        col[u] = num_color;
        ins[u] = false;
        S.pop();
    }
}

/*
 * Ⱦ��ɫ��, ͬһ��color��������һ����, ����ͼ�൱��һ��DAG
 * DAG��������, ����ɶ��, Ϊ����Ϊ
 */