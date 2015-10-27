#include <iostream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <vector>
#include <ctime>
#include <map>
#include <sys/time.h>

using namespace std;

const int Max_M = 11000;        //最大文档数
const int Max_K = 100;          //最大主题数
const int Max_V = 11000;        //最大单词数
const int Max_Len = 100;        //文档最多包含单词数
const int Line_Len = 11000;     //读入一行最大长度
const int Count = 1000;         //迭代次数
const double alpha = 0.01;      //先验参数
const double beta = 0.1;        //先验参数
int M, V, K = 3;                //分别表示文档数、单词数、主题数
int N[Max_M];                   //N[i]表示第i篇文档的单词数

map<string, int> ID;            //单词到序号的映射
map<int, string> Word;          //序号到单词的映射
vector<int> W[Max_M];           //W[i]为第i篇文档全部单词构成的向量

int n1[Max_M][Max_K], n2[Max_K][Max_V], Sum[Max_K]; //分别表示Nmk. , N.kt 和 sigma N.kt
int Z[Max_M][Max_Len];          //表示Zmn
double theta[Max_M][Max_K], phi[Max_K][Max_V];      //结果参数theta和phi
double likelihood = 0;          //log－likelihood

struct timeval start, end;      //记录程序的开始和结束时间
struct A {                      //topic explanation时对单词进行排序的结构体
    int id, value;
};

bool cmp(const A &a, const A &b) {
    return a.value > b.value;
}

int getId(string s) { //获取单词的序号
    if (ID[s] == 0) {
        V++;
        ID[s] = V;
        Word[V] = s;
    }
    return ID[s];
}

bool valid(char ch) { //判断字符是否属于单词
    if (ch >= 'a' && ch <= 'z') return true;
    if (ch >= 'A' && ch <= 'Z') return true;
    if (ch >= '0' && ch <= '9') return true;
    if (ch == '-') return true;
    return false;
}

void read(char s[], int m) { //从s中读入第m个文档的所有单词
    string word = "";
    W[m].push_back(0);
    for (int p = 0; s[p]; p++) {
        if (!valid(s[p])) {
            if (word != "") {
                N[m]++;
                W[m].push_back(getId(word));
                word = "";
            }
        } else {
            word += s[p];
        }
    }
    if (word != "") {
        N[m]++;
        W[m].push_back(getId(word));
        word = "";
    }
}

void init() { //初始化，指定z[m][n] = 0
    srand(time(NULL));
    memset(n1, 0, sizeof(n1));
    memset(n2, 0, sizeof(n2));
    memset(Sum, 0, sizeof(Sum));
    for (int m = 1; m <= M; m++) {
        for (int n = 1; n <= N[m]; n++) {
            Z[m][n] = 0;
            int t = W[m][n];
            int k = Z[m][n];
            n1[m][k]++;
            n2[k][t]++;
            Sum[k]++;
        }
    }
}

int sampling(int m, int t) { //gibbs sampling
    double p[Max_K], sum = 0;
    for (int k = 1; k <= K; k++) {
        p[k] = (n2[k][t] + beta) * (n1[m][k] + alpha) / (Sum[k] + beta * V);
        sum += p[k];
    }
    for (int k = 1; k <= K; k++) {
        p[k] /= sum;
    }
    double r = rand() * 1.0/ RAND_MAX, add = 0;
    for (int k = 1; k <= K; k++) {
        add += p[k];
        if (r < add) return k;
    }
    return K;
}

void lda_gibbs() { //lda
    for (int count = 1; count <= Count; count++) {
        for (int m = 1; m <= M; m++) {
            for (int n = 1; n <= N[m]; n++) {
                int t = W[m][n];
                int k = Z[m][n];
                n1[m][k]--;
                n2[k][t]--;
                Sum[k]--;
                k = sampling(m, t);
                Z[m][n] = k;
                n1[m][k]++;
                n2[k][t]++;
                Sum[k]++;
            }
        }
    }
}

void cal_parameter() { //计算theta 和 phi
    for (int m = 1; m <= M; m++) {
        int sum = 0;
        for (int k = 1; k <= K; k++)
            sum += n1[m][k];
        for (int k = 1; k <= K; k++)
            theta[m][k] = (n1[m][k] + alpha) / (sum + K * alpha);
    }
    for (int k = 1; k <= K; k++) {
        int sum = 0;
        for (int t = 1; t <= V; t++)
            sum += n2[k][t];
        for (int t = 1; t <= V; t++)
            phi[k][t] = (n2[k][t] + beta) / (sum + V * beta);
    }

}

void cal_likelihood() { //计算相似度
    for (int m = 1; m <= M; m++) {
        for (int n = 1; n <= N[m]; n++) {
            double tmp = 0;
            int t = W[m][n];
            for (int k = 1; k <= K; k++) {
                tmp += theta[m][k] * phi[k][t];
            }
            likelihood += log(tmp);
        }
    }
}

void output() { //输出结果
    gettimeofday(&end, NULL);
    printf("Topic number: %d\n", K);
    printf("log-likelihood: %.6lf\n", likelihood);
    printf("time: %lld microseconds\n",(end.tv_sec - start.tv_sec) * 1000ll + (end.tv_usec-start.tv_usec) / 1000);
}

void topic_explanation() { //主题解释
    vector<A> vec;
    for (int k = 1; k <= K; k++) {
        printf("\nTopic %d:\n", k);
        vec.clear();
        for (int t = 1; t <= V; t++) {
            A tmp;
            tmp.id = t;
            tmp.value = n2[k][t];
            vec.push_back(tmp);
        }
        sort(vec.begin(), vec.end(), cmp);
        for (int i = 0; i < 10; i++) {
            int t = vec[i].id;
            cout<<Word[t];
            double tmp = vec[i].value * 1.0 / Sum[k];
            printf(" : %.10lf\n", tmp);
        }
    }
}

int main() {
    gettimeofday(&start, NULL);
    char s[Line_Len];
    cin>>M;
    gets(s);
    for (int i = 1; i <= M; i++) {
        gets(s);
        read(s, i);
    }
    init();
    lda_gibbs();
    cal_parameter();
    cal_likelihood();
    output();
    topic_explanation();
    return 0;
}