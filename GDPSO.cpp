#include<bits/stdc++.h>
using namespace std;

#define rep(i,a,b,x)  for(int i=a;i<=b;i+=x)

random_device rd;   
mt19937 gen(rd());

const int N=100;
const double c1=1.4961,c2=1.4961;
const double w=0.7298;

int T=100;

int n,m;
vector<vector<int>> E;
vector<vector<int>> P(N+1);
vector<vector<int>> Pb(N+1);
vector<int> Pg;
vector<double> Q(N+1,0);
vector<double> Qb(N+1,0);
double Qg=0;

vector<vector<double>> V;
vector<double> k;
vector<vector<double>> A;

double modularity(vector<int> l){
    int S=0;
    double Q=0;
    for (int label:l)
        S=max(S,label);
    
    vector<int> cs[S+1];

    for (int i=1;i<=n;i++)
        cs[l[i]].push_back(i);

    for (auto c:cs){
        double sumd=0;
        for (int u:c){
            sumd+=k[u];
            for (int v:c)
                if (u<v)
                    Q+=double(A[u][v])/double(m);
        }
        Q-=pow(sumd/(2*m),2);
    }
    return Q;
}

// double modularity(vector<int> l){
//     double Q=0;

//     for (int i=1;i<=n;i++)
//         for (int j=1;j<=n;j++)
//                 Q+=double(A[i][j]-k[i]*k[j]/double(2*m))*double(l[i]==l[j]);
//     return Q/double(2*m);
// }

void LAR_rand(vector<vector<int>> &a){
    rep(u,1,n,1){
        if (!E[u].size()) continue;
        uniform_int_distribution<int> disv(0,E[u].size()-1);
        int v=E[u][disv(gen)];
        a[u].push_back(v);
        a[v].push_back(u);
    }
}

vector<int> decoding(vector<vector<int>> a){
    bool dd[n+1]={};
    vector<int> l(n+1);
    int cnt=0;

    rep(i,1,n,1)
        if (!dd[i]){
            ++cnt;
            queue<int> q;
            q.push(i);
            while (!q.empty()){
                int u=q.front();
                q.pop();
                l[u]=cnt;
                for (int v:a[u])
                    if (!dd[v]){
                        dd[v]=true;
                        q.push(v);
                    }
            }
        }
    return l;
}

void initialization(){
    rep(i,1,N,1){
        vector<vector<int>> a(n+1);
        LAR_rand(a);
        P[i]=decoding(a);
        Q[i]=modularity(P[i]);
        Pb[i]=P[i];
        Qb[i]=Q[i];
        if (Qb[i]>Qg){
            Qg=Qb[i];
            Pg=Pb[i];
        }
    } 
}
 

void standardization(vector<int> &p){
    map<int,int> mp;
    int cnt=0;
    rep(i,1,n,1)
        if (!mp[p[i]])
            mp[p[i]]=++cnt;
    rep(i,1,n,1)
        p[i]=mp[p[i]];
}

double Div(vector<int> Pbi){
    int res=0;
    rep(k,1,n,1){
        res+=(Pbi[k]!=Pg[k]);
    }
    return double(res)/double(n);
}

vector<int> subtraction(vector<int> p1,vector<int> p2){
    vector<int> resV(n+1,0);
    rep(i,1,n,1)
        resV[i]=(p1[i]!=p2[i]);
    
    return resV;
}

vector<double> multiplicationRcAndV(vector<int> v,double rc){
    vector<double> res(v.begin(),v.end());
    rep(i,1,n,1)
        res[i]*=rc;
    
    return res;
}

vector<double> merge(vector<double> v1,vector<double> v2){
    vector<double> res(n+1,0);
    rep(i,1,n,1)
        res[i]=((v1[i]+v2[i])>=1);
    
    return res;
}

void keyOperator(int t){

    rep(k,1,n,1) {
        
        if (V[t][k]==1) {
            bool dd[n+1]={};
            // Find the community with largest modularity increment among neighbors
            int best_community = P[t][k];
            double best_increment = 0;
            
            // Calculate original modularity
            double orig_mod = modularity(P[t]);
            dd[P[t][k]] = 1;
            for (int neighbor : E[k]) {
                int l_neigh = P[t][neighbor];
                if (dd[l_neigh]) continue;
                dd[l_neigh] = 1;
                
                // Try moving node k to community l_neigh
                int original_community = P[t][k];
                P[t][k] = l_neigh;
                
                // Calculate new modularity
                double new_mod = modularity(P[t]);
                
                // Calculate increment
                double increment = new_mod - orig_mod;
                
                if (increment > best_increment) {
                    best_increment = increment;
                    best_community = l_neigh;
                }
                
                // Restore original community
                P[t][k] = original_community;
            }
            
            // Apply mutation to the best community found
            if (best_community != P[t][k]) {
                P[t][k] = best_community;
            }
        }
    }
}

void GDPSO(){
    initialization();

    //standardization
    rep(i,1,N,1){
            standardization(P[i]);
            standardization(Pb[i]);
    }
    standardization(Pg);
    cout<<Qg<<"\n";

    while (T--){
        //update each particle
        uniform_real_distribution<double> randR(0.0,1);
        rep(i,1,N,1){
            vector<int> Vtmp(V[i].begin(),V[i].end());
            vector<double> wVi=multiplicationRcAndV(Vtmp,w);

            double r1=randR(gen),r2=randR(gen);
            vector<double> V1=multiplicationRcAndV(subtraction(Pb[i],P[i]),r1*c1);  
            vector<double> V2=multiplicationRcAndV(subtraction(Pg,P[i]),r2*c2);

            vector<double> Vtotal(n+1,0.0);
            rep(j,1,n,1)
                Vtotal[j] = V1[j] + V2[j];


            V[i]=merge(wVi,Vtotal);
            keyOperator(i);
        }
        

        // standardization
        rep(i,1,N,1){
            standardization(P[i]);
            standardization(Pb[i]);
        }
        standardization(Pg);

        //calculate fitness and update the personal best and global best positions
        rep(i,1,N,1){
            Q[i]=modularity(P[i]);
            if (Qb[i]<Q[i]){
                Qb[i]=Q[i];
                Pb[i]=P[i];
            }

            if (Qb[i]>Qg)
                Qg=Qb[i],Pg=Pb[i];
        }
    }

    cout<<Qg<<"\n";
    rep(i,1,n,1)
        cout<<Pg[i]<<" ";
}

int main(){
    freopen("input.txt","r",stdin);
    // freopen("output.txt","w",stdout);

    cin>>n>>m;

    E.resize(n+1);
    k.resize(n+1,0);  
    A.resize(n+1,vector<double>(n+1,0));
    V.resize(N+1,vector<double>(n+1,0));

    int u,v;
    rep(i,1,m,1){
        cin>>u>>v;
        E[u].push_back(v);
        E[v].push_back(u);
        A[u][v]=1,A[v][u]=1;
        k[u]++,k[v]++;
    }

    GDPSO();
}