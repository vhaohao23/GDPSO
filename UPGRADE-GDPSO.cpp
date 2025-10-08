#include<bits/stdc++.h>
using namespace std;

#define rep(i,a,b,x)  for(int i=a;i<=b;i+=x)

random_device rd;   
mt19937 gen(rd());

int N=100;
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

void mutation(vector<int> &l,double u){
    uniform_real_distribution<double> dis(0,1);
    
    int S = *max_element(l.begin(), l.end());
    vector<int> ltmp;
    for (int i=1;i<=n;i++){
        double x=dis(gen);
        if (x<u){
            ltmp=l;
            ++S;
            double y=dis(gen);
            if (y<0.5){                
                l[i]=S;
            }
            else{
                l[i]=S;
                for (int neigbor:E[i])
                    l[neigbor]=S;
            }
            
            if (modularity(l)<=modularity(ltmp)){
                --S;
                l=ltmp;
            }
        }
    }
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

void consolidation(vector<int> &l,int l1,int l2){
    for (int i=1;i<=N;i++)
        if (l[i]==l1)
            l[i]=l2;
}

void SecondaryCommunityConsolidation(vector<int>&l){
    map<int,bool> mp;
    map<int,int> cntNode;
    int numC=0;
    for (int i=1;i<=n;i++){
        if (mp.find(l[i])==mp.end()){
            ++numC;
            mp[l[i]]=1;
        }
        cntNode[l[i]]++;
    }

    
    vector<pair<int,int>> decComminities;
    for (auto [x,_]:mp)
        decComminities.push_back({x,cntNode[x]});
    
        sort(decComminities.begin(),decComminities.end(),[](pair<int,int> a,pair<int,int> b){
        return a.second>b.second;
    });

    int i=numC-1;
    vector<int> xtmp;
    while (i>0){
        int j=0;
        bool check=0;
        while (i>j){
            xtmp=l;
            consolidation(l,decComminities[i].first,decComminities[j].first);

            if (modularity(l)>modularity(xtmp)){
                --i;
                check=1;
                break;
            }
            else l=xtmp;
            ++j;
        }

        if (!check) --i;
    }
    
}

void EPD(){
    if (P.size()<10) return;

    vector<pair<double, int>> modularityValues;
    for (int i = 1; i <= N; i++) {
        modularityValues.push_back({Q[i], i});
    }

    sort(modularityValues.begin(), modularityValues.end());

    vector<vector<int>> sortedP(N + 1);
    vector<vector<int>> sortedPb(N + 1);
    vector<double> sortedQ(N + 1);
    vector<double> sortedQb(N + 1);

    for (int i = 0; i < N; i++) {
        sortedP[i + 1] = P[modularityValues[i].second];
        sortedPb[i + 1] = Pb[modularityValues[i].second];
        sortedQ[i + 1] = Q[modularityValues[i].second];
        sortedQb[i + 1] = Qb[modularityValues[i].second];
    }

    P=sortedP,Pb=sortedPb,Q=sortedQ,Qb=sortedQb;

    double N_nor=N-(N/2+1)+1;
    uniform_real_distribution<double> dis(0.0,1.0);
    for (int i=N/2+1;i<=N;i++){
        double C=1.0-exp(-double(i)/N_nor);
        double rand=dis(gen);
        if (rand<=C){
            P.erase(P.begin() + i);
            Pb.erase(Pb.begin() + i);
            Q.erase(Q.begin() + i);
            Qb.erase(Qb.begin() + i);
            
            --N;
            --i;
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
            vector<int> testVb=subtraction(Pb[i],P[i]);
            rep(i,1,n,1)
                if (testVb[i]==1) cout<<"siuuu"; 
            vector<double> V2=multiplicationRcAndV(subtraction(Pg,P[i]),r2*c2);

            vector<double> Vtotal(n+1,0.0);
            rep(j,1,n,1)
                Vtotal[j] = V1[j] + V2[j];


            V[i]=merge(wVi,Vtotal);
            keyOperator(i);
            mutation(P[i],0.3);
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
        EPD();
    }

    SecondaryCommunityConsolidation(Pg);
    Qg=modularity(Pg);
    for (int i=1;i<=N;i++){
        SecondaryCommunityConsolidation(Pb[i]);
        Qb[i]=modularity(Pb[i]);

        SecondaryCommunityConsolidation(P[i]);
        Q[i]=modularity(P[i]);
        if (Q[i]>Qb[i]){
            Qb[i]=Q[i];
            Pb[i]=P[i];
        }

        if (Qb[i]>Qg){
            Qg=Qb[i];
            Pg=Pb[i];
        }
    }

    cout<<Qg<<"\n";
    rep(i,1,n,1)
        cout<<Pg[i]<<" ";
}

int main(){
    clock_t tStart = clock();

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

    printf("\nTime taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
}