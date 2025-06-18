//0 based holee child 2*ind+1,2*ind+2
//1 based holee child 2*ind, 2*ind+1
#include <bits/stdc++.h>
using namespace std;typedef long long int ll;typedef unsigned long long int ull;typedef long double ld;
#define fastio ios_base::sync_with_stdio(false); cin.tie(NULL);
#define all(x) x.begin(),x.end()
#define endl "\n"
#ifndef ONLINE_JUDGE
#define dbg(...) cout<<"["<<#__VA_ARGS__<<"] = ",debugg(__VA_ARGS__)
#else
#define dbg(...)
#endif
void debugg(){cout<<endl;}template<typename Head,typename...Tail>void debugg(Head H,Tail...T){cout<<H; if(sizeof...(T))cout<<", "; debugg(T...);}
ll modExp(ll base,ll exp,ll modVal){base%=modVal;ll result=1;while(exp){if(exp&1)result=(result*base)%modVal;base=(base*base)%modVal;exp>>=1;}return result;}
template<typename T>using minHeap=priority_queue<T, vector<T>, greater<T>>;
struct SegTree{
    vector<array<ll,2>>tree;
    ll n;
    SegTree(ll arr[], ll sz){
        n=sz;
        tree.resize(4*n);
        build(arr,0,0,n-1);
    }
    void build(ll arr[],ll ind,ll st,ll ed){
        if(st==ed){tree[ind][0]=arr[st];tree[ind][1]=1;return;}
        ll mid=(st+ed)/2;
        build(arr,2*ind+1,st,mid);build(arr,2*ind+2,mid+1,ed);
        ll f=tree[2*ind+1][1],mn=tree[2*ind+1][0];
        if(mn>tree[2*ind+2][0])f=tree[2*ind+2][1],mn=tree[2*ind+2][0];else if(mn==tree[2*ind+2][0])f+=tree[2*ind+2][1];
        tree[ind][0]=mn;tree[ind][1]=f;
    }
    array<ll,2> mini(ll ind,ll l, ll r, ll st, ll ed){
        if(r<st or l>ed)return {LLONG_MAX,0};
        if(l==st and r==ed)return {tree[ind][0],tree[ind][1]};
        ll mid=(st+ed)/2;
        array<ll,2>tmp=mini(2*ind+1,l,min(r,mid),st,mid);
        array<ll,2>tmp1=mini(2*ind+2,max(l,mid+1),r,mid+1,ed);
        if(tmp[0]<tmp1[0])return tmp;
        else if(tmp[0]>tmp1[0])return tmp1;
        return {tmp[0],tmp[1]+tmp1[1]};
    }
    void update(ll pos,ll val,ll ind,ll st,ll ed){
        if(pos>ed or pos<st)return;
        if(st==ed){tree[ind][0]=val;tree[ind][1]=1;return;}
        ll mid=(st+ed)/2;
        update(pos,val,2*ind+1,st,mid);update(pos,val,2*ind+2,mid+1,ed);
        ll f=tree[2*ind+1][1],mn=tree[2*ind+1][0];
        if(mn>tree[2*ind+2][0])f=tree[2*ind+2][1],mn=tree[2*ind+2][0];else if(mn==tree[2*ind+2][0])f+=tree[2*ind+2][1];
        tree[ind][0]=mn;tree[ind][1]=f;
    }
};

void solve(){
    ll n,q;cin>>n>>q;
    ll arr[n];for(ll i=0;i<n;i++)cin>>arr[i];
    SegTree st(arr,n);
    while(q--){
        ll x,y,z;cin>>x>>y>>z;
        if(x==1)st.update(y,z,0,0,n-1);
        else cout<<st.mini(0,y,z-1,0,n-1)[0]<<" "<<st.mini(0,y,z-1,0,n-1)[1]<<endl;
    }  
}
 
int main(){
    fastio
    ll t=1; 
    while(t--) solve(); 
    return 0;
}
//https://codeforces.com/edu/course/2/lesson/4/1/practice/contest/273169/problem/C
