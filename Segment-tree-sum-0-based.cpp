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
    vector<ll>tree;
    ll n;
    SegTree(ll arr[], ll sz){
        n=sz;
        tree.resize(4*n);
        build(arr,0,0,n-1);
    }
    void build(ll arr[],ll ind,ll st,ll ed){
        if(st==ed){tree[ind]=arr[st];return;}
        ll mid=(st+ed)/2;
        build(arr,2*ind+1,st,mid);build(arr,2*ind+2,mid+1,ed);
        tree[ind]=tree[2*ind+1]+tree[2*ind+2];
    }
    //l theke r er sum
    ll sum(ll ind,ll l, ll r, ll st, ll ed){
        if(r<st or l>ed)return 0;
        if(l<=st and r>=ed)return tree[ind];
        ll mid=(st+ed)/2;
        return sum(2*ind+1,l,min(r,mid),st,mid)+sum(2*ind+2,max(l,mid+1),r,mid+1,ed);
    }
    void update(ll pos,ll val,ll ind,ll st,ll ed){
        if(st==ed){tree[ind]=val;return;}
        ll mid=(st+ed)/2;
        if(pos<=mid)update(pos,val,2*ind+1,st,mid);
        else update(pos,val,2*ind+2,mid+1,ed);
        tree[ind]=tree[2*ind+1]+tree[2*ind+2];
    }
};

void solve(){
    ll n,q;cin>>n>>q;
    ll arr[n];for(ll i=0;i<n;i++)cin>>arr[i];
    SegTree st(arr,n);
    while(q--){
        ll x,y,z;cin>>x>>y>>z;
        if(x==1)st.update(y,z,0,0,n-1);
        else cout<<st.sum(0,y,z-1,0,n-1)<<endl;
    }  
}
 
int main(){
    fastio
    ll t=1; 
    while(t--) solve(); 
    return 0;
}
//https://codeforces.com/edu/course/2/lesson/4/1/practice/contest/273169/problem/A
