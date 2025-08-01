struct LazySegTree{
    vector<ll>tree,lazy;
    ll n;
    LazySegTree(ll arr[], ll sz){
        n=sz;
        tree.resize(4*n);
        lazy.resize(4*n);
        build(arr,0,0,n-1);
    }
    void push(ll ind, ll st, ll ed){
        if(lazy[ind]==0)return;
        tree[ind]+=lazy[ind]*(ed-st+1);
        if(st!=ed){
            lazy[2*ind+1]+=lazy[ind];
            lazy[2*ind+2]+=lazy[ind];
        }
        lazy[ind]=0;
    }
    void build(ll arr[],ll ind,ll st,ll ed){
        if(st==ed){tree[ind]=arr[st];return;}
        ll mid=(st+ed)/2;
        build(arr,2*ind+1,st,mid);build(arr,2*ind+2,mid+1,ed);
        tree[ind]=tree[2*ind+1]+tree[2*ind+2];
    }
    //l theke r er sum
    ll sum(ll ind,ll l, ll r, ll st, ll ed){
        push(ind,st,ed);
        if(r<st or l>ed)return 0;
        if(l<=st and r>=ed)return tree[ind];
        ll mid=(st+ed)/2;
        return sum(2*ind+1,l,min(r,mid),st,mid)+sum(2*ind+2,max(l,mid+1),r,mid+1,ed);
    }
    //l theke r e value add
    void updateRange(ll val,ll ind,ll l,ll r,ll st,ll ed){
        push(ind,st,ed);
        if(r<st or ed<l)return;
        if(l<=st and ed<=r){
            lazy[ind]+=val;
            push(ind,st,ed);return;
        }
        ll mid=(st+ed)/2;
        updateRange(val,2*ind+1,l,min(r,mid),st,mid);
        updateRange(val,2*ind+2,max(l,mid+1),r,mid+1,ed);
        tree[ind]=tree[2*ind+1]+tree[2*ind+2];
    }
};
