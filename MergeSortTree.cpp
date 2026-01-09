struct MergeSortTree{
    vector<vector<ll>> tree;
    ll n;
    MergeSortTree(ll arr[], ll sz){
        n=sz;
        tree.resize(4*n);
        build(arr,0,0,n-1);
    }
    void build(ll arr[], ll ind, ll st, ll ed){
        if(st==ed){tree[ind]={arr[st]};return;}
        ll mid=(st+ed)/2;
        build(arr,2*ind+1,st,mid);
        build(arr,2*ind+2,mid+1,ed);
        merge(tree[2*ind+1].begin(),tree[2*ind+1].end(),tree[2*ind+2].begin(),tree[2*ind+2].end(),back_inserter(tree[ind]));
    }
    // l theke r er moddhe k er <= koyta
    ll cnt(ll ind,ll l,ll r,ll st,ll ed,ll k){
        if(r<st || l>ed) return 0;
        if(l<=st && r>=ed)return upper_bound(tree[ind].begin(),tree[ind].end(),k)-tree[ind].begin();
        ll mid=(st+ed)/2;
        return cnt(2*ind+1,l,r,st,mid,k)+cnt(2*ind+2,l,r,mid+1,ed,k);
    }
};
