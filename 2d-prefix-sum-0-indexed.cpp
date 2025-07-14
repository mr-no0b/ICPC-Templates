ll pre[501][501];
string arr[501];
ll n,m;
void build(){
    for(ll i=0;i<n;i++){
        for(ll j=0;j<m;j++){
            pre[i][j]=(arr[i][j]=='g');//jetar upre prefix calabo
            pre[i][j]+=(i?pre[i-1][j]:0)+(j?pre[i][j-1]:0)-(i&&j?pre[i-1][j-1]:0);
        }
    }
}
ll sum(ll x1,ll y1,ll x2,ll y2){
    return pre[x2][y2]-(x1?pre[x1-1][y2]:0)-(y1?pre[x2][y1-1]:0)+(x1&&y1?pre[x1-1][y1-1]:0);
}
