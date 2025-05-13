ll cross(const pair<ll,ll>& O, const pair<ll,ll>& A, const pair<ll,ll>& B) {
    return (A.first - O.first)*(B.second - O.second)- (A.second - O.second)*(B.first - O.first);
}
struct CH {
    vector<pair<ll,ll>>& p;
    CH(vector<pair<ll,ll>>& v): p(v) {}
    vector<int> build() {
        int n = p.size(), k = 0;
        vector<int> ord(n), H(2*n);
        iota(ord.begin(), ord.end(), 0);
        sort(ord.begin(), ord.end(), [&](int i, int j){ return p[i] < p[j]; });
        for (int i : ord) {
            while (k >= 2 && cross(p[H[k-2]], p[H[k-1]], p[i]) <= 0) --k;
            H[k++] = i;
        }
        for (int t = k+1, i = n-2; i >= 0; --i) {
            int idx = ord[i];
            while (k >= t && cross(p[H[k-2]], p[H[k-1]], p[idx]) <= 0) --k;
            H[k++] = idx;
        }
        H.resize(k-1);
        return H;
    }
};
