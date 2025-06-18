//Multiply returns ans[k]=∑ ai*bj [such that i+j==k]
//Frequency arrays are often used
const double PI=acos(-1);
struct base{
    double a,b;
    base(double a=0,double b=0):a(a),b(b){}
    base operator+(const base&c)const{return base(a+c.a,b+c.b);}
    base operator-(const base&c)const{return base(a-c.a,b-c.b);}
    base operator*(const base&c)const{return base(a*c.a-b*c.b,a*c.b+b*c.a);}
};
void fft(vector<base>&p,bool inv=false){
    int n=p.size(),i=0;
    for(int j=1;j<n-1;++j){
        for(int k=n>>1;k>(i^=k);k>>=1);
        if(j<i)swap(p[i],p[j]);
    }
    for(int len=1;len<<1<=n;len<<=1){
        int m=len<<1;double ang=2*PI/m;
        base wn=base(cos(ang),(inv?1:-1)*sin(ang)),w;
        for(int i=0;i<n;i+=m){
            w=base(1,0);
            for(int j=i;j<i+len;++j){
                base t=w*p[j+len];
                p[j+len]=p[j]-t;
                p[j]=p[j]+t;
                w=w*wn;
            }
        }
    }
    if(inv)for(int i=0;i<n;++i){p[i].a/=n;p[i].b/=n;}
}
vector<ll> multiply(const vector<int>&a,const vector<int>&b){
    int n=a.size(),m=b.size(),t=n+m-1,sz=1;
    while(sz<t)sz<<=1;
    vector<base>x(sz),y(sz),z(sz);
    for(int i=0;i<sz;++i){
        x[i]=i<n?base(a[i],0):base();
        y[i]=i<m?base(b[i],0):base();
    }
    fft(x);fft(y);
    for(int i=0;i<sz;++i)z[i]=x[i]*y[i];
    fft(z,true);
    vector<ll>r(sz);
    for(int i=0;i<sz;++i)r[i]=llround(z[i].a);
    while(r.size()>1&&r.back()==0)r.pop_back();
    return r;
}
