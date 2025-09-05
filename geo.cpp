#include <bits/stdc++.h>
using namespace std;
using ld = long double;
const ld EPS = 1e-12L;
const ld PI = acosl(-1.0L);

inline int sgn(ld x){ if (x > EPS) return 1; if (x < -EPS) return -1; return 0; }
inline ld sqr(ld x){ return x*x; }

struct Point {
    ld x, y;
    Point(): x(0), y(0) {}
    Point(ld _x, ld _y): x(_x), y(_y) {}
};

static inline ostream& operator<<(ostream &os, const Point &p){
    os << fixed << setprecision(12) << (double)p.x << " " << (double)p.y;
    return os;
}

Point operator+(const Point &a, const Point &b){ return Point(a.x + b.x, a.y + b.y); }
Point operator-(const Point &a, const Point &b){ return Point(a.x - b.x, a.y - b.y); }
Point operator*(const Point &a, ld k){ return Point(a.x * k, a.y * k); }
Point operator*(ld k, const Point &a){ return Point(a.x * k, a.y * k); }
Point operator/(const Point &a, ld k){ return Point(a.x / k, a.y / k); }

ld dot(const Point &a, const Point &b){ return a.x*b.x + a.y*b.y; }
ld cross(const Point &a, const Point &b){ return a.x*b.y - a.y*b.x; }
ld cross(const Point &o, const Point &a, const Point &b){ return cross(a - o, b - o); }
ld norm2(const Point &a){ return dot(a,a); }
ld absval(const Point &a){ return sqrtl(max((ld)0.0L, norm2(a))); }

ld dist(const Point &a, const Point &b){ return absval(a - b); }

// angle and normalization
ld angle(const Point &p){ return atan2l(p.y, p.x); }
ld normalizeAngle(ld a){ while (a <= -PI) a += 2*PI; while (a > PI) a -= 2*PI; return a; }

bool pointEqual(const Point &a, const Point &b){ return sgn(a.x - b.x) == 0 && sgn(a.y - b.y) == 0; }

// Lexicographic sort
bool pointCmp(const Point &a, const Point &b){ if (sgn(a.x - b.x) != 0) return a.x < b.x; return a.y < b.y; }

Point rotate(const Point &p, ld ang){ ld c = cosl(ang), s = sinl(ang); return Point(p.x*c - p.y*s, p.x*s + p.y*c); }

// Line in parametric form p + t * dir
struct Line {
    Point p, dir;
    Line(){}
    Line(Point _p, Point _dir): p(_p), dir(_dir){}
};

// Segment representation
struct Segment { Point a, b; Segment(){} Segment(Point _a, Point _b): a(_a), b(_b){} };

// Ray representation
struct Ray { Point p, dir; Ray(){} Ray(Point _p, Point _dir): p(_p), dir(_dir){} };

// Projection of point p onto line ab (infinite)
Point projectOnLine(const Point &a, const Point &b, const Point &p){
    Point ab = b - a; ld t = dot(p - a, ab) / norm2(ab); return a + ab * t;
}

// Distance point to line (infinite)
// FIX: use fabsl for scalar absolute of cross(...)
ld distPointLine(const Point &a, const Point &b, const Point &p){ return fabsl(cross(b - a, p - a)) / absval(b - a); }

// Distance point to segment
ld distPointSegment(const Point &a, const Point &b, const Point &p){
    Point ab = b - a; ld r = dot(p - a, ab) / norm2(ab);
    if (r < 0) return absval(p - a);
    if (r > 1) return absval(p - b);
    Point proj = a + ab * r; return absval(p - proj);
}

// Check if point p is on segment ab
bool onSegment(const Point &a, const Point &b, const Point &p){
    if (sgn(cross(a,b,p)) != 0) return false;
    return sgn(dot(p - a, p - b)) <= 0;
}

// Segment intersection test (proper or not)
bool segmentsIntersect(const Point &a, const Point &b, const Point &c, const Point &d){
    ld c1 = cross(a,b,c), c2 = cross(a,b,d), c3 = cross(c,d,a), c4 = cross(c,d,b);
    int s1 = sgn(c1), s2 = sgn(c2), s3 = sgn(c3), s4 = sgn(c4);
    if (s1 == 0 && onSegment(a,b,c)) return true;
    if (s2 == 0 && onSegment(a,b,d)) return true;
    if (s3 == 0 && onSegment(c,d,a)) return true;
    if (s4 == 0 && onSegment(c,d,b)) return true;
    return (s1 * s2 < 0) && (s3 * s4 < 0);
}

// Line-line intersection (infinite). returns pair(hasIntersection, point)
pair<bool, Point> lineLineIntersection(const Point &p1, const Point &p2, const Point &q1, const Point &q2){
    Point r = p2 - p1, s = q2 - q1;
    ld rxs = cross(r, s);
    if (sgn(rxs) == 0) return {false, Point()};
    ld t = cross(q1 - p1, s) / rxs;
    return {true, p1 + r * t};
}

// Distance between segments
ld distSegments(const Point &a, const Point &b, const Point &c, const Point &d){
    if (segmentsIntersect(a,b,c,d)) return 0;
    ld d1 = distPointSegment(a,b,c);
    ld d2 = distPointSegment(a,b,d);
    ld d3 = distPointSegment(c,d,a);
    ld d4 = distPointSegment(c,d,b);
    return min(min(d1,d2), min(d3,d4));
}

// Polygon utilities
ld polygonArea(const vector<Point> &P){
    int n = P.size(); if (n < 3) return 0;
    ld a = 0; for (int i = 0; i < n; ++i) a += cross(P[i], P[(i+1)%n]); return a/2.0L;
}

ld polygonPerimeter(const vector<Point> &P){ ld s = 0; for (int i = 0; i < (int)P.size(); ++i) s += dist(P[i], P[(i+1)%P.size()]); return s; }

Point polygonCentroid(const vector<Point> &P){
    ld A = polygonArea(P);
    if (sgn(A) == 0) return Point(0,0);
    ld cx = 0, cy = 0; int n = P.size();
    for (int i = 0; i < n; ++i){
        int j = (i+1)%n; ld cr = cross(P[i], P[j]);
        cx += (P[i].x + P[j].x) * cr;
        cy += (P[i].y + P[j].y) * cr;
    }
    cx /= (6*A); cy /= (6*A); return Point(cx, cy);
}

bool isConvexPolygon(const vector<Point> &P){
    int n = P.size(); if (n < 3) return false;
    int sign = 0;
    for (int i = 0; i < n; ++i){
        Point a = P[i], b = P[(i+1)%n], c = P[(i+2)%n];
        int s = sgn(cross(a,b,c));
        if (s == 0) continue;
        if (sign == 0) sign = s; else if (sign != s) return false;
    }
    return true;
}

// Monotone chain convex hull (returns CCW, no duplicate last element)
vector<Point> convexHull(vector<Point> pts){
    sort(pts.begin(), pts.end(), pointCmp);
    pts.erase(unique(pts.begin(), pts.end(), [](const Point &a, const Point &b){ return pointEqual(a,b); }), pts.end());
    int n = pts.size(); if (n <= 1) return pts;
    vector<Point> lo, hi;
    for (int i = 0; i < n; ++i){ while (lo.size() >= 2 && sgn(cross(lo[lo.size()-2], lo.back(), pts[i])) <= 0) lo.pop_back(); lo.push_back(pts[i]); }
    for (int i = n-1; i >= 0; --i){ while (hi.size() >= 2 && sgn(cross(hi[hi.size()-2], hi.back(), pts[i])) <= 0) hi.pop_back(); hi.push_back(pts[i]); }
    lo.pop_back(); hi.pop_back(); lo.insert(lo.end(), hi.begin(), hi.end());
    return lo;
}

// Sutherland–Hodgman polygon clipping (clip subject by a single directed line a->b) - keeps left side
vector<Point> polygonCut(const vector<Point> &poly, const Point &a, const Point &b){
    vector<Point> res; int n = poly.size();
    for (int i = 0; i < n; ++i){
        Point cur = poly[i], prev = poly[(i + n - 1) % n];
        ld curSide = cross(a,b,cur), prevSide = cross(a,b,prev);
        if (sgn(curSide) >= 0){
            if (sgn(prevSide) < 0){ auto inter = lineLineIntersection(prev, cur, a, b); if (inter.first) res.push_back(inter.second); }
            res.push_back(cur);
        } else if (sgn(prevSide) >= 0){ auto inter = lineLineIntersection(prev, cur, a, b); if (inter.first) res.push_back(inter.second); }
    }
    return res;
}

// Convex polygon intersection (O(n+m)), inputs must be counter-clockwise convex polygons. Returns intersection polygon CCW.
vector<Point> convexPolyIntersection(vector<Point> A, vector<Point> B){
    for (int i = 0; i < (int)B.size(); ++i){
        Point b1 = B[i], b2 = B[(i+1)%B.size()];
        A = polygonCut(A, b1, b2);
        if (A.empty()) return {};
    }
    return A;
}

// Point in polygon: 0 outside, 1 on boundary, 2 inside (ray casting)
int pointInPolygon(const vector<Point> &poly, const Point &p){
    bool in = false; int n = poly.size();
    for (int i = 0, j = n-1; i < n; j = i++){
        const Point &a = poly[i], &b = poly[j];
        if (onSegment(a,b,p)) return 1;
        bool intersect = ((a.y > p.y) != (b.y > p.y)) && (p.x < (b.x - a.x) * (p.y - a.y) / (b.y - a.y) + a.x);
        if (intersect) in = !in;
    }
    return in ? 2 : 0;
}

// Rotating calipers: diameter (max pairwise distance) for convex polygon
pair<ld, pair<Point,Point>> convexDiameter(const vector<Point> &P){
    int n = P.size(); if (n == 0) return {0, {Point(),Point()}}; if (n == 1) return {0, {P[0], P[0]}};
    int j = 1; ld best = 0; pair<Point,Point> bp = {P[0], P[0]};
    for (int i = 0; i < n; ++i){
        int ni = (i+1)%n;
        while (abs(cross(P[i], P[ni], P[(j+1)%n])) > abs(cross(P[i], P[ni], P[j]))) j = (j+1)%n;
        ld d1 = dist(P[i], P[j]); if (d1 > best){ best = d1; bp = {P[i], P[j]}; }
        ld d2 = dist(P[ni], P[j]); if (d2 > best){ best = d2; bp = {P[ni], P[j]}; }
    }
    return {best, bp};
}

// Minimum width of convex polygon
ld convexMinWidth(const vector<Point> &P){
    int n = P.size(); if (n < 3) return 0;
    int k = 1; ld best = 1e300L;
    for (int i = 0; i < n; ++i){
        int ni = (i+1)%n;
        while (abs(cross(P[i], P[ni], P[(k+1)%n])) > abs(cross(P[i], P[ni], P[k]))) k = (k+1)%n;
        ld width = abs(cross(P[ni] - P[i], P[k] - P[i])) / absval(P[ni] - P[i]);
        best = min(best, width);
    }
    return best;
}

// Minimum Enclosing Circle (Welzl) - randomized expected linear time
struct Circle { Point c; ld r; Circle(): c(), r(0) {} Circle(Point _c, ld _r): c(_c), r(_r) {} };

Circle circleFrom2(const Point &a, const Point &b){ Point c = (a + b) / 2.0L; ld r = dist(a,b)/2.0L; return Circle(c,r); }

// forward declaration to fix usage before definition
Point circumcenter(const Point &a, const Point &b, const Point &c);

Circle circleFrom3(const Point &a, const Point &b, const Point &c){
    Point o = circumcenter(a,b,c); if (!isfinite((double)o.x)) return Circle(Point(0,0), -1); ld R = dist(o,a); return Circle(o,R);
}

Point circumcenter(const Point &a, const Point &b, const Point &c){
    ld d = 2.0L * cross(a,b,c);
    if (sgn(d) == 0) return Point(nan(""), nan(""));
    ld ax = a.x, ay = a.y, bx = b.x, by = b.y, cx = c.x, cy = c.y;
    ld A = ax*ax + ay*ay;
    ld B = bx*bx + by*by;
    ld C = cx*cx + cy*cy;
    ld ux = (A*(by - cy) + B*(cy - ay) + C*(ay - by)) / d;
    ld uy = (A*(cx - bx) + B*(ax - cx) + C*(bx - ax)) / d;
    return Point(ux, uy);
}

bool circleContains(const Circle &C, const Point &p){ return sgn(dist(C.c, p) - C.r) <= 0; }

Circle mec(vector<Point> pts){
    random_device rd; mt19937 gen(rd()); shuffle(pts.begin(), pts.end(), gen);
    Circle C = Circle(Point(0,0), -1);
    int n = pts.size();
    for (int i = 0; i < n; ++i){
        if (C.r < 0 || !circleContains(C, pts[i])){
            C = Circle(pts[i], 0);
            for (int j = 0; j < i; ++j) if (!circleContains(C, pts[j])){
                C = circleFrom2(pts[i], pts[j]);
                for (int k = 0; k < j; ++k) if (!circleContains(C, pts[k])){
                    C = circleFrom3(pts[i], pts[j], pts[k]);
                }
            }
        }
    }
    return C;
}

// Half-plane intersection (HPI) - lines given as directed lines (keep left side). Returns convex polygon CCW or empty if no intersection.
struct HLine { Point p; Point dir; ld ang; HLine(){} HLine(Point _p, Point _dir): p(_p), dir(_dir){ ang = atan2l(dir.y, dir.x); } };

bool hlineOnLeft(const HLine &L, const Point &q){ return sgn(cross(L.p, L.p + L.dir, q)) >= 0; }

Point intersectH(const HLine &a, const HLine &b){ auto pr = lineLineIntersection(a.p, a.p + a.dir, b.p, b.p + b.dir); return pr.second; }

vector<Point> halfPlaneIntersection(vector<HLine> &lines){
    sort(lines.begin(), lines.end(), [](const HLine &A, const HLine &B){ if (sgn(A.ang - B.ang) == 0) return cross(A.p, A.p + A.dir, B.p) < 0; return A.ang < B.ang; });
    deque<HLine> dq;
    deque<Point> pts;
    for (int i = 0; i < (int)lines.size(); ++i){
        if (i > 0 && sgn(lines[i].ang - lines[i-1].ang) == 0) continue; // parallel and later
        while (dq.size() >= 2 && !hlineOnLeft(lines[i], pts.back())){ dq.pop_back(); pts.pop_back(); }
        while (dq.size() >= 2 && !hlineOnLeft(lines[i], pts.front())){ dq.pop_front(); pts.pop_front(); }
        if (!dq.empty()){
            Point ip = intersectH(dq.back(), lines[i]); pts.push_back(ip);
        }
        dq.push_back(lines[i]);
    }
    while (dq.size() > 2 && !hlineOnLeft(dq.front(), pts.back())){ dq.pop_back(); pts.pop_back(); }
    while (dq.size() > 2 && !hlineOnLeft(dq.back(), pts.front())){ dq.pop_front(); pts.pop_front(); }
    if (dq.size() < 3) return {};
    vector<Point> res; res.reserve(dq.size());
    for (int i = 0; i < (int)dq.size(); ++i){ int j = (i+1)%dq.size(); res.push_back(intersectH(dq[i], dq[j])); }
    return res;
}

// Minkowski sum of two convex polygons (CCW, no duplicate last). Returns convex polygon CCW.
vector<Point> minkowskiSum(const vector<Point> &A, const vector<Point> &B){
    if (A.empty() || B.empty()) return {};
    auto idx = [](const vector<Point> &P){ int n = P.size(); int id = 0; for (int i = 1; i < n; ++i) if (pointCmp(P[i], P[id])) id = i; return id; };
    int ia = idx(A), ib = idx(B);
    int na = A.size(), nb = B.size();
    vector<Point> C; C.reserve(na + nb);
    int i = ia, j = ib;
    do {
        C.push_back(A[i] + B[j]);
        int ni = (i+1)%na; int nj = (j+1)%nb;
        Point va = A[ni] - A[i]; Point vb = B[nj] - B[j];
        ld cr = cross(va, vb);
        if (sgn(cr) >= 0) i = ni; if (sgn(cr) <= 0) j = nj;
    } while (i != ia || j != ib);
    auto hull = convexHull(C);
    return hull;
}

// Circle-line intersection
vector<Point> circleLineIntersection(const Circle &C, const Point &a, const Point &b){
    Point proj = projectOnLine(a,b,C.c);
    ld d = absval(proj - C.c);
    if (sgn(d - C.r) > 0) return {};
    if (sgn(fabsl(d - C.r)) == 0) return {proj};
    Point dir = (b - a);
    ld len = absval(dir);
    dir = dir / len;
    ld h = sqrtl(max((ld)0.0L, C.r*C.r - d*d));
    return { proj + dir * h, proj - dir * h };
}

// Circle-circle intersection
vector<Point> circleCircleIntersection(const Circle &A, const Circle &B){
    ld d = absval(B.c - A.c);
    if (sgn(d) == 0) return {};
    if (sgn(d - (A.r + B.r)) > 0) return {};
    if (sgn(d - fabsl(A.r - B.r)) < 0) return {};
    ld x = (d*d - B.r*B.r + A.r*A.r) / (2*d);
    ld h2 = A.r*A.r - x*x;
    Point v = (B.c - A.c) / d;
    Point p = A.c + v * x;
    if (sgn(h2) == 0) return {p};
    ld h = sqrtl(max((ld)0.0L, h2)); Point perp = Point(-v.y, v.x);
    return { p + perp * h, p - perp * h };
}

// Tangents from point to circle (returns tangent points on circle)
vector<Point> tangentsFromPointToCircle(const Circle &C, const Point &p){
    ld d2 = norm2(p - C.c);
    if (sgn(d2 - C.r*C.r) < 0) return {};
    if (sgn(d2 - C.r*C.r) == 0) return {p};
    ld d = sqrtl(d2); ld l = C.r*C.r / d2;
    Point m = C.c + (p - C.c) * l;
    ld h = sqrtl(max((ld)0.0L, C.r*C.r - l*l*d2));
    Point v = Point(-(p.y - C.c.y) / d, (p.x - C.c.x) / d);
    return { m + v * (h/d), m - v * (h/d) };
}

// Polar sort by angle around origin or pivot
void polarSort(vector<Point> &pts, const Point &pivot = Point(0,0)){
    sort(pts.begin(), pts.end(), [&](const Point &a, const Point &b){
        Point A = a - pivot, B = b - pivot;
        int ha = (A.y > 0 || (A.y == 0 && A.x >= 0));
        int hb = (B.y > 0 || (B.y == 0 && B.x >= 0));
        if (ha != hb) return ha > hb;
        ld cr = cross(A, B); if (sgn(cr) != 0) return cr > 0; return norm2(A) < norm2(B);
    });
}

// Input helpers for fast reading
inline Point readPoint(){ double x,y; if(!(cin>>x>>y)) return Point(0,0); return Point(x,y); }

// Example main: shows usage and available functions
int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n; if (!(cin >> n)) return 0;
    vector<Point> pts(n); for (int i = 0; i < n; ++i) pts[i] = readPoint();
    auto hull = convexHull(pts);
    cout << (int)hull.size() << ' ';
    for (auto &p : hull) cout << p << ' ';
    auto C = mec(pts);
    cout << "MEC center: " << C.c << " radius: " << fixed << setprecision(12) << (double)C.r << ' ';
    return 0;
}
