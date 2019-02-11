#include <bits/stdc++.h>
using namespace std;

#define debug(...) do {\
    cerr << __PRETTY_FUNCTION__ << ":" << __LINE__ << " - ";\
    _debug(#__VA_ARGS__, __VA_ARGS__);\
} while(0)

template<typename T>
void _debug(const char* name, T&& head) {
    cerr << name << "=" << head << endl;
}

template<typename T, typename... K>
void _debug(const char* name, T&& head, K&&... tail) {
    int c = 0;
    while (*name != ',' or c != 0) {
        if (*name == '(' or *name == '[' or *name == '{') {
            c++;
        }
        if (*name == ')' or *name == ']' or *name == '}') {
            c--;
        }
        cerr << *name++;
    }
    cerr << "=" << head << ", ";
    _debug(name + 1, tail...);
}

template<typename I>
ostream& _out(ostream& os, I b, I e, char ob, char cb){
    os << ob;
    for (auto it = b; it != e; it++) {
        os << (it == b ? "" : " ") << *it;
    }
    return os << cb;
}

template<typename A, typename B>
ostream& operator<<(ostream &os, const pair<A, B> &p) {
    return os << "(" << p.first << ", " << p.second << ")";
}

template<typename T>
ostream& operator<<(ostream &os, const vector<T> &v) {
    return _out(os, v.begin(), v.end(), '[', ']');
}

template<typename T, size_t N>
ostream& operator<<(ostream &os, const array<T, N> &a) {
    return _out(os, a.begin(), a.end(), '[', ']');
}

template<typename T>
ostream& operator<<(ostream &os, const set<T> &s) {
    return _out(os, s.begin(), s.end(), '{', '}');
}

template<typename A, typename B>
ostream& operator<<(ostream &os, const map<A, B> &m) {
    return _out(os, m.begin(), m.end(), '{', '}');
}

template<typename T>
ostream& operator<<(ostream &os, const unordered_set<T> &s) {
    return _out(os, s.begin(), s.end(), '{', '}');
}

template<typename A, typename B>
ostream& operator<<(ostream &os, const unordered_map<A, B> &m) {
    return _out(os, m.begin(), m.end(), '{', '}');
}
