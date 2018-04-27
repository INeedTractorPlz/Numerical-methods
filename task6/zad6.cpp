#include <boost/math/special_functions/legendre.hpp>
#include<iostream>
#include <boost/numeric/ublas/vector.hpp>
#include<boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include<boost/shared_ptr.hpp>
#include<algorithm>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include<fstream>
#include<ios>



using namespace boost::math;
using namespace boost::numeric::ublas;
using namespace boost;

typedef std::vector<std::vector<long double> > state_type;

template<typename T, typename U>
std::vector<T> operator/(const std::vector<T>& a, U b){
    std::vector<T> c(a.size());
    for(unsigned i=0;i<a.size();++i)
        c[i]=std::move(a[i]/b);
    return c;
}
template<typename T, typename U>
std::vector<T> operator*(const std::vector<T>& a, U b){
    std::vector<T> c(a.size());
    for(unsigned i=0;i<a.size();++i)
        c[i]=std::move(a[i]*b);
    return c;
}
template<typename T, typename U>
std::vector<T> operator*(U b, const std::vector<T>& a){
    std::vector<T> c(a.size());
    for(unsigned i=0;i<a.size();++i)
        c[i]=std::move(a[i]*b);
    return c;
}

template<typename T>
std::vector<T> operator+(const std::vector<T>& a,const std::vector<T>& b){
    std::vector<T> c(a.size());
    for(unsigned i=0;i<a.size();++i)
        c[i]=std::move(a[i]+b[i]);
    return c;
}

template<typename T, typename U>
std::vector<T> operator+(U b, const std::vector<T>& a){
    std::vector<T> c(a.size());
    for(unsigned i=0;i<a.size();++i)
        c[i]=std::move(a[i]+b);
    return c;
}

template<typename T, typename U>
std::vector<T> operator+(const std::vector<T>& a,U b){
    std::vector<T> c(a.size());
    for(unsigned i=0;i<a.size();++i)
        c[i]=std::move(a[i]+b);
    return c;
}
struct Sweep{
    std::vector<long double> alpha, betta;
    std::vector<long double>& a, &b, &c, &d;
    int n;
    Sweep(std::vector<long double>& a, std::vector<long double>& b,
    std::vector<long double>& c, std::vector<long double>& d) : a(a), b(b), c(c), d(d){
        n=a.size()-1;
        alpha.resize(n+1); betta.resize(n+1);
        alpha[0]=c[0]/b[0];
	    betta[0]=-d[0]/b[0];
	    for(unsigned i=1;i<n+1;i++)
	    {
	        alpha[i]=c[i]/(b[i]-a[i]*(alpha[i-1]));
	        betta[i]=(a[i]*(betta[i-1])-d[i])/(b[i]-a[i]*(alpha[i-1])) ;
	    }
    }
    template<typename type>
    void GetY(type y)
    {
        y[n]=betta[n]; //i=i-1; <=> i--;
        for(int i=n-1;i>=0;i--)
	        y[i]=(y[i+1])*alpha[i]+betta[i];
    }
};

struct HeatEquation{
    state_type &u;
    int n,m,M,N;
    std::vector<long double> x,t;
    long double h,tau,sigma,tau_stable;
    HeatEquation(state_type &result, int n, int m, long double a, long double b, long double T) : 
    u(result), n(n) , m(m){
        h=(b-a)/n;
        tau=T/m;
        N=ceil(100*T*tau/(m*h*h));
        tau_stable=tau/N;
        std::cout << "tau= " << tau << std::endl;
        std::cout << "tau_stable= " << tau_stable << std::endl;
        M=N*m;
        sigma=tau_stable/(h*h);
        x.resize(n+1);
        x[0]=a;
        for(unsigned i=1;i<x.size();++i){
            x[i]=x[i-1]+h;
        }
        t.resize(M+1);
        t[0]=0;
        for(unsigned i=1;i<t.size();++i){
            t[i]=t[i-1]+tau_stable;
        }
    }

    static long double  psi1(long double t){
        return exp(-0.25*t);
    }
    static long double psi2(long double t){
        return exp(-0.25*t)*cos(0.5);
    }
    static long double phi(long double x){
        return cos(0.5*x)+(1-x)*x;
    }
    static long double f(long double x, long double t){
        return exp(-t)*(x*x-x+2);
    }
    void GetImplicitResult(){
        u.resize(0);
        for(unsigned k=0;k<m+1;++k){
            u.push_back(std::vector<long double>(n+1));
        }
        for(unsigned i=0;i<n+1;++i)
            u[0][i]=phi(x[i]);
        for(unsigned k=1;k<m+1;++k){
            u[k][0]=psi1(t[k*N]); u[k][n]=psi2(t[k*N]);
        }
        
        sigma=tau/(h*h);
        std::vector<long double> a,b,c,d;
            for(unsigned i=0;i<n-1;++i){
                a.push_back(sigma);
                b.push_back(2*sigma+1);
                c.push_back(sigma);
            }
            a[0]=0;c[n-2]=0;
            d.resize(n-1);
        for(unsigned k=1;k<=m;++k){
            for(unsigned i=1;i<n;++i){
                d[i-1]=-u[k-1][i]-tau*f(x[i],t[k*N]);
            }
            d[0]+=-sigma*psi1(t[k*N]);d[n-2]+=-sigma*psi2(t[k*N]);
            Sweep run(a,b,c,d);
            run.GetY(u[k].begin()+1);
        }
    }
    void GetExplicitResult(){
        u.resize(0);
        for(unsigned k=0;k<M+1;++k){
            u.push_back(std::vector<long double>(n+1));
        }
        for(unsigned i=0;i<n+1;++i)
            u[0][i]=phi(x[i]);
        for(unsigned k=1;k<M+1;++k){
            u[k][0]=psi1(t[k]); u[k][n]=psi2(t[k]);
        }
        for(unsigned k=0;k<M;++k){
            for(unsigned i=1;i<n;++i){
                u[k+1][i]=sigma*u[k][i+1]-(2*sigma-1)*u[k][i]+sigma*u[k][i-1]+tau_stable*f(x[i],t[k]);
            }
        }
    }
};

long double Solution(long double x, long double t){
    return exp(-0.25*t)*cos(0.5*x)+exp(-t)*(1-x)*x;
}

int main(){
    std::ios_base::sync_with_stdio(false);

    std::ifstream file_in;
    std::ofstream explicit_result, implicit_result, solution_result;
    int n, m; 
    state_type result;
    long double a, b, T;
    file_in.open("file_in.dat",std::ios_base::in);
        file_in >> n >> m >> a >> b >> T;
    file_in.close();
    
    HeatEquation he(result, n, m, a, b, T);


    he.GetExplicitResult();
    explicit_result.open("explicit.dat",std::ios_base::trunc); 
        for(unsigned i=0;i<result.size();i=i+he.N){
            explicit_result << he.t[i] << " ";
            for(unsigned j=0;j<result[0].size();++j)
                explicit_result << result[i][j] << " ";
            explicit_result << std::endl;
        }
    explicit_result.close();

    solution_result.open("solution.dat",std::ios_base::trunc); 
        for(unsigned i=0;i<result.size();i=i+he.N){
            solution_result << he.t[i] << " ";
            for(unsigned j=0;j<result[0].size();++j)
                solution_result << fabs(result[i][j] - Solution(he.x[j],he.t[i])) << " ";
            solution_result << std::endl;
        }
    solution_result.close();

    he.GetImplicitResult();
    implicit_result.open("implicit.dat",std::ios_base::trunc); 
        for(unsigned i=0;i<result.size();++i){
            implicit_result << he.t[i*(he.N)] << " ";
            for(unsigned j=0;j<result[0].size();++j)
                implicit_result << result[i][j] << " ";
            implicit_result << std::endl;
        }
    implicit_result.close();
    return 0;
}
