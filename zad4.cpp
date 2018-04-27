
#define BOOST_UBLAS_TYPE_CHECK 0
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

typedef matrix<long double, row_major, std::vector<long double>> Matrix;
typedef vector<long double, std::vector<long double>> Vector;
struct IntegralEquation{
    std::vector<long double> a,rootLeg;
    Vector U;
    int N;
    IntegralEquation(int N) : N(N){}
    static long double f(const long double& x);    
    static long double K(const long double& x,const long double& t);
    Matrix Invert(const Matrix& _a); 
    template<typename type1, typename type2>
    void getU(type1 sysf=IntegralEquation::f,type2 sysK=IntegralEquation::K,long double alpha=1.);
    template<typename type1, typename type2>
    long double finally(const long double& x, type1 sysf=IntegralEquation::f,
                        type2 sysK=IntegralEquation::K,long double alpha=1.);   
};

struct first_kind{
    std::vector<long double>& a,&rootLeg;
    Vector &U;
    int &N;
    first_kind(IntegralEquation& IEq) : a(IEq.a), rootLeg(IEq.rootLeg),U(IEq.U),N(IEq.N) {}
    long double operator()(const long double& x){
        long double result=0;
        for(unsigned i=0;i<N;++i){
            result+=a[i]*IntegralEquation::K(x,rootLeg[i])*IntegralEquation::f(rootLeg[i]);
        }
        return result;
    }    
    long double operator()(const long double& x,const long double& t){
        long double result=0;
        for(unsigned i=0;i<N;++i){
            result-=a[i]*IntegralEquation::K(x,rootLeg[i])*IntegralEquation::K(t,rootLeg[i]);
        }
        return result;
    }
};

long double IntegralEquation::f(const long double& x){
    return 1/x-log(1+x)/(x*x);
}
long double IntegralEquation::K(const long double& x,const long double& t){
    return 1/(1+x*t);
}


template<typename type1, typename type2>
long double IntegralEquation::finally(const long double& x,type1 sysf, type2 sysK, 
                                        long double alpha){
    long double result=0;
    for(unsigned i=0;i<N;++i){
        result+=a[i]*sysK(x,rootLeg[i])*U(i);
    }
    result+=sysf(x);
    return result/alpha;
}
//Украдено с https://gist.github.com/lilac/2464434
Matrix IntegralEquation::Invert(const Matrix& _a)
{
    Matrix A(_a);
    permutation_matrix<decltype(A)::size_type> P(N);
    int res = lu_factorize(A, P);

    if( res != 0 )
        throw std::logic_error("lu_factorize error");
    
    Matrix inverse(identity_matrix<long double>(A.size1())); // identity_matrix - единичная матрица
    lu_substitute(A, P, inverse);
    //std::cout << prod(inverse,_a) << std::endl;
    
    return inverse;
}

template<typename type1, typename type2>
void IntegralEquation::getU(type1 sysf, type2 sysK,long double alpha){
    bool first=0; 
    if(rootLeg.size()==0){
        first=1;
        rootLeg=legendre_p_zeros<long double>(N);
        int j=rootLeg.size();
        for(unsigned i=0;i<j;++i){
            rootLeg.push_back(-rootLeg[i]);
        }
        std::sort(rootLeg.begin(),rootLeg.end());
    }
    
    Vector frootLeg(N);
    Matrix A(N,N);
    std::vector<long double> LegDer(N);
    for(unsigned i=0;i<N;++i)
        LegDer[i]= legendre_p_prime<long double>(N,rootLeg[i]);
    if(a.size()==0){
        for(unsigned i=0;i<N;++i){
            a.push_back(1/(1-rootLeg[i]*rootLeg[i])/(LegDer[i]*LegDer[i]));
        }
    }
    if(first){
        rootLeg=0.5+rootLeg/2.;
    }
    
    for(unsigned i=0;i<N;++i){
        frootLeg(i)=sysf(rootLeg[i]);
        for(unsigned j=0;j<N;++j){
            A(j,i)=-a[i]*sysK(rootLeg[i],rootLeg[j]);
        }
        A(i,i)+=alpha;
    }
    frootLeg(4)+=10e-9;
    /*for(unsigned i=0;i<N;++i){
        std::cout << "A(" << i << "): ";
        for(unsigned j=0;j<N;++j){
            std::cout <<  A(i,j) << " ";
        }
        std::cout << std::endl;
    }*/
    
    U=prod(Invert(A),frootLeg);
}
int main(){
    std::ios_base::sync_with_stdio(false);

    int N;
    long double alpha=0.1;
    std::ifstream file_in;
    std::ofstream file_alpha;
    file_in.open("file_in.dat",std::ios_base::in);
        file_in >> N;
    file_in.close();
    IntegralEquation IEq(N);
    first_kind IEq_first(IEq);
    std::vector<std::vector<long double> > data(16); 
    
    IEq.getU(IEq.f,IEq.K);
    for(unsigned i=0;i<IEq.U.size();++i){
        std::cout << "U[" << i << "]= " << (IEq.U)[i] << std::endl;
    }
    for(unsigned i=0;i<IEq.a.size();++i){
        std::cout << "a[" << i << "]= " << (IEq.a)[i] << std::endl;
    }
    for(long double i=0.;i<=10;++i)
        std::cout << "U(" << i/10 << ")= " << IEq.finally(i/10,IEq.f,IEq.K) << std::endl;

    for(long double i=0.;i<=50;++i)
        data[0].push_back(i/50);
    for(unsigned k=1;k<16;k++){
        alpha/=10;
        IEq.getU(IEq_first,IEq_first,alpha);
        for(long double i=0.;i<=50;++i)
            data[k].push_back(IEq.finally(i/50,IEq_first,IEq_first,alpha));
    }
    file_alpha.open("file_alpha.dat",std::ios_base::trunc); 
        for(unsigned i=0;i<data[0].size();++i){
            for(unsigned j=0;j<data.size();++j)
                file_alpha << data[j][i] << " ";
            file_alpha << std::endl;
        }
    file_alpha.close();
}