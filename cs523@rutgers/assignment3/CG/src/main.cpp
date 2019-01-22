#include <iostream>
#include <CG_System.h>
#include <Arrays/Array.h>
#include <Krylov_System_Base.h>
#include <CG_Vector.h>
#include <Conjugate_Gradient.h>
double Phi(double x_, double y_)
{
    return x_ * x_ + y_ * y_ - R * R;
}

int Index(int i_, int j_)
{
    return j_ * (RES + 1) + i_;
}
using namespace Nova;
int main()
{
    double f;
    double x_loc, y_loc;
    double cellwith = 1.0 / RES;
    // const int tmp;
    // std::cout << "Enter (RES+1)olution: " << std::endl;
    // std::cin >> tmp;

    //Vector_Base& x,const Vector_Base& b,Vector_Base& q,Vector_Base& s,Vector_Base& r,Vector_Base& k,Vector_Base& z,
    //const T tolerance,const int min_iterations,const int max_iterations
    Array<Vector<double, 2>> xarray = Array<Vector<double, 2>>((RES + 1) * (RES + 1), 0);
    CG_Vector<double, 2> x = CG_Vector<double, 2>(xarray);
    Array<Vector<double, 2>> barray = Array<Vector<double, 2>>((RES + 1) * (RES + 1), 0);
    for (size_t j = 0; j < RES + 1; j++)
    {
        y_loc = -0.5 + j * cellwith;
        for (size_t i = 0; i < RES + 1; i++)
        {
            x_loc = -0.5 + i * cellwith;
            f = Phi(x_loc, y_loc);
            if (f < 0)
            {
                //std::cout << "inside: " << Index(i, j) << std::endl;
                barray(Index(i, j))(0) = 4;
            }
        }
    }
    CG_Vector<double, 2> b = CG_Vector<double, 2>(barray);
    Array<Vector<double, 2>> qarray = Array<Vector<double, 2>>((RES + 1) * (RES + 1), 0);
    CG_Vector<double, 2> q = CG_Vector<double, 2>(qarray);
    Array<Vector<double, 2>> sarray = Array<Vector<double, 2>>((RES + 1) * (RES + 1), 0);
    CG_Vector<double, 2> s = CG_Vector<double, 2>(sarray);
    Array<Vector<double, 2>> rarray = Array<Vector<double, 2>>((RES + 1) * (RES + 1), 0);
    CG_Vector<double, 2> r = CG_Vector<double, 2>(rarray);
    Array<Vector<double, 2>> karray = Array<Vector<double, 2>>((RES + 1) * (RES + 1), 0);
    CG_Vector<double, 2> k = CG_Vector<double, 2>(karray);
    Array<Vector<double, 2>> zarray = Array<Vector<double, 2>>((RES + 1) * (RES + 1), 0);
    CG_Vector<double, 2> z = CG_Vector<double, 2>(zarray);

    Conjugate_Gradient<double> cg;
    CG_System<double, 2> cgs;
    //TEST
    //Array<Vector<double, 2>> IParray1 = Array<Vector<double, 2>>((RES + 1) * (RES + 1), -1);
    //CG_Vector<double, 2> IP1 = CG_Vector<double, 2>(IParray1);
    //std::cout<<cgs.ComputeError(IP1)<<std::endl;
    // Array<Vector<double, 2>> IParray2 = Array<Vector<double, 2>>((RES + 1) * (RES + 1), 2);
    // CG_Vector<double, 2> IP2 = CG_Vector<double, 2>(IParray2);
    //std::cout << cgs.Convergence_Norm(IP1) << std::endl;
    //std::cout << x.array << std::endl;
    cg.Solve(cgs, x, b, q, s, r, k, z, (double)1e-6, 1, 200);
    //std::cout << "x:" << CG_Vector<double, 2>::CG_Array(x) << std::endl;
    std::cout << "x error:" << cgs.ComputeError(x) << std::endl;
    //cgs.Multiply(x, b);
    //std::cout << "b: " << CG_Vector<double, 2>::CG_Array(b) << std::endl;
    return 0;
}
