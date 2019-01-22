//!#####################################################################
//! \file CG_System.h
//!#####################################################################
// Class CG_System
//######################################################################
#ifndef __CG_System__
#define __CG_System__

#include <Krylov_System_Base.h>
#include <CG_Vector.h>
#define R 0.25
#define RES 16
namespace Nova
{
template <class T, int d>
class CG_System : public Krylov_System_Base<T>
{
    using TV = Vector<T, d>;
    using Base = Krylov_System_Base<T>;
    using Vector_Base = Krylov_Vector_Base<T>;

  public:
    std::vector<int> flag;
    T cellwidth;
    T error;
    int res;
    CG_System(int res_)
    {
        res = res_;
        cellwidth = 1.0 / res;
        Init();
    }
    CG_System()
    {
        res = RES;
        cellwidth = 1.0 / res;
        Init();
    }

    void Init()
    {
        T x, y;
        T f;
        for (size_t i = 0; i < res + 1; i++)
        {
            x = -0.5 + cellwidth * i;
            for (size_t j = 0; j < res + 1; j++)
            {
                y = -0.5 + cellwidth * j;
                f = Phi(x, y);

                if (f < 0)
                    flag.push_back(-1);
                else
                    flag.push_back(0);
            }
        }
    }
    void Multiply(const Vector_Base &x, Vector_Base &b) const
    {
        const Array<TV> &x_array = CG_Vector<T, d>::CG_Array(x);
        Array<TV> &b_array = CG_Vector<T, d>::CG_Array(b);
        int index = 0;
        int top, bottom, left, right;
        T laplace;
        for (size_t j = 0; j < res + 1; j++)
        {
            for (size_t i = 0; i < res + 1; i++)
            {
                index = Index(i, j);
                if (flag[index] == -1)
                {
                    laplace = 0.0;
                    left = Index(i - 1, j);
                    right = Index(i + 1, j);
                    top = Index(i, j + 1);
                    bottom = Index(i, j - 1);
                    laplace -= 4 * x_array(index)(0);
                    if (left >= 0 && flag[left] == -1)
                    {
                        laplace += x_array(left)(0);
                    }
                    if (bottom >= 0 && flag[bottom] == -1)
                    {
                        laplace += x_array(bottom)(0);
                    }
                    if (right < (res + 1) * (res + 1) && flag[right] == -1)
                    {
                        laplace += x_array(right)(0);
                    }
                    if (top < (res + 1) * (res + 1) && flag[top] == -1)
                    {
                        laplace += x_array(top)(0);
                    }
                    //std::cout << index << ": " << laplace << std::endl;
                    b_array(index)(0) = laplace / (cellwidth * cellwidth);
                }
            }
        }
        //std::cout << b_array << std::endl;
        b = CG_Vector<T, d>(b_array);
        //std::cout << CG_Vector<T, d>::CG_Array(x) << std::endl;
    }

    T Inner_Product(const Vector_Base &x, const Vector_Base &y) const
    {
        const Array<TV> &x_array = CG_Vector<T, d>::CG_Array(x);
        const Array<TV> &y_array = CG_Vector<T, d>::CG_Array(y);

        T result = (T)0.;
        for (size_t i = 0; i < (res + 1) * (res + 1); ++i)
        {
            if (flag[i] == -1)
                result += x_array(i)(0) * y_array(i)(0);
        }
        return result;
    }

    T Convergence_Norm(const Vector_Base &x) const
    {
        //const Array<TV> &x_array = CG_Vector<T, d>::CG_Array(x);
        const Array<TV> &x_array = CG_Vector<T, d>::CG_Array(x);
        T max = (T)-100.0;
        for (size_t i = 0; i < (res + 1) * (res + 1); ++i)
            if (flag[i] == -1)
            {
                //std::cout << i << ": " << x_array(i)(0) << std::endl;
                max = fabs(x_array(i)(0)) > max ? fabs(x_array(i)(0)) : max;
            }

        return max;
    }

    T ComputeError(const Vector_Base &x)
    {
        const Array<TV> &x_array = CG_Vector<T, d>::CG_Array(x);
        int index;
        T err = (T)-100.0;
        T result;
        T f;
        T x_loc, y_loc;
        for (size_t j = 0; j < res + 1; j++)
        {
            y_loc = -0.5 + j * cellwidth;
            for (size_t i = 0; i < res + 1; i++)
            {
                index = Index(i, j);

                if (flag[index] == -1)
                {
                    x_loc = -0.5 + i * cellwidth;
                    f = Phi(x_loc, y_loc);
                    //std::cout << "phi(" << index << "): " << f << std::endl;
                    //std::cout << "x(" << index << "): " << x_array(index)(0) << std::endl;
                    result = x_array(index)(0) - f;
                    //std::cout << "result(" << index << "): " << result << std::endl;
                    err = fabs(result) > err ? fabs(result) : err;
                }
            }
        }

        return err;
    }
    /*
    T ComputeError(Vector_Base &x)
    {

    }
*/
    // void ComputeError(Vector_Base &x) const
    // {
    //     Array<TV> &x_array = CG_Vector<T, d>::CG_Array(x);
    //     T x_loc, y_loc;
    //     for (size_t i = 0; i < RES + 1; i++)
    //     {
    //         x_loc = -0.5 + i * cellwidth;
    //         for (size_t j = 0; j < RES + 1; j++)
    //         {
    //             y_loc = -0.5 + j * cellwidth;
    //             x_array(Index(i, j))(0) -= Phi(x_loc, y_loc);
    //         }
    //     }
    //     error = Convergence_Norm(CG_Vector<T, d>(x_array));
    // }
    void Project(Vector_Base &x) const
    {
        return;
    }

    void Set_Boundary_Conditions(Vector_Base &v) const { return; }
    void Project_Nullspace(Vector_Base &x) const { return; }
    int Index(int i_, int j_) const
    {
        return j_ * (res + 1) + i_;
    }
    T Phi(T x_, T y_)
    {
        return x_ * x_ + y_ * y_ - R * R;
    }
    //const Vector_Base& Precondition(const Vector_Base& r,Vector_Base& z) const{return ();}
    //void Apply_Preconditioner(const Vector_Base& r,Vector_Base& z) const{return;}
};
} // namespace Nova
#endif
