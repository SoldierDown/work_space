//!#####################################################################
//! \file Krylov_System_Base.h
//!#####################################################################
// Class Krylov_System_Base
//######################################################################
#ifndef __Krylov_System_Base__
#define __Krylov_System_Base__

#include <Krylov_Vector_Base.h>

namespace Nova{
template<class T>
class Krylov_System_Base
{
  public:
    bool use_preconditioner,preconditioner_commutes_with_projection;

    //constructor
    Krylov_System_Base()
        :use_preconditioner(false),preconditioner_commutes_with_projection(false)
    {}
    Krylov_System_Base(const bool use_preconditioner_input,const bool preconditioner_commutes_with_projection_input)
        :use_preconditioner(use_preconditioner_input),preconditioner_commutes_with_projection(preconditioner_commutes_with_projection_input)
    {}

    virtual ~Krylov_System_Base() {}

//###################################################################### 
    virtual void Multiply(const Krylov_Vector_Base<T>& x,Krylov_Vector_Base<T>& result) const=0;
    virtual double Inner_Product(const Krylov_Vector_Base<T>& x,const Krylov_Vector_Base<T>& y) const=0;
    virtual T Convergence_Norm(const Krylov_Vector_Base<T>& x) const=0;
    // modifies system Ax=b to PAPx=Pb
    virtual void Project(Krylov_Vector_Base<T>& x) const=0;
    // implements Sx=Px+x0, where x0 is the desired component of x in the nullspace of P
    virtual void Set_Boundary_Conditions(Krylov_Vector_Base<T>& x) const=0;
    //virtual T ComputeError(Krylov_Vector_Base<T>& x) =0;
//###################################################################### 
    // removes component of x in nullspace of A (used to project residual for stopping conditions)
    virtual void Project_Nullspace(Krylov_Vector_Base<T>& x) const;
    //virtual T ComputeError(const Krylov_Vector_Base<T> &x);
    //virtual T Compute_Error(const Krylov_Vector_Base<T>&x);
    const Krylov_Vector_Base<T>& Precondition(const Krylov_Vector_Base<T>& r,Krylov_Vector_Base<T>& z) const;
    //virtual void ComputeError(const Krylov_Vector_Base<T>& x);
  protected:
    virtual void Apply_Preconditioner(const Krylov_Vector_Base<T>& r,Krylov_Vector_Base<T>& z) const;
//###################################################################### 
};
}
#endif