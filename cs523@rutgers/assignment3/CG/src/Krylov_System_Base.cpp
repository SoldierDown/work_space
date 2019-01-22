//!#####################################################################
//! \file Krylov_System_Base.cpp
//!#####################################################################
#include <Krylov_System_Base.h>
#include <Log/Debug_Utilities.h>

using namespace Nova;
//######################################################################
// Project_Nullspace
//######################################################################
template<class T> void Krylov_System_Base<T>::
Project_Nullspace(Krylov_Vector_Base<T>& x) const
{printf("function Project_Nullspace not define\n");}
//######################################################################
// Apply_Preconditioner
//######################################################################
template<class T> void Krylov_System_Base<T>::
Apply_Preconditioner(const Krylov_Vector_Base<T>& r,Krylov_Vector_Base<T>& z) const
{printf("function Apply_Preconditioner not define\n");}
//######################################################################
// Precondition
//######################################################################
template<class T> const Krylov_Vector_Base<T>& Krylov_System_Base<T>::
Precondition(const Krylov_Vector_Base<T>& r,Krylov_Vector_Base<T>& z) const
{
    if(!use_preconditioner) return r;
    Apply_Preconditioner(r,z);
    if(!preconditioner_commutes_with_projection){Project(z);Project_Nullspace(z);}
    return z;
}

template class Nova::Krylov_System_Base<float>;
#ifdef COMPILE_WITH_DOUBLE_SUPPORT
template class Nova::Krylov_System_Base<double>;
#endif