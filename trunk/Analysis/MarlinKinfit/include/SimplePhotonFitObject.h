 
#ifndef __SIMPLEPHOTONFITOBJECT_H
#define __SIMPLEPHOTONFITOBJECT_H

#include "ParticleFitObject.h"

// Class SimplePhotonFitObject
/// Class for ISR/Beamstrahlung photons with (p_x,p_y (both fix),  p_z (free)) in kinematic fits
/**
 *
 * Author: Moritz Beckmann
 * $Date: 2010/06/11 20:32:51 $
 * $Author: mbeckman $
 *
 * \b Changelog:
 * - 
 *
 * \b CVS Log messages:
 * - $Log: SimplePhotonFitObject.h,v $
 * - Revision 1.1  2010/06/11 20:32:51  mbeckman
 * - Renamed PhotonFitObjects, cleaned them up for usage
 * -
 * - Revision 1.3  2009/02/23 12:03:18  mbeckman
 * - - PhotonFitObject:     bug fix (1/0), removed dispensable variables
 * - - PhotonFitObjectPxyg: bug fixes (1/0, order of computing variables), modified parametrization
 * -
 * - Revision 1.2  2009/02/18 11:53:42  mbeckman
 * - documentation, debug output
 *
 */ 
class SimplePhotonFitObject : public ParticleFitObject {
  public:
    SimplePhotonFitObject(double px, double py, double pz,   /// initial values for photon (p_x,p_y fix)
                 double Dpz);                          /// sigma for pz
    virtual ~SimplePhotonFitObject();
    
    /// Get name of parameter ilocal
    virtual const char *getParamName (int ilocal     ///< Local parameter number
                                     ) const;
 
    /// Set value of parameter ilocal; return: significant change
    virtual bool setParam (int ilocal,    ///< Local parameter number
                           double par_    ///< New parameter value
                          );  
    /// Set value and measured flag of parameter ilocal; return=significant change
    virtual bool   setParam (int ilocal,         ///< Local parameter number
                             double par_,        ///< New parameter value
                             bool measured_,     ///< New "measured" flag
                             bool fixed_ = false ///< New "fixed" flag
                            );  
    /// Read values from global vector, readjust vector; return: significant change
    virtual bool   updateParams (double p[],   ///< The parameter vector
                                 int idim      ///< Length of the vector                         
                                );  
    
    
    // these depend on actual parametrisation!
    virtual double getPx() const;
    virtual double getPy() const;
    virtual double getPz() const;
    virtual double getE() const;
    
    virtual double getP() const;
    virtual double getP2() const;  
    virtual double getPt() const;
    virtual double getPt2() const;
     
    virtual double getDPx(int ilocal) const;
    virtual double getDPy(int ilocal) const;
    virtual double getDPz(int ilocal) const;
    virtual double getDE(int ilocal) const;
    
    /// add derivatives to vector der of size idim
    /// pxfact*dpx/dx_i + pyfact*dpy/dx_i + pzfact*dpz/dx_i + efact*dE/dx_i 
    virtual void   addToDerivatives (double der[],      ///< Derivatives vector, length idim
                                     int idim,          ///< Length of derivatives vector
                                     double efact=0,    ///< Factor for dE/dx_i
                                     double pxfact=0,   ///< Factor for dpx/dx_i
                                     double pyfact=0,   ///< Factor for dpy/dx_i 
                                     double pzfact=0    ///< Factor for dpz/dx_i
                                     ) const;

    /// add second order derivatives to matrix der2 of size idim x idim
    /// pxfact*d^2px/(dx_i dx_j) + pyfact...
    virtual void   addTo2ndDerivatives (double der2[],  ///< Derivatives vector, size idim x idim
                                        int    idim,    ///< First dimension of derivatives matrix
                                        double efact,   ///< Factor for d^2E/dx_i dx_j
                                        double pxfact,  ///< Factor for d^2px/dx_i dx_j
                                        double pyfact,  ///< Factor for d^2py/dx_i dx_j
                                        double pzfact   ///< Factor for d^2pz/dx_i dx_j
                                       ) const;
    /// Add second order derivatives to matrix M of size idim x idim
    /// der[0]*d^2E/(dx_i dx_j) + der[1]*d^2px/(dx_i dx_j) + ...
    virtual void   addTo2ndDerivatives (double M[],     ///< Global derivatives matrix size idim x idim
                                        int    idim,    ///< First dimension of derivatives matrix
                                        double lamda,   ///< Global factor
                                        double der[]    ///< Factors for d^2(E,px,py,pz)/dx_i dx_j
                                       ) const;
    /// Add first order derivatives to matrix M of size idim x idim
    /// der[0]*dE/dx_i + der[1]*dpx/dx_i + ...
    virtual void   addTo1stDerivatives (double M[],     ///< Global derivatives matrix size idim x idim
                                        int    idim,    ///< First dimension of derivatives matrix
                                        double der[],   ///< Factors for d^2(E,px,py,pz)/dx_i dx_j
                                        int kglobal     ///< Global parameter number of constraint
                                       ) const;
    
    virtual void addToGlobalChi2DerVector (double *y,     ///< Vector of chi2 derivatives
                                           int idim,      ///< Vector size 
                                           double lambda, ///< The lambda value
                                           double der[]   ///< 4-vector with dg/dE, dg/dpx, dg/dpy, dg/dpz
                                           ) const;
    /// Calculates the squared error for a quantity with derivatives w.r.t. E, dx, dy, dz
    virtual double getError2 (double der[]    ///< Factors for d(E,px,py,pz)/dx_i
                             ) const;

    virtual void invalidateCache() const;
  
  protected:
    
    virtual void initCov();
    
    void updateCache() const;
  
    mutable bool cachevalid;
    
    mutable double pt2, p2, p,
                   dE0, dE1, dE2, 
                   chi2;

//     mutable double dpx0, dpy0, dpz0, dpx1, dpy1, dpz1, dpx2, dpy2, dpz2;
  
};



#endif // __SIMPLEPHOTONFITOBJECT_H
