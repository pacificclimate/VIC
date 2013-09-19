
#ifndef NEWT_RAPH_FUNC_FAST_H_
#define NEWT_RAPH_FUNC_FAST_H_

class NewtonRaphsonMethod {
public:
  NewtonRaphsonMethod(double T_2[], double res[], int n, double deltat, int FS_ACTIVE,
      int NOFLUX, int EXP_TRANS, double *T0, double *moist, double *ice,
      double *kappa, double *Cs, double *max_moist, double *bubble, double *expt,
      double *porosity, double *effective_porosity, double *alpha, double *beta,
      double *gamma, double *Zsum, double Dp, double *bulk_dens_min,
      double *soil_dens_min, double *quartz, double *bulk_density,
      double *soil_density, double *organic, double *depth, int Nlayers);
  virtual ~NewtonRaphsonMethod() {}
  int compute(double x[], int n);
private:
  virtual void fda_heat_eqn(double [], double [], int n, int focus);
  void fdjac3(double x[], double fvec[], double a[], double b[], double c[], int n);

  double deltat;
  int FS_ACTIVE;
  int NOFLUX;
  int EXP_TRANS;
  double *T0;
  double *moist;
  double *ice;
  double *kappa;
  double *Cs;
  double *max_moist;
  double *bubble;
  double *expt;
  double *porosity;
  double *effective_porosity;
  double *alpha;
  double *beta;
  double *gamma;
  double *Zsum;
  double Dp;
  double *bulk_dens_min;
  double *soil_dens_min;
  double *quartz;
  double *bulk_density;
  double *soil_density;
  double *organic;
  double *depth;
  int Nlayers;

  double Bexp;
  // variables used to calculate residual of the heat equation
  // defined here
  double Ts;
  double Tb;
};



#endif /* NEWT_RAPH_FUNC_FAST_H_ */
