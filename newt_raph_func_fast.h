
#ifndef NEWT_RAPH_FUNC_FAST_H_
#define NEWT_RAPH_FUNC_FAST_H_

class NewtonRaphsonMethod {
public:
  NewtonRaphsonMethod(double T_2[], double res[], int n, double deltat, int FS_ACTIVE,
      int NOFLUX, int EXP_TRANS, double *T0, double *moist, double *ice,
      double *kappa, double *Cs, const double *max_moist, const double *bubble, const double *expt,
      const double *porosity, const double *effective_porosity, const double *alpha, const double *beta,
      const double *gamma, const double *Zsum, double Dp, const double *bulk_dens_min,
      const double *soil_dens_min, const double *quartz, const double *bulk_density,
      const double *soil_density, const double *organic, const double *depth, int Nlayers);
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
  const double *max_moist;
  const double *bubble;
  const double *expt;
  const double *porosity;
  const double *effective_porosity;
  const double *alpha;
  const double *beta;
  const double *gamma;
  const double *Zsum;
  double Dp;
  const double *bulk_dens_min;
  const double *soil_dens_min;
  const double *quartz;
  const double *bulk_density;
  const double *soil_density;
  const double *organic;
  const double *depth;
  int Nlayers;

  double Bexp;
  // variables used to calculate residual of the heat equation
  // defined here
  double Ts;
  double Tb;
};



#endif /* NEWT_RAPH_FUNC_FAST_H_ */
