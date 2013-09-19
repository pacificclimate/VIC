
#ifndef ROOT_BRENT_H_
#define ROOT_BRENT_H_

class RootBrent {
public:
  RootBrent() {}
  virtual ~RootBrent() {}
  double root_brent(double LowerBound, double UpperBound, char* ErrorString);
  virtual double calculate(double) = 0;
};

#endif /* ROOT_BRENT_H_ */
