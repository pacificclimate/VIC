
#ifndef ROOT_BRENT_H_
#define ROOT_BRENT_H_

class RootBrent {
public:
  RootBrent() {}
  virtual ~RootBrent() {}
  double root_brent(double LowerBound, double UpperBound, char* ErrorString);
  virtual double calculate(double) = 0;
  //if the result of any calculation is less than -998 (ie -999) then there is an error
  static bool resultIsError(double result) { return result <= -998; }
};

#endif /* ROOT_BRENT_H_ */
