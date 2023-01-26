#include <R.h>

char *fcn_message(char *msg, int info, int n, int nit)
{
    if      (info == 1)
      snprintf(msg, 256, "Relative error in the sum of squares is at most `ftol'.");
    else if (info == 2)
      snprintf(msg, 256, "Relative error between `par' and the solution is at most `ptol'.");
    else if (info == 3)
      snprintf(msg, 256, "Conditions for `info = 1' and `info = 2' both hold.");
    else if (info == 4)
      snprintf(msg, 256, "The cosine of the angle between `fvec' and any column of the Jacobian is at most `gtol' in absolute value.");
    else if (info == 5)
      snprintf(msg, 256, "Number of calls to `fcn' has reached or exceeded `maxfev' == %d.", n);
    else if (info == 6)
      snprintf(msg, 256, "`ftol' is too small. No further reduction in the sum of squares is possible.");
    else if (info == 7)
      snprintf(msg, 256, "`ptol' is too small. No further improvement in the approximate solution `par' is possible.");
    else if (info == 8)
      snprintf(msg, 256, "`gtol' is too small. `fvec' is orthogonal to the columns of the Jacobian to machine precision.");
    else if (info < 0)
      snprintf(msg, 256, "Number of iterations has reached `maxiter' == %d.", nit);
    else if (info == 0)
      snprintf(msg, 256, "Improper input parameters.");
    return msg;
}
