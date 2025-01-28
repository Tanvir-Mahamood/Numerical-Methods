#include <bits/stdc++.h>
using namespace std;
#define ll long long
#define nn  "\n"

                                   // Basic Globally declared variables & container
double p, m, hh = 0.0001;          // hh is a small value used for differentiation
ll n;                              // number of data
vector<vector<long double>> table; // Difference table for Interpolation
vector<long double> xx, yy;        // Year and Population holder vector
double Error = 0.0001;
double a0, a1, a2;                 // coefficients of polynomial of degree 2

// In order to use curve fitting algorithm (Least Square Curve fitting), we assume that the data may fit in a polynomial of order 2.
// In this case, we will be needed to solve a system of liner equations.
// In case of order 2 polynomial, there will be 3 linear equations that can be solved using Cramer's Rule.
// Here, I have assumed that the data can be fitted into a polynomial of order 2

ll determinant(ll mat[3][3]) // This function is used to calculate determinant of a (3 x 3) matrix
{
    ll ans;
    ans = + mat[0][0] * (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2]) // Formula of determinant of matrix
          - mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0])
          + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
    return ans;
}

void findSolution(ll coeff[3][4]) // This whole function is used for applying Cramer's Rule.
{
    ll del[3][3] =
    {
        { coeff[0][0], coeff[0][1], coeff[0][2] },
        { coeff[1][0], coeff[1][1], coeff[1][2] },
        { coeff[2][0], coeff[2][1], coeff[2][2] },
    };
    ll delx[3][3] =
    {
        { coeff[0][3], coeff[0][1], coeff[0][2] },
        { coeff[1][3], coeff[1][1], coeff[1][2] },
        { coeff[2][3], coeff[2][1], coeff[2][2] },
    };
    ll dely[3][3] =
    {
        { coeff[0][0], coeff[0][3], coeff[0][2] },
        { coeff[1][0], coeff[1][3], coeff[1][2] },
        { coeff[2][0], coeff[2][3], coeff[2][2] },
    };
    ll delz[3][3] =
    {
        { coeff[0][0], coeff[0][1], coeff[0][3] },
        { coeff[1][0], coeff[1][1], coeff[1][3] },
        { coeff[2][0], coeff[2][1], coeff[2][3] },
    };

    ll D  = determinant(del);
    ll Dx = determinant(delx);
    ll Dy = determinant(dely);
    ll Dz = determinant(delz);

    if (D != 0) // Coefficients have a unique solution. Apply Cramer's Rule
    {
        a0 = (double)Dx / D;
        a1 = (double)Dy / D;
        a2 = (double)Dz / D;
    }
    else
    {
        if (Dx == 0 && Dy == 0 && Dz == 0) cout << "Infinite solutions\n";
        else if (Dx != 0 || Dy != 0 || Dz != 0) cout << "No solution\n";
        return; // Program has stopped, cannot be solved.
    }
}

double get_p(ll k, int sig)             // This function is used for the construction of Newton's Interpolation formula
{
    k--;
    double ans = m*(p + k * sig);
    m = ans;
    return ans;
}

double Newton_For(double p, ll n) // Newton's Forward Interpolation formulation
{
    long long fact = 1;
    double y = table[0][0];
    for(ll i=1; i<n; i++) fact *= i, y += ((long double)get_p(i,-1) / fact) * table[0][i];
    return y;
}

double Newton_Bac(double p, ll n) // Newton's Backward Interpolation formulation
{
    long long fact = 1;
    double y = table[n-1][0];
    for(ll i=1; i<n; i++) fact *= i, y += ((double)get_p(i, 1) / fact ) * table[n-1-i][i];
    return y;
}

// This function decides which interpolation is to be applied.
// Accroding to the direction of target data, Newton's Forward or Backward Interpolation is applied.
long double Function(double x)
{
    ll t;
    double h = xx[1] - xx[0];

    if(abs(x - xx[0]) >= abs(x - xx[n-1])) t = 1; // Forward Interpolation should be called
    else t = 2;                                   // Backward Interpolation should be called

    if(t == 1)
    {
        p = (double)(x - xx[0]) / h;
        m = 1;
        return Newton_For(p, n); // Calling Newton's Forward Interpolation
    }
    else if(t == 2)
    {
        p = (double)(x - xx[n-1]) / h;
        m = 1;
        return Newton_Bac(p, n); // Calling Newton's Backward Interpolation
    }
}

// It is easy to run a Newton-Raphson method on an equation.
// But a tricky modification can made this method applicable for tabulated data of an unknown equation.
// So, we don't need to know the original equation to apply Newton-Raphson method.
// Tabular data is very enough to do the task.
// Modification:
//      Iterative Formula: x_(n+1) = x_n - f(x_n) / f'(x_n)
//      For example: x_1 = x_0 - f(x_0) / f'(x_1)
//      as we don't know the equation (don't know about the nature of f(x)),
//          we may use interpolation(x) instead of f(x). I have just renamed interpolation() as Function()
//          we may use differential formula instead of f'(x). Here, f'(x) = lim h->0 (f(x+h) - f(x-h)) / 2h

// Well, our equation is: population(x) - target = 0
// So, f(x) is replaced with Function(x) - target = 0, and the differential equation is simply f'(x).
// As 'target' is a constant, it does not hamper the differential equation. Thus, we can find our modified formula.
long double Newton_Raphson(long double target)
{
    double past, present, x = xx[0]; //initial guess
    past = x;

    x = x - ((double)(Function(x) - target)) / ((Function(x + hh) - Function(x - hh)) / (2 * hh));

    while(1)
    {
        past = x;
        x = x - ((double)(Function(x) - target)) / ((Function(x + hh) - Function(x - hh)) / (2 * hh));
        present = x;
        if(abs(present - past) < Error) break;
    }
    return x;
}

int main()
{
    cout << fixed << setprecision(5);
    ifstream inputFile("population_data.txt");

    ll i, j, t = 100, target_year, target_pop;
    ll sum_x, sum_x2, sum_x3, sum_x4, sum_y, sum_xy, sum_x2y;
    long double  val_x, val_y, h, x, y;

    cout << "Number of data: ";
    inputFile >> n;

    xx.resize(n);
    yy.resize(n);

    vector<vector<long double>> grid(n,vector<long double>(n));

    sum_x = sum_x2 = sum_x3 = sum_x4 = sum_y = sum_xy = sum_x2y = 0;
    for(i=0; i<n; i++) // Input Loop
    {
        inputFile >> val_x >> val_y;
        xx[i] = val_x;
        yy[i] = val_y;
        grid[i][0] = val_y;

//      Here, calculating some parameters for curve fitting algorithm.

        sum_x  += xx[i];                           // Summation of x
        sum_x2 += (xx[i] * xx[i]);                 // Summation of x^2
        sum_x3 += (xx[i] * xx[i] * xx[i]);         // Summation of x^3
        sum_x4 += (xx[i] * xx[i] * xx[i] * xx[i]); // Summation of x^4

        sum_y   +=  yy[i];                         // Summation of y
        sum_xy  += (xx[i] * yy[i]);                // Summation of x*y
        sum_x2y += (xx[i] * xx[i] * yy[i]);        // Summation of x^2*y
    }

    // System of linear equations
    // As mentioned before, I have assumed the data fits in a polynomial equation of order 2
    // Thus 3 linear equations will be found, those are to be solved using Cramer's Rule.

    ll coeff[3][4] =
    {
        { n,       sum_x, sum_x2, sum_y },
        { sum_x,  sum_x2, sum_x3, sum_xy},
        { sum_x2, sum_x3, sum_x4, sum_x2y },
    };

    for(j=1; j<n; j++) // Difference table for Interpolation
    {
        for(i=0; i<n-1; i++)
        {
            grid[i][j] = grid[i+1][j-1] - grid[i][j-1];
        }
    }
    table = grid;
    // Basic Implementation done. Its time for menu program.

    vector<string> menu {"0. Exit", "1. Newton's Interpolation", "2. Newton-Raphson", "3. Curve Fitting"};
    for(i=0; i<menu.size(); i++) cout << "\t" << menu[i] << nn;

    while(t)
    {
        cout << "Enter a choice between 0 to 3 :";
        cin >> t;
        cout << menu[t] << nn;
        if(t == 1) //Interpolation
        {
            cout << "Enter Target Year: ";
            cin >> target_year;
            cout << "Interpolated Population: " << Function(target_year) << nn << nn;
        }
        else if(t == 2) // Root finding (Newton Raphson Method)
        {
            cout << "Enter Target Population: ";
            cin >> target_pop;
            cout << "Estimated Year: " << (ll)Newton_Raphson(target_pop) << nn << nn;
        }
        else if(t == 3) // Least Square Curve fitting
        {
            findSolution(coeff);
            cout << "Nearly matched polynomial: y= " << a0 << " + " << a1 << "x + " << a2 << "x^2\n\n";

            cout << "Error analysis:\n";
            for(i=0; i<n; i++)
            {
                double cal_y = a0 + a1 * xx[i] + a2 * xx[i] * xx[i];
                double er = (double)(abs(cal_y - yy[i])) / yy[i];
                er *= 100;
                cout << xx[i] << " " << yy[i] << " " << cal_y << " " << er << " % " << nn;
            }

            cout << "\n Next 10 years:\n\n";
            cout << "Year       Least Square Fitting        Interpolation\n";
            ll last = xx[n-1] + 1;
            for(i=0; i<10; i++)
            {
                double cal_y = a0 + a1 * last + a2 * last * last;
                cout << last << "         " << cal_y << "         " << Function(last) << nn;
                last++;
            }
            cout << nn;
        }
    }

    return 0;
}

// Real life data (i.e Population) are not followed by mathematical function.
// Thus, Newton's Interpolation Methods cannot provide appropriate answer. It includes huge error.
// But Least Square curve fitting algorithm reduces error and provide result with less error.













