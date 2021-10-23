/***************************************************
Code written for the optimization exercises purposes
by Lukasz Sztangret, PhD
Department of Applied Computer Science and Modelling
AGH University of Science and Technology
***************************************************/

#include <iomanip>
#include "opt_alg.h"

int main() {
    try {
        cout << "LAB NUMBER " << LAB_NO << endl;
        cout << "LAB PART " << LAB_PART << endl << endl;
#if LAB_NO == 0

#elif LAB_NO == 1 && LAB_PART == 1
        double t0 = 0, dt = 0.1, tend = 50;
        matrix Y0 = matrix(2, new double[2]{ 0,0 });
        matrix* Y = solve_ode(t0, dt, tend, Y0);
        ofstream S("out.csv");
        matrix OUT = hcat(Y[0], Y[1]);
        S << OUT;
        S.close();

#elif LAB_NO == 1 && LAB_PART == 2

#elif LAB_NO == 2 && LAB_PART == 1
        double x0 = -20, d = 1, alpha = 2, epsilon = 1e-5, gamma = 1e-200;
        int Nmax = 1000;
        double *p = expansion(x0, d, alpha, Nmax);
        std::cout << p[0] << "\t" << p[1] << std::endl;
        solution::clear_calls();

        solution opt_f = fib(p[0], p[1], epsilon);
        cout << opt_f << endl;
        solution::clear_calls();

        solution opt_l = lag(p[0], p[1], epsilon, gamma, Nmax);
        cout << opt_l << endl;

#elif LAB_NO == 2 && LAB_PART == 2

        ofstream res_1;
        res_1.open("res_1.txt", ofstream::out);

        cout << fixed << setprecision(5);

        std::default_random_engine random_engine(static_cast<long unsigned int>(time(nullptr)));
        std::uniform_real_distribution<double> uniform_dist(-100.0, 100.0);

        array<double, 3> alphas = {2,3,4};

        double d = 1, epsilon = 1e-5, gamma = 1e-200;
        int Nmax = 1000;

        for (double alpha : alphas) {
            for (size_t j = 0; j < 100; ++j) {
                double x0 = uniform_dist(random_engine);

                double *p = expansion(x0, d, alpha, Nmax);

                res_1 << x0 << "\t"
                      << p[0] << "\t"
                      << p[1] << "\t"
                      << solution::f_calls << "\t";

                solution::clear_calls();

                matrix ab_F(1, 1, 200);
                solution opt_f = fib(p[0], p[1], epsilon, &ab_F);

                res_1 << opt_f.x() << "\t"
                      << opt_f.y() << "\t"
                      << solution::f_calls << "\t"
                      << "???" << "\t";

                solution::clear_calls();

                matrix ab_L(1, 1, 200);
                solution opt_l = lag(p[0], p[1], epsilon, gamma, Nmax, &ab_L);

                res_1 << opt_l.x() << "\t"
                      << opt_l.y() << "\t"
                      << solution::f_calls << "\t"
                      << "???" << endl;

                solution::clear_calls();
            }
        }
        res_1.close();
#if 0
        matrix ab_F(1, 1, 200);
        solution opt_f = fib(-100, 100, epsilon, &ab_F);
        cout << "fib:" << endl;
        cout << "opt_f = " << opt_f << endl;
        cout << "ab_F = " << ab_F << endl;
        cout << endl;

        matrix ab_L(1, 1, 200);
        solution opt_l = lag(-100, 100, epsilon, gamma, Nmax, &ab_L);
        cout << "lag: " << endl;
        cout << "opt_f = " << opt_l << endl;
        cout << "ab_F = " << ab_L << endl;
#endif
#elif LAB_NO == 2 && LAB_PART == 3
        solution test(0.001);
        test.fit_fun();
        cout << test << endl;
#elif LAB_NO==3 && LAB_PART==1

#elif LAB_NO==3 && LAB_PART==2

#elif LAB_NO==3 && LAB_PART==3

#elif LAB_NO==4 && LAB_PART==1

#elif LAB_NO==4 && LAB_PART==2

#elif LAB_NO==5 && LAB_PART==1

#elif LAB_NO==5 && LAB_PART==2

#elif LAB_NO==5 && LAB_PART==3

#elif LAB_NO==6 && LAB_PART==1

#elif LAB_NO==6 && LAB_PART==2

#elif LAB_NO==7 && LAB_PART==1

#elif LAB_NO==7 && LAB_PART==2

#endif
    }
    catch(const exception &e) {
        cout << e.what() << endl;
    }
//    catch (char *EX_INFO) {
//        cout << EX_INFO << endl;
//    }
//    system("pause");
    return 0;
}
