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
#define tab1 1
#define tab2 1
        double epsilon = 1e-5, gamma = 1e-200;
        int Nmax = 1000;

#if tab1
        ofstream res_1;
        res_1.open("res_1.txt", ofstream::out);

        cout << fixed << setprecision(5);

        std::default_random_engine random_engine(static_cast<long unsigned int>(time(nullptr)));
        std::uniform_real_distribution<double> uniform_dist(-100.0, 100.0);

        array<double, 3> alphas = {2, 3, 4};
        double d = 1;

        for (double alpha: alphas) {
            for (size_t j = 0; j < 100; ++j) {
                double x0 = uniform_dist(random_engine);

                double *p = expansion(x0, d, alpha, Nmax);

                res_1 << x0 << "\t" << p[0] << "\t" << p[1] << "\t" << solution::f_calls << "\t";
                solution::clear_calls();

                matrix ab_F(1, 1, 200);
                solution opt_f = fib(p[0], p[1], epsilon, &ab_F);

                res_1 << opt_f.x() << "\t" << opt_f.y() << "\t" << solution::f_calls << "\t" << "???" << "\t";
                solution::clear_calls();

                matrix ab_L(1, 1, 200);
                solution opt_l = lag(p[0], p[1], epsilon, gamma, Nmax, &ab_L);

                res_1 << opt_l.x() << "\t" << opt_l.y() << "\t" << solution::f_calls << "\t" << "???" << endl;
                solution::clear_calls();

            }
        }
        res_1.close();
#endif
#if tab2
        matrix ab_F(1, 1, 200);
        solution opt_f = fib(-100, 100, epsilon, &ab_F);
        cout << "fib:" << endl;
        cout << "opt_f = " << opt_f << endl;
        cout << "ab_F = " << endl << ab_F << endl;
        cout << endl;
        solution::clear_calls();

        matrix ab_L(1, 1, 200);
        solution opt_l = lag(-100, 100, epsilon, gamma, Nmax, &ab_L);
        cout << "lag: " << endl;
        cout << "opt_f = " << opt_l << endl;
        cout << "ab_F = " << endl << ab_L << endl;
        solution::clear_calls();
#endif

#elif LAB_NO == 2 && LAB_PART == 3
        cout << "fib\n";
        solution opt_f = fib(0.0001, 0.01, 1e-5);
        cout << opt_f << '\n';

        solution::clear_calls();

        cout << "lag\n";
        solution opt_g = lag(0.0001, 0.01, 1e-5,1e-200,300);
        cout << opt_g << '\n';

        solution::clear_calls();

        solution sim_f(0.00241025);
        sim_f.fit_fun();
        cout << sim_f;

        solution::clear_calls();

        solution sim_l(0.00241257);
        sim_l.fit_fun();
        cout << sim_l;
#elif LAB_NO == 3 && LAB_PART == 1
        double s = 0.1, epsilon = 1e-3, alpha_HJ = 0.5, alpha_Rosen = 2, beta = 0.5;
        int Nmax = 5000;
        matrix x0 = 2 * rand_mat(2,1) -1; //lowosa macierz z przedzialu 0,1 * 2 - 1 = los matrix z przeizialu -1, 1
        matrix s0(2,1,5);
        cout << x0 << endl << endl;
        solution opt_HJ = HJ(x0, s, alpha_HJ, epsilon, Nmax);
        cout << "Opt_HJ = " << endl;
        cout << opt_HJ << endl;
        solution::clear_calls();

        solution opt_R = Rosen(x0, s0, alpha_Rosen, beta, epsilon, Nmax);
        cout << "Opt_R = " << endl;
        cout << opt_R << endl;
        solution::clear_calls();


#elif LAB_NO == 3 && LAB_PART == 2
        double s = 0.1, epsilon = 1e-3, alpha_HJ = 0.5, alpha_Rosen = 2, beta = 0.5;
        int Nmax = 5000;
        matrix x0 = 2 * rand_mat(2, 1) - 1; //lowosa macierz z przedzialu 0,1 * 2 - 1 = los matrix z przeizialu -1, 1
        matrix s0(2, 1, 5);
        cout << x0 << endl << endl;
        matrix Xs_HJ = trans(x0);
        matrix Xs_R = trans(x0);
        solution opt_HJ = HJ(x0, s, alpha_HJ, epsilon, Nmax, &Xs_HJ);
        cout << "Opt_HJ = " << endl;
        cout << opt_HJ << endl;
        cout << "Xs_HJ = " << endl;
        cout << Xs_HJ << endl;
        solution::clear_calls();

        solution opt_R = Rosen(x0, s0, alpha_Rosen, beta, epsilon, Nmax, &Xs_R);
        cout << "Opt_R = " << endl;
        cout << opt_R << endl;
        cout << "Xs_R = " << endl;
        cout << Xs_R << endl;
        solution::clear_calls();
#elif LAB_NO == 3 && LAB_PART == 3
        // X = [K1 / K2], y = Q(k1,k2) = calka od 0 t_end ... dt
        // x = [0,10], t0 = 0s, t_end = 100s, dt = 0.1s
        matrix a(2,1,1);
        solution s(a);
        s.fit_fun();
        cout << s << endl;
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
    catch (const exception &e) {
        cout << e.what() << endl;
    }
//    catch (char *EX_INFO) {
//        cout << EX_INFO << endl;
//    }
//    system("pause");
    return 0;
}
