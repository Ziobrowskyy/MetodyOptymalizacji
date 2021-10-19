/***************************************************
Code written for the optimization exercises purposes
by Lukasz Sztangret, PhD
Department of Applied Computer Science and Modelling
AGH University of Science and Technology
***************************************************/

#include "opt_alg.h"

int main()
{
    try
    {
        cout << "LAB NUMBER " << LAB_NO << endl;
        cout << "LAB PART " << LAB_PART << endl << endl;
#if LAB_NO==0

#elif LAB_NO==1 && LAB_PART==1
        double t0 = 0, dt = 0.1, tend = 50;
		matrix Y0 = matrix(2, new double[2]{ 0,0 });
		matrix* Y = solve_ode(t0, dt, tend, Y0);
		ofstream S("out.csv");
		matrix OUT = hcat(Y[0], Y[1]);
		S << OUT;
		S.close();

#elif LAB_NO==1 && LAB_PART==2

#elif LAB_NO==2 && LAB_PART==1

#elif LAB_NO==2 && LAB_PART==2

#elif LAB_NO==2 && LAB_PART==3

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
    catch (char * EX_INFO)
    {
        cout << EX_INFO << endl;
    }
    system("pause");
    return 0;
}
