#include"opt_alg.h"

#if LAB_NO > 1

double *expansion(double x0, double d, double alpha, int Nmax, matrix *ud, matrix *ad) {
    auto *p = new double[2];
    solution X0(x0), X1(x0 + d);
    X0.fit_fun(ud, ad);
    X1.fit_fun(ud, ad);
    if (X0.y == X1.y) {
        p[0] = X0.x(); // Operator () - zwraca wartosc macierzy
        p[1] = X1.x();
        return p;
    }
    if (X0.y < X1.y) {
        d *= -1;
        X1.x = x0 + d;
        X1.fit_fun(ud, ad);
        if (X0.y <= X1.y) {
            p[0] = X1.x();
            p[1] = X0.x() - d; //bo d jest ujemne dlatego minus
            return p;
        }
    }
    solution X2;
    int i = 1;
    while (true) {
        X2.x = x0 + pow(alpha, i) * d; //x0 - punkt startowy
        X2.fit_fun(ud, ad);
        if (X1.y <= X2.y || solution::f_calls > Nmax)
            break;
        X0 = X1;
        X1 = X2;
        ++i;
    }
    d > 0 ? (p[0] = X0.x(), p[1] = X2.x()) : (p[0] = X2.x(), p[1] = X0.x());
    return p;
}

/*
 *   if (d > 0) {
        p[0] = X0.x();
        p[1] = X2.x();
    } else {
        p[0] = X2.x();
        p[1] = X0.x();
    }
 */
solution fib(double a, double b, double epsilon, matrix *ud, matrix *ad) {
    int n = static_cast<int>(
            ceil(log2(sqrt(5) * (b - a) / epsilon) / log2((1 + sqrt(5)) * 0.5))
    );
    int *F = new int[n]{1, 1};
    for (int i = 2; i < n; ++i)
        F[i] = F[i - 1] + F[i - 2];
    solution A(a), B(b), C, D;
    //mnozymy przez długość przedzialu
    C.x = B.x - 1.0 * F[n - 2] / F[n - 1] * (B.x - A.x);
    D.x = A.x + B.x - C.x;
    C.fit_fun(ud, ad);
    D.fit_fun(ud, ad);
    for (int i = 0; i <= n - 3; ++i) {
        if (C.y < D.y)
            B = D;
        else
            A = C;
        C.x = B.x - 1.0 * F[n - i - 2] / F[n - i - 1] * (B.x - A.x);
        D.x = A.x + B.x - C.x;
        C.fit_fun(ud, ad);
        D.fit_fun(ud, ad);

#if LAB_NO == 2 && LAB_PART == 2
        ud->add_row((B.x- A.x)());
#endif
    }
    return C;
}

solution lag(double a, double b, double epsilon, double gamma, int Nmax, matrix *ud, matrix *ad) {
    solution A(a), B(b), C(0.5 * (a + b)), D, D_old(a);
    A.fit_fun(ud, ad);
    B.fit_fun(ud, ad);
    C.fit_fun(ud, ad);
    double l, m;

    while (true) {
        l = A.y(0) * (pow(B.x(0), 2) - pow(C.x(0), 2)) + B.y(0) * (pow(C.x(0), 2) - pow(A.x(0), 2)) +
            C.y(0) * (pow(A.x(0), 2) - pow(B.x(0), 2));
        m = A.y(0) * (B.x(0) - C.x(0)) + B.y(0) * (C.x(0) - A.x(0)) + C.y(0) * (A.x(0) - B.x(0));
        if (m <= 0) {
            C.x = NAN;
            C.y = NAN;
            return C;
        }
        D.x = 0.5 * l / m;
        D.fit_fun(ud, ad);
        if (A.x <= D.x && D.x <= C.x) {
            if (D.y < C.y) {
                B = C;
                C = D;
            } else
                A = D;
        } else if (C.x <= D.x && D.x <= B.x) {
            if (D.y < C.y) {
                A = C;
                C = D;
            } else
                B = D;
        } else {
            C.x = NAN;
            C.y = NAN;
            return C;
        }
#if LAB_NO == 2 && LAB_PART == 2
        ud->add_row((B.x - A.x)());
#endif
        if (B.x - A.x < epsilon || abs(D.x() - D_old.x()) < gamma || solution::f_calls > Nmax)
            return C;
        D_old = D;
    }
}

#endif

#if LAB_NO > 2

solution HJ(matrix x0, double s, double alpha, double epsilon, int Nmax, matrix *ud, matrix *ad) {
    //xb - pnkt bazowy, old - stara baza, x
    solution XB(x0), XB_old, X;
    /* f: R2 => R
     * przesuwamy sie wzdluz x1 i x2,
     * jesli krok poprawia wartosc to ok, jest nie to w druga strone
     * - przesuniecie w pierwszej plaszczyznie
     * 1 < 2 (wartosc celu w p1 jest mniejsa w p1 niz w p2)
     * jesli nie to dodajemy 3
     * jesli w 3 < 1 - przesuwamy sie do 3
     * - potem przesuniecie w drugiej plaszczyznie
     * - jesli kroki sa gorsze to sie nie przesuwamy
     * - kazdy z krokow zwraca jeden w 9 punktów otaczajacych p startowy
     * p4    #   #
     * p3    p1  p2
     * p5    #   #
    */
    XB.fit_fun(ud, ad);
    while (true) {

        // etap probny moze sie sie zakonczyc zwroceniem tego samego lub lepszego punktu
        // jesli zakoczy sie porazka to zmniejszamy dl kroku (s)
        // musimy szukac punktu jescze blizej
        // kryterium jest znalezienie s < epsilon - to znaczy ze wszedzie jest gorzej
        // jesli f(xB) < f(x)
        //jeden z tych 9 pkt otaczahjacych (8 obok lub xb startowe)
        X = HJ_trial(XB, s, ud, ad);

        if (X.y < XB.y) { //zakonczylo sie sukcesem
            //etap roboczy
            while (true) {
                XB_old = XB;
                XB = X;
#if LAB_NO == 3 && LAB_PART == 2
                ud->add_row(trans(XB.x));
#endif
                X.x = 2 * XB.x - XB_old.x; //odbicie symetryczne
                X.fit_fun(ud, ad);
                //ocena nowego punktu po odbicuu przesuwamy go
                X = HJ_trial(X, s, ud, ad);
                if (X.y >= XB.y) //nowy X jest gorszy, przerywamy etap roboczy
                    break;
//                if (s < epsilon)
                if (solution::f_calls > Nmax)
                    return XB;
            }
        } else //nie zwrocilo lepszego punktu
            s *= alpha; //zmniejszenie kroku
        if (s < epsilon || solution::f_calls > Nmax) //sprawdz czy krok nie jest mnejszy od kroku granicznego
            return XB;
    }
}

// s - dlugosc kroku, ud (user data), ad (alg data) - chwoliwowo niewykorzystywane
solution HJ_trial(solution XB, double s, matrix *ud, matrix *ad) {
    int n = get_dim(XB); //sprawdzamy rozmiar problemu u nas 2D
    matrix D = ident_mat(n); //macierz zawierajaca kierunki d1 = [1/0] d2 = [0/1], D = [1/0 0/1]
    /* D:
     * [1 0]
     * [0 1]
     */
    solution X;
    for (int i = 0; i < n; ++i) { //dla kazdego z kierunkow
        X.x = XB.x + s * D[i]; //D[i] - i-ty kierunek
        X.fit_fun(ud, ad); //wartosc funkcji celu w tym punkcie
        if (X.y < XB.y) { // jesli pkt jest lepszy to przesuwamy
            XB = X;
        } else { // jesli nie jest lepszy to:
            X.x = XB.x - s * D[i]; //krok w drugim kierunku
            X.fit_fun(ud, ad);
            if (X.y < XB.y) //znowu sprawdzawmy czy lepszy
                XB = X;
        }
    }
    return XB;
}

/*
 * bezgradientowa metoda
 * wystepuje jedynie etap probny - nie ma odbic symetrycznych - poruszamy sie jedynie wzdluz osi
 * zmieniamy kierunki w etapie probnym w tej metodzie w przeciwienstwie do metoty HJ
 * f: R2 -> R - dwa kierunki x1 i x2
 * zaczynamy od x0
 * zaczynamy od bazy kierunkow D = [1/0 0/1]
 * wykonujemy krok i sprawdzamy czy jest lepiej w pierwszym a potem w drugim kierunku
 *         x2
 *         ^    //jesli x2 lepszy od x1
 *         |
 * x0 ->  x1 //jesli x1 lepszy od x0
 * jesli jest gorzej to zostajemy a nie robimy kroku w drugim kierunku
 * s[0] = [1/1]
 * s[1] = [2/2] // jesli obydwa kroku byly lepsze w obydwu kierunkach, zakladajac ze alpha = 2
 * w nastepnej interacji wydluzamy dl kroku jesli jest lepiej
 * w tej metodzie kazdy kierunke ma wlasna dl kroku, s startowe = [1/1]
 * jesli krok jest nieudany to zmieniamy krok przemnazajac przez -beta, nie ruszajac sie
 * s[2] = [-1/???] // s[i] = -beta * s(i), jesli beta = 0.5
 *
 * alg:
 * x(i) -> x(i) + s(i) * d(i)
 * jesli f(xi) > f(x(i) + s(i) * d(i))
 *      x(x+1) = x(i) + s(i) * d(i) // przesuwm sie
 *      s(i+1) = alpha * s(i) // zwiekszam dl kroku
 *      L(i + 1 )  = L(i) + s(i) //lambda - sumaryczne przesuniecie
 * else
 *      x(i+1) = x(i) // nie przesuwam sie
 *      s(i+1) = -beta * s(i) //zmniejszam dl kroku
 *      p(i+1) = p(i) + 1 //ile razy byl nieudany krok w danym kierunku - p - licznik porazek
 */
// s0 - macierz dl kroku dla kazego kierunku
// alpha/beta - zmniejszanie/zwiekszanie dl kroku
#if 1

solution Rosen(matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix *ud, matrix *ad) {
    solution X(x0), Xt; //pkt startowy w x0 , xt - pkt tymaczasowy do sprawdzania nowego kroku
    int n = get_dim(X); //wymiar problemu u nas = 2
    // l - macierz lamba, sumaryczne przesuniecie
    // p - wartosci (0, 1, 2, 3) - ile porazek w danym kierunku
    // s - macierz dl krokow, poniewaz koniczne bedzie wracanie do poczatkowej dl kroku
    // d - macierz kierunkow - poczatkowa macierz jednostrowa D = [1/0 0/1]
    matrix l(n, 1), p(n, 1), s(s0), D = ident_mat(n);
    X.fit_fun(ud, ad);
    while (true) {
        // w kazdym kierunku
        for (int i = 0; i < n; ++i) {
            Xt.x = X.x + s(i) * D[i]; // nowy krok do sprawdzenia [wektor + skalar * wektror]
            Xt.fit_fun(ud, ad);
            if (Xt.y < X.y) { // pkt jest lepszy to sie przesuwamy
                X = Xt;
                l(i) += s(i); // dodanie do lamddy dl kroku
                s(i) *= alpha; //zwiekszenie dl kroku
            } else { // jesli krok jest nieudany
                ++p(i);
                s(i) *= -beta; //zmniejszenie dl kroku
            }
        }
#if LAB_NO == 3 && LAB_PART == 2
        ud->add_row(trans(X.x));
#endif
        /*
         * zmiana bazy kierunkow, czyli zmiana macierzy D
         * obrot wykonujemy jesli w kazdym kierunku jest kazdym kierunku jest porazka
         * wszyskite lambdy sa rozne od 0 - nowe kierunki generujemy z wykorzystaniem przesuniec w kazdym z kieruno
         * kierunek nie mozem byc zerowy !!!
         */
        bool change = true; // czy zrobic obrot
        for (int i = 0; i < n; ++i) {
            if (l(i) == 0 || p(i) == 0) { // ktores sie zeruje
                change = false;
                break;
            }
        }
        /*
         * nowe kierunki zaleza od lambd
         * nowe kierunki sa prostopadle i wektory te musza byc znormalizowane
         * uzywamy ortonormalizacji Grahma-Schmidta
         * 1. wyzmaczamy macierz Q = D * [l1  0  0   0 ...   0]
         *                               [l2  l2 0   0 ...   0]
         *                               [l3  l3 l3  0 ...   0]
         *                               [...                 ]
         *                               [ln  ln ln  ln ... ln]
         *
         * Vj - V1 = Q[1]
         * dj = Vj/ || Vj || // dzielimy przez norme z Vj
         * Vj = Q[j] - E (k =1, k -1) Q[j]^T * d_k * d_k
         */
        if (change) {
            matrix Q(n, n), v(n, 1); // macierz kwadratowa oraz wektor pionowy
            for (int i = 0; i < n; ++i) //wypelniay macierz Q lambdami
                for (int j = 0; j <= i; ++j) {
                    Q(i, j) = l(i); // ity wiersz w Q to ...
                }
            Q = D * Q; // mnozenie przez kierunek
            v = Q[0] / norm(Q[0]);
            D.set_col(v, 0); // mamy ustawony pierwszy kierunek
            // liczymy kolejne kierunki
            for (int i = 1; i < n; ++i) {
                matrix temp(n, 1);
                for (int j = 0; j < i; ++j) {
                    temp = temp + trans(Q[i]) * D[j] * D[j];
                }
                v = (Q[i] - temp) / norm(Q[i] - temp);
                D.set_col(v, i);
            } // koniec obrotu, mamy wszystkie kierunki w macierzy D
            s = s0;
            l = matrix(n, 1); //zerujemy lambdy
            p = matrix(n, 1); // wektory porazek rowniez
        }
        double max_s = abs(s(0));
        //znajdz najmniejsza dl kroku, uwzgledniajac fakt ze dl kroku moga byc ujemne dlatego abs(s(i))
        for (int i = 1; i < n; ++i)
            if (max_s < abs(s(i)))
                max_s = abs(s(i));
        if (max_s < epsilon || solution::f_calls > Nmax)
            return X;
    }
}

#else
solution Rosen(matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
    // Xt - tymczasowy
    solution X(x0), Xt;

    // wymiar problemu
    int n = get_dim(X);

    // l - macierz lambda (sumaryczne przesuniecie dla kazdego pierunku), p - wektor porazek dla kazdego kierunku,
    // s - macierz dlugosci kroku (kazdy kierunek ma swoja), D - kierunki (poczatkowo zgodne z osiami)
    matrix l(n, 1), p(n, 1), s(s0), D = ident_mat(n);
    X.fit_fun(ud, ad);
    while (true)
    {
        // dla kazdego kierunku
        for (int i = 0; i < n; ++i)
        {
            // s(i) - liczba, D[i] - kolumna
            Xt.x = X.x + s(i) * D[i];

            // liczmy wartosc f celu
            Xt.fit_fun(ud, ad);

            // punkt Xt jest lepszy
            if (Xt.y < X.y)
            {

                // przesuwamy sie
                X = Xt;

                // przesuwamy sie o s, lambda - suma dlugosci przesuniec -> zwiekszamy ja o s(i)
                l(i) += s(i);

                // zwiekszamy dlugosc kroku
                s(i) *= alpha;
            }

                // krok nieudany
            else
            {
                // zwiekszamy licznik porazek
                ++p(i);
                // zmniejszamy krok i zmieniamy kierunek
                s(i) *= -beta;
            }
        }
#if LAB_NO==3 && LAB_PART==2
        (*ud).add_row(trans(X.x));
#endif
        // zmiany kierunkow
        // obrot wykonujemy gdy w kazdym kierunku jest przesuniecie i w kazdym kierunku jest porazka
        // p(i) > 0 i l(i) > 0
        bool change = true;
        for (int i = 0; i < n; ++i)
            if (l(i) == 0 || p(i) == 0)
            {
                change = false;
                break;
            }

        // wykonujemy obrot
        if (change)
        {
            // elementy sa wyzerowane
            matrix Q(n,n), v(n,1);

            // wypleniamy macierz q lambdami
            for (int i = 0; i<n; ++i)
                for (int j = 0; j<=i; ++j)
                    Q(i, j) = l(i);

            // teraz q jest macierza trojkatna

            // wyznaczamy nowe Q na podstawie starych kierunkow
            Q = D * Q;

            // pierwszy kierunek
            v = Q[1] / norm(Q[0]);

            // wstawiamy nowa kolumne -> pierwszy kierunek zostal zmieniony
            D.set_col(v,0);

            // kazdy kolejny kierunek wyznczamy tu w petli
            for (int i = 1; i<n; ++i)
            {
                // elementy wyzerowane
                matrix temp(n,1);

                for (int j = 0; j < i; ++j)
                    temp = temp + trans(Q[i]) * D[j] * D[j];

                // normalizacja
                v = (Q[i] - temp) / norm(Q[i] - temp);

                // wstawiamy nowy kierunek
                D.set_col(v,i);
            }
            // poczatkowa dlugosc kroku
            s = s0;
            // zerujemy lambdy
            l = matrix(n, 1);
            // zerujemy porazki
            p = matrix(n, 1);
        }

        // szukamy modulu najdluzeszego kroku (bo dlugosc kroku moze byc ujemna)
        double max_s = abs(s(0));
        for (int i = 1; i < n; ++i)
            if (max_s < abs(s(i)))
                max_s = abs(s(i));
        if (max_s < epsilon || solution::f_calls > Nmax)
            return X;
    }
}
#endif
#endif

#if LAB_NO > 3
/*
 * f(x) = ...
 * F_i = f(x) + c_i + S(x) - S - funkcja kary
 *
 * c_i+1 = dc * c_i - dc > 1 - kara zew, dc < 1 - kara wew
 *
 * dla dc > 1 - kara zew
 *  s(x) = Suma od i = 1, do 3 =  max(0, g_i(x))^2 - zawsze wieksze od 0
 *  dla dc < 1 - kara wew
 *  s(x) = - suma od 1 do 3,  = 1 / g_i
 */
// c0 - pocz wartosc c, sily dary, dc - zmiana  kary, reszta standard
solution pen(matrix x0, double c0, double dc, double epsilon, int Nmax, matrix *ud, matrix *ad) {
    double alpha = 1, beta = 0.5, gamma = 2, delta = 0.5, s = 0.5;
    solution X(x0), X1;
    //dane przekazywane do funkcji, aby wiedziala ktora funkcje kary wybrac
    matrix c(2, new double[2]{c0, dc});
    while (true) {
        X1 = sym_NM(X.x, s, alpha, beta, gamma, delta, epsilon, Nmax, ud, &c);
        //punkt sie nie zmienia
        if (norm(X.x - X1.x) < epsilon || solution::f_calls > Nmax)
            return X1;
        //zmiana wspolczynnika c
        c(0) *= dc;
        //podmiana punktu
        X = X1;
    }
}

solution
sym_NM(matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix *ud,
       matrix *ad) {
    // sympleks ma n+1 wierzcholkow gdzie n to wymiar problemu - u nas 2d -> trojkat
    int n = get_len(x0); //rozmiar problemu
    // wektory w macierzy d - kierunki poruszania sie po osi
    matrix D = ident_mat(n);
    int N = n + 1;
    // nasz sympleks
    solution *S = new solution[N];
    S[0].x = x0;
    S[0].fit_fun(ud, ad);
    for (int i = 1; i < N; ++i) {
        S[i].x = S[0].x + s * D[i - 1];
        S[i].fit_fun(ud, ad);
    }
    //ropzwiazania r - reflecnted, e - ekspansja , n - zawezone
    solution PR, PE, PN;
    // punkt ciezkosci
    matrix pc;
    int i_min, i_max;
    while (true) {
        // zerowanie indeksow wierzcholkow
        i_min = i_max = 0;
        for (int i = 1; i < N; ++i) {
            if (S[i_min].y < S[i].y)
                i_min = i;
            if (S[i_max].y > S[i].y)
                i_max = i;
        }
        // liczumy srodek ciezkoci przez srednia figuyr
        pc = matrix(n, 1);
        for (int i = 0; i < N; ++i)
            if (i != i_max)
                pc = pc + S[i].x;
        pc = pc / (N - 1);
        // odbicie
        PR.x = pc + alpha * (pc - S[i_max].x);
        // ocena punktu odbitego
        PR.fit_fun(ud, ad);
        //warunek akceptacji odbicia
        if (S[i_min].y <= PR.y && PR.y < S[i_max].y) {
            S[i_max] = PR;
        }
            // kiedy robimly ekspansje
        else if (PR.y < S[i_min].y) {
            PE.x = pc + gamma * (PR.x - pc);
            PE.fit_fun(ud, ad);
            if (PR.y <= PE.y)
//            if (PE.y < PR.y)
                S[i_max] = PR;
            else
                S[i_max] = PE;
        } else {
            //zawezenie
            PN.x = pc + beta * (S[i_max].x - pc);
            PN.fit_fun(ud, ad);
            // warunek akceptacji zawezenia
            if (PN.y < S[i_max].y)
                S[i_max] = PN;
            else {
                // robimy redukcje, ktora polega na przysuniecu wszystkich wierzcholko do wierzchoka
                // najlepszego
                for (int i = 0; i < N; ++i)
                    if (i != i_min) { // wszystkie oprocz najlepszego
                        S[i].x = delta * (S[i].x + S[i_min].x);
                        S[i].fit_fun(ud, ad);
                    }
            }
        }
        // sprawdzanie warunku stopu
        double max_s = norm(S[0].x - S[i_min].x); // odlegl pomiedzy wierz zerowym a minumalnym
        for (int i = 1; i < N; ++i)
            if (max_s < norm(S[i].x - S[i_min].x))
                max_s = norm(S[i].x - S[i_min].x);
        if (max_s < epsilon || solution::f_calls > Nmax)
            return S[i_min];
    }
}

#endif

#if LAB_NO > 4
solution SD(matrix x0, double h0, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
    int n = get_len(x0);
    solution X, X1;
    X.x = x0;
    matrix d(n, 1), *P = new matrix[2];
    solution h;
    double *ab;
    while (true)
    {
        X.grad();
        d = ???
        if (h0<0)
        {
            P[0] = ???
            P[1] = ???
            ab = ???
            h = ???
            X1.x = ???
        }
        else
            X1.x = ???
#if LAB_NO==5 && LAB_PART==2
        ???
#endif
        if (???)
        {
            X1.fit_fun(ud, ad);
            return X1;
        }
        ???
    }
}

solution CG(matrix x0, double h0, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
    int n = get_len(x0);
    solution X, X1;
    X.x = x0;
    matrix d(n, 1), *P = new matrix[2];
    solution h;
    double *ab, beta;
    X.grad();
    d = ???
    while (true)
    {
        if (h0<0)
        {
            P[0] = ???
            P[1] = ???
            ab = ???
            h = ???
            X1.x = ???
        }
        else
            X1.x = ???
#if LAB_NO==5 && LAB_PART==2
        ???
#endif
        if (???)
        {
            X1.fit_fun(ud);
            return X1;
        }
        X1.grad();
        beta = ???
        d = ???
        ???
    }
}

solution Newton(matrix x0, double h0, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
    int n = get_len(x0);
    solution X, X1;
    X.x = x0;
    matrix d(n, 1), *P = new matrix[2];
    solution h;
    double *ab;
    while (true)
    {
        X.grad();
        X.hess();
        d = ???
        if (h0<0)
        {
            P[0] = ???
            P[1] = ???
            ab = ???
            h = ???
            X1.x = ???
        }
        else
            X1.x = ???
#if LAB_NO==5 && LAB_PART==2
        ???
#endif
        if (???)
        {
            X1.fit_fun(ud);
            return X1;
        }
        ???
    }
}

solution golden(double a, double b, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
    double alfa = ???
    solution A, B, C, D;
    A.x = a;
    B.x = b;
    C.x = ???
    C.fit_fun(ud, ad);
    D.x = ???
    D.fit_fun(ud, ad);
    while (true)
    {
        if (???)
        {
            ???
            ???
            ???
            C.fit_fun(ud, ad);
        }
        else
        {
            ???
            ???
            ???
            D.fit_fun(ud, ad);
        }
        if (???)
        {
            A.x = ???
            A.fit_fun(ud, ad);
            return A;
        }
    }
}

#endif

#if LAB_NO > 5
solution Powell(matrix x0, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
    int n = get_len(x0);
    matrix D = ident_mat(n), *A = new matrix[2];
    solution X, P, h;
    X.x = x0;
    double *ab;
    while (true)
    {
        P = ???
        for (int i = 0; ???; ++i)
        {
            A[0] = ???
            A[1] = ???
            ab = ???
            h = ???
            P.x = ???
        }
        if (???)
        {
            P.fit_fun(ud);
            return P;
        }
        for (int i = 0; i < n - 1; ++i)
            D.set_col(???);
        D.set_col(???);
        A[0] = ???
        A[1] = ???
        ab = ???
        h = ???
        X.x = ???
    }
}
#endif

#if LAB_NO > 6
solution EA(int N, matrix limits, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
    solution *P = new solution[mi + lambda];
    solution *Pm = new solution[mi];
    random_device rd;
    default_random_engine gen;
    gen.seed(static_cast<unsigned int>(chrono::system_clock::now().time_since_epoch().count()));
    normal_distribution<double> distr(0.0, 1.0);
    matrix IFF(mi, 1), temp(N, 2);
    double r, s, s_IFF;
    double tau = pow(2 * N, -0.5), tau1 = pow(2 * pow(N, 0.5), -0.5);
    int j_min;
    for (int i = 0; i < mi; ++i)
    {
        P[i].x = matrix(N, 2);
        for (int j = 0; j < N; ++j)
        {
            P[i].x(j, 0) = (limits(j, 1) - limits(j, 0))*rand_mat(1, 1)() + limits(j, 0);
            P[i].x(j, 1) = sigma0(j);
        }
        P[i].fit_fun(ud, ad);
        if (P[i].y < epsilon)
            return P[i];
    }
    while (true)
    {
        s_IFF = 0;
        for (int i = 0; ???; ++i)
        {
            IFF(i) = ???
            s_IFF += ???
        }
        for (int i = 0; ???; ++i)
        {
            r = ???
            s = ???
            for (int j = 0; ???; ++j)
            {
                s += ???
                if (???)
                {
                    P[mi + i] = ???
                    break;
                }
            }
        }
        for (int i = 0; ???; ++i)
        {
            r = ???
            for (int j = 0; ???; ++j)
            {
                ???
                ???
            }
        }
        for (int i = 0; ???; i += 2)
        {
            r = ???
            temp = P[mi + i].x;
            P[mi + i].x = ???
            P[mi + i + 1].x = ???
        }
        for (int i = 0; ???; ++i)
        {
            P[mi + i].fit_fun(ud, ad);
            if (???)
                return P[mi + i];
        }
        for (int i = 0; ???; ++i)
        {
            j_min = 0;
            for (int j = 1; ???; ++j)
                if (???)
                    j_min = j;
            Pm[i] = ???
            P[j_min].y = ???
        }
        for (int i = 0; i < mi; ++i)
            P[i] = Pm[i];
        if (???)
            return P[0];
    }
}
#endif
