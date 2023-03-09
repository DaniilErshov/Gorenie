// ConsoleApplication14.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

// ConsoleApplication13.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.

#include <iostream>
#include <fstream>
using namespace std;

double Lambda(double T) {
   return 1 + 0.01 * T + pow(10, -5) * pow(T, 2);
   //return 10;
}

double Lambda_prime(double T) {
    return 0.01 / 2. + pow(10, -5) *  T / 2.;
    //return 0;
}

void Jacobian(double** jac, double* T, const int N1, double T_l, double T_r, const double* x, const double x_l, const double x_r)
{
    jac[0][0] = (2. / (x[1] - x_l)) * (Lambda_prime((T[1] + T[0]) / 2.) * (T[1] - T[0]) / (x[1] - x[0]) - Lambda((T[1] + T[0]) / 2.) / (x[1] - x[0])
        - Lambda_prime((T[0] + T_l) / 2.) * (T[0] - T_l) / (x[0] - x_l) - Lambda((T_l + T[0]) / 2.) / (x[0] - x_l));
    jac[0][1] = (2. / (x[1] - x_l)) *  (Lambda_prime((T[1] + T[0]) / 2.) * (T[1] - T[0]) / (x[1] - x[0]) + Lambda((T[1] + T[0]) / 2.) / (x[1] - x[0]));
    jac[0][2] = 0;
    for (int i = 1; i < N1; i++) {
        jac[i][i - 1] = (2. / (x[i + 1] - x[i - 1])) * (- Lambda_prime((T[i] + T[i - 1]) / 2.) * (T[i] - T[i - 1]) / (x[i] - x[i - 1]) 
            + Lambda((T[i] + T[i - 1]) / 2.) / (x[i] - x[i - 1]));
        jac[i][i] = (2. / (x[i + 1] - x[i - 1])) * (Lambda_prime((T[i] + T[i + 1]) / 2.) * (T[i + 1] - T[i]) / (x[i + 1] - x[i])
            - Lambda((T[i] + T[i + 1]) / 2.) / (x[i + 1] - x[i]) - Lambda_prime((T[i] + T[i - 1]) / 2.) * (T[i] - T[i - 1]) / (x[i] - x[i - 1]) 
            - Lambda((T[i] + T[i - 1]) / 2.) / (x[i] - x[i - 1])) ;
        jac[i][i + 1] = (2. / (x[i + 1] - x[i - 1])) *(Lambda_prime((T[i + 1] + T[i]) / 2.) * (T[i + 1] - T[i]) / (x[i + 1] - x[i])
            + Lambda((T[i + 1] + T[i]) / 2.) / (x[i + 1] - x[i]));
    }
    jac[N1][N1 - 2] = 0;
    jac[N1][N1 - 1] = (2. / (x_r - x[N1 - 1])) * (-Lambda_prime((T[N1] + T[N1 - 1]) / 2.) * (T[N1] - T[N1 - 1]) / (x[N1] - x[N1 - 1])
        + Lambda((T[N1] + T[N1 - 1]) / 2.) / (x[N1] - x[N1 - 1]));
    jac[N1][N1] = (2. / (x_r - x[N1 - 1])) * (Lambda_prime((T[N1] + T_r) / 2.) * (T_r - T[N1]) / (x_r - x[N1]) 
        - Lambda((T_r + T[N1]) / 2.) / (x_r - x[N1])
        - Lambda_prime((T[N1] + T[N1 - 1]) / 2.) * (T[N1] - T[N1 - 1]) / (x[N1] - x[N1 - 1]) 
        - Lambda((T[N1] + T[N1 - 1]) / 2.) / (x[N1] - x[N1 - 1]));
}

void F_vector(double* f, double* T, const int N1, double T_l, double T_r, const double* x, const double x_l, const double x_r) {
    f[0] = -(2. / (x[1] - x_l)) * (Lambda((T[1] + T[0]) / 2.) * (T[1] - T[0]) / (x[1] - x[0]) - Lambda((T[0] + T_l) / 2.) * (T[0] - T_l) / (x[0] - x_l));
    for (int j = 1; j < N1; j++)
    {
        f[j] = -(2. / (x[j + 1] - x[j - 1])) * (Lambda((T[j + 1] + T[j]) / 2.) * (T[j + 1] - T[j]) / (x[j + 1] - x[j]) - Lambda((T[j] + T[j - 1]) / 2.) * (T[j] - T[j - 1]) / (x[j] - x[j - 1]));
    }
    f[N1] = -(2. / (x_r - x[N1 - 1])) * (Lambda((T_r + T[N1]) / 2.) * (T_r - T[N1]) / (x_r - x[N1]) - Lambda((T[N1] + T[N1 - 1]) / 2.) * (T[N1] - T[N1 - 1]) / (x[N1] - x[N1 - 1]));
}

void Solve(const int N, const double T_l, const double T_r, const double x_l, const double x_r, ofstream &fout, const double* x)
{
    ofstream file_h;
    file_h.open("file_h.dat");
    const int N1 = N - 1;
    double y;
    double * a, * B, * T_vect, * dT, * f_nevaz;
    a = new double[N]; B = new double[N];
    T_vect = new double[N]; dT = new double[N]; f_nevaz = new double[N];
    double** jac = new double* [N];
    for (int i = 0; i < N; i++) jac[i] = new double[N];
    double f_mod = 100;
    for (int j = 0; j < N; j++) {
        T_vect[j] = T_l + (double)j * (T_r - T_l) / N;
        dT[j] = 0;
    }
    int num_iter = 0;
    file_h << "TITLE=\"" << "Graphics" << "\"" << endl;
    file_h << R"(VARIABLES= "i", "F")" << endl;
    F_vector(f_nevaz, T_vect, N1, T_l, T_r, x, x_l, x_r);
    while (f_mod > pow(10, -6))
    {
        f_mod = 0;
        Jacobian(jac, T_vect, N1, T_l, T_r, x, x_l, x_r);
        cout << "iteration N " << num_iter << endl;
        y = jac[0][0];
        y = jac[0][0];
        a[0] = -jac[0][1] / y;
        B[0] = f_nevaz[0] / y;
        for (int i = 1; i < N1; i++) {
            y = jac[i][i] + jac[i][i - 1] * a[i - 1];
            a[i] = -jac[i][i + 1] / y;
            B[i] = (f_nevaz[i] - jac[i][i - 1] * B[i - 1]) / y;
        }
        dT[N1] = (f_nevaz[N1] - jac[N1][N1 - 1] * B[N1 - 1]) / (jac[N1][N1] + jac[N1][N1 - 1] * a[N1 - 1]);
        T_vect[N1] += dT[N1];
        //std::cout << "T" << N1 + 2 << " =" << T_vect[N1] << " ";
        for (int i = N1 - 1; i >= 0; i--) {
            dT[i] = a[i] * dT[i + 1] + B[i];
            T_vect[i] += dT[i];
            //std::cout << "T" << i + 2 << " =" << T_vect[i] << " ";
        }
        num_iter += 1;
        //std::cout << std::endl;
        F_vector(f_nevaz, T_vect, N1, T_l, T_r, x, x_l, x_r);
        for (int i = 0; i <= N1; i++) {
            f_mod += f_nevaz[i] * f_nevaz[i];
        }
        f_mod = pow(f_mod, 0.5);
        file_h << num_iter << " " << f_mod << endl;
        cout << endl << "modul F = " << f_mod << " ";
        std::cout << std::endl << endl;;
    }
    fout << "TITLE=\"" << "Graphics" << "\"" << endl;
    fout << R"(VARIABLES= "x", "T", "lambda")" << endl;
    fout << 0 << " " << T_l << " " << Lambda(T_l) << endl;
    for (int i = 0; i <= N1; i++) {
        fout << x[i] << " " << T_vect[i] << " " << Lambda(T_vect[i]) << endl;
    }
    fout << x[N1] << " " << T_r << " " << Lambda(T_r);
    std::cout << "Hello World!\n";
    file_h.close();
}

int main()
{
    const int N = 100;
    const double T_l = 1000;
    const double T_r = 5000;
    const double x_l = 0;
    const double x_r = 10;
    double* x;
    double h = (x_r - x_l) * 2. / (1. + N) / N;
    x = new double[N];
    x[0] = x_l + h;
    ofstream fout;
    fout.open("file8.dat");

    for (int i = 1; i < N; i++)
    {
        x[i] = x[i - 1] + h * i;
        cout << "i = " << i << "  x[i] = " << x[i] << endl;

    }
    Solve(N, T_l, T_r, x_l, x_r, fout, x);
    fout.close();
}
// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
