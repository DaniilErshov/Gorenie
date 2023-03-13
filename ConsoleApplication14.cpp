// ConsoleApplication14.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

// ConsoleApplication13.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

double Lambda(double T) {
   return 1 + 0.01 * T + pow(10, -5) * pow(T, 2);
   //return 10;
}

double Lambda_prime(double T) {
    return 0.01 / 2. + pow(10, -5) *  T / 2.;
    //return 0;
}

double FprimeLeft(vector<double>& res_vect, vector<double>& x, const int i)
{
    return (2. / (x[i + 1] - x[i - 1])) * (-Lambda_prime((res_vect[i] + res_vect[i - 1]) / 2.) * (res_vect[i] - res_vect[i - 1]) / (x[i] - x[i - 1])
        + Lambda((res_vect[i] + res_vect[i - 1]) / 2.) / (x[i] - x[i - 1]));
}

double FprimeCenter(vector<double>& res_vect, vector<double>& x, const int i)
{
    return (2. / (x[i + 1] - x[i - 1])) * (Lambda_prime((res_vect[i] + res_vect[i + 1]) / 2.) * (res_vect[i + 1] - res_vect[i]) / (x[i + 1] - x[i])
        - Lambda((res_vect[i] + res_vect[i + 1]) / 2.) / (x[i + 1] - x[i]) - Lambda_prime((res_vect[i] + res_vect[i - 1]) / 2.) * (res_vect[i] - res_vect[i - 1]) / (x[i] - x[i - 1])
        - Lambda((res_vect[i] + res_vect[i - 1]) / 2.) / (x[i] - x[i - 1]));
}

double FprimeRight(vector<double>& res_vect, vector<double>& x, const int i)
{
    return (2. / (x[i + 1] - x[i - 1])) * (Lambda_prime((res_vect[i + 1] + res_vect[i]) / 2.) * (res_vect[i + 1] - res_vect[i]) / (x[i + 1] - x[i])
        + Lambda((res_vect[i + 1] + res_vect[i]) / 2.) / (x[i + 1] - x[i]));
}

void Jacobian(vector<vector<double>>& jac, vector<double>& T, const int N1, vector<double>& x)
{
    jac[1][0] = FprimeCenter(T, x, 1);
    jac[1][1] = FprimeRight(T, x, 1);
    jac[1][2] = 0;
    for (int i = 2; i < N1; i++) {
        jac[i][i - 1] = FprimeLeft(T, x, i);
        jac[i][i] = FprimeCenter(T, x, i);
        jac[i][i + 1] = FprimeRight(T, x, i);
    }
    jac[N1][N1 - 2] = 0;
    jac[N1][N1 - 1] = FprimeLeft(T, x, N1);
    jac[N1][N1] = FprimeCenter(T, x, N1);
}

void F_vector(vector<double>& f, vector<double>& T, const int N1, vector<double>& x) {
    for (int j = 1; j <= N1; j++)
    {
        f[j] = -(2. / (x[j + 1] - x[j - 1])) * (Lambda((T[j + 1] + T[j]) / 2.) * (T[j + 1] - T[j]) / (x[j + 1] - x[j]) - Lambda((T[j] + T[j - 1]) / 2.) * (T[j] - T[j - 1]) / (x[j] - x[j - 1]));
    }

}

void Solve(const int N,  ofstream &fout, vector<double>& x, vector <double>& T_vect)
{
    ofstream file_h;
    file_h.open("file_h.dat");
    const int N1 = N - 1;
    double y;
    vector<double> a(N), B(N), dT(N + 1), f_nevaz(N);
    vector < vector <double> > jac(N);
    for (int i = 0; i < N; i++) {
        jac[i].resize(N);
    }
    double f_mod = 100;
    int num_iter = 0;
    file_h << "TITLE=\"" << "Graphics" << "\"" << endl;
    file_h << R"(VARIABLES= "i", "F")" << endl;
    F_vector(f_nevaz, T_vect, N1, x);
    while (f_mod > pow(10, -6))
    {
        f_mod = 0;
        Jacobian(jac, T_vect, N1, x);
        cout << "iteration N " << num_iter << endl;
        y = jac[1][0];
        a[1] = -jac[1][1] / y;
        B[1] = f_nevaz[1] / y;
        for (int i = 2; i < N1; i++) {
            y = jac[i][i] + jac[i][i - 1] * a[i - 1];
            a[i] = -jac[i][i + 1] / y;
            B[i] = (f_nevaz[i] - jac[i][i - 1] * B[i - 1]) / y;
        }
        dT[N1] = (f_nevaz[N1] - jac[N1][N1 - 1] * B[N1 - 1]) / (jac[N1][N1] + jac[N1][N1 - 1] * a[N1 - 1]);
        T_vect[N1] += dT[N1];
        //std::cout << "T" << N1 + 2 << " =" << T_vect[N1] << " ";
        for (int i = N1 - 1; i >= 1; i--) {
            dT[i] = a[i] * dT[i + 1] + B[i];
            T_vect[i] += dT[i];
            //std::cout << "T" << i + 2 << " =" << T_vect[i] << " ";
        }
        num_iter += 1;
        //std::cout << std::endl;
        F_vector(f_nevaz, T_vect, N1, x);
        for (int i = 1; i < N1; i++) {
            f_mod += f_nevaz[i] * f_nevaz[i];
        }
        f_mod = pow(f_mod, 0.5);
        file_h << num_iter << " " << f_mod << endl;
        cout << endl << "modul F = " << f_mod << " ";
        std::cout << std::endl << endl;;
    }
    fout << "TITLE=\"" << "Graphics" << "\"" << endl;
    fout << R"(VARIABLES= "x", "T", "lambda")" << endl;
    for (int i = 0; i <= N; i++) {
        fout << x[i] << " " << T_vect[i] << " " << Lambda(T_vect[i]) << endl;
    }
    std::cout << "Hello World!\n";
    file_h.close();
}

void Init(const int N, double x_l, double x_r, const double T_l, const double T_r, vector<double> &x, vector<double>& T_vect)
{
    T_vect.resize(N + 1);
    x.resize(N + 1);
    x[0] = x_l;
    double h = (x_r - x_l) * 2. / (N + 1) / (N);
    T_vect[0] = T_l;
    for (int i = 1; i <= N; i++)
    {
        T_vect[i] = T_l + (double)i * (T_r - T_l) / N;
        x[i] = x[i - 1] + h * i;
        cout << "i = " << i << "  x[i] = " << x[i] << endl;

    }

}

int main()
{
    const int N = 50;
    const double T_l = 1000;
    const double T_r = 5000;
    const double x_l = 0;
    const double x_r = 10;
    vector<double> x;
    vector<double> T_vect;
    Init(N, x_l, x_r, T_l, T_r, x, T_vect);
    ofstream fout;
    fout.open("file10.dat");
    Solve(N,fout, x, T_vect);
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
