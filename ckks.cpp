#include <stdio.h>
#include <iomanip>
#include <cmath>
#include <array>
#include <cstddef>
#include <iostream>
#include <random>

constexpr std::size_t herfN = 2;
constexpr std::size_t N = 2 * herfN;
constexpr std::size_t M = 2 * N;
constexpr double PI = 3.14159265358979323846;

using Complex = std::array<double, 2>;

// 加算
void complex_add(const Complex &a, const Complex &b, Complex &result)
{
    double r_r = a[0] + b[0]; // 実部の加算
    double r_i = a[1] + b[1]; // 虚部の加算
    result[0] = r_r;
    result[1] = r_i;
}

// 減算
void complex_sub(const Complex &a, const Complex &b, Complex &result)
{
    double r_r = a[0] - b[0]; // 実部の減算
    double r_i = a[1] - b[1]; // 虚部の減算
    result[0] = r_r;
    result[1] = r_i;
}

// 乗算
void complex_mul(const Complex &a, const Complex &b, Complex &result)
{
    double r_r = a[0] * b[0] - a[1] * b[1]; // 実部の乗算
    double r_i = a[0] * b[1] + a[1] * b[0]; // 虚部の乗算
    result[0] = r_r;
    result[1] = r_i;
}

// 除算
void complex_div(const Complex &a, const Complex &b, Complex &result)
{
    double denominator = b[0] * b[0] + b[1] * b[1];         // 分母
    double r_r = (a[0] * b[0] + a[1] * b[1]) / denominator; // 実部の除算
    double r_i = (a[1] * b[0] - a[0] * b[1]) / denominator; // 虚部の除算
    result[0] = r_r;
    result[1] = r_i;
}

void mod_q(const __int128_t a, const long long q, __int128_t& b)
{
    __int128_t r = a % q;
    if (r < 0){
        r += q;
    }
    b = r;
}

template <std::size_t N>
void poly_mult(const std::array<__int128_t, N> &a, const std::array<__int128_t, N> &b, std::array<__int128_t,N>& c){
    //mult　polys

    std::array<__int128_t,N> temp1 = {};//c_0 ~ c_N-1
    std::array<__int128_t,N> temp2 = {};//c_N ~ c_2N-1

    for(std::size_t i = 0; i < N; ++i){
        for (std::size_t j = 0; j <= i; ++j){
            temp1[i] += a[j] * b[i - j];
        }
    }

    for(std::size_t i = 0; i < N - 1; ++i){
        for(std::size_t j = 0; j < N - i - 1; ++ j){
            temp2[i] += a[j + i +1] * b[N - 1 -j];
        }
    }

    for(std::size_t i = 0; i < N; ++i){
        c[i] = temp1[i] - temp2[i];
    }
}

template <std::size_t herfN>
void expand(const std::array<Complex, herfN> &data,
            std::array<Complex, 2 * herfN> &expanded_data)
{
    // 前半：コピー
    for (std::size_t i = 0; i < herfN; ++i)
    {
        expanded_data[i] = data[i];
        // std::cout << expanded_data[i][0] << " + " << expanded_data[i][1] << "i\n";
    }

    // 後半：反転 + 虚部の符号を反転
    for (std::size_t i = 0; i < herfN; ++i)
    {
        expanded_data[herfN + i][0] = data[herfN - i - 1][0];  // 実部
        expanded_data[herfN + i][1] = -data[herfN - i - 1][1]; // 虚部反転
        // std::cout << expanded_data[herfN+i][0] << " + " << expanded_data[herfN+i][1] << "i\n";
    }
}

template <std::size_t herfN, std::size_t N>
void shorten(const std::array<Complex, N> &input,
             std::array<Complex, herfN> &output)
{
    for (std::size_t i = 0; i < herfN; ++i)
    {
        output[i] = input[i];
    }
}

template <std::size_t N>
void calc_sigmaR(
    const std::array<Complex, N> &inputarray,
    const std::array<std::array<Complex, N>, N> &inputmatrix,
    std::array<Complex, N> &output)
{
    std::array<Complex, N> z_i;

    for (std::size_t i = 0; i < N; ++i)
    {
        Complex dotnum = {0.0, 0.0};
        double norm_sq = 0.0;

        // ノルム
        for (std::size_t j = 0; j < N; ++j)
        {
            const Complex &bji = inputmatrix[j][i];
            norm_sq += bji[0] * bji[0] + bji[1] * bji[1];
        }

        // 内積<bi,zi>
        for (std::size_t j = 0; j < N; ++j)
        {
            Complex conj_bji = {inputmatrix[j][i][0], -inputmatrix[j][i][1]};
            Complex tmp;
            complex_mul(inputarray[j], conj_bji, tmp);
            complex_add(dotnum, tmp, dotnum);
        }

        Complex abs_b = {norm_sq, 0.0};
        complex_div(dotnum, abs_b, z_i[i]);
        // std::cout << z_i[i][0] << " + " << z_i[i][1] << "i\n";
    }

    // z_i を最近傍の整数に丸める（実部のみ）
    for (std::size_t i = 0; i < N; ++i)
    {
        z_i[i][0] = std::round(z_i[i][0]);
        z_i[i][1] = 0.0;
        // std::cout << z_i[i][0] << " + " << z_i[i][1] << "i\n";
    }

    // A^T * z_i
    for (std::size_t i = 0; i < N; ++i)
    {
        Complex sum = {0.0, 0.0};
        for (std::size_t j = 0; j < N; ++j)
        {
            Complex tmp;
            complex_mul(inputmatrix[i][j], z_i[j], tmp);
            complex_add(sum, tmp, sum);
        }
        output[i] = sum;
        // std::cout << output[i][0] << " + " << output[i][1] << "i\n";
    }
}

template <std::size_t N>
void invsigma(const std::array<Complex, N> &inputarray,
              const std::array<std::array<Complex, N>, N> &inputmatrix,
              std::array<Complex, N> &output)
{
    // calculate alpha= A^-1 * z
    for (std::size_t i = 0; i < N; ++i)
    {
        Complex sum = {0.0, 0.0};
        for (std::size_t j = 0; j < N; ++j)
        {
            Complex temp = {0, 0};
            complex_mul(inputmatrix[i][j], inputarray[j], temp);
            complex_add(sum, temp, sum);
        }
        output[i] = sum;

        // 出力の確認（デバッグ用）
        // std::cout << sum[0] << " + " << sum[1] << "i\n";
    }
}

template <std::size_t N>
void sigma(const std::array<std::array<Complex, N>, N> &inputmatrix,
           const std::array<Complex, N> &input,
           std::array<Complex, N> &output)
{
    /*apply sigma. m'(xi),m'(xi^3),m'(xi^5),...*/
    for (size_t i = 0; i < N; ++i)
    {
        Complex base = inputmatrix[i][1]; // is the same as xi^(2i+1)
        Complex sum = {0, 0};
        Complex val = {1.0, 0.0};

        /*m'(xi^{2i+1})*/
        for (size_t j = 0; j < N; ++j)
        {
            Complex temp;
            complex_mul(input[j], val, temp);

            // next val
            complex_add(temp, sum, sum);
            complex_mul(val, base, val);
        }

        output[i] = sum;
    }
}

template <std::size_t herfN, std::size_t N>
void encode(
    double delta,
    const std::array<Complex, herfN> &inputdata,
    const std::array<std::array<Complex, N>, N> &A,
    const std::array<std::array<Complex, N>, N> &invA,
    std::array<__int128_t, N> &output)
{
    static_assert(N == 2 * herfN, "N must be 2 * herfN");

    // Step 1: expand
    std::array<Complex, N> expanded;
    expand<herfN>(inputdata, expanded);

    // Step 2: scale
    for (auto &x : expanded)
    {
        x[0] *= delta;
        x[1] *= delta;
        //std::cout << x[0] << " + " << x[1] << "i\n";
    }

    // Step 3: project to Z^N
    std::array<Complex, N> projected_real;
    calc_sigmaR<N>(expanded, A, projected_real);
    std::cout << "FT" << std::endl;
    for (const auto &v : projected_real)
    {
        std::cout << v[0]; // 実部
        if (v[1] >= 0)
            std::cout << '+';     // 符号を明示
        std::cout << v[1] << ' '; // 虚部
    }

    // Step 4: apply invsigma
    std::array<Complex, N> real_coeffs;
    invsigma<N>(projected_real, invA, real_coeffs);
    std::cout << "invsigma'" << std::endl;
    for (const auto &v : real_coeffs)
    {
        std::cout << v[0]; // 実部
        if (v[1] >= 0)
            std::cout << '+';     // 符号を明示
        std::cout << v[1] << ' '; // 虚部
    }

    // Step 5: round and output
    for (std::size_t i = 0; i < N; ++i)
    {
        output[i] = static_cast<__int128_t>(std::round(real_coeffs[i][0]));
        //std::cout << "+" << output[i] << " ";
    }
    // std::cout << std::endl;
}

template <std::size_t herfN, std::size_t N>
void decode(
    std::array<std::array<Complex, N>, N> A,
    double delta,
    const std::array<__int128_t, N> &input,
    std::array<Complex, herfN> &output)
{
    static_assert(N == 2 * herfN, "N must be 2 * herfN");

    // Step 1: input / delta -> Complex[N] (虚部は0)
    std::array<Complex, N> temp1;
    for (std::size_t i = 0; i < N; ++i)
    {
        temp1[i][0] = static_cast<double>(input[i]) / delta;
        temp1[i][1] = 0.0;
    }

    // Step 2: apply sigma
    std::array<Complex, N> temp2;
    sigma<N>(A, temp1, temp2);

    // Step 3: shorten to herfN
    shorten<herfN, N>(temp2, output);
}

template <std::size_t N>
void HWT(int h, std::array<__int128_t, N> &output)
{
    static std::mt19937 rng{std::random_device{}()};
    static std::uniform_int_distribution<int> dist(0, N - 1);

    // set +-1 and 0
    for (std::size_t i = 0; i < h / 2; ++i)
    {
        output[2 * i] = 1;
        output[2 * i + 1] = -1;
    }
    for (std::size_t i = h; i < N; ++i)
        output[i] = 0;

    for (std::size_t i = 0; i < 1000; ++i)
    {
        int num = dist(rng);
        std::swap(output[0], output[num]);
    }
}

template <std::size_t N>
void DG(std::array<int, N> &output)
{
    static std::mt19937 rng{std::random_device{}()};
    const double sigma = 8.0 / std::sqrt(2.0 * PI); // HE.org recommendation
    const double cutoff = 6.0 * sigma;
    std::normal_distribution<double> dist(0.0, sigma);

    for (std::size_t i = 0; i < N; ++i)
    {
        double x;
        do
        {
            x = dist(rng);
        } while (std::abs(x) > cutoff);

        output[i] = static_cast<int>(std::round(x));
    }
}

template <std::size_t N>
void ZO(std::array<__int128_t, N> &output)
{
    static std::mt19937 rng{std::random_device{}()};
    static std::uniform_int_distribution<int> dist(0, N-1);

    // set +-1 and 0
    for (std::size_t i = 0; i < N / 4; ++i)
    {
        output[2 * i] = 1;
        output[2 * i + 1] = -1;
    }
    for (std::size_t i = N / 2; i < N; ++i)
        output[i] = 0;

    for (std::size_t i = 0; i < 1000; ++i)
    {
        int num = dist(rng);
        std::swap(output[0], output[num]);
    }
}

template <std::size_t N>
void sample_vector(std::array<__int128_t, N> &output, long long qL)
{
    static std::mt19937 rng(std::random_device{}());
    long long min = -(qL >> 1);
    long long max = qL >> 1;
    static std::uniform_int_distribution<long long> dist(min, max-1);

    for (std::size_t i = 0; i < N; ++i)
        output[i] = static_cast<__int128_t>(dist(rng));

    /*std::cout << "a = [";
    for (std::size_t i = 0; i < N; ++i) {
        std::cout << output[i];
        if (i + 1 != N) std::cout << ", ";
    }
    std::cout << "]\n";*/
}

template <std::size_t N>
void KeyGen(int h, std::array<std::array<__int128_t, N>, 2> &pk, std::array<std::array<__int128_t, N>, 2> &sk, long long qL)
{
    std::array<__int128_t, N> s={};
    std::array<int, N> e={};
    std::array<__int128_t,N> a={};
    /*Generate secret key = (1,s)*/
    HWT<N>(h, s);
    for (std::size_t i = 0; i < N; ++i)
    {
        sk[0][i] = 1;
        sk[1][i] = s[i];
    }

    /*std::cout << "sk = [\n";
    for (std::size_t row = 0; row < 2; ++row)
    {
        std::cout << "  [";
        for (std::size_t col = 0; col < N; ++col)
        {
            std::cout << sk[row][col];
            if (col + 1 != N)
                std::cout << ", ";
        }
        std::cout << "]";
        if (row == 0)
            std::cout << ','; // 1 行目の後にカンマ
        std::cout << '\n';
    }
    std::cout << "]\n";
    */

    /*don't keep the error vector*/
    DG<N>(e);

    /*Generate public key = (b,a)*/
    sample_vector<N>(a, qL);

    std::array<__int128_t,N> temp={};
    poly_mult(a,s,temp);
    for (std::size_t i = 0; i < N; ++i)
    {
        pk[1][i] = a[i];
        __int128_t b = -temp[i] + e[i];
        pk[0][i] = b;
    }

    /*
    std::cout << "pk = [\n";
    for (std::size_t row = 0; row < 2; ++row)
    {
        std::cout << "  [";
        for (std::size_t col = 0; col < N; ++col)
        {
            std::cout << pk[row][col];
            if (col + 1 != N)
                std::cout << ", ";
        }
        std::cout << "]";
        if (row == 0)
            std::cout << ','; // 1 行目の後にカンマ
        std::cout << '\n';
    }
    std::cout << "]\n";
    */
}

template <std::size_t N>
void Encrypt(std::array<__int128_t, N> &m, std::array<std::array<__int128_t, N>, 2> &pk, long long qL, std::array<std::array<__int128_t, N>, 2> &output)
{
    /*generate encryption parameters*/
    std::array<__int128_t, N> v;
    std::array<int, N> e_0,e_1;
    ZO(v);
    DG<N>(e_0);
    DG<N>(e_1);

    /*std::cout << "v = [";
    for (std::size_t i = 0; i < N; ++i)
    {
        std::cout << v[i];
        if (i + 1 != N)
            std::cout << ", ";
    }
    std::cout << "]\n";

    std::cout << "e_0 = [";
    for (std::size_t i = 0; i < N; ++i)
    {
        std::cout << e_0[i];
        if (i + 1 != N)
            std::cout << ", ";
    }
    std::cout << "]\n";

    std::cout << "e_1 = [";
    for (std::size_t i = 0; i < N; ++i)
    {
        std::cout << e_1[i];
        if (i + 1 != N)
            std::cout << ", ";
    }
    std::cout << "]\n";
    */
    /*Encryption*/
    /*v · pk + (m + e0, e1) (mod qL)*/
    std::array<__int128_t,N> temp1 = {};
    std::array<__int128_t,N> temp2 = {};
    poly_mult(v,pk[0],temp1);
    poly_mult(v,pk[1],temp2);
    for (std::size_t i = 0; i < N; ++i) {
        mod_q(temp1[i], qL, temp1[i]);
        mod_q(temp2[i], qL, temp2[i]);
    }
    for (std::size_t i = 0; i < N; ++i)
    {
        __int128_t t0 = temp1[i] + m[i] + e_0[i];
        mod_q(t0,qL,output[0][i]);

        __int128_t t1 = temp2[i] + e_1[i];
        mod_q(t1,qL,output[1][i]);
    }
}

template <std::size_t N>
void Decrypt(std::array<std::array<__int128_t, N>, 2> &sk, std::array<std::array<__int128_t, N>, 2> &ct, long long ql, std::array<__int128_t, N> &output)
{
    /*ct・sk (mod ql) = c0 +c1・s (mod ql)*/
    std::array<__int128_t,N> temp={};
    poly_mult(ct[1],sk[1],temp);
    for (std::size_t i = 0; i < N; ++i)
    {
        __int128_t t =  ct[0][i] + temp[i];

        __int128_t r;
        mod_q(t, ql, r);

        if (r > ql / 2) r -= ql;

        output[i] = r;
    }
}

template <std::size_t N>
void generate_vandermonde(double theta,
                          std::array<std::array<Complex, N>, N> &A,
                          std::array<std::array<Complex, N>, N> &InvA)
{
    // Step 1: A の生成
    for (std::size_t i = 0; i < N; ++i)
    {
        double row_theta = theta * (2 * i + 1);
        for (std::size_t j = 0; j < N; ++j)
        {
            A[i][j][0] = std::cos(row_theta * j);
            A[i][j][1] = std::sin(row_theta * j);
        }
    }
    // std::cout << "matrix A" << std::endl;

    // Step 2: InvA を単位行列で初期化
    for (std::size_t i = 0; i < N; ++i)
    {
        for (std::size_t j = 0; j < N; ++j)
        {
            InvA[i][j] = (i == j) ? Complex{1.0, 0.0} : Complex{0.0, 0.0};
        }
    }

    // Step 3: Gauss-Jordan で A を単位行列化、InvA を逆行列に
    std::array<std::array<Complex, N>, N> temp_A = A;

    for (std::size_t i = 0; i < N; ++i)
    {
        // ピボットが0でないか確認
        if (std::abs(temp_A[i][i][0]) < 1e-12 && std::abs(temp_A[i][i][1]) < 1e-12)
        {
            std::cerr << "Error: Singular matrix.\n";
            return;
        }

        // 対角成分を1に正規化
        Complex ratio;
        complex_div({1.0, 0.0}, temp_A[i][i], ratio);

        for (std::size_t j = 0; j < N; ++j)
        {
            complex_mul(temp_A[i][j], ratio, temp_A[i][j]);
            complex_mul(InvA[i][j], ratio, InvA[i][j]);
        }

        // 他の行から i列 を消去
        for (std::size_t k = 0; k < N; ++k)
        {
            if (k == i)
                continue;

            Complex factor;
            complex_mul(temp_A[k][i], {-1.0, 0.0}, factor);

            for (std::size_t j = 0; j < N; ++j)
            {
                Complex t1, t2;
                complex_mul(factor, temp_A[i][j], t1);
                complex_add(temp_A[k][j], t1, temp_A[k][j]);

                complex_mul(factor, InvA[i][j], t2);
                complex_add(InvA[k][j], t2, InvA[k][j]);
            }
        }
    }
}

int main()
{
    double theta = 2 * PI / M;
    double delta =64;//pow(2.0, 40);

    std::array<std::array<Complex, N>, N> A, InvA;
    generate_vandermonde<N>(theta, A, InvA);
    // std::cout << "matrix" << std::endl;
    // std::array<Complex, herfN> data = {{{1.0, 2.0},
    //                                    {3.0, 4.0}}};
    std::array<Complex, herfN> data;
    for (std::size_t i = 0; i < herfN; ++i)
    {
        data[i] = {0, 0};
    }

    data[0] = {1.0, 2.0};
    data[1] = {3.0, 4.0};

    std::cout << std::fixed << std::setprecision(10); // ← ここで設定

    std::cout << "z'" << std::endl;
    for (const auto &v : data)
    {
        std::cout << v[0]; // 実部
        if (v[1] >= 0)
            std::cout << '+';     // 符号を明示
        std::cout << v[1] << ' '; // 虚部
    }
    std::cout << std::endl;

    std::cout << std::defaultfloat;

    std::array<__int128_t, N> encoded, decrypted;
    encode<herfN, N>(delta, data, A, InvA, encoded);

    std::cout << "m" << std::endl;
    for (auto v : encoded)
    {
        std::cout << static_cast<long long>(v) << ' ';
    }
    std::cout << std::endl;

    std::array<std::array<__int128_t,N>,2> ct;
    std::array<std::array<__int128_t, N>, 2> sk, pk;
    long long qL = 36028797019488260; // log q = 55
    int hamming_weight = 2;
    KeyGen<N>(hamming_weight, pk, sk, qL);
    Encrypt<N>(encoded, pk, qL, ct);

    std::cout << "ct = [\n";
    for (std::size_t row = 0; row < 2; ++row)
    {
        std::cout << "  [";
        for (std::size_t col = 0; col < N; ++col)
        {
            std::cout << static_cast<long long>(ct[row][col]);
            if (col + 1 != N)
                std::cout << ", ";
        }
        std::cout << "]";
        if (row == 0)
            std::cout << ','; // 1 行目の後にカンマ
        std::cout << '\n';
    }
    std::cout << "]\n";

    Decrypt<N>(sk, ct, qL, decrypted);
    std::cout << "m'" << std::endl;
    for (auto v : decrypted)
    {
        std::cout << static_cast<long long>(v) << ' ';
    }
    std::cout << std::endl;

    std::array<Complex, herfN> decoded;
    decode<herfN, N>(A, delta, decrypted, decoded);
    std::cout << std::fixed << std::setprecision(15); 

    std::cout << "z'" << std::endl;
    for (const auto &v : decoded)
    {
        std::cout << v[0]; // 実部
        if (v[1] >= 0)
            std::cout << '+';     // 符号を明示
        std::cout << v[1] << ' '; // 虚部
    }
    std::cout << std::endl;

    std::cout << std::defaultfloat;

    return 0;
}