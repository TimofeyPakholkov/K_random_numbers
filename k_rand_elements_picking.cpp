#include <iostream>
#include <random>
#include <string>


double random_summation_normal_dist(const double* const arr, size_t N, size_t K);
double random_summation_residual(const double* const arr, size_t N, size_t K, size_t items_amount, int* multipliers_arr, double accumulated_deviation, double accumulated_average);
double random_summation(const double* const arr, size_t N, size_t K);

int main()
{
    //const size_t const N = 1000;
    //size_t K = 10000;
    const int const N = 1000000;
    int K = 100000000;
    double* arr = new double[N];

    std::default_random_engine re;
    //std::uniform_real_distribution<double> values_dist(-DBL_MAX / 2, DBL_MAX / 2);
    std::uniform_real_distribution<double> values_dist(-INT_MAX / 2, INT_MAX / 2);
    for (size_t i = 0; i < N; i++)
    {
        arr[i] = values_dist(re);
    }

    try
    {
        std::cout << random_summation_normal_dist(arr, N, K) << std::endl;
    }
    catch (std::string error_message)
    {
        std::cout << error_message << std::endl;
    }

    delete[] arr;
    return 0;
}

/*
    !!! Used only in case when K > N !!!
    In case K <= N execution of the task redirects to the function 'random_summation'

    This function calculates sum of K > 0 random values from a given array of size N > 0.
    Iterates over a range of N, on each step generates random multiplier using normal distribuion,
    which depends on ratio K / N. Afterwards redirects to the function 'random_summation_residual'
    to eliminate the difference between K and total amount of items.
*/
double random_summation_normal_dist(const double* const arr, size_t N, size_t K)
{
    if (K == 0 || N == 0) return 0;
    if (K <= N)
    {
        double sum;
        try
        {
            sum = random_summation(arr, N, K);
        }
        catch (std::string error_message)
        {
            throw error_message;
        }
        return sum;
    }

    std::random_device dev;
    std::mt19937 rand_gen(dev());
    double ratio = (double)K / N;
    std::normal_distribution<double>normal_dist(ratio, sqrt(ratio));
    //std::cout << ratio << '\t' << sqrt(ratio) << std::endl;

    size_t items_amount = 0;
    double sum = 0;
    double accumulated_average = 0;
    double accumulated_deviation = 0;
    double current_average = 0;
    int* multipliers_arr = new int[N];

    //std::cout << 0 << ": " << picked_index << "\t" << arr[picked_index] << std::endl;
    for (size_t i = 0; i < N; i++)
    {
        multipliers_arr[i] = round(normal_dist(rand_gen));
        if (multipliers_arr[i] < 0)
        {
            multipliers_arr[i] = 0;
            continue;
        }
        items_amount += multipliers_arr[i];
        //std::cout << i << ": " << picked_index << "\t" << arr[picked_index] << std::endl;
        current_average = (accumulated_average / items_amount) * (items_amount - multipliers_arr[i]) + (arr[i] / items_amount) * multipliers_arr[i];
        accumulated_deviation += (arr[i] - current_average) * multipliers_arr[i] + (accumulated_average - current_average) * (items_amount - multipliers_arr[i]);
        accumulated_average = current_average;
        sum += multipliers_arr[i] * arr[i];
        //std::cout << i << ": " << accumulated_deviation << " + " << items_amount << " * " << accumulated_average << '\t' << arr[i] << std::endl;
    }

    //std::cout << "total items amount: " << items_amount << std::endl;
    //std::cout << accumulated_deviation << " + " << K << " * " << accumulated_average << std::endl;

    try
    {
        sum = random_summation_residual(arr, N, K, items_amount, multipliers_arr, accumulated_deviation, accumulated_average);
    }
    catch (std::string error_message)
    {
        delete[] multipliers_arr;
        throw error_message;
    }
    delete[] multipliers_arr;
    return sum;
}

/*
    This function varies the sum in accordance to difference between K and total amount of items.
    The function handles the case when summation variable is overflowed - 
    in this circumstance throws the exception with error message, 
    containing summation result in form: deviation + average * K.
    If any of those components are overflowed redirects the execution of the task to the function 'random_summation'
*/
double random_summation_residual(const double* const arr, size_t N, size_t K, size_t items_amount, int* multipliers_arr, double accumulated_deviation, double accumulated_average)
{
    std::random_device dev;
    std::mt19937 rand_gen(dev());
    std::uniform_int_distribution<std::mt19937::result_type> indices_dist(0, N - 1);

    double sum = accumulated_deviation + accumulated_average * items_amount;
    size_t picked_index;
    double current_average = 0;
    if (items_amount < K)
    {
        for (size_t i = items_amount; i < K; i++)
        {
            picked_index = indices_dist(rand_gen);
            current_average = (accumulated_average / (i + 1)) * i + arr[picked_index] / (i + 1);
            accumulated_deviation += (arr[picked_index] - current_average) + (accumulated_average - current_average) * i;
            accumulated_average = current_average;
            sum += arr[picked_index];
        }
    }
    if (items_amount > K)
    {
        for (size_t i = items_amount; i > K; i--)
        {
            picked_index = indices_dist(rand_gen);
            if (multipliers_arr[picked_index] == 0)
            {
                i++;
                continue;
            }
            accumulated_deviation -= arr[picked_index] - accumulated_average;
            multipliers_arr[picked_index]--;
            sum -= arr[picked_index];
        }
    }

    //std::cout << sum << '\t' << accumulated_deviation + accumulated_average * K << std::endl;

    std::string sum_str;
    if (isinf(sum) || isnan(sum))
    {
        if (isinf(accumulated_average) || isnan(accumulated_average) || isinf(accumulated_deviation) || isnan(accumulated_deviation))
        {
            try
            {
                sum = random_summation(arr, N, K);
            }
            catch (std::string error_message)
            {
                throw error_message;
            }
            return sum;
        }
        sum_str = std::to_string(accumulated_deviation) + " + " + std::to_string(K) + " * " + std::to_string(accumulated_average);
        throw std::string{ "Variable is overflowed. Sum value:\n" + sum_str };
    }
    else
    {
        return sum;
    }
}

/*
    !!! Recommended to use in case when K <= N !!!
    This function calculates sum of K > 0 random values from a given array of size N > 0
    by generating K random indices and sequentially adding corresponging values.
    The function handles the case when summation variable is overflowed -
    in this circumstance throws the exception with error message,
    containing summation result in form: deviation + average * K.
*/
double random_summation(const double* const arr, size_t N, size_t K)
{
    if (K == 0 || N == 0) return 0;

    std::random_device dev;
    std::mt19937 rand_gen(dev());
    std::uniform_int_distribution<std::mt19937::result_type> indices_dist(0, N - 1);

    size_t picked_index = indices_dist(rand_gen);
    double accumulated_average = arr[picked_index];
    double accumulated_deviation = 0;
    double current_average = arr[picked_index];
    //std::cout << 0 << ": " << picked_index << "\t" << arr[picked_index] << std::endl;
    for (size_t i = 1; i < K; i++)
    {
        picked_index = indices_dist(rand_gen);
        //std::cout << i << ": " << picked_index << "\t" << arr[picked_index] << std::endl;
        current_average = (accumulated_average / (i + 1)) * i + arr[picked_index] / (i + 1);
        accumulated_deviation += (arr[picked_index] - current_average) + (accumulated_average - current_average) * i;
        accumulated_average = current_average;
        //std::cout << i << ": " << accumulated_deviation << " + " << i + 1 << " * " << accumulated_average << std::endl;
    }
    double sum = accumulated_average * K + accumulated_deviation;
    std::string sum_str;
    //std::cout << accumulated_deviation << " + " << K << " * " << accumulated_average << std::endl;
    if (isinf(sum) || isnan(sum))
    {
        sum_str = std::to_string(accumulated_deviation) + " + " + std::to_string(K) + " * " + std::to_string(accumulated_average);
        throw std::string{ "Variable is overflowed. Sum value:\n" + sum_str };
    }
    else
    {
        return sum;
    }
}
