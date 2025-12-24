#include <vector>
#include <random>
#include <chrono>
#include <iostream>
#include <thread>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <set>

#include "../CEC/2010demo/Header.h"
// #include "../CEC/CEC2013/Header.h"

using namespace std;

extern const int DIMENSION = 1000; // 个体维度
// extern const int DIMENSION = 2000; // 个体维度
const int POPULATION_SIZE = 700; // 种群大小
// const int POPULATION_SIZE = 1400;   // 种群大小
int FEs = 3000 * DIMENSION; // 最大评估次数
double fai;
double beita = 0.2;
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); // 随机数种子
mt19937 generator_start(seed);                                               // 随机数生成器
uniform_int_distribution<int> distribution_ind(0, POPULATION_SIZE - 1);      // 均匀分布
uniform_real_distribution<double> distribution_rand(0.0, 1.0);               // 均匀分布
// std::vector<std::array<int, 4>> itercountall;

int count3 = 0;
int count4 = 0;
int count5 = 0;

int count_iter = 0;
double rate_in=0;

// 个体结构体
struct individual
{
    int id;
    vector<double> x; // 个体的位置
    vector<double> v; // 个体的速度
    double fitness;   // 个体适应度
};

vector<individual> temp_ind;
set<int> selected_indices;
vector<int> temp_vec;
int r1, r2, r3, r4, r5;
double randw, rand1, rand2, rand3, rand4, rand5, fai1, fai2, fai3, fai4, fai5;
individual ind3, ind4, ind5;

// vector2array
double *vector2array(vector<double> v)
{
    double *a = new double[v.size()];
    for (int i = 0; i < v.size(); i++)
    {
        a[i] = v[i];
    }
    return a;
}

// 随机初始化种群
void init_population(vector<individual> &population, Benchmarks *benchtool, double X_Max, double X_Min, double &best)
{
    // 设置随机数生成器
    uniform_real_distribution<double> distribution(X_Min, X_Max); // 均匀分布

    for (int j = 0; j < POPULATION_SIZE; j++)
    {
        // 维度
        individual ind;
        vector<double> p(DIMENSION);
        vector<double> v(DIMENSION);
        for (int k = 0; k < DIMENSION; k++)
        {
            p[k] = distribution(generator_start);       // 随机初始化
            v[k] = distribution(generator_start) / 5.0; // 初始速度
        }
        ind.id = j + 1;
        ind.x = p;
        ind.v = v;

        double *z = vector2array(p);
        ind.fitness = benchtool->compute(z);
        free(z);

        if (ind.fitness < best)
        {
            best = ind.fitness;
        }

        population.push_back(ind);
    }
}

// 计算fai
double calculate_fai(double x, int FEs, int numevals)
{
    double fai = 1 / x - (1 / x - (1 / (x + 1))) * ((double)numevals / FEs);
    return fai;
}

// 更新
void optimization(vector<individual> &population, Benchmarks *benchtool, double X_Max, double X_Min, double &best)
{
    double out=0;
    double in=0;
    // 排序，小的在前
    // sort(population.begin(), population.end(), [](individual a, individual b)
    //      { return a.fitness < b.fitness; });

    // 从前往后更新    
    for (int i = 0; i < POPULATION_SIZE; i++)
    // 从后往前更新
    // for (int i = POPULATION_SIZE - 1; i > 0; i--)
    {
        out++;
        // vector<individual> temp_ind;
        // // 随机选5个个体
        // int r1 = distribution_ind(generator_start);
        // int r2 = distribution_ind(generator_start);
        // int r3 = distribution_ind(generator_start);
        // int r4 = distribution_ind(generator_start);
        // int r5 = distribution_ind(generator_start);
        // while (i == r1)
        //     r1 = distribution_ind(generator_start);
        // while (i == r2 || r2 == r1)
        //     r2 = distribution_ind(generator_start);
        // while (i == r3 || r3 == r1 || r3 == r2)
        //     r3 = distribution_ind(generator_start);
        // while (i == r4 || r4 == r1 || r4 == r2 || r4 == r3)
        //     r4 = distribution_ind(generator_start);
        // while (i == r5 || r5 == r1 || r5 == r2 || r5 == r3 || r5 == r4)
        //     r5 = distribution_ind(generator_start);

        // // 5个个体排序
        // int array[5] = {r1, r2, r3, r4, r5};
        // sort(array, array + 5);
        // // sort(array, array + 5, [&](int a, int b)
        // //      { return population[a].fitness < population[b].fitness; });
        // r1 = array[0];
        // r2 = array[1];
        // r3 = array[2];
        // r4 = array[3];
        // r5 = array[4];

        // set<int> selected_indices;
        selected_indices.clear();
        // 使用set确保选择5个不同的个体（排除当前索引i）
        while (selected_indices.size() < 5)
        {
            int random_index = distribution_ind(generator_start);
            if (random_index != i)
            {
                selected_indices.insert(random_index);
            }
        }

        // 将set转换为vector以便排序
        // vector<int> temp_vec(selected_indices.begin(), selected_indices.end());
        temp_vec.clear();
        temp_vec.assign(selected_indices.begin(), selected_indices.end());

        sort(temp_vec.begin(), temp_vec.end(), [&](int a, int b)
             { return population[a].fitness < population[b].fitness; });

        // int r1 = temp_vec[0];
        // int r2 = temp_vec[1];
        // int r3 = temp_vec[2];
        // int r4 = temp_vec[3];
        // int r5 = temp_vec[4];
        r1 = temp_vec[0];
        r2 = temp_vec[1];
        r3 = temp_vec[2];
        r4 = temp_vec[3];
        r5 = temp_vec[4];

        // if (r5 < i)
        if(population[r5].fitness <= population[i].fitness)
        {
            in++;

            temp_ind.clear();
            fai1 = 1;
            fai2 = calculate_fai(2, FEs, benchtool->getNumevals());
            fai3 = calculate_fai(3, FEs, benchtool->getNumevals());
            fai4 = calculate_fai(4, FEs, benchtool->getNumevals());
            fai5 = calculate_fai(5, FEs, benchtool->getNumevals());
            // double randw, rand1, rand2, rand3, rand4, rand5;
            // double fai1 = 1;
            // double fai2 = calculate_fai(2, FEs, benchtool->getNumevals());
            // double fai3 = calculate_fai(3, FEs, benchtool->getNumevals());
            // double fai4 = calculate_fai(4, FEs, benchtool->getNumevals());
            // double fai5 = calculate_fai(5, FEs, benchtool->getNumevals());

            // individual ind2 = population[i];
            // // 更新速度
            // // v=w*v+rand*(r1-x)+rand*fai*(r3-x)
            // for (int j = 0; j < DIMENSION; j++)
            // {
            //     randw = distribution_rand(generator_start);
            //     rand1 = distribution_rand(generator_start);
            //     rand2 = distribution_rand(generator_start);
            //     ind2.v[j] = randw * ind2.v[j] + rand1 * fai1 * (population[r1].x[j] - ind2.x[j]) + rand2 * fai2 * (population[r2].x[j] - ind2.x[j]);
            // }

            // // 更新位置
            // for (int j = 0; j < DIMENSION; j++)
            // {
            //     ind2.x[j] = ind2.x[j] + ind2.v[j];
            //     // 边界处理
            //     if (ind2.x[j] > X_Max)
            //     {
            //         ind2.x[j] = X_Max;
            //     }
            //     else if (ind2.x[j] < X_Min)
            //     {
            //         ind2.x[j] = X_Min;
            //     }
            // }

            // // 更新适应度
            // double *z = vector2array(ind2.x);
            // ind2.fitness = benchtool->compute(z);
            // free(z);

            // // 更好则更新
            // if (ind2.fitness < population[i].fitness)
            // {
            //     population[i] = ind2;
            //     // 更新最优
            //     if (population[i].fitness < best)
            //     {
            //         best = population[i].fitness;
            //     }
            // }
            // else // 如果不好于当前个体,则一个学习样本
            // {
            //     temp_ind.push_back(ind2);
            ind3 = population[i];
            // individual ind3 = population[i];
            // 更新速度
            // v=w*v+rand*(r1-x)+rand*fai*(r2-x)+rand*beita*(r3-x)
            for (int j = 0; j < DIMENSION; j++)
            {
                randw = distribution_rand(generator_start);
                rand1 = distribution_rand(generator_start);
                rand2 = distribution_rand(generator_start);
                rand3 = distribution_rand(generator_start);
                ind3.v[j] = randw * ind3.v[j] + rand1 * fai1 * (population[r1].x[j] - ind3.x[j]) + rand2 * fai2 * (population[r2].x[j] - ind3.x[j]) + rand3 * fai3 * (population[r3].x[j] - ind3.x[j]);
            }

            // 更新位置
            for (int j = 0; j < DIMENSION; j++)
            {
                ind3.x[j] = ind3.x[j] + ind3.v[j];
                // 边界处理
                if (ind3.x[j] > X_Max)
                {
                    ind3.x[j] = X_Max;
                }
                else if (ind3.x[j] < X_Min)
                {
                    ind3.x[j] = X_Min;
                }
            }

            // 更新适应度
            double *z = vector2array(ind3.x);
            ind3.fitness = benchtool->compute(z);
            free(z);

            if (ind3.fitness < population[i].fitness)
            {
                count3++;
                population[i] = ind3;
                // 更新最优
                if (population[i].fitness < best)
                {
                    best = population[i].fitness;
                }
            }
            else // 如果不好于当前个体,则两个学习样本
            {
                temp_ind.push_back(ind3);
                ind4 = population[i];
                // individual ind4 = population[i];
                // 更新速度
                // v=w*v+rand*(r1-x)+rand*fai*(r2-x)+rand*beita*(r3-x)
                for (int j = 0; j < DIMENSION; j++)
                {
                    randw = distribution_rand(generator_start);
                    rand1 = distribution_rand(generator_start);
                    rand2 = distribution_rand(generator_start);
                    rand3 = distribution_rand(generator_start);
                    rand4 = distribution_rand(generator_start);
                    ind4.v[j] = randw * ind4.v[j] + rand1 * fai1 * (population[r1].x[j] - ind4.x[j]) + rand2 * fai2 * (population[r2].x[j] - ind4.x[j]) + rand3 * fai3 * (population[r3].x[j] - ind4.x[j]) + rand4 * fai4 * (population[r4].x[j] - ind4.x[j]);
                }

                // 更新位置
                for (int j = 0; j < DIMENSION; j++)
                {
                    ind4.x[j] = ind4.x[j] + ind4.v[j];
                    // 边界处理
                    if (ind4.x[j] > X_Max)
                    {
                        ind4.x[j] = X_Max;
                    }
                    else if (ind4.x[j] < X_Min)
                    {
                        ind4.x[j] = X_Min;
                    }
                }

                // 更新适应度
                double *z = vector2array(ind4.x);
                ind4.fitness = benchtool->compute(z);
                free(z);

                if (ind4.fitness < population[i].fitness)
                {
                    count4++;
                    population[i] = ind4;
                    // 更新最优
                    if (population[i].fitness < best)
                    {
                        best = population[i].fitness;
                    }
                }
                else // 如果不好于当前个体,则3个学习样本
                {
                    temp_ind.push_back(ind4);
                    ind5 = population[i];
                    // individual ind5 = population[i];
                    // 更新速度
                    // // v=w*v+rand*(r1-x)+rand*fai*(r2-x)+rand*beita*(r3-x)
                    for (int j = 0; j < DIMENSION; j++)
                    {
                        randw = distribution_rand(generator_start);
                        rand1 = distribution_rand(generator_start);
                        rand2 = distribution_rand(generator_start);
                        rand3 = distribution_rand(generator_start);
                        rand4 = distribution_rand(generator_start);
                        rand5 = distribution_rand(generator_start);
                        ind5.v[j] = randw * ind5.v[j] + rand1 * fai1 * (population[r1].x[j] - ind5.x[j]) + rand2 * fai2 * (population[r2].x[j] - ind5.x[j]) + rand3 * fai3 * (population[r3].x[j] - ind5.x[j]) + rand4 * fai4 * (population[r4].x[j] - ind5.x[j]) + rand5 * fai5 * (population[r5].x[j] - ind5.x[j]);
                    }

                    // 更新位置
                    for (int j = 0; j < DIMENSION; j++)
                    {
                        ind5.x[j] = ind5.x[j] + ind5.v[j];
                        // 边界处理
                        if (ind5.x[j] > X_Max)
                        {
                            ind5.x[j] = X_Max;
                        }
                        else if (ind5.x[j] < X_Min)
                        {
                            ind5.x[j] = X_Min;
                        }
                    }

                    // 更新适应度
                    double *z = vector2array(ind5.x);
                    ind5.fitness = benchtool->compute(z);
                    free(z);

                    if (ind5.fitness < population[i].fitness)
                    {
                        count5++;
                        population[i] = ind5;
                        // 更新最优
                        if (population[i].fitness < best)
                        {
                            best = population[i].fitness;
                        }
                    }
                    else // 如果不好于当前个体,则4个学习样本
                    {
                        temp_ind.push_back(ind5);
                        // 比较ind1 2 3 4 5选择更好的
                        sort(temp_ind.begin(), temp_ind.end(), [](individual a, individual b)
                             { return a.fitness < b.fitness; });
                        population[i] = temp_ind[0];
                        // 对应的count++
                        // if (temp_ind[0].fitness == ind3.fitness)
                        // {
                        //     count3++;
                        // }
                        // else if (temp_ind[0].fitness == ind4.fitness)
                        // if (temp_ind[0].fitness == ind4.fitness)
                        // {
                        //     count4++;
                        // }
                        // else if (temp_ind[0].fitness == ind5.fitness)
                        // {
                        //     count5++;
                        // }
                    // }
                }
            }
            }
        }
    }

    rate_in+=in/out;
}

// 迭代器
void iteration(int FID, int thread)
{
    temp_ind.reserve(3);   // 预设容量为3，但size仍为0
    temp_vec.reserve(5);   // 预设容量为5，但size仍为0

    // 记录开始时间
    auto start_time = std::chrono::high_resolution_clock::now();

    // 初始化评估函数
    Benchmarks *benchtool;
    benchtool = generateFuncObj(FID);
    benchtool->setThreadNum(thread);
    // benchtool->setDimension(DIMENSION);
    benchtool->nextRun();

    double X_Max = benchtool->getMaxX();
    double X_Min = -X_Max;

    vector<individual> population;

    // vector<double> x_avg(DIMENSION);
    double best = INFINITY;

    // 初始化种群
    init_population(population, benchtool, X_Max, X_Min, best);
    cout << "F" << FID << "\tinit compilte" << endl;

    // 写出初始化数据，文件在"./result_2010/results_f" + std::to_string(ID) + ".csv"
    ofstream outfile;
    outfile.open("./result_2010/results_f" + std::to_string(FID) + ".csv", ios::app);
    outfile << 0 << ", " << FID << ", " << thread << "," << std::scientific << best << endl;
    outfile.close();

    // 迭代
    while (benchtool->getNumevals() < FEs)
    {
        // optimition
        optimization(population, benchtool, X_Max, X_Min, best);

        count_iter++;
        // cout << "F" << FID << "\tFEs=" << benchtool->getNumevals() << "\tBest:" << best
        //      << endl
        //      << endl
        //      << endl;
    }

    // // 写出初始化数据，文件在"./result_2010/results_f" + std::to_string(ID) + ".csv"
    // ofstream out;
    // out.open("./result_2010/count/results_f" + std::to_string(FID)+"_count" + ".csv", ios::app);
    // out << count3 << ", " << count4 << ", " << count5 << endl;
    // out.close();
    // count3 = 0;
    // count4 = 0;
    // count5 = 0;

    // //输出itercountall,每一个元素是一个数组，包含5个元素，输出一个元素换一行
    // outfile.open("./result_2010/count/f" + std::to_string(FID) + ".csv", ios::app);
    // for (int i = 0; i < itercountall.size(); i++)
    // {
    //     // outfile << itercountall[i][0] << ", " << itercountall[i][1] << ", " << itercountall[i][2] << ", " << itercountall[i][3] << ", " << itercountall[i][4] << endl;
    //     outfile << itercountall[i][0] << ", " << itercountall[i][1] << ", " << itercountall[i][2] << ", " << itercountall[i][3] << endl;
    // }
    // outfile.close();

    free(benchtool);

    // 记录结束时间
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    double total_time_seconds = duration.count() / 1000.0;

    // // 写入时间文件
    // ofstream time_file;
    // time_file.open("./time_2010.csv", ios::app);
    // time_file << FID << ", " << thread << ", " << std::fixed << std::setprecision(6) << total_time_seconds << endl;
    // time_file.close();

    // 写入平均更新率
    double average_rate_in = rate_in / count_iter;
    ofstream rate_file;
    rate_file.open("./rate_2010.csv", ios::app);
    rate_file << FID << ", " << average_rate_in << endl;
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        cout << "Missing input argument." << std::endl;
        return 1;
    }
    int FID = stoi(argv[1]);

    unsigned run = 1; // 每个函数运行次数

    for (int j = 0; j < run; j++)
    {
        iteration(FID, j); // 启动一个线程并加入到向量中
    }
}
