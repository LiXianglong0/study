using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PSOTest2
{
    public class MyPSO
    {
        const int NUM = 100;//粒子数
        const int DIM = 30;//维度
        const double c1 = 1.8; //参数
        const double c2 = 1.8; //参数

        static double xmin = -100;//位置下限
        static double xmax = 100;//位置上限
        static double[] gbestx = new double[DIM];//全局最有接位置
        static double gbestf;//全局优适应度
        static Random rand = new Random();//随机数
        static int cont = 0;
        class particle
        {
            public double[] x = new double[DIM]; //当前位置矢量
            public double[] bestx = new double[DIM]; //历史最优秀位置
            public double f; //当前适应度
            public double bestf; //历史适应度
        }
        static particle[] swarm = new particle[NUM];
        private static double f1(double[] x)
        {
            cont++;
            return x.Sum(o => o * o);
        }
        public void run()
        {
            for (int i = 0; i < DIM; i++)
            {
                gbestx[i] = rand.NextDouble() * (xmax - xmin) + xmin;
            }
            gbestf = double.MaxValue;
            for (int i = 0; i < NUM; i++)
            {
                particle p1 = new particle();
                for (int j = 0; j < DIM; j++)
                {
                    p1.x[j] = rand.NextDouble() * (xmax - xmin) + xmin;
                }
                p1.f = f1(p1.x);
                p1.bestf = double.MaxValue;
                swarm[i] = p1;
            }
            for (int t = 0; t < 5000; t++)
            {
                for (int i = 0; i < NUM; i++)
                {
                    particle p1 = swarm[i];
                    for (int j = 0; j < DIM; j++)
                    {
                        p1.x[j] += c1 * rand.NextDouble() * (p1.bestx[j] - p1.x[j])
                            + c2 * rand.NextDouble() * (gbestx[j] - p1.x[j]);

                    }

                    p1.f = f1(p1.x);
                    if (p1.f < p1.bestf)
                    {
                        p1.x.CopyTo(p1.bestx, 0);
                        p1.bestf = p1.f;
                    }
                    if (p1.f < gbestf)
                    {
                        p1.x.CopyTo(gbestx, 0);
                        for (int j = 0; j < DIM; j++)
                        {
                            p1.x[j] = rand.NextDouble() * (xmax - xmin) + xmin;
                        }
                        gbestf = p1.f;
                    }
                }
            }
            Console.WriteLine("最优适应度：{0},{1}", gbestf, cont);
            foreach (var item in gbestx)
            {
                Console.WriteLine("最优解：{0}", item);
            }
            Console.ReadLine();
        }
    }
}



//    class Program
//    {
//        /* 已知函数 𝑦=𝑓(𝑥1,𝑥2)=𝑥1²+𝑥2² ，其中−10≤𝑥1,𝑥2≤10
//        请用粒子群优化算法求解y的最小值 */

//        static int N = 10;//初始化种群大小
//        int min = -10;
//        int max = 10;//函数自变量取值范围
//        double v_max = 5;
//        double v_min = -5;//速度范围
//        double w = 0.5;//惯量权重
//        double c1 = 2;
//        double c2 = 2;//加速系数
//        Random r = new Random();//随机数
//        double[,] particle = new double[N, 2];//初始化3个粒子
//        double[,] current_v = new double[N, 2];//粒子当前速度
//        double[,] pBest = new double[N, 2];//局部最优解
//        double[] gBest = new double[2];//全局最优解
//        double[] current_Fitness = new double[N];//粒子当前适应度（函数值）
//        double[] particle_local_best_Fitness = new double[N];//粒子局部最优适应度（函数值）
//        double particle_global_Fitness;//粒子全局最优适应度（函数值）
//        public void Initial()//粒子初始化
//        {
//            int i, j;
//            for (i = 0; i < N; i++)
//            {
//                for (j = 0; j < 2; j++)
//                {
//                    particle[i, j] = min + (max - min) * r.NextDouble();//初始化群体
//                    pBest[i, j] = particle[i, j];//将当前最优结果写入局部最优
//                    current_v[i, j] = v_min + 2 * v_max * r.NextDouble(); //初始化粒子的速度
//                }

//            }

//            for (i = 0; i < N; i++) //计算每个粒子的适应度
//            {
//                current_Fitness[i] = Fitness(particle[i, 0], particle[i, 1]);
//                particle_local_best_Fitness[i] = current_Fitness[i];
//            }

//            particle_global_Fitness = particle_local_best_Fitness[0];//找出全局最优的适应度(函数值)
//            j = 0;
//            for (i = 0; i < N; i++)
//            {
//                if (particle_local_best_Fitness[i] < particle_global_Fitness)
//                {
//                    particle_global_Fitness = particle_local_best_Fitness[i];
//                    j = i;
//                }
//            }

//            for (i = 0; i < 2; i++) //更新全局最优向量
//            {
//                gBest[i] = pBest[j, i];
//            }
//            Console.WriteLine("初始化粒子完成!");
//            Console.WriteLine("正在进行迭代...");
//        }

//        public void Renew_location()//迭代更新粒子所在的位置和速度
//        {
//            for (int i = 0; i < N; i++)
//            {
//                for (int j = 0; j < 2; j++)
//                {
//                    double r1 = r.NextDouble();
//                    double r2 = r.NextDouble();
//                    current_v[i, j] += (w * current_v[i, j] + c1 * r1 * (pBest[i, j] - particle[i, j]) +
//                        c2 * r2 * (gBest[j] - particle[i, j]));//更新速度
//                    particle[i, j] += current_v[i, j];//更新位置
//                    if (particle[i, j] > max)//越界的位置，合法性调整(定义域范围内)
//                    {
//                        particle[i, j] = max;
//                    }
//                    if (particle[i, j] < min)
//                    {
//                        particle[i, j] = min;
//                    }
//                }
//            }

//        }
//        static double Fitness(double x1, double x2)            //返回y=f(x1,x2)=x1²+x2²的函数值
//        {
//            return Math.Pow(x1, 2.0) + Math.Pow(x2, 2.0);
//        }


//        public void Renew_Fitness()//评估例子的适应度函数值
//        {
//            int j = -1;
//            for (int i = 0; i < N; i++)
//            {
//                if (Fitness(particle[i, 0], particle[i, 1]) < current_Fitness[i])//更新局部最优解pBest
//                {
//                    pBest[i, 0] = particle[i, 0];
//                    pBest[i, 1] = particle[i, 1];
//                }
//                particle_local_best_Fitness[i] = Fitness(pBest[i, 0], pBest[i, 1]);//更新局部最优适应值

//                if (particle_local_best_Fitness[i] < particle_global_Fitness)//更新全局最优适应度
//                {
//                    particle_global_Fitness = particle_local_best_Fitness[i];
//                    j = i;
//                }
//            }


//        }
//        static void Main(string[] args)
//        {
//            Program PSO = new Program();
//            PSO.Initial();
//            for (int i = 0; i < 1; i++)//迭代10000次
//            {
//                PSO.Renew_location();
//                PSO.Renew_Fitness();
//            }
//            Console.WriteLine("迭代完成，最小值y={0}", PSO.particle_global_Fitness);
//            Console.WriteLine("最优解X1={0}，X2={1}", PSO.gBest[0], PSO.gBest[1]);
//            Console.ReadLine();
//        }
//    }
//}

