using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PSOTest2
{
    class Program
    {
        static void Main(string[] args)
        {
            MyPSO pso = new MyPSO();
            pso.run();
            //GeneticAlgorithm gen = new GeneticAlgorithm();
            //gen.Run();
        }
    }
}

