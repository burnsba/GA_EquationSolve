using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Timers;

namespace GA
{
    /*
     * Notes:
     * 
     * 2014.07.30
     * Initially had 64 registers, and had slow or no convergence after 50,000 iterations. Changed the number of registers to 4,
     * max blocks to 3, max statements to 4 and have very fast convergence now.
     * 
     * Seems like there should be only a small number of mutations per generation.
     * 
     * Increasing the number of blocks and statements on both initial seed and on NewRand dramatically lowers the number
     * of generations to converge.
     * 
     * It seems convergence requires much more death than I originally imagined.
     * 
     * If the population is too small, there doesn't seem to be enough variety, but if the population is large, each generation
     * takes a long time to calculate.
     * 
     * 2014.07.31
     * Changed to Parallel.ForEach for massive speed increase. Introduced several threading bugs which are now fixed.
     * 
     * Fitness is now the average of a number of runs, instead of a single run. This was also a huge increase in the speed of convergence.
     * 
     * It seems the number of mutations is a magic number. Too few or too many, and it takes too long, if ever, to converge.
     * 
     * Added more death. There's an ApocalypseGeneration count, where everyone dies. Seems to resolve getting out of a sticking point.
     * Also added heroes, when a new bestFitness is found, which are recalled every RecallHeroesGeneration. This also seems to
     * help convergence.
     * 
     * Changed the fitness formula to weight correct output slots much higher than anything else, which is weighted higher than the 
     * difference of remaining outputs, which is weighted higher than the number of statements used. I also changed the mutate function
     * to be more likely to change a number by +/- 10%. I thought constants would be more likely to appear then. Here's a recent simulation
     * showing the results (about 5200 generations):
     * 
     * 
            _expectedOutput[0] = _input[0] + _input[1];
            _expectedOutput[1] = _input[2];
            _expectedOutput[2] = _input[3] + 2121;

            {r0 = 1
            r4 = r0 * in2

            r1 = r2 ^ in1

            r3 = r0 + in3

            r1 = r1 + in0


            r0 = 8
            r2 = r0 << 8
            r2 = r2 + r0

            r0 = r0 * r0

            r2 = r2 | r0


            out1 = r4 | r4

            out0 = r1 | r1
            out2 = r3 + r2


            }
     * which calculates 2120 from ((8 << 8) + 8) | (8 * 8) = 2056 | 64. Not quite what I had in mind.
     * */

    // mutations can occur at the block level
    // // loop counter can increment or decrement
    // + block can be duplicated or deleted
    // + block can swap position with neighbor
    // + a new block can be randomly generated
    // + (creation should be less likely the more there are)

    // mutations can occur at the statment level
    // + statement can be duplicated or deleted
    // + statement can swap position with neighbor
    // + a new statement can be randomly generated
    // + RHS variable index can be incremented or decremented, or const can be incremented or decremented
    // + operator can be "incremented" or "decremented"
    // + (creation should be less likely the more there are)

    // crossover can occur at the block level
    // crossover can occur at the chromosome level

    class Program
    {
        static void Main(string[] args)
        {
            Simulation sim = new Simulation();

            sim.SolveForEquation = (input, expectedOutput) =>
            {
                // solve-for equation
                expectedOutput[0] = input[0] + input[1];
                expectedOutput[1] = input[2];
                expectedOutput[2] = input[3] + 21;
            };

            sim.InitPopulation();

            while (!sim.HasConverged)
            {
                sim.CalculateFitness();
                sim.MoveNext();
            }
        }
    }

    /// <summary>
    /// Thread safe random number provider.
    /// </summary>
    /// <remarks>
    /// http://csharpindepth.com/Articles/Chapter12/Random.aspx
    /// </remarks>
    public static class RandomProvider
    {
        private static int seed = Environment.TickCount;

        private static ThreadLocal<Random> randomWrapper = new ThreadLocal<Random>(() =>
            new Random(Interlocked.Increment(ref seed))
        );

        /// <summary>
        /// Gets a thread safe instance of <see cref="Random/>.
        /// </summary>
        /// <returns>A thread safe instance of <see cref="Random/>.</returns>
        public static Random GetThreadRandom()
        {
            return randomWrapper.Value;
        }
    }

    public static class Utilities
    {
        public static string OpToString(this Op o)
        {
            switch (o)
            {
                case Op.Add: return "+";
                case Op.Subtract: return "-";
                case Op.Multiply: return "*";
                case Op.Divide: return "/";
                case Op.Modulus: return "%";
                case Op.LeftShift: return "<<";
                case Op.RightShift: return ">>";
                case Op.BitwiseAnd: return "&";
                case Op.BitwiseOr: return "|";
                case Op.BitwiseXor: return "^";

                default: return "unknown";
            }
        }

        /// <summary>
        /// Preconditions: Op is not unknown, and no exceptions will be thrown (e.g., divide by zero)
        /// </summary>
        /// <param name="o"></param>
        /// <param name="left"></param>
        /// <param name="right"></param>
        /// <returns></returns>
        public static int OpCalc(Op o, int left, int right)
        {
            switch (o)
            {
                case Op.Add: return left + right;
                case Op.Subtract: return left - right;
                case Op.Multiply: return left * right;
                case Op.Divide: return left / right;
                case Op.Modulus: return left % right;
                case Op.LeftShift: return left << right;
                case Op.RightShift: return left >> right;
                case Op.BitwiseAnd: return left & right;
                case Op.BitwiseOr: return left | right;
                case Op.BitwiseXor: return left ^ right;

                default: return -1;
            }
        }

        public static Op Next(this Op o)
        {
            switch (o)
            {
                case Op.Add: return Op.Subtract;
                case Op.Subtract: return Op.Multiply;
                case Op.Multiply: return Op.Divide;
                case Op.Divide: return Op.Modulus;
                case Op.Modulus: return Op.LeftShift;
                case Op.LeftShift: return Op.RightShift;
                case Op.RightShift: return Op.BitwiseAnd;
                case Op.BitwiseAnd: return Op.BitwiseOr;
                case Op.BitwiseOr: return Op.BitwiseXor;
                case Op.BitwiseXor: return Op.Add;

                default: return Op.Unknown;
            }
        }

        public static Op Prev(this Op o)
        {
            switch (o)
            {
                case Op.Add: return Op.BitwiseXor;
                case Op.Subtract: return Op.Add;
                case Op.Multiply: return Op.Subtract;
                case Op.Divide: return Op.Multiply;
                case Op.Modulus: return Op.Divide;
                case Op.LeftShift: return Op.Modulus;
                case Op.RightShift: return Op.LeftShift;
                case Op.BitwiseAnd: return Op.RightShift;
                case Op.BitwiseOr: return Op.BitwiseAnd;
                case Op.BitwiseXor: return Op.BitwiseOr;

                default: return Op.Unknown;
            }
        }
    }

    public enum Op
    {
        Unknown,

        Add,

        Subtract,

        Multiply,

        Divide,

        Modulus,

        LeftShift,

        RightShift,

        BitwiseAnd,

        BitwiseOr,

        BitwiseXor
    }

    public class Simulation
    {
        public const int MaxRegisters = 5;
        public const int MaxInputs = 4;
        public const int MaxOutputs = 3;

        public const int PopulationSize = 1200;

        // when creatures are ranked by fitness, everyone after this index will die.
        public const int GenerationalDeathPopulationStartIndex = PopulationSize / 4;

        // Creatures ranked by fitness into groups. Number of groups including the first to keep.
        public const int OtherFitnessRanksCount = 5;

        // Number of creatures per fitness rank group to keep alive. The first group is superceded
        // by GenerationalDeathPopulationStartIndex. Example, GenerationalDeathPopulationStartIndex is 250,
        // OtherFitnessCount is 20, all members of the best fitness group up to 250 will be alive for the
        // next generation. Then the first 20 members of the best fitness group will be marked as alive
        // which is redundant.
        public const int OtherFitnessCount = 20;

        // Number of random new creatures (without parents) to generate each generation. This always
        // occurs whether this is an apocalypse generation or not.
        public const int SpontaneousAdditionsPerGeneration = 25;

        // When a creature is created, chance it will have a single parent.
        public const int ChanceToClone = 80;

        public const int MutationsBetweenGenerations = 5;

        public const int MaxChromosomesPerCreature = 3;
        public const int MaxBlocksPerChromosome = 4;
        public const int MaxStatementsPerBlock = 5;

        //public const int ChanceToMutateElement = 30;

        public const int ChanceToCreate = 4;
        public const int ChanceToDelete = 8;

        public const int FitnessAverageCount = 12;

        /// <summary>
        /// Getting an answer correct/wrong in each position is worth more than the difference of other solutions.
        /// </summary>
        public const int CorrectCountResultWeight = 10000;

        /// <summary>
        /// Difference between measured and expected solution weight.
        /// </summary>
        public const int ResultDifferenceWeight = 100;

        /// <summary>
        /// Wait for a best fitness less than _bestMaxFitness this many times in a row before halting.
        /// </summary>
        private const int BestBeforeConvergence = 50;

        /// <summary>
        /// If the bestFitness variable hasn't improved in this many generations, everyone dies.
        /// </summary>
        private const int ApocalypseGeneration = 250;

        /// <summary>
        /// 
        /// </summary>
        private const int RecallHeroesGeneration = ApocalypseGeneration / 2;

        private int _generation = 0;

        private List<Creature> _population = new List<Creature>();

        private const int _hallOfFameMaxLength = PopulationSize / 10;
        private List<Creature> _hallOfFame = new List<Creature>();
        private UInt64 _bestFitness = UInt64.MaxValue;
        private UInt64 _lastBestFitnessGeneration = 0;

        private int _convergenceCount = 0;
        private int _bestMaxFitness = 0;
        private bool _hasConverged = false;

        public int Generation
        {
            get
            {
                return _generation;
            }

            set
            {
                _generation = value;
            }
        }

        public bool HasConverged
        {
            get
            {
                return _hasConverged;
            }

            private set
            {
                _hasConverged = value;
            }
        }

        public static Random Rand
        {
            get
            {
                return RandomProvider.GetThreadRandom();
            }
        }

        public List<Creature> Population
        {
            get
            {
                return _population;
            }
        }

        public Simulation()
        {
            _bestMaxFitness = MaxChromosomesPerCreature * MaxBlocksPerChromosome * MaxStatementsPerBlock;
        }

        public void InitPopulation()
        {
            for (int i=0; i<PopulationSize; i++)
            {
                _population.Add(Creature.NewSeedCreature());
            }
        }

        public Action<int[], int[]> SolveForEquation
        {
            get;
            set;
        }

        public void CalculateFitness()
        {
            Parallel.ForEach(_population, creature =>
            {
                creature.Fitness = int.MaxValue;

                bool success = false;

                int[] registers = new int[MaxRegisters];
                int[] input = new int[MaxInputs];
                int[] expectedOutput = new int[MaxOutputs];
                int[] measuredOutput = new int[MaxOutputs];

                UInt64[] fitness = new UInt64[FitnessAverageCount];

                for (int i = 0; i < FitnessAverageCount; i++)
                {
                    for (int counter = 0; counter < MaxInputs; counter++)
                    {
                        input[counter] = Simulation.Rand.Next(0, 10000);
                    }

                    Array.Clear(registers, 0, registers.Length);
                    Array.Clear(expectedOutput, 0, expectedOutput.Length);
                    Array.Clear(measuredOutput, 0, measuredOutput.Length);

                    SolveForEquation(input, expectedOutput);

                    try
                    {
                        success = creature.Eval(registers, input, measuredOutput);
                    }
                    catch
                    {
                        success = false;
                        creature.Fitness = int.MaxValue;
                        creature.Health = Creature.LifeState.Dead;
                    }

                    if (success)
                    {
                        try
                        {
                            // a solution with 3 incorrect outputs is penalized more than a solution
                            // with 2 correct outputs and 1 incorrect output with the same difference

                            UInt64 weight = 0;
                            UInt64 t = 0;
                            UInt64 solutionFitness = 0;

                            for (int counter = 0; counter < MaxOutputs; counter++)
                            {
                                t = (UInt64)Math.Abs(measuredOutput[counter] - expectedOutput[counter]);

                                if (t != 0)
                                {
                                    weight += 1;
                                }

                                solutionFitness += t;
                            }

                            // fewer number of total elements will rank better
                            UInt64 otherFitness = creature.CountTotalChildElements();

                            fitness[i] = (solutionFitness * Simulation.ResultDifferenceWeight) + (Simulation.CorrectCountResultWeight * weight) + otherFitness;
                        }
                        catch
                        {
                            success = false;
                            creature.Fitness = int.MaxValue;
                            creature.Health = Creature.LifeState.Dead;
                        }
                    }
                    else
                    {
                        success = false;
                        creature.Fitness = int.MaxValue;
                        creature.Health = Creature.LifeState.Dead;
                    }

                    // next creature
                    if (success == false)
                    {
                        break;
                    }
                }

                if (success)
                {
                    BigInteger sum = 0;
                    BigInteger average = 1;

                    for (int i = 0; i < Simulation.FitnessAverageCount; i++)
                    {
                        sum += fitness[i];
                    }

                    average = sum / Simulation.FitnessAverageCount;

                    // should not cause an overflow as each value should be less than UInt64.maxvalue
                    creature.Fitness = (UInt64)average;
                }
            });
        }

        public void MoveNext()
        {
            // sort
            _population = _population.OrderBy(x => x.Fitness).ToList();

            // update lifespan
            Parallel.ForEach(_population, creature =>
            {
                creature.Lifespan++;
            });

            // calculate best fitness
            BigInteger sum = 0;
            BigInteger average = 1;

            for (int i = 0; i < Simulation.FitnessAverageCount; i++)
            {
                sum += _population[i].Fitness;
            }

            average = sum / Simulation.FitnessAverageCount;

            Console.WriteLine(String.Format("{1,10}: top (avg): {0:n0}", average, Generation.ToString("n0")));
            
            // save best fitness to track convergence
            if (_population[0].Fitness <= (UInt64)_bestMaxFitness)
            {
                _convergenceCount++;

                if (_convergenceCount > BestBeforeConvergence)
                {
                    _hasConverged = true;

                    _population = _population.OrderBy(x => x.Fitness).ToList();

                    return;
                }
            }
            else
            {
                _convergenceCount = 0;
            }

            // save record breaking creatures to the hall of fame
            if (_population[0].Fitness < (UInt64)(_bestFitness / 2) && _hallOfFame.Count < _hallOfFameMaxLength)
            {
                _hallOfFame.Add(_population[0].Clone());

                _bestFitness = _population[0].Fitness;

                _lastBestFitnessGeneration = (UInt64)Generation;
            }

            int replaceCount = 0;

            int startIndex = GenerationalDeathPopulationStartIndex;
            bool death = false;

            // is this an apocalypse generation?
            if ((UInt64)Generation - _lastBestFitnessGeneration > ApocalypseGeneration)
            {
                startIndex = 0;
                _bestFitness = UInt64.MaxValue;

                death = true;
            }
            // is this a hero generation?
            else if ((UInt64)Generation - _lastBestFitnessGeneration > RecallHeroesGeneration)
            {
                // add some heros
                foreach (var h in _hallOfFame)
                {
                    _population.Insert(0, h.Clone());
                }
            }
            
            Parallel.For(startIndex, _population.Count, i =>
            {
                _population[i].Health = Creature.LifeState.Dead;
            });

            // if this is an apocalypse generation, don't ressurect creatures.
            if (!death)
            {
                _population.GroupBy(x => x.Fitness).Take(OtherFitnessRanksCount).ToList().ForEach(x =>
                {
                    int notPrimeCount = OtherFitnessCount;

                    x.ToList().ForEach(y =>
                    {
                        if (notPrimeCount-- > 0)
                        {
                            y.Health = Creature.LifeState.Alive;
                        }
                    });
                });
            }

            // add some spontaneous additions
            for (int i = 0; i < SpontaneousAdditionsPerGeneration; i++)
            {
                _population.Add(Creature.NewSeedCreature());
            }
            
            _population.RemoveAll(x => x.Health == Creature.LifeState.Dead);

            // now increase population back from the remaining population
            replaceCount = PopulationSize - _population.Count;
            
            List<Creature> nextGeneration = new List<Creature>();
            int currentPopulationSize = _population.Count;
            int bestFitnessCutoff = (int)(replaceCount / 10);

            // tried replacing this with a parallel for, but that was much much much slower.
            while (replaceCount-- > 0)
            {
                bool clone = Simulation.Rand.Next(0, 100) < ChanceToClone;

                // more fit creatures have more offspring
                int parent1 = replaceCount > bestFitnessCutoff ? Simulation.Rand.Next(0, currentPopulationSize) : 0;
                int parent2 = replaceCount > 5 ? Simulation.Rand.Next(0, currentPopulationSize) : 0;

                if (clone)
                {
                    nextGeneration.Add(Creature.NewCreature(_population[parent1]));
                }
                else
                {
                    nextGeneration.Add(Creature.NewCreature(_population[parent1], _population[parent2]));
                }
            }
             
             nextGeneration.ForEach(x => _population.Add(x));

            _generation++;
        }
    }

    public class Creature
    {
        private UInt64 _fitnessLevel = int.MaxValue;
        private LifeState _health = LifeState.Alive;

        /// <summary>
        /// Number of generations survived.
        /// </summary>
        private int _lifespan = 0;

        List<Chromosome> _chromosomes = new List<Chromosome>();

        public LifeState Health
        {
            get
            {
                return _health;
            }

            set
            {
                _health = value;
            }
        }

        public List<Chromosome> Chromosomes
        {
            get
            {
                return _chromosomes;
            }

            set
            {
                _chromosomes = value;
            }
        }

        public UInt64 Fitness
        {
            get
            {
                return _fitnessLevel;
            }

            set
            {
                _fitnessLevel = value;
            }
        }

        public int Lifespan
        {
            get
            {
                return _lifespan;
            }

            set
            {
                _lifespan = value;
            }
        }

        [Flags]
        public enum LifeState
        {
            Unknown,

            Alive,

            Zombie,

            Dead
        }

        public UInt64 CountTotalChildElements()
        {
            UInt64 sum = 0;

            foreach (var c in _chromosomes)
            {
                sum += c.CountTotalChildElements();
            }

            return sum;
        }

        public static Creature NewCreature(Creature parent)
        {
            Random r = Simulation.Rand;
            Creature c = parent.Clone();

            c._fitnessLevel = int.MaxValue;

            int maxMutations = r.Next(0, Simulation.MutationsBetweenGenerations);

            // mutate!
            for (int i = 0; i < maxMutations; i++)
            {
                int child = r.Next(0, c._chromosomes.Count);

                c.Chromosomes[child].Mutate();
            }

            c.Chromosomes[0].SetPreserve(Statement.IOFlags.RightSideInput);
            c.Chromosomes[1].SetPreserve(Statement.IOFlags.None);
            c.Chromosomes[2].SetPreserve(Statement.IOFlags.LeftSideOutput);

            return c;
        }

        public static Creature NewCreature(Creature parent1, Creature parent2)
        {
            Random r = Simulation.Rand;
            Creature c = new Creature();

            c.Health = LifeState.Alive;
            c._fitnessLevel = int.MaxValue;

            for (int i = 0; i < Simulation.MaxChromosomesPerCreature; i++)
            {
                if (r.Next(0, 2) == 0)
                {
                    c.Chromosomes.Add(parent1.Chromosomes[i].Clone());
                }
                else
                {
                    c.Chromosomes.Add(parent2.Chromosomes[i].Clone());
                }
            }

            int maxMutations = r.Next(0, Simulation.MutationsBetweenGenerations);

            // mutate!
            for (int i = 0; i < maxMutations; i++)
            {
                int child = r.Next(0, c._chromosomes.Count);

                c.Chromosomes[child].Mutate();
            }

            c.Chromosomes[0].SetPreserve(Statement.IOFlags.RightSideInput);
            c.Chromosomes[1].SetPreserve(Statement.IOFlags.None);
            c.Chromosomes[2].SetPreserve(Statement.IOFlags.LeftSideOutput);

            return c;
        }

        public static Creature NewSeedCreature()
        {
            Chromosome c = null;
            Creature creature = new Creature();

            // input
            c = new Chromosome();
            c.Blocks.Add(Block.NewRand(Statement.IOFlags.RightSideInput));
            //c.Blocks.Add(Block.NewRand(Statement.Unmutateable.RightSideEternal));
            creature.Chromosomes.Add(c);

            // normal
            c = new Chromosome();
            c.Blocks.Add(Block.NewRand(Statement.IOFlags.None));
            //c.Blocks.Add(Block.NewRand(Statement.Unmutateable.None));
            creature.Chromosomes.Add(c);

            // output
            c = new Chromosome();
            c.Blocks.Add(Block.NewRand(Statement.IOFlags.LeftSideOutput));
            //c.Blocks.Add(Block.NewRand(Statement.Unmutateable.LeftSideEternal));
            creature.Chromosomes.Add(c);

            return creature;
        }

        public bool Eval(int[] registers, int[] input, int[] output)
        {
            bool success = false;

            foreach (var c in _chromosomes)
            {
                success = c.Eval(registers, input, output);

                if (success == false)
                {
                    return false;
                }
            }

            return success;
        }

        public Creature Clone()
        {
            Creature c = new Creature();

            c._fitnessLevel = _fitnessLevel;
            c.Health = Health;

            _chromosomes.ForEach(x => c._chromosomes.Add(x.Clone()));

            return c;
        }

        public override string ToString()
        {
            string s = String.Empty;

            foreach (var c in _chromosomes)
            {
                s += c.ToString() + "\r\n";
            }

            return s;
        }
    }

    public class Chromosome
    {
        List<Block> _blocks = new List<Block>();

        public List<Block> Blocks
        {
            get
            {
                return _blocks;
            }

            set
            {
                _blocks = value;
            }
        }

        public Chromosome Clone()
        {
            Chromosome c = new Chromosome();

            _blocks.ForEach(x => c._blocks.Add(x.Clone()));

            return c;
        }

        public UInt64 CountTotalChildElements()
        {
            UInt64 sum = 0;

            foreach (var b in _blocks)
            {
                sum += b.CountTotalChildElements();
            }

            return sum;
        }

        public void SetPreserve(Statement.IOFlags state)
        {
            _blocks.ForEach(x => x.SetPreserve(state));
        }

        public void Mutate()
        {
            Random r = Simulation.Rand;

            // duplicate statement
            // create new statement
            if (_blocks.Count < Simulation.MaxStatementsPerBlock && r.Next(0, 100) < Simulation.ChanceToCreate)
            {
                // duplicate
                if (r.Next(0, 2) == 0)
                {
                    int child = r.Next(0, _blocks.Count);

                    Block b = _blocks[child].Clone();

                    // chance to mutate new child
                    if (r.Next(0, 2) == 0)
                    {
                        b.Mutate();
                    }

                    _blocks.Insert(child, b);
                }
                // create
                else
                {
                    _blocks.Add(Block.NewRand(Statement.IOFlags.None));
                }

                return;
            }

            if (_blocks.Count > 1 && r.Next(0, 100) < Simulation.ChanceToDelete)
            {
                // delete statement
                int removeOne = r.Next(0, Simulation.MaxBlocksPerChromosome - 1);
                if (_blocks.Count > removeOne)
                {
                    int child = r.Next(0, _blocks.Count);

                    _blocks.RemoveAt(child);

                    return;
                }
            }

            // otherwise, swap two items
            int indexA = r.Next(0, _blocks.Count);
            int indexB = r.Next(0, _blocks.Count);

            Block temp = _blocks[indexA];
            _blocks[indexA] = _blocks[indexB];
            _blocks[indexB] = temp;

            if (r.Next(0, 2) == 0)
            {
                _blocks[indexA].Mutate();
            }

            if (r.Next(0, 2) == 0)
            {
                _blocks[indexB].Mutate();
            }
        }

        public bool Eval(int[] registers, int[] input, int[] output)
        {
            bool success = false;

            foreach (var b in _blocks)
            {
                success = b.Eval(registers, input, output);

                if (success == false)
                {
                    return false;
                }
            }

            return success;
        }

        public override string ToString()
        {
            string s = String.Empty;

            foreach (var b in _blocks)
            {
                s += b.ToString() + "\r\n";
            }

            return s;
        }
    }

    public class Block
    {
        List<Statement> _statements = new List<Statement>();

        public List<Statement> Statements
        {
            get
            {
                return _statements;
            }

            set
            {
                _statements = value;
            }
        }

        public Block Clone()
        {
            Block b = new Block();

            _statements.ForEach(x => b._statements.Add(x.Clone()));

            return b;
        }

        public static Block NewRand(Statement.IOFlags mutateFlags)
        {
            Random r = Simulation.Rand;
            Block b = new Block();

            b._statements.Add(Statement.NewRand(mutateFlags));

            if (r.Next(0, 2) == 0)
            {
                b._statements.Add(Statement.NewRand(mutateFlags));
            }

            if (r.Next(0, 3) == 0)
            {
                b._statements.Add(Statement.NewRand(mutateFlags));
            }

            if (r.Next(0, 4) == 0)
            {
                b._statements.Add(Statement.NewRand(mutateFlags));
            }

            return b;
        }

        public void SetPreserve(Statement.IOFlags state)
        {
            _statements.ForEach(x => x.SetPreserve(state));
        }

        public UInt64 CountTotalChildElements()
        {
            return (UInt64)_statements.Count;
        }

        public void Mutate()
        {
            Random r = Simulation.Rand;

            // duplicate statement
            // create new statement
            if (_statements.Count < Simulation.MaxStatementsPerBlock && r.Next(0, 100) < Simulation.ChanceToCreate)
            {
                // duplicate
                if (r.Next(0, 2) == 0)
                {
                    int child = r.Next(0, _statements.Count);

                    Statement s = _statements[child].Clone();

                    // chance to mutate new child
                    if (r.Next(0, 2) == 0)
                    {
                        s.Mutate();
                    }

                    _statements.Insert(child, s);
                }
                // create
                else
                {
                    _statements.Add(Statement.NewRand(Statement.IOFlags.None));
                }

                return;
            }

            if (_statements.Count > 1 && r.Next(0, 100) < Simulation.ChanceToDelete)
            {
                // delete statement
                int removeOne = r.Next(0, Simulation.MaxStatementsPerBlock - 1);
                if (_statements.Count > removeOne)
                {
                    int child = r.Next(0, _statements.Count);

                    _statements.RemoveAt(child);

                    return;
                }
            }

            // otherwise, swap two items
            int indexA = r.Next(0, _statements.Count);
            int indexB = r.Next(0, _statements.Count);

            Statement temp = _statements[indexA];
            _statements[indexA] = _statements[indexB];
            _statements[indexB] = temp;

            if (r.Next(0, 2) == 0)
            {
                _statements[indexA].Mutate();
            }

            if (r.Next(0, 2) == 0)
            {
                _statements[indexB].Mutate();
            }
        }

        public bool Eval(int[] registers, int[] input, int[] output)
        {
            bool success = false;

            foreach (var s in _statements)
            {
                success = s.Eval(registers, input, output);

                if (success == false)
                {
                    return false;
                }
            }

            return success;
        }

        public override string ToString()
        {
            string s = String.Empty;

            foreach (var statement in _statements)
            {
                s += statement.ToString() + "\r\n";
            }

            return s;
        }
    }

    public class Statement
    {
        private StatementClass _family;
        private IOFlags _preserve;

        private int _lhsIndex;

        private int _rhsIndex1;

        /// <summary>
        /// Second rhs var or constant
        /// </summary>
        private int _rhs2;
        private Op _op;

        public int LhsIndex
        {
            get
            {
                return _lhsIndex;
            }

            set
            {
                _lhsIndex = value;
            }
        }

        public int RhsIndex1
        {
            get
            {
                switch (_family)
                {
                    case StatementClass.VarConst: /* fall through */
                    case StatementClass.VarVar:
                        return _rhsIndex1;
                    case StatementClass.Const:
                    case StatementClass.Unknown:
                    default:
                        return -1;
                }
            }

            set
            {
                switch (_family)
                {
                    case StatementClass.VarConst:
                        _rhsIndex1 = value;
                        break;
                    case StatementClass.VarVar: /* can't fall through */
                        _rhsIndex1 = value;
                        break;
                    case StatementClass.Const:
                    case StatementClass.Unknown:
                    default:
                        return;
                }
            }
        }

        public int RhsIndex2
        {
            get
            {
                switch (_family)
                {
                    case StatementClass.VarVar:
                        return _rhs2;
                    case StatementClass.VarConst: /* fall through */
                    case StatementClass.Const:
                    case StatementClass.Unknown:
                    default:
                        return -1;
                }
            }

            set
            {
                switch (_family)
                {
                    case StatementClass.VarVar:
                        _rhs2 = value;
                        break;
                    case StatementClass.VarConst: /* fall through */
                    case StatementClass.Const:
                    case StatementClass.Unknown:
                    default:
                        return;
                }
            }
        }

        public int Constant
        {
            get
            {
                switch (_family)
                {
                    case StatementClass.VarConst: /* fall through */
                    case StatementClass.Const:
                        return _rhs2;
                    case StatementClass.VarVar:
                    case StatementClass.Unknown:
                    default:
                        return -1;
                }
            }

            set
            {
                switch (_family)
                {
                    case StatementClass.VarConst: /* fall through */
                    case StatementClass.Const:
                        _rhs2 = value;
                        break;
                    case StatementClass.VarVar:
                    case StatementClass.Unknown:
                    default:
                        return;
                }
            }
        }

        public Op Operator
        {
            get
            {
                return _op;
            }

            set
            {
                _op = value;
            }
        }

        public enum StatementClass
        {
            Unknown = 0,

            /// <summary>
            /// var = var op var
            /// </summary>
            VarVar = 1,

            /// <summary>
            /// var = var op const
            /// </summary>
            VarConst = 2,

            /// <summary>
            /// var = const
            /// </summary>
            Const = 3
        }

        [Flags]
        public enum IOFlags
        {
            None = 0,

            /// <summary>
            /// Left hand side of statement is output.
            /// </summary>
            LeftSideOutput = 1,

            /// <summary>
            /// Right hand side of statement is input.
            /// </summary>
            RightSideInput = 2
        }

        public static Statement NewRand(IOFlags mutateFlags)
        {
            Statement s = new Statement();

            s._preserve = mutateFlags;

            s._family = (StatementClass)Simulation.Rand.Next(1, 4);

            s._lhsIndex = Simulation.Rand.Next(0, 3);
            s._rhsIndex1 = Simulation.Rand.Next(0, 3);

            s._op = (Op)Simulation.Rand.Next(1, 11);

            if (mutateFlags != IOFlags.RightSideInput && (s._family == StatementClass.VarConst || s._family == StatementClass.Const))
            {
                s._rhs2 = NumberChange(true, 0);
            }
            else
            {
                s._rhs2 = Simulation.Rand.Next(0, 3);
            }

            return s;
        }

        private static int NumberChange(bool init, int original)
        {
            int rand = Simulation.Rand.Next(0, 100);

            // increment 50% of the time for small changes
            bool increment = Simulation.Rand.Next(0, 2) == 0;

            int newNumber = 0;

            if (init)
            {
                // usually a small number
                if (rand < 80)
                {
                    newNumber = Simulation.Rand.Next(0, 10);
                }
                else
                {
                    newNumber = Simulation.Rand.Next(0, int.MaxValue);
                }
            }
            else
            {
                newNumber = original;

                if (rand < 60)
                {
                    newNumber += increment ? 1 : -1;
                }
                else if (rand < 80)
                {
                    newNumber += increment ? (int)(newNumber * 0.1) : -(int)(newNumber * 0.1);
                }
                else
                {
                    newNumber = Simulation.Rand.Next(0, int.MaxValue);
                }
            }

            return newNumber;
        }

        public void Mutate()
        {
            int state;

            // increment or decrement
            int increment = Simulation.Rand.Next(0, 2);

            switch (_family)
            {
                // var = var op var
                case StatementClass.VarVar:
                    state = Simulation.Rand.Next(0, 4);
                    switch (state)
                    {
                        case 0:
                            _lhsIndex += increment == 1 ? 1 : -1;
                            break;
                        case 1:
                            _op = increment == 1 ? _op.Next() : _op.Prev();
                            break;
                        case 2:
                            _rhsIndex1 += increment == 1 ? 1 : -1;
                            break;
                        case 3:
                            _rhs2 += increment == 1 ? 1 : -1;
                            break;
                        default:
                            break;
                    }
                    break;
                // var = var op const
                case StatementClass.VarConst:
                    state = Simulation.Rand.Next(0, 6);
                    switch (state)
                    {
                        case 0:
                            _lhsIndex += increment == 1 ? 1 : -1;
                            break;
                        case 1:
                            _op = increment == 1 ? _op.Next() : _op.Prev();
                            break;
                        case 2:
                            _rhsIndex1 += increment == 1 ? 1 : -1;
                            break;
                        case 3: // fall through //
                        case 4:
                        case 5:
                            _rhs2 = NumberChange(false, _rhs2);
                            break;
                        default:
                            break;
                    }
                    break;
                // var = const
                case StatementClass.Const:
                    state = Simulation.Rand.Next(0, 5);
                    switch (state)
                    {
                        case 0:
                            _lhsIndex += increment == 1 ? 1 : -1;
                            break;
                        case 1: // fall through //
                        case 2: 
                        case 3:
                        case 4:
                            _rhs2 = NumberChange(false, _rhs2);
                            break;
                        default:
                            break;
                    }
                    break;
                case StatementClass.Unknown:
                default:
                    break;
            }
        }

        public void SetPreserve(Statement.IOFlags state)
        {
            _preserve = state;
        }

        public Statement Clone()
        {
            Statement s = new Statement();

            s._family = _family;
            s._lhsIndex = _lhsIndex;
            s._op = _op;
            s._preserve = _preserve;
            s._rhs2 = _rhs2;
            s._rhsIndex1 = _rhsIndex1;

            return s;
        }

        public override string ToString()
        {
            string s = String.Empty;

            string lhs = String.Empty;

            if (_preserve == IOFlags.LeftSideOutput)
            {
                lhs = String.Format("out{0} = ", _lhsIndex);
            }
            else
            {
                lhs = String.Format("r{0} = ", _lhsIndex);
            }
            
            string rhs = String.Empty;

            if (_preserve == IOFlags.RightSideInput)
            {
                switch (_family)
                {
                    case StatementClass.VarVar:
                        rhs = String.Format("r{0} {1} in{2}", _rhsIndex1, _op.OpToString(), _rhs2);
                        break;
                    case StatementClass.VarConst:
                        rhs = String.Format("r{0} {1} {2}", _rhsIndex1, _op.OpToString(), _rhs2);
                        break;
                    case StatementClass.Const:
                        rhs = String.Format("{0}", _rhs2);
                        break;
                    case StatementClass.Unknown:
                    default:
                        rhs = "unknown";
                        break;
                }
            }
            else
            {
                switch (_family)
                {
                    case StatementClass.VarVar:
                        rhs = String.Format("r{0} {1} r{2}", _rhsIndex1, _op.OpToString(), _rhs2);
                        break;
                    case StatementClass.VarConst:
                        rhs = String.Format("r{0} {1} {2}", _rhsIndex1, _op.OpToString(), _rhs2);
                        break;
                    case StatementClass.Const:
                        rhs = String.Format("{0}", _rhs2);
                        break;
                    case StatementClass.Unknown:
                    default:
                        rhs = "unknown";
                        break;
                }
            }

            s = lhs + rhs;

            return s;
        }

        public bool Eval(int[] registers, int[] consumeSource, int[] produceDest)
        {
            bool success = false;

            int registerLen = registers.Length - 1;
            int consumeLen = consumeSource.Length - 1;
            int produceLen = produceDest.Length - 1;

            if (_op == Op.Unknown)
            {
                return false;
            }

            switch (_family)
            {
                // var = var op var
                case StatementClass.VarVar:
                    
                    // Left side doesn't change, therefore writing to final output
                    if (_preserve == IOFlags.LeftSideOutput)
                    {
                        if (_lhsIndex < 0 || _lhsIndex > produceLen ||
                            _rhsIndex1 < 0 || _rhsIndex1 > registerLen ||
                            _rhs2 < 0 || _rhs2 > registerLen)
                        {
                            return false;
                        }

                        if ((_op == Op.Divide || _op == Op.Modulus) && (registers[_rhs2] == 0 || (registers[_rhsIndex1] == -2147483648 && registers[_rhs2] == -1)))
                        {
                            return false;
                        }

                        produceDest[_lhsIndex] = Utilities.OpCalc(_op, registers[_rhsIndex1], registers[_rhs2]);
                    }
                    // right side doesn't change, therefore reading from initial input
                    else if (_preserve == IOFlags.RightSideInput)
                    {
                        if (_lhsIndex < 0 || _lhsIndex > registerLen ||
                            _rhsIndex1 < 0 || _rhsIndex1 > registerLen ||
                            _rhs2 < 0 || _rhs2 > consumeLen)
                        {
                            return false;
                        }

                        if ((_op == Op.Divide || _op == Op.Modulus) && (consumeSource[_rhs2] == 0 || (registers[_rhsIndex1] == -2147483648 && consumeSource[_rhs2] == -1)))
                        {
                            return false;
                        }

                        registers[_lhsIndex] = Utilities.OpCalc(_op, registers[_rhsIndex1], consumeSource[_rhs2]);
                    }
                    else
                    {
                        if (_lhsIndex < 0 || _lhsIndex > registerLen ||
                            _rhsIndex1 < 0 || _rhsIndex1 > registerLen ||
                            _rhs2 < 0 || _rhs2 > registerLen)
                        {
                            return false;
                        }

                        if ((_op == Op.Divide || _op == Op.Modulus) && (registers[_rhs2] == 0 || (registers[_rhsIndex1] == -2147483648 && registers[_rhs2] == -1)))
                        {
                            return false;
                        }

                        registers[_lhsIndex] = Utilities.OpCalc(_op, registers[_rhsIndex1], registers[_rhs2]);
                    }

                    success = true;
                    break;

                // var = var op const
                case StatementClass.VarConst:

                    // Left side doesn't change, therefore writing to final output
                    if (_preserve == IOFlags.LeftSideOutput)
                    {
                        if (_lhsIndex < 0 || _lhsIndex > produceLen ||
                            _rhsIndex1 < 0 || _rhsIndex1 > registerLen)
                        {
                            return false;
                        }
  
                        if ((_op == Op.Divide || _op == Op.Modulus) && (_rhs2 == 0 || (registers[_rhsIndex1] == -2147483648 && _rhs2 == -1)))
                        {
                            return false;
                        }

                        produceDest[_lhsIndex] = Utilities.OpCalc(_op, registers[_rhsIndex1], _rhs2);
                    }
                    // right side doesn't change, therefore reading from initial input, but there isn't an input register to read from
                    else 
                    {
                        if (_lhsIndex < 0 || _lhsIndex > registerLen ||
                            _rhsIndex1 < 0 || _rhsIndex1 > registerLen)
                        {
                            return false;
                        }

                        if ((_op == Op.Divide || _op == Op.Modulus) && (_rhs2 == 0 || (registers[_rhsIndex1] == -2147483648 && _rhs2 == -1)))
                        {
                            return false;
                        }

                        registers[_lhsIndex] = Utilities.OpCalc(_op, registers[_rhsIndex1], _rhs2);
                    }

                    
                    success = true;
                    break;

                // var = const
                case StatementClass.Const:

                    // Left side doesn't change, therefore writing to final output, but this is illegal.
                    // Or right side doesn't change, therefore reading from initial input, but that's illegal too.
                    if (_preserve == (IOFlags.LeftSideOutput | IOFlags.RightSideInput))
                    {
                        return false;
                    }

                    if (_lhsIndex < 0 || _lhsIndex > registerLen)
                    {
                        return false;
                    }

                    registers[_lhsIndex] = _rhs2;
                    success = true;
                    break;

                case StatementClass.Unknown:
                default:
                    return false;
            }

            return success;
        }
    }
}
