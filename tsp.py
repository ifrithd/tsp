import math, random
import matplotlib.pyplot as plt
from copy import deepcopy

DEBUG = False
POPULATION = 150
GENERATION = 1000
MUTATE_RATE = (0.1, 0.2)
CROSS_RATE = (0.9, 0.9)
INF = 0x3f3f3f


class Gene:
    cities = None

    def __init__(self, dna):
        self.dna = dna
        self.len = len(dna)
        self._distance = None
        self._fitness = None

    def mutate(self):
        r1 = random.randrange(self.len - 2)
        r2 = random.randrange(r1 + 1, self.len - 1)
        r3 = random.randrange(r2 + 1, self.len)
        self.dna = self.dna[:r1] + self.dna[r2+1:r3+1] + self.dna[r1:r2+1] + self.dna[r3+1:]
        self._distance = self._fitness = None

    @property
    def fitness(self):
        if self._fitness:
            return self._fitness
        self._fitness = 1 / self.distance
        return self._fitness

    @property
    def distance(self):
        if self._distance:
            return self._distance
        self._distance = self.calc_city_distance(self.dna[0], self.dna[-1])
        for i in range(len(self.dna) - 1):
            self._distance += self.calc_city_distance(self.dna[i], self.dna[i+1])
        return self._distance

    @staticmethod
    def calc_city_distance(a, b):
        a, b = a % len(Gene.cities), b % len(Gene.cities)
        city_a, city_b  = Gene.cities[a], Gene.cities[b]
        return math.sqrt((city_a.x - city_b.x) ** 2 + (city_a.y - city_b.y) ** 2)

    @staticmethod
    def cross(gene1, gene2):
        r1 = random.randrange(gene1.len - 1)
        r2 = random.randrange(r1 + 1, gene1.len)
        parent1, parent2 = gene1.dna, gene2.dna
        child1 = parent1[:r1] + parent2[r1:r2] + parent1[r2:]
        child2 = parent2[:r1] + parent1[r1:r2] + parent2[r2:]

        dup1 = []
        dup2 = []
        for j in range(gene1.len):
            try:
                dup1.append(child1.index(child1[j], j + 1))
                if r1 <= dup1[-1] < r2:
                    dup1[-1] = j
            except ValueError:
                pass
        for j in range(gene1.len):
            try:
                dup2.append(child2.index(child2[j], j + 1))
                if r1 <= dup2[-1] < r2:
                    dup2[-1] = j
            except ValueError:
                pass
        for j in range(len(dup1)):
            child1[dup1[j]], child2[dup2[j]] = child2[dup2[j]], child1[dup1[j]]
        return Gene(child1), Gene(child2)

    def __str__(self):
        ret = 'gene: {} distance:{}'.format(self.dna, self.distance)
        return ret


def init_genes(city_num):
    genes = []
    for i in range(POPULATION):
        gene = list(range(city_num))
        random.shuffle(gene)
        dis = Gene.calc_city_distance
        for u in range(len(gene)-1):
            for v in range(u+1, len(gene)):
                if dis(u-1,v) + dis(u,v+1) - dis(u-1,u) - dis(v,v+1):
                    gene = gene[:u] + gene[u:v+1][::-1] + gene[v+1:]
        genes.append(Gene(gene))
    return genes


def roulette_wheel(genes):
    # 轮盘赌法选择新一代
    distances = [gene.distance for gene in genes]
    minn = 0x3f3f3f
    for d in distances:
        minn = min(minn, d)
    for i in range(len(genes)):
        distances[i] = distances[i] - minn + 100
    fitness = [1/d for d in distances]
    s_fitness = sum(fitness)
    new_genes = []

    for t in range(POPULATION):
        p = random.random()
        for i in range(len(genes)):
            if fitness[i]/s_fitness > p:
                new_genes.append(deepcopy(genes[i]))
                break
            p -= fitness[i] / s_fitness
    return new_genes


def tournament(genes):
    # 锦标赛法选择新一代
    new_genes = []

    while len(new_genes) < POPULATION:
        r1 = random.randrange(len(genes))
        r2 = random.randrange(len(genes))
        #r3 = random.randrange(len(genes))
        #x = compete(compete(r1, r2), r3)
        x = r1 if genes[r1].fitness > genes[r2].fitness else r2
        new_genes.append(deepcopy(genes[x]))
    return new_genes


def elitism(genes):
    new_genes = sorted(genes, key=lambda gene:gene.fitness, reverse=True)
    new_genes = new_genes[:POPULATION]
    return new_genes


def genes_mutate(genes, pm):
    for gene in genes:
        p = random.random()
        if p > pm:
            continue
        gene.mutate()


def genes_cross(genes, pc):
    genes_num = len(genes) - 1 if len(genes)&1 else len(genes)
    for i in range(0, genes_num, 2):
        p = random.random()
        if p > pc:
            continue
        child1, child2 = Gene.cross(genes[i],genes[i+1])
        genes.append(child1)
        genes.append(child2)


def find_best_gene(genes):
    minn = INF
    ret = None
    for gene in genes:
        if minn > gene.distance:
            minn = gene.distance
            ret = gene
    return deepcopy(ret)


def trans_ans(start_city, cities_num, ans):
    # 指定路径（环）的起始城市
    try:
        st = ans.dna.index(start_city)
    except ValueError:
        return ans
    ret = []
    for i in range(cities_num):
        ret.append(ans.dna[(st+i) % cities_num])
    return ret


def tsp(cities):
    """计算经过所有城市的最短路径

    Args:
        cities: 参与计算的城市的列表

    Returns:
        返回一个Gene的实例，其为最佳个体，代表最短路径。
    """
    Gene.cities = cities
    cities_num = len(cities)
    genes = init_genes(cities_num)
    best_of_gen = []
    best_of_all = []
    ans_gen = 0
    ans = None
    #genes_select = roulette_wheel
    #genes_select = tournament
    genes_select = elitism

    t = 0
    while t<GENERATION:
        pc = CROSS_RATE[1] - t * (CROSS_RATE[1] - CROSS_RATE[0]) // GENERATION
        genes_cross(genes, pc)
        pm = MUTATE_RATE[1] - t * (MUTATE_RATE[1] - MUTATE_RATE[0]) // GENERATION
        genes_mutate(genes, pm)
        genes = genes_select(genes)

        tmp_ans = find_best_gene(genes)
        best_of_gen.append(tmp_ans.distance)
        try:
            best_of_all.append(best_of_all[-1])
        except IndexError:
            best_of_all.append(INF)
        if best_of_gen[-1] < best_of_all[-1]:
            best_of_all[-1] = best_of_gen[-1]
            ans_gen = t
            ans = tmp_ans
            print('update at Gen{}'.format(t), best_of_all[-1])

        t += 1

    print('shortest distance: {:.3f}'.format(best_of_all[-1]))
    print(ans.dna)
    print('found it at Gen{}'.format(ans_gen))

    if DEBUG:
        fig1 = plt.figure('figure')
        ax1 = fig1.add_subplot(2, 1, 1)
        ax1.plot(best_of_all)
        ax1.set_title('best ans')

        ax2 = fig1.add_subplot(2, 1, 2)
        ax2.set_title('best ans of generation')
        ax2.plot(best_of_gen)

        plt.show()
    return ans