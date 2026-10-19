#include "driver.hpp"
#include "common.hpp"
#include "caster.hpp"

using namespace std::string_literals;
using namespace std;

using Caster = caster::Color<caster::StepwiseColorDefaultAttributes<unsigned short, unsigned short> >;
using AnnotatedBinaryTree = common::AnnotatedBinaryTree;

class NJ {
    size_t n, nextID;
    vector<vector<double>> dist;
    vector<double> R;
    vector<bool> active;
    vector<size_t> activeList;
    vector<AnnotatedBinaryTree> clusters;

private:
    pair<size_t,size_t> minQPair() {
        double minQ = 1e100;
        pair<size_t,size_t> minQpair{0, 0};
        double m = static_cast<double>(activeList.size());
        for (size_t active_i = 0; active_i < activeList.size(); ++active_i) {
            size_t i = activeList[active_i];
            for (size_t active_j = 0; active_j < active_i; ++active_j) {
                size_t j = activeList[active_j];
                double Q = (m - 2) * dist[i][j] - R[i] - R[j];
                if (Q < minQ) {
                    minQ = Q;
                    minQpair = pair<size_t,size_t>(i, j);
                }
            }
        }
        return minQpair;
    }

    void merge(size_t i, size_t j) {
        double m = static_cast<double>(activeList.size());
        double dij = dist[i][j];
        double Li = dij / 2 + (R[i] - R[j]) / (2 * (m - 2));
        double Lj = dij / 2 - (R[i] - R[j]) / (2 * (m - 2));

        clusters[i].root()->set("LENGTH", Li);
        clusters[j].root()->set("LENGTH", Lj);

        clusters[nextID] = std::move(clusters[i]);
        clusters[nextID].merge(clusters[j]);

        active[i] = false;
        active[j] = false;

        vector<size_t> newactiveList;
        newactiveList.reserve(activeList.size());
        dist[nextID][nextID] = 0.0;

        for (size_t k : activeList) {
            if (k == i || k == j) continue;
            newactiveList.push_back(k);
            dist[k][nextID] = dist[nextID][k] = (dist[i][k] + dist[j][k] - dij) / 2;
            R[k] = R[k] - dist[i][k] - dist[j][k] + dist[k][nextID];
        }

        newactiveList.push_back(nextID);
        R[nextID] = (R[i] + R[j] - m * dij) / 2;
        active[nextID] = true;
        activeList = std::move(newactiveList);
        ++nextID;
    }

public:
    NJ(vector<vector<double>> const& d):
        n(d.size()),
        nextID(d.size()),
        dist((d.empty()) ? 0 : d.size() * 2 - 1,
             vector<double>((d.empty()) ? 0 : d.size() * 2 - 1, 0.0)),
        R((d.empty()) ? 0 : d.size() * 2 - 1),
        active((d.empty()) ? 0 : 2 * d.size() - 1, false),
        clusters((d.empty()) ? 0 : d.size() * 2 - 1) {
        if (n == 0) return;
        size_t N = 2 * n - 1;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                dist[i][j] = d[i][j];
                if (i != j) R[i] += d[i][j];
            }
            dist[i][i] = 0.0;
        }
        for (size_t i = n; i < N; ++i) {
            dist[i][i] = 0.0;
        }
        for (size_t i = 0; i < n; ++i) {
            active[i] = true;
            activeList.push_back(i);
            clusters[i].emplaceRoot();
            clusters[i].root()->set(AnnotatedBinaryTree::LEAF_ID, i);
        }
    }

    AnnotatedBinaryTree buildTree() {
        if (n == 0) return AnnotatedBinaryTree();

        while (activeList.size() > 3) {
            auto [i, j] = minQPair();
            merge(i, j);
        }

        if (activeList.size() == 1) {
            return std::move(clusters[activeList[0]]);
        }

        if (activeList.size() == 2) {
            size_t i = activeList[0], j = activeList[1];
            double dij = dist[i][j];
            clusters[i].root()->set("LENGTH", dij / 2);
            clusters[j].root()->set("LENGTH", dij / 2);
            size_t rootID = nextID++;
            clusters[rootID] = std::move(clusters[i]);
            clusters[rootID].merge(clusters[j]);
            return std::move(clusters[rootID]);
        }

        // activeList.size() == 3
        size_t i = activeList[0], j = activeList[1], k = activeList[2];
        double dij = dist[i][j], dik = dist[i][k], djk = dist[j][k];
        double Li = (dij + dik - djk) / 2;
        double Lj = (dij + djk - dik) / 2;
        double Lk = (dik + djk - dij) / 2;

        clusters[i].root()->set("LENGTH", Li);
        clusters[j].root()->set("LENGTH", Lj);
        size_t t = nextID++;
        clusters[t] = std::move(clusters[i]);
        clusters[t].merge(clusters[j]);

        clusters[t].root()->set("LENGTH", Lk / 2);
        clusters[k].root()->set("LENGTH", Lk / 2);
        size_t rootID = nextID++;
        clusters[rootID] = std::move(clusters[t]);
        clusters[rootID].merge(clusters[k]);

        return std::move(clusters[rootID]);
    }
};

void run(Caster::SharedConstData &data) {
	Caster caster(data);
	size_t n = common::taxonName2ID.nTaxa();
	size_t k = data.nElements();
	// scores[k][i][j]
	vector<vector<vector<double> > > scores(k, vector<vector<double> >(n, vector<double>(n, 0.0)));
	
	for (int i = 0; i < n; ++i) {
		for (int x = 0; x < k; ++x) {
			caster.elementSetTaxonColor(x, i, 0);
			caster.elementSetTaxonColor(x, i, 3);
		}
	}

	for (int i = 0; i < n; ++i) {
		for (int x = 0; x < k; ++x) {
			caster.elementClearTaxonColor(x, i, 0);
			caster.elementClearTaxonColor(x, i, 3);
			caster.elementSetTaxonColor(x, i, 1);
		}

		for (int j = 0; j < i; ++j) {
			for (int x = 0; x < k; ++x) {
				caster.elementClearTaxonColor(x, j, 0);
				caster.elementClearTaxonColor(x, j, 3);
				caster.elementSetTaxonColor(x, j, 2);
				scores[x][i][j] = scores[x][j][i] = caster.elementQuadripartitionScores(x)[0];
				array<double, 3> s = caster.elementQuadripartitionScores(x);
				caster.elementClearTaxonColor(x, j, 2);
				caster.elementSetTaxonColor(x, j, 0);
				caster.elementSetTaxonColor(x, j, 3);
			}
		}

		for (int x = 0; x < k; ++x) {
			caster.elementClearTaxonColor(x, i, 1);
			caster.elementSetTaxonColor(x, i, 0);
			caster.elementSetTaxonColor(x, i, 3);
		}
	}

	for (int x = 0; x < k; ++x) {
		if (x > 0) cout << endl;
		NJ nj(scores[x]);
		AnnotatedBinaryTree tree = nj.buildTree();
		tree.displaySimpleNewick<int, double>(
        std::cout, "__NO_INTERNAL_LABEL__", "LENGTH");
    	std::cout << "\n";
	}
}

int main(int argc, char* argv[]) {
	using std::string;
	using std::size_t;
	using Clock = std::chrono::high_resolution_clock;
	using Driver = caster::Driver<true>;

	std::pair<string, string> programNames = Driver::programNames();
	ARG.set("SHORT_NAME", programNames.first);
	ARG.set("FULL_NAME", programNames.second);
	
	ARG.addArgument('h', "help", "flag", "Display help message", 6, true);
	ARG.addArgument('i', "input", "string", "Input file path", 5);
	//ARG.addArgument('o', "output", "string", "Output file path, print to stdout if not provided", 5, true);
	ARG.addArgument('a', "mapping", "string", "Mapping file path, a list of gene/speicesman name to taxon name maps, each line contains one gene/speicesman name followed by one taxon name separated by a space or tab", 3, true, false);
	ARG.addArgument('\0', "verbose", "integer", "Verbose level", 0, true, true, "4");
	ARG.addArgument('\0', "no-log", "flag", "Don't generate log file", 1, true);
	ARG.addArgument('\0', "log", "string", "Log file path", 0, true, true, "log.txt");
	ARG.addArgument('\0', "log-verbose", "integer", "Verbose level in log file", 0, true, true, "5");

	ARG.addArgument('\0', "chunk", "integer", "The maximum number of sites in each local aligment block for parameter estimation", 0, true, true, "10000");

	ARG.parse(argc, argv);
	common::LogInfo::setVerbose(ARG.has("no-log") ? nullptr : new std::ofstream(ARG.get<string>("log")), ARG.get<size_t>("log-verbose"), ARG.get<size_t>("verbose"));
	ARG.print();

	if (ARG.has("root")) common::taxonName2ID[ARG.get<string>("root")];
	ARG.log() << "Parsing input file(s)..." << endl;
	Caster::SharedConstData data = caster::DriverHelper::read<Caster::SharedConstData>();
	ARG.log() << "Start running ..." << endl;
	run(data);

	return 0;
}