#include "SingleLinkage.h"

using namespace std;

void quickprobs::SingleLinkage::build()
{
	int i, j;
	int next;
	int n_seq = numSeqs;

	vector<int> pi(n_seq);
	vector<float> lambda(n_seq);
	vector<float> sim_vector(n_seq);

	for (i = 0; i < n_seq; ++i) {
		pi[i] = i;
		lambda[i] = std::numeric_limits<float>::min() ;

		#pragma omp parallel for schedule(static)
		for(j = 0; j < i; ++j) {
			sim_vector[j] = distances[i][j];
		}

		for (j = 0; j < i; ++j) {
			next = pi[j];
			if(lambda[j] > sim_vector[j]) {
				sim_vector[next] = max(sim_vector[next], sim_vector[j]);
			} else {
				sim_vector[next] = max(lambda[j], sim_vector[next]);
				pi[j] = i;
				lambda[j] = sim_vector[j];
			}
		}

		for (j = 0; j < i; ++j) {
			next = pi[j];
			if(lambda[next] >= lambda[j]) {
				pi[j] = i;
			}
		}
	}

	vector<int> elements(n_seq-1);
	for (i = 0; i < n_seq-1; ++i) {
		elements[i] = i;
	}

	sort(elements.begin(), elements.end(), [&](int x, int y){
		return lambda[x] > lambda[y];
	});

	vector<int> index(n_seq);
	for (i = 0; i < n_seq; ++i) {
		index[i] = i;
	}
	
	for (i = 0; i < n_seq-1; ++i) {
		j = elements[i];
		next = pi[j];
		
		Node &n1 = nodes[next];
		Node &n2 = nodes[j];

		float branchLength = 0.5f;
		this->connectNodes(&nodes[next], next, &nodes[j], branchLength ,&nodes[next], branchLength);
		
		
		auto p = make_pair(index[j], index[next]);
		
		index[next] = n_seq+i;
	} 
}

