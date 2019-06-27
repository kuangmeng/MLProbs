////////////////////////////////////////////////////////////////
// AlignGraph.h
//
// Utilities for constructing the alignment graph.
/////////////////////////////////////////////////////////////////

#ifndef ALIGNGRAPH_H
#define ALIGNGRAPH_H

#include "SafeVector.h"
#include "MultiSequence.h"
#include "Sequence.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/////////////////////////////////////////////////////////////////
// AlignGraph
//
// Class for multiple sequence alignment input/output.
/////////////////////////////////////////////////////////////////

class AlignGraph {

	//Sequences in the alignment
	MultiSequence *sequences;

	//Deduced alignment
	MultiSequence *alignment;

	//The alignment graph
	// G[i][j] is the j'th child of the i'th node
	VVI G;

	//Columns(nodes) in Graph
	// cols[i]: vector of residues in node i
	VVVI cols;

	// Matrix of residues in each node:
	// IsPresent[i][j] is -1 if  j'th residue in i'th sequence
	// has not appeared in any node and is m (>=1) if it is
	// present in the  m'th node
	VVI IsPresent;

	//Matrices of Ancestors and Descendents:
	// Ancs[i][j]=true if node j is ancestor of node i;  otherwise it is false
	// Descs[i][j]=true if node j is descendant of node i;  otherwise it is false
	SafeVector<SafeVector<bool> > Ancs, Descs;

	//Zero Vector
	SafeVector<bool> ZZ;

private:

	/////////////////////////////////////////////////////////////////
	// Partition()
	//
	// Partition the set for Quick Sort
	/////////////////////////////////////////////////////////////////

	int Partition(int low, int high, float arr[], int ind[]) {
		float high_vac, low_vac, pivot;
		int il, ih, ip;
		ip = ind[low];
		il = ind[low];
		ih = ind[high];
		ih = ind[high];
		pivot = arr[low];

		while (high > low) {
			high_vac = arr[high];
			ih = ind[high];
			while (pivot <= high_vac) {
				if (high <= low)
					break;
				high--;
				high_vac = arr[high];
				ih = ind[high];
			}
			arr[low] = high_vac;
			ind[low] = ih;
			low_vac = arr[low];
			il = ind[low];
			while (pivot >= low_vac) {
				if (high <= low)
					break;
				low++;
				low_vac = arr[low];
				il = ind[low];
			}
			arr[high] = low_vac;
			ind[high] = il;
		}
		arr[low] = pivot;
		ind[low] = ip;
		return low;
	}
	/////////////////////////////////////////////////////////////////
	// Quick_sort()
	//
	// Quick Sort Algorithm
	/////////////////////////////////////////////////////////////////

	void Quick_sort(int low, int high, float arr[], int ind[]) {
		int Piv_index;
		if (low < high) {
			Piv_index = Partition(low, high, arr, ind);
			Quick_sort(low, Piv_index - 1, arr, ind);
			Quick_sort(Piv_index + 1, high, arr, ind);
		}
	}

	//////////////////////////////////////////////////////////////
	// merge()
	// Merge sorted lists A and B into list A.  A must have dim >= m+n
	// Added by Yongtao Ye, May 2014
	/////////////////////////////////////////////////////////////
	void merge(float A[], int Aind[], int m, float B[], int Bind[], int n) 
	{
		int i=0, j=0, k=0;
		int size = m+n;
		float *C = (float *)malloc(size*sizeof(float));
		int *Cind = (int *)malloc(size*sizeof(int));
		while (i < m && j < n) {
            if (A[i] <= B[j]) {
            	C[k] = A[i];
            	Cind[k] = Aind[i];
            	i++;
            }
            else {
            	C[k] = B[j];
            	Cind[k] = Bind[j];
            	j++;
            }
            k++;
		}
		if (i < m) 
			for (int p = i; p < m; p++,k++) {
				C[k] = A[p];
				Cind[k] = Aind[p];
			}
		else 
			for (int p = j; p < n; p++,k++) {
				C[k] = B[p];
				Cind[k] = Bind[p];
			}
		for( i=0; i<size; i++ ) {
			A[i] = C[i];
			Aind[i] = Cind[i];
		}
		free(C);
		free(Cind);
	}

	//////////////////////////////////////////////////////////////
	// arraymerge()
	// Merges N sorted sub-sections of array a into final, fully sorted array a
	// Added by Yongtao Ye, May 2014
	/////////////////////////////////////////////////////////////
	void arraymerge(float *a, int *aind, int size, int *index, int N)
	{
		int i;
		while ( N>1 ) {
	    	for( i=0; i<N; i++ ) index[i]=i*size/N; index[N]=size;
//#ifdef _OPENMP	    		
//pragma omp parallel for private(i) 
//#endif	    		
	    	for( i=0; i<N; i+=2 ) {
			//if( VERBOSE ) fprintf(stderr,"merging %d and %d, index %d and %d (up to %d)\n",i,i+1,index[i],index[i+1],index[i+2]);
				merge(a+index[i], aind+index[i], index[i+1]-index[i], 
					  a+index[i+1], aind+index[i+1], index[i+2]-index[i+1]);
			//if( VERBOSE ) for(int i=0; i<size; i++) fprintf(stderr,"after: %d %d\n",i,a[i]);
	    	}
	   		N /= 2;
		}
	}


	/////////////////////////////////////////////////////////////////
	// FindCloseNodes()
	//
	// Find the preceding and succeeding nodes in the graph G that
	// contain residues in the sequence of residue x
	/////////////////////////////////////////////////////////////////

	VVI FindCloseNodes(VI x) {

		VVI Sx;

		// The preceding node
		VI Parent;

		// The succeeding node
		VI Child;
		int pi;
		int ci;
		int clx = -1;
		int crx = 10000; //inf
		for (int i = x[1] - 1; i >= (int) 0; i--)
			if (IsPresent[x[0]][i] != -1) {
				pi = IsPresent[x[0]][i];
				clx = i;
				break;
			}
		for (int i = x[1] + 1; i < (int) IsPresent[x[0]].size(); i++)
			if (IsPresent[x[0]][i] != -1) {
				ci = IsPresent[x[0]][i];
				crx = i;
				break;
			}
		if (clx != -1)
			Parent.push_back(pi);
		if (crx != 10000)
			Child.push_back(ci);

		Sx.push_back(Parent);
		Sx.push_back(Child);

		return Sx;
	}

	/////////////////////////////////////////////////////////////////
	// GiveNodes()
	//
	// Give the ancestors or descendents of a given node
	/////////////////////////////////////////////////////////////////

	VI GiveNodes(SafeVector<bool> A, int msize) {

		VI B;
		for (int i = 0; i < msize; i++)
			if (A[i])
				B.push_back(i);
		return B;

	}

	/////////////////////////////////////////////////////////////////
	// Union_VBVI()
	//
	// Perform the union of set a (in the boolean format)
	// with set b  (in integer format)
	/////////////////////////////////////////////////////////////////

	SafeVector<bool> Union_VBVI(SafeVector<bool> A, VI B) {
		SafeVector<bool> C(A);
		for (int i = 0; i < (int) B.size(); i++) {
			C[B[i]] = true;
		}
		return C;
	}

	/////////////////////////////////////////////////////////////////
	// Union()
	//
	// Perform the union of sets a and b ( both in the boolean format)
	/////////////////////////////////////////////////////////////////

	SafeVector<bool> Union(SafeVector<bool> A, SafeVector<bool> B, int msize) {
		SafeVector<bool> C(A);
		for (int i = 0; i < msize; i++) {
			C[i] = A[i] | B[i];
		}
		return C;
	}

	/////////////////////////////////////////////////////////////////
	// Union_VI()
	//
	// Perform the union of sets a and b ( both in the integer format)
	/////////////////////////////////////////////////////////////////

	VI Union_VI(VI A, VI B) {

		for (int j = 0; j < (int) B.size(); j++) {
			bool flag = false;
			for (int k = 0; k < (int) A.size(); k++)
				if (A[k] == B[j]) {
					flag = true;
					break;
				}
			if (!flag)
				A.push_back(B[j]);
		}
		return A;
	}

	/////////////////////////////////////////////////////////////////
	// Ismember()
	//
	// Check whether i is in set A (boolean format)
	/////////////////////////////////////////////////////////////////

	bool Ismember(SafeVector<bool> A, int i) {
		return (A[i]);
	}

	/////////////////////////////////////////////////////////////////
	// Ismember_VI()
	//
	// Check whether i is in set A (integer format)
	/////////////////////////////////////////////////////////////////

	bool Ismember_VI(VI A, int i) {
		for (int j = 0; j < (int) A.size(); j++)
			if (A[j] == i)
				return true;
		return false;
	}

	/////////////////////////////////////////////////////////////////
	// Remove()
	//
	// Remove i from set A (integer format)
	/////////////////////////////////////////////////////////////////

	VI Remove(VI A, int i) {
		VI B;
		for (int j = 0; j < (int) A.size(); j++)
			if (A[j] != i)
				B.push_back(A[j]);
		return B;
	}

	/////////////////////////////////////////////////////////////////
	// Update()
	//
	// Update the set A after we merges to components cx and cy
	/////////////////////////////////////////////////////////////////

	SafeVector<bool> Update(SafeVector<bool> A, int cy, int msize) {

		for (int j = cy; j < msize - 1; j++)
			A[j] = A[j + 1];
		A[msize - 1] = false;
		return A;
	}

	/////////////////////////////////////////////////////////////
	// Update()
	//
	// Update the value of i after we merges to components cx and cy
	/////////////////////////////////////////////////////////////////

	int Update(int i, int cx, int cy) {
		if (i < cy)
			return i;
		else if (i == cy)
			return cx;
		else
			return i - 1;
	}

	/////////////////////////////////////////////////////////////
	// GiveParent()
	//
	// Find the parents of the node i in graph G
	/////////////////////////////////////////////////////////////////

	VI GiveParent(VVI G, int i) {
		VI a;
		for (int j = 0; j < (int) G.size(); j++) {
			for (int k = 0; k < (int) G[j].size(); k++) {
				if (G[j][k] == i) {
					a.push_back(j);
					break;
				}
			}
		}
		return a;
	}

	/////////////////////////////////////////////////////////////
	// CheckAddNewNode()
	//
	// Procedure to insert the new node in the graph provided that
	// it does not introduce any cycle
	/////////////////////////////////////////////////////////////////

	bool CheckAddNewNode(VI x, VI y) {

		//Develop the temporary Graph to check for cycles
		VVI Temp_G = G;

		//Find the closest nodes for residues x and y
		VVI Sx = FindCloseNodes(x);
		VVI Sy = FindCloseNodes(y);

		// Find the set of Parents and Children
		VI Parent = Union_VI(Sx[0], Sy[0]);
		VI Child = Union_VI(Sx[1], Sy[1]);

		//Add the Parents and Children to the graph
		Temp_G.push_back(Child);

		for (int j = 0; j < (int) Parent.size(); j++)
			Temp_G[Parent[j]].push_back(G.size());

		//Check for cycle
		bool Check_Cycle = true;
		//Check whether one residue's parent is the descendent of the other one
		if ((Sx[0].size() == 1) & (Sy[1].size() == 1))
			Check_Cycle = Check_Cycle & !Descs[Sy[1][0]][Sx[0][0]] & !(Sx[0][0]
					== Sy[1][0]);
		if ((Sy[0].size() == 1) & (Sx[1].size() == 1))
			Check_Cycle = Check_Cycle & !Descs[Sx[1][0]][Sy[0][0]] & !(Sy[0][0]
					== Sx[1][0]);

		// If no cycle appears update the graph
		if (Check_Cycle) {

			// Remove redundant edges
			if ((Sy[0].size() == 1) & (Sx[0].size() == 1)) {
				if (Ismember(Descs[Sx[0][0]], Sy[0][0]))
					Temp_G[Sx[0][0]] = Remove(Temp_G[Sx[0][0]], G.size());
				if (Ismember(Descs[Sy[0][0]], Sx[0][0]))
					Temp_G[Sy[0][0]] = Remove(Temp_G[Sy[0][0]], G.size());
			}
			if ((Sy[1].size()) & (Sx[1].size())) {
				if (Ismember(Descs[Sx[1][0]], Sy[1][0]))
					Temp_G[G.size()] = Remove(Temp_G[G.size()], Sy[1][0]);
				if (Ismember(Descs[Sy[1][0]], Sx[1][0]))
					Temp_G[G.size()] = Remove(Temp_G[G.size()], Sx[1][0]);
			}
			for (int j = 0; j < (int) Parent.size(); j++)
				for (int k = 0; k < (int) Child.size(); k++)
					Temp_G[Parent[j]] = Remove(Temp_G[Parent[j]], Child[k]);

			//Update IsPresent
			IsPresent[x[0]][x[1]] = G.size();
			IsPresent[y[0]][y[1]] = G.size();

			//Update G
			G = Temp_G;

			/*	SafeVector <VI*>  tt;
			 SafeVector<VI* >  tt2;
			 VI *t;
			 tt.push_back(t);
			 delete tt[0];
			 tt=tt2;
			 */

			int Gsz = G.size();

			//New Node Ancestors
			SafeVector<bool> curr_Anc(ZZ);
			if (Parent.size() > 0)
				curr_Anc = Ancs[Parent[0]];
			if (Parent.size() == 2)
				curr_Anc = Union(curr_Anc, Ancs[Parent[1]], Gsz - 1);
			curr_Anc = Union_VBVI(curr_Anc, Parent);

			//Update Ancestors for the new node
			Ancs.push_back(curr_Anc);

			//New Node descendants
			SafeVector<bool> curr_Desc(ZZ);
			if (Child.size() > 0)
				curr_Desc = Descs[Child[0]];
			if (Child.size() == 2)
				curr_Desc = Union(curr_Desc, Descs[Child[1]], Gsz - 1);
			curr_Desc = Union_VBVI(curr_Desc, Child);

			//Update Descendants for the new node
			Descs.push_back(curr_Desc);

			//Update Ancestors and Descendants for other nodes
			VI DD, AA;
			for (int j = 0; j < Gsz; j++) {
				if (curr_Anc[j])
					AA.push_back(j);
				if (curr_Desc[j])
					DD.push_back(j);
			}
			for (int j = 0; j < (int) DD.size(); j++) {
				Ancs[DD[j]][Gsz - 1] = true;
				for (int k = 0; k < (int) AA.size(); k++) {
					Ancs[DD[j]][AA[k]] = true;
					Descs[AA[k]][DD[j]] = true;
				}
			}
			for (int k = 0; k < (int) AA.size(); k++)
				Descs[AA[k]][Gsz - 1] = true;
		}
		return Check_Cycle;
	}

	/////////////////////////////////////////////////////////////
	// CheckAddColumnEx()
	//
	// Procedure to extend an existing column in the graph provided that
	// it does not introduce any cycle
	/////////////////////////////////////////////////////////////////

	bool CheckAddColumnEx(VI y, int cx) {

		//Develop the temporary Graph to check for cycles
		VVI Temp_G = G;

		//Find the closest nodes for residue y
		VVI Sy = FindCloseNodes(y);

		// Find the set of Parents and Children
		VI Parent = Sy[0];
		VI Child = Sy[1];

		// If x is not already the child of parents of y insert
		// that as a child of them
		for (int j = 0; j < (int) Parent.size(); j++) {
			bool flag = false;
			for (int k = 0; k < (int) Temp_G[Parent[j]].size(); k++)
				if (Temp_G[Parent[j]][k] == cx) {
					flag = true;
					break;
				}
			if (!flag)
				Temp_G[Parent[j]].push_back(cx);
		}

		// If x is not already the parent of children of y insert
		// that as a parent of them
		for (int j = 0; j < (int) Child.size(); j++) {
			bool flag = false;
			for (int k = 0; k < (int) Temp_G[cx].size(); k++)
				if (Temp_G[cx][k] == Child[j]) {
					flag = true;
					break;
				}
			if (!flag)
				Temp_G[cx].push_back(Child[j]);
		}

		//Check for cycle
		bool Check_Cycle = true;
		//Check whether x is the descendant of any of y's children
		//or if any Descendants of x is a parent of y
		if (Child.size() > 0)
			Check_Cycle = !Descs[Child[0]][cx] & !(Child[0] == cx);
		if (Parent.size() > 0)
			Check_Cycle = Check_Cycle & !(Descs[cx][Parent[0]]) & !(Parent[0]
					== cx);

		// If no cycle appears update the graph
		if (Check_Cycle) {

			// Remove redundant edges
			if (Sy[0].size() == 1)
				if ((Ismember(Descs[Sy[0][0]], cx)) & (!Ismember_VI(
						G[Sy[0][0]], cx)))
					Temp_G[Sy[0][0]] = Remove(Temp_G[Sy[0][0]], cx);
			if (Sy[1].size() == 1)
				if ((Ismember(Descs[cx], Sy[1][0])) && (!Ismember_VI(G[cx],
						Sy[1][0])))
					Temp_G[cx] = Remove(Temp_G[cx], Sy[1][0]);
			if ((Sy[1].size() == 1) & (Sy[0].size() == 1))
				Temp_G[Sy[0][0]] = Remove(Temp_G[Sy[0][0]], Sy[1][0]);

			//Update IsPresent
			IsPresent[y[0]][y[1]] = cx;

			//Update G
			G = Temp_G;
			int Gsz = G.size();

			//cx Ancestors
			SafeVector<bool> curr_Anc(ZZ);
			if (Parent.size() > 0)
				curr_Anc = Ancs[Parent[0]];
			curr_Anc = Union_VBVI(curr_Anc, Parent);

			//Update cx's Ancestors
			Ancs[cx] = Union(Ancs[cx], curr_Anc, Gsz);

			//cx Descendants
			SafeVector<bool> curr_Desc(ZZ);
			if (Child.size() > 0)
				curr_Desc = Descs[Child[0]];
			curr_Desc = Union_VBVI(curr_Desc, Child);

			//Update cx's Descendants
			Descs[cx] = Union(Descs[cx], curr_Desc, Gsz);

			//Update Ancestors and Descendants for other nodes
			VI DD, AA;
			for (int j = 0; j < Gsz; j++) {
				if (Ancs[cx][j])
					AA.push_back(j);
				if (Descs[cx][j])
					DD.push_back(j);
			}
			for (int j = 0; j < (int) DD.size(); j++) {
				Ancs[DD[j]][cx] = true;
				for (int k = 0; k < (int) AA.size(); k++) {
					Ancs[DD[j]][AA[k]] = true;
					Descs[AA[k]][DD[j]] = true;
				}
			}
			for (int k = 0; k < (int) AA.size(); k++)
				Descs[AA[k]][cx] = true;

		}
		return Check_Cycle;
	}

	/////////////////////////////////////////////////////////////
	// CheckAddColumnMrg()
	//
	// Procedure to Merge two columns in the graph provided that
	// it does not introduce any cycle
	/////////////////////////////////////////////////////////////////

	bool CheckAddColumnMrg(int cx, int cy) {

		//Develop the temporary Graph to check for cycles
		VVI Temp_G;

		// Find the set of Children
		VI Child_x = G[cx];
		VI Child_y = G[cy];
		VI Child;
		Child_x = Union_VI(Child_x, Child_y);
		for (int j = 0; j < (int) Child_x.size(); j++)
			Child.push_back(Update(Child_x[j], cx, cy));

		//Construct the updated temporary Graph
		for (int j = 0; j < (int) G.size(); j++) {
			bool flag = false;
			VI nodes;
			for (int k = 0; k < (int) G[j].size(); k++)
				if ((G[j][k] == cx) | (G[j][k] == cy)) {
					if (flag == false) {
						nodes.push_back(cx);
						flag = true;
					}
				} else if (G[j][k] < cy)
					nodes.push_back(G[j][k]);
				else if (G[j][k] > cy)
					nodes.push_back(G[j][k] - 1);
			if (j == cx)
				Temp_G.push_back(Child);
			else if (j != cy)
				Temp_G.push_back(nodes);

		}

		//Check for cycle
		bool Check_Cycle;
		//Check whether any descendant of x is y or visevers
		Check_Cycle = !(Descs[cx][cy]) & !(Descs[cy][cx]);

		// If no cycle appears update the graph
		if (Check_Cycle) {

			// Remove redundant edges
			VI ax = GiveNodes(Ancs[cx], G.size());
			VI dy = GiveNodes(Descs[cy], G.size());
			for (int j = 0; j < (int) ax.size(); j++)
				for (int k = 0; k < (int) dy.size(); k++)
					if (Ismember_VI(G[ax[j]], dy[k]))
						Temp_G[Update(ax[j], cx, cy)] = Remove(Temp_G[Update(
								ax[j], cx, cy)], Update(dy[k], cx, cy));
			for (int j = 0; j < (int) ax.size(); j++)
				if (Ismember_VI(G[ax[j]], cy) & !Ismember_VI(G[ax[j]], cx))
					Temp_G[Update(ax[j], cx, cy)] = Remove(Temp_G[Update(ax[j],
							cx, cy)], cx);
			VI ay = GiveNodes(Ancs[cy], G.size());
			VI dx = GiveNodes(Descs[cx], G.size());
			for (int j = 0; j < (int) ay.size(); j++)
				for (int k = 0; k < (int) dx.size(); k++)
					if (Ismember_VI(G[ay[j]], dx[k])) {
						Temp_G[Update(ay[j], cx, cy)] = Remove(Temp_G[Update(
								ay[j], cx, cy)], Update(dx[k], cx, cy));
					}
			for (int j = 0; j < (int) ay.size(); j++)
				if (Ismember_VI(G[ay[j]], cx) & !Ismember_VI(G[ay[j]], cy)) {
					Temp_G[Update(ay[j], cx, cy)] = Remove(Temp_G[Update(ay[j],
							cx, cy)], cx);
				}
			VI pax = GiveParent(G, cx);
			VI chx = G[cx];
			VI pay = GiveParent(G, cy);
			VI chy = G[cy];
			for (int j = 0; j < (int) pax.size(); j++)
				if ((Ismember_VI(ay, pax[j])) & !(Ismember_VI(G[pax[j]], cy)))
					Temp_G[Update(pax[j], cx, cy)] = Remove(Temp_G[Update(
							pax[j], cx, cy)], cx);
			for (int j = 0; j < (int) pay.size(); j++)
				if ((Ismember_VI(ax, pay[j])) & !(Ismember_VI(G[pay[j]], cx)))
					Temp_G[Update(pay[j], cx, cy)] = Remove(Temp_G[Update(
							pay[j], cx, cy)], cx);
			for (int j = 0; j < (int) chx.size(); j++)
				if ((Ismember_VI(dy, chx[j])) & !(Ismember_VI(G[cy], chx[j])))
					Temp_G[cx] = Remove(Temp_G[cx], Update(chx[j], cx, cy));
			for (int j = 0; j < (int) chy.size(); j++)
				if ((Ismember_VI(dx, chy[j])) & !(Ismember_VI(G[cx], chy[j])))
					Temp_G[cx] = Remove(Temp_G[cx], Update(chy[j], cx, cy));

			//Update IsPresent
			for (int j = 0; j < (int) IsPresent.size(); j++)
				for (int k = 0; k < (int) IsPresent[j].size(); k++)
					IsPresent[j][k] = Update(IsPresent[j][k], cx, cy);

			//Update G
			G = Temp_G;
			int Gsz = G.size();

			//Merged node Ancestors and Descendants
			SafeVector<bool> curr_Anc(ZZ), curr_Desc(ZZ);
			curr_Anc = Union(Ancs[cx], Ancs[cy], Gsz + 1);
			curr_Desc = Union(Descs[cx], Descs[cy], Gsz + 1);

			//Update Ancestors and Descendants for other nodes
			SafeVector<SafeVector<bool> > new_Ancs, new_Descs;
			for (int j = 0; j < (int) Ancs.size(); j++) {
				if (j == cx) {
					curr_Anc = Update(curr_Anc, cy, Gsz + 1);
					new_Ancs.push_back(curr_Anc);
				} else if (j != cy) {
					SafeVector<bool> n_anc = Update(Ancs[j], cy, Gsz + 1);
					new_Ancs.push_back(n_anc);
				}
			}
			Ancs = new_Ancs;
			for (int j = 0; j < (int) Descs.size(); j++) {
				if (j == cx) {
					curr_Desc = Update(curr_Desc, cy, Gsz + 1);
					new_Descs.push_back(curr_Desc);
				} else if (j != cy) {
					SafeVector<bool> n_Desc = Update(Descs[j], cy, Gsz + 1);
					new_Descs.push_back(n_Desc);
				}
			}
			Descs = new_Descs;
			VI DD, AA;
			for (int j = 0; j < Gsz; j++) {
				if (curr_Anc[j])
					AA.push_back(j);
				if (curr_Desc[j])
					DD.push_back(j);
			}
			for (int j = 0; j < (int) DD.size(); j++) {
				Ancs[DD[j]][cx] = true;
				for (int k = 0; k < (int) AA.size(); k++) {
					Ancs[DD[j]][AA[k]] = true;
					Descs[AA[k]][DD[j]] = true;
				}
			}
			for (int k = 0; k < (int) AA.size(); k++)
				Descs[AA[k]][cx] = true;

		}
		return Check_Cycle;
	}

	/////////////////////////////////////////////////////////////
	// AddtoPath()
	//
	// Add node N2 after node N1 in the path
	/////////////////////////////////////////////////////////////////

	void AddtoPath(VI &Path, int N1, int N2) {

		int Ps = Path.size();
		int h;
		if (N1 == -1) {
			// If N2 is a root node
			h = -1;
		} else {
			for (h = 0; h < Ps; h++)
				if (Path[h] == N1)
					break;
		}
		if (h == (Ps - 1))
			Path.push_back(N2);
		else {
			Path.push_back(Path[Ps - 1]);
			for (int k = Ps - 1; k > (h + 1); k--)
				Path[k] = Path[k - 1];
			Path[h + 1] = N2;
		}
	}

	/////////////////////////////////////////////////////////////
	// FindPath()
	//
	// Update the path starting from the node N1 and considering the
	// passed nodes in "marked"
	/////////////////////////////////////////////////////////////////

	void FindPath(int N1, SafeVector<bool> &marked, VI &Path) {

		for (int j = 0; j < (int) G[N1].size(); j++) {
			if (!marked[G[N1][j]]) {
				marked[G[N1][j]] = true;
				AddtoPath(Path, N1, G[N1][j]);
				FindPath(G[N1][j], marked, Path);
			}
		}
	}

	/////////////////////////////////////////////////////////////
	// Path2Align()
	//
	// Convert path to a valid alignment
	/////////////////////////////////////////////////////////////////
	void Path2Align(VI Path, VVVI SRC, VVI ZeroPos, int SRCsize) {

		//Construct the alignment structure
		alignment = new MultiSequence();
		VI lengths(sequences->GetNumSequences());
		for (int i = 0; i < sequences->GetNumSequences(); i++) {

			Sequence *seq = sequences->GetSequence(i);
			Sequence *ret = new Sequence();
			assert (ret);

			ret->isValid = seq->isValid;
			ret->header = seq->header;
			ret->data = new SafeVector<char> ;
			assert (ret->data);
			ret->length = (int) Path.size() + SRCsize;
			ret->sequenceLabel = seq->sequenceLabel;
			ret->inputLabel = seq->inputLabel;
			ret->data->push_back('@');
			alignment->AddSequence(ret);
		}

		VI temp_seqs(sequences->GetNumSequences());
		for (int i = 0; i < sequences->GetNumSequences(); i++)
			temp_seqs[i] = i;

		//insert single residue columns which should appear
		//at position zero
		for (int i = 0; i < (int) ZeroPos.size(); i++) {
			int seqnum = ZeroPos[i][0];
			int resnum = ZeroPos[i][1];
			alignment->GetSequence(seqnum)->data->push_back(
					sequences->GetSequence(seqnum)->GetPosition(resnum + 1));
			for (int k = 0; k < sequences->GetNumSequences(); k++)
				if (k != seqnum)
					alignment->GetSequence(k)->data->push_back('-');
		}

		//insert other columns
		for (int i = 0; i < (int) Path.size(); i++) {
			VI b(temp_seqs);
			//insert residues
			for (int j = 0; j < (int) cols[Path[i]].size(); j++) {
				int seqnum = cols[Path[i]][j][0];
				int resnum = cols[Path[i]][j][1];
				alignment->GetSequence(seqnum)->data->push_back(
						sequences->GetSequence(seqnum)->GetPosition(resnum + 1));
				b = Remove(b, seqnum);
			}
			//insert gaps
			for (int j = 0; j < (int) b.size(); j++) {
				alignment->GetSequence(b[j])->data->push_back('-');
			}

			//insert single residue columns which should
			// appear after i'th rgular column
			for (int j = 0; j < (int) SRC[i].size(); j++) {
				int seqnum = SRC[i][j][0];
				int resnum = SRC[i][j][1];
				alignment->GetSequence(seqnum)->data->push_back(
						sequences->GetSequence(seqnum)->GetPosition(resnum + 1));
				for (int k = 0; k < sequences->GetNumSequences(); k++)
					if (k != seqnum)
						alignment->GetSequence(k)->data->push_back('-');
			}
		}
	}

public:

	/////////////////////////////////////////////////////////////////
	// AlignGraph::AlignGraph()
	//
	// Default constructor.
	/////////////////////////////////////////////////////////////////

	AlignGraph() {
	}

	/////////////////////////////////////////////////////////////////
	// AlignGraph::ConstGraph()
	//
	// Construct the Graph From Postrior Probabilities.
	// Modified by Yongtao Ye, May 2014
	/////////////////////////////////////////////////////////////////

	AlignGraph(VVI aligns, VF alignp, MultiSequence *sequences) :
		sequences(sequences) {

		int numSeqs = sequences->GetNumSequences();

		//sort the posterior probabilities using quick sort algorithm
		int n, low, high;
		float *a;
		int *ind;
		n = alignp.size();
		a = new float[n];
		ind = new int[n];
		for (int i = 0; i < n; i++) {
			a[i] = alignp[i];
			ind[i] = i;
		}
		high = n - 1;
		low = 0;

		int threads = 1;
#ifdef _OPENMP
		threads = omp_get_num_threads();
#endif
		if(threads < 2) Quick_sort(low, high, a, ind);
		else{
			int *index = (int *)malloc((threads+1)*sizeof(int));
			for(int i=0; i<threads; i++) index[i]=i*n/threads; 
			index[threads]=n;
//#ifdef _OPENMP
//#pragma omp parallel for private(i)
//#endif
			//sub-sort 
			for(int i=0; i<threads; i++) 
				Quick_sort(index[i], index[i+1]-1, a, ind);
			//merge
			arraymerge(a, ind, n, index, threads);
			free(index);
		}	

		// Length of the longest sequence
		int maxlength = 0;
		for (int i = 0; i < sequences->GetNumSequences(); i++)
			if (sequences->GetSequence(i)->GetLength() > maxlength)
				maxlength = sequences->GetSequence(i)->GetLength();

		// initialize the zero vector:
		// initial guess of number of columns in the alignment is 1.5*max_length
		for (int i = 0; i < 1.5 * maxlength; i++)
			ZZ.push_back(false);

		//initialize IsPresent (Matrix of residues in each node)
		for (int i = 0; i < numSeqs; i++) {
			VI s;
			for (int j = 0; j < maxlength; j++) {
				s.push_back(-1);
			}
			IsPresent.push_back(s);
		}

		// initialize counters
		int cnt_NN = 0; //number of new nodes addition calls
		int cnt_CE = 0; //number of column extension calls
		int cnt_CM = 0; //number of column merging calls
		int cnt_tot = 0; // totla number of calls
		int cnt_NoCycle = 0;// number of calls with no cycle introduction

		//Start Graph construction from the largest posteriro probabilities
		for (int i = 0; i < (int) aligns.size(); i++) {
			//residue pair to be aligned
			VI res_pair(aligns[ind[aligns.size() - i - 1]]);

			//residues x and y
			// x[0] sequence number of residue x
			// x[1] residue position of residue x
			VI x(2);
			x[0] = res_pair[0];
			x[1] = res_pair[1];
			VI y(2);
			y[0] = res_pair[2];
			y[1] = res_pair[3];

			//Check whether x and/or y have appeared in any node
			int cx = -1;
			int cy = -1;
			bool fx = false;
			bool fy = false;
			if (IsPresent[x[0]][x[1]] != -1) {
				fx = true;
				cx = IsPresent[x[0]][x[1]];
			}
			if (IsPresent[y[0]][y[1]] != -1) {
				fy = true;
				cy = IsPresent[y[0]][y[1]];
			}

			bool Check_Cycle = false;

			//if niether x nor y have appeared:  Add a new node
			if ((!fx) & (!fy)) {
				cnt_tot++;
				cnt_NN++;
				Check_Cycle = CheckAddNewNode(x, y);

				if (Check_Cycle)
					cnt_NoCycle++;
			}

			//if only one of x or y have appeared:  Extend the column
			else if (fx ^ fy) {
				//make x the residue which is already in the graph
				if (fy) {
					swap(x, y);
					swap(cx, cy);
					swap(fx, fy);
				}

				//immediate check for cycle : check whether cx
				// has common residue in sequence of residue y
				bool iflag = false;
				for (int j = 0; j < (int) IsPresent[y[0]].size(); j++)
					if (IsPresent[y[0]][j] == cx) {
						iflag = true;
						break;
					}
				if (iflag == false) {
					cnt_tot++;
					cnt_CE++;
					Check_Cycle = CheckAddColumnEx(y, cx);
					if (Check_Cycle)
						cnt_NoCycle++;
				}
			}

			//if both of x and y have appeared: Merge the columns
			else if ((fx & fy) & (cx != cy)) { //Column Merge
				//immediate check for cycle : check whether cx
				// has common residue in sequence of residue y or vice versa
				bool iflag = false;
				for (int j = 0; j < (int) IsPresent[y[0]].size(); j++)
					if (IsPresent[y[0]][j] == cx) {
						iflag = true;
						break;
					}
				for (int j = 0; j < (int) IsPresent[x[0]].size(); j++)
					if (IsPresent[x[0]][j] == cy) {
						iflag = true;
						break;
					}
				if (iflag == false) {
					if (cx > cy) {
						swap(x, y);
						swap(cx, cy);
						swap(fx, fy);
					}
					cnt_tot++;
					cnt_CM++;
					Check_Cycle = CheckAddColumnMrg(cx, cy);
					if (Check_Cycle)
						cnt_NoCycle++;
				}
			}

			// If the initial guess of the number of nodes
			// is close to reach extend number of possible nodes
			if (!((fx & fy) & (cx == cy))) {
				if (G.size() > ZZ.size() - 10) {
					for (int j = 0; j < 100; j++)
						ZZ.push_back(false);
					for (int k = 0; k < (int) G.size(); k++)
						for (int j = 0; j < 100; j++) {
							Descs[k].push_back(false);
							Ancs[k].push_back(false);
						}
				}
			}
		}

		//Construct the columns
		// cols[i]: vector of residues in node i

		VVI temp;
		for (int i = 0; i < (int) G.size(); i++)
			cols.push_back(temp);

		for (int i = 0; i < sequences->GetNumSequences(); i++)
			for (int j = 0; j < sequences->GetSequence(i)->GetLength(); j++)
				if (IsPresent[i][j] != -1) {
					VI x;
					x.push_back(i);
					x.push_back(j);
					cols[IsPresent[i][j]].push_back(x);
				}

		delete []a;
		delete []ind;
	}

	/////////////////////////////////////////////////////////////////
	// AlignGraph::Graph2Align()
	//
	// Convert the graph to a valid alignment
	/////////////////////////////////////////////////////////////////
	void Graph2Align() {

		//vector of roots in G
		VI roots;
		for (int i = 0; i < (int) G.size(); i++)
			if (GiveParent(G, i).size() == 0)
				roots.push_back(i);

		//vector of passed nodes
		SafeVector<bool> marked(G.size());

		//Alignment Path deduced from G
		VI Path;
		for (int i = 0; i < (int) roots.size(); i++) {
			AddtoPath(Path, -1, roots[i]);
			FindPath(roots[i], marked, Path);
		}

		//Mapping the place of each node in Path
		VI PathMap(Path.size());
		for (int i = 0; i < (int) Path.size(); i++) {
			PathMap[Path[i]] = i;
		}

		// Single residue columns which should appear after
		// each node in Path
		VVVI SRC(Path.size());

		// single residue columns which should appear
		//at position zero
		VVI ZeroPos;

		//figure out the single residue columns
		int SRCsize = 0;
		for (int i = 0; i < sequences->GetNumSequences(); i++)
			for (int j = 0; j < sequences->GetSequence(i)->GetLength(); j++)
				if (IsPresent[i][j] == -1) {
					SRCsize++;
					VI x(2);
					x[0] = i;
					x[1] = j;
					int ct = j - 1;
					while (ct >= 0) {
						if (IsPresent[i][ct] != -1) {
							SRC[PathMap[IsPresent[i][ct]]].push_back(x);
							break;
						}
						ct--;
					}
					if (ct == -1)
						ZeroPos.push_back(x);
				}

		//Convert the path to valid alignment
		Path2Align(Path, SRC, ZeroPos, SRCsize);

	}

	/////////////////////////////////////////////////////////////////
	// AlignGraph::GetAlignment()
	//
	// Retrieve the alignment from the AlignGraph object.
	/////////////////////////////////////////////////////////////////

	MultiSequence* GetAlignment() {
		return alignment;
	}

	/////////////////////////////////////////////////////////////////
	// AlignGraph::~AlignGraph()
	//
	// Destructor.  Gets rid of objects contained in the graph.
	/////////////////////////////////////////////////////////////////

	~AlignGraph() {
		/*
		if (sequences) {
			delete sequences;
			sequences = NULL;
		}
		
		if (alignment) {
			delete alignment;
			alignment = NULL;
		}
		*/
	}
};
#endif
