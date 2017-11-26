#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>
#include <string.h>

//#define __DEBUG__	// to create files for debugging

#define EPSILON 0.00000000001		// error threshold of computing floating point

#define POPSIZE 100			// population size in GA
#define OFFSIZE 50			// offsping size in GA

using namespace std;

typedef struct {			// subfragments
	char *str;
	double *qStr;
	double *qStr_star;
	int start_pos;
	int end_pos;
	int length;
}SubFragType;

typedef struct {			// fragment(read)
	SubFragType * subFragment;
	int num_subFrag;
	int start_pos;
	int end_pos;
	int length;
}FragType;

typedef struct {			// fragment inforamtion for each position
	int start_frag;
	int end_frag;
	int endPos;
}FragsForPosType;

typedef struct {			// block
	int start_frag;
	int end_frag;
	int length;
}BlockType;

typedef struct {			// Weighted MEC (distance) for each fragment
	double D_hf;			// It is used in toggling stage
	double D_hfstar;		// in order to speed up calculation of weighted MEC
}DFragType;

typedef struct {			// Haplotype sequeucne
	char *seq, *seq2;
	double sum_D;
}HaploType;



/////////////////////////////////////////////////////////////////

char *MatrixData;				// raw matrix data
FragType *Fragment;
FragsForPosType *FragsForPos;
DFragType *DFrag;

vector <BlockType> Block;		// block of fragments

HaploType *Population;			// GA population
int *GAcnt;						// counter used in GA
int *GAcnt2;

HaploType BestHaplo, TmpHaplo;				// haplotype sequences
char *Covered;						// positions covered by fragments

int NumPos;						// number of positions
int NumFrag;						// number of fragments
int MaxBlkLen;					// maximum block length

/////////////////////////////////////////////////////////////////

int calc_MEC(HaploType &, BlockType &);
double calc_sumD(HaploType &, BlockType &, int);
double calc_sumD_tog(HaploType &, BlockType &, int, int);
double calc_sumD_2seq(HaploType &, BlockType &);

void find_covered_sites(BlockType &);
void GA(HaploType &, BlockType &);
int range_switch_toggling(HaploType &, BlockType &);
int single_switch_toggling(HaploType &, BlockType &);
int create_complement(HaploType &, BlockType &);
void find_homoSNP(HaploType &, BlockType &);

int haplotype_phasing(BlockType &, int);

void load_matrixFile(char []);
void build_aux_struct();
void allocate();
void deallocate();

void procedure(char [], char [], int);

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

// compare function for sorting
bool compare_sumD_val(HaploType left, HaploType right) {
	return left.sum_D < right.sum_D;
}

bool compare_frag_pos(FragType left, FragType right) {
	return left.start_pos < right.start_pos;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

//
// calculating non-weighted MEC
// This function is used only for analysis, but not used in algorithm
//
int calc_MEC(HaploType& haplo, BlockType& block) {

	int startFrag, endFrag;
	int sum_D = 0;
	int pos_k;

	startFrag = block.start_frag;
	endFrag = block.end_frag;

	int realStartIndex = Fragment[block.start_frag].start_pos;

	for (int i = startFrag; i <= endFrag; i++) {	// for each fragment

		int D_hf = 0, D_hstarf = 0;

		for (int j = 0; j<Fragment[i].num_subFrag; j++) {		// for each subfragment
			for (int k = 0; k<Fragment[i].subFragment[j].length; k++) {		// for each position

				pos_k = k + Fragment[i].subFragment[j].start_pos - realStartIndex;

				if (haplo.seq[pos_k] == '-' || haplo.seq2[pos_k] == '-'
					|| Fragment[i].subFragment[j].str[k] == '-')
					continue;

				if (haplo.seq[pos_k] != Fragment[i].subFragment[j].str[k])
					++D_hf;
				if (haplo.seq2[pos_k] != Fragment[i].subFragment[j].str[k])
					++D_hstarf;

			}  // end_for k (each position)
		}  // end_for j (each subfragment)

		if (D_hf < D_hstarf)
			sum_D += D_hf;
		else
			sum_D += D_hstarf;

	}// end_for i (each fragment)

	return sum_D;
}


//
// calculating fitness values (normal version)
// : if 'update' is true, DFrag array and haplo (except for seq/seq_star) are updated
//
double calc_sumD(HaploType& haplo, BlockType& block, int update) {

	int startFrag, endFrag;
	double sum_D = 0.0;
	int pos_k;

	startFrag = block.start_frag;
	endFrag = block.end_frag;

	int realStartIndex = Fragment[block.start_frag].start_pos;

	for (int i = startFrag; i <= endFrag; i++) {	// for each fragment

		double D_hf = 0.0, D_hfstar = 0.0;

		for (int j = 0; j<Fragment[i].num_subFrag; j++) {	// for each subfragment
			for (int k = 0; k<Fragment[i].subFragment[j].length; k++) {	// for each position

				pos_k = k + Fragment[i].subFragment[j].start_pos - realStartIndex;

				double q_j = Fragment[i].subFragment[j].qStr[k];
				double q_j_star = Fragment[i].subFragment[j].qStr_star[k];

				// calculating distance for a position
				if (haplo.seq[pos_k]	!= Fragment[i].subFragment[j].str[k]) {
					D_hf += q_j_star;
					D_hfstar += q_j;
				}
				else {
					D_hf += q_j;
					D_hfstar += q_j_star;
				}

			} // end_for k (each position)

		} // end_for j (each subfragment)

		if (D_hf < D_hfstar)	// select min( D_h, D_h*)
			sum_D += D_hf;
		else
			sum_D += D_hfstar;
		
		if (update) {					// *** if update is tree  ***
			DFrag[i].D_hf = D_hf;			// *** the calculated values are stored in DFrag ***
			DFrag[i].D_hfstar = D_hfstar;
		}
	}// end_for i (each fragment)

	haplo.sum_D = sum_D;

	return sum_D;
}

//
// calculating fitness values (position version) used only in xxx_switch_procedure( )
// : calculating weighted MEC with assumption that seq[i] is toggled
// : recalculating the distance only at position pos
// : if 'update' is true, DFrag array and haplotype (except for seq/seq_star) are updated
//
double calc_sumD_tog(HaploType& haplo, BlockType& block, int pos, int update) {

	int startFrag, endFrag;
	double sum_D = haplo.sum_D;

	if (FragsForPos[pos + Fragment[block.start_frag].start_pos].start_frag
		> FragsForPos[pos + Fragment[block.start_frag].start_pos].end_frag)	// no fragment is located at pos
		return sum_D;

	pos = pos + Fragment[block.start_frag].start_pos; // convert pos to real position (matrix positino)
	startFrag = FragsForPos[pos].start_frag;	// first fragment located at pos
	endFrag = FragsForPos[pos].end_frag;		// last fragment located at pos

	for (int i = startFrag; i <= endFrag; i++) {	// for each fragment located at pos

		// 1. finding subfragment located at pos
		int j = 0;
		while (j < Fragment[i].num_subFrag && pos > Fragment[i].subFragment[j].end_pos)
			++j;								// skip subfragments before pos
		if (j >= Fragment[i].num_subFrag || pos < Fragment[i].subFragment[j].start_pos)
			continue;							// no subfragment is located at pos

		// 2. update sum_D : subFragment[j] is located at pos
		double D_hf = DFrag[i].D_hf, D_hfstar = DFrag[i].D_hfstar;  // previous DFrag values

		if (D_hf < D_hfstar) sum_D -= D_hf;		// substract the previous DFrag value
		else sum_D -= D_hfstar;

		int k = pos - Fragment[i].subFragment[j].start_pos;
		double q_j = Fragment[i].subFragment[j].qStr[k];
		double q_j_star = Fragment[i].subFragment[j].qStr_star[k];

		// computing under the assumption that the the bit is toggled. So, != -> ==	
		if (haplo.seq[pos - Fragment[block.start_frag].start_pos]
			== Fragment[i].subFragment[j].str[k]) {
			D_hf += q_j_star - q_j;			// + new value - previous value
			D_hfstar += q_j - q_j_star;
		}
		else {
			D_hf += q_j - q_j_star;			// + new value - previous value
			D_hfstar += q_j_star - q_j;
		}

		if (D_hf < D_hfstar)		// select min( D_h, D_h*)
			sum_D += D_hf;			// add the new DFrag value
		else
			sum_D += D_hfstar;

		if ( update ) {				// *** if update is true  ***
			DFrag[i].D_hf = D_hf;			// *** the calculated values are stored in DFrag ***
			DFrag[i].D_hfstar = D_hfstar;
		}

	}// end_for i (each fragment)

	if (update) 					// *** if update is true  ***
		haplo.sum_D = sum_D;

	return sum_D;
}

//
// calculate fitness values (2seq version)
//		using two sequences (normal and complement)
// This function is used only in find_homoSNP
//

double calc_sumD_2seq(HaploType& haplo, BlockType& block) {

	int startFrag, endFrag;
	double sum_D = 0.0;
	int pos_k;

	startFrag = block.start_frag;
	endFrag = block.end_frag;

	int realStartIndex = Fragment[block.start_frag].start_pos;

	for (int i = startFrag; i <= endFrag; i++) {	// for each fragment

		double D_hf = 0.0, D_hfstar = 0.0;

		for (int j = 0; j<Fragment[i].num_subFrag; j++) {		// for each subfragment
			for (int k = 0; k<Fragment[i].subFragment[j].length; k++) {		// for each position

				pos_k = k + Fragment[i].subFragment[j].start_pos - realStartIndex;

				double q_j = Fragment[i].subFragment[j].qStr[k];
				double q_j_star = Fragment[i].subFragment[j].qStr_star[k];

				if (haplo.seq[pos_k] != Fragment[i].subFragment[j].str[k]) 
					D_hf += q_j_star;
				else 
					D_hf += q_j;

				if (haplo.seq2[pos_k] != Fragment[i].subFragment[j].str[k]) 
					D_hfstar += q_j_star;
				else 
					D_hfstar += q_j;

			}  // end_for k (each position)

		}  // end_for j (each subfragment)

		if (D_hf < D_hfstar)
			sum_D += D_hf;
		else
			sum_D += D_hfstar;

	}// end_for i (each fragment)

	haplo.sum_D = sum_D;

	return sum_D;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// 
// finding position uncovered by matrix (Covered[])
// return : the number of uncovered position by matrix
//
void find_covered_sites(BlockType &block) {

	int startFrag, endFrag;
	int pos;

	for (int i = 0; i < block.length; i++)
		Covered[i] = 0;

	startFrag = block.start_frag;
	endFrag = block.end_frag;
	int realStartIndex = Fragment[block.start_frag].start_pos;

	// finding uncovered positions
	for (int i = startFrag; i <= endFrag; i++) 	// for each fragment
		for (int j = 0; j < Fragment[i].num_subFrag; j++) // for each subfragment
			for (int k = 0; k < Fragment[i].subFragment[j].length; k++) {	// for each position
				pos = k + Fragment[i].subFragment[j].start_pos - realStartIndex;
				Covered[pos] = 1;
			} // end_for k (each position)
}

//
// Genetic algorithm (EDA)
//
void GA(HaploType &besthaplo, BlockType &block) {

	bool isStop = true;
	double bestSum = -10;
	int stopCnt = 0;

	// uniform distribution function
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);
	uniform_real_distribution<double> distribution(0.0, 1.0);

	for (int i = 0; i < block.length; ++i)
		GAcnt[i] = GAcnt2[i] = 0;  	// array for storing the frequency of 1's for each position

	for (int i = 0; i<POPSIZE; i++) {	// initialize population randomly
		Population[i].sum_D = 0;
		Population[i].seq[block.length] = '\0';		// .seq2 is not unsed in GA

		for (int j = 0; j<block.length; j++) {

			if (distribution(generator)<0.5)
				Population[i].seq[j] = '0';
			else {
				Population[i].seq[j] = '1';

				if (i<OFFSIZE)
					GAcnt[j]++;
				else
					GAcnt2[j]++;
			}
		}

		calc_sumD(Population[i], block, 0);	// compute population[i].sum_D
	} // end_for

	sort(Population, Population+POPSIZE, compare_sumD_val);	// sorting with fitness values

	int g_iter = 50;		// generation number of GA
	while (g_iter--) {

		double prob;
		for (int i = 0; i<block.length; i++)
			GAcnt[i] += GAcnt2[i];

		for (int i = OFFSIZE; i<POPSIZE - 1; i++) {		// replace lower individual with new chromosomes
			for (int j = 0; j<block.length; j++) {

				prob = (double)GAcnt[j] / POPSIZE;

				if (Population[i].seq[j] == '1')
					GAcnt2[j]--;

				if (distribution(generator)<prob) {		// according to the frequency of 1 at the position
					Population[i].seq[j] = '1';
					GAcnt2[j]++;
				}
				else
					Population[i].seq[j] = '0';
			}
		}

		for (int i = OFFSIZE; i<POPSIZE - 1; i++)		// calculate the fitness values of new individuals
			calc_sumD(Population[i], block, 0);
		sort(Population, Population+POPSIZE, compare_sumD_val);


		if (bestSum - Population[0].sum_D > EPSILON || bestSum < -1 ) {	// check convergence 
			bestSum = Population[0].sum_D;
			stopCnt = 0;
		}
		else
			stopCnt++;

		if (stopCnt >= 10)		// if population is converged, stop
			break;
	}

	besthaplo.sum_D = Population[0].sum_D;
	strcpy(besthaplo.seq, Population[0].seq);		// storing best individual
	besthaplo.seq2[0] = '\0';						// empty sequence
}

// 
// Toggling for range switch
// return 1 if better soultion is found
//
int range_switch_toggling(HaploType & haplo, BlockType &block) {

	int rval=0, bestIdx = -1;
	bool isImp = true;
	double bestSum = haplo.sum_D;
	double tmpSum;

	calc_sumD(haplo, block, 1);	// update DFrag array for the current haplotype sequence

	while (isImp) {		// while finding better solution

		isImp = false;

		for (int i = 0; i<block.length; i++) {	// for each position in block

			// calculating the fitness value assuming that haplo.seq[i] is toggled
			tmpSum = calc_sumD_tog(haplo, block, i, 1); // DFrag array & haplo are updated

			if (bestSum - tmpSum > EPSILON) {	// if finding better solution
				isImp = true;
				bestSum = tmpSum;
				bestIdx = i;
			}
		} // end_for

		if (isImp) {	// update haplotype and DFrag array for the next iteration,.
			for (int i = 0; i <= bestIdx; i++)
				haplo.seq[i] = (haplo.seq[i] == '0') ? '1' : '0';
			calc_sumD(haplo, block, 1);
			rval = 1;
		}
	} // end_while

	return rval;
}

// 
// Toggling for single switch
// return 1 if better soultion is found
//
int single_switch_toggling(HaploType & haplo, BlockType &block) {

	int rval = 0, bestIdx = -1;
	bool isImp = true;
	double bestSum = haplo.sum_D;
	double tmpSum;

	calc_sumD(haplo, block, 1);	// update DFrag array for the current haplotype sequence

	while (isImp) {		// while finding better solution

		isImp = false;

		for (int i = 0; i<block.length; i++) {	// for each position in block

			// calculating the fitness value assuming that haplo.seq[i] is toggled
			tmpSum = calc_sumD_tog(haplo, block, i, 0);	// DFrag array & haplo are NOT updated

			if (bestSum - tmpSum > EPSILON) {
				isImp = true;
				bestSum = tmpSum;
				bestIdx = i;
			}
		} // end_for

		if (isImp) {	// update haplotype and DFrag array for the next iteration,.
			haplo.seq[bestIdx] = (haplo.seq[bestIdx] == '0') ? '1' : '0';
			calc_sumD(haplo, block, 1);
			rval = 1;
		}
	} // end_while

	return rval;
}

// 
// Creating complement sequence
// return : the number of phased position
//
int create_complement(HaploType & haplo, BlockType &block) {

	int phased = 0;

	for (int i = 0; i < block.length; i++)
		if (Covered[i]) {
			haplo.seq2[i] = (haplo.seq[i] == '0') ? '1' : '0';
			++phased;
		}
		else
			haplo.seq[i] = haplo.seq2[i] = '-';

	haplo.seq2[block.length] = '\0';

	return phased;
}

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
// test
//
// Toggling for range switch (toggling all [is..ie]'s)
// return 1 if better soultion is found
//
int range_switch_toggling2(HaploType & haplo, BlockType &block) {

	int rval = 0, bestsIdx = -1, besteIdx = -1;
	bool isImp = true;
	double bestSum = haplo.sum_D;
	double tmpSum;

	calc_sumD(haplo, block, 1);	// update DFrag array for the current haplotype sequence

	while (isImp) {		// while finding better solution

		isImp = false;

		for (int is = 0; is<block.length; is++) {	// for each position in block

			calc_sumD(haplo, block, 1);

			for (int ie = is; ie < block.length; ie++) {	// for each position in block

															// calculating the fitness value assuming that haplo.seq[ie] is toggled
				tmpSum = calc_sumD_tog(haplo, block, ie, 1); // DFrag array & haplo are updated

				if (bestSum - tmpSum > EPSILON) {	// if finding better solution
					isImp = true;
					bestSum = tmpSum;
					bestsIdx = is;
					besteIdx = ie;
				}
			} // end_for ie

		} // end_for is

		if (isImp) {	// update haplotype and DFrag array for the next iteration,.
			for (int i = bestsIdx; i <= besteIdx; i++)
				haplo.seq[i] = (haplo.seq[i] == '0') ? '1' : '0';
			calc_sumD(haplo, block, 1);
			rval = 1;
		}
	} // end_while

	return rval;
}


//
// finding homoSNP sites
//
void find_homoSNP(HaploType & haplo, BlockType &block) {

	double prevSum, tmpSum;

	//calc_sumD(haplo, block, 1);	// update DFrag array for the current haplotype sequence

	//range_switch_toggling2(haplo, block);	// Toggling for range switch [1..i]

	//create_complement(haplo, block);

	for (int i = 0; i<block.length; i++) {	// for each position in block

		if (haplo.seq[i] == '-') continue;

		haplo.seq[i] = (haplo.seq[i] == '0') ? '1' : '0';
		prevSum = haplo.sum_D;
		tmpSum = calc_sumD_2seq(haplo, block);

		if (tmpSum > prevSum) {	// rollback
			haplo.seq[i] = (haplo.seq[i] == '0') ? '1' : '0';
			haplo.sum_D = prevSum;
		}

		haplo.seq2[i] = (haplo.seq2[i] == '0') ? '1' : '0';
		prevSum = haplo.sum_D;
		tmpSum = calc_sumD_2seq(haplo, block);

		if (tmpSum > prevSum) {	// rollback
			haplo.seq2[i] = (haplo.seq2[i] == '0') ? '1' : '0';
			haplo.sum_D = prevSum;
		}
	} // end_for

}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

//
// haplotype phasing procedure
// return number of phased positions
//
int haplotype_phasing(BlockType &block, int phasing_iter) {

	find_covered_sites(block);

	BestHaplo.sum_D = -10;	// init for while loop

	while (phasing_iter--) {

		for (int j = 0; j < block.length; j++)	// init with "000.....000"
			TmpHaplo.seq[j] = '0';
		TmpHaplo.seq[block.length] = '\0';

		GA(TmpHaplo, block);							// Step 1: Genetic algorithm (EDA)

		int tog_iter = 10;			// toggling iteration number
		while (tog_iter--) {
			int imp = range_switch_toggling(TmpHaplo, block);	// Toggling for range switch [1..i]
//			int imp = range_switch_toggling2(TmpHaplo, block);	// Toggling for range switch [1..i]
			imp += single_switch_toggling(TmpHaplo, block);	// Toggling for single switch [1..i]
			if (!imp) break;				// if not improved, stop
		}

		if (BestHaplo.sum_D > TmpHaplo.sum_D || BestHaplo.sum_D < -1 ) {		// update best solution
			BestHaplo.sum_D = TmpHaplo.sum_D;
			strcpy(BestHaplo.seq, TmpHaplo.seq);
		}

		if (calc_MEC(BestHaplo, block) == 0)		// stop loop if mec value is 0
			break;
	}

	int phasedpos = create_complement(BestHaplo, block);

	/////////////////////////////////////////////////////////
	// if you want to find a sequence with homoSNP sites
	// comment out the following procedure;
	/////////////////////////////////////////////////////////

//	find_homoSNP(BestHaplo, block);


	return phasedpos;
}


///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

//
// load input matrix file
//
void load_matrixFile(char matrixFileName[])
{
	FILE * matrixFile = fopen(matrixFileName, "r");  	// open input file

	if (!matrixFile) {
		cout << "Inputfile \"" << matrixFileName << "\" does not exist." << endl;
		exit(1);
	}

	fseek(matrixFile, 0, SEEK_END);
	int fileSize = ftell(matrixFile);
	fseek(matrixFile, 0, SEEK_SET);

	// 1. load naive matrix data

	MatrixData = new char[fileSize + 10];
	int fragDataSize;;

	fragDataSize = fread(MatrixData, sizeof(char), fileSize + 10, matrixFile);	// read raw data in input file
	MatrixData[fragDataSize] = '\0';

	fclose(matrixFile);


	// 2. get # of fragments and # of positions (columns in input matrix)

	const char *token = " \t\n";	// token used in strtok
	char * p;

	p = strtok(MatrixData, token);
	NumFrag = atoi(p);				// number of fragments

	p = strtok(NULL, token);
	NumPos = atoi(p);				// number of positions

	// 3. build structure for matrix data

	Fragment = new FragType[NumFrag];

	for (int i = 0; i < NumFrag; i++) {	// storing fragment data in data structures

		if (!(p = strtok(NULL, token))) {
			NumFrag = i;			// set NumFrag to real number of frags
			break;
		}

		int num_subFrag = atoi(p); // number of subfragments in fragment[i]

		Fragment[i].num_subFrag = num_subFrag;
		Fragment[i].subFragment = new SubFragType[num_subFrag]; // array for storing subfragments

		p = strtok(NULL, token);		// skip fragment name

		// storing subfragment metadata
		for (int j = 0; j < num_subFrag; j++) {

			p = strtok(NULL, token);
			Fragment[i].subFragment[j].start_pos = atoi(p) - 1;

			p = strtok(NULL, token);
			Fragment[i].subFragment[j].str = p;
			Fragment[i].subFragment[j].length = strlen(Fragment[i].subFragment[j].str);
			Fragment[i].subFragment[j].end_pos = Fragment[i].subFragment[j].start_pos + Fragment[i].subFragment[j].length - 1;
		} // end_for j

		  // storing fragment metadata 
		Fragment[i].start_pos = Fragment[i].subFragment[0].start_pos;
		Fragment[i].end_pos = Fragment[i].subFragment[num_subFrag - 1].end_pos;
		Fragment[i].length = Fragment[i].end_pos - Fragment[i].start_pos + 1;

		p = strtok(NULL, "\n");

		// if fragment contains only one site, exclude the fragment in phasing
		if (Fragment[i].length == 1) {
			delete[] Fragment[i--].subFragment;
			continue;
		}

		// storing quality scores
		double * qStr = new double[Fragment[i].length];
		double * qStr_star = new double[Fragment[i].length];

		for (int j = 0; j < num_subFrag; j++)
		{
			// precomputing error probability
			for (int k = 0; k < Fragment[i].subFragment[j].length; k++) {
				int qual_j = (int)(*(p++));
				qual_j -= 33;
				qStr[k] = pow(10.0, ((-1 * (double)qual_j) / 10.0));
				qStr_star[k] = 1 - qStr[k];
			}

			Fragment[i].subFragment[j].qStr = qStr;
			Fragment[i].subFragment[j].qStr_star = qStr_star;

			qStr += (Fragment[i].subFragment[j].length);		// for the next subfragment
			qStr_star += (Fragment[i].subFragment[j].length);	// for the next subfragment

		} // end_for j

	} // end i

	sort(Fragment, Fragment + NumFrag, compare_frag_pos);	// sorting with starting positions

#ifdef __DEBUG__
															// recording fragment data in file for checking if fragment data is well loaded
	ofstream fragmentOutFile;
	fragmentOutFile.open("ZfragmentOut.txt");

	for (int i = 0; i < NumFrag; i++) {
		fragmentOutFile << "start: " << Fragment[i].start_pos
			<< " end: " << Fragment[i].end_pos
			<< " length: " << Fragment[i].length
			<< " num_subFrag: " << Fragment[i].num_subFrag << endl;

		for (int j = 0; j < Fragment[i].num_subFrag; j++) {
			fragmentOutFile << "start: " << Fragment[i].subFragment[j].start_pos
				<< " end: " << Fragment[i].subFragment[j].end_pos
				<< " length: " << Fragment[i].subFragment[j].length
				<< " str:" << Fragment[i].subFragment[j].str
				<< " qual:" << Fragment[i].subFragment[j].qStr;
		}
		fragmentOutFile << "==============================================================" << endl;
	}
	fragmentOutFile.close();
#endif
}

//
// build auxilary data structures
// 1. fragments info. for each position
// 2. blocks of fragments
//
void build_aux_struct(){

	// 1. build an array of fragment numbers locating at each position

	FragsForPos = new FragsForPosType[NumPos];
	for (int i = 0; i<NumPos; i++)
		FragsForPos[i] = { -1,-1 };

	int cutStart = 0;
	for (int i = 0; i<NumFrag; i++) 		// storing starting fragment #
		for (; cutStart <= Fragment[i].end_pos; cutStart++)
			if (FragsForPos[cutStart].start_frag == -1)
				FragsForPos[cutStart].start_frag = i;

	int cutEnd = NumPos - 1;
	for (int i = NumFrag - 1; i >= 0; i--) 	// storing ending fragment #
		for (; cutEnd >= Fragment[i].start_pos; cutEnd--)
			if (FragsForPos[cutEnd].end_frag == -1)
				FragsForPos[cutEnd].end_frag = i;

#ifdef __DEBUG__
	// recording range positions in file
	ofstream rangePosFile;
	rangePosFile.open("ZrangePos.txt");
	for (int i = 0; i<NumPos; i++) {
		rangePosFile << "[" << i << "] " << FragsForPos[i].start_frag << " " << FragsForPos[i].end_frag << endl;
	}
	rangePosFile.close();
#endif


	// 2. divide fragments into several blocks (sets of fragments whose covered positions were overlapped)

	int cutStartFrag = 0;	// set starting fragment
	int maxEndIdx = Fragment[0].end_pos;
	int block_len;

	MaxBlkLen = 0;

	for (int i = 1; i<NumFrag; i++) {
		if (Fragment[i].start_pos > maxEndIdx) {	// if frag[i] is NOT overlapped with the previous frags
			block_len = maxEndIdx - Fragment[cutStartFrag].start_pos + 1;
			Block.push_back({ cutStartFrag,i - 1, block_len }); // stroing current block info
			cutStartFrag = i;						// set the new starting fragment

			if (MaxBlkLen < block_len)	MaxBlkLen = block_len;
		}
		maxEndIdx = max(maxEndIdx, Fragment[i].end_pos); // update end position
	}

	block_len = maxEndIdx - Fragment[cutStartFrag].start_pos + 1;
	Block.push_back({ cutStartFrag, NumFrag - 1, block_len }); // for the last block
	if (MaxBlkLen < block_len)	MaxBlkLen = block_len;

#ifdef __DEBUG__
	// recording fragment numbers in each block
	ofstream blockFile;
	blockFile.open("ZBlockPos.txt");
	for (unsigned int i = 0; i<Block.size(); i++) {
		blockFile << "[" << i + 1 << "] " << Block[i].start_frag << " " << Block[i].end_frag << " " << Block[i].length << endl;

	}
	blockFile.close();
#endif

#ifdef __DEBUG__
	// recording fragment data in matrix form
	ofstream fragmentForm("ZFragmentForm.SORTED");

	for (unsigned int i = 0; i<Block.size(); i++) {

		fragmentForm << "block # : " << i + 1 << endl;

		for (int j = Block[i].start_frag; j <= Block[i].end_frag; j++) {

			int startIdx = Fragment[Block[i].start_frag].start_pos;

			char * frag = new char[(Block[i].length + 1) * 2];
			for (int k = 0; k<Block[i].length; k++)
				frag[2 * k] = ' ', frag[2 * k + 1] = '-';
			frag[Block[i].length * 2] = '\0';

			for (int k = 0; k<Fragment[j].num_subFrag; k++) {
				int idx = Fragment[j].subFragment[k].start_pos - startIdx;
				for (int l = 0; l<Fragment[j].subFragment[k].length; l++) {
					frag[2 * (idx + l) + 1] = Fragment[j].subFragment[k].str[l];
				}
			}
			fragmentForm << frag << endl;

			delete[] frag;
		}
	}

	fragmentForm.close();
#endif
	
}

void allocate() {

	DFrag = new DFragType[NumFrag];	// array for storing weighted MEC for each frament

	Population = new HaploType[POPSIZE];		// Population array

	for (int i = 0; i<POPSIZE; i++)
		Population[i].seq = new char[MaxBlkLen + 1];

	GAcnt = new int[MaxBlkLen + 1];		// array for storing the frequency of 1's for each position
	GAcnt2 = new int[MaxBlkLen + 1];

	BestHaplo.seq = new char[MaxBlkLen + 1];
	BestHaplo.seq2 = new char[MaxBlkLen + 1];
	TmpHaplo.seq = new char[MaxBlkLen + 1];
	TmpHaplo.seq2 = new char[MaxBlkLen + 1];

	Covered = new char[MaxBlkLen + 1];
}

void deallocate() {

	for (int i = 0; i < NumFrag; ++i) {
		delete[] Fragment[i].subFragment[0].qStr;	// qStr's is allocated with one array
		delete[] Fragment[i].subFragment[0].qStr_star;
		delete[] Fragment[i].subFragment;
	}
	delete[] Fragment, MatrixData, FragsForPos, DFrag;

	for (int i = 0; i<POPSIZE; i++)
		delete[] Population[i].seq;
	delete[] Population;

	delete[] GAcnt, GAcnt2;

	delete[] BestHaplo.seq, BestHaplo.seq2, TmpHaplo.seq, TmpHaplo.seq2;
	delete[] Covered;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

//
// main procedure
//
void procedure(char matrixFileName[], char outputFileName[], int phasing_iter) {

	std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

	load_matrixFile(matrixFileName); // load matrix data in input file into data structures

	build_aux_struct( );		// build auxilary data structures

	allocate();		// memory allocation

	int totalBlockSize = 0;
	double totalWMEC = 0.0;		// weighted mec value
	int totalMEC = 0;			// mec value
	int totalPhased = 0;		// total number of phased positions
	int totalReads = 0;			// total number of reads used for phasing
	int mec, phased;

	ofstream phasingResultFile;   // results file
	phasingResultFile.open(outputFileName);
	//	phasingResultFile.precision(5);

	for (unsigned int i = 0; i< Block.size() ; i++) {	// for each block, performing phasing

		//std::cout << ".";
		phased = haplotype_phasing(Block[i], phasing_iter);		// haplotyp phasing procedure

		mec = calc_MEC(BestHaplo, Block[i]);

		totalBlockSize += Block[i].length;
		totalPhased += phased;
		totalReads += Block[i].end_frag - Block[i].start_frag + 1;
		totalWMEC += BestHaplo.sum_D;
		totalMEC += mec;

		// ---------------------
		// 1. Printing Block Header
		// ---------------------
		phasingResultFile << "Block Number: " << i + 1 ;
		phasingResultFile << "  Block Length: " << Block[i].length;
		phasingResultFile << "  Phased Length: " << phased;
		phasingResultFile << "  Number of Reads: " << Block[i].end_frag - Block[i].start_frag + 1;
		phasingResultFile << "  Start position: " << Fragment[Block[i].start_frag].start_pos + 1;
		phasingResultFile << "  Weighted MEC: " << BestHaplo.sum_D;
		phasingResultFile << "  MEC: " << mec << endl;

		// ---------------------
		// 2. Printg Haloptype (Format 1) : printing each phased site in each line
		// ---------------------
		int offset = Fragment[Block[i].start_frag].start_pos + 1;	// 1st index is 1 
		for (int k = 0; k < Block[i].length; ++k)
			if(BestHaplo.seq[k] != '-')
				phasingResultFile << offset+k << "\t " << BestHaplo.seq[k]
					<< "\t " << BestHaplo.seq2[k] << endl;

		// ---------------------
		// 2. Printg Haloptype (Format 2) : printing phased sequence in one line (without blank)
		// ---------------------
		//phasingResultFile << BestHaplo.seq1 << endl;
		//phasingResultFile << BestHaplo.seq2 << endl << endl;

		// ---------------------
		// 2. Printg Haloptype (Format 3) : printing phased sequence in one line (with blank)
		// ---------------------
		//for (int k = 0; k < Block[i].length; ++k)
		//	phasingResultFile << " " << BestHaplo.seq[k];
		//phasingResultFile << endl;
		//for (int k = 0; k < Block[i].length; ++k)
		//	phasingResultFile << " " << BestHaplo.seq2[k];
		//phasingResultFile << endl << endl;

		// ---------------------
		// 3. Printg Block Closer
		// ---------------------
		phasingResultFile << "********" << endl;
		
	}

	deallocate();

	std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
	std::chrono::duration<double> DefaultSec = end - start;

	// ==============================================
	// Printing Total Results
	// ==============================================
	//phasingResultFile << "-------------------------------------" << endl;
	//phasingResultFile << "Total Block Length : " << totalBlockSize << endl;
	//phasingResultFile << "Total Phased Length : " << totalPhased << endl;
	//phasingResultFile << "Total Number of used Reads : " << totalReads << endl;
	//phasingResultFile << "Total Weighted MEC : " << totalWMEC << endl;
	//phasingResultFile << "Total MEC : " << totalMEC << endl;
	//phasingResultFile << "Total Time : " << DefaultSec.count() << " seconds" << std::endl;
	//phasingResultFile << "-------------------------------------" << endl;
	//----------------------------------------------------------------
	phasingResultFile.close();

	//--------------------
	// Printing in stdout
	//--------------------
	// std::cout << "Block_len: "<<totalBlockSize << "   Phased_len: " << totalPhased
	//	<< "   Used_Reads: " << totalReads 
	//	<< "   w_MEC: " << totalWMEC << "   MEC: " << totalMEC
	//	<< "   Time: " << DefaultSec.count() << endl;
	std::cout << totalBlockSize << "   " << totalPhased << "   " << totalReads << "  "
		<< totalWMEC << "   " << totalMEC << "   " << DefaultSec.count() << endl;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

//
// main function
// 1st argument : input file name
// 2nd argument : output filne name

int main(int argc, char ** argv) {

	if (argc != 3 && argc != 4) {
		cout << "usage: " << argv[0] << " <input_file> <output_file> <param>" << endl;
		cout << "<param> is a positive integer (optional, default : 50)" << endl;
		exit(1);
	}

	char matrixFileName[100], outputFileName[100];

	strcpy(matrixFileName, argv[1]);
	strcpy(outputFileName, argv[2]);

	int phasing_iter = 50;				// default : 50

	if (argc == 4) {
		phasing_iter = atoi(argv[3]);
		if (phasing_iter < 1) {
			cout << "usage: " << argv[0] << " <input_file> <output_file> <param>" << endl;
			cout << "<param> is a positive integer (optional, default : 50)" << endl;
			exit(1);
		}
	}

	procedure(matrixFileName, outputFileName, phasing_iter);		//// main procedure

	return 0;
}

