#include <iostream>
#include <fstream>
#include <algorithm>
#include<functional>
#include <string.h>

using namespace std;

#define MAX_PNUM 200000
#define MAX_BLK 10000

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

typedef struct {			// block
	int start_pos;
	int end_pos;
	int phased;
	char *seq1, *seq2;
}BlockType2;

///////////////////////////////////////////////
///////////////////////////////////////////////

static int NumBLK;
static int NumPos;
static int NumFrag;

///////////////////////////////////////////////
///////////////////////////////////////////////

// using two haplotype sequences
int calc_MEC(BlockType2 block, FragType fragment[]) {

	int sum_D = 0, pos_k;

	char *seq1 = block.seq1;
	char *seq2 = block.seq2;

	int pos_s = block.start_pos;
	int pos_e = block.end_pos;

	int i;

	for ( i = 0; i < strlen(seq1); ++i)			// for sdhap evaluation
		if (seq1[i] == '0' || seq2[i] == '0')	break;
	if (i == strlen(seq1))
		return 0;

	for ( i = 0; i < NumFrag; i++) {	// for each fragment

		if (fragment[i].start_pos > pos_e || fragment[i].end_pos < pos_s)  continue;	// skip

		int D_hf = 0, D_hstarf = 0;

		for (int j = 0; j<fragment[i].num_subFrag; j++) {		// for each subfragment

			if (fragment[i].subFragment[j].start_pos > pos_e || fragment[i].subFragment[j].end_pos < pos_s)  continue;	// skip

			for (int k = 0; k<fragment[i].subFragment[j].length; k++) {		// for each position

				pos_k = fragment[i].subFragment[j].start_pos + k;		// convert to real positions

				if (pos_k < block.start_pos || pos_k > block.end_pos)		// outside of block
					continue;

				pos_k -= pos_s;					// position in block 

				if (fragment[i].subFragment[j].str[k] == '-')
					continue;
				if (seq1[pos_k] == '-' && seq2[pos_k] == '-')
					continue;

				if (seq1[pos_k] != '-' && seq1[pos_k] != fragment[i].subFragment[j].str[k])
					++D_hf;
				if (seq2[pos_k] != '-' && seq2[pos_k] != fragment[i].subFragment[j].str[k])
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


void build_pattern(char sw[], char ansseq1[], char ansseq2[], BlockType2 block, ofstream &LogFile){

	int k, i = 0;
	char *hapseq1 = block.seq1;
	char *hapseq2 = block.seq2;
	int pos_s = block.start_pos;

	LogFile << "seq1: " << hapseq1 << endl;
	LogFile << "seq2: " << hapseq2 << endl;
	LogFile << "      " ;
	for (k = block.start_pos; k <= block.end_pos; ++k) {	// conversion
		if (hapseq1[k-pos_s] == '-') {
			LogFile << '.';
			continue;
		}
		if (ansseq1[k] == '-') {
			LogFile << '+';
			continue;
		}

		if (ansseq1[k] == hapseq1[k-pos_s]) LogFile << '0';
		else LogFile << '1';

		if (ansseq1[k] == hapseq1[k-pos_s]) sw[i] = '0';
		else sw[i] = '1';

		i++;
	}
	sw[i] = '\0';
	LogFile << endl;

	LogFile << "ans1: ";
	for (k = block.start_pos; k <= block.end_pos; ++k)
		LogFile << ansseq1[k];
	LogFile << endl;
	LogFile << "ans2: ";
	for (k = block.start_pos; k <= block.end_pos; ++k)
		LogFile << ansseq2[k];
	LogFile << endl << endl;

}

// from Probhap code
void correct_pattern(char sw[], ofstream &LogFile) {

	int len = strlen(sw);

	if (len < 2) return;

	if (sw[0] != sw[1] )	sw[0] = sw[1];		// 01... -> 11... , 10... -> 00...
	if (sw[len - 2] != sw[len - 1])		sw[len - 1] = sw[len - 2];  // ...01 -> ...00 , ...10 -> ...11 

	for (int i = 0; i < len; ++i)
		if (!strncmp(sw + i, "01010", 5))
			sw[i + 1] = sw[i + 3] = '0';
		else if (!strncmp(sw + i, "10101", 5))
			sw[i + 1] = sw[i + 3] = '1';

	for (int i = 0; i < len; ++i)
		if (!strncmp(sw + i, "010", 3))
			sw[i + 1] = '0';
		else if (!strncmp(sw + i, "101", 3))
			sw[i + 1] = '1';
}

int count_swer(char sw[]) {

	int len = strlen(sw);
	int swer = 0;

	int cnt = 1;
	for (int j = 1; j < len; ++j) 
		if (sw[j - 1] == sw[j])	cnt++;
		else ++swer, cnt = 1;

	return swer;
}

// compare function for sorting
bool compare_phased(BlockType2 left, BlockType2 right) {
	return left.phased > right.phased;		// reverse
}

int N50(BlockType2 block[], int tot_len) {
	int i, sum;

	sort(block, block+NumBLK, compare_phased);

	tot_len /= 2;
	sum = i = 0;
	while (sum < tot_len)
		sum += block[i++].phased;

	return block[i - 1].phased;
}

///////////////////////////////////////////////
///////////////////////////////////////////////

void evaluate_haplo(char *EvalOutFileName, char ansseq1[], char ansseq2[], BlockType2 block[], FragType frag[]){

	int tot_blk_len, tot_phased, tot_sw_len, tot_Rswer, tot_Sswer, tot_mec;

	ofstream LogFile;
	LogFile.open(EvalOutFileName);

	tot_blk_len = tot_phased = tot_sw_len = tot_Rswer = tot_Sswer = tot_mec = 0;

	for (int i = 0; i < NumBLK; ++i) {

		LogFile << "================================================" << endl;
		LogFile << "Block: " << i + 1 << "  ( "<< block[i].start_pos+1 <<" , " << block[i].end_pos+1 << " )" <<endl;

		int phased = block[i].phased;
		int blk_len = block[i].end_pos - block[i].start_pos + 1;
		char *sw = new char[blk_len + 10];

		build_pattern(sw, ansseq1, ansseq2, block[i], LogFile);
		LogFile << "sw: " << sw << endl;

		if (phased < 2) continue;	// skip for block of length 1

		int swer = count_swer(sw);

		correct_pattern(sw, LogFile);
		LogFile << "sw: " << sw << endl<<endl;

		int Rswer = count_swer(sw);

		int sw_len = strlen(sw);
		if (sw_len < 1) sw_len = 1;

		delete[] sw;

		int mec = calc_MEC(block[i], frag);

		tot_blk_len += blk_len;
		tot_phased += phased;
		tot_sw_len += (sw_len - 1);
		tot_Rswer += Rswer;
		tot_Sswer += swer - Rswer;
		tot_mec += mec;

		LogFile << "block_len (tot): " << blk_len << "  (" << tot_blk_len << ")" << endl
			<< "phased_len (tot): " << phased << "  (" << tot_phased << ")" << endl
			<< "swcheck (tot): " << sw_len-1 << "  (" << tot_sw_len << ")" << endl
			<< "swer (tot): " << swer << "  ("<< tot_Rswer+tot_Sswer << ")"<< endl
			<< "rswer (tot): " << Rswer << "  (" << tot_Rswer << ")" << endl
			<< "sswer (tot): " << swer - Rswer << "  (" << tot_Sswer << ")" << endl
			<< "mec (tot): " << mec << "  (" << tot_mec << ")" << endl <<endl;
	} // end_for each block

	int n50 = N50(block, tot_phased);

	LogFile << "------------------------------------------------" << endl;
	LogFile << "chr_len: " << NumPos
		<< "  block_len: " << tot_blk_len << "  phased_len: " << tot_phased
		<< "  swcheck: " << tot_sw_len	<< "  swer:  " << tot_Sswer+tot_Rswer 
		<< "  rswer: " << tot_Rswer << "  sswer: " << tot_Sswer
		<< "  mec: " << tot_mec << "  N50: " << n50 << endl;

	//cout << "chr_len: " << NumPos
	//	<< "  block_len: " << tot_blk_len << "  phased_len: " << tot_phased
	//	<< "  swcheck: " << tot_sw_len << "  swer:  " << tot_Sswer + tot_Rswer
	//	<< "  rswer: " << tot_Rswer << "  sswer: " << tot_Sswer
	//	<< "  mec: " << tot_mec << "  N50: " << n50 << endl;

	cout << NumPos
		<< "  " << tot_blk_len << "  " << tot_phased
		<< "  " << tot_sw_len << "  " << tot_Sswer + tot_Rswer
		<< "  " << tot_Rswer << "  " << tot_Sswer
		<< "  " << tot_mec << "  " << n50 << endl;
	
	LogFile.close();
}

///////////////////////////////////////////////
///////////////////////////////////////////////

char *load_to_mem(char inFileName[] ) {

	FILE * inFile = fopen(inFileName, "r");  	// open input file

	if (!inFile) {
		cout << "Inputfile \"" << inFileName << "\" does not exist." << endl;
		exit(1);
	}

	fseek(inFile, 0, SEEK_END);
	int fileSize = ftell(inFile);
	fseek(inFile, 0, SEEK_SET);

	char *buf = new char[fileSize + 10];
	int buf_size = fread(buf, sizeof(char), fileSize + 10, inFile);
	buf[buf_size] = '\0';

	fclose(inFile);

	return buf;
}

void load_answer_vcf_file(char answerFileName[], char ansseq1[], char ansseq2[]) {

	const char *token = " \t\n";	// token used in strtok
	char *ptr;
	int k;

	for (k = 0; k < MAX_PNUM; ++k)		// initialization
		ansseq1[k] = ansseq2[k] = '-';

	char *buf = load_to_mem(answerFileName);

	ptr = strtok(buf, token);	// 1st column

	k = 0;			// position (1 - origin in file, 0 - origin in program)
					// line number is position in .vcf file
	while (ptr) {
		ptr = strtok(NULL, token);	// 2nd column
		ansseq1[k] = *ptr;
		ptr = strtok(NULL, token);	// 3rd column
		ansseq2[k] = *ptr;

		if (++k > MAX_PNUM - 10) {
			cout << "too small array size" << endl;
			exit(1);
		}

		ptr = strtok(NULL, token);	// 1st column
	}
	delete[] buf;
}

void load_answer_valid_file(char answerFileName[], char ansseq1[], char ansseq2[]) {

	const char *token = " \t\n";	// token used in strtok
	char *ptr;
	int k;

	for (k = 0; k < MAX_PNUM; ++k)		// initialization
		ansseq1[k] = ansseq2[k] = '-';

	char *buf = load_to_mem(answerFileName);

	ptr = strtok(buf, token);	// 1st column

	while (ptr) {
		k = atoi(ptr) -1 ;			// position (1-origin in file, 0-origin in program)
		ptr = strtok(NULL, token);	// skip 2nd column
		ptr = strtok(NULL, token);	// 3rd column
		ansseq1[k] = *ptr;
		ptr = strtok(NULL, token);	// 4th column
		ansseq2[k] = *ptr;

		if (k >= MAX_PNUM) {
			cout << "too small array size" << endl;
			exit(1);
		}

		ptr = strtok(NULL, token);	// 1st column
	}

	delete[] buf;
}

char *get_line(char buf[], int init) {
	static int offset;
	if (init) offset = 0;
	char *line = strtok(buf + offset, "\n");
	if (line)	offset += strlen(line) + 1;
	return line;
}
void load_haplo_file(char haploFileName[], BlockType2 block[], char maskseq[]) {

	const char *token = " \t";	// token used in strtok
	char *ptr;
	int k, blk_idx, pos_s, offset = 0, cnt, phased;

	char *hapseq1 = new char[MAX_PNUM];
	char *hapseq2 = new char[MAX_PNUM];

	for (k = 0; k < MAX_PNUM; ++k)		// initialization
		hapseq1[k] = hapseq2[k] = '-';

	char *buf = load_to_mem(haploFileName);

	blk_idx = 0;	// block index;

	ptr = get_line(buf, 1);	// read 1st line (block header)

	while (ptr) {					// line by line

		ptr = strtok(ptr, token);		// 1st column

		if (*ptr == '*') {		// end of block  "*******"

			block[blk_idx].end_pos = k;
			block[blk_idx].phased = phased;

			int blk_len = k - pos_s + 1;

			block[blk_idx].seq1 = new char[blk_len + 1];
			block[blk_idx].seq2 = new char[blk_len + 1];

			strncpy(block[blk_idx].seq1, hapseq1, blk_len);
			strncpy(block[blk_idx].seq2, hapseq2, blk_len);

			block[blk_idx].seq1[blk_len] = block[blk_idx].seq2[blk_len] = '\0';

			for (k = 0; k < blk_len; ++k)		// initialization
				hapseq1[k] = hapseq2[k] = '-';

			blk_idx++;
			if (blk_idx >= MAX_BLK) {
				cout << "too small number of blocks" << endl;
				exit(1);
			}
		}
		else if (isdigit(*ptr)) {		// if haplotype line

			k = atoi(ptr) - 1;			// position (1-origin in file, 0-origin in program)

			if (cnt++ == 0) 	// first element of block
				block[blk_idx].start_pos = pos_s = k;

			ptr = strtok(NULL, token);	// 2nd column
			if( !maskseq[k] )	hapseq1[k - pos_s] = *ptr;
			ptr = strtok(NULL, token);	// 3rd column
			if (!maskseq[k])	hapseq2[k - pos_s] = *ptr;

			if (hapseq1[k - pos_s] != '-' || hapseq2[k - pos_s] != '-')
				++phased;				// counting phased positions

			if (k >= MAX_PNUM - 10) {
				cout << "too small array size" << endl;
				exit(1);
			}
		}
		else		// if header line
			cnt = phased = 0;

		ptr = get_line(buf, 0);	// read 1st line (block header)
	}
	NumBLK = blk_idx;

	delete[] buf;
	delete[] hapseq1;
	delete[] hapseq2;
}

char *get_line_mask(char buf[], int init) {
	static int offset;
	if (init) offset = 0;

	char *line = strtok(buf + offset, "\n");
	if (line)	offset += strlen(line) + 1;
	return line;
}

void load_mask_haplo_file(char maskhaploFileName[], char maskseq[]) {

	const char *token = " \t";	// token used in strtok
	char *ptr;
	int k;

	for (k = 0; k < MAX_PNUM; ++k)		// initialization
		maskseq[k] = '-';

	char *buf = load_to_mem(maskhaploFileName);

	ptr = get_line_mask(buf, 1);	// read 1st line (block header)

	while (ptr) {					// line by line

		ptr = strtok(ptr, token);		// 1st column

		if (isdigit(*ptr)) {		// if haplotype line
			k = atoi(ptr) - 1;			// position (1-origin in file, 0-origin in program)
			if (k >= MAX_PNUM - 10) {
				cout << "too small array size" << endl;
				exit(1);
			}
			ptr = strtok(NULL, token);	// 2nd column
			if (*ptr != '-') {
				ptr = strtok(NULL, token);	// 3rd column
				if (*ptr != '-')
					maskseq[k] = 0;
			}
		}

		ptr = get_line_mask(buf, 0);	// read 1st line (block header)
	}

	delete[] buf;
}

///////////////////////////

// compare function for sorting
bool compare_frag_pos2(FragType left, FragType right) {
	return left.start_pos < right.start_pos;
}

//
// load input matrix file
//
void load_matrixFile2(char matrixFileName[], FragType *frag[], char *fragdata[]) {

	FILE * matrixFile = fopen(matrixFileName, "r");  	// open input file

	if (!matrixFile) {
		cout << "Inputfile \"" << matrixFileName << "\" does not exist." << endl;
		exit(1);
	}

	fseek(matrixFile, 0, SEEK_END);
	int fileSize = ftell(matrixFile);
	fseek(matrixFile, 0, SEEK_SET);

	// 1. load naive matrix data

	char *fragData = new char[fileSize + 10];
	int fragDataSize;;
	*fragdata = fragData;

	fragDataSize = fread(fragData, sizeof(char), fileSize + 10, matrixFile);	// read raw data in input file
	fragData[fragDataSize] = '\0';

	fclose(matrixFile);


	// 2. get # of fragments and # of positions (columns in input matrix)

	const char *token = " \t\n";	// token used in strtok
	char * p;

	p = strtok(fragData, token);
	NumFrag = atoi(p);				// number of fragments

	p = strtok(NULL, token);
	NumPos = atoi(p);				// number of positions

	// 3. build structure for matrix data

	FragType *fragment = new FragType[NumFrag];
	*frag = fragment;

	for (int i = 0; i < NumFrag; i++) {	// storing fragment data in data structures

		if (!(p = strtok(NULL, token))) {
			NumFrag = i;			// set NumFrag to real number of frags
			break;
		}

		int num_subFrag = atoi(p); // number of subfragments in fragment[i]

		fragment[i].num_subFrag = num_subFrag;
		fragment[i].subFragment = new SubFragType[num_subFrag]; // array for storing subfragments

		p = strtok(NULL, token);		// skip fragment name

										// storing subfragment metadata
		for (int j = 0; j < num_subFrag; j++) {

			p = strtok(NULL, token);
			fragment[i].subFragment[j].start_pos = atoi(p) - 1;

			p = strtok(NULL, token);
			fragment[i].subFragment[j].str = p;
			fragment[i].subFragment[j].length = strlen(fragment[i].subFragment[j].str);
			fragment[i].subFragment[j].end_pos = fragment[i].subFragment[j].start_pos + fragment[i].subFragment[j].length - 1;
		} // end_for j

		  // storing fragment metadata 
		fragment[i].start_pos = fragment[i].subFragment[0].start_pos;
		fragment[i].end_pos = fragment[i].subFragment[num_subFrag - 1].end_pos;
		fragment[i].length = fragment[i].end_pos - fragment[i].start_pos + 1;

		p = strtok(NULL, "\n");		// skip quality information

		// if fragment contains only one site, exclude the fragment in phasing
		if (fragment[i].length == 1) {
			delete[] fragment[i--].subFragment;
			continue;
		}
	} // end i

	sort(fragment, fragment + NumFrag, compare_frag_pos2);	// sorting with starting positions

}

void deallocate_fragment(FragType *fragment, char *fragData, BlockType2 block[]) {

	for (int i = 0; i < NumFrag; ++i) {
		delete[] fragment[i].subFragment;
	}
	delete[] fragment, fragData;

	for (int i = 0; i < NumBLK; ++i) {
		delete[] block[i].seq1;
		delete[] block[i].seq2;
	}
}

/////////////////////////

//
// main function
// 1st argument : matrix file
// 2nd argument : answer file
// 3rd argument : haplotype file
// 4th argument : output file (evaluation results)
// 5th argument : mask_file (optional)

int main(int argc, char ** argv) {

	if (argc < 5) {
		cout << "usage: " << argv[0] << " <matrix_file> <answer_file> <haplo_file> <output_file> <mask_file>" << endl;
		exit(1);
	}

	// matrix file loading (argv[1])
	FragType *fragm = NULL;
	char *fragdata = NULL;				// raw matrix data

	load_matrixFile2(argv[1], &fragm, &fragdata);

	// answer file loading (argv[2])
	char ansseq1[MAX_PNUM];
	char ansseq2[MAX_PNUM];

	char *f;
	char filename[100];

	strcpy(filename, argv[2]);
	strtok(filename, ".");
	f = strtok(NULL, ".");
	
	if( !strcmp (f, "vcf") )
		load_answer_vcf_file(argv[2], ansseq1, ansseq2 );
	else if (!strcmp(f, "valid"))
		load_answer_valid_file(argv[2], ansseq1, ansseq2);
	else {
		cout << "Invalid Answer filename" << endl;
		exit(1);
	}
	
	// mask file loading  (argv[5])
	char maskseq[MAX_PNUM] = { 0 };

	if( argc == 6)
		load_mask_haplo_file(argv[5], maskseq);	// MixSIH, PEATH, ProbHap, SDhaP
	
	// haplo file loading  (argv[3])
	BlockType2 block[MAX_BLK];

	load_haplo_file(argv[3], block, maskseq);	// MixSIH, PEATH, ProbHap, SDhaP

	// evaluation (argv[4])
	evaluate_haplo(argv[4], ansseq1, ansseq2, block, fragm);

	deallocate_fragment(fragm, fragdata, block);

	return 0;
}

