#include <omp.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include <sys/types.h> 
#include <sys/stat.h> 
#include <math.h>
#include <limits.h>
#include <bitset>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <ctype.h> 
#include <sstream>
#include <set>
#include <memory>
#include <typeinfo>
#include <pthread.h>
#include <vector>
#include <functional>
#include <iterator>
#include <string.h>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/seeds.h>
#include <seqan/score.h>
#include <seqan/modifier.h>
#include "seed.hpp"
#include "score.hpp"

#include "xcl2.hpp"

//HBM Banks requirements
#define MAX_HBM_BANKCOUNT 32
#define BANK_NAME(n) n | XCL_MEM_TOPOLOGY
const int bank[MAX_HBM_BANKCOUNT] = {
        BANK_NAME(0),  BANK_NAME(1),  BANK_NAME(2),  BANK_NAME(3),  BANK_NAME(4),
        BANK_NAME(5),  BANK_NAME(6),  BANK_NAME(7),  BANK_NAME(8),  BANK_NAME(9),
        BANK_NAME(10), BANK_NAME(11), BANK_NAME(12), BANK_NAME(13), BANK_NAME(14),
        BANK_NAME(15), BANK_NAME(16), BANK_NAME(17), BANK_NAME(18), BANK_NAME(19),
        BANK_NAME(20), BANK_NAME(21), BANK_NAME(22), BANK_NAME(23), BANK_NAME(24),
        BANK_NAME(25), BANK_NAME(26), BANK_NAME(27), BANK_NAME(28), BANK_NAME(29),
        BANK_NAME(30), BANK_NAME(31)};

#define NOW std::chrono::high_resolution_clock::now()

using namespace std;


#define MIN -32768
#define BYTES_INT 4
//#define N_STREAMS 60
//#define MAX_SIZE_ANTIDIAG 8000
#define MAX_DEVICES 1
#define NUM_KERNEL 2
//trying to see if the scoring scheme is a bottleneck in some way
#define MATCH     1
#define MISMATCH -1
#define GAP_EXT  -1
#define GAP_OPEN -1
#define UNDEF -32767
#define WARP_DIM 32 
#define EXTEND_NONEL  0
#define EXTEND_LEFTL  1
#define EXTEND_RIGHTL 2
#define EXTEND_BOTHL  3

typedef int ExtensionDirectionL;


using namespace chrono;


//=======================================================================
//
// Common functions
//
//=======================================================================

typedef seqan::Seed<seqan::Simple> TSeed;
typedef std::tuple< int, int, int, int, int, double > myinfo;	// score, start seedV, end seedV, start seedH, end seedH, runtime

char dummycomplement (char n)
{
	switch(n)
	{
	case 'A':
		return 'T';
	case 'T':
		return 'A';
	case 'G':
		return 'C';
	case 'C':
		return 'G';
	}
	assert(false);
	return ' ';
}

vector<std::string> split (const std::string &s, char delim)
{
	std::vector<std::string> result;
	std::stringstream ss (s);
	std::string item;

	while (std::getline (ss, item, delim))
	{
		result.push_back (item);
	}

	return result;
}

//=======================================================================
//
// SeqAn and FPGA function calls
//
//=======================================================================


inline void extendSeedFPGA(vector<SeedL> &seeds,
			ExtensionDirectionL direction,
			vector<string> &target,
			vector<string> &query,
			vector<ScoringSchemeL> &penalties,
			int const& XDrop,
			int const& kmer_length,
			int *res,
			int numAlignments,
			std::string binary_file
			)
{

	if(scoreGapExtend(penalties[0]) >= 0){

		cout<<"Error: Does not support gap extension penalty >= 0\n";
		exit(-1);
	}
	if(scoreGapOpen(penalties[0]) >= 0){

		cout<<"Error: Does not support gap opening penalty >= 0\n";
		exit(-1);
	}

	// NB nSequences is correlated to the number of FPGASs that we have
	
	int nSequences = numAlignments; //single device
	//int nSequencesLast = nSequences+numAlignments;

	//final result of the alignment
	std::vector<short, aligned_allocator<short>> scoreLeft(numAlignments);
	std::vector<short, aligned_allocator<short>> scoreRight(numAlignments);
	
	//create two sets of seeds
	//copy seeds
	std::vector<SeedL, aligned_allocator<SeedL>>  seeds_r;
	std::vector<SeedL, aligned_allocator<SeedL>>  seeds_l;
	seeds_r.reserve(numAlignments);
	//seeds_l.reserve(numAlignments);
	for (int i=0; i<seeds.size(); i++){
			seeds_r.push_back(seeds[i]);
			seeds_l.push_back(seeds[i]);
	}

	//sequences offsets	 		
	std::vector<int, aligned_allocator<int>> offsetLeftQ;
	std::vector<int, aligned_allocator<int>> offsetLeftT;	
	std::vector<int, aligned_allocator<int>> offsetRightQ;	
	std::vector<int, aligned_allocator<int>> offsetRightT;

	//shared_mem_size per block per FPGA
	int ant_len_left;
	int ant_len_right;

	//total lenght of the sequences
	int totalLengthQPref;
	int totalLengthTPref;
	int totalLengthQSuff;
	int totalLengthTSuff;
	
	int dim = nSequences;
	// if(i==nfpgas-1)
	// 	dim = nSequencesLast;
	//compute offsets and shared memory per block
	auto start_setup_ithread = NOW;
	ant_len_left=0;
	ant_len_right=0;
	for(int j = 0; j < dim; j++){

		offsetLeftQ.push_back(getBeginPositionV(seeds[j]));
		offsetLeftT.push_back(getBeginPositionH(seeds[j]));
		ant_len_left = std::max(std::min(offsetLeftQ[j],offsetLeftT[j]), ant_len_left);
		
		offsetRightQ.push_back(query[j].size()-getEndPositionV(seeds[j]));
		offsetRightT.push_back(target[j].size()-getEndPositionH(seeds[j]));
		ant_len_right = std::max(std::min(offsetRightQ[j], offsetRightT[j]), ant_len_right);
	}
	
	//compute antidiagonal offsets
	partial_sum(offsetLeftQ.begin(),offsetLeftQ.end(),offsetLeftQ.begin());	
	partial_sum(offsetLeftT.begin(),offsetLeftT.end(),offsetLeftT.begin());
	partial_sum(offsetRightQ.begin(),offsetRightQ.end(),offsetRightQ.begin());
	partial_sum(offsetRightT.begin(),offsetRightT.end(),offsetRightT.begin());
	//set total length of the sequences
	totalLengthQPref = offsetLeftQ[dim-1];
	totalLengthTPref = offsetLeftT[dim-1];
	totalLengthQSuff = offsetRightQ[dim-1];
	totalLengthTSuff = offsetRightT[dim-1];
	//allocate sequences prefix and suffix on the CPU
	//declare and allocate sequences prefixes and suffixes
	std::vector<char, aligned_allocator<char>> prefQ(totalLengthQPref);
	std::vector<char, aligned_allocator<char>> suffQ(totalLengthQSuff); 
	std::vector<char, aligned_allocator<char>> prefT(totalLengthTPref);
	std::vector<char, aligned_allocator<char>> suffT(totalLengthTSuff);
	// prefQ = (char*)malloc(sizeof(char)*totalLengthQPref);
	// prefT = (char*)malloc(sizeof(char)*totalLengthTPref);
	// suffQ = (char*)malloc(sizeof(char)*totalLengthQSuff);
	// suffT = (char*)malloc(sizeof(char)*totalLengthTSuff);
	//generate prefix and suffix on the CPU
	//std::cout << "SETTING UP PREF/SUFF" << std::endl;
	reverse_copy(query[0].c_str(),query[0].c_str()+offsetLeftQ[0],prefQ.data());
	//memcpy(prefQ[i], query[0+i*nSequences].c_str(), offsetLeftQ[i][0]);
	memcpy(prefT.data(), target[0].c_str(), offsetLeftT[0]);
	memcpy(suffQ.data(), query[0].c_str()+getEndPositionV(seeds[0]), offsetRightQ[0]);
	reverse_copy(target[0].c_str()+getEndPositionH(seeds[0]),target[0].c_str()+getEndPositionH(seeds[0])+offsetRightT[0],suffT.data());
	//memcpy(suffT[i], target[0+i*nSequences].c_str()+getEndPositionH(seeds[0+i*nSequences]), offsetRightT[i][0]);
	for(int j = 1; j<dim; j++){
		char *seqptr = prefQ.data() + offsetLeftQ[j-1];
		reverse_copy(query[j].c_str(),query[j].c_str()+(offsetLeftQ[j]-offsetLeftQ[j-1]),seqptr);
		//memcpy(seqptr, query[j+i*nSequences].c_str(), offsetLeftQ[i][j]-offsetLeftQ[i][j-1]);
		seqptr = prefT.data() + offsetLeftT[j-1];
		memcpy(seqptr, target[j].c_str(), offsetLeftT[j]-offsetLeftT[j-1]);
		seqptr = suffQ.data() + offsetRightQ[j-1];
		memcpy(seqptr, query[j].c_str()+getEndPositionV(seeds[j]), offsetRightQ[j]-offsetRightQ[j-1]);
		seqptr = suffT.data() + offsetRightT[j-1];
		reverse_copy(target[j].c_str()+getEndPositionH(seeds[j]),target[j].c_str()+getEndPositionH(seeds[j])+(offsetRightT[j]-offsetRightT[j-1]),seqptr);
		//memcpy(seqptr, target[j+i*nSequences].c_str()+getEndPositionH(seeds[j+i*nSequences]), offsetRightT[i][j]-offsetRightT[i][j-1]);

	}
	
	std::vector<short, aligned_allocator<short>> ant_l(sizeof(short)*ant_len_left*3*dim);
	std::vector<short, aligned_allocator<short>> ant_r(sizeof(short)*ant_len_right*3*dim);



	cl_int err;
    cl::CommandQueue q;
    std::string krnl_name = "xdrop";
	std::vector<cl::Kernel> krnls(NUM_KERNEL);
	cl::Context context;
	auto devices = xcl::get_xil_devices();

    // read_binary_file() command will find the OpenCL binary file created using the
    // V++ compiler load into OpenCL Binary and return pointer to file buffer.
    auto fileBuf = xcl::read_binary_file(binary_file);

    cl::Program::Binaries bins{{fileBuf.data(), fileBuf.size()}};
    int valid_device = 0;
    for (unsigned int i = 0; i < devices.size(); i++) {
        auto device = devices[i];
        // Creating Context and Command Queue for selected Device
        OCL_CHECK(err, context = cl::Context(device, NULL, NULL, NULL, &err));
        OCL_CHECK(err,
                q = cl::CommandQueue(context,
                                     device,
                                     CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE |
                                     CL_QUEUE_PROFILING_ENABLE,
                                     &err));

        // std::cout << "Trying to program device[" << i
                            // << "]: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
        cl::Program program(context, {device}, bins, NULL, &err);
        if (err != CL_SUCCESS) {
            // std::cout << "Failed to program device[" << i
                    // << "] with xclbin file!\n";
        } else {
            // std::cout << "Device[" << i << "]: program successful!\n";
            // Creating Kernel object using Compute unit names

            for (int i = 0; i < NUM_KERNEL; i++) {
                std::string cu_id = std::to_string(i + 1);
                std::string krnl_name_full =
                        krnl_name + ":{" + "xdrop_" + cu_id + "}";

                // printf("Creating a kernel [%s] for CU(%d)\n",
                             // krnl_name_full.c_str(),
                             // i + 1);

                //Here Kernel object is created by specifying kernel name along with compute unit.
                //For such case, this kernel object can only access the specific Compute unit

                OCL_CHECK(err,
                krnls[i] = cl::Kernel(
                program, krnl_name_full.c_str(), &err));
            }
            valid_device++;
            break; // we break because we found a valid device
        }
    }
    if (valid_device == 0) {
        std::cout << "Failed to program any device found, exit!\n";
        exit(EXIT_FAILURE);
    }

	//SETUP DATA TRANSFER
	std::vector<cl_mem_ext_ptr_t> inBufExt1(NUM_KERNEL);
    std::vector<cl_mem_ext_ptr_t> inBufExt2(NUM_KERNEL);
    std::vector<cl_mem_ext_ptr_t> inBufExt3(NUM_KERNEL);
    std::vector<cl_mem_ext_ptr_t> inBufExt4(NUM_KERNEL);
    std::vector<cl_mem_ext_ptr_t> inBufExt5(NUM_KERNEL);
    std::vector<cl_mem_ext_ptr_t> inBufExt6(NUM_KERNEL);
    std::vector<cl_mem_ext_ptr_t> outBufExt(NUM_KERNEL);

    std::vector<cl::Buffer> buffer_input1(NUM_KERNEL);
    std::vector<cl::Buffer> buffer_input2(NUM_KERNEL);
    std::vector<cl::Buffer> buffer_input3(NUM_KERNEL);
    std::vector<cl::Buffer> buffer_input4(NUM_KERNEL);
    std::vector<cl::Buffer> buffer_input5(NUM_KERNEL);
    std::vector<cl::Buffer> buffer_input6(NUM_KERNEL);
    std::vector<cl::Buffer> buffer_output(NUM_KERNEL);

    for (int i = 0; i < NUM_KERNEL; i++) {

    	if(i==0){
	        inBufExt1[i].obj = seeds_l.data();
	        inBufExt1[i].param = 0;
	        inBufExt1[i].flags = bank[i * 7];

	        inBufExt2[i].obj = prefQ.data();
	        inBufExt2[i].param = 0;
	        inBufExt2[i].flags = bank[(i * 7) + 1];

	        inBufExt3[i].obj = prefT.data();
	        inBufExt3[i].param = 0;
	        inBufExt3[i].flags = bank[(i * 7) + 2];

	        inBufExt4[i].obj = offsetLeftQ.data();
	        inBufExt4[i].param = 0;
	        inBufExt4[i].flags = bank[(i * 7) + 3];

	        inBufExt5[i].obj = offsetLeftT.data();
	        inBufExt5[i].param = 0;
	        inBufExt5[i].flags = bank[(i * 7) + 4];

	        inBufExt6[i].obj = ant_l.data();
	        inBufExt6[i].param = 0;
	       	inBufExt6[i].flags = bank[(i * 7) + 5];

	       	outBufExt[i].obj = scoreLeft.data();     
	        outBufExt[i].param = 0;
	        outBufExt[i].flags = bank[(i * 7) + 6];
    	}
    	else{
    		inBufExt1[i].obj = seeds_r.data();
	        inBufExt1[i].param = 0;
	        inBufExt1[i].flags = bank[i * 7];

	        inBufExt2[i].obj = suffQ.data();
	        inBufExt2[i].param = 0;
	        inBufExt2[i].flags = bank[(i * 7) + 1];

	        inBufExt3[i].obj = suffT.data();
	        inBufExt3[i].param = 0;
	        inBufExt3[i].flags = bank[(i * 7) + 2];

	        inBufExt4[i].obj = offsetRightQ.data();
	        inBufExt4[i].param = 0;
	        inBufExt4[i].flags = bank[(i * 7) + 3];

	        inBufExt5[i].obj = offsetRightT.data();
	        inBufExt5[i].param = 0;
	        inBufExt5[i].flags = bank[(i * 7) + 4];

	        inBufExt6[i].obj = ant_r.data();
	        inBufExt6[i].param = 0;
	       	inBufExt6[i].flags = bank[(i * 7) + 5];

	       	outBufExt[i].obj = scoreRight.data();     
	        outBufExt[i].param = 0;
	        outBufExt[i].flags = bank[(i * 7) + 6];
    	}

    }

	// These commands will allocate memory on the FPGA. The cl::Buffer objects can
    // be used to reference the memory locations on the device.
    //Creating Buffers
    for (int i = 0; i < NUM_KERNEL; i++) {
	        OCL_CHECK(err,
	                    buffer_input1[i] =            
	                    cl::Buffer(context,
	                    CL_MEM_READ_WRITE | CL_MEM_EXT_PTR_XILINX |
	                    CL_MEM_USE_HOST_PTR,
	                    dim*sizeof(SeedL),
	                    &inBufExt1[i],
	                    &err));
	        if(i==0){
		        OCL_CHECK(err,
		                    buffer_input2[i] =
		                    cl::Buffer(context,
		                    CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX |
		                    CL_MEM_USE_HOST_PTR,
		                    totalLengthQPref*sizeof(char),
		                    &inBufExt2[i],
		                    &err));
		        OCL_CHECK(err,
		                    buffer_input3[i] =
		                    cl::Buffer(context,
		                    CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX |
		                    CL_MEM_USE_HOST_PTR,
		                    totalLengthTPref*sizeof(char),
		                    &inBufExt3[i],
		                    &err));
		        OCL_CHECK(err,
		                    buffer_input4[i] =
		                    cl::Buffer(context,
		                    CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX |
		                    CL_MEM_USE_HOST_PTR,
		                    sizeof(int) * dim,
		                    &inBufExt4[i],
		                    &err));
		        OCL_CHECK(err,
		                    buffer_input5[i] =
		                    cl::Buffer(context,
		                    CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX |
		                    CL_MEM_USE_HOST_PTR,
		                    sizeof(int) * dim,
		                    &inBufExt5[i],
		                    &err));
		        OCL_CHECK(err,
		                    buffer_input6[i] =
		                    cl::Buffer(context,
		                    CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX |
		                    CL_MEM_USE_HOST_PTR,
		                    sizeof(short)*ant_len_left*3*dim,
		                    &inBufExt6[i],
		                    &err));
		        OCL_CHECK(err,
		                    buffer_output[i] =
		                    cl::Buffer(context,
		                    CL_MEM_WRITE_ONLY | CL_MEM_EXT_PTR_XILINX |
		                    CL_MEM_USE_HOST_PTR,
		                    sizeof(short) * dim,
		                    &outBufExt[i],
		                    &err));
    	}else{
    		OCL_CHECK(err,
	                    buffer_input2[i] =
	                    cl::Buffer(context,
	                    CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX |
	                    CL_MEM_USE_HOST_PTR,
	                    totalLengthQSuff*sizeof(char),
	                    &inBufExt2[i],
	                    &err));
	        OCL_CHECK(err,
	                    buffer_input3[i] =
	                    cl::Buffer(context,
	                    CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX |
	                    CL_MEM_USE_HOST_PTR,
	                    totalLengthTSuff*sizeof(char),
	                    &inBufExt3[i],
	                    &err));
	        OCL_CHECK(err,
	                    buffer_input4[i] =
	                    cl::Buffer(context,
	                    CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX |
	                    CL_MEM_USE_HOST_PTR,
	                    sizeof(int) * dim,
	                    &inBufExt4[i],
	                    &err));
	        OCL_CHECK(err,
	                    buffer_input5[i] =
	                    cl::Buffer(context,
	                    CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX |
	                    CL_MEM_USE_HOST_PTR,
	                    sizeof(int) * dim,
	                    &inBufExt5[i],
	                    &err));
	        OCL_CHECK(err,
	                    buffer_input6[i] =
	                    cl::Buffer(context,
	                    CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX |
	                    CL_MEM_USE_HOST_PTR,
	                    sizeof(short)*ant_len_right*3*dim,
	                    &inBufExt6[i],
	                    &err));
	        OCL_CHECK(err,
	                    buffer_output[i] =
	                    cl::Buffer(context,
	                    CL_MEM_WRITE_ONLY | CL_MEM_EXT_PTR_XILINX |
	                    CL_MEM_USE_HOST_PTR,
	                    sizeof(short) * dim,
	                    &outBufExt[i],
	                    &err));

    	}
    }

    // Copy input data to Device Global Memory
    for (int i = 0; i < NUM_KERNEL; i++) {
        OCL_CHECK(err,
        err = q.enqueueMigrateMemObjects(
        {buffer_input1[i], buffer_input2[i], buffer_input3[i], buffer_input4[i], buffer_input5[i], buffer_input6[i]},
        0 /* 0 means from host*/));
    }
    q.finish();
    // double kernel_time_in_sec = 0, result = 0;
    // std::chrono::duration<double> kernel_time(0);
	// auto kernel_start = NOW;

    for (int i = 0; i < NUM_KERNEL; i++) {
    	//Setting kernel arguments
		OCL_CHECK(err, err = krnls[i].setArg(0, buffer_input1[i]));
		OCL_CHECK(err, err = krnls[i].setArg(1, buffer_input2[i]));
		OCL_CHECK(err, err = krnls[i].setArg(2, buffer_input3[i]));
		if(i == 0){
			OCL_CHECK(err, err = krnls[i].setArg(3, EXTEND_LEFTL)); 
		}
		else{
			OCL_CHECK(err, err = krnls[i].setArg(3, EXTEND_RIGHTL)); 
		}
		OCL_CHECK(err, err = krnls[i].setArg(4, XDrop));
		OCL_CHECK(err, err = krnls[i].setArg(5, dim));
		OCL_CHECK(err, err = krnls[i].setArg(6, buffer_input4[i]));  
		OCL_CHECK(err, err = krnls[i].setArg(7, buffer_input5[i]));
		if(i == 0){
			OCL_CHECK(err, err = krnls[i].setArg(8, ant_len_left)); 
		}
		else{
			OCL_CHECK(err, err = krnls[i].setArg(8, ant_len_right)); 
		}
		OCL_CHECK(err, err = krnls[i].setArg(9, buffer_input6[i]));
		OCL_CHECK(err, err = krnls[i].setArg(10, buffer_output[i]));
		//Invoking the kernel
        OCL_CHECK(err, err = q.enqueueTask(krnls[i]));	
	}

	q.finish();
    
    //auto kernel_end = NOW;	 
	// kernel_time = std::chrono::duration<double>(kernel_end - kernel_start);
    // kernel_time_in_sec = kernel_time.count();   
    // std::cout << "TIME = " << kernel_time_in_sec << "s" << std::endl;       

    for (int i = 0; i < NUM_KERNEL; i++) {
        OCL_CHECK(err,
            err = q.enqueueMigrateMemObjects(
            {buffer_input1[i],buffer_output[i]},
            CL_MIGRATE_MEM_OBJECT_HOST));
    }
    q.finish();
	
	for(int i = 0; i < numAlignments; i++){
		res[i] = scoreLeft[i]+scoreRight[i]+kmer_length;
		// std::cout<<"Align: "<<i<<" res: "<<scoreLeft[i] << " " << scoreRight[i]<<std::endl;
		setEndPositionH(seeds[i], getEndPositionH(seeds_r[i]));    
		setEndPositionV(seeds[i], getEndPositionV(seeds_r[i])); 
		//cout << res[i] <<endl;
	}
	
	// //host mem reset
	// free(scoreLeft);
	// free(scoreRight);
}

void computeXdrop(std::vector< std::vector<std::string> > &v, int mat, int mis, int gap, int kmerLen, int xdrop, int numpair, std::string binary_file)
{

	std::cout << "EXECUTING HOST CODE FOR PALADIN\n";
	//Result result(kmerLen);
	int n_align = v.size();
	//int result;
	//myinfo FPGAresult;
	vector<ScoringSchemeL> penalties(n_align);
	vector<int> posV(n_align);
	vector<int> posH(n_align);
	vector<string> seqV(n_align);
	vector<string> seqH(n_align);
	vector<SeedL> seeds(n_align);
	for(int i = 0; i < v.size(); i++){
		ScoringSchemeL tmp_sscheme(mat, mis, -1, gap);
		penalties[i]=tmp_sscheme;
		posV[i]=stoi(v[i][1]);
		posH[i]=stoi(v[i][3]);
		seqV[i]=v[i][0];
		seqH[i]=v[i][2];
		std::string strand = v[i][4];

		if(strand == "c"){
			std::transform(
				std::begin(seqH[i]),
				std::end(seqH[i]),
				std::begin(seqH[i]),
				dummycomplement);
				posH[i] = seqH[i].length()-posH[i]-kmerLen;
		}
		
		SeedL tmp_seed(posH[i], posV[i], kmerLen);
		seeds[i] = tmp_seed;
	}
	//seqan testbench
	seqan::Score<int, seqan::Simple> scoringScheme_s(mat, mis, -1, gap);
	cout<< "PERFORMING "<< numpair << " ALIGNMENTS"<<endl;
	int *scoreSeqAn;
	scoreSeqAn = (int *)malloc(sizeof(int)*numpair);
	std::cout << "STARTING CPU" << std::endl;
	std::chrono::duration<double>  diff_s;
	vector<seqan::Dna5String> seqV_s_arr(numpair);
	vector<seqan::Dna5String> seqH_s_arr(numpair);
	vector<TSeed> seed(numpair);
	for(int i = 0; i < numpair; i++){
		seqan::Dna5String seqV_s(seqV[i]);
		seqan::Dna5String seqH_s(seqH[i]);
		seqV_s_arr[i]=seqV_s;
		seqH_s_arr[i]=seqH_s;
		TSeed tmp(posH[i], posV[i], kmerLen);
		seed[i]=tmp;
	}
	auto start_s = std::chrono::high_resolution_clock::now();
	#pragma omp parallel for
	for(int i = 0; i < numpair; i++){
		scoreSeqAn[i] = seqan::extendSeed(seed[i], seqH_s_arr[i], seqV_s_arr[i], seqan::EXTEND_LEFT, scoringScheme_s, xdrop, seqan::GappedXDrop(), kmerLen);
	}
	auto end_s = std::chrono::high_resolution_clock::now();
	diff_s = end_s-start_s;
	cout << "SEQAN TIME:\t" <<  diff_s.count() <<endl;

	int* scoreFPGA; // assuming this is 'short* res' in the xdrop() function
	scoreFPGA = (int *)malloc(sizeof(int)*numpair);
	std::chrono::duration<double>  diff_l;
	std::cout << "STARTING FPGA" << std::endl;
	auto start_l = NOW;
	extendSeedFPGA(seeds, EXTEND_BOTHL, seqH, seqV, penalties, xdrop, kmerLen, scoreFPGA, numpair, binary_file);
	auto end_l = NOW;
	diff_l = end_l-start_l;

	cout << "FPGA TIME:\t" <<  diff_l.count() <<endl;
	
	cout << "SPEEDUP " << diff_s.count()/diff_l.count()<<"X"<< endl;
}

//=======================================================================
//
// Function call main
//
//=======================================================================

int main(int argc, char **argv)
{
	// add optlist library later
	ifstream input(argv[1]);       // file name with sequences and seed positions
	int kmerLen = atoi(argv[2]);   // kmerLen
	int xdrop = atoi(argv[3]);     // xdrop
	std::string binary_file = argv[4]; //xcl binary

	int mat = 1, mis = -1, gap = -1;	
	const char* filename =  (char*) malloc(20 * sizeof(char));
	//filename = temp.c_str();
	std::cout << "STARTING BENCHMARK" << std::endl;

	//setting up the fpga environment

	uint64_t numpair = std::count(std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>(), '\n');
	input.seekg(0, std::ios_base::beg);

	vector<std::string> entries;

	// read input file
	if(input)
		for (int i = 0; i < numpair; ++i){

			std::string line;
			std::getline(input, line);
			entries.push_back(line);
		}
	input.close();
		// compute pairwise alignments
	vector< vector<string> > v(numpair);
	for(uint64_t i = 0; i < numpair; i++){
		//int ithread = i;//omp_get_thread_num();
		vector<string> temp = split(entries[i], '\t');
		// format: seqV, posV, seqH, posH, strand -- GGGG: generate this input with BELLA
		v[i]=temp;
	}
	
	computeXdrop(v, mat, mis, gap, kmerLen, xdrop, numpair, binary_file);

	return 0;
}
