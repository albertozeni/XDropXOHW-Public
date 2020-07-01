#include "seed.hpp"
#include "score.hpp"
#include <iostream>

#define MIN -32768
#define BYTES_INT 4

//trying to see if the scoring scheme is a bottleneck in some way
#define MATCH     1
#define MISMATCH -1
#define GAP_EXT  -1
#define GAP_OPEN -1
#define UNDEF -32767
#define WARP_DIM 32 
#define UNROLL_FACTOR 64 //since we are using chars 
#define EXTEND_NONEL  0
#define EXTEND_LEFTL  1
#define EXTEND_RIGHTL 2
#define EXTEND_BOTHL  3
#define NOW std::chrono::high_resolution_clock::now()
//dim = number of sequences

typedef int ExtensionDirectionL;


void updateExtendedSeedL(SeedL &seed,
                    ExtensionDirectionL direction, //as there are only 4 directions we may consider even smaller data types
                    int cols,
                    int rows,
                    int lowerDiag,
                    int upperDiag)
{
    if (direction == EXTEND_LEFTL)
    {
        int beginDiag = seed.beginDiagonal;
        // Set lower and upper diagonals.
        
        if (getLowerDiagonal(seed) > beginDiag + lowerDiag)
            setLowerDiagonal(seed, beginDiag + lowerDiag);
        if (getUpperDiagonal(seed) < beginDiag + upperDiag)
            setUpperDiagonal(seed, beginDiag + upperDiag);

        // Set new start position of seed.
        seed.beginPositionH -= rows;
        seed.beginPositionV -= cols;
    } else {  // direction == EXTEND_RIGHTL
        // Set new lower and upper diagonals.
        int endDiag = seed.endDiagonal;
        if (getUpperDiagonal(seed) < endDiag - lowerDiag)
            setUpperDiagonal(seed, (endDiag - lowerDiag));
        if (getLowerDiagonal(seed) > (endDiag - upperDiag))
            setLowerDiagonal(seed, endDiag - upperDiag);

        // Set new end position of seed.
        seed.endPositionH += rows;
        seed.endPositionV += cols;
        
    }
}
void calcExtendedLowerDiag(int &lowerDiag,
              int const &minCol,
              int const &antiDiagNo)
{
    int minRow = antiDiagNo - minCol;
    if (minCol - minRow < lowerDiag)
        lowerDiag = minCol - minRow;
}

void calcExtendedUpperDiag(int &upperDiag,
              int const &maxCol,
              int const &antiDiagNo)
{
    int maxRow = antiDiagNo + 1 - maxCol;
    if (maxCol - 1 - maxRow > upperDiag)
        upperDiag = maxCol - 1 - maxRow;
}

short reduce_max(short *input, int dim, int n_cus){
    
    reduce_loop:for(int i = n_cus/2; i > 0; i>>=1){
        max_unroll:for(int cu = 0; cu < UNROLL_FACTOR; cu++){
            #pragma HLS unroll
            if(cu < i)
                input[cu] = (input[cu] > input[cu + i]) ? input[cu] : input[cu + i];
        }
    }
    return input[0];
}

void initAntiDiags(
               short *antiDiag1,
               short *antiDiag2,
               short *antiDiag3,
               int &a2size,
               int &a3size,
               int const dropOff,
               int const gapCost,
               int const undefined){

    a2size = 1;

    antiDiag2[0] = 0;

    a3size = 2;

    antiDiag3[0] = gapCost;
    antiDiag3[1] = gapCost;

}

void initAntiDiag3(short *antiDiag3,
                    int &a3size,
                    int const offset,
                    int const maxCol,
                    int const antiDiagNo,
                    int const minScore,
                    int const gapCost,
                    int const undefined)
{
    a3size = maxCol + 1 - offset;

    antiDiag3[0] = undefined;
    antiDiag3[maxCol - offset] = undefined;

    if (antiDiagNo * gapCost > minScore)
    {
        if (offset == 0) // init first column
            antiDiag3[0] = antiDiagNo * gapCost;
        if (antiDiagNo - maxCol == 0) // init first row
            antiDiag3[maxCol - offset] = antiDiagNo * gapCost;
    }
}

void computeAntidiag(short *antiDiag1,
						short *antiDiag2,
						short *antiDiag3,
						char* querySeg,
						char* databaseSeg,
						int best,
						int scoreDropOff,
						int cols,
						int rows,
						int minCol,
						int maxCol,
						int antiDiagNo,
						int offset1,
						int offset2,
						ExtensionDirectionL direction
						//int n_threads
                                    ){

    
    for(int i = 0; i < maxCol; i+=UNROLL_FACTOR){
    #pragma HLS pipeline II=1 rewind
        for(int j = 0; j < UNROLL_FACTOR; j++){
            #pragma HLS unroll
            int col = j + minCol + i;
            int queryPos, dbPos;
    
            queryPos = col - 1;
            dbPos = col + rows - antiDiagNo - 1;
            
            if(col < maxCol){
    
                int tmp = (antiDiag2[col-offset2] > antiDiag2[col-offset2-1]) ? antiDiag2[col-offset2] : antiDiag2[col-offset2-1];
    
                tmp += GAP_EXT;

                int score = (querySeg[queryPos] == databaseSeg[dbPos]) ? MATCH : MISMATCH;
                

                tmp = (antiDiag1[col-offset1-1]+score > tmp) ? antiDiag1[col-offset1-1]+score : tmp;
                
                antiDiag3[j+1+i] = (tmp < best - scoreDropOff) ? UNDEF : tmp;
            }
        }
    }
}

extern "C" {
    void xdrop(SeedL *seed, char *querySegArray, char *databaseSegArray, ExtensionDirectionL direction, 
        int scoreDropOff, int numberOfAlignments, int *offsetQuery, int *offsetTarget, int offAntidiag,
        short *antidiag, short *res) {

        #pragma HLS interface m_axi port=seed offset=slave bundle=gmem0
        #pragma HLS interface m_axi port=querySegArray offset=slave bundle=gmem1
        #pragma HLS interface m_axi port=databaseSegArray offset=slave bundle=gmem2
        #pragma HLS interface m_axi port=offsetQuery offset=slave bundle=gmem3
        #pragma HLS interface m_axi port=offsetTarget offset=slave bundle=gmem4
        #pragma HLS interface m_axi port=antidiag offset=slave bundle=gmem5
        #pragma HLS interface m_axi port=res offset=slave bundle=gmem6

        #pragma HLS interface s_axilite port=seed bundle=control
        #pragma HLS interface s_axilite port=querySegArray bundle=control
        #pragma HLS interface s_axilite port=databaseSegArray bundle=control
        #pragma HLS interface s_axilite port=direction bundle=control
        #pragma HLS interface s_axilite port=scoreDropOff bundle=control
        #pragma HLS interface s_axilite port=numberOfAlignments bundle=control
        #pragma HLS interface s_axilite port=offsetQuery bundle=control
        #pragma HLS interface s_axilite port=offsetTarget bundle=control
        #pragma HLS interface s_axilite port=offAntidiag bundle=control
        #pragma HLS interface s_axilite port=antidiag bundle=control
        #pragma HLS interface s_axilite port=res bundle=control

        #pragma HLS interface s_axilite port=return bundle=control

        #pragma HLS DATA_PACK variable = seed //? forse non so se serve

        #pragma HLS INTERFACE ap_ctrl_chain port=return bundle=control
        #pragma HLS dataflow
        #pragma HLS inline recursive

        for(int alignmentId = 0; alignmentId < numberOfAlignments; alignmentId++){

            char *querySeg;
            char *databaseSeg;

            if(alignmentId==0){
                querySeg = querySegArray;
                databaseSeg = databaseSegArray;
            }
            else{
                querySeg = querySegArray + offsetQuery[alignmentId-1];
                databaseSeg = databaseSegArray + offsetTarget[alignmentId-1];
            }

            // for(int i = 0; i < 128; i++){
            //     std::cout << querySeg[i];
            // } std::cout << std::endl;

            // for(int i = 0; i < 128; i++){
            //     std::cout << databaseSeg[i];
            // } std::cout << std::endl;

            short *antiDiag1 = &antidiag[alignmentId*offAntidiag*3]; 
            short *antiDiag2 = &antiDiag1[offAntidiag];
            short *antiDiag3 = &antiDiag2[offAntidiag];

            SeedL mySeed(seed[alignmentId]);

            int a1size = 0, a2size = 0, a3size = 0;
            int cols, rows;

            if(alignmentId == 0){
            
                cols = offsetQuery[alignmentId]+1;
                rows = offsetTarget[alignmentId]+1;
            
            }else{
            
                cols = offsetQuery[alignmentId]-offsetQuery[alignmentId-1]+1;
                rows = offsetTarget[alignmentId]-offsetTarget[alignmentId-1]+1;
            
            }

            if (rows == 1 || cols == 1)
                return; // if the alignment is too short return

            int minCol = 1;
            int maxCol = 2;

            int offset1 = 0; // number of leading columns that need not be calculated in antiDiag1
            int offset2 = 0; // in antiDiag2
            int offset3 = 0; // in antiDiag3

            initAntiDiags(antiDiag1,antiDiag2, antiDiag3, a2size, a3size, scoreDropOff, GAP_EXT, UNDEF);
            int antiDiagNo = 1; // the currently calculated anti-diagonal

            int best = 0; // maximal score value in the DP matrix (for drop-off calculation)

            int lowerDiag = 0;
            int upperDiag = 0;

            short temp[64];

            while (minCol < maxCol){

                ++antiDiagNo;
                //antidiagswap
                //antiDiag2 -> antiDiag1
                //antiDiag3 -> antiDiag2
                //antiDiag1 -> antiDiag3
                short *t = antiDiag1;
                antiDiag1 = antiDiag2;
                antiDiag2 = antiDiag3;
                antiDiag3 = t;
                int t_l = a1size;
                a1size = a2size;
                a2size = a3size;
                a3size = t_l;
                offset1 = offset2;
                offset2 = offset3;
                offset3 = minCol-1;

                initAntiDiag3(antiDiag3, a3size, offset3, maxCol, antiDiagNo, best - scoreDropOff, GAP_EXT, UNDEF);
            
                computeAntidiag(antiDiag1, antiDiag2, antiDiag3, querySeg, databaseSeg, best, scoreDropOff, cols, rows, minCol, maxCol, antiDiagNo, offset1, offset2, direction);
            
                int tmp, antiDiagBest = UNDEF;  

                // for(int i = 0; i < a3size; i ++){
                //     antiDiagBest = (antiDiag3[i] > antiDiagBest) ? antiDiag3[i] : antiDiagBest;
                // }
                max:for(int i=0; i<a3size; i+=UNROLL_FACTOR){

                    int size = a3size-i;
                    
                    for(int j = 0; j < UNROLL_FACTOR; j++){ //might be able to unroll everything
                    #pragma HLS unroll
                        temp[j] = (j<size) ? antiDiag3[j+i]:UNDEF;              
                    }
                    
                    tmp = reduce_max(temp,size, UNROLL_FACTOR);
                    antiDiagBest = (tmp>antiDiagBest) ? tmp:antiDiagBest;

                }

                best = (best > antiDiagBest) ? best : antiDiagBest;
                
                mincol:while (minCol - offset3 < a3size && antiDiag3[minCol - offset3] == UNDEF &&
                       minCol - offset2 - 1 < a2size && antiDiag2[minCol - offset2 - 1] == UNDEF){
                    #pragma HLS pipeline
                    ++minCol;
                }

                // Calculate new maxCol
                maxcol:while (maxCol - offset3 > 0 && (antiDiag3[maxCol - offset3 - 1] == UNDEF) &&
                                               (antiDiag2[maxCol - offset2 - 1] == UNDEF)){
                    #pragma HLS pipeline
                    --maxCol;
                }
                ++maxCol;

                // Calculate new lowerDiag and upperDiag of extended seed
                calcExtendedLowerDiag(lowerDiag, minCol, antiDiagNo);
                calcExtendedUpperDiag(upperDiag, maxCol - 1, antiDiagNo);
                
                // end of databaseSeg reached?
                minCol = (minCol > (antiDiagNo + 2 - rows)) ? minCol : (antiDiagNo + 2 - rows);
                // end of querySeg reached?
                maxCol = (maxCol < cols) ? maxCol : cols;
                //}

            }
            int longestExtensionCol = a3size + offset3 - 2;
            int longestExtensionRow = antiDiagNo - longestExtensionCol;
            int longestExtensionScore = antiDiag3[longestExtensionCol - offset3];
            
            if (longestExtensionScore == UNDEF){
                if (antiDiag2[a2size -2] != UNDEF){
                    // reached end of query segment
                    longestExtensionCol = a2size + offset2 - 2;
                    longestExtensionRow = antiDiagNo - 1 - longestExtensionCol;
                    longestExtensionScore = antiDiag2[longestExtensionCol - offset2];
                    
                }
                else if (a2size > 2 && antiDiag2[a2size-3] != UNDEF){
                    // reached end of database segment
                    longestExtensionCol = a2size + offset2 - 3;
                    longestExtensionRow = antiDiagNo - 1 - longestExtensionCol;
                    longestExtensionScore = antiDiag2[longestExtensionCol - offset2];
                    
                }
            }


            //could be parallelized in some way
            if (longestExtensionScore == UNDEF){

                // general case
                for (int i = 0; i < a1size; ++i){
                    
                    if (antiDiag1[i] > longestExtensionScore){

                        longestExtensionScore = antiDiag1[i];
                        longestExtensionCol = i + offset1;
                        longestExtensionRow = antiDiagNo - 2 - longestExtensionCol;
                    
                    }
                }
            }
            if (longestExtensionScore != UNDEF)
                updateExtendedSeedL(mySeed, direction, longestExtensionCol, longestExtensionRow, lowerDiag, upperDiag);
            
            seed[alignmentId] = mySeed;
            res[alignmentId] = longestExtensionScore;
        }
    }
}
