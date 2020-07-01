<p align="center">
  <img width="682" height="170" src="https://github.com/albertozeni/XdropXOHW/blob/master/media/PALADIN.png">
</p>

## Project Information

Team number: xohw20_405  	<br />
Project name: X-drop on FPGA		<br />
Date: 30/06/2020			<br />
Version of uploaded archive:1	<br />
													<br />
University name: Politecnico di Milano				<br />
Supervisor name: Marco Domenico Santambrogio		<br />
Supervisor e-mail: marco.santambrogio@polimi.it		<br />
Participant: Alberto Zeni							<br />	
Email: alberto.zeni@mail.polimi.it					<br />
Participant: Guido Walter Di Donato							<br />
Email: guidowalter.didonato@mail.polimi.it					<br />
Participant: Francesco Peverelli					<br />
Email: francesco1.peverelli@mail.polimi.it			<br />
													<br />
Board used: Xilinx Alveo U280     <br />
Vitis Version: 2019.2							<br />

## Project Description
<p align="justify">
Pairwise sequence alignment is one of the most computationally intensive kernels in genomic data analysis, accounting for more than 90% of the run time for key bioinformatics applications. This method is particularly expensive for third-generation sequences due to the high computational expense of analyzing these long read lengths (1Kb-1Mb). Given the quadratic overhead of exact pairwise algorithms such as Smith-Waterman, for long alignments the community primarily relies on approximate algorithms that search only for high-quality alignments and stop early when one is not found. In this work, we present the first FPGA implementation of the popular X-drop alignment algorithm, named PALADIN.
</p>

## Usage

The repo already includes a host and a bitstream for the Alveo U280.
Before executing PALADIN be sure so source XRT.

### Compilation
Before executing PALADIN be sure so source XRT and Vitis.
PALADIN requires Vitis 2019.2 and C++14. To build PALADIN simply type:
```
make all TARGET=hw
```
PALADIN has been written to run on the Xilinx Alveo U280.
PALADIN will an executable called **host** and a bitstream file for the Alveo U280 called **xdrop.xclbin**.

### Demo

To check everything works properly type:
```
./host inputs/example.txt 17 21 build_dir.hw.xilinx_u280_xdma_201920_3/xdrop.xclbin 
```
This command executes PALADIN on our example dataset with a k-mer length of 17, an X-drop value of 21.
If everything executes correctly you can start using PALADIN with any input, and any X-drop.

The command line inputs are:
```
./host [input] [k-mer-length] [X-drop] [bitstream]
```
The input format for this demo is:
```
[seqV] [posV] [seqH] [posH] [strand]
```
**Each line of the input contains a pair of sequences to align**: the query sequence (seqV), the starting position of the seed on the query sequence (posV), the target sequence (seqH), the starting position of the seed on the target sequence (posH), and the relative strand ("c" if on opposite strands, "n" otherwise). Tab separated.
